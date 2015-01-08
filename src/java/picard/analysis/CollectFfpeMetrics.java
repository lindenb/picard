package picard.analysis;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.NotPrimaryAlignmentFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.util.DbSnpBitSetUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static java.lang.Math.log10;

/**
 * TODO javadoc
 */
@CommandLineProgramProperties(
        usage = CollectFfpeMetrics.USAGE,
        usageShort = CollectFfpeMetrics.USAGE,
        programGroup = Metrics.class
)
public class CollectFfpeMetrics extends CommandLineProgram {
    static final String USAGE = "TODO";

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Input BAM file for analysis.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Base path of output files to write.")
    public File OUTPUT;

    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME,
            doc = "Reference sequence to which BAM is aligned.")
    public File REFERENCE_SEQUENCE;

    @Option(doc = "An optional list of intervals to restrict analysis to.",
            optional = true)
    public File INTERVALS;

    @Option(doc = "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis.",
            optional = true)
    public File DB_SNP;

    @Option(shortName = "Q",
            doc = "The minimum base quality score for a base to be included in analysis.")
    public int MINIMUM_QUALITY_SCORE = 20;

    @Option(shortName = "MQ",
            doc = "The minimum mapping quality score for a base to be included in analysis.")
    public int MINIMUM_MAPPING_QUALITY = 30;

    @Option(shortName = "MIN_INS",
            doc = "The minimum insert size for a read to be included in analysis. Set of 0 to allow unpaired reads.")
    public int MINIMUM_INSERT_SIZE = 60;

    @Option(shortName = "MAX_INS",
            doc = "The maximum insert size for a read to be included in analysis. Set of 0 to allow unpaired reads.")
    public int MAXIMUM_INSERT_SIZE = 600;

    @Option(doc = "When available, use original quality scores for filtering.")
    public boolean USE_OQ = true;

    @Option(doc = "The number of context bases to include on each side of the assayed base.")
    public int CONTEXT_SIZE = 1;

    @Option(doc = "The optional set of sequence contexts to restrict analysis to. Their reverse complements will also be analyzed, " +
            "even if not specified here. If not supplied all contexts are analyzed.")
    public Set<String> CONTEXTS = new HashSet<String>();

    @Option(doc = "For debugging purposes: stop after visiting this many sites with at least 1X coverage.")
    public int STOP_AFTER = Integer.MAX_VALUE;

    private final Log log = Log.getInstance(CollectFfpeMetrics.class);
    private static final String UNKNOWN_LIBRARY = "UnknownLibrary";
    private static final String UNKNOWN_SAMPLE = "UnknownSample";
    private static final String ALL_CONTEXTS = "***";
    private static final byte[] BASES = {'A', 'C', 'G', 'T'};

    public static class FfpeSummaryMetrics extends MetricBase {
        public String SAMPLE_ALIAS;
        public String LIBRARY;

        public double A_TO_C_QSCORE;
        public double A_TO_G_QSCORE;
        public double A_TO_T_QSCORE;

        public double C_TO_A_QSCORE;
        public double C_TO_G_QSCORE;
        public double C_TO_T_QSCORE;

        public double G_TO_A_QSCORE;
        public double G_TO_C_QSCORE;
        public double G_TO_T_QSCORE;

        public double T_TO_A_QSCORE;
        public double T_TO_C_QSCORE;
        public double T_TO_G_QSCORE;
    }

    /**
     * Metrics for a specific mutation type, possibly limited to a specific context.
     */
    public static class FfpeDetailMetrics extends MetricBase {
        public String SAMPLE_ALIAS;
        public String LIBRARY;

        public byte REF_BASE;
        public byte ALT_BASE;

        public String CONTEXT;

        public long FWD_REF_REF_BASES;
        public long FWD_REF_ALT_BASES;
        public long REV_REF_REF_BASES;
        public long REV_REF_ALT_BASES;

        public double ERROR_RATE;
        public double QSCORE;
        // TODO p-value?

        /**
         * Compute error rate and Phred-scaled quality score given the observed counts.
         * Due to the symmetry of the problem, half of the error rates will be negative.
         * Instead of reporting these as-is, we effectively set them to zero, suggesting
         * that there is no evidence for that particular artifact.
         *
         * What's really happening is that each artifact is confounded by the presence
         * of another (its reverse complement). However from a chemistry perspective it
         * seems safe to assume that both artifacts won't occur significantly within the
         * same prepped library.
         */
        private void calculateDerivedStatistics() {
            // TODO fix divide by zero edge case
            final double rawErrorRate = (this.FWD_REF_ALT_BASES - this.REV_REF_ALT_BASES)
                    / (double) (this.FWD_REF_REF_BASES + this.FWD_REF_ALT_BASES + this.REV_REF_REF_BASES + this.REV_REF_ALT_BASES);
            this.ERROR_RATE = Math.max(rawErrorRate, 1e-10);
            this.QSCORE = -10 * log10(this.ERROR_RATE);
        }
    }

    /**
     * SAM filter for insert size range.
     */
    static class InsertSizeFilter implements SamRecordFilter {
        final int minInsertSize;
        final int maxInsertSize;

        InsertSizeFilter(final int minInsertSize, final int maxInsertSize) {
            this.minInsertSize = minInsertSize;
            this.maxInsertSize = maxInsertSize;
        }

        @Override
        public boolean filterOut(final SAMRecord rec) {
            // Treat both parameters == 0 as not filtering
            if (minInsertSize == 0 && maxInsertSize == 0) return false;

            if (rec.getReadPairedFlag()) {
                final int ins = Math.abs(rec.getInferredInsertSize());
                return ins < minInsertSize || ins > maxInsertSize;
            }

            // If the read isn't paired and either min or max is specified filter it out
            return minInsertSize != 0 || maxInsertSize != 0;
        }

        @Override
        public boolean filterOut(final SAMRecord r1, final SAMRecord r2) {
            return filterOut(r1) || filterOut(r2);
        }
    }

    // Stock main method
    public static void main(final String[] args) {
        new CollectFfpeMetrics().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        final int size = 1 + 2 * CONTEXT_SIZE;
        final List<String> messages = new ArrayList<String>();

        for (final String ctx : CONTEXTS) {
            if (ctx.length() != size) {
                messages.add("Context " + ctx + " is not " + size + " long as implied by CONTEXT_SIZE=" + CONTEXT_SIZE);
            }
        }

        if (CONTEXT_SIZE < 0) messages.add("CONTEXT_SIZE cannot be negative");

        return messages.isEmpty() ? null : messages.toArray(new String[messages.size()]);
    }

    /** Mimic of Oracle's nvl() - returns the first value if not null, otherwise the second value. */
    private <T> T nvl(final T value1, final T value2) {
        if (value1 != null) return value1;
        else return value2;
    }

    @Override
    protected int doWork() {
        final File SUMMARY_OUT = new File(OUTPUT + ".ffpe_summary_metrics");
        final File DETAILS_OUT = new File(OUTPUT + ".ffpe_detail_metrics");

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(SUMMARY_OUT);
        IOUtil.assertFileIsWritable(DETAILS_OUT);

        if (INTERVALS != null) IOUtil.assertFileIsReadable(INTERVALS);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
        final SamReader in = SamReaderFactory.makeDefault().open(INPUT);

        final Set<String> samples = new HashSet<String>();
        final Set<String> libraries = new HashSet<String>();
        for (final SAMReadGroupRecord rec : in.getFileHeader().getReadGroups()) {
            samples.add(nvl(rec.getSample(), UNKNOWN_SAMPLE));
            libraries.add(nvl(rec.getLibrary(), UNKNOWN_LIBRARY));
        }

        // Setup the calculator
        final int contextLength = (2 * CONTEXT_SIZE) + 1;
        final Set<String> contexts = CONTEXTS.isEmpty() ? makeContextStrings(contextLength) : includeReverseComplements(CONTEXTS);
        final FfpeCalculator counts = new FfpeCalculator(libraries, contexts);

        // Load up dbSNP if available
        log.info("Loading dbSNP File: " + DB_SNP);
        final DbSnpBitSetUtil dbSnp;
        if (DB_SNP != null) dbSnp = new DbSnpBitSetUtil(DB_SNP, in.getFileHeader().getSequenceDictionary());
        else dbSnp = null;

        // Make an iterator that will filter out funny looking things
        final SamLocusIterator iterator;
        if (INTERVALS != null) {
            final IntervalList intervals = IntervalList.fromFile(INTERVALS);
            iterator = new SamLocusIterator(in, intervals.uniqued(), false);
        } else {
            iterator = new SamLocusIterator(in);
        }
        iterator.setEmitUncoveredLoci(false);
        iterator.setMappingQualityScoreCutoff(MINIMUM_MAPPING_QUALITY);
        iterator.setSamFilters(Arrays.asList(
                new NotPrimaryAlignmentFilter(),
                new DuplicateReadFilter(),
                new InsertSizeFilter(MINIMUM_INSERT_SIZE, MAXIMUM_INSERT_SIZE)
        ));

        log.info("Starting iteration.");
        long nextLogTime = 0;
        int sites = 0;

        // Iterate over reference loci
        for (final SamLocusIterator.LocusInfo info : iterator) {
            // Skip dbSNP sites
            final String chrom = info.getSequenceName();
            final int pos = info.getPosition();
            final int index = pos - 1;
            if (dbSnp != null && dbSnp.isDbSnpSite(chrom, pos)) continue;

            // Skip sites at the end of chromosomes
            final byte[] bases = refWalker.get(info.getSequenceIndex()).getBases();
            if (pos < contextLength || pos > bases.length - contextLength) continue;

            // Get the reference context string and perform counting
            final String context = StringUtil.bytesToString(bases, index - CONTEXT_SIZE, contextLength).toUpperCase();
            if (contexts.contains(context)) counts.countAlleles(info, context);

            // See if we need to stop
            if (++sites % 100 == 0) {
                final long now = System.currentTimeMillis();
                if (now > nextLogTime) {
                    log.info("Visited " + sites + " sites of interest. Last site: " + chrom + ":" + pos);
                    nextLogTime = now + 60000;
                }
            }
            if (sites >= STOP_AFTER) break;
        }

        // Finish up and write metrics
        final MetricsFile<FfpeSummaryMetrics, Integer> summaryMetricsFile = getMetricsFile();
        final MetricsFile<FfpeDetailMetrics, Integer> detailMetricsFile = getMetricsFile();
        final String sampleAlias = StringUtil.join(",", new ArrayList<String>(samples));
        final List<LibraryLevelMetrics> allMetrics = counts.finalize(sampleAlias);
        for (final LibraryLevelMetrics llm : allMetrics) {
            for (final FfpeDetailMetrics clm : llm.contextLevelMetrics) detailMetricsFile.addMetric(clm);
            for (final FfpeDetailMetrics alm : llm.artifactLevelMetrics) detailMetricsFile.addMetric(alm);
            summaryMetricsFile.addMetric(llm.summaryMetrics);
        }
        summaryMetricsFile.write(SUMMARY_OUT);
        detailMetricsFile.write(DETAILS_OUT);
        CloserUtil.close(in);
        return 0;
    }

    /**
     * Little method to expand a set of sequences to include all reverse complements.
     * Necessary if a user manually specifies a set of contexts, due to symmetries that the analysis depends on.
     */
    private Set<String> includeReverseComplements(final Set<String> sequences) {
        final Set<String> all = new HashSet<String>();
        for (final String seq : sequences) {
            all.add(seq);
            all.add(SequenceUtil.reverseComplement(seq));
        }
        return all;
    }

    private Set<String> makeContextStrings(final int length) {
        final Set<String> contexts = new HashSet<String>();

        for (final byte[] kmer : generateAllKmers(length)) {
           contexts.add(StringUtil.bytesToString(kmer));
        }

        log.info("Generated " + contexts.size() + " context strings.");
        return contexts;
    }

    /** Generates all possible unambiguous kmers of length and returns them as byte[]s. */
    private List<byte[]> generateAllKmers(final int length) {
        final List<byte[]> sofar = new LinkedList<byte[]>();

        if (sofar.size() == 0) {
            sofar.add(new byte[length]);
        }

        while (true) {
            final byte[] bs = sofar.remove(0);
            final int indexOfNextBase = findIndexOfNextBase(bs);

            if (indexOfNextBase == -1) {
                sofar.add(bs);
                break;
            } else {
                for (final byte b : BASES) {
                    final byte[] next = Arrays.copyOf(bs, bs.length);
                    next[indexOfNextBase] = b;
                    sofar.add(next);
                }
            }
        }

        return sofar;
    }

    /** Finds the first non-zero character in the array, or returns -1 if all are non-zero. */
    private int findIndexOfNextBase(final byte[] bs) {
        for (int i = 0; i < bs.length; ++i) {
            if (bs[i] == 0) return i;
        }

        return -1;
    }

    /** A container for metrics at the per-library level. */
    private class LibraryLevelMetrics {
        protected final List<FfpeDetailMetrics> contextLevelMetrics;
        protected final List<FfpeDetailMetrics> artifactLevelMetrics;
        protected final FfpeSummaryMetrics summaryMetrics;

        private LibraryLevelMetrics(final List<FfpeDetailMetrics> contextLevelMetrics,
                                    final List<FfpeDetailMetrics> artifactLevelMetrics,
                                    final FfpeSummaryMetrics summaryMetrics) {
            this.contextLevelMetrics = contextLevelMetrics;
            this.artifactLevelMetrics = artifactLevelMetrics;
            this.summaryMetrics = summaryMetrics;
        }
    }

    private class AlleleCounter {
        private final Map<String, Map<Byte, Long>> alleleCountsPerContext;

        private AlleleCounter(final Set<String> contexts) {
            // populate empty counts
            this.alleleCountsPerContext = new HashMap<String, Map<Byte, Long>>();
            for (final String context : contexts) {
                final Map<Byte, Long> baseCounts = new HashMap<Byte, Long>();
                for (final byte b : BASES) baseCounts.put(b, 0l);
                this.alleleCountsPerContext.put(context, baseCounts);
            }
        }

        private void addCallToContext(final String context, final byte calledBase) {
            final Map<Byte, Long> baseCounts = this.alleleCountsPerContext.get(context);
            baseCounts.put(calledBase, baseCounts.get(calledBase) + 1);
        }

        private LibraryLevelMetrics finalizeAndComputeMetrics(final String sampleAlias, final String library) {

            /**
             * 1. compute context-level metrics.
             */
            final List<FfpeDetailMetrics> contextLevelMetrics = new ArrayList<FfpeDetailMetrics>();
            for (final String context : this.alleleCountsPerContext.keySet()) {
                final byte refBase = (byte) context.charAt(CONTEXT_SIZE);
                for (final byte altBase : BASES) {
                    if (altBase != refBase) {
                        final FfpeDetailMetrics clm = new FfpeDetailMetrics();
                        clm.SAMPLE_ALIAS = sampleAlias;
                        clm.LIBRARY = library;
                        clm.CONTEXT = context;
                        clm.REF_BASE = refBase;
                        clm.ALT_BASE = altBase;

                        /**
                         * TODO explain what we're counting here
                         */
                        clm.FWD_REF_REF_BASES = this.alleleCountsPerContext.get(context).get(refBase);
                        clm.FWD_REF_ALT_BASES = this.alleleCountsPerContext.get(context).get(altBase);
                        clm.REV_REF_REF_BASES = this.alleleCountsPerContext.get(SequenceUtil.reverseComplement(context)).get(SequenceUtil.complement(refBase));
                        clm.REV_REF_ALT_BASES = this.alleleCountsPerContext.get(SequenceUtil.reverseComplement(context)).get(SequenceUtil.complement(altBase));

                        clm.calculateDerivedStatistics();
                        contextLevelMetrics.add(clm);
                    }
                }
            }

            /**
             * 2. collapse context-level metrics into artifact-level metrics.
             * TODO refactor this to be less cumbersome?
             */
            final CollectionUtil.TwoKeyHashMap<Byte, Byte, FfpeDetailMetrics> artifactLevelMetricsMap =
                    new CollectionUtil.TwoKeyHashMap<Byte, Byte, FfpeDetailMetrics>();
            for (final byte refBase : BASES) {
                for (final byte altBase : BASES) {
                    if (altBase != refBase) {
                        final FfpeDetailMetrics alm = new FfpeDetailMetrics();
                        alm.SAMPLE_ALIAS = sampleAlias;
                        alm.LIBRARY = library;
                        alm.CONTEXT = ALL_CONTEXTS;
                        alm.REF_BASE = refBase;
                        alm.ALT_BASE = altBase;
                        alm.FWD_REF_REF_BASES = 0;
                        alm.FWD_REF_ALT_BASES = 0;
                        alm.REV_REF_REF_BASES = 0;
                        alm.REV_REF_ALT_BASES = 0;
                        artifactLevelMetricsMap.put(refBase, altBase, alm);
                    }
                }
            }
            for (final FfpeDetailMetrics clm : contextLevelMetrics) {
                final FfpeDetailMetrics alm = artifactLevelMetricsMap.get(clm.REF_BASE, clm.ALT_BASE);
                alm.FWD_REF_REF_BASES += clm.FWD_REF_REF_BASES;
                alm.FWD_REF_ALT_BASES += clm.FWD_REF_ALT_BASES;
                alm.REV_REF_REF_BASES += clm.REV_REF_REF_BASES;
                alm.REV_REF_ALT_BASES += clm.REV_REF_ALT_BASES;
            }
            for (final FfpeDetailMetrics alm : artifactLevelMetricsMap.values()) {
                alm.calculateDerivedStatistics();
            }

            /**
             * 3. collapse artifact-level metrics into a single summary metric.
             */
            final FfpeSummaryMetrics summaryMetrics = new FfpeSummaryMetrics();
            summaryMetrics.SAMPLE_ALIAS = sampleAlias;
            summaryMetrics.LIBRARY = library;
            summaryMetrics.A_TO_C_QSCORE = artifactLevelMetricsMap.get((byte)'A', (byte)'C').QSCORE;
            summaryMetrics.A_TO_G_QSCORE = artifactLevelMetricsMap.get((byte)'A', (byte)'G').QSCORE;
            summaryMetrics.A_TO_T_QSCORE = artifactLevelMetricsMap.get((byte)'A', (byte)'T').QSCORE;
            summaryMetrics.C_TO_A_QSCORE = artifactLevelMetricsMap.get((byte)'C', (byte)'A').QSCORE;
            summaryMetrics.C_TO_G_QSCORE = artifactLevelMetricsMap.get((byte)'C', (byte)'G').QSCORE;
            summaryMetrics.C_TO_T_QSCORE = artifactLevelMetricsMap.get((byte)'C', (byte)'T').QSCORE;
            summaryMetrics.G_TO_A_QSCORE = artifactLevelMetricsMap.get((byte)'G', (byte)'A').QSCORE;
            summaryMetrics.G_TO_C_QSCORE = artifactLevelMetricsMap.get((byte)'G', (byte)'C').QSCORE;
            summaryMetrics.G_TO_T_QSCORE = artifactLevelMetricsMap.get((byte)'G', (byte)'T').QSCORE;
            summaryMetrics.T_TO_A_QSCORE = artifactLevelMetricsMap.get((byte)'T', (byte)'A').QSCORE;
            summaryMetrics.T_TO_C_QSCORE = artifactLevelMetricsMap.get((byte)'T', (byte)'C').QSCORE;
            summaryMetrics.T_TO_G_QSCORE = artifactLevelMetricsMap.get((byte)'T', (byte)'G').QSCORE;

            /** 4. bring them all and in the darkness bind them */
            final List<FfpeDetailMetrics> artifactLevelMetrics = new ArrayList<FfpeDetailMetrics>(artifactLevelMetricsMap.values());
            return new LibraryLevelMetrics(contextLevelMetrics, artifactLevelMetrics, summaryMetrics);
        }
    }

    /**
     *
     */
    private class FfpeCalculator {
        private final Set<String> acceptedContexts;
        private final Set<String> acceptedLibraries;
        private final Map<String, AlleleCounter> libraryMap;

        private FfpeCalculator(final Set<String> libraries, final Set<String> contexts) {
            this.acceptedLibraries = libraries;
            this.acceptedContexts = contexts;
            this.libraryMap = new HashMap<String, AlleleCounter>();
            for (final String library : libraries) {
                this.libraryMap.put(library, new AlleleCounter(contexts));
            }
        }

        private void countAlleles(final SamLocusIterator.LocusInfo info, final String refContext) {
            for (final SamLocusIterator.RecordAndOffset rec : info.getRecordAndPositions()) {
                final byte qual;
                final SAMRecord samrec = rec.getRecord();

                if (USE_OQ) {
                    final byte[] oqs = samrec.getOriginalBaseQualities();
                    if (oqs != null) qual = oqs[rec.getOffset()];
                    else qual = rec.getBaseQuality();
                } else {
                    qual = rec.getBaseQuality();
                }

                // Skip if below qual
                if (qual < MINIMUM_QUALITY_SCORE) continue;

                // Skip if context or library is unknown
                final String library = nvl(samrec.getReadGroup().getLibrary(), UNKNOWN_LIBRARY);
                if (!acceptedLibraries.contains(library)) continue;
                if (!acceptedContexts.contains(refContext)) continue;

                final String forwardTemplateContext;
                final byte forwardCalledBase;
                /**
                 * Remember that if the read is aligned to the negative strand, the reported ref + call bases must be reverse complemented
                 * to reflect the original sequence. Likewise if the read is second in a pair (and thus derived from the reverse template),
                 * the ref + call bases must be revc'd to reflect the forward template. If both of these conditions are true, they cancel
                 * each other out (hence the XOR).
                 */
                final boolean isNegativeStrand = samrec.getReadNegativeStrandFlag();
                final boolean isReadTwo = samrec.getReadPairedFlag() && samrec.getSecondOfPairFlag();
                if (isNegativeStrand ^ isReadTwo) {
                    forwardTemplateContext = SequenceUtil.reverseComplement(refContext);
                    forwardCalledBase = SequenceUtil.complement(rec.getReadBase());
                } else {
                    forwardTemplateContext = refContext;
                    forwardCalledBase = rec.getReadBase();
                }

                // Count the base
                this.libraryMap.get(library).addCallToContext(forwardTemplateContext, forwardCalledBase);
            }
        }

        private List<LibraryLevelMetrics> finalize(final String sampleAlias) {
            final List<LibraryLevelMetrics> allMetrics = new ArrayList<LibraryLevelMetrics>();
            for (final String library : this.libraryMap.keySet()) {
                allMetrics.add(this.libraryMap.get(library).finalizeAndComputeMetrics(sampleAlias, library));
            }
            return allMetrics;
        }
    }
}
