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
import htsjdk.samtools.util.ListMap;
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
    private static final byte[] BASES = {'A', 'C', 'G', 'T'};
    private static final double MIN_ERROR = 1e-10;

    private static class FfpeArtifact {
        // TODO
    }

    public static class FfpeSummaryMetrics extends MetricBase {
        public String SAMPLE_ALIAS;
        public String LIBRARY;

        // TODO include error rates?
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
     *
     * Here's an example: the characteristic OxoG artifact is that a G on the reference strand gets oxidated
     * and binds to an A, rather than a C. Thus the negative strand says 'A', and when it gets reverse transcribed
     * we get a 'T' where the original strand had a 'G'. The resulting alignment would call a G>T variant on read1,
     * or a C>A variant on read2 (or the opposite, if the read aligns to the negative strand of the reference).
     * The error rate is then computed as:
     *
     * ((G>T read 1 + C>A read 2) - (G>T read 2 + C>A read 1)) / (G>T read 1 + G>T read 2 + G>G read 1 + G>G read 2 + C>A read 1 + C>A read 2 + C>C read 1 + C>C read 2)
     * = (alt alignments supporting OxoG - alt alignments refuting OxoG) / (alt alignments relating to OxoG + corresponding ref alignments)
     *
     * The key here is that IF the damage occurs on the reference strand before PCR, read 1 will preferentially see one kind of variant
     * and read 2 will preferentially see the reverse complement of that variant (whereas actual mutations should occur equally on both strands,
     * so their contributions will cancel out in the calculation).
     *
     */
    public static class FfpeDetailMetrics extends MetricBase {
        public String SAMPLE_ALIAS;
        public String LIBRARY;
        public String CONTEXT;

        public char REF_BASE;
        public char ALT_BASE;

        public long PRO_REF_BASES;
        public long PRO_ALT_BASES;
        public long CON_REF_BASES;
        public long CON_ALT_BASES;

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
            this.ERROR_RATE = MIN_ERROR;
            final long totalBases = this.PRO_REF_BASES + this.PRO_ALT_BASES + this.CON_REF_BASES + this.CON_ALT_BASES;
            if (totalBases > 0) {
                final double rawErrorRate = (this.PRO_ALT_BASES - this.CON_ALT_BASES) / (double) totalBases;
                this.ERROR_RATE = Math.max(MIN_ERROR, rawErrorRate);
            }
            this.QSCORE = -10 * log10(this.ERROR_RATE);
        }
    }

    /**
     * TODO explain
     *
     * Continuing the OxoG example above:
     * (G>T read 1 + G>T read 2) / (G>T read 1 + G>T read 2 + G>G read 1 + G>G read 2) = G>T alignments / (G>T alignments + G>G alignments)
     * (C>A read 1 + C>A read 2) / (C>A read 1 + C>A read 2 + C>C read 1 + C>C read 2) = C>A alignments / (C>A alignments + C>C alignments)
     *
     */
    public static class ReferenceBiasMetrics extends MetricBase {
        public String SAMPLE_ALIAS;
        public String LIBRARY;
        public String CONTEXT;

        public char REF_BASE;
        public char ALT_BASE;

        public long FWD_CXT_REF_BASES;
        public long FWD_CXT_ALT_BASES;
        public long REV_CXT_REF_BASES;
        public long REV_CXT_ALT_BASES;

        double FWD_ERROR_RATE;
        double REV_ERROR_RATE;

        double REF_BIAS;

        /**
         * This calculation follows from the C_REF and G_REF error rates in CollectOxoGMetrics.
         */
        private void calculateDerivedStatistics() {
            this.FWD_ERROR_RATE = MIN_ERROR;
            final long fwdBases = this.FWD_CXT_REF_BASES + this.FWD_CXT_ALT_BASES;
            if (fwdBases > 0) {
                final double fwdErr = this.FWD_CXT_ALT_BASES / (double) fwdBases;
                this.FWD_ERROR_RATE = Math.max(MIN_ERROR, fwdErr);
            }

            this.REV_ERROR_RATE = MIN_ERROR;
            final long revBases = this.REV_CXT_REF_BASES + this.REV_CXT_ALT_BASES;
            if (revBases > 0) {
                final double revErr = this.REV_CXT_ALT_BASES / (double) revBases;
                this.REV_ERROR_RATE = Math.max(MIN_ERROR, revErr);
            }

            this.REF_BIAS = this.FWD_ERROR_RATE - this.REV_ERROR_RATE;
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
        final File REFBIAS_OUT = new File(OUTPUT + ".reference_bias_metrics");

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(SUMMARY_OUT);
        IOUtil.assertFileIsWritable(DETAILS_OUT);
        IOUtil.assertFileIsWritable(REFBIAS_OUT);

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
            // TODO what should the buffer be? In OxoG it's always 3, but that could crash if CONTEXT_SIZE is greater than 1
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
        final MetricsFile<ReferenceBiasMetrics, Integer> refBiasMetricsFile = getMetricsFile();

        final String sampleAlias = StringUtil.join(",", new ArrayList<String>(samples));
        final List<LibraryLevelMetrics> allMetrics = counts.finalize(sampleAlias);

        for (final LibraryLevelMetrics llm : allMetrics) {
            for (final FfpeDetailMetrics clm : llm.contextLevelMetrics) detailMetricsFile.addMetric(clm);
            for (final FfpeDetailMetrics alm : llm.artifactLevelMetrics) detailMetricsFile.addMetric(alm);
            for (final ReferenceBiasMetrics rbm : llm.referenceBiasMetrics) refBiasMetricsFile.addMetric(rbm);
            summaryMetricsFile.addMetric(llm.summaryMetrics);
        }

        summaryMetricsFile.write(SUMMARY_OUT);
        detailMetricsFile.write(DETAILS_OUT);
        refBiasMetricsFile.write(REFBIAS_OUT);
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

    /**
     * Breaks down alignments by read1/read2 and positive/negative strand.
     */
    private static class AlignmentAccumulator {
        private long R1_POS = 0;
        private long R1_NEG = 0;
        private long R2_POS = 0;
        private long R2_NEG = 0;

        private void addRecord(final SAMRecord rec) {
            final boolean isNegativeStrand = rec.getReadNegativeStrandFlag();
            final boolean isReadTwo = rec.getReadPairedFlag() && rec.getSecondOfPairFlag();

            if (isReadTwo) {
                if (isNegativeStrand) this.R2_NEG++;
                else this.R2_POS++;
            } else {
                if (isNegativeStrand) this.R1_NEG++;
                else this.R1_POS++;
            }
        }

        private static AlignmentAccumulator combine(final Iterable<AlignmentAccumulator> aas) {
            final AlignmentAccumulator combined = new AlignmentAccumulator();
            for (final AlignmentAccumulator aa : aas) {
                combined.R1_POS += aa.R1_POS;
                combined.R1_NEG += aa.R1_NEG;
                combined.R2_POS += aa.R2_POS;
                combined.R2_NEG += aa.R2_NEG;
            }

            return combined;
        }
    }

    private static class CalledBaseAccumulator extends HashMap<Byte, AlignmentAccumulator> {
        private CalledBaseAccumulator() {
            for (final byte b : BASES) {
                this.put(b, new AlignmentAccumulator());
            }
        }

        private void addRecord(final SAMRecord rec, final byte readBase) {
            this.get(readBase).addRecord(rec);
        }

        private static CalledBaseAccumulator combine(final Iterable<CalledBaseAccumulator> cbas) {
            final ListMap<Byte, AlignmentAccumulator> invertedCbas = new ListMap<Byte, AlignmentAccumulator>();
            for (final CalledBaseAccumulator cba : cbas) {
                for (final byte b : BASES) {
                    invertedCbas.add(b, cba.get(b));
                }
            }

            final CalledBaseAccumulator combined = new CalledBaseAccumulator();
            for (final byte b : BASES) {
                combined.put(b, AlignmentAccumulator.combine(invertedCbas.get(b)));
            }

            return combined;
        }
    }

    private static class ContextAccumulator extends HashMap<String, CalledBaseAccumulator> {
        private ContextAccumulator(final Iterable<String> allowedContexts) {
            for (final String context : allowedContexts) {
                this.put(context, new CalledBaseAccumulator());
            }
        }

        private void addRecord(final SAMRecord rec, final byte readBase, final String context) {
            this.get(context).addRecord(rec, readBase);
        }
    }

    /** A container for metrics at the per-library level. */
    private static class LibraryLevelMetrics {
        private final List<FfpeDetailMetrics> contextLevelMetrics;
        private final List<FfpeDetailMetrics> artifactLevelMetrics;
        private final List<ReferenceBiasMetrics> referenceBiasMetrics;
        private final FfpeSummaryMetrics summaryMetrics;

        private LibraryLevelMetrics(final List<FfpeDetailMetrics> contextLevelMetrics,
                                    final List<FfpeDetailMetrics> artifactLevelMetrics,
                                    final List<ReferenceBiasMetrics> referenceBiasMetrics,
                                    final FfpeSummaryMetrics summaryMetrics) {
            this.contextLevelMetrics = contextLevelMetrics;
            this.artifactLevelMetrics = artifactLevelMetrics;
            this.referenceBiasMetrics = referenceBiasMetrics;
            this.summaryMetrics = summaryMetrics;
        }
    }

    private LibraryLevelMetrics extractMetrics(final String sampleAlias, final String library, final ContextAccumulator contextAccumulator) {
        final List<FfpeDetailMetrics> contextLevelMetrics = new ArrayList<FfpeDetailMetrics>();
        final List<FfpeDetailMetrics> artifactLevelMetrics = new ArrayList<FfpeDetailMetrics>();
        final List<ReferenceBiasMetrics> referenceBiasMetrics = new ArrayList<ReferenceBiasMetrics>();
        final FfpeSummaryMetrics summaryMetrics = new FfpeSummaryMetrics();

        for (final String context : contextAccumulator.keySet()) {
            final byte refBase = (byte) context.charAt(CONTEXT_SIZE);
            for (final byte altBase : BASES) {
                if (altBase != refBase) {
                    final FfpeDetailMetrics detailMetric = new FfpeDetailMetrics();
                    final ReferenceBiasMetrics refBiasMetric = new ReferenceBiasMetrics();

                    // retrieve all the necessary alignment counters
                    final AlignmentAccumulator fwdRefAlignments = contextAccumulator.get(context).get(refBase);
                    final AlignmentAccumulator fwdAltAlignments = contextAccumulator.get(context).get(altBase);
                    final AlignmentAccumulator revRefAlignments = contextAccumulator.get(SequenceUtil.reverseComplement(context)).get(SequenceUtil.complement(refBase));
                    final AlignmentAccumulator revAltAlignments = contextAccumulator.get(SequenceUtil.reverseComplement(context)).get(SequenceUtil.complement(altBase));

                    // populate basic fields
                    detailMetric.SAMPLE_ALIAS = sampleAlias;
                    detailMetric.LIBRARY = library;
                    detailMetric.CONTEXT = context;
                    detailMetric.REF_BASE = (char) refBase;
                    detailMetric.ALT_BASE = (char) altBase;

                    refBiasMetric.SAMPLE_ALIAS = sampleAlias;
                    refBiasMetric.LIBRARY = library;
                    refBiasMetric.CONTEXT = context;
                    refBiasMetric.REF_BASE = (char) refBase;
                    refBiasMetric.ALT_BASE = (char) altBase;

                    // do the actual counting, as explained in the metric definitions
                    detailMetric.PRO_REF_BASES = fwdRefAlignments.R1_POS + fwdRefAlignments.R2_NEG + revRefAlignments.R1_NEG + revRefAlignments.R2_POS;
                    detailMetric.PRO_ALT_BASES = fwdAltAlignments.R1_POS + fwdAltAlignments.R2_NEG + revAltAlignments.R1_NEG + revAltAlignments.R2_POS;
                    detailMetric.CON_REF_BASES = fwdRefAlignments.R1_NEG + fwdRefAlignments.R2_POS + revRefAlignments.R1_POS + revRefAlignments.R2_NEG;
                    detailMetric.CON_ALT_BASES = fwdAltAlignments.R1_NEG + fwdAltAlignments.R2_POS + revAltAlignments.R1_POS + revAltAlignments.R2_NEG;

                    refBiasMetric.FWD_CXT_REF_BASES = fwdRefAlignments.R1_POS + fwdRefAlignments.R1_NEG + fwdRefAlignments.R2_POS + fwdRefAlignments.R2_NEG;
                    refBiasMetric.FWD_CXT_ALT_BASES = fwdAltAlignments.R1_POS + fwdAltAlignments.R1_NEG + fwdAltAlignments.R2_POS + fwdAltAlignments.R2_NEG;
                    refBiasMetric.REV_CXT_REF_BASES = revRefAlignments.R1_POS + revRefAlignments.R1_NEG + revRefAlignments.R2_POS + revRefAlignments.R2_NEG;
                    refBiasMetric.REV_CXT_ALT_BASES = revAltAlignments.R1_POS + revAltAlignments.R1_NEG + revAltAlignments.R2_POS + revAltAlignments.R2_NEG;

                    // calculate derived stats (conveniently in a separate method)
                    detailMetric.calculateDerivedStatistics();
                    refBiasMetric.calculateDerivedStatistics();

                    // add to list
                    contextLevelMetrics.add(detailMetric);
                    referenceBiasMetrics.add(refBiasMetric);
                }
            }
        }

        // TODO compute higher-level metrics

        return new LibraryLevelMetrics(contextLevelMetrics, artifactLevelMetrics, referenceBiasMetrics, summaryMetrics);
    }

    /**
     * TODO remove this
     */
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

                        clm.PRO_REF_BASES = this.alleleCountsPerContext.get(context).get(refBase);
                        clm.PRO_ALT_BASES = this.alleleCountsPerContext.get(context).get(altBase);
                        clm.CON_REF_BASES = this.alleleCountsPerContext.get(SequenceUtil.reverseComplement(context)).get(SequenceUtil.complement(refBase));
                        clm.CON_ALT_BASES = this.alleleCountsPerContext.get(SequenceUtil.reverseComplement(context)).get(SequenceUtil.complement(altBase));

                        clm.calculateDerivedStatistics();
                        contextLevelMetrics.add(clm);
                    }
                }
            }

            /**
             * 2. collapse context-level metrics into artifact-level metrics.
             */
            final CollectionUtil.TwoKeyHashMap<Byte, Byte, FfpeDetailMetrics> artifactLevelMetricsMap =
                    new CollectionUtil.TwoKeyHashMap<Byte, Byte, FfpeDetailMetrics>();
            for (final byte refBase : BASES) {
                for (final byte altBase : BASES) {
                    if (altBase != refBase) {
                        final FfpeDetailMetrics alm = new FfpeDetailMetrics();
                        alm.SAMPLE_ALIAS = sampleAlias;
                        alm.LIBRARY = library;
                        alm.PRO_REF_BASES = 0;
                        alm.PRO_ALT_BASES = 0;
                        alm.CON_REF_BASES = 0;
                        alm.CON_ALT_BASES = 0;
                        artifactLevelMetricsMap.put(refBase, altBase, alm);
                    }
                }
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
            //return new LibraryLevelMetrics(contextLevelMetrics, artifactLevelMetrics, summaryMetrics);
            return null;
        }
    }

    /**
     *
     */
    private class FfpeCalculator {
        private final Set<String> acceptedContexts;
        private final Set<String> acceptedLibraries;
        private final Map<String, ContextAccumulator> libraryMap;

        private FfpeCalculator(final Set<String> libraries, final Set<String> contexts) {
            this.acceptedLibraries = libraries;
            this.acceptedContexts = contexts;
            this.libraryMap = new HashMap<String, ContextAccumulator>();
            for (final String library : libraries) {
                this.libraryMap.put(library, new ContextAccumulator(contexts));
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

                // Count the base
                final byte readBase = rec.getReadBase();
                this.libraryMap.get(library).addRecord(samrec, readBase, refContext);
            }
        }

        private List<LibraryLevelMetrics> finalize(final String sampleAlias) {
            final List<LibraryLevelMetrics> allMetrics = new ArrayList<LibraryLevelMetrics>();
            for (final String library : this.libraryMap.keySet()) {
                final ContextAccumulator accumulator = this.libraryMap.get(library);
                final LibraryLevelMetrics metrics = extractMetrics(sampleAlias, library, accumulator);
                allMetrics.add(metrics);
            }
            return allMetrics;
        }
    }
}
