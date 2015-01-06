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

    @Option(doc = "The optional set of sequence contexts to restrict analysis to. If not supplied all contexts are analyzed.")
    public Set<String> CONTEXTS = new HashSet<String>();

    @Option(doc = "For debugging purposes: stop after visiting this many sites with at least 1X coverage.")
    public int STOP_AFTER = Integer.MAX_VALUE;

    private final Log log = Log.getInstance(CollectFfpeMetrics.class);
    private static final String UNKNOWN_LIBRARY = "UnknownLibrary";
    private static final String UNKNOWN_SAMPLE = "UnknownSample";

    public static class FfpeSummaryMetrics extends MetricBase {
        public String SAMPLE_ALIAS;
        public String LIBRARY;

        public double A_TO_C_ERROR_RATE;
        public double A_TO_C_QSCORE;
        public double A_TO_G_ERROR_RATE;
        public double A_TO_G_QSCORE;
        public double A_TO_T_ERROR_RATE;
        public double A_TO_T_QSCORE;

        public double C_TO_A_ERROR_RATE;
        public double C_TO_A_QSCORE;
        public double C_TO_G_ERROR_RATE;
        public double C_TO_G_QSCORE;
        public double C_TO_T_ERROR_RATE;
        public double C_TO_T_QSCORE;

        public double G_TO_A_ERROR_RATE;
        public double G_TO_A_QSCORE;
        public double G_TO_C_ERROR_RATE;
        public double G_TO_C_QSCORE;
        public double G_TO_T_ERROR_RATE;
        public double G_TO_T_QSCORE;

        public double T_TO_A_ERROR_RATE;
        public double T_TO_A_QSCORE;
        public double T_TO_C_ERROR_RATE;
        public double T_TO_C_QSCORE;
        public double T_TO_G_ERROR_RATE;
        public double T_TO_G_QSCORE;
    }

    public static class FfpeDetailMetrics extends MetricBase {
        public String SAMPLE_ALIAS;
        public String LIBRARY;
        public String CONTEXT;

        public long TOTAL_SITES;
        public long TOTAL_BASES;

        public long ALT_A_BASES;
        public long ALT_C_BASES;
        public long ALT_G_BASES;
        public long ALT_T_BASES;

        public double A_ERROR_RATE;
        public double C_ERROR_RATE;
        public double G_ERROR_RATE;
        public double T_ERROR_RATE;

        public double A_QSCORE;
        public double C_QSCORE;
        public double G_QSCORE;
        public double T_QSCORE;

        // TODO p-values?
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

        return messages.isEmpty() ? null : messages.toArray(new String[messages.size()]);
    }

    /** Mimic of Oracle's nvl() - returns the first value if not null, otherwise the second value. */
    private <T> T nvl(final T value1, final T value2) {
        if (value1 != null) return value1;
        else return value2;
    }

    @Override
    protected int doWork() {
        final File SUMMARY_OUT = new File(OUTPUT + "ffpe_summary_metrics");
        final File DETAILS_OUT = new File(OUTPUT + "ffpe_detail_metrics");

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

        // Setup the calculators
        final Set<String> contexts = CONTEXTS.isEmpty() ? makeContextStrings(CONTEXT_SIZE) : CONTEXTS;
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

        for (final SamLocusIterator.LocusInfo info : iterator) {
            // Skip dbSNP sites
            final String chrom = info.getSequenceName();
            final int pos = info.getPosition();
            final int index = pos - 1;
            if (dbSnp != null && dbSnp.isDbSnpSite(chrom, pos)) continue;

            // Skip sites at the end of chromosomes
            final byte[] bases = refWalker.get(info.getSequenceIndex()).getBases();
            if (pos < 3 || pos > bases.length - 3) continue;

            // Get the reference context string and perform counting
            final String context = StringUtil.bytesToString(bases, index - CONTEXT_SIZE, 1 + (2 * CONTEXT_SIZE)).toUpperCase();
            counts.countAlleles(info, context);

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

        final MetricsFile<FfpeSummaryMetrics, Integer> summaryMetricsFile = getMetricsFile();
        final MetricsFile<FfpeDetailMetrics, Integer> detailMetricsFile = getMetricsFile();
        // TODO add metrics

        summaryMetricsFile.write(SUMMARY_OUT);
        detailMetricsFile.write(DETAILS_OUT);
        CloserUtil.close(in);
        return 0;
    }

    private Set<String> makeContextStrings(final int contextSize) {
        final Set<String> contexts = new HashSet<String>();

        for (final byte[] kmer : generateAllKmers(2 * contextSize + 1)) {
           contexts.add(StringUtil.bytesToString(kmer));
        }

        log.info("Generated " + contexts.size() + " context strings.");
        return contexts;
    }

    /** Generates all possible unambiguous kmers of length and returns them as byte[]s. */
    private List<byte[]> generateAllKmers(final int length) {
        final List<byte[]> sofar = new LinkedList<byte[]>();
        final byte[] bases = {'A', 'C', 'G', 'T'};

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
                for (final byte b : bases) {
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

    private static class AlleleCounter {
        private final byte refBase;
        private long A = 0;
        private long C = 0;
        private long G = 0;
        private long T = 0;

        private AlleleCounter(final byte refBase) {
            this.refBase = refBase;
        }

        private void add(byte base) {
            switch (base) {
                case 'A': this.A++; break;
                case 'C': this.C++; break;
                case 'G': this.G++; break;
                case 'T': this.T++; break;
                default: throw new IllegalStateException(base + " is not a valid DNA base");
            }
        }

        private long total() { return A + C + G + T; }

        private long altA() { return (refBase == 'A') ? -1 : A; }
        private long altC() { return (refBase == 'C') ? -1 : C; }
        private long altG() { return (refBase == 'G') ? -1 : G; }
        private long altT() { return (refBase == 'T') ? -1 : T; }
    }

    /**
     * A library-level accumulator for state transitions.
     */
    private class FfpeCalculator {
        private final Set<String> acceptedContexts;
        private final Set<String> acceptedLibraries;
        private final Map<String, Map<String, AlleleCounter>> countsPerLibrary;

        private FfpeCalculator(final Set<String> libraries, final Set<String> contexts) {
            this.acceptedLibraries = libraries;
            this.acceptedContexts = contexts;

            // TODO make this less ugly?
            this.countsPerLibrary = new HashMap<String, Map<String, AlleleCounter>>();
            for (final String library : libraries) {
                final Map<String, AlleleCounter> countsPerContext = new HashMap<String, AlleleCounter>();
                for (final String context : contexts) {
                    final byte refBase = (byte) context.charAt(CONTEXT_SIZE);
                    countsPerContext.put(context, new AlleleCounter(refBase));
                }
                this.countsPerLibrary.put(library, countsPerContext);
            }
        }

        /**
         * This is where we actually do the counting - examine all reads overlapping the given reference locus, and
         * count up the observed genotypes.
         */
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

                final String originalTemplateContext;
                final byte calledBase;
                /**
                 * Remember that if the read is aligned to the negative strand, the reported ref + call bases must be reverse complemented
                 * to reflect the original sequence. Likewise if the read is second in a pair (and thus derived from the reverse template),
                 * the ref + call bases must be revc'd to reflect the forward template. If both of these conditions are true, they cancel
                 * each other out (hence the XOR).
                 */
                final boolean isNegativeStrand = samrec.getReadNegativeStrandFlag();
                final boolean isReadTwo = samrec.getReadPairedFlag() && samrec.getSecondOfPairFlag();
                if (isNegativeStrand ^ isReadTwo) {
                    originalTemplateContext = SequenceUtil.reverseComplement(refContext);
                    calledBase = SequenceUtil.complement(rec.getReadBase());
                } else {
                    originalTemplateContext = refContext;
                    calledBase = rec.getReadBase();
                }

                final String library = nvl(samrec.getReadGroup().getLibrary(), UNKNOWN_LIBRARY);

                // Skip if context or library is unknown
                if (!acceptedLibraries.contains(library)) continue;
                if (!acceptedContexts.contains(originalTemplateContext)) continue;

                // Count the base
                this.countsPerLibrary.get(library).get(originalTemplateContext).add(calledBase);
            }
        }

        private void writeMetrics(final String sampleAlias) {
            List<FfpeSummaryMetrics> summaryMetrics = new ArrayList<FfpeSummaryMetrics>();
            List<FfpeDetailMetrics> detailMetrics = new ArrayList<FfpeDetailMetrics>();


            /**
             * for each library:
             *     for each base:
             *         for each ref context centered on that base:
             *             (1) counts of each alt base
             *             (2) error rates / q-values for each alt base based on (1)
             *         (3) error rates / q-values for each alt base based on sum of (1) across contexts
             *
             */
            for (String library : this.acceptedLibraries) {

                for (String context : this.acceptedContexts) {
                    FfpeDetailMetrics detail = new FfpeDetailMetrics();
                    detail.SAMPLE_ALIAS = sampleAlias;
                    detail.LIBRARY = library;
                    detail.CONTEXT = context;

                    // TODO
                }

                FfpeSummaryMetrics summary = new FfpeSummaryMetrics();
                summary.SAMPLE_ALIAS = sampleAlias;
                summary.LIBRARY = library;
                // TODO
            }
        }
    }
}
