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
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
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
            doc = "Location of output metrics file to write.")
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

    @Option(doc = "The number of context bases to include on each side of the assayed G/C base.")
    public int CONTEXT_SIZE = 1;

    @Option(doc = "The optional set of sequence contexts to restrict analysis to. If not supplied all contexts are analyzed.")
    public Set<String> CONTEXTS = new HashSet<String>();

    @Option(doc = "For debugging purposes: stop after visiting this many sites with at least 1X coverage.")
    public int STOP_AFTER = Integer.MAX_VALUE;

    private final Log log = Log.getInstance(CollectFfpeMetrics.class);
    private static final String UNKNOWN_LIBRARY = "UnknownLibrary";
    private static final String UNKNOWN_SAMPLE = "UnknownSample";

    /**
     * Metrics class for outputs. TODO expand
     */
    public static final class CpcgMetrics extends MetricBase {
        /** The name of the sample being assayed. */
        public String SAMPLE_ALIAS;
        /** The name of the library being assayed. */
        public String LIBRARY;
        /** The sequence context being reported on. */
        public String CONTEXT;
        /** The total number of sites that had at least one base covering them. */
        public int TOTAL_SITES;
        /** The total number of basecalls observed at all sites. */
        public long TOTAL_BASES;

        /** The number of reference alleles observed as C in read 1 and G in read 2. */
        public long REF_NONOXO_BASES;
        /** The number of reference alleles observed as G in read 1 and C in read 2. */
        public long REF_OXO_BASES;
        /** The total number of reference alleles observed */
        public long REF_TOTAL_BASES;

        /**
         * The count of observed A basecalls at C reference positions and T basecalls
         * at G reference bases that are correlated to instrument read number in a way
         * that rules out oxidation as the cause
         */
        public long ALT_NONOXO_BASES;

        /**
         * The count of observed A basecalls at C reference positions and T basecalls
         * at G reference bases that are correlated to instrument read number in a way
         * that is consistent with oxidative damage.
         */
        public long ALT_OXO_BASES;

        /** The oxo error rate, calculated as max(ALT_OXO_BASES - ALT_NONOXO_BASES, 1) / TOTAL_BASES */
        public double OXIDATION_ERROR_RATE;

        /** -10 * log10(OXIDATION_ERROR_RATE) */
        public double OXIDATION_Q;

        // Fields below this point are metrics that are calculated to see if there is oxidative damage that is
        // biased toward the reference base - i.e. occurs more where the reference base is a C vs. a G or vice
        // versa.

        /** The number of ref basecalls observed at sites where the genome reference == C. */
        public long C_REF_REF_BASES;
        /** The number of ref basecalls observed at sites where the genome reference == G. */
        public long G_REF_REF_BASES;
        /** The number of alt (A/T) basecalls observed at sites where the genome reference == C. */
        public long C_REF_ALT_BASES;
        /** The number of alt (A/T) basecalls observed at sites where the genome reference == G. */
        public long G_REF_ALT_BASES;

        /**
         * The rate at which C>A and G>T substitutions are observed at C reference sites above the expected rate if there
         * were no bias between sites with a C reference base vs. a G reference base.
         */
        public double C_REF_OXO_ERROR_RATE;
        /** C_REF_OXO_ERROR_RATE expressed as a phred-scaled quality score. */
        public double C_REF_OXO_Q;
        /**
         * The rate at which C>A and G>T substitutions are observed at G reference sites above the expected rate if there
         * were no bias between sites with a C reference base vs. a G reference base.
         */
        public double G_REF_OXO_ERROR_RATE;
        /** G_REF_OXO_ERROR_RATE expressed as a phred-scaled quality score. */
        public double G_REF_OXO_Q;




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
    private final <T> T nvl(final T value1, final T value2) {
        if (value1 != null) return value1;
        else return value2;
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
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
        final ListMap<String, Calculator> calculators = new ListMap<String, Calculator>();
        for (final String context : contexts) {
            for (final String library : libraries) {
                calculators.add(context, new Calculator(library, context));
            }
        }

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

            // Get the reference base at this locus
            final byte base = StringUtil.toUpperCase(bases[index]);

            // Get the reference context string
            // TODO expand
            final String context;
            {
                // the point of the reverse complement here is merely because makeContextStrings only generated C-contexts.
                // however that's no longer the case, so we don't need it?
                final String tmp = StringUtil.bytesToString(bases, index - CONTEXT_SIZE, 1 + (2 * CONTEXT_SIZE)).toUpperCase();
                if (base == 'C') context = tmp;
                else /* if (base == 'G') */ context = SequenceUtil.reverseComplement(tmp);
            }

            final List<Calculator> calculatorsForContext = calculators.get(context);
            if (calculatorsForContext == null) continue; // happens if we get ambiguous (e.g. N) bases in the reference
            for (final Calculator calc : calculatorsForContext) calc.accept(info, base);

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

        final MetricsFile<CpcgMetrics, Integer> file = getMetricsFile();
        for (final List<Calculator> calcs : calculators.values()) {
            for (final Calculator calc : calcs) {
                final CpcgMetrics m = calc.finish();
                m.SAMPLE_ALIAS = StringUtil.join(",", new ArrayList<String>(samples));
                file.addMetric(m);
            }
        }

        file.write(OUTPUT);
        CloserUtil.close(in);
        return 0;
    }

    private Set<String> makeContextStrings(final int contextSize) {
        final Set<String> contexts = new HashSet<String>();

        // TODO should we limit this to say, C and T contexts?
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

    /**
     * A little class for counting alleles. TODO expand
     */
    private static class Counts {
        int controlA;
        int oxidatedA;
        int controlC;
        int oxidatedC;
        int controlG;
        int oxidatedG;
        int controlT;
        int oxidatedT;

        int total() { return controlC + oxidatedC + controlA + oxidatedA + controlG + oxidatedG + controlT + oxidatedT; }
    }




    /**
     * Class that calculated CpCG metrics for a specific library. TODO expand
     */
    private class Calculator {
        private final String library;
        private final String context;

        // Things to be accumulated
        int sites = 0;
        long refCcontrolA = 0;
        long refCoxidatedA = 0;
        long refCcontrolC = 0;
        long refCoxidatedC = 0;
        long refGcontrolA = 0;
        long refGoxidatedA = 0;
        long refGcontrolC = 0;
        long refGoxidatedC = 0;

        Calculator(final String library, final String context) {
            this.library = library;
            this.context = context;
        }

        void accept(final SamLocusIterator.LocusInfo info, final byte refBase) {
            final Counts counts = computeAlleleFraction(info, refBase);

            if (counts.total() > 0) {
                // Things calculated on all sites with coverage
                this.sites++;
                if (refBase == 'C') {
                    this.refCcontrolA += counts.controlA;
                    this.refCoxidatedA += counts.oxidatedA;
                    this.refCcontrolC += counts.controlC;
                    this.refCoxidatedC += counts.oxidatedC;
                } else if (refBase == 'G') {
                    this.refGcontrolA += counts.controlA;
                    this.refGoxidatedA += counts.oxidatedA;
                    this.refGcontrolC += counts.controlC;
                    this.refGoxidatedC += counts.oxidatedC;
                } else {
                    throw new IllegalStateException("Reference bases other than G and C not supported.");
                }
            }
        }

        CpcgMetrics finish() {
            final CpcgMetrics m = new CpcgMetrics();
            m.LIBRARY = this.library;
            m.CONTEXT = this.context;
            m.TOTAL_SITES = this.sites;
            m.TOTAL_BASES = this.refCcontrolC + this.refCoxidatedC + this.refCcontrolA + this.refCoxidatedA +
                    this.refGcontrolC + this.refGoxidatedC + this.refGcontrolA + this.refGoxidatedA;
            m.REF_OXO_BASES = this.refCoxidatedC + refGoxidatedC;
            m.REF_NONOXO_BASES = this.refCcontrolC + this.refGcontrolC;
            m.REF_TOTAL_BASES = m.REF_OXO_BASES + m.REF_NONOXO_BASES;
            m.ALT_NONOXO_BASES = this.refCcontrolA + this.refGcontrolA;
            m.ALT_OXO_BASES = this.refCoxidatedA + this.refGoxidatedA;

            /**
             * Why do we calculate the oxo error rate using oxidatedA - controlA you ask?  We know that all the
             * bases counted in oxidatedA are consistent with 8-oxo-G damage during shearing, but not all of them
             * will have been caused by this. If we assume that C>A errors caused by other factors will occur randomly
             * with respect to read1/read2, then we should see as many in the 8-oxo-G consistent state as not.  So we
             * assume that controlA is half the story, and remove the other half from oxidatedA.
             */
            m.OXIDATION_ERROR_RATE = Math.max(m.ALT_OXO_BASES - m.ALT_NONOXO_BASES, 1) / (double) m.TOTAL_BASES;
            m.OXIDATION_Q = -10 * log10(m.OXIDATION_ERROR_RATE);

            /** Now look for things that have a reference base bias! */
            m.C_REF_REF_BASES = this.refCcontrolC + this.refCoxidatedC;
            m.G_REF_REF_BASES = this.refGcontrolC + this.refGoxidatedC;
            m.C_REF_ALT_BASES = this.refCcontrolA + this.refCoxidatedA;
            m.G_REF_ALT_BASES = this.refGcontrolA + this.refGoxidatedA;

            final double cRefErrorRate = m.C_REF_ALT_BASES / (double) (m.C_REF_ALT_BASES + m.C_REF_REF_BASES);
            final double gRefErrorRate = m.G_REF_ALT_BASES / (double) (m.G_REF_ALT_BASES + m.G_REF_REF_BASES);

            m.C_REF_OXO_ERROR_RATE = Math.max(cRefErrorRate - gRefErrorRate, 1e-10);
            m.G_REF_OXO_ERROR_RATE = Math.max(gRefErrorRate - cRefErrorRate, 1e-10);
            m.C_REF_OXO_Q = -10 * log10(m.C_REF_OXO_ERROR_RATE);
            m.G_REF_OXO_Q = -10 * log10(m.G_REF_OXO_ERROR_RATE);

            return m;
        }

        /**
         * This is where we actually do the counting - examine all reads overlapping the given reference locus, and
         * count up the observed genotypes.
         */
        private Counts computeAlleleFraction(final SamLocusIterator.LocusInfo info, final byte refBase) {
            final Counts counts = new Counts();

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

                // Skip if below qual, or if library isn't a match
                if (qual < MINIMUM_QUALITY_SCORE) continue;
                if (!this.library.equals(nvl(samrec.getReadGroup().getLibrary(), UNKNOWN_LIBRARY))) continue;

                // Get the read base, and get it in "as read" orientation
                final byte base = rec.getReadBase();
                final byte baseAsRead = samrec.getReadNegativeStrandFlag() ? SequenceUtil.complement(base) : base;
                final int read = samrec.getReadPairedFlag() && samrec.getSecondOfPairFlag() ? 2 : 1;

                // TODO this is where the magic happens
                // Figure out how to count the alternative allele. If the damage is caused by oxidation of G
                // during shearing (in non-rnaseq data), then we know that:
                //     G>T observation is always in read 1
                //     C>A observation is always in read 2
                // But if the substitution is from other causes the distribution of A/T across R1/R2 will be
                // random.
                if (refBase == 'C') {
                    if (base == 'C') {
                        if (baseAsRead == 'G' && read == 1) ++counts.oxidatedC;
                        else if (baseAsRead == 'G' && read == 2) ++counts.controlC;
                        else if (baseAsRead == 'C' && read == 1) ++counts.controlC;
                        else if (baseAsRead == 'C' && read == 2) ++counts.oxidatedC;
                    } else if (base == 'A') {
                        if (baseAsRead == 'T' && read == 1) ++counts.oxidatedA;
                        else if (baseAsRead == 'T' && read == 2) ++counts.controlA;
                        else if (baseAsRead == 'A' && read == 1) ++counts.controlA;
                        else if (baseAsRead == 'A' && read == 2) ++counts.oxidatedA;
                    } else if (base == 'G') {
                        if (baseAsRead == 'G' && read == 1) ++counts.oxidatedG;
                        else if (baseAsRead == 'G' && read == 2) ++counts.controlG;
                        else if (baseAsRead == 'C' && read == 1) ++counts.controlG;
                        else if (baseAsRead == 'C' && read == 2) ++counts.oxidatedG;
                    } else if (base == 'T') {
                        if (baseAsRead == 'T' && read == 1) ++counts.oxidatedT;
                        else if (baseAsRead == 'T' && read == 2) ++counts.controlT;
                        else if (baseAsRead == 'A' && read == 1) ++counts.controlT;
                        else if (baseAsRead == 'A' && read == 2) ++counts.oxidatedT;
                    }
                }
            } // TODO check if above is correct; then handle other refBases

            return counts;
        }
    }
}
