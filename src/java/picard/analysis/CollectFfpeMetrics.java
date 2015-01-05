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

    @Option(doc = "The number of context bases to include on each side of the assayed base.")
    public int CONTEXT_SIZE = 1;

    @Option(doc = "The optional set of sequence contexts to restrict analysis to. If not supplied all contexts are analyzed.")
    public Set<String> CONTEXTS = new HashSet<String>();

    @Option(doc = "For debugging purposes: stop after visiting this many sites with at least 1X coverage.")
    public int STOP_AFTER = Integer.MAX_VALUE;

    private final Log log = Log.getInstance(CollectFfpeMetrics.class);
    private static final String UNKNOWN_LIBRARY = "UnknownLibrary";
    private static final String UNKNOWN_SAMPLE = "UnknownSample";

    /**
     * Metrics class for outputs.
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

        // TODO

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
            final String context = StringUtil.bytesToString(bases, index - CONTEXT_SIZE, 1 + (2 * CONTEXT_SIZE)).toUpperCase();

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
     * A little class for counting alleles.
     */
    private static class Counts {

        // TODO

        int total() {
            return 0;
        }
    }

    /**
     * Class that calculated CpCG metrics for a specific library.
     */
    private class Calculator {
        private final String library;
        private final String context;

        // Things to be accumulated
        int sites = 0;
        // TODO

        Calculator(final String library, final String context) {
            this.library = library;
            this.context = context;
        }

        void accept(final SamLocusIterator.LocusInfo info, final byte refBase) {
            final Counts counts = computeAlleleFraction(info, refBase);

            if (counts.total() > 0) {
                // Things calculated on all sites with coverage
                this.sites++;

                // TODO
            }
        }

        CpcgMetrics finish() {
            final CpcgMetrics m = new CpcgMetrics();
            m.LIBRARY = this.library;
            m.CONTEXT = this.context;
            m.TOTAL_SITES = this.sites;
            // TODO
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

                // TODO magic
            }

            return counts;
        }
    }
}
