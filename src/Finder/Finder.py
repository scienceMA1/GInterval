from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from specificsite.interval import Interval
import pysam


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


# The base-class of finders which extract information from alignment results.
class Finder(object):
    def __init__(self, file=None):
        self.fw = None
        if file is not None:
            self.fw = open(file, "w+")

    def __del__(self):
        if self.fw is not None:
            self.fw.close()
            self.fw = None

    def close(self):
        if self.fw is not None:
            self.fw.close()
            self.fw = None

    # Processing the AlignedSegment object and
    def process(self, **kwargs):
        assert False

    def write_line(self, s):
        if self.fw is not None:
            self.fw.write(s + "\n")


# Extracting the five-end site intervals from alignment results.
class FiveFinder(Finder):
    def process(self, chrom, read, **kwargs):
        if not read.is_reverse:
            position = read.reference_start
            forward = True
        else:
            position = read.reference_end - 1
            forward = False
        interval = Interval(chrom=chrom, chrom_start=position, chrom_end=position, name=read.query_name + "_FE",
                            forward=forward)
        self.write_line(interval.to_bed_format_string(6))
        return interval


# Extracting the three-end site intervals from alignment results.
class ThreeFinder(Finder):
    def process(self, chrom, read, **kwargs):
        if not read.is_reverse:
            position = read.reference_end - 1
            forward = True
        else:
            position = read.reference_start
            forward = False
        interval = Interval(chrom=chrom, chrom_start=position, chrom_end=position, name=read.query_name + "_TE",
                            forward=forward)
        self.write_line(interval.to_bed_format_string(6))
        return interval


# Extracting the deletion site intervals from alignment results.
class DeletionFinder(Finder):
    def process(self, chrom, read, fa, **kwargs):
        start = read.reference_start
        intervals = list()
        for operations, number in read.cigartuples:
            if operations == pysam.CMATCH or operations == pysam.CREF_SKIP:
                start += number
            if operations == pysam.CDEL:
                end = start + number - 1
                interval = Interval(chrom=chrom, chrom_start=start, chrom_end=end, name=read.query_name + "_D",
                                    forward=not read.is_reverse)
                start += number
                seq = fa.get_sequence(interval)     # The deleted bases.
                interval.name = interval.name + "_" + str(seq)
                self.write_line(interval.to_bed_format_string(6))
                intervals.append(interval)
        return intervals


# Extracting the insertion site intervals from alignment results.
class InsertionFinder(Finder):
    def process(self, chrom, read, **kwargs):
        start = read.reference_start
        index = 0
        for operations, number in read.cigartuples:
            if operations == pysam.CMATCH or operations == pysam.CDEL or operations == pysam.CREF_SKIP:
                start += number
            if operations == pysam.CMATCH or operations == pysam.CINS or operations == pysam.CSOFT_CLIP:
                index += number
            if operations == pysam.CINS:
                if not read.is_reverse:
                    position = start
                    forward = True
                else:
                    position = start - 1
                    forward = False
                name = read.query_name + "_I_" + read.query_sequence[index-number:index]
                interval = Interval(chrom=chrom, chrom_start=position, chrom_end=position, name=name, forward=forward)
                self.write_line(interval.to_bed_format_string(6))


class FiveClippingFinder(Finder):
    def process(self, chrom, read, fasta, **kwargs):
        interval = None
        if len(read.cigartuples) >= 3:
            qseq = None
            if not read.is_reverse and read.cigartuples[0][0] == pysam.CSOFT_CLIP:
                number = read.cigartuples[0][1]
                start = read.reference_start - number
                end = read.reference_start - 1
                qseq = read.query_sequence[:number]
                interval = Interval(chrom=chrom, chrom_start=start, chrom_end=end, name=read.query_name, forward=True)
            elif read.is_reverse and read.cigartuples[-1][0] == pysam.CSOFT_CLIP:
                number = read.cigartuples[-1][1]
                start = read.reference_end
                end = read.reference_end + number - 1
                qseq = reverse_complement(read.query_sequence[-number:])
                interval = Interval(chrom=chrom, chrom_start=start, chrom_end=end, name=read.query_name, forward=False)
            if interval is not None:
                rseq = fasta.get_sequence(interval)
                # [read name]_[query sequence]_[reference sequence]
                interval.name += "_" + qseq + "_" + rseq
                self.write_line(interval.to_bed_format_string(6))
        return interval


class ThreeClippingFinder(Finder):
    def process(self, chrom, read, fasta, **kwargs):
        interval = None
        if len(read.cigartuples) >= 3:
            qseq = None
            if not read.is_reverse and read.cigartuples[-1][0] == pysam.CSOFT_CLIP:
                number = read.cigartuples[-1][1]
                start = read.reference_end
                end = read.reference_end + number - 1
                qseq = read.query_sequence[-number:]
                interval = Interval(chrom=chrom, chrom_start=start, chrom_end=end, name=read.query_name, forward=True)
            elif read.is_reverse and read.cigartuples[0][0] == pysam.CSOFT_CLIP:
                number = read.cigartuples[0][1]
                start = read.reference_start - number
                end = read.reference_start - 1
                qseq = reverse_complement(read.query_sequence[:number])
                interval = Interval(chrom=chrom, chrom_start=start, chrom_end=end, name=read.query_name, forward=False)
            if interval is not None:
                rseq = fasta.get_sequence(interval)
                # [read name]_[query sequence]_[reference sequence]
                interval.name += "_" + qseq + "_" + rseq
                self.write_line(interval.to_bed_format_string(6))
        return interval


# Extracting the query sequences from the alignment results.
class QueryFinder(Finder):
    def process(self, chrom, read, **kwargs):
        name = ">%s %s:%d-%d(%s)" % (
            read.query_name, chrom,
            read.reference_start,
            read.reference_end,
            "-" if read.is_reverse else "+"
        )

        seq = ""
        start = 0
        for operations, number in read.cigartuples:
            if operations == 0:
                seq += read.query_alignment_sequence[start:start + number]
                start += number
            if operations == 1:
                start += number
        if read.is_reverse:
            seq = str(Seq(seq).reverse_complement())

        self.write_line(name)
        self.write_line(seq)
        return SeqRecord(seq=seq, name=name)


class MismatchFinder(Finder):
    MISMATCH_TYPES = [
        "A-C", "A-T", "A-G", "A-N",
        "C-A", "C-T", "C-G", "C-N",
        "T-A", "T-C", "T-G", "T-N",
        "G-A", "G-C", "G-T", "G-N",
        "N-A", "N-C", "N-T", "N-G"
    ]

    def __init__(self, file=None):
        super(MismatchFinder, self).__init__(file)
        self.summary = {mt: 0 for mt in self.MISMATCH_TYPES}
        self.mismatch_number = dict()

    def process(self, rseq=None, qseq=None, interval=None, **kwargs):
        intervals = []
        if len(rseq) != len(qseq):
            raise Exception("The length of reference sequence and query sequence is not identical.")

        index = 0
        rseq = rseq.upper()
        for a, b in zip(rseq, qseq):
            if a != b:
                position = interval.get_position_by_index(index)
                if position is None:
                    position = interval.get_position_by_index(index)
                mt = "%s-%s" % (a, b)
                interval0 = Interval(chrom=interval.chrom, chrom_start=position, chrom_end=position,
                                     name="%s_%s" % (interval.name, mt), forward=interval.strand == "+")
                intervals.append(interval0)
                self.summary[mt] += 1
                self.write_line(interval0.to_bed_format_string(6))
            index += 1

        number = len(intervals)
        self.mismatch_number.setdefault(number, 0)
        self.mismatch_number[number] += 1

        return intervals


if __name__ == '__main__':
    bam = pysam.AlignmentFile("")
    for record in bam.fetch(until_eof=True):
        #record = pysam.AlignedSegment()
        if record.is_unmapped:
            continue
        number = None
        if record.is_reverse:
            if record.cigartuples[-1][0] == pysam.CSOFT_CLIP:
                number = record.cigartuples[-1][0]
        else:
            if record.cigartuples[0][0] == pysam.CSOFT_CLIP:
                number = record.cigartuples[0][0]
        if number is not None:
            print(number)

    fcf = FiveClippingFinder("../five.clipping.bed")
    exit(0)

