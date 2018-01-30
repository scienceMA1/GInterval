# Author: Zonggui Chen
# Date: 2017-10-11

# This script is a small module for bed records manipulation.
# Usage: import this script on your scripts.

from Bio.Seq import Seq
from pysam import AlignedSegment
from pysam import FastaFile
import os


class Interval(object):
    BED4 = 4
    BED6 = 6
    BED8 = 8
    BED12 = 12

    def __init__(self, chrom, chrom_start, chrom_end, name, score=None, forward=True, rgb="255,0,0",
                 chrom_thick_start=None, chrom_thick_end=None, chrom_blocks=None):
        """

        An Interval object present a interval of genome.

        Including all the information stored in the BED12 format files.

        The coordinate of the genome is 0-base.

        For more detail definition, please reference to the definition of BED12 format by UCSC.

        :param chrom:           The chromosome location of this interval. For example, chr1
        :param chrom_start:     The start position (included) of the interval, strand is not considered.
        :param chrom_end:       The end position (included) of the interval, strand is not considered.
        :param name:            The name of the interval.
        :param score:           The score of the interval. "." represent None, like the BED12 definition of UCSC.
        :param forward:         The True when the interval is located on the forward strand.
        :param rgb:             The RGB color of the interval for drawing.
        :param chrom_thick_start:   The start position of thick interval, usually represent CDS.
        :param chrom_thick_end:     The end position of thick interval.
        :param chrom_blocks:        The blocks store the pair of [start, end], strand is not considered.
        """

        """
        Assignment the essential attributes.
        """
        self.chrom = chrom
        self.name = name
        self.score = score
        self.rgb = rgb
        self.forward = forward

        self.__chrom_start = chrom_start
        self.__chrom_end = chrom_end
        self.__chrom_thick_start = chrom_thick_start
        self.__chrom_thick_end = chrom_thick_end

        """
        Default blocks.
        """
        if chrom_blocks is None:
            self.__chrom_blocks = [[self.__chrom_start, self.__chrom_end]]
        else:
            self.__chrom_blocks = chrom_blocks  # [[12, 14], [16, 17], [19,20]]

        self.block_count = len(self.__chrom_blocks)

    """
    The properties start with "chrom_" will return the position without considering the strand.
    """

    @property
    def chrom_start(self):
        return self.__chrom_start

    @property
    def chrom_end(self):
        return self.__chrom_end

    @property
    def chrom_thick_start(self):
        return self.__chrom_thick_start

    @property
    def chrom_thick_end(self):
        return self.__chrom_thick_end

    @property
    def chrom_blocks(self):
        return self.__chrom_blocks

    """
    The properties below will consider the strand.
    """

    @property
    def start(self):
        return self.__chrom_start if self.forward else self.__chrom_end

    @property
    def start_index(self):
        return 0

    @property
    def end(self):
        return self.__chrom_end if self.forward else self.__chrom_start

    @property
    def end_index(self):
        return self.length - 1

    @property
    def thick_start(self):
        return self.__chrom_thick_start if self.forward else self.__chrom_thick_end

    @property
    def thick_end(self):
        return self.__chrom_thick_end if self.forward else self.__chrom_thick_start

    @property
    def strand(self):
        return "+" if self.forward else "-"

    @property
    def blocks(self):
        if self.forward:
            blocks = self.__chrom_blocks
        else:
            blocks = [[end, start] for start, end in self.__chrom_blocks]
            blocks.reverse()
        return blocks

    def __len__(self):
        """
        The length of an interval is the total length of all blocks.
        :return: Total length of blocks.
        """
        return sum([abs(block[1] - block[0]) + 1 for block in self.__chrom_blocks])

    def __str__(self):
        bstr = "\n".join(["%d - %d" % (block[0], block[1]) for block in self.__chrom_blocks])
        rstr = """Chromosome: %s
Name: %s
Strand: %s
Position: %d - %d
Block Count: %d
%s""" % (self.chrom, self.name, self.strand, self.__chrom_start, self.__chrom_end, self.block_count, bstr)
        return rstr

    def to_bed_format_string(self, fmt=BED6):
        """

        This function will convert interval object to BED format string.
        The string can be writen to the BED file directory.

        :param fmt: Valid values: BED4,6,8,12.
        :return: BED format string.
        """
        rstr = ""
        if fmt >= self.BED4:
            rstr = "\t".join(map(str, [self.chrom, self.__chrom_start, self.__chrom_end + 1, self.name]))
        if fmt >= self.BED6:
            rstr = "\t".join(map(str, [rstr, self.score if self.score is not None else ".", self.strand]))
        if fmt >= self.BED8:
            if self.__chrom_thick_start is None or self.__chrom_thick_end is None:
                rstr = "\t".join(map(str, [rstr, self.__chrom_start, self.__chrom_start]))
            else:
                rstr = "\t".join(map(str, [rstr, self.__chrom_thick_start, self.__chrom_thick_end + 1]))
        if fmt >= self.BED12:
            sizes = ",".join(map(str, [block[1] - block[0] + 1 for block in self.__chrom_blocks]))
            starts = ",".join(map(str, [block[0] - self.__chrom_start for block in self.__chrom_blocks]))
            rstr = "\t".join(map(str, [rstr, self.rgb, self.block_count, sizes, starts]))
        return rstr

    """
    The transform between position and index.
    """

    def get_index_by_position(self, position):
        """

        The index of specific position on the reference interval. Strand is considered.

        :param position: The position on genome.
        :return: If the position is not hit on the interval (blocks), will return None.
        """

        index = 0

        if not self.__chrom_start <= position <= self.__chrom_end:
            return None

        for start, end in self.__chrom_blocks:
            if start <= position <= end:
                index += position - start
                break
            elif end < position:
                index += end - start + 1
            elif position < start:
                return None

        return (len(self) - 1 - index) if self.strand == "-" else index

    def get_position_by_index(self, index):
        """

        The index on the reference interval.

        :param index: The position of specific index.
        :return:
        """

        index_count = len(self)

        if not -index_count <= index < index_count:
            raise IndexError("Invalid index.")

        if index < 0:
            index = index + index_count
        if not self.forward:
            index = index_count - 1 - index

        max_index = -1
        for start, end in self.__chrom_blocks:
            min_index = max_index + 1
            max_index = min_index + end - start
            if min_index <= index <= max_index:
                position = (index - min_index) + start
                return position

    """
    Some utility functions.
    """

    def get_sub_interval(self, start, end):

        if max(start, self.__chrom_start) > min(min(end, self.__chrom_end)):
            return None

        blocks = self.chrom_blocks()
        new_blocks = []

        for bstart, bend in blocks:
            block_start = max(bstart, start)
            block_end = min(bend, end)
            if block_end >= block_start:
                new_blocks.append([bstart, bend])

        if len(new_blocks) <= 0:
            return None

        pars = dict()
        pars["chrom"] = self.chrom
        pars["chrom_start"] = new_blocks[0][0]
        pars["chrom_end"] = new_blocks[-1][1]
        pars["name"] = self.name
        pars["score"] = self.score
        pars["forward"] = self.forward
        if self.__chrom_thick_start is not None and self.__chrom_thick_end is not None:
            thick_start = max(start, self.__chrom_thick_start)
            thick_end = min(end, self.__chrom_thick_end)
            if thick_end >= thick_start:
                pars["chrom_thick_start"] = thick_start
                pars["chrom_thick_end"] = thick_end
        pars["rgb"] = self.rgb
        pars["chrom_blocks"] = new_blocks

        return Interval(**pars)

    def cal_distance(self, shift, reference=None):
        """

        Calculate the distance of two intervals inside genome or reference interval.

        :param shift: The shift interval.
        :param reference: The reference interval.
        :return: The distance.
        """

        position1 = self.center_position()
        position2 = shift.center_position()

        if reference is None:
            if self.forward:
                return position2 - position1
            else:
                return position1 - position2
        else:
            index1 = reference.get_index_by_position(position1)
            index2 = reference.get_index_by_position(position2)
            if index1 is not None and index2 is not None:
                return index2 - index1
            else:
                return None

    @property
    def center_index(self):
        return self.length // 2

    @property
    def center_position(self):
        return self.get_position_by_index(self.center_index)

    @property
    def thick_start_index(self):
        return self.get_index_by_position(self.thick_start)

    @property
    def thick_end_index(self):
        return self.get_index_by_position(self.thick_end)

    @property
    def length(self):
        return sum([block[1] - block[0] + 1 for block in self.__chrom_blocks])

    @property
    def cds_interval(self):
        return self.get_sub_interval(self.__chrom_thick_start, self.__chrom_thick_end)

    @property
    def utr5_interval(self):
        if self.forward:
            start = self.__chrom_start
            end = self.get_position_by_index(self.thick_start_index - 1)
        else:
            start = self.get_position_by_index(self.thick_start_index - 1)
            end = self.__chrom_end
        return self.get_sub_interval(start, end)

    @property
    def utr3_interval(self):
        if self.forward:
            start = self.get_position_by_index(self.thick_end_index + 1)
            end = self.__chrom_end
        else:
            start = self.__chrom_start
            end = self.get_position_by_index(self.thick_end_index + 1)
        return self.get_sub_interval(start, end)

    @property
    def has_thick_interval(self):
        if self.__chrom_thick_start is None or self.__chrom_thick_end is None:
            return False
        return self.__chrom_thick_start <= self.__chrom_thick_end

    @property
    def is_mrna(self):
        return self.has_thick_interval

    @property
    def is_protein_coding(self):
        return self.has_thick_interval

    @property
    def utr5_length(self):
        return self.thick_start_index - 0

    @property
    def cds_length(self):
        return self.thick_end_index - self.thick_start_index + 1

    @property
    def utr3_length(self):
        return self.length - self.thick_end_index - 1


class IntervalFactory(object):
    @classmethod
    def from_bed_format_string(cls, record):
        """

        Constructing the interval object from BED format string.

        :param record:
        :return:
        """
        cols = record.split("\t")
        cnum = len(cols)
        if cnum < 4:
            raise Exception("At least 4 column is necessary to construct interval object.")
        pars = dict()
        pars["chrom"] = cols[0]
        chrom_start = int(cols[1])
        pars["chrom_start"] = chrom_start
        pars["chrom_end"] = int(cols[2]) - 1
        pars["name"] = cols[3]
        if cnum >= 6:
            pars["score"] = int(cols[4]) if cols[4] != "." else None
            pars["forward"] = cols[5] == "+"
        if cnum >= 8:
            pars["chrom_thick_start"] = int(cols[6])
            pars["chrom_thick_end"] = int(cols[7]) - 1
        if cnum >= 12:
            pars["rgb"] = cols[8]
            sizes = map(int, filter(lambda x: x != "", cols[10].split(",")))
            starts = map(int, filter(lambda x: x != "", cols[11].split(",")))
            pars["chrom_blocks"] = [[start + chrom_start, start + chrom_start + size - 1]
                                    for start, size in zip(starts, sizes)]
        return Interval(**pars)

    @classmethod
    def from_bed_file(cls, infile):
        """

        Generator for constructing interval objects from BED format file.

        :param infile:
        :return:
        """
        with open(infile) as f:
            for line in f:
                interval = IntervalFactory.from_bed_format_string(line.strip("\n"))
                yield interval

    @classmethod
    def from_aligned_segment(cls, chrom, record):
        """

        Transform the pysam.AlignedSegment object to interval object.

        :param chrom:
        :param record:
        :return:
        """
        assert isinstance(record, AlignedSegment)
        pars = dict()
        pars["chrom"] = chrom
        pars["chrom_start"] = record.reference_start
        pars["chrom_end"] = record.reference_end - 1
        pars["name"] = record.query_name
        pars["score"] = record.mapq
        pars["forward"] = not record.is_reverse
        pars["chrom_thick_start"] = record.reference_start
        pars["chrom_thick_end"] = record.reference_end - 1
        pars["chrom_blocks"] = [[block[0], block[1] - 1] for block in record.get_blocks()]
        return Interval(**pars)

    @classmethod
    def from_gtf_file(cls, infile):
        assert False


class SequenceExtractor(object):
    """

    Extracting sequences from FASTA file by interval objects.

    """
    def __init__(self, path):
        assert os.path.exists(path)
        self.fa = FastaFile(path)

    def get_sequence(self, interval):
        seq = ""
        for start, end in interval.chrom_blocks:
            seq += self.fa.fetch(interval.chrom, start, end + 1)
        seq = Seq(seq)
        if not interval.forward:
            seq = seq.reverse_complement()
        return seq

    def __del__(self):
        if not self.fa.closed:
            self.fa.close()


class ShiftLoader(object):
    """

    The reference and shift bed files must bed sorted by bedtools sort.

    Or use command: sort -k1,1 -k2,2n

    The sorted results of igvtools is not supported.

    """

    def __init__(self, path):
        super(ShiftLoader, self).__init__()

        assert os.path.exists(path)

        self.__generator = IntervalFactory.from_bed_file(path)

        # The list store the loaded intervals. The list should not to be accessed directly.
        self.__buffer = []
        self.__buffer_chrom = None
        self.__buffer_index = 0
        self.__eoc = False  # The end of current chromosome.

        # It store the result of next(generator). It will be None if reach the end of the file.
        self.__interval = None

        # At the beginning, we need to load the first interval.
        self.__load_one()
        if self.__interval is not None:
            self.__buffer.append(self.__interval)
            self.__buffer_chrom = self.__interval.chrom

    def __load_one(self):
        # Load one interval by the generator, and store it in the last loaded record.
        try:
            interval = next(self.__generator)
            if self.__interval is not None:
                if interval.chrom == self.__interval.chrom:
                    if interval.chrom_start < self.__interval.chrom_start:
                        print(self.__interval)
                        print(interval)
                        raise Exception("The shift intervals must be sorted!")
                elif interval.chrom < self.__interval.chrom:
                    print(self.__interval)
                    print(interval)
                    raise Exception("The shift intervals must be sorted!")
            self.__interval = interval
        except StopIteration as e:
            self.__interval = None
        return self.__interval

    def __buffer_expand(self):
        if not self.__eoc:
            self.__load_one()
            if self.__interval is None or self.__interval.chrom != self.__buffer_chrom:
                self.__eoc = True
                return None
            elif self.__interval.chrom == self.__buffer_chrom:
                self.__buffer.append(self.__interval)
                return self.__next_interval()
        else:
            return None

    def __next_interval(self):
        # Iterated the intervals on the buffer list. The list will be extended automatically if necessary.
        if self.__buffer_index >= len(self.__buffer):  # Dynamically load next one from file.
            return self.__buffer_expand()
        else:
            interval = self.__buffer[self.__buffer_index]
            self.__buffer_index += 1
            return interval

    def __next_chrom(self):  # Jump to next chromosome.
        while True:

            if self.__interval is None:
                self.__buffer_chrom = None
                self.__buffer = []
                break
            if self.__buffer_chrom != self.__interval.chrom:
                self.__buffer_chrom = self.__interval.chrom
                self.__buffer = [self.__interval]
                self.__eoc = False
                break
            self.__load_one()

    def get_shifts(self, reference, extend=0):
        # 1.
        if reference.chrom == self.__buffer_chrom:

            # Init some variables.
            valid_start = reference.chrom_start - extend
            valid_end = reference.chrom_end + extend
            intervals = []
            removable = True
            invalid_number = 0

            # Reset the index.
            self.__buffer_index = 0

            while True:
                interval = self.__next_interval()
                if interval is None:
                    break
                if interval.chrom_end < valid_start and removable:
                    invalid_number += 1
                else:
                    removable = False
                    if interval.chrom_start > valid_end:
                        break
                    else:
                        intervals.append(interval)

            # Remove the invalid intervals.
            if invalid_number > 0:
                self.__buffer = self.__buffer[invalid_number:]

            return intervals

        # 2. The shift interval reach the end of the file or reference interval has not reach the buffer chromosome.
        if self.__buffer_chrom is None or reference.chrom < self.__buffer_chrom:
            return None

        # 3. The buffer should jump to next chromosome and rerun this function.
        if reference.chrom > self.__buffer_chrom:
            self.__next_chrom()
            return self.get_shifts(reference, extend)


if __name__ == '__main__':
    pass
