from ginterval.parser.factory import ParserFactory


class Writer(object):
    pass


class ShiftLoader(object):
    """

    The reference and shift bed files must bed sorted by bedtools sort.

    Or use command: sort -k1,1 -k2,2n

    The sorted results of igvtools is not supported.

    """

    def __init__(self, path, fmt=None):

        self.__generator = ParserFactory(path=path, fmt=fmt).parser.__iter__()

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
                # Check the record has been sorted.
                try:
                    assert (interval.chrom == self.__interval.chrom and interval._x >= self.__interval.x) or \
                           (interval.chrom > self.__interval.chrom)
                except AssertionError as e1:
                    print(self.__interval)
                    print(interval)
                    raise Exception("The shift intervals must be sorted!")
            self.__interval = interval
        except StopIteration as e2:
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
            valid_start = reference.x - extend
            valid_end = reference.y + extend
            intervals = []
            removable = True
            invalid_number = 0

            # Reset the index.
            self.__buffer_index = 0

            while True:
                interval = self.__next_interval()
                if interval is None:
                    break
                if interval.y < valid_start and removable:
                    invalid_number += 1
                else:
                    removable = False
                    if interval.x > valid_end:
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
    '''fw = open('../data/combine.bed', 'w+')
    parser1 = ParserFactory("../data/no_cut8_1M_R1_uniq.100000.bed", fmt='BED').parser
    parser2 = ParserFactory("../data/no_cut8_1M_R2_uniq.100000.bed", fmt='BED').parser
    for g1, g2 in zip(parser1, parser2):
        if g1.name == 'E00509:219:HGGJGCCXY:2:1216:7669:29261/1':
            print(g1.to_bed_format_string())
            print(g2.to_bed_format_string())
            print((g1 + g2).to_bed_format_string())

        if g1.chrom != g2.chrom:
            continue

        if g1.forward == g2.forward:
            continue

        fw.write((g1 + g2).to_bed_format_string() + "\n")
    fw.close()'''

    parser1 = ParserFactory("../../test/data/test.bed", fmt='BED').parser
    for gi in parser1:
        print(gi.to_bed_format_string())

    pass


    # for interval in ParserFactory("../data/genecode.v27lift37", fmt='GTF').parser:
    #    print(interval.to_bed_format_string())
    #    break
