from GInterval import GComplexInterval, Interval
import os


class Parser(object):
    def __init__(self, path):
        self.path = path


class ParserFactory(object):
    format_map = {
        '.BED': "BED",
        '.FQ': 'FASTQ',
        '.FASTQ': 'FASTQ',
        '.FA': 'FASTA',
        '.FASTA': 'FASTA',
        '.GTF': 'GTF'
    }

    def __init__(self, path, fmt=None):
        self.path = path
        self.fmt = fmt

    @property
    def parser(self):
        if self.fmt is None:
            self.fmt = ParserFactory.infer_format_by_suffix(self.path)
        if self.fmt is None:
            raise ValueError('Format is None and we can not infer the format!')
        fmt = self.fmt
        parser = None
        if fmt == 'BED':
            parser = BEDParser(self.path)
        if parser is None:
            raise ValueError('We can not find proper Parser for %s format.' % fmt)
        return parser

    @classmethod
    def infer_format_by_suffix(cls, path):
        fmt = None
        x = os.path.splitext(path)[1]
        if x != '':
            fmt = cls.format_map.get(x.upper())
        return fmt

    @classmethod
    def infer_format_by_context(cls, path):
        return None


class BEDParser(Parser):
    def __init__(self, path):
        super(BEDParser, self).__init__(path)

    def __iter__(self):
        with open(self.path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cols = line.strip("\n").split("\t")
                cnum = len(cols)
                if cnum < 4:
                    raise Exception("Record least than 4 columns")
                parameters = dict()
                parameters['chrom'] = cols[0]
                start = int(cols[1])
                end = int(cols[2])
                parameters['x'] = start
                parameters['y'] = end
                parameters['name'] = cols[3]
                if cnum > 4:
                    if cols[4] != '.':
                        parameters['score'] = float(cols[4])
                    else:
                        parameters['score'] = None
                    parameters['strand'] = cols[5]
                if cnum > 6:
                    thick_x = int(cols[6])
                    thick_y = int(cols[7])
                    parameters['thick'] = None
                    try:
                        parameters['thick'] = Interval(thick_x, thick_y)
                    except AssertionError as err:
                        pass
                    parameters['rgb'] = cols[8]
                    sizes = map(int, filter(lambda obj: obj != '', cols[10].split(',')))
                    starts = map(int, filter(lambda obj: obj != '', cols[11].split(',')))
                    intervals = []
                    for bstart, bsize in zip(starts, sizes):
                        x = start + bstart
                        y = x + bsize
                        intervals.append(Interval(x, y))
                    parameters['blocks'] = intervals

                yield GComplexInterval(**parameters)


class GTFParser(Parser):
    def __init__(self, path):
        super(GTFParser, self).__init__(path)


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
                    assert (interval.chrom == self.__interval.chrom and interval.x >= self.__interval.x) or \
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
    for interval in ParserFactory("../data/P-Body-Enrichment_Class1.sorted.bed").parser:
        print(interval.to_bed_format_string())
