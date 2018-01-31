from GInterval import GInterval, Interval
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
                        interval = Interval(x, y)
                        intervals.append(interval)
                    parameters['blocks'] = intervals

                yield GInterval(**parameters)


class GTFParser(Parser):
    def __init__(self, path):
        super(GTFParser, self).__init__(path)


class Writer(object):
    pass

if __name__ == '__main__':
    for x in ParserFactory("../data/P-Body-Enrichment_Class1.sorted.bed").parser:
        print(GInterval.to_bed_format_string(x))
