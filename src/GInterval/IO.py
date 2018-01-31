from GInterval import GInterval, Interval
#from GInterval.Interval import Interval, IntervalList


def parser(infile, fmt):
    fmt = fmt.upper()
    if fmt == 'BED':
        return __parser_for_bed(infile)
    elif fmt == 'GTF':
        pass
    elif fmt == 'GFF':
        pass


def __parser_for_bed(infile):
    with open(infile) as f:
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


def __parser_for_gtf(infile):
    pass


def __parser_for_gff(infile):
    pass


def writer(f, interval, fmt):
    pass


if __name__ == '__main__':
    for record in parser("../data/P-Body-Enrichment_Class1.sorted.bed", "bed"):
        print(GInterval.to_bed_format_string(record))
