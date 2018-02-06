from ginterval.parser import Parser
from ginterval import GInterval

class BEDParser(Parser):
    """
    BEDParser is a iterable parser which create :class:`ginterval.GInterval` instance
    from BED format file. The lines in file start with "#" will be
    ignored.

    BED4, BED6 and BED12 are supported. The parser will infers
    which is the real format according to columns number and decides
    how to create GInterval instance. If the provided ncol is less
    than real column number, the residual columns will be ignored.

    Attributes:
        path (str): The path of BED file.
        ncol (int, optional): The max column number for each record.

    Examples:
        This example will demonstrates how to obtain GInterval
        instance from BED format file. If the real column number
        is less than ncol, parser will create GInterval instance
        by real column number.

        >>> parser = BEDParser("test/test.bed", ncol=12)
        >>> for gi in parser:
        >>>     print(gi.to_bed_format_string())

        This example will demonstrates that if the real column number
        is larger than ncol, the extra will be ignored.

        >>> parser = BEDParser("test/test.bed", ncol=4)
        >>> for gi in parser:
        >>>     print(gi.to_bed_format_string())

    """

    BED4 = 4
    BED6 = 6
    BED12 = 12

    def __init__(self, path, ncol=BED12):
        super(BEDParser, self).__init__(path)
        self.__ncol = ncol

    def __iter__(self):
        is_bed6 = self.__ncol >= BEDParser.BED6
        is_bed12 = self.__ncol >= BEDParser.BED12

        with open(self.path) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                cols = line.strip("\n").split("\t")
                ncol = len(cols)
                if ncol < 4:
                    raise TypeError("Invalid BED format file. Record least than 4 columns!")

                parameters = dict()
                parameters["chrom"] = cols[0]
                start = int(cols[1])
                end = int(cols[2])
                parameters["x"] = start
                parameters["y"] = end
                parameters["name"] = cols[3]

                if ncol > 4 and is_bed6:
                    if cols[4] != ".":
                        parameters["score"] = float(cols[4])
                    else:
                        parameters["score"] = None
                    parameters['strand'] = cols[5]

                if ncol > 6 and is_bed12:
                    thick_x = int(cols[6])
                    thick_y = int(cols[7])
                    parameters["thick"] = None
                    if thick_y > thick_x:
                        parameters["thick"] = (thick_x, thick_y)
                    parameters["rgb"] = cols[8]
                    sizes = map(int, filter(lambda obj: obj != '', cols[10].split(',')))
                    starts = map(int, filter(lambda obj: obj != '', cols[11].split(',')))
                    blocks = []
                    for bstart, bsize in zip(starts, sizes):
                        x = start + bstart
                        y = x + bsize
                        blocks.append((x, y))
                    parameters["blocks"] = blocks

                yield GInterval(**parameters)
