import os

from gparser.bed import BEDParser
from gparser.gtf import GTFParser


class ParserFactory(object):
    _FORMATS_MAP = {
        '.BED': "BED",
        '.FQ': 'FASTQ',
        '.FASTQ': 'FASTQ',
        '.FA': 'FASTA',
        '.FASTA': 'FASTA',
        '.GTF': 'GTF',
        ".BED12": "BED",
        ".BED6": "BED",
        ".BED8": "BED"
    }

    def __init__(self, path, fmt=None):
        self.path = path
        self.fmt = fmt

    @property
    def parser(self):
        """Parser: The iterable gparser which creates :class:`GInterval` instance
        for format-specific file.

        Raises:
            TypeError: The format is unknown or not supported.
        """

        if self.fmt is None:
            self.fmt = ParserFactory.infer_format_by_suffix(self.path)
        if self.fmt is None:
            raise TypeError('Format is None and we can not infer the format!')

        fmt = self.fmt
        parser = None
        if fmt == 'BED':
            parser = BEDParser(self.path)
        elif fmt == 'GTF':
            parser = GTFParser(self.path)

        if parser is None:
            raise TypeError('We can not find proper Parser for %s format.' % fmt)

        return parser

    @staticmethod
    def infer_format_by_suffix(path):
        """
        Inferring the format of file by their file name suffix.
        :param path:
        :return:
        """
        fmt = None
        x = os.path.splitext(path)[1]
        if x != '':
            fmt = ParserFactory._FORMATS_MAP.get(x.upper())
        return fmt

    @staticmethod
    def infer_format_by_context(path):
        """
        Inferring the format of file by their context.
        :param path:
        :return:
        """
        return None
