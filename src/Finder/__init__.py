from GInterval import GInterval
from abc import abstractmethod


class SpecificSiteInterval(GInterval):
    def __init__(self, chrom=None, x=0, y=0, strand='+', *args, **kwargs):
        super(SpecificSiteInterval, self).__init__(chrom, x, y, strand, *args, **kwargs)


class SpecificSiteFinder(object):
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

    @abstractmethod
    def process(self, *args, **kwargs):
        raise NotImplementedError("This method must be implemented by sub-class.")

    def write(self, specific_site_interval):
        if self.fw is not None:
            self.fw.write(specific_site_interval + "\n")