
class Interval(object):
    def __init__(self, start=0, end=0):
        self._start = start
        self._end = end

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def length(self):
        return self._end - self._start + 1

    def __len__(self):
        return self.length


class IntervalList(object):
    def __init__(self, intervals=None):
        self._intervals = intervals
        if self._intervals is None:
            self._intervals = []

    def __getitem__(self, item):
        return self._intervals[item]

    def __iter__(self):
        return self._intervals.__iter__()

    def __len__(self):
        return len(self._intervals)

    @property
    def length(self):
        return sum([interval.length for interval in self._intervals])


class GInterval(Interval):
    def __init__(self, chrom=None, start=0, end=0, strand='+', *args, **kwargs):
        super(GInterval, self).__init__(start, end)
        self._chrom = chrom
        self._strand = strand
        self._attributes = kwargs

    def __getitem__(self, item):
        return self._attributes.get(item)

    @property
    def chrom(self):
        return self._chrom

    @property
    def strand(self):
        return self._strand

    @property
    def forward(self):
        return self._strand == '+'

    @property
    def reverse(self):
        return self._strand == '-'

    @property
    def score(self):
        return self._attributes.get('score', None)

    def get_start(self, strand=True):
        if strand and self._strand == '-':
            return self._end
        else:
            return self._start

    def get_end(self, strand=True):
        if strand and self._strand == '-':
            return self._start
        else:
            return self._end

    def get_index(self, position, strand=True):
        if not self._start <= position <= self._end:
            raise Exception("Position [%d] out of valid position range [%d:%d]." % (position, self._start, self._end))
        if strand and self._strand == '-':
            index = self._end - position
        else:
            index = position - self._start
        return index

    def get_position(self, index, strand=True):
        length = self.length
        max_index = length - 1
        min_index = -length
        if not min_index <= index <= max_index:
            raise IndexError("Index [%d] out of valid index range [%d:%d]." % (index, min_index, max_index))
        if index < 0:
            index += length
        if strand and self._strand == '-':
            position = self._end - index
        else:
            position = index + self._start
        return position


class GMultiInterval(GInterval):
    def __init__(self, chrom=None, start=0, end=0, strand=None, interval_list=None, *args, **kwargs):
        super(GMultiInterval, self).__init__(chrom, start, end, strand, *args, **kwargs)
        self._interval_list = interval_list
        if self._interval_list is None:
            self._interval_list = IntervalList()

    @property
    def length(self):
        return self._intervals.length

    def get_index(self, position, strand=True):
        index = 0
        for interval in self._interval_list:
            if interval.end < position:
                index += interval.length
                continue
            elif interval.start <= position <= interval.end:
                index += position - interval.start
            else:
                raise Exception("Position [%d] is not inside this GMultiInterval" % position)
        if strand and self.reverse:
            return self.length - index

    def get_position(self, index, strand=True):
        pass


class GComplexInterval(GMultiInterval):
    def __init__(self):
        super(GComplexInterval, self).__init__()
        self._thick_interval = None


class BedReader(object):
    pass


if __name__ == '__main__':
    gi = GInterval(start=2, end=5, strand='-')
    print(gi.get_position(5, strand=False))


