
import unittest
from GInterval.Interval import Interval, IntervalList


class GInterval(Interval):
    def __init__(self, chrom=None, x=0, y=0, strand='+', *args, **kwargs):
        super(GInterval, self).__init__(x, y)
        self._chrom = chrom
        self._strand = strand
        self._attributes = kwargs

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

    @property
    def attributes(self):
        return self._attributes

    def get_sub_interval(self, x, y):
        interval = super(GInterval, self).get_sub_interval(x, y)
        if interval is None:
            return None
        return GInterval(self._chrom, interval.x, interval.y, self.strand, **self._attributes)

    def get_start_position(self, strand=True):
        if strand and self._strand == '-':
            return self._y - 1
        else:
            return self._x

    def get_end_position(self, strand=True):
        if strand and self._strand == '-':
            return self._x
        else:
            return self._y - 1

    def get_index(self, position, strand=True):
        if not self._x <= position < self._y:
            raise Exception("Position [%d] out of valid position range [%d:%d]." % (position, self._x, self._y))
        index = abs(position - self.get_start_position(strand))
        return index

    def get_position(self, index, strand=True):
        length = self.length
        max_index = length - 1
        min_index = -length
        if not min_index <= index <= max_index:
            raise IndexError("Index [%d] out of valid index range [%d:%d]." % (index, min_index, max_index))
        if index < 0:
            index += length
        start_position = self.get_start_position(strand)
        if strand and self._strand == '-':
            position = start_position - index
        else:
            position = start_position + index
        return position

    def __getitem__(self, item):
        return self._attributes.get(item)


class GMultiInterval(GInterval):
    def __init__(self, chrom=None, x=0, y=0, strand=None, interval_list=None, *args, **kwargs):
        super(GMultiInterval, self).__init__(chrom, x, y, strand, *args, **kwargs)
        self._interval_list = interval_list
        if self._interval_list is None:
            self._interval_list = IntervalList([Interval(x, y)])

    @property
    def length(self):
        return self._interval_list.length

    def get_index(self, position, strand=True):
        index = 0
        for interval in self._interval_list:
            if interval.y <= position:
                index += interval.length
                continue
            elif interval.x <= position < interval.y:
                index += position - interval.x
                break
            else:
                raise Exception("Position [%d] is not inside this GMultiInterval" % position)
        if strand and self.reverse:
            index = self.length - index - 1
        return index

    def get_position(self, index, strand=True):
        length = self.length
        min_index = -length
        max_index = length - 1
        if not min_index <= index <= max_index:
            raise IndexError("Index [%d] out of valid index range [%d:%d]." % (index, min_index, max_index))
        if index < 0:
            index += length
        if strand and self.reverse:
            index = length - index - 1
        index0 = 0
        index1 = 0
        for interval in self._interval_list:
            index0 = index1
            index1 = index0 + interval.length
            if index0 <= index < index1:
                return index - index0 + interval.x
        raise Exception('Invalid index value [%d]' % index)

    def get_sub_interval(self, x, y):
        intervals = list(filter(lambda a: a is not None,
                                [interval.get_sub_interval(x, y) for interval in self._interval_list]))
        if len(intervals) == 0:
            return None
        x0 = intervals[0].x
        y0 = intervals[-1].y
        return GMultiInterval(chrom=self._chrom, x=x0, y=y0,
                              strand=self._strand, interval_list=IntervalList(intervals), **self._attributes)


class GComplexInterval(GMultiInterval):
    def __init__(self, chrom=None, x=0, y=0, strand=None, interval_list=None, thick_interval=None, *args, **kwargs):
        super(GComplexInterval, self).__init__(chrom=chrom, x=x, y=y, strand=strand,
                                               interval_list=interval_list, *args, **kwargs)
        self._thick_interval = thick_interval

    def get_sub_interval(self, x, y):
        intervals = list(filter(lambda a: a is not None,
                                [interval.get_sub_interval(x, y) for interval in self._interval_list]))
        if len(intervals) == 0:
            return None
        x0 = intervals[0].x
        y0 = intervals[-1].y

        thick_interval = self._thick_interval.get_sub_interval(x, y)

        return GComplexInterval(chrom=self._chrom, x=x0, y=y0, strand=self._strand, thick_interval=thick_interval,
                                interval_list=IntervalList(intervals), **self._attributes)


class GIntervalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ginterval = GInterval('chr1', 20, 50, '-', score=999, color='red')
        cls.interval_list = IntervalList([
            Interval(20, 30),
            Interval(40, 50),
            Interval(60, 70),
        ])
        cls.gminterval = GMultiInterval('chr2', 20, 50, '-', interval_list=cls.interval_list, score=999, color='red')

    def test_ginterval(self):
        self.assertEqual(self.ginterval.chrom, 'chr1')
        self.assertEqual(self.ginterval.x, 20)
        self.assertEqual(self.ginterval.y, 50)
        self.assertEqual(self.ginterval.get_start_position(True), 49)
        self.assertEqual(self.ginterval.length, 30)
        self.assertFalse(self.ginterval.forward)
        self.assertTrue(self.ginterval.reverse)
        self.assertEqual(self.ginterval.get_index(22), 27)
        self.assertEqual(self.ginterval.get_index(49), 0)
        self.assertEqual(self.ginterval.get_index(22, False), 2)
        self.assertEqual(self.ginterval.get_index(49, False), 29)
        self.assertEqual(self.ginterval.get_position(0), 49)
        self.assertEqual(self.ginterval.get_position(29), 20)
        self.assertEqual(self.ginterval.get_position(0, False), 20)
        self.assertEqual(self.ginterval.get_position(29, False), 49)
        self.assertRaises(Exception, self.ginterval.get_index, 50)
        self.assertEqual(self.ginterval.score, 999)
        self.assertEqual(self.ginterval['color'], 'red')

    def test_gminterval(self):
        self.assertEqual(self.gminterval.chrom, 'chr2')
        self.assertEqual(self.gminterval.length, 30)
        self.assertEqual(self.gminterval['color'], 'red')
        self.assertRaises(Exception, self.gminterval.get_index, 30)

        self.assertEqual(self.gminterval.get_index(29, False), 9)
        self.assertEqual(self.gminterval.get_index(40, False), 10)
        self.assertEqual(self.gminterval.get_index(49, False), 19)
        self.assertEqual(self.gminterval.get_index(69, False), 29)

        self.assertEqual(self.gminterval.get_index(29), 20)
        self.assertEqual(self.gminterval.get_index(40), 19)
        self.assertEqual(self.gminterval.get_index(49), 10)
        self.assertEqual(self.gminterval.get_index(69), 0)

        self.assertEqual(self.gminterval.get_position(0), 69)
        self.assertEqual(self.gminterval.get_position(9), 60)
        self.assertEqual(self.gminterval.get_position(10), 49)
        self.assertEqual(self.gminterval.get_position(20), 29)
        self.assertEqual(self.gminterval.get_position(29), 20)
        self.assertEqual(self.gminterval.get_position(-15), 44)

        self.assertEqual(self.gminterval.get_position(0, False), 20)
        self.assertEqual(self.gminterval.get_position(9, False), 29)
        self.assertEqual(self.gminterval.get_position(10, False), 40)
        self.assertEqual(self.gminterval.get_position(20, False), 60)
        self.assertEqual(self.gminterval.get_position(29, False), 69)
        self.assertEqual(self.gminterval.get_position(-15, False), 45)

    def test_gcinterval(self):
        pass


if __name__ == '__main__':
    unittest.main()
