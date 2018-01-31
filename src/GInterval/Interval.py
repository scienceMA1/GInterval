import unittest


class Interval(object):
    """
    Interval class.

    0-base coordination.
    half-interval [x, y)

    """

    def __init__(self, x=0, y=0):
        assert x is not None
        assert y is not None
        assert y >= x >= 0

        self._x = x
        self._y = y

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def length(self):
        return self._y - self._x

    def get_sub_interval(self, x, y):
        x0 = max(self._x, x)
        y0 = min(self._y, y)
        if x0 >= y0:
            return None
        return Interval(x0, y0)

    def __len__(self):
        return self.length

    def __str__(self):
        return 'Start: %d, End: %d, Length: %d' % (self._x, self._y, self.length)


class IntervalList(object):
    def __init__(self, intervals=None):
        self._intervals = intervals
        if self._intervals is None:
            self._intervals = []

    @property
    def length(self):
        return sum([interval.length for interval in self._intervals])

    def append(self, interval):
        self._intervals.append(interval)

    def __getitem__(self, item):
        return self._intervals[item]

    def __iter__(self):
        return self._intervals.__iter__()

    def __len__(self):
        return len(self._intervals)

    def __str__(self):
        pass


class IntervalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.interval = Interval(2, 7)
        cls.intervals = [
            Interval(2, 3),
            Interval(7, 10),
            Interval(20, 30),
            Interval(200, 300)
        ]
        cls.interval_list = IntervalList(cls.intervals)

    def test_interval(self):
        self.assertRaises(AssertionError, Interval, x=2, y=1)
        self.assertRaises(AssertionError, Interval, x=-1, y=1)
        self.assertEqual(self.interval.x, 2)
        self.assertEqual(self.interval.y, 7)
        self.assertEqual(self.interval.length, 5)

    def test_interval_list(self):
        self.assertEqual(self.interval_list.length, 114)
        self.assertEqual(len(self.interval_list), 4)
        for a, b in zip(self.intervals, self.interval_list):
            self.assertEqual(a, b)


if __name__ == '__main__':
    unittest.main()
