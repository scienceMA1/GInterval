import unittest


class Interval(object):
    def __init__(self, x=0, y=0):
        assert x is not None
        assert y is not None
        assert y > x >= 0

        self.x = x
        self.y = y

    def get_start_position(self, forward=True):
        position = None
        if forward:
            position = self.x
        else:
            position = self.y - 1
        return position

    def get_end_position(self, forward=True):
        position = None
        if forward:
            position = self.y - 1
        else:
            position = self.x
        return position

    @property
    def length(self):
        return self.y - self.x

    def __len__(self):
        return self.length

    def __str__(self):
        return 'Start: %d, End: %d, Length: %d' % (self.x, self.y, self.length)

    def __getitem__(self, item):
        if isinstance(item, slice):
            x = max(self.x, slice.start)
            y = min(self.y, slice.stop)
            if x >= y:
                return None
            return Interval(x, y)
        raise ValueError("Only slice supported!")


class GInterval(Interval):
    def __init__(self, x, y, chrom=None, name=None, strand=None, blocks=None, thick=None, *args, **kwargs):
        super(GInterval, self).__init__(x, y)
        self.chrom = chrom
        self.name = name
        self.strand = strand
        self.blocks = blocks
        self.thick = thick
        self.attributes = kwargs

    @property
    def forward(self):
        return self.strand == '+'

    @property
    def reverse(self):
        return self.strand == '-'

    def get_start_position(self, strand_sensitive=True):
        if strand_sensitive and self.reverse:
            return self.y - 1
        return self.x

    def get_end_position(self, strand_sensitive=True):
        if strand_sensitive and self.reverse:
            return self.x
        return self.y - 1

    def get_index(self, position, strand_sensitive=True, block_sensitive=True):
        if self.blocks is None:
            block_sensitive = False
        if block_sensitive:
            index = 0
            for block in self.blocks:
                if block.y <= position:
                    index += len(block)
                    continue
                elif block.x <= position < block.y:
                    index += position - block.x
                    break
                else:
                    raise ValueError("Position [%d] is not hit in block!" % position)
            if strand_sensitive and self.reverse:
                index = self.get_length(block_sensitive) - 1 - index
            return index
        else:
            if not self.x <= position < self.y:
                raise ValueError("Invalid position: %d, valid position range: [%d:%d]" % (position, self.x, self.y - 1))
            return abs(position - self.get_start_position(strand_sensitive))

    def get_position(self, index, strand_sensitive=True, block_sensitive=True):
        if self.blocks is None:
            block_sensitive = False
        length = self.get_length(block_sensitive)
        min_index = -length
        max_index = length - 1
        if not min_index <= index <= max_index:
            raise ValueError("Invalid index [%d], valid index range: [%d:%d]" % (index, min_index, max_index))
        if index < 0:
            index += length
        if block_sensitive:
            if strand_sensitive and self.reverse:
                index = length - 1 - index
            index0 = 0
            index1 = 0
            position = None
            for block in self.blocks:
                index0 += index1
                index1 = index0 + len(block) - 1
                if index > index1:
                    continue
                elif index0 <= index <= index1:
                    position = block.x + index - index0
                    break
                else:
                    break
            if position is None:
                raise Exception('Unknown Exception!')
            return position




        else:
            start = self.get_start_position(strand_sensitive)
            if strand_sensitive and self.reverse:
                return start - index
            else:
                return start + index

    def get_length(self, block_sensitive=True):
        if block_sensitive and self.blocks is not None:
            return sum([len(block) for block in self.blocks])
        return len(self)

    def get_max_index(self, block_sensitive=True):
        return self.get_length(block_sensitive=block_sensitive) - 1

    @staticmethod
    def get_min_index():
        return 0

    def get_thick_length(self, block_sensitive=True):
        index1 = self.get_index(self.thick.get_start_position(), block_sensitive=block_sensitive)
        index2 = self.get_index(self.thick.get_end_position(), block_sensitive=block_sensitive)
        return abs(index1 - index2) + 1

    def get_thin_length_1(self, strand_sensitive=True, block_sensitive=True):
        position = self.thick.get_start_position(self.forward | (not strand_sensitive))
        index = self.get_index(self.thick.get_start_position(self.forward | (not strand_sensitive)),
                               strand_sensitive=strand_sensitive, block_sensitive=block_sensitive)
        return index

    def get_thin_length_2(self, strand_sensitive=True, block_sensitive=True):
        index = self.get_index(self.thick.get_end_position(self.forward | (not strand_sensitive)),
                               strand_sensitive=strand_sensitive, block_sensitive=block_sensitive)
        return self.get_max_index(block_sensitive=block_sensitive) - index

    def get_sub_interval(self, x, y):
        pass

    def __getitem__(self, item):
        if isinstance(item, slice):
            return self.get_sub_interval(slice.start, slice.stop)
        else:
            return self.attributes[item]

    @classmethod
    def to_bed_format_string(cls, gi, ncol=12):
        cols = []
        if ncol >= 4:
            cols.extend([gi.chrom, gi.x, gi.y, gi.name])
        if ncol >= 6:
            score = gi['score']
            if score is None:
                score = '.'
            cols.extend([gi['score'], gi.strand])
        if ncol >= 12:
            tx = '.'
            ty = '.'
            rgb = gi['rgb']
            if rgb is None:
                rgb = '.'
            bnum = 0
            sizes = ''
            starts = ''
            if gi.thick is not None:
                tx = gi.thick.x
                ty = gi.thick.y
            if gi.blocks is not None:
                bnum = len(gi.blocks)
                for block in gi.blocks:
                    sizes += '%d,' % len(block)
                    starts += '%d,' % (block.x - gi.x)
            cols.extend([tx, ty, rgb, bnum, sizes, starts])
        return "\t".join(map(str, cols))

    @classmethod
    def to_gtf_format_string(cls, gi):
        pass


class IntervalTest(unittest.TestCase):
    def test_Interval(self):
        interval = Interval(10, 20)
        self.assertEqual(interval.length, 10)
        self.assertEqual(interval.get_start_position(), 10)
        self.assertEqual(interval.get_end_position(), 19)
        self.assertEqual(interval.get_start_position(False), 19)
        self.assertEqual(interval.get_end_position(False), 10)

    def test_GInterval(self):
        intervals = [Interval(10, 20), Interval(30, 40), Interval(45, 50)]
        gi1 = GInterval(10, 50, 'chr1', 'GInterval1', '-')
        gi2 = GInterval(10, 50, 'chr2', 'GInterval2', '-', blocks=intervals)

        self.assertEqual(gi1.get_length(True), 40)
        self.assertEqual(gi1.get_length(False), 40)
        self.assertEqual(gi2.get_length(True), 25)
        self.assertEqual(gi2.get_length(False), 40)

        self.assertEqual(gi1.get_index(10), 39)
        self.assertEqual(gi1.get_index(20), 29)
        self.assertEqual(gi1.get_index(40), 9)
        self.assertEqual(gi1.get_index(49), 0)

        self.assertEqual(gi2.get_index(10), 24)
        self.assertEqual(gi2.get_index(19), 15)
        self.assertEqual(gi2.get_index(39), 5)
        self.assertEqual(gi2.get_index(49), 0)

        self.assertEqual(gi2.get_index(10, strand_sensitive=False), 0)
        self.assertEqual(gi2.get_index(19, strand_sensitive=False), 9)
        self.assertEqual(gi2.get_index(39, strand_sensitive=False), 19)
        self.assertEqual(gi2.get_index(49, strand_sensitive=False), 24)

        self.assertEqual(gi2.get_index(10, strand_sensitive=False, block_sensitive=False), 0)
        self.assertEqual(gi2.get_index(19, strand_sensitive=False, block_sensitive=False), 9)
        self.assertEqual(gi2.get_index(39, strand_sensitive=False, block_sensitive=False), 29)
        self.assertEqual(gi2.get_index(49, strand_sensitive=False, block_sensitive=False), 39)

    def test_GInterval_Thick(self):
        intervals = [Interval(10, 20), Interval(30, 40), Interval(45, 50)]
        thick = Interval(15, 35)
        gi = GInterval(10, 50, 'chr2', 'GInterval2', '-', blocks=intervals, thick=thick)

        self.assertEqual(gi.get_thick_length(block_sensitive=True), 10)
        self.assertEqual(gi.get_thick_length(block_sensitive=False), 20)

        self.assertEqual(gi.get_thin_length_1(), 10)
        self.assertEqual(gi.get_thin_length_2(), 5)
        self.assertEqual(gi.get_thin_length_1(block_sensitive=False), 15)
        self.assertEqual(gi.get_thin_length_2(block_sensitive=False), 5)

        self.assertEqual(gi.get_thin_length_1(strand_sensitive=False), 5)
        self.assertEqual(gi.get_thin_length_2(strand_sensitive=False), 10)


if __name__ == '__main__':
    unittest.main()
