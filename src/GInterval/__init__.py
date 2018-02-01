import unittest


class Interval(object):
    def __init__(self, x, y):
        assert y > x >= 0
        self.x = x
        self.y = y

    def __len__(self):
        return self.y - self.x

    def __str__(self):
        return 'x: %d, y: %d, len: %d' % (self.x, self.y, len(self))

    def __getitem__(self, item):
        if isinstance(item, slice):
            x = max(self.x, slice.start)
            y = min(self.y, slice.stop)
            if x >= y:
                return None
            return Interval(x, y)
        raise ValueError("Only slice supported!")


class GInterval(Interval):
    def __init__(self, x, y, chrom, name, strand='+', **kwargs):
        super(GInterval, self).__init__(x, y)
        self.chrom = chrom
        self.name = name
        self.strand = strand
        self._strand_backup = self.strand
        self.attribute = kwargs

    def __getitem__(self, item):
        return self.attribute.get(item)

    def set_strand_sensitivity(self, sensitive=True):
        if sensitive:
            self.strand = self._strand_backup
        else:
            self.strand = '+'

    @property
    def forward(self):
        return self.strand == '+'

    @property
    def reverse(self):
        return self.strand == '-'

    @property
    def start_position(self):
        if self.forward:
            return self.x
        else:
            return self.y - 1

    @property
    def end_position(self):
        if self.forward:
            return self.y - 1
        else:
            return self.x

    @property
    def start_index(self):
        return 0

    @property
    def end_index(self):
        return len(self) - 1

    @property
    def center_index(self):
        return len(self) // 2

    @property
    def center_position(self):
        return self.position(self.center_index)

    def index(self, position):
        return abs(position - self.start_position)

    def position(self, index):
        # start_index = self.start_index
        end_index = self.end_index
        min_index = -(end_index + 1)
        if not min_index <= index <= end_index:
            raise ValueError('Invalid index %d' % index)
        if index < 0:
            index = index - min_index
        if self.reverse:
            index = end_index - index
        return self.x + index


class GMultiInterval(GInterval):
    def __init__(self, x, y, chrom, name, strand='+', blocks=None, **kwargs):
        super(GMultiInterval, self).__init__(x, y, chrom, name, strand, **kwargs)
        self._blocks = None
        if blocks is None:
            self._blocks = [y - x]
        elif isinstance(blocks, list):
            self._blocks = list()
            last_block = None
            for b in blocks:
                if last_block is not None:
                    self._blocks.append(b.x - last_block.y)
                self._blocks.append(b.y - b.x)
                last_block = b

    @property
    def blocks(self):
        x = self.x
        # y = 0
        for i, v in enumerate(self._blocks):
            if i % 2 == 0:
                y = x + v
                yield Interval(x, y)
                x = y
            else:
                x += v

    @property
    def block_count(self):
        return (len(self._blocks) + 1) // 2

    def __len__(self):
        length = 0
        for i in range(0, len(self._blocks), 2):
            length += self._blocks[i]
        return length

    def index(self, position):
        if not self.x <= position < self.y:
            raise ValueError('Invalid position %d' % position)
        index = 0
        for block in self.blocks:
            if block.y <= position:
                index += len(block)
                continue
            elif block.x <= position < block.y:
                index += position - block.x
                break
            else:
                raise ValueError('')
        if self.reverse:
            index = len(self) - index - 1
        return index

    def position(self, index):
        end_index = self.end_index
        min_index = -(end_index + 1)
        if not min_index <= index <= end_index:
            raise ValueError('Invalid index %d' % index)
        if index < 0:
            index -= min_index
        if self.reverse:
            index = end_index - index
        for block in self.blocks:
            index -= len(block)
            if index < 0:
                return block.y + index
        raise Exception('Unknown Exception!')


class GComplexInterval(GMultiInterval):
    def __init__(self, x, y, chrom, name, strand='+', blocks=None, thick=None, **kwargs):
        super(GComplexInterval, self).__init__(x, y, chrom, name, strand, blocks, **kwargs)
        self.thick = None
        if thick is not None:
            self.thick = thick

    @property
    def thick_start_position(self):
        if self.forward:
            return self.thick.x
        else:
            return self.thick.y - 1

    @property
    def thick_end_position(self):
        if self.forward:
            return self.thick.y - 1
        else:
            return self.thick.x

    @property
    def thick_start_index(self):
        return self.index(self.thick_start_position)

    @property
    def thick_end_index(self):
        return self.index(self.thick_end_position)

    def has_thick(self):
        return self.thick is not None

    def to_bed_format_string(self, ncol=12):
        cols = []
        if ncol >= 4:
            cols.extend([self.chrom, self.x, self.y, self.name])
        if ncol >= 6:
            score = self['score']
            if score is None:
                score = '.'
            cols.extend([score, self.strand])
        if ncol >= 12:
            tx = '.'
            ty = '.'
            rgb = self['rgb']
            if rgb is None:
                rgb = '.'
            bnum = 0
            sizes = ''
            starts = ''
            if self.thick is not None:
                tx = self.thick.x
                ty = self.thick.y
            for block in self.blocks:
                bnum += 1
                sizes += '%d,' % len(block)
                starts += '%d,' % (block.x - self.x)
            cols.extend([tx, ty, rgb, bnum, sizes, starts])
        return "\t".join(map(str, cols))


class IntervalTest(unittest.TestCase):
    def test_Interval(self):
        i = Interval(10, 20)
        self.assertEqual(len(i), 10)
        self.assertEqual(i.x, 10)
        self.assertEqual(i.y, 20)

    def test_GInterval(self):
        gi = GInterval(10, 50, 'chr1', 'GInterval1', '-')

        self.assertEqual(len(gi), 40)

        gi.set_strand_sensitivity(False)
        self.assertEqual(gi.start_position, 10)
        self.assertEqual(gi.end_position, 49)
        self.assertEqual(gi.start_index, 0)
        self.assertEqual(gi.end_index, 39)

        gi.set_strand_sensitivity(True)
        self.assertEqual(gi.start_position, 49)
        self.assertEqual(gi.end_position, 10)
        self.assertEqual(gi.start_index, 0)
        self.assertEqual(gi.end_index, 39)

    def test_GMultiInterval(self):
        intervals = [Interval(10, 20), Interval(30, 35), Interval(50, 65)]
        gi = GMultiInterval(10, 65, 'chr2', 'GInterval2', '-', blocks=intervals)

        self.assertEqual(gi.start_index, 0)
        self.assertEqual(gi.end_index, 29)
        self.assertEqual(len(gi), 30)

        self.assertEqual(gi.index(10), 29)
        gi.set_strand_sensitivity(False)
        self.assertEqual(gi.index(15), 5)
        gi.set_strand_sensitivity(True)
        self.assertEqual(gi.index(15), 24)

        self.assertEqual(gi.position(29), 10)

    def test_GComplexInterval(self):
        intervals = [Interval(10, 20), Interval(30, 35), Interval(50, 65)]
        thick = Interval(15, 60)
        gi = GComplexInterval(10, 65, 'chr2', 'GInterval2', '-', blocks=intervals, thick=thick)

        self.assertEqual(gi.thick_start_position, 59)
        self.assertEqual(gi.thick_end_position, 15)
        self.assertEqual(gi.thick_start_index, 5)
        self.assertEqual(gi.thick_end_index, 24)


if __name__ == '__main__':
    unittest.main()
