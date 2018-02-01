import unittest


class Interval(object):
    def __init__(self, x, y):
        assert y > x >= 0
        self._x = x
        self._y = y

    '''def get_start_position(self, forward=True):
        position = None
        if forward:
            position = self._x
        else:
            position = self._y - 1
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
        return self.y - self.x'''

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    def __len__(self):
        return self._y - self._x

    def __str__(self):
        return 'x: %d, y: %d, len: %d' % (self._x, self._y, len(self))

    def __getitem__(self, item):
        if isinstance(item, slice):
            x = max(self._x, slice.start)
            y = min(self._y, slice.stop)
            if x >= y:
                return None
            return Interval(x, y)
        raise ValueError("Only slice supported!")


class GInterval(Interval):
    def __init__(self, x, y, chrom, name, strand='+', **kwargs):
        super(GInterval, self).__init__(x, y)
        self._chrom = chrom
        self._name = name
        self._strand = strand
        self._strand_backup = self._strand
        self._attribute = kwargs

    def __getitem__(self, item):
        return self._attribute.get(item)

    def set_strand_sensitivity(self, sensitive=True):
        if sensitive:
            self._strand = self._strand_backup
        else:
            self._strand = '+'

    @property
    def strand(self):
        return self._strand

    @strand.setter
    def strand(self, strand):
        self._strand = strand

    @property
    def forward(self):
        return self.strand == '+'

    @property
    def reverse(self):
        return self.strand == '-'

    @property
    def chrom(self):
        return self._chrom

    @property
    def name(self):
        return self._name

    @property
    def start_position(self):
        if self.forward:
            return self._x
        else:
            return self._y - 1

    @property
    def end_position(self):
        if self.forward:
            return self._y - 1
        else:
            return self._x

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
        return self._x + index


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
        x = self._x
        # y = 0
        for i, v in enumerate(self._blocks):
            if i % 2 == 0:
                y = x + v
                yield Interval(x, y)
                x = y
            else:
                x += v

    def __len__(self):
        length = 0
        for i in range(0, len(self._blocks), 2):
            length += self._blocks[i]
        return length

    def index(self, position):
        if not self._x <= position < self._y:
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
        self._thick = None
        if thick is not None:
            self._thick = thick

    @property
    def thick_start_position(self):
        if self.forward:
            return self._thick.x
        else:
            return self._thick.y - 1

    @property
    def thick_end_position(self):
        if self.forward:
            return self._thick.y - 1
        else:
            return self._thick.x

    @property
    def thick_start_index(self):
        return self.index(self.thick_start_position)

    @property
    def thick_end_index(self):
        return self.index(self.thick_end_position)

    @property
    def thick(self):
        return self._thick

    def has_thick(self):
        return self._thick is not None

    def to_bed_format_string(self, ncol=12):
        cols = []
        if ncol >= 4:
            cols.extend([self._chrom, self._x, self._y, self._name])
        if ncol >= 6:
            score = self['score']
            if score is None:
                score = '.'
            cols.extend([score, self._strand])
        if ncol >= 12:
            tx = '.'
            ty = '.'
            rgb = self['rgb']
            if rgb is None:
                rgb = '.'
            bnum = 0
            sizes = ''
            starts = ''
            if self._thick is not None:
                tx = self._thick.x
                ty = self._thick.y
            for block in self.blocks:
                bnum += 1
                sizes += '%d,' % len(block)
                starts += '%d,' % (block.x - self.x)
            cols.extend([tx, ty, rgb, bnum, sizes, starts])
        return "\t".join(map(str, cols))


'''
class GInterval0(Interval):
    def __init__(self, x, y, chrom=None, name=None, strand=None, blocks=None, thick=None, *args, **kwargs):
        super(GInterval0, self).__init__(x, y)
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

    def get_thick_start_position(self, strand_sensitive=True):
        return self.thick.get_start_position(self.forward | (not strand_sensitive))

    def get_thick_end_position(self, strand_sensitive=True):
        return self.thick.get_end_position(self.forward | (not strand_sensitive))

    def get_thick_start_index(self, strand_sensitive=True, block_sensitive=True):
        return self.get_index(self.get_thick_start_position(strand_sensitive), block_sensitive=block_sensitive)

    def get_thick_end_index(self, strand_sensitive=True, block_sensitive=True):
        return self.get_index(self.get_thick_end_position(strand_sensitive), block_sensitive=block_sensitive)

    def get_thick_length(self, block_sensitive=True):
        # index1 = self.get_index(self.thick.get_start_position(), block_sensitive=block_sensitive)
        # index2 = self.get_index(self.thick.get_end_position(), block_sensitive=block_sensitive)
        # return abs(index1 - index2) + 1
        return self.get_thick_end_index(block_sensitive=block_sensitive) - \
               self.get_thick_start_index(block_sensitive=block_sensitive) + 1

    def get_thin_length_1(self, strand_sensitive=True, block_sensitive=True):
        index = self.get_index(self.thick.get_start_position(self.forward | (not strand_sensitive)),
                               strand_sensitive=strand_sensitive, block_sensitive=block_sensitive)
        return index

    def get_thin_length_2(self, strand_sensitive=True, block_sensitive=True):
        index = self.get_index(self.thick.get_end_position(self.forward | (not strand_sensitive)),
                               strand_sensitive=strand_sensitive, block_sensitive=block_sensitive)
        return self.get_max_index(block_sensitive=block_sensitive) - index

    def get_center_position(self, strand_sensitive=True, block_sensitive=True):
        return self.get_position(self.get_center_index(block_sensitive), strand_sensitive, block_sensitive)

    def get_center_index(self, block_sensitive=True):
        return self.get_length(block_sensitive) // 2

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
'''


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
