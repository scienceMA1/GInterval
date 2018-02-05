import unittest


class GInterval(object):
    def __init__(self, x=None, y=None, blocks=None, thick=None, **kwargs):
        self._x = x
        self._y = y
        self.attribute = kwargs
        self._blocks = blocks
        self.thick = thick

        if self._x is None and self._y is None:
            self._x = self._blocks[0].x
            self._y = self._blocks[-1].y

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def name(self):
        return self['name']

    @property
    def chrom(self):
        return self['chrom']

    @property
    def strand(self):
        return self.attribute.get('strand', '+')

    @property
    def forward(self):
        return self.strand == '+'

    @property
    def reverse(self):
        return self.strand == '-'

    @property
    def blocks(self):
        if self._blocks is not None:
            return self._blocks
        else:
            return [self]

    @property
    def gaps(self):
        last_block = None
        results = list()
        blocks = self.combine_adjacent(self.blocks)
        for block in blocks:
            if last_block is None:
                last_block = block
            else:
                results.append(GInterval(last_block.y, block.x))
                last_block = block
        return results

    def set_strand_sensitivity(self, sensitivity=True):
        self.attribute.setdefault('strand_backup', self.attribute.get('strand'))
        if sensitivity:
            self.attribute['strand'] = self.attribute.get('strand_backup')
        else:
            self.attribute['strand'] = '+'

    def __getitem__(self, item):
        if isinstance(item, slice):
            x = max(self._x, item.start)
            y = min(self._y, item.stop)

            if x >= y:
                return None

            thick = None
            if self.thick is not None:
                thick = self.thick[x:y]

            blocks = self.blocks
            bnum = len(blocks)
            if bnum == 0:
                return None
            elif bnum == 1:
                return GInterval(x, y, thick=thick, **self.attribute)
            blocks = list(filter(lambda z: z is not None, [block[x:y] for block in blocks]))

            return GInterval(blocks=blocks, thick=thick, **self.attribute)
        return self.attribute.get(item)

    def __add__(self, other):
        """

        The strand of other will be ignored.

        :param other:
        :return:
        """
        last_block = None
        integrated = []
        for block in sorted(self.blocks + other.blocks, key=lambda z: z.x):
            if last_block is None:
                last_block = block
                continue
            else:
                x = max(last_block._x, block.x)
                y = min(last_block._y, block.y)
                if y > x:
                    last_block = GInterval(x, y)
                else:
                    integrated.append(last_block)
                    last_block = block
        if last_block is not None:
            integrated.append(last_block)
        thicks = list(filter(lambda z: z is not None, [self.thick, other.thick]))
        thick = None
        if len(thicks) > 0:
            thick = GInterval(min(thicks[0].x, thicks[-1].x), max(thicks[0].y, thicks[-1].y))
        return GInterval(blocks=integrated, thick=thick, **self.attribute)

    def __len__(self):
        bs = self.blocks
        if len(bs) == 1:
            b = bs[0]
            return b._y - b._x
        else:
            return sum([len(b) for b in bs])

    def overlap(self, other):
        return max(self._x, other.x) < min(self._y.other.y)
        '''
        blocks1 = list(other.blocks)
        blocks2 = list(self.blocks)
        i1 = 0
        i2 = 0
        mi1 = len(blocks1)
        mi2 = len(blocks2)
        while i1 < mi1 and i2 < mi2:
            b1 = blocks1[i1]
            b2 = blocks2[i2]
            if b1.y <= b2.x:
                i1 += 1
            elif b1.x >= b2.y:
                i2 += 1
            else:
                return True
        return False
        '''

    def contain(self, other):
        blocks1 = self.combine_adjacent(self.blocks)
        blocks2 = self.combine_adjacent(other.blocks)

        i1 = 0
        mi1 = len(blocks1)
        for b2 in blocks2:
            while True:
                if i1 >= mi1:
                    return False
                b1 = blocks1[i1]
                if b2.x >= b1.y:
                    i1 += 1
                    continue
                elif b2.x >= b1.x and b2.y <= b1.y:
                    break
                else:
                    return False
        return True

    def coincide(self, other):
        if not (self._x <= other.x and self._y >= other.y):
            return False

        if len(other.blocks) == 1:
            return self.contain(other)

        gaps1 = self.gaps
        gaps2 = other.gaps
        if len(gaps1) == len(gaps2) == 0:
            return True

        if len(gaps1) < len(gaps2):
            return False

        start = None
        g2 = gaps2[0]
        for i, g1 in enumerate(gaps1):
            if g1.x == g2.x and g1.y == g2.y:
                start = i
                break
        if start is None:
            return False

        residual = len(gaps2) - 1
        if len(gaps1) - start - 1 < residual:
            return False

        for i in range(residual):
            g1 = gaps1[start + 1 + i]
            g2 = gaps2[i + 1]
            if not (g1.x == g2.x and g1.y == g2.y):
                return False
        return True

    def index(self, position):
        if not self._x <= position < self._y:
            raise ValueError('Invalid position %d' % position)
        index = 0
        for block in self.blocks:
            if block._y <= position:
                index += len(block)
                continue
            elif block._x <= position < block._y:
                index += position - block._x
                break
            else:
                raise ValueError('')
        if self.reverse:
            index = len(self) - index - 1
        return index

    def position(self, index):
        end_index = len(self) - 1
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
                return block._y + index
        raise Exception('Unknown Exception!')

    # Additional methods.

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

    # Thick operation.

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

    # Output format

    def to_bed_format_string(self, ncol=12):
        cols = []
        if ncol >= 4:
            cols.extend([self.chrom, self._x, self._y, self.name])
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
                starts += '%d,' % (block._x - self._x)
            cols.extend([tx, ty, rgb, bnum, sizes, starts])
        return "\t".join(map(str, cols))

    @classmethod
    def combine_adjacent(cls, intervals):
        intervals = sorted(intervals, key=lambda z: z.x)
        results = []
        last_interval = None
        for interval in intervals:
            if last_interval is None:
                last_interval = interval
            else:
                if interval.x <= last_interval.y:
                    last_interval = GInterval(last_interval.x, max(last_interval.y, interval.y))
                else:
                    results.append(last_interval)
                    last_interval = interval
        if last_interval is not None:
            results.append(last_interval)
        return results


'''
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
            x = max(self.x, item.start)
            y = min(self.y, item.stop)
            if x >= y:
                return None
            return Interval(x, y)
        raise ValueError("Only slice supported!")

    def __contains__(self, item):
        print(item)


class GInterval(Interval):
    def __init__(self, x, y, chrom, name, strand='+', **kwargs):
        super(GInterval, self).__init__(x, y)
        self.chrom = chrom
        self.name = name
        self.strand = strand
        self._strand_backup = self.strand
        self.attribute = kwargs

    def __getitem__(self, item):
        if isinstance(item, slice):
            x = max(self.x, item.start)
            y = min(self.y, item.stop)
            if x >= y:
                return None
            return GInterval(x=x, y=y, chrom=self.chrom, name=self.name, strand=self.strand, **self.attribute)
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

    def __getitem__(self, item):
        if isinstance(item, slice):
            x = max(self.x, item.start)
            y = min(self.y, item.stop)
            if x >= y:
                return None
            blocks = list(filter(lambda z: z is not None, [block[x:y] for block in self.blocks]))
            if len(blocks) == 0:
                return None
            x = min([block.x for block in blocks])
            y = max([block.y for block in blocks])
            return GMultiInterval(x=x, y=y, chrom=self.chrom, name=self.name, strand=self.strand, blocks=blocks,
                                  **self.attribute)
        return super(GMultiInterval, self).__getitem__(item)

    def __len__(self):
        length = 0
        for i in range(0, len(self._blocks), 2):
            length += self._blocks[i]
        return length

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

    def __getitem__(self, item):
        if isinstance(item, slice):
            x = max(self.x, item.start)
            y = min(self.y, item.stop)
            if x >= y:
                return None
            blocks = list(filter(lambda z: z is not None, [block[x:y] for block in self.blocks]))
            if len(blocks) == 0:
                return None
            thick = None
            if self.thick is not None:
                thick = self.thick[x:y]
            x = min([block.x for block in blocks])
            y = max([block.y for block in blocks])
            return GComplexInterval(x=x, y=y, chrom=self.chrom, name=self.name, strand=self.strand, blocks=blocks,
                                    thick=thick, **self.attribute)
        return super(GMultiInterval, self).__getitem__(item)

    def __add__(self, other):
        blocks = list(self.blocks)
        blocks.extend(list(other.blocks))
        blocks = sorted(blocks, key=lambda z: z.x)
        last_block = None
        integrated = []
        for block in blocks:
            if last_block is None:
                last_block = block
                continue
            else:
                x = max(last_block.x, block.x)
                y = min(last_block.y, block.y)
                if y > x:
                    last_block = Interval(x, y)
                else:
                    integrated.append(last_block)
                    last_block = block
        if last_block is not None:
            integrated.append(last_block)
            # last_block = None
        x = min([block.x for block in integrated])
        y = max([block.y for block in integrated])
        return GComplexInterval(x=x, y=y, chrom=self.chrom, name=self.name, strand=self.strand, blocks=integrated,
                                **self.attribute)

    def overlap(self, other):
        blocks1 = list(other.blocks)
        blocks2 = list(self.blocks)
        i1 = 0
        i2 = 0
        mi1 = len(blocks1)
        mi2 = len(blocks2)
        while i1 < mi1 and i2 < mi2:
            b1 = blocks1[i1]
            b2 = blocks2[i2]
            if b1.y <= b2.x:
                i1 += 1
            elif b1.x >= b2.y:
                i2 += 1
            else:
                return True
        return False

    def contain(self, other):
        blocks1 = list(other.blocks)
        blocks2 = list(self.blocks)
        i2 = 0
        mi2 = len(blocks2)
        for b1 in blocks1:
            while True:
                if i2 >= mi2:
                    return False
                b2 = blocks2[i2]
                if b1.x >= b2.y:
                    i2 += 1
                    continue
                elif not (b1.x >= b2.x and b1.y <= b2.y):
                    return False
        return True

    def coincide(self, other):
        pass

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
'''


class IntervalTest(unittest.TestCase):
    def test_GInterval(self):
        '''i = GInterval(10, 20)
        self.assertEqual(len(i), 10)
        self.assertEqual(i.x, 10)
        self.assertEqual(i.y, 20)

        gi = GInterval(10, 50, chrom='chr1', name='GInterval1', strand='-')

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

        intervals = [GInterval(10, 20), GInterval(30, 35), GInterval(50, 65)]
        gi = GInterval(chrom='chr2', name='GInterval2', strand='-', blocks=intervals)

        self.assertEqual(gi.start_index, 0)
        self.assertEqual(gi.end_index, 29)
        self.assertEqual(len(gi), 30)

        self.assertEqual(gi.index(10), 29)
        gi.set_strand_sensitivity(False)
        self.assertEqual(gi.index(15), 5)
        gi.set_strand_sensitivity(True)
        self.assertEqual(gi.index(15), 24)

        self.assertEqual(gi.position(29), 10)

        intervals = [GInterval(10, 20), GInterval(30, 35), GInterval(50, 65)]
        thick = GInterval(15, 60)
        gi = GInterval(chrom='chr2', name='GInterval2', strand='-', blocks=intervals, thick=thick, rgb='read')

        self.assertEqual(gi.thick_start_position, 59)
        self.assertEqual(gi.thick_end_position, 15)
        self.assertEqual(gi.thick_start_index, 5)
        self.assertEqual(gi.thick_end_index, 24)
        self.assertEqual(gi['red'], None)

        g1 = gi[10:31]
        self.assertEqual(len(g1), 11)

        g1 = GInterval(blocks=[GInterval(10, 20), GInterval(20, 22), GInterval(100, 200)])
        g2 = GInterval(25, 27, chrom='chr2', name='GInterval2', strand='-', rgb='read')
        g3 = g1 + g2
        self.assertEqual(len(g3), 114)'''

        g1 = GInterval(blocks=[GInterval(10, 20), GInterval(30, 50), GInterval(90, 100)])
        g2 = GInterval(blocks=[GInterval(40, 50), GInterval(90, 95)])
        self.assertTrue(g1.coincide(g2))

        g1 = GInterval(blocks=[GInterval(10, 20), GInterval(30, 50), GInterval(90, 100)])
        g2 = GInterval(blocks=[GInterval(40, 52), GInterval(90, 95)])
        self.assertFalse(g1.coincide(g2))


if __name__ == '__main__':
    unittest.main()
