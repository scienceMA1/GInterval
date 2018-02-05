"""GInterval object represent a interval on the genomic coordination.
It can consist many blocks (smaller interval which only contain x and y value).
Additional information, such as chromosome, name and strand, are also contained in the GInterval objects.

"""
import numpy as np
from operator import itemgetter


class GInterval(object):
    """GInterval is the basic element interval represent a range of genomic.

    Attributes:
        x (int): The start position (0-base, included) of the interval.
        y (int): The end position (0-base, excluded) of the interval.
        chrom (str): The chromosome on which the interval is located.
        name (str): The name of the interval.
        strand (str): The value is "+" or "-".
        forward (boolean): True when the strand is "+".
        reverse (boolean): True when the strand is "-".

    Examples:
        Create a GInterval instance before any operation:

        g = GInterval(10, 20)

        g = GInterval(blocks=[10, 20, 30, 40])



    """

    def __init__(self, x=None, y=None, blocks=None, thick=None, **kwargs):
        self._x = x
        self._y = y
        self._annotations = kwargs

        self._blocks = None
        if blocks is None:
            if self._x is not None and self._y is not None:
                self._blocks = [self._x, self._y]
            else:
                raise ValueError("The x and y or blocks must be provided!")
        else:
            if isinstance(blocks, list) or isinstance(blocks, np.ndarray):
                self._blocks = []
                for x0, y0 in blocks:
                    if not x0 < y0:
                        raise ValueError("The block.x value must be smaller than block.y value.")
                    self._blocks.append(x0)
                    self._blocks.append(y0)
                lp = None
                for p in self._blocks:
                    if lp is not None and p < lp:
                        raise ValueError("The position values of blocks are not sorted by oder.")
                    lp = p
                if len(self._blocks) < 2:
                    raise ValueError("Too little blocks position!")

                if self._x is None or self._y is None:
                    self._x = self._blocks[0]
                    self._y = self._blocks[-1]
                else:
                    if self._x != self._blocks[0] or self._y != self._blocks[-1]:
                        raise ValueError("The x, y position of GInterval is not consistent with blocks range!")
            else:
                raise ValueError("Invalid blocks values!")

        if self._x >= self._y:
            raise ValueError("The x position of GInterval must smaller than y position.")
        elif self._x < 0:
            raise ValueError("The position of GInterval can not be smaller than 0!")

        self._thick = None
        if thick is not None:
            tx, ty = thick
            tx0 = max(tx, self._x)
            ty0 = min(ty, self._y)
            if ty0 > tx0:
                flag_x = True
                flag_y = True
                for bx, by in self.gaps:
                    if flag_x and bx <= tx0 < by:
                        tx0 = by
                        flag_x = False
                    if flag_y and bx <= ty0 < by:
                        ty0 = bx
                        flag_y = False
                if ty0 > tx0:
                    self._thick = [tx0, ty0]

    @property
    def x(self):
        """

        :return: The x position (0-base included) of the interval instance.
        """
        return self._x

    @property
    def y(self):
        """

        :return: The y position (0-base excluded) of the interval instance.
        """
        return self._y

    @property
    def strand(self):
        """

        :return: The strand information of GInterval instance "+" or "-". Default "+"
        """
        return self._annotations.get('strand', '+')

    @property
    def forward(self):
        """

        :return: True when the strand equal to "+"
        """
        return self.strand == '+'

    @property
    def reverse(self):
        """

        :return: True when the strand equal to "-".
        """
        return self.strand == '-'

    @property
    def thick(self):
        return self._thick

    @property
    def blocks(self):
        for i in range(0, len(self._blocks), 2):
            yield (self._blocks[i], self._blocks[i + 1])

    @property
    def block_count(self):
        return len(self._blocks) // 2

    @property
    def gaps(self):
        if self._blocks is not None:
            for i in range(1, len(self._blocks) - 2, 2):
                yield (self._blocks[i], self._blocks[i + 1])

    def set_strand_sensitivity(self, sensitivity=True):
        self._annotations.setdefault('strand_backup', self._annotations.get('strand'))
        if sensitivity:
            self._annotations['strand'] = self._annotations.get('strand_backup')
        else:
            self._annotations['strand'] = '+'

    def __getitem__(self, item):
        """
        instance[1]
        instance[2:8]
        :param item:
        :return:
        """
        if isinstance(item, slice):
            if item.start >= item.stop:
                return None

            x = max(self._x, item.start)
            y = min(self._y, item.stop)

            if x >= y:
                return None

            blocks = []
            for x1, y1 in self.blocks:
                x2 = max(x1, x)
                y2 = min(y1, y)
                if y2 > x2:
                    blocks.append([x2, y2])
            if len(blocks) == 0:
                return None

            thick = None
            if self.thick is not None:
                tx, ty = self.thick
                tx1 = max(tx, x)
                ty1 = min(ty, y)
                if ty1 > tx1:
                    thick = (tx1, ty1)

            return GInterval(blocks=blocks, thick=thick, **self._annotations)

        elif isinstance(item, int):
            return self.index(item)
        else:
            raise TypeError("The type of item is not supported!")

    def __add__(self, other):
        """

        The strand of other will be ignored.

        :param other:
        :return:
        """
        integrated = []
        blocks = list(self.blocks) + list(other.blocks)
        blocks = sorted(blocks, key=itemgetter(0))
        x0 = None
        y0 = None
        for x, y in blocks:
            if x0 is None:
                x0, y0 = x, y
            else:
                x1 = max(x0, x)
                y1 = min(y0, y)
                if x1 < y1:
                    x0 = min(x0, x)
                    y0 = max(y0, y)
                else:
                    integrated.append([x0, y0])
                    x0, y0 = x, y

        if x0 is not None:
            integrated.append([x0, y0])

        for x, y in integrated:
            print(x, y)

        thick = None
        thick1 = self._thick
        thick2 = other.thick
        if thick1 is not None or thick2 is not None:
            if thick1 is None:
                thick = thick2
            elif thick2 is None:
                thick = thick1
            else:
                x1, y1 = thick1
                x2, y2 = thick2
                thick = [min(x1, x2), max(y1, y2)]

        return GInterval(blocks=integrated, thick=thick, **self._annotations)

    def __len__(self):
        """
        In generally, the base number is the total base number of all blocks.
        :return: The number of base consist the GInterval.
        """
        return sum((y - x for x, y in self.blocks))

    def __getattr__(self, item):
        return self._annotations.get(item)

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
        return self._thick is not None

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


if __name__ == '__main__':
    gi1 = GInterval(blocks=[[10, 20], [30, 40], [50, 60]], thick=(25, 55))
    gi2 = GInterval(blocks=[[10, 20], [30, 40], [45, 65]], thick=(25, 55))
    # gi3 = gi1 + gi2
    # print(len(gi3))
