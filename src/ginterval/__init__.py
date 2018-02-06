"""GInterval object represent a interval on the genomic coordination.
It can consist many blocks (smaller interval which only contain x and y value).
Additional information, such as chromosome, name and strand, are also contained in the GInterval objects.

"""
import numpy as np
from operator import itemgetter


class GInterval(object):
    """GInterval is the basic element interval represent a range
    of genomic.

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

        self._block_info = None
        if blocks is None:
            if self._x is not None and self._y is not None:
                self._block_info = [self._x, self._y]
            else:
                raise ValueError("The x and y or blocks must be provided!")
        else:
            if isinstance(blocks, list) or isinstance(blocks, np.ndarray):
                self._block_info = []
                for x0, y0 in blocks:
                    if not x0 < y0:
                        raise ValueError("The block.x value must be smaller than block.y value.")
                    self._block_info.append(x0)
                    self._block_info.append(y0)
                lp = None
                for p in self._block_info:
                    if lp is not None and p < lp:
                        print(kwargs["name"])
                        print(self._block_info)
                        print(p,lp)
                        raise ValueError("The position values of blocks are not sorted by oder.")
                    lp = p
                if len(self._block_info) < 2:
                    raise ValueError("Too little blocks position!")

                if self._x is None or self._y is None:
                    self._x = self._block_info[0]
                    self._y = self._block_info[-1]
                else:
                    if self._x != self._block_info[0] or self._y != self._block_info[-1]:
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
        for i in range(0, len(self._block_info), 2):
            yield (self._block_info[i], self._block_info[i + 1])

    @property
    def block_count(self):
        return len(self._block_info) // 2

    @property
    def gaps(self):
        if self._block_info is not None:
            for i in range(1, len(self._block_info) - 2, 2):
                yield (self._block_info[i], self._block_info[i + 1])

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
            return self.position(item)
        else:
            raise TypeError("The type of item is not supported!")

    def __add__(self, other):
        """Rewrite the operation of GInterval plus GInterval.
        In fact, it is a process that all blocks merge together.
        The overlapped blocks will be combined to one.
        The thick will be combine if exist.
        The annotation information will be the same as left one.
        The annotation information of `other` will be ignored.

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
        """
        The overlap of GInterval without blocks considered.
        :param other:
        :return:
        """
        return max(self._x, other.x) < min(self._y.other.y)

    def contain(self, other):
        """
        The bases of other are all in self block intervals.
        :param other:
        :return:
        """
        blocks1 = list(GInterval.combine_adjacent(self.blocks))
        i = 0
        i_max = len(blocks1)
        for x, y in other.blocks:
            find = False
            while True:
                x1, y1 = blocks1[i]
                if x1 <= x and y1 >= y:
                    find = True
                    break
                i += 1
                if i == i_max:
                    break
            if not find:
                return False

        return True

    def coincide(self, other, adjacent1=False, adjacent2=True):
        """
        The other block interval is a continuous subset of self.
        :param adjacent1:
        :param adjacent2:
        :param other:
        :return:
        """
        if not (self._x <= other.x and self._y >= other.y):
            return False

        blocks1 = self.blocks
        blocks2 = other.blocks
        if not adjacent1:
            blocks1 = GInterval.combine_adjacent(blocks1)
        if not adjacent2:
            blocks2 = GInterval.combine_adjacent(blocks2)

        x1, y1 = next(blocks1)
        head_allowed = True
        tail_allowed = True
        try:
            for x2, y2 in blocks2:
                while True:
                    if x2 >= y1:
                        x1, y1 = next(blocks1)
                        continue
                    elif x1 >= y2:
                        return False
                    else:
                        if x1 < x2 and y1 == y2:
                            if head_allowed:
                                head_allowed = False
                                break
                            else:
                                return False
                        elif x1 == x2 and y1 > y2:
                            if tail_allowed:
                                head_allowed = False
                                break
                            else:
                                return False
                        elif x1 < x2 and y1 > y2:
                            if head_allowed and tail_allowed:
                                head_allowed = False
                                tail_allowed = False
                                break
                            else:
                                return False
                        elif x1 == x2 and y1 == y2:
                            continue
                        else:
                            return False
        except StopIteration as err:
            return False

        return True

    def index(self, position):
        """
        The 0-base and strand-specific index correspond to specific position.
        :param position:
        :return:
        """
        if not self._x <= position < self._y:
            raise ValueError('Invalid position %d' % position)
        index = 0
        for x, y in self.blocks:
            if y <= position:
                index += (y - x)
                continue
            elif x <= position < y:
                index += position - x
                break
            else:
                raise ValueError('')
        if self.reverse:
            index = len(self) - index - 1
        return index

    def position(self, index):
        """
        The genomic position of specific strand-specific index.
        :param index:
        :return:
        """
        end_index = len(self) - 1
        min_index = -(end_index + 1)
        if not min_index <= index <= end_index:
            raise ValueError('Invalid index %d. The valid index range is [%d:%d]' %
                             (index, min_index, end_index))
        if index < 0:
            index -= min_index
        if self.reverse:
            index = end_index - index
        for x, y in self.blocks:
            index -= (y - x)
            if index < 0:
                return y + index
        raise Exception('Unknown Exception!')

    # Additional methods.

    @property
    def start_position(self):
        """

        :return: The strand-specific start genomic position (0-base, included).
        """
        if self.forward:
            return self._x
        else:
            return self._y - 1

    @property
    def end_position(self):
        """

        :return: The strand-specific end genomic position (0-base, included).
        """
        if self.forward:
            return self._y - 1
        else:
            return self._x

    @property
    def start_index(self):
        """

        :return: The strand-specific index of start genomic position (0-base).
        """
        return 0

    @property
    def end_index(self):
        """

        :return: The strand-specific index of end genomic position (0-base).
        """
        return len(self) - 1

    @property
    def center_index(self):
        """

        :return: The index of the strand-specific center genomic position. Blocks and Gaps are considered.
        """
        return len(self) // 2

    @property
    def center_position(self):
        """

        :return: The center genomic position of this GInterval instance. Blocks and Gaps are considered.
        """
        return self.position(self.center_index)

    # Thick operation.

    @property
    def thick_start_position(self):
        """

        :return: The genomic start position of thick interval. (0-base, included)
        """
        x, y = self._thick
        if self.forward:
            return x
        else:
            return y - 1

    @property
    def thick_end_position(self):
        """

        :return: The genomic end position of thick interval. (0-base, included)
        """
        x, y = self._thick
        if self.forward:
            return y - 1
        else:
            return x

    @property
    def thick_start_index(self):
        """

        :return: The index of the genomic start position of thick relative to GInterval instance. (0-base)
        """
        return self.index(self.thick_start_position)

    @property
    def thick_end_index(self):
        """

        :return: The index of the genomic end position of thick relative to GInterval instance. (0-base)
        """
        return self.index(self.thick_end_position)

    def has_thick(self):
        """

        :return: True when there is a thick interval in this GInterval instance.
        """
        return self._thick is not None

    # Output format

    def to_bed_format_string(self, ncol=12):
        cols = []
        if ncol >= 4:
            cols.extend([self.chrom, self._x, self._y, self.name])
        if ncol >= 6:
            score = self.score
            if score is None:
                score = '.'
            cols.extend([score, self.strand])
        if ncol >= 12:
            tx = self._x
            ty = self._y
            rgb = self.rgb
            if rgb is None:
                rgb = "."
            nblock = 0
            sizes = ""
            starts = ""
            if self.thick is not None:
                tx, ty = self.thick
            for bx, by in self.blocks:
                nblock += 1
                sizes += "%d," % (by - bx)
                starts += "%d," % (bx - self._x)
            cols.extend([tx, ty, rgb, nblock, sizes, starts])
        return "\t".join(map(str, cols))

    def to_gtf_format_string(self):
        chrom = self.chrom
        source = self.source
        feature = self.feature
        start = self._x + 1
        end = self._y
        score = self.score
        strand = self.strand
        frame = self.frame
        texts = []
        for k, v in self.group.items():
            texts.append("%s \"%s\";" % (k, v))
        groups = " ".join(texts)
        return "\t".join(map(str, [chrom, source, feature, start, end, score, strand, frame, groups]))

    @staticmethod
    def sorted_blocks(blocks):
        """
        Input a iterator of blocks and return the sorted iterator.
        Sorted by x position.
        :param blocks:
        :return:
        """
        return sorted(blocks, key=itemgetter(0))

    @staticmethod
    def combine_blocks(blocks1, blocks2, adjacent=True):
        """

        :param adjacent: If the adjacent blocks is allowed.
        :param blocks1:
        :param blocks2:
        :return:
        """
        blocks1 = list(blocks1)
        blocks2 = list(blocks2)
        blocks = blocks1 + blocks2
        blocks = GInterval.sorted_blocks(blocks)
        x0 = None
        y0 = None
        for x, y in blocks:
            if x0 is None:
                x0, y0 = x, y
            else:
                if x <= y0:
                    if adjacent and x == y0:
                        yield (x0, y0)
                        x0, y0 = x, y
                    else:
                        y0 = max(y0, y)
                else:
                    yield (x0, y0)
                    x0, y0 = x, y
        if x0 is not None:
            yield (x0, y0)

    @staticmethod
    def combine_adjacent(blocks):
        """
        Combine the adjacent blocks to one block.
        :param blocks: the list of sorted blocks.
        :return: 
        """
        x0 = None
        y0 = None
        for x, y in blocks:
            if x0 is None:
                x0, y0 = x, y
            else:
                if x == y0:
                    y0 = y
                else:
                    yield (x0, y0)
                    x0, y0 = x, y
        if x0 is not None:
            yield (x0, y0)


if __name__ == '__main__':
    gi1 = GInterval(20, 30)
    gi2 = GInterval(50, 60)

    pass
