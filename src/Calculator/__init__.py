from GInterval.IO import ShiftLoader, ParserFactory
from GInterval import GInterval

def cal_shift_around_anchor():
    pass


def cal_shift_inside_reference(spath, rpath, opath, sfmt=None, rfmt=None, extend=0, strand_sensitive=True):
    number = 0
    with open(opath, "w+") as fw:
        sl = ShiftLoader(path=spath, fmt=sfmt)
        for ri in ParserFactory(path=rpath, fmt=rfmt).parser:
            if ri.thick is None:
                continue

            index0 = ri.start_index
            index1 = ri.thick_start_index
            index2 = ri.thick_end_index
            index3 = ri.end_index

            flen = index1
            clen = index2 - index1 + 1
            tlen = index3 - index2

            sfts = sl.get_shifts(ri, extend=extend)
            if sfts is None:
                continue

            for shift in sfts:
                if not strand_sensitive:
                    shift.strand = ri.strand
                if shift.strand != ri.strand:
                    continue
                try:
                    index = ri.index(shift.center_position)
                except ValueError as e:
                    continue

                # 1. Calculate the relative distance.
                if index0 <= index < index1:
                    relative = (index - index0) / flen
                elif index1 <= index <= index2:
                    relative = (index - index1) / clen + 1
                else:
                    relative = (index - index2 - 1) / tlen + 2

                # 2. Calculate the absolute distance.
                dis0 = index - index0
                dis1 = index - index1
                dis2 = index - index2
                dis3 = index - index3

                # 3. Save the result.
                fw.write("\t".join([
                    ri.to_bed_format_string(6),
                    shift.to_bed_format_string(6),
                    "\t".join(map(str, [relative, dis0, dis1, dis2, dis3]))
                ]) + "\n")

                number += 1

    return number

if __name__ == '__main__':
    '''cal_shift_inside_reference(spath='../data/GSM2300428_m6A_peak_region_HeLa_CY_hg19.sorted.bed',
                               rpath='../data/RFP_inputRNA_TE_Class1.sorted.bed',
                               opath='../data/RFP_inputRNA_TE_Class1_m6A.bed',
                               strand_sensitive=False)
    cal_shift_inside_reference(spath='../data/GSM2300428_m6A_peak_region_HeLa_CY_hg19.sorted.bed',
                               rpath='../data/RFP_inputRNA_TE_Class2.sorted.bed',
                               opath='../data/RFP_inputRNA_TE_Class2_m6A.bed',
                               strand_sensitive=False)
    cal_shift_inside_reference(spath='../data/GSM2300428_m6A_peak_region_HeLa_CY_hg19.sorted.bed',
                               rpath='../data/RFP_inputRNA_TE_Class3.sorted.bed', 
                               opath='../data/RFP_inputRNA_TE_Class3_m6A.bed',
                               strand_sensitive=False)'''

    '''cal_shift_inside_reference(spath='../data/GSM2300428_m6A_peak_region_HeLa_CY_hg19.sorted.bed',
                               rpath='../data/P-Body-Enrichment_Class1.sorted.bed',
                               opath='../data/P-Body-Enrichment_Class1_m6A.bed',
                               strand_sensitive=False)
    cal_shift_inside_reference(spath='../data/GSM2300428_m6A_peak_region_HeLa_CY_hg19.sorted.bed',
                               rpath='../data/P-Body-Enrichment_Class2.sorted.bed',
                               opath='../data/P-Body-Enrichment_Class2_m6A.bed',
                               strand_sensitive=False)
    cal_shift_inside_reference(spath='../data/GSM2300428_m6A_peak_region_HeLa_CY_hg19.sorted.bed',
                               rpath='../data/P-Body-Enrichment_Class3.sorted.bed',
                               opath='../data/P-Body-Enrichment_Class3_m6A.bed',
                               strand_sensitive=False)'''

    '''cal_shift_inside_reference(spath='../data/GSM2300428_m6A_peak_region_HeLa_CY_hg19.sorted.bed',
                               rpath='../data/hg19_knownGene.bed',
                               opath='../data/hg19_knownGene_m6A.bed',
                               strand_sensitive=False)'''


    '''cal_shift_inside_reference(spath='../data/random.bed',
                               rpath='../data/hg19_knownGene.bed',
                               opath='../data/hg19_knownGene_random.bed')'''

    '''cal_shift_inside_reference(spath='../data/random.bed',
                               rpath='../data/RFP_inputRNA_TE_Class1.sorted.bed',
                               opath='../data/RFP_inputRNA_TE_Class1_random.bed')
    cal_shift_inside_reference(spath='../data/random.bed',
                               rpath='../data/RFP_inputRNA_TE_Class2.sorted.bed',
                               opath='../data/RFP_inputRNA_TE_Class2_random.bed')
    cal_shift_inside_reference(spath='../data/random.bed',
                               rpath='../data/RFP_inputRNA_TE_Class3.sorted.bed',
                               opath='../data/RFP_inputRNA_TE_Class3_random.bed')'''

class Calculator(object):
    pass
    ''''@classmethod
    def cal_shift_around_anchor(cls, apath, spath, opath, width=30):
        """

        Calculate the distance of shift intervals around anchor intervals.

        :param apath: The path of file which contain anchor intervals
        :param spath: The path of file which contain shift intervals.
        :param opath: The path of output file which stores the results.
        :param width: The valid distance for shift intervals to anchor intervals.
        :return: The number of valid results.
        """
        assert apath is not None
        assert spath is not None
        assert opath is not None

        number = 0

        with open(opath, "w+") as fw:
            sl = ShiftLoader(spath)
            for anchor in IntervalFactory.from_bed_file(apath):
                sfts = sl.get_shifts(anchor, width)
                if sfts is None:
                    continue
                for shift in sfts:
                    if shift.forward != anchor.forward:
                        continue
                    dis = anchor.cal_distance(shift)
                    if abs(dis) > width:
                        continue
                    fw.write("\t".join([
                        anchor.to_bed_format_string(6),
                        shift.to_bed_format_string(6),
                        str(dis)
                    ]) + "\n")
                    number += 1

        return number

    @classmethod
    def cal_shift_on_reference(cls, spath, rpath, opath):
        """

        Calculate the distribution of shift intervals inside reference intervals.

        :param spath: The path of file which contain shift intervals.
        :param rpath: The path of file which contain reference intervals.
        :param opath: The path of output file which stores the results.
        :return: The number of valid results.
        """
        assert spath is not None
        assert rpath is not None
        assert opath is not None

        number = 0

        with open(opath, "w+") as fw:
            sl = ShiftLoader(spath)
            for reference in IntervalFactory.from_bed_file(rpath):
                if not reference.has_thick_interval:
                    continue

                five_length = reference.utr5_length
                cds_length = reference.cds_length
                three_length = reference.utr3_length

                index0 = 0
                index1 = reference.thick_start_index
                index2 = reference.thick_end_index
                index3 = reference.length - 1

                sfts = sl.get_shifts(reference, 0)
                if sfts is None:
                    continue

                for shift in sfts:
                    if shift.forward != reference.forward:
                        continue

                    index = reference.get_index_by_position(shift.center_position)
                    if index is None:
                        continue

                    # 1. Calculate the relative distance.
                    if index0 <= index < index1:
                        relative = (index - index0) / five_length
                    elif index1 <= index <= index2:
                        relative = (index - index1) / cds_length + 1
                    else:
                        relative = (index - index2 - 1) / three_length + 2

                    # 2. Calculate the absolute distance.
                    dis0 = index - index0
                    dis1 = index - index1
                    dis2 = index - index2
                    dis3 = index - index3

                    # 3. Save the result.
                    fw.write("\t".join([
                        reference.to_bed_format_string(6),
                        shift.to_bed_format_string(6),
                        "\t".join(map(str, [relative, dis0, dis1, dis2, dis3]))
                    ]) + "\n")

                    number += 1

        return number

    @classmethod
    def cal_shift_around_anchor_on_reference(cls, apath, spath, rpath, opath, width=30):
        """

        Calculate the relative distances of shift intervals around anchor intervals on reference.

        :param apath: The path of file which contain the anchor intervals.
        :param spath: The path of file which contain the shift intervals.
        :param rpath: The path of file which contain the reference intervals.
        :param width: The valid distance for shift intervals to anchor intervals.
        :param opath: The path of output file which stores the results (distances).
        :return: The number of valid results.
        """
        assert apath is not None
        assert spath is not None
        assert rpath is not None
        assert opath is not None

        number = 0

        with open(opath, "w+") as fw:
            shift_loader = ShiftLoader(spath)
            anchor_loader = ShiftLoader(apath)
            for reference in IntervalFactory.from_bed_file(rpath):
                if not reference.has_thick_interval:
                    continue

                shifts = shift_loader.get_shifts(reference, 0)
                anchors = anchor_loader.get_shifts(reference, 0)
                if shifts is None or anchors is None:
                    continue

                for shift in shifts:
                    if shift.forward != reference.forward:
                        continue
                    sindex = reference.get_index_by_position(shift.center_position)
                    if sindex is None:
                        continue
                    for anchor in anchors:
                        if anchor.forward != reference.forward:
                            continue
                        aindex = reference.get_index_by_position(anchor.center_position)
                        if aindex is None:
                            continue
                        distance = anchor.cal_distance(shift)
                        if abs(distance) > width:
                            continue
                        fw.write("\t".join([
                            anchor.to_bed_format_string(6),
                            shift.to_bed_format_string(6),
                            str(distance)
                        ]) + "\n")
                        number += 1

        return number'''
