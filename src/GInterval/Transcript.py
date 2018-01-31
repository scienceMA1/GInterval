from GInterval import GInterval


def has_thick(gi):
    return gi.thick is not None


class TranscriptAnalysis(object):
    def __init__(self):
        self.five_utr_lengths = []
        self.cds_lengths = []
        self.three_utr_lengths = []
        self.has_thick_count = 0
        self.total_count = 0
        self.without_thick_count = 0

    def process(self, gi):
        self.total_count += 1
        if has_thick(gi):
            if gi.name == 'uc009ski.1':
                for block in gi.blocks:
                    print(block)
            self.has_thick_count += 1
            length1 = gi.get_thin_length_1(strand_sensitive=True, block_sensitive=True)
            length2 = gi.get_thick_length(block_sensitive=True)
            length3 = gi.get_thin_length_2(strand_sensitive=True, block_sensitive=True)
            self.five_utr_lengths.append(length1)
            self.cds_lengths.append(length2)
            self.three_utr_lengths.append(length3)

        else:
            self.without_thick_count += 1

if __name__ == '__main__':
    from GInterval.IO import parser
    import numpy as np
    ta = TranscriptAnalysis()
    for record in parser("../data/mm9_knownGene.bed", "bed"):
        #print(GInterval.to_bed_format_string(record))
        ta.process(record)
    print(np.mean(ta.five_utr_lengths))
    print(np.mean(ta.cds_lengths))
    print(np.mean(ta.three_utr_lengths))
