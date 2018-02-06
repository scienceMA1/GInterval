from ginterval.parser import Parser
from ginterval import GInterval
from operator import itemgetter

class GTFParser(Parser):
    @staticmethod
    def _process_features(features):
        cdss = []
        exons = []
        for feature in features:
            feature_name = feature.feature
            if feature_name == 'exon':
                exons.append((feature.x, feature.y))
            elif feature_name == 'CDS':
                cdss.append((feature.x, feature.y))
        thick = None
        if len(cdss) > 0:
            tx = min((b[0] for b in cdss))
            ty = max((b[1] for b in cdss))
            thick = (tx, ty)
        feature = features[0]
        exons = GInterval.sorted_blocks(exons)
        exons = list(GInterval.combine_blocks(exons, exons))

        return GInterval(blocks=exons, thick=thick, chrom=feature.chrom,
                         name=feature.group["transcript_id"], strand=feature.strand, score=feature.score)

    @staticmethod
    def __parse_one_row(row):
        cols = row.split("\t")
        if len(cols) < 9:
            raise TypeError("The file is not standard GTF format file.")

        parameters = dict()
        parameters["chrom"] = cols[0]
        parameters["source"] = cols[1]
        parameters["feature"] = cols[2]
        parameters["x"] = int(cols[3]) - 1
        parameters["y"] = int(cols[4])
        parameters["score"] = cols[5]
        parameters["strand"] = cols[6]
        parameters["frame"] = cols[7]

        group = dict()
        a = filter(lambda x: x != "", map(lambda y: y.strip(), cols[8].split(";")))
        for b in a:
            b = b.strip()
            k, v = b.split()
            v = v.strip("\"")
            group[k] = v

        parameters["group"] = group

        return GInterval(**parameters)

    def __iter__(self):
        with open(self.path) as f:
            features = list()
            for line in f:
                if line.startswith('#'):
                    continue

                gi = GTFParser.__parse_one_row(line.strip("\n"))

                if len(features) == 0:
                    features.append(gi)
                else:
                    if gi.group["transcript_id"] == features[-1].group["transcript_id"]:
                        features.append(gi)
                    else:
                        yield self._process_features(features)
                        features = [gi]

            if len(features) > 0:
                yield self._process_features(features)

    @staticmethod
    def sort_gtf_by_transcript_id(infile, outfile):
        def sort_by_attribute(obj):
            return obj.group["transcript_id"]

        with open(infile) as f:
            with open(outfile, "w+") as fw:
                for gi in sorted((GTFParser.__parse_one_row(line.strip("\n")) for line in f), key=sort_by_attribute):
                    fw.write(gi.to_gtf_format_string() + "\n")

    @staticmethod
    def convert_gtf_to_bed(infile, outfile):
        with open(outfile, "w+") as fw:
            for gi in GTFParser(infile):
                fw.write(gi.to_bed_format_string() + "\n")


if __name__ == '__main__':
    '''parser = GTFParser("../../../test/data/test.gtf")
    with open("../../../test/data/test.gtf.bed", "w+") as fw:
        for gi in parser:
            fw.write(gi.to_bed_format_string() + "\n")'''

    file1 = "../../../data/genecode.v27lift37.sorted.gtf"
    file2 = "../../../data/genecode.v27lift37.sorted_by_transcript_id.gtf"
    file3 = "../../../data/genecode.v27lift37.sorted_by_transcript_id.bed"
    # GTFParser.sort_gtf_by_transcript_id(file1, file2)
    GTFParser.convert_gtf_to_bed(file2, file3)
