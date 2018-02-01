from GInterval.IO import ParserFactory
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


class TranscriptAnalysis(object):
    def __init__(self, raw_data=None, fmt=None, sfile=None):

        self.summary = None

        self.raw_data = raw_data
        self.fmt = fmt
        self.sfile = sfile
        self.header = ['Name', '5UTR', 'CDS', '3UTR']

        if self.raw_data is not None:
            self.__load_from_raw_data()
        if self.sfile is not None:
            self.__load_from_summary_file()

    def __load_from_raw_data(self):
        names = []
        lengths1 = []
        lengths2 = []
        lengths3 = []

        for record in ParserFactory(self.raw_data, self.fmt).parser:
            name = record.name
            length1 = None
            length2 = None
            length3 = None

            if record.thick is not None:
                index0 = record.start_index
                index1 = record.thick_start_index
                index2 = record.thick_end_index
                index3 = record.end_index
                length1 = index1 - index0
                length2 = index2 - index1 + 1
                length3 = index3 - index2

            names.append(name)
            lengths1.append(length1)
            lengths2.append(length2)
            lengths3.append(length3)

        self.summary = pd.DataFrame(index=names, data={
            self.header[1]: lengths1,
            self.header[2]: lengths2,
            self.header[3]: lengths3
        })

        self.summary.index.name = self.header[0]

    def __load_from_summary_file(self):
        self.summary = pd.read_csv(self.sfile, delimiter='\t', index_col=0)

    @property
    def total_count(self):
        return len(self.summary)

    @property
    def thick_count(self):
        return sum(self.summary['CDS'].notnull())

    @property
    def smooth_count(self):
        return sum(self.summary['CDS'].isnull())

    @property
    def average_5UTR(self):
        return np.mean(self.summary['5UTR'])

    @property
    def average_CDS(self):
        return np.mean(self.summary['CDS'])

    @property
    def average_3UTR(self):
        return np.mean(self.summary['3UTR'])

    def __str__(self):
        lines = list()
        total_count = self.total_count
        thick_count = self.thick_count
        smooth_count = self.smooth_count
        lines.append('Total input record: %d' % total_count)
        lines.append('Record with thick interval: %d (%.2f%%)' %
                     (thick_count, thick_count * 100 / total_count))
        lines.append('Record without thick interval: %d (%.2f%%)' %
                     (smooth_count, smooth_count * 100 / total_count))
        lines.append('Average 5UTR length: %.2f' % self.average_5UTR)
        lines.append('Average CDS length: %.2f' % self.average_CDS)
        lines.append('Average 3UTR length: %.2f' % self.average_3UTR)
        return "\n".join(lines)

    def plot(self, output=None):
        plt.figure(figsize=(6, 3))

        data = self.summary[self.summary['CDS'].notnull()]
        plt.title('The length distribution of transcript components')
        plt.boxplot([data['5UTR'].values, data['CDS'].values, data['3UTR'].values], showfliers=False)
        plt.xticks(range(1, 4), self.header[1:])
        plt.xlabel('Transcript Components')
        plt.ylabel('Length Distribution (nt)')

        plt.tight_layout()
        if output is None:
            plt.show()
        else:
            plt.savefig(output, dpi=300)
        plt.close()

    def save(self, path):
        self.summary.to_csv(path, sep="\t", na_rep='NA')

    @classmethod
    def get_summary(cls, infile, fmt=None, summary=None):
        if summary is not None:
            if os.path.exists(summary):
                return TranscriptAnalysis(sfile=summary)
        ta = TranscriptAnalysis(infile, fmt=fmt)
        if summary is not None:
            ta.save(summary)
        return ta


if __name__ == '__main__':
    pass
