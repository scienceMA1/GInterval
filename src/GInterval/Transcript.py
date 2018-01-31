from GInterval.IO import ParserFactory
import numpy as np


class TranscriptAnalysis(object):
    def __init__(self, raw_data=None, fmt=None, processed_data=None):
        self.data = None

        self.__length1 = []
        self.__length2 = []
        self.__length3 = []
        self.__names = []
        self.header = ['Name', '5UTR', 'CDS', '3UTR']

        self.count_total = 0
        self.count_thick = 0
        self.count_without_thick = 0

        if raw_data is not None:
            for record in ParserFactory(raw_data, fmt).parser:
                self.__process(record)
            self.__after_process()
        if processed_data is not None:
            pass

    def __process(self, record):
        self.count_total += 1
        name = record.name
        length1 = None
        length2 = None
        length3 = None

        if record.thick is None:
            self.count_without_thick += 1

        else:
            self.count_thick += 1
            length1 = record.get_thin_length_1()
            length2 = record.get_thick_length()
            length3 = record.get_thin_length_2()
        self.__names.append(name)
        self.__length1.append(length1)
        self.__length2.append(length2)
        self.__length3.append(length3)

    def __after_process(self):
        import pandas as pd
        self.data = pd.DataFrame(index=self.__names, data={
            self.header[1]: self.__length1,
            self.header[2]: self.__length2,
            self.header[3]: self.__length3
        })
        self.data.index.name = self.header[0]
        self.__names = None
        self.__length1 = None
        self.__length2 = None
        self.__length3 = None

    def __str__(self):
        lines = []
        lines.append('Total input record: %d' % self.count_total)
        lines.append('Record with thick interval: %d (%.2f%%)' %
                     (self.count_thick, self.count_thick * 100 / self.count_total))
        lines.append('Record without thick interval: %d (%.2f%%)' %
                     (self.count_without_thick, self.count_without_thick * 100 / self.count_total))
        lines.append('Average 5UTR length: %.2f' % np.mean(self.data['5UTR']))
        lines.append('Average CDS length: %.2f' % np.mean(self.data['CDS']))
        lines.append('Average 3UTR length: %.2f' % np.mean(self.data['3UTR']))
        return "\n".join(lines)

    def plot(self, output=None):
        import matplotlib.pyplot as plt
        plt.figure(figsize=(4, 3))
        data = self.data[self.data[self.header[1]].notnull()]
        plt.boxplot([data['5UTR'].values, data['CDS'].values, data['3UTR'].values], showfliers=False)
        plt.xticks(range(1, 4), self.header[1:])
        plt.tight_layout()
        if output is None:
            plt.show()
        else:
            plt.savefig(output, dpi=300)
        plt.close()


if __name__ == '__main__':
    ta = TranscriptAnalysis("../data/mm9_knownGene.bed")
    print(ta)
    ta.plot("../data/transcript.png")
