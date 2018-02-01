import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def draw_shift_on_reference(data, output=None, tas=None,
                            bsize=20, bcount=100, title=None, ticks=None, labels=None, ralen=None):
    """

    Drawing the distribution of shift intervals inside reference intervals.

    :param output: The path of output file.
    :param data:  The path of files contains the information of shift distribution.
                   It should be tab-delimited and only the last five columns is available.
    :param bsize: The bin-size (number of base) of absolute distribution.
    :param bcount: The bin-count for both relative distribution and absolute distribution.
    :param title: The title of figure.
    :param ticks: The ticks for gene body classifiers.
    :param labels: The labels of samples.
    :param ralen: The average length of gene body.
    :return: None
    """

    """
    1. Checking and init the input parameters. 
    """

    if not isinstance(data, list):
        data = [data]
    if title is None:
        title = "The distribution of shift intervals inside reference intervals"
    if ticks is None:
        ticks = ["5Cap", "Start Codon", "Stop Codon", "PolyA"]
    if labels is None:
        labels = list(map(str, range(len(data))))
    if not isinstance(labels, list):
        labels = [labels]
    if ralen is None:
        ralen = (238.93, 1530.46, 1122.19)
    width = bsize * bcount
    offset = bsize // 2

    """
    2. Init the canvas for drawing.
    """

    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(2, 6)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1:3])
    ax4 = fig.add_subplot(gs[1, 3:5])
    ax5 = fig.add_subplot(gs[1, 5])
    ls = {"color": "gray", "ls": "--", "lw": 0.5}

    xs1 = np.array(list(range(0, 3 * bcount))) + 0.5
    ax1.set_title(title)
    ax1.set_xlim(0, 3 * bcount)
    ax1.set_ylabel("Normalized Number (Number/kb)")
    ax1.set_xticks([0, bcount, 2 * bcount, 3 * bcount])
    ax1.set_xticklabels(ticks)
    ax1.axvline(bcount, **ls)
    ax1.axvline(2 * bcount, **ls)

    xs2 = list(range(offset, width, bsize))
    ax2.set_ylabel("Number")
    ax2.set_xlim(0, width)
    ax2.set_xlabel("Distance to " + ticks[0])

    xs3 = list(range(-width + offset, width, bsize))
    ax3.axvline(0, **ls)
    ax3.set_xlim(-width, width)
    ax3.set_xlabel("Distance to " + ticks[1])

    xs4 = xs3
    ax4.axvline(0, **ls)
    ax4.set_xlim(-width, width)
    ax4.set_xlabel("Distance to " + ticks[2])

    xs5 = list(range(-width + offset, 0, bsize))
    ax5.set_xlim(-width, 0)
    ax5.set_xlabel("Distance to " + ticks[3])

    """
    3. Loading the data for drawing and plot them.
    """

    relative_max = 0
    number_max = 0

    for shift, label, ta in zip(data, labels, tas):
        relatives = np.zeros(3 * bcount)
        dis0s = np.zeros(bcount)
        dis1s = np.zeros(2 * bcount)
        dis2s = np.zeros(2 * bcount)
        dis3s = np.zeros(bcount)

        with open(shift) as f:
            for line in f:
                relative, dis0, dis1, dis2, dis3 = line.strip("\n").split("\t")[-5:]
                relatives[int(float(relative) * bcount)] += 1
                dis0, dis1, dis2, dis3 = map(int, [dis0, dis1, dis2, dis3])

                if dis0 < width:
                    dis0s[dis0 // bsize] += 1
                if -width <= dis1 < width:
                    dis1s[(dis1 + width) // bsize] += 1
                if -width <= dis2 < width:
                    dis2s[(dis2 + width) // bsize] += 1
                if -width < dis3 <= 0:
                    dis3s[(dis3 + width - 1) // bsize] += 1

        dis0s = dis0s / bsize
        dis1s = dis1s / bsize
        dis2s = dis2s / bsize
        dis3s = dis3s / bsize

        number_max = max(number_max, max(map(max, [dis0s, dis1s, dis2s, dis3s])))

        relatives = relatives.reshape((3, -1))

        relative0 = relatives[0] * 1000 / ta.average_5UTR / ta.thick_count
        relative1 = relatives[1] * 1000 / ta.average_CDS / ta.thick_count
        relative2 = relatives[2] * 1000 / ta.average_3UTR / ta.thick_count

        relatives = np.append(relative0, np.append(relative1, relative2))
        relative_max = max(relative_max, max(relatives))

        ax1.plot(xs1, relatives, label=label, lw=3, ls="-")
        ax2.plot(xs2, dis0s)
        ax3.plot(xs3, dis1s)
        ax4.plot(xs4, dis2s)
        ax5.plot(xs5, dis3s)

    """
    4. Adjusting location and save the figure.
    """

    ylim_max = number_max * 1.2
    ax1.set_ylim(0, relative_max * 1.2)
    ax2.set_ylim(0, ylim_max)
    ax3.set_ylim(0, ylim_max)
    ax4.set_ylim(0, ylim_max)
    ax5.set_ylim(0, ylim_max)

    if len(data) > 1:
        ax1.legend()

    fig.tight_layout()
    if output is None:
        plt.show()
    else:
        plt.savefig(output, dip=300)
    plt.close()

if __name__ == '__main__':
    from GInterval.Transcript import TranscriptAnalysis
    from Calculator import cal_shift_inside_reference

    if True:
        references = ['../data/RFP_inputRNA_TE_Class%d.sorted.bed' % i for i in range(1, 4)]
        shift = '../data/GSE49339_A-PARCLIP.sorted.bed'
        transcripts = [TranscriptAnalysis.get_summary(reference, fmt='BED',
                                                      summary='../data/RFP_inputRNA_TE_Class%d.txt' % (i + 1))
                       for i, reference in enumerate(references)]
        distributions = ["../data/RFP_inputRNA_TE_Class%d_YTHDF2.txt" % i for i in range(1, 4)]
        # for distribution, reference in zip(distributions, references):
        # cal_shift_inside_reference(shift, reference, distribution, strand_sensitive=False)
        draw_shift_on_reference(data=distributions, tas=transcripts, bcount=10,
                                output='../data/RFP_inputRNA_TE_YTHDF2.png')

    if False:
        references = ['../data/RFP_inputRNA_TE_Class%d.sorted.bed' % i for i in range(1, 4)]
        shift = '../data/GSM2300428_m6A_peak_region_HeLa_CY_hg19.sorted.bed'
        transcripts = [TranscriptAnalysis.get_summary(reference, fmt='BED',
                                                      summary='../data/RFP_inputRNA_TE_Class%d.txt' % (i + 1))
                       for i, reference in enumerate(references)]
        distributions = ["../data/RFP_inputRNA_TE_Class%d_m6A.txt" % i for i in range(1, 4)]
        # for distribution, reference in zip(distributions, references):
        #    cal_shift_inside_reference(shift, reference, distribution, strand_sensitive=False)
        draw_shift_on_reference(data=distributions, tas=transcripts, bcount=10,
                                output='../data/RFP_inputRNA_TE_m6A.png')

    if False:
        references = ['../data/P-Body-Enrichment_Class%d.sorted.bed' % i for i in range(1, 4)]
        shift = '../data/GSE49339_A-PARCLIP.sorted.bed'
        transcripts = [TranscriptAnalysis.get_summary(reference, fmt='BED',
                                                      summary='../data/P-Body-Enrichment_Class%d.txt' % (i + 1))
                       for i, reference in enumerate(references)]
        distributions = ["../data/P-Body-Enrichment_Class%d_YTHDF2.txt" % i for i in range(1, 4)]
        # for distribution, reference in zip(distributions, references):
        #    cal_shift_inside_reference(shift, reference, distribution, strand_sensitive=False)
        draw_shift_on_reference(data=distributions, tas=transcripts, bcount=10,
                                output='../data/P-Body-Enrichment_YTHDF2.png')

    if False:
        references = ['../data/P-Body-Enrichment_Class%d.sorted.bed' % i for i in range(1, 4)]
        shift = '../data/GSM2300428_m6A_peak_region_HeLa_CY_hg19.sorted.bed'
        transcripts = [TranscriptAnalysis.get_summary(reference, fmt='BED',
                                                      summary='../data/P-Body-Enrichment_Class%d.txt' % (i + 1))
                       for i, reference in enumerate(references)]
        distributions = ["../data/P-Body-Enrichment_Class%d_m6A.txt" % i for i in range(1, 4)]
        # for distribution, reference in zip(distributions, references):
        #    cal_shift_inside_reference(shift, reference, distribution, strand_sensitive=False)
        draw_shift_on_reference(data=distributions, tas=transcripts, bcount=10,
                                output='../data/P-Body-Enrichment_m6A.png')

    if False:
        shift = "../data/random.bed"
        reference = "../data/hg19_knownGene.bed"
        cpath = "../data/hg19_knownGene_random.bed"
        tpath = "../data/hg19_knownGene.transcript"
        ta = TranscriptAnalysis.get_summary(reference, summary=tpath)
        dat = cal_shift_inside_reference(spath=shift, rpath=reference, opath=cpath, strand_sensitive=False)
        draw_shift_on_reference(dat, tas=[ta])
