from collections import Counter
import os
import numpy as np

from mvlib.common import NTCOLORS

import matplotlib
if not os.environ.get('DISPLAY'):
    # Use non-interactive Agg backend
    matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plotCoverage(minorVariantInfo):
    """
    Plot the coverage of a MinorVariantInfo instance.
    """
    fig, ax = plt.subplot(1, 1, figsize=(15, 5))

    ax.plot(list(range(len(minorVariantInfo.coveragePerBase))),
            minorVariantInfo.coveragePerBase, '-')

    ax.set_title(minorVariantInfo.name)
    ax.set_xlabel('Position')
    ax.set_ylabel('Coverage')


def plotFrequencies(minorVariantInfo, positions, outFilename=None, ax=None,
                    outFormat='png', title=None, plotFurinAAs=True):
    """
    Make a single bar chart showing minor variants at specified positions.
    Note that the graph will produce have the 1-based positions on the x-axis.

    @param minorVariantInfo: A C{MinorVariantInfo} instance.
    @param positions: a C{list} of 1-based positions that should be plotted.
    @param ax: If not C{None}, use this as the subplot for plotting the base
        frequencies.
    @param outFilename: If not C{None}, a C{str} filename where the plot is
        saved to.
    @param outFormat: If not C{None}, a C{str} of the file format of the
        figure.
    @param title: If not C{None}, a C{str} that is added to the filename as
        title.
    @param plotFurinAAs: If C{True} plot the furin cleavage site amino acids
        at the bottom of the plot.
    """
    ind = np.arange(len(positions))
    width = 0.5

    furinSeq = 'GGTATATGCGCTAGTTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGCT'

    if not ax:
        fig, ax = plt.subplots(1, 1, figsize=(16, 3.75))

    # aFreqs = []
    # tFreqs = []
    # cFreqs = []
    # gFreqs = []

    coverage = []

    for i, position in enumerate(sorted(positions)):
        baseInfo = Counter(minorVariantInfo.countsPerBase[position - 1])
        total = sum(list(baseInfo.values()))
        coverage.append(total)
        aFreqs = []
        tFreqs = []
        cFreqs = []
        gFreqs = []
        originalBase = furinSeq[i]
        aAlpha = 1.0
        tAlpha = 1.0
        cAlpha = 1.0
        gAlpha = 1.0
        if originalBase == 'A':
            aAlpha = 0.3
        if originalBase == 'T':
            tAlpha = 0.3
        if originalBase == 'C':
            cAlpha = 0.3
        if originalBase == 'G':
            gAlpha = 0.3
        try:
            aFreqs.append(baseInfo['A'] / total)
        except ZeroDivisionError:
            aFreqs.append(0)
        try:
            tFreqs.append(baseInfo['T'] / total)
        except ZeroDivisionError:
            tFreqs.append(0)
        try:
            cFreqs.append(baseInfo['C'] / total)
        except ZeroDivisionError:
            cFreqs.append(0)
        try:
            gFreqs.append(baseInfo['G'] / total)
        except ZeroDivisionError:
            gFreqs.append(0)
        ax.bar(i, aFreqs, width, facecolor=NTCOLORS['A'],
               edgecolor=NTCOLORS['A'], alpha=aAlpha)
        ax.bar(i, tFreqs, width, bottom=aFreqs,
               facecolor=NTCOLORS['T'], edgecolor=NTCOLORS['T'], alpha=tAlpha)
        ax.bar(i, cFreqs, width,
               bottom=[i + j for i, j in zip(aFreqs, tFreqs)],
               facecolor=NTCOLORS['C'], edgecolor=NTCOLORS['C'], alpha=cAlpha)
        ax.bar(i, gFreqs, width,
               bottom=[i + j + k for i, j, k in zip(aFreqs, tFreqs, cFreqs)],
               facecolor=NTCOLORS['G'], edgecolor=NTCOLORS['G'], alpha=gAlpha)

    # ax.bar(ind, aFreqs, width, facecolor=NTCOLORS['A'],
    #        edgecolor=NTCOLORS['A'])
    # ax.bar(ind, tFreqs, width, bottom=aFreqs,
    #        facecolor=NTCOLORS['T'], edgecolor=NTCOLORS['T'])
    # ax.bar(ind, cFreqs, width,
    #        bottom=[i + j for i, j in zip(aFreqs, tFreqs)],
    #        facecolor=NTCOLORS['C'], edgecolor=NTCOLORS['C'])
    # ax.bar(ind, gFreqs, width,
    #        bottom=[i + j + k for i, j, k in zip(aFreqs, tFreqs, cFreqs)],
    #        facecolor=NTCOLORS['G'], edgecolor=NTCOLORS['G'])

    ax.hlines(y=0.5, xmin=-0.5, xmax=len(positions), linestyle=':')

    ax.set_xlim(-0.5, len(positions) - 0.5)

    ax.set_ylim(0.0, 1.0)

    titleAddition = '' if not title else ', %s' % title
    ax.set_title('%s%s' % (minorVariantInfo.name, titleAddition), fontsize=16)

    ax.xaxis.set_ticks(
        [i + 0.2 for i in range(len(positions))])

    xLabels = []
    for i, pos in enumerate(positions):
        xLabels.append('%d (%d)' % (pos, coverage[i]))
    ax.set_xticklabels(xLabels, fontsize=13, rotation=90)

    if plotFurinAAs:
        secax = ax.secondary_xaxis('bottom')
        secax.set_xlim(-0.5, len(positions) - 0.5)
        secax.tick_params(axis='both', which='major', pad=100, bottom=False)
        secax.xaxis.set_ticks([1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37,
                               40, 43, 46, 49, 52, 55, 58, 61])
        secax.set_xticklabels(['G', 'I', 'C', 'A', 'S', 'Y', 'Q', 'T', 'Q',
                               'T', 'N', 'S', 'P', 'R/W', 'R', 'A', 'R', 'S',
                               'V', 'A'], fontsize=13)

    if outFilename:
        plt.savefig(outFilename, format=outFormat, bbox_inches='tight')
