import json
import numpy as np
import sys

import allel

from mvlib.functions import getBaseFrequencies, isMinorVariantPosition


class MinorVariantInfo():
    """
    Hold information about the minor variants in one sample.

    @param bamFile: The C{str} filename of the bam file with the reads aligned
        to a reference.
    @param name: A C{str} name of this instance.
    @param jsonFile: If not C{None}, a C{str} filename of a json file that
        contains the result of getBaseFrequencies().
    @param frequenciesDict: A C{dict} as returned by getBaseFrequencies().
    """
    def __init__(self, bamFile=None, jsonFile=None, frequenciesDict=False):
        if frequenciesDict:
            self.countsPerBase = frequenciesDict
            self.name = None
        elif jsonFile:
            with open(jsonFile, 'r') as fp:
                xx = json.load(fp)
                self.countsPerBase = {int(k): v for k, v in xx.items()}
            self.name = jsonFile.split('/')[-1].split('.')[0]
        elif bamFile:
            self.countsPerBase = getBaseFrequencies(bamFile)
            self.name = bamFile.split('/')[-1].split('.')[0]
        else:
            raise ('At least one out of bamFile, jsonFile, or frequenciesDict '
                   'must be specified.')

        self.coveragePerBase = [sum(self.countsPerBase[i].values()) for i in
                                range(max(self.countsPerBase)) if
                                len(self.countsPerBase[i].values()) > 0]

        self.length = len(self.coveragePerBase)

        self.maxFreqPerBase = []
        for i in range(max(self.countsPerBase) + 1):
            bases = list(self.countsPerBase[i].values())
            if bases:
                try:
                    frequency = max(bases) / sum(bases)
                except ZeroDivisionError:
                    frequency = 0.0
            else:
                frequency = 0.0
            self.maxFreqPerBase.append(frequency)

    def save(self, outFilename=False):
        """
        Save self.countsPerBase to a json file.

        @param outFilename: A C{str} filename of the json file where
            self.countsPerBase should be written to.
        """
        if outFilename:
            with open(outFilename, 'w') as fp:
                json.dump(self.countsPerBase, fp)
        else:
            json.dump(self.countsPerBase, sys.stdout)

    def meanCoverage(self):
        """
        Return the mean coverage of the entire file.
        """
        return np.mean(self.coveragePerBase)

    def richness(self, minCoverage=50, minFrequency=0.3):
        """
        Calculate the richness.

        @param minCoverage: The C{int} number of read coverage that needs to be
            present at a position for it to be considered a minor variant.
        @param minFrequency: A C{float} minimum frequency with which at least
            two nucleotides need to be present at a position for it to be
            considered variable.
        """
        richness = 0
        for position in self.countsPerBase:
            if isMinorVariantPosition(self.countsPerBase[position],
                                      minCoverage, minFrequency):
                richness += 1
        return richness

    def complexity(self, minCoverage=50, minFrequency=0.3):
        """
        Calculate the complexity.

        @param minCoverage: The C{int} number of read coverage that needs to be
            present at a position for it to be considered a minor variant.
        @param minFrequency: A C{float} minimum frequency with which at least
            two nucleotides need to be present at a position for it to be
            considered variable.
        """
        shannonEntropies = []
        for position in self.countsPerBase:
            if isMinorVariantPosition(self.countsPerBase[position],
                                      minCoverage, minFrequency):
                freq = max(self.maxFreqPerBase[position])
                shannonEntropy = (-(freq * np.log(freq)) +
                                  (1 - freq) * np.log(1 - freq)) / (np.log(2))
                shannonEntropies.append(shannonEntropy)
        return np.mean(shannonEntropies)

    def distance(self, minCoverage=50, minFrequency=0.3):
        """
        Calculate the distance.

        @param minCoverage: The C{int} number of read coverage that needs to be
            present at a position for it to be considered a minor variant.
        @param minFrequency: A C{float} minimum frequency with which at least
            two nucleotides need to be present at a position for it to be
            considered variable.
        """
        distance = 0
        for position in self.countsPerBase:
            if isMinorVariantPosition(self.countsPerBase[position],
                                      minCoverage, minFrequency):
                distance += max(self.maxFreqPerBase[position])

        return distance

    def allelArray(self):
        """
        Return an array as used by scikit-allel.
        """
        P2C = {'A': 0, 'C': 1, 'T': 2, 'G': 3}
        C2P = {0: 'A', 1: 'C', 2: 'T', 3: 'G'}
        alleleCounts = []
        align_pos = 1
        for position in ref_record:
            alleles = [0, 0, 0, 0]
            if not np.isnan(align_index[align_pos]):
                for seq in seqs:
                    if seq[align_pos-1] in P2C:
                        alleles[P2C[seq[align_pos-1]]] += 1


        allele_counts.append(alleles)
    align_pos += 1
allel1 = allel.AlleleCountsArray(allele_counts)

    def nucleotideDiversity(self, perPosition=False, overall=False,
                            offsets=False):
        """
        Calculate the nucleotide diversity pi.

        @param perPosition: if C{True} calculate the nucleotide diversity per
            position and return a list of the nucleotide diversity at each
            position.
        @param overall: if C{True}, calculate the nucleotide diversity across
            the entire sequence.
        @param offsets: if not C{False} a C{tuple} of 0-based start, stop
            offsets between which the nucleotide diversity should be computed.
        """
        if perPosition:
            result = []
            for position in self.countsPerBase:
                # calculate dif.
                result.append(pi)

