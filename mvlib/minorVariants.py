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
    @param minBaseQuality: Minimum base quality. Bases below the minimum
        quality will not be output.
    @param minMappingQuality: Only use reads above a minimum mapping quality.
    @param sequencingTech: The sequencing technology that was used to create
        the bam file. Should be one of 'Unknown', 'MiSeq', 'NextSeq',
        'NovaSeq'.
    @param referenceId: The Id of the reference sequence to use. In case a BAM
        file contains multiple references.
    """
    def __init__(self, bamFile=None, jsonFile=None, frequenciesDict=False,
                 minBaseQuality=None, minMappingQuality=None,
                 sequencingTech=None, referenceId=False):

        self.minBaseQuality = minBaseQuality if minBaseQuality else 0
        self.minMappingQuality = minMappingQuality if minMappingQuality else 0

        if frequenciesDict:
            self.countsPerBase = frequenciesDict
            self.name = None
            self.sequencingTech = sequencingTech

        elif jsonFile:
            with open(jsonFile, 'r') as fp:
                openJson = json.load(fp)
                self.countsPerBase = {int(k): v for k, v in
                                      openJson['countsPerBase'].items()}
                params = openJson['parameters']
                self.sequencingTech = params['sequencingTech']
                self.minBaseQuality = params['minBaseQuality']
                self.minMappingQuality = params['minMappingQuality']
            self.name = jsonFile.split('/')[-1].split('.')[0]
        elif bamFile:
            self.countsPerBase = getBaseFrequencies(
                bamFile, self.minBaseQuality, self.minMappingQuality,
                referenceId=referenceId)
            self.name = bamFile.split('/')[-1].split('.')[0]
            self.sequencingTech = sequencingTech
        else:
            raise ('At least one out of bamFile, jsonFile, or frequenciesDict '
                   'must be specified.')

        self.coveragePerBase = []
        for position in range(len(self.countsPerBase)):
            self.coveragePerBase.append(self.countsPerBase[position])
        self.coveragePerBase = [sum(self.countsPerBase[i].values()) for i in
                                range(len(self.countsPerBase))]

        self.length = len(self.coveragePerBase)

        self.maxFreqPerBase = {}
        for position in self.countsPerBase:
            bases = list(self.countsPerBase[position].values())
            if bases:
                try:
                    frequency = max(bases) / sum(bases)
                except ZeroDivisionError:
                    frequency = 0.0
            else:
                frequency = 0.0
            self.maxFreqPerBase[position] = frequency

        assert self.length == len(self.countsPerBase)

    def save(self, outFilename=False):
        """
        Save self.countsPerBase to a json file.

        @param outFilename: A C{str} filename of the json file where
            self.countsPerBase should be written to.
        """
        data = {}
        data['parameters'] = {
            'sequencingTech': self.sequencingTech,
            'minBaseQuality': self.minBaseQuality,
            'minMappingQuality': self.minMappingQuality,
        }
        data['countsPerBase'] = self.countsPerBase

        if outFilename:
            with open(outFilename, 'w') as fp:
                json.dump(data, fp)
        else:
            json.dump(data, sys.stdout, indent=4, sort_keys=True)

    def meanCoverage(self):
        """
        Return the mean coverage of the entire file.
        """
        return np.mean(self.coveragePerBase)

    def richness(self, minCoverage=50, minFrequency=0.03,
                 printFrequencies=False):
        """
        Calculate the richness.

        @param minCoverage: The C{int} number of read coverage that needs to be
            present at a position for it to be considered a minor variant.
        @param minFrequency: A C{float} minimum frequency with which at least
            two nucleotides need to be present at a position for it to be
            considered variable.
        @param printFrequencies: if C{True}, print the base frequencies for the
            position that is a minor variant.
        """
        richnessCount = 0
        for position in self.countsPerBase:
            if isMinorVariantPosition(self.countsPerBase[position],
                                      minCoverage, minFrequency):
                richnessCount += 1
                if printFrequencies:
                    print('\t', position, self.countsPerBase[position],
                          richnessCount)
        return richnessCount

    def complexity(self, minCoverage=50, minFrequency=0.03):
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
                totalBases = sum(self.countsPerBase[position].values())
                shannonEntropy = 0
                for base in self.countsPerBase[position]:
                    freq = self.countsPerBase[position][base] / totalBases
                    shannonEntropy += (-freq * np.log(freq))
                shannonEntropies.append(shannonEntropy)
        return np.mean(shannonEntropies)

    def distance(self, minCoverage=50, minFrequency=0.03):
        """
        Calculate the distance.

        @param minCoverage: The C{int} number of read coverage that needs to be
            present at a position for it to be considered a minor variant.
        @param minFrequency: A C{float} minimum frequency with which at least
            two nucleotides need to be present at a position for it to be
            considered variable.
        """
        distance = 0
        for p in self.countsPerBase:
            if isMinorVariantPosition(self.countsPerBase[p],
                                      minCoverage, minFrequency):
                totalBases = sum(list(self.countsPerBase[p].values()))
                frequencies = sorted([self.countsPerBase[p][b] / totalBases for
                                      b in self.countsPerBase[p]],
                                     reverse=True)
                distance += frequencies[1]

        return distance

    def allelArray(self, minCoverage):
        """
        Return an array as used by scikit-allel.
        """
        bToC = {'A': 0, 'C': 1, 'T': 2, 'G': 3}

        allelArray = []
        for position in range(len(self.countsPerBase)):
            alleles = [0, 0, 0, 0]
            if (not minCoverage or (minCoverage and sum(
                    self.countsPerBase[position].values()) >= minCoverage)):
                for base, count in self.countsPerBase[position].items():
                    try:
                        alleles[bToC[base]] += count
                    except KeyError:
                        continue
                allelArray.append(alleles)

        if allelArray:
            return allel.AlleleCountsArray(allelArray)

    def nucleotideDiversity(self, perPosition=False, offsets=False,
                            minCoverage=False):
        """
        Calculate the nucleotide diversity pi.

        @param perPosition: if C{True} calculate the nucleotide diversity per
            position and return a list of the nucleotide diversity at each
            position.
        @param offsets: if not C{False} a C{tuple} of 0-based start, stop
            offsets between which the nucleotide diversity should be computed.
        @param minCoverage: if not C{False}, the minimum number of reads
            covering a position for it to be included in the calculation.
        """
        assert perPosition != offsets, ("Don't use 'perPosition' and "
                                        "'offsets' at the same time.")

        alleleCountsArray = self.allelArray(minCoverage=minCoverage)

        if alleleCountsArray:
            if perPosition:
                result = []
                for position in range(self.length):
                    pi = allel.sequence_diversity(list(range(self.length)),
                                                  alleleCountsArray,
                                                  start=position,
                                                  stop=position)
                    result.append(pi)

            if offsets:
                result = allel.sequence_diversity(list(range(self.length)),
                                                  alleleCountsArray,
                                                  start=offsets[0],
                                                  stop=offsets[1])

            return result
