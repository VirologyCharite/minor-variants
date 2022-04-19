from collections import defaultdict, Counter

from dark.reads import Reads
from dark.sam import SAMFilter, samfile


def isMinorVariantPosition(bases, minDepth, minFrequency):
    """
    Identify positions with minority variants, given some criteria.

    @param bases: A C{collections.Counter} of number of bases.
    @param minDepth: The C{int} number of read coverage that needs to be
        present at a position for it to be considered a minor variant.
    @param minFrequency: A C{float} minimum frequency with which at least two
        nucleotides need to be present at a position for it to be considered
        variable.
    """
    baseCount = sum(list(bases.values()))

    if baseCount == 0:
        return False

    if baseCount >= minDepth:
        if len(bases.keys()) == 1:
            return False

        frequencies = defaultdict(int)

        for base in bases:
            if bases[base] / baseCount > minFrequency:
                frequencies['above'] += 1
            else:
                frequencies['below'] += 1

        if frequencies['above'] >= 2:
            return True
    else:
        return False


def getBaseFrequencies(bamFile, minBaseQuality=0, minMappingQuality=0):
    """
    Takes a bam file and returns a dictionary where the key maps to a position
    and the values map to a Counter with the number of each base at that
    position.

    @param bamFile: A C{str} filename of a bam file.
    @param minBaseQuality: Minimum base quality. Bases below the minimum
        quality will not be output.
    @param minMappingQuality: Only use reads above a minimum mapping quality.
    """
    reads = Reads().filter()

    samFilter = SAMFilter(bamFile, filterRead=reads.filterRead)

    referenceLengths = samFilter.referenceLengths()

    result = {}

    with samfile(bamFile) as sam:
        if samFilter.referenceIds:
            # No need to check if the given reference id is in referenceLengths
            # because the samFilter.referenceLengths call above caught that.
            referenceId = samFilter.referenceIds.pop()
        else:
            if len(referenceLengths) == 1:
                referenceId = list(referenceLengths)[0]
            else:
                print('SAM file %r contains %d references (%s). Only one '
                      'reference id can be analyzed at a time.' % (
                          bamFile, len(referenceLengths),
                          ', '.join(sorted(referenceLengths))))

        for i, column in enumerate(sam.pileup(
                                   reference=referenceId,
                                   min_base_quality=minBaseQuality,
                                   min_mapping_quality=minMappingQuality,
                                   ignore_overlap=False)):
            bases = Counter()
            for read in column.pileups:
                if (not read.is_del and not read.is_refskip):
                    base = read.alignment.query_sequence[read.query_position]
                    bases[base] += 1
                elif read.is_del:
                    bases['-'] += 1

            result[column.reference_pos] = bases

        for position in range(referenceLengths[referenceId]):
            if position not in result:
                result[position] = Counter({'A': 0, 'T': 0, 'G': 0, 'C': 0})

    return result
