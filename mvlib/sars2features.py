from mvlib.common import (SARS2GENOME, SARS2OFFSETS, CODONSTOAA, YFVGENOME,
                          YFVOFFSETS, WNVGENOME, WNVOFFSETS)

VI = {
    'SARS2': {
        'g': SARS2GENOME,
        'o': SARS2OFFSETS,
    },
    'WNV': {
        'g': WNVGENOME,
        'o': WNVOFFSETS,
    },
    'YFV': {
        'g': YFVGENOME,
        'o': YFVOFFSETS,
    },
}


def getGene(position, virus):
    """
    Get the gene a particular position is in.

    @param position: The C{int} 0-based offset off a position in the SARS2
        genome.
    """
    if virus == 'SARS2' and position in (13467, 13468):
        return 'ORF1ab'
    else:
        for gene in VI[virus]['o']:
            if position in range(VI[virus]['o'][gene][0],
                                 VI[virus]['o'][gene][1]):
                return gene
        return 'UTR'


def getCodonAtPosition(position, nt, virus):
    """
    Get the codon that a particular position is in.
    """
    gene = getGene(position, virus)

    if gene == 'UTR':
        return 'non-coding'
    elif gene == 'ORF1ab':
        raise "Hoping this never happens"
        # orf1aSequence = SARS2GENOME.sequence[265:13468]
        # orf1aCodons = [orf1aSequence[i:i + 3]
        #                for i in range(0, len(orf1aSequence), 3)]
        # orf1aPos = position - 265
        # orf1bSequence = SARS2GENOME.sequence[13467:21555]
        # corf1bCodons = [orf1bSequence[i:i + 3]
        #                 for i in range(0, len(orf1bSequence), 3)]
        # orf1bPos = position - 13467
    else:
        ntSequence = VI[virus]['g'].sequence[VI[virus]['o'][gene][0]:
                                             VI[virus]['o'][gene][1]]
        codons = [ntSequence[i:i + 3] for i in range(0, len(ntSequence), 3)]
        ntPos = position - VI[virus]['o'][gene][0]
        newNtSequence = ntSequence[:ntPos] + nt + ntSequence[ntPos + 1:]
        newCodons = [newNtSequence[i:i + 3] for i in
                     range(0, len(newNtSequence), 3)]
        codonPosition = int(ntPos / 3)
        return codons[codonPosition], newCodons[codonPosition]


def isNS(position, nt, virus):
    """
    Return true if a change at a position is non-synonymous.
    """
    assert virus in {'SARS2', 'WNV', 'YFV'}, ('"virus" must be one of '
                                              '"SARS2", "WNV", "YFV".')
    try:
        oldCodon, newCodon = getCodonAtPosition(position, nt, virus)
    except ValueError:
        return False

    if CODONSTOAA[oldCodon] != CODONSTOAA[newCodon]:
        # print('\t', CODONSTOAA[oldCodon], CODONSTOAA[newCodon])
        return True
    else:
        return False
