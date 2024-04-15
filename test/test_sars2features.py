from unittest import TestCase

from mvlib.sars2features import getGene, getCodonAtPosition, isNS


class TestSARS2FeaturesGetGene(TestCase):
    """
    Tests for the sars2features.getGene function.
    """
    def testGetGeneUTR(self):
        """
        Position in ORF1 must be returned correctly.
        """
        self.assertEqual('UTR', getGene(50, 'SARS2'))

    def testGetGeneORF1aStart(self):
        """
        Position in ORF1 must be returned correctly.
        """
        self.assertEqual('ORF1a', getGene(265, 'SARS2'))

    def testGetGeneORF1a(self):
        """
        Position in ORF1 must be returned correctly.
        """
        self.assertEqual('ORF1a', getGene(300, 'SARS2'))

    def testGetGeneORF1b(self):
        """
        Position in ORF1 must be returned correctly.
        """
        self.assertEqual('ORF1b', getGene(20000, 'SARS2'))

    def testGetGeneORF1bEnd(self):
        """
        Position in ORF1 must be returned correctly.
        """
        self.assertEqual('ORF1b', getGene(21554, 'SARS2'))

    def testGetGeneORF1ab(self):
        """
        Position in ORF1 must be returned correctly.
        """
        self.assertEqual('ORF1ab', getGene(13467, 'SARS2'))

    def testGetGeneWNVPoly(self):
        """
        Position in WNV polyprotein must be returned correctly.
        """
        self.assertEqual('Polyprotein', getGene(300, 'WNV'))

    def testGetGeneWNVUTR(self):
        """
        Position in WNV UTR must be returned correctly.
        """
        self.assertEqual('UTR', getGene(1, 'WNV'))

    def testGetGeneYFVPoly(self):
        """
        Position in YFV Polyprotein must be returned correctly.
        """
        self.assertEqual('Polyprotein', getGene(300, 'YFV'))

    def testGetGeneYFVUTR(self):
        """
        Position in YFV UTR must be returned correctly.
        """
        self.assertEqual('UTR', getGene(1, 'YFV'))


class TestSARS2FeaturesGetCodonAtPosition(TestCase):
    """
    Tests for the sars2features.getCodonAtPosition function.
    """
    def testCorrectCodon1(self):
        """
        The correct codon must be returned.
        """
        oldCodon, newCodon, position = getCodonAtPosition(265, 'A', 'SARS2')
        self.assertEqual('ATG', oldCodon)
        self.assertEqual('ATG', newCodon)

    def testCorrectCodon2(self):
        """
        The correct codon must be returned.
        """
        oldCodon, newCodon, position = getCodonAtPosition(266, 'T', 'SARS2')
        self.assertEqual('ATG', oldCodon)
        self.assertEqual('ATG', newCodon)

    def testCorrectCodon3(self):
        """
        The correct codon must be returned.
        """
        oldCodon, newCodon, position = getCodonAtPosition(267, 'G', 'SARS2')
        self.assertEqual('ATG', oldCodon)
        self.assertEqual('ATG', newCodon)

    def testCorrectCodon4(self):
        """
        The correct codon must be returned.
        """
        oldCodon, newCodon, position = getCodonAtPosition(268, 'G', 'SARS2')
        self.assertEqual('GAG', oldCodon)
        self.assertEqual('GAG', newCodon)

    def testCorrectCodonWNV(self):
        """
        The correct codon must be returned in WNV.
        """
        oldCodon, newCodon, position = getCodonAtPosition(96, 'A', 'WNV')
        self.assertEqual('ATG', oldCodon)
        self.assertEqual('ATG', newCodon)

    def testCorrectCodonYFV(self):
        """
        The correct codon must be returned.
        """
        oldCodon, newCodon, position = getCodonAtPosition(118, 'A', 'YFV')
        self.assertEqual('ATG', oldCodon)
        self.assertEqual('ATG', newCodon)


class TestSARS2FeaturesIsNS(TestCase):
    """
    Tests for the ars2features.isNS function.
    """
    def testIsNSTrue(self):
        """
        If a mutation is non-synonymous, return correctly
        """
        self.assertEqual(['start', 'L', 1], isNS(265, 'C', 'SARS2'))

    def testIsNSFalse(self):
        """
        If a mutation is synonymous, return correctly
        """
        self.assertEqual([False, False, 1], isNS(265, 'A', 'SARS2'))
