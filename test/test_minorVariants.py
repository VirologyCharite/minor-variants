from os.path import join
from unittest import TestCase

from mvlib.minorVariants import MinorVariantInfo
from mvlib.common import DATADIR


class TestMinorVariantInfo(TestCase):
    """
    Tests for the MinorVariantInfo class.
    """
    def testCoveragePerBaseCompleteCoverage(self):
        """
        coveragePerBase must be correct when the reference is completely
        covered.
        """
        bamFile = join(DATADIR, 'complete-coverage-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        coveragePerBase = [2, 8, 19, 30, 37, 38, 46, 48, 50, 54, 56, 56, 65,
                           72, 73, 76, 107, 136, 143, 150, 161, 176, 177, 185,
                           196, 196, 196, 199, 207, 207, 208, 208, 209, 210,
                           211, 218, 219, 224, 228, 238, 241, 241, 241, 246,
                           246, 246, 251, 252, 257, 257, 260, 264, 264, 264,
                           264, 264, 266, 264, 265, 265, 265, 265, 265, 264,
                           264, 262, 263, 263, 264, 264, 265, 264, 263, 268,
                           270, 273, 271, 272, 266, 258, 258, 272, 267, 276,
                           272, 273, 271, 275, 268, 267, 266, 263, 232, 205,
                           192, 192, 190, 175, 174, 171]
        self.assertEqual(coveragePerBase, mvi.coveragePerBase)

    def testCoveragePerBasePartialCoverage(self):
        """
        coveragePerBase must be correct when the reference is only partially
        covered.
        """
        bamFile = join(DATADIR, 'partial-coverage-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        coveragePerBase = [5, 5, 5, 5, 5, 5, 5, 5, 4, 3, 3, 3, 3, 3, 3, 3, 3,
                           3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                           2, 2, 2, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        self.assertEqual(coveragePerBase, mvi.coveragePerBase)

    def testCoveragePerBaseDeletion(self):
        """
        coveragePerBase must be correct when the reads have a deletion relative
        to the reference.
        """
        bamFile = join(DATADIR, 'complete-coverage-deletion-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        coveragePerBase = [2, 8, 19, 30, 37, 38, 46, 48, 0, 54, 56, 56, 65, 72,
                           73, 76, 107, 136, 143, 150, 161, 176, 177, 185, 196,
                           196, 196, 199, 207, 207, 208, 208, 209, 210, 211,
                           218, 219, 224, 228, 238, 241, 241, 241, 246, 246,
                           246, 251, 252, 257, 257, 260, 264, 264, 264, 264,
                           264, 266, 264, 265, 265, 265, 265, 265, 264, 264,
                           262, 263, 263, 264, 264, 265, 264, 263, 268, 270,
                           273, 271, 272, 266, 258, 258, 272, 267, 276, 272,
                           273, 271, 275, 268, 267, 266, 263, 232, 205, 192,
                           192, 190, 175, 174, 171]
        self.assertEqual(coveragePerBase, mvi.coveragePerBase)

    def testCoveragePerBaseInsertion(self):
        """
        coveragePerBase must be correct when the reads have a insertion
        relative to the reference.
        NOTE THAT THE INSERTED POSITION ISN'T PART OF WHAT IS RETURNED!
        """
        bamFile = join(DATADIR, 'complete-coverage-insertion-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        coveragePerBase = [2, 8, 19, 30, 37, 38, 46, 48, 50, 54, 56, 56, 65,
                           72, 73, 76, 107, 136, 143, 150, 161, 176, 177, 185,
                           196, 196, 196, 199, 207, 207, 208, 208, 209, 210,
                           211, 218, 219, 224, 228, 238, 241, 241, 241, 246,
                           246, 246, 251, 252, 257, 257, 260, 264, 264, 264,
                           264, 264, 266, 264, 265, 265, 265, 265, 265, 264,
                           264, 262, 263, 263, 264, 264, 265, 264, 263, 264,
                           265, 265, 262, 257, 247, 238, 227, 224, 218, 217,
                           211, 211, 209, 213, 206, 205, 12, 10, 10, 10, 10,
                           10, 10, 10, 10, 10]
        self.assertEqual(coveragePerBase, mvi.coveragePerBase)

    def testMaxFreqPerBaseCompleteCoverage(self):
        """
        maxFreqPerBase must be correct when the reference is completely
        covered.
        """
        bamFile = join(DATADIR, 'complete-coverage-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        maxFreqPerBase = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          0.9933333333333333, 1.0, 1.0, 1.0, 1.0,
                          0.9948979591836735, 1.0, 1.0, 1.0, 1.0, 1.0,
                          0.9951923076923077, 0.9951923076923077, 1.0, 1.0,
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          0.995850622406639, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          1.0, 1.0, 1.0, 1.0, 0.9962121212121212, 1.0, 1.0,
                          1.0, 1.0, 1.0, 0.9962264150943396, 1.0, 1.0, 1.0,
                          0.9886363636363636, 1.0, 0.9961832061068703,
                          0.9961977186311787, 1.0, 1.0, 0.9962121212121212,
                          1.0, 1.0, 0.9961977186311787, 1.0, 1.0, 1.0,
                          0.996309963099631, 0.9926470588235294,
                          0.9962406015037594, 0.9961240310077519, 1.0, 1.0,
                          1.0, 1.0, 1.0, 1.0, 1.0, 0.9963636363636363, 1.0,
                          1.0, 1.0, 0.9961977186311787, 1.0,
                          0.9951219512195122, 1.0, 0.9947916666666666, 1.0,
                          0.9942857142857143, 1.0, 1.0]
        self.assertEqual(maxFreqPerBase, mvi.maxFreqPerBase)

    def testMaxFreqPerBasePartialCoverage(self):
        """
        maxFreqPerBase must be correct when the reference is only partially
        covered.
        """
        bamFile = join(DATADIR, 'partial-coverage-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        maxFreqPerBase = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        self.assertEqual(maxFreqPerBase, mvi.maxFreqPerBase)

    def testMeanCoverageCompleteCoverage(self):
        """
        The correct mean coverage must be returned for the complete coverage
        bam file.
        """
        bamFile = join(DATADIR, 'complete-coverage-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        self.assertEqual(205.25, mvi.meanCoverage())

    def testMeanCoveragePartialCoverage(self):
        """
        The correct mean coverage must be returned for the partial coverage
        bam file.
        """
        bamFile = join(DATADIR, 'partial-coverage-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        self.assertEqual(1.0692307692307692, mvi.meanCoverage())
