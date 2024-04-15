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
        coveragePerBase = [
            2, 8, 19, 30, 37, 38, 46, 48, 50, 54, 56, 56, 65, 72, 73, 76, 107,
            136, 143, 150, 161, 176, 177, 185, 196, 196, 196, 199, 207, 207,
            208, 208, 209, 210, 211, 218, 219, 224, 228, 238, 241, 241, 241,
            246, 246, 246, 251, 252, 257, 257, 260, 264, 264, 264, 264, 265,
            266, 264, 265, 265, 265, 265, 265, 265, 264, 262, 263, 263, 264,
            264, 265, 264, 264, 269, 270, 273, 271, 272, 266, 258, 258, 272,
            267, 276, 272, 273, 271, 275, 268, 267, 266, 263, 232, 205, 192,
            192, 190, 175, 174, 171
        ]
        self.assertEqual(coveragePerBase, mvi.coveragePerBase)

    def testCoveragePerBasePartialCoverage(self):
        """
        coveragePerBase must be correct when the reference is only partially
        covered.
        """
        bamFile = join(DATADIR, 'partial-coverage-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        coveragePerBase = [
            5, 5, 5, 5, 5, 5, 5, 5, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
        ]
        self.assertEqual(coveragePerBase, mvi.coveragePerBase)

    def testCoveragePerBaseDeletion(self):
        """
        coveragePerBase must be correct when the reads have a deletion relative
        to the reference.
        """
        bamFile = join(DATADIR, 'complete-coverage-deletion-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        coveragePerBase = [
            2, 8, 19, 30, 37, 38, 46, 48, 48, 54, 56, 56, 65, 72, 73, 76, 107,
            136, 143, 150, 161, 176, 177, 185, 196, 196, 196, 199, 207, 207,
            208, 208, 209, 210, 211, 218, 219, 224, 228, 238, 241, 241, 241,
            246, 246, 246, 251, 252, 257, 257, 260, 264, 264, 264, 264, 265,
            266, 264, 265, 265, 265, 265, 265, 265, 264, 262, 263, 263, 264,
            264, 265, 264, 264, 269, 270, 273, 271, 272, 266, 258, 258, 272,
            267, 276, 272, 273, 271, 275, 268, 267, 266, 263, 232, 205, 192,
            192, 190, 175, 174, 171
        ]
        self.assertEqual(coveragePerBase, mvi.coveragePerBase)

    def testCoveragePerBaseInsertion(self):
        """
        coveragePerBase must be correct when the reads have a insertion
        relative to the reference.
        NOTE THAT THE INSERTED POSITION ISN'T PART OF WHAT IS RETURNED!
        """
        bamFile = join(DATADIR, 'complete-coverage-insertion-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        coveragePerBase = [
            2, 8, 19, 30, 37, 38, 46, 48, 50, 54, 56, 56, 65, 72, 73, 76, 107,
            136, 143, 150, 161, 176, 177, 185, 196, 196, 196, 199, 207, 207,
            208, 208, 209, 210, 211, 218, 219, 224, 228, 238, 241, 241, 241,
            246, 246, 246, 251, 252, 257, 257, 260, 264, 264, 264, 264, 265,
            266, 264, 265, 265, 265, 265, 265, 265, 264, 262, 263, 263, 264,
            264, 265, 264, 264, 265, 265, 265, 262, 257, 247, 238, 227, 224,
            218, 217, 211, 211, 209, 213, 206, 205, 12, 10, 10, 10, 10, 10,
            10, 10, 10, 10
        ]
        self.assertEqual(coveragePerBase, mvi.coveragePerBase)

    def testMaxFreqPerBaseCompleteCoverage(self):
        """
        maxFreqPerBase must be correct when the reference is completely
        covered.
        """
        bamFile = join(DATADIR, 'complete-coverage-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        maxFreqPerBase = {
            0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0, 6: 1.0, 7: 1.0,
            8: 1.0, 9: 1.0, 10: 1.0, 11: 1.0, 12: 1.0, 13: 1.0, 14: 1.0,
            15: 1.0, 16: 1.0, 17: 1.0, 18: 1.0, 19: 0.9933333333333333,
            20: 1.0, 21: 1.0, 22: 1.0, 23: 1.0, 24: 0.9948979591836735,
            25: 1.0, 26: 1.0, 27: 1.0, 28: 1.0, 29: 1.0,
            30: 0.9951923076923077, 31: 0.9951923076923077, 32: 1.0, 33: 1.0,
            34: 1.0, 35: 1.0, 36: 1.0, 37: 1.0, 38: 1.0, 39: 1.0, 40: 1.0,
            41: 1.0, 42: 0.995850622406639, 43: 1.0, 44: 1.0, 45: 1.0,
            46: 1.0, 47: 1.0, 48: 1.0, 49: 1.0, 50: 1.0, 51: 1.0, 52: 1.0,
            53: 0.9962121212121212, 54: 1.0, 55: 0.9962264150943396,
            56: 1.0, 57: 1.0, 58: 1.0, 59: 0.9962264150943396, 60: 1.0,
            61: 1.0, 62: 1.0, 63: 0.9849056603773585, 64: 1.0,
            65: 0.9961832061068703, 66: 0.9961977186311787, 67: 1.0, 68: 1.0,
            69: 0.9962121212121212, 70: 1.0, 71: 1.0, 72: 0.9924242424242424,
            73: 0.9962825278810409, 74: 1.0, 75: 1.0, 76: 0.996309963099631,
            77: 0.9926470588235294, 78: 0.9962406015037594,
            79: 0.9961240310077519, 80: 1.0, 81: 1.0, 82: 1.0, 83: 1.0,
            84: 1.0, 85: 1.0, 86: 1.0, 87: 0.9963636363636363, 88: 1.0,
            89: 1.0, 90: 1.0, 91: 0.9961977186311787, 92: 1.0,
            93: 0.9951219512195122, 94: 1.0, 95: 0.9947916666666666, 96: 1.0,
            97: 0.9942857142857143, 98: 1.0, 99: 1.0
        }
        self.assertEqual(maxFreqPerBase, mvi.maxFreqPerBase)

    def testMaxFreqPerBasePartialCoverage(self):
        """
        maxFreqPerBase must be correct when the reference is only partially
        covered.
        """
        bamFile = join(DATADIR, 'partial-coverage-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        maxFreqPerBase = {
            0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0, 6: 1.0, 7: 1.0,
            8: 1.0, 9: 1.0, 10: 1.0, 11: 1.0, 12: 1.0, 13: 1.0, 14: 1.0,
            15: 1.0, 16: 1.0, 17: 1.0, 18: 1.0, 19: 1.0, 20: 1.0, 21: 1.0,
            22: 1.0, 23: 1.0, 24: 1.0, 25: 1.0, 26: 1.0, 27: 1.0, 28: 1.0,
            29: 1.0, 30: 1.0, 31: 1.0, 32: 1.0, 33: 1.0, 34: 1.0, 35: 1.0,
            36: 1.0, 37: 1.0, 38: 1.0, 39: 1.0, 40: 1.0, 108: 1.0, 109: 1.0,
            110: 1.0, 111: 1.0, 112: 1.0, 113: 1.0, 114: 1.0, 115: 1.0,
            116: 1.0, 117: 1.0, 118: 1.0, 119: 1.0, 120: 1.0, 121: 1.0,
            122: 1.0, 123: 1.0, 124: 1.0, 125: 1.0, 126: 1.0, 127: 1.0,
            128: 1.0, 129: 1.0, 41: 0.0, 42: 0.0, 43: 0.0, 44: 0.0, 45: 0.0,
            46: 0.0, 47: 0.0, 48: 0.0, 49: 0.0, 50: 0.0, 51: 0.0, 52: 0.0,
            53: 0.0, 54: 0.0, 55: 0.0, 56: 0.0, 57: 0.0, 58: 0.0, 59: 0.0,
            60: 0.0, 61: 0.0, 62: 0.0, 63: 0.0, 64: 0.0, 65: 0.0, 66: 0.0,
            67: 0.0, 68: 0.0, 69: 0.0, 70: 0.0, 71: 0.0, 72: 0.0, 73: 0.0,
            74: 0.0, 75: 0.0, 76: 0.0, 77: 0.0, 78: 0.0, 79: 0.0, 80: 0.0,
            81: 0.0, 82: 0.0, 83: 0.0, 84: 0.0, 85: 0.0, 86: 0.0, 87: 0.0,
            88: 0.0, 89: 0.0, 90: 0.0, 91: 0.0, 92: 0.0, 93: 0.0, 94: 0.0,
            95: 0.0, 96: 0.0, 97: 0.0, 98: 0.0, 99: 0.0, 100: 0.0, 101: 0.0,
            102: 0.0, 103: 0.0, 104: 0.0, 105: 0.0, 106: 0.0, 107: 0.0
        }
        self.assertEqual(maxFreqPerBase, mvi.maxFreqPerBase)

    def testMeanCoverageCompleteCoverage(self):
        """
        The correct mean coverage must be returned for the complete coverage
        bam file.
        """
        bamFile = join(DATADIR, 'complete-coverage-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        self.assertEqual(205.29, mvi.meanCoverage())

    def testMeanCoveragePartialCoverage(self):
        """
        The correct mean coverage must be returned for the partial coverage
        bam file.
        """
        bamFile = join(DATADIR, 'partial-coverage-sorted.bam')
        mvi = MinorVariantInfo(bamFile=bamFile)
        self.assertEqual(1.0692307692307692, mvi.meanCoverage())
