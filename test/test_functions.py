from unittest import TestCase
from collections import Counter

from mvlib.functions import isMinorVariantPosition


class TestIsMinorVariantPosition(TestCase):
    """
    Tests for the isMinorVariantPosition function.
    """
    def testOneBaseAboveFreqCutoff(self):
        """
        If one out of four bases is above the function must return False.
        """
        bases = Counter(['A', 'C', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G',
                         'G', 'G', 'G', 'G', 'G', 'G'])
        self.assertFalse(isMinorVariantPosition(bases, 10, 0.2))

    def testTwoBasesAboveFreqCutoff(self):
        """
        If two out of four bases is above the function must return True.
        """
        bases = Counter(['A', 'C', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G',
                         'G', 'G', 'G', 'G', 'G', 'G', 'C', 'C', 'C', 'C', 'C',
                         'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'])
        self.assertTrue(isMinorVariantPosition(bases, 10, 0.2))

    def testCoverageDepthTooLow(self):
        """
        If two out of four bases is above the function must return True.
        """
        bases = Counter(['A', 'T', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G',
                         'G', 'G', 'G', 'G', 'G', 'G', 'C', 'C', 'C', 'C', 'C',
                         'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'])
        self.assertFalse(isMinorVariantPosition(bases, 33, 0.2))

    def testThreeBasesAboveFreqCutoff(self):
        """
        If one out of four bases is above the function must return False.
        """
        bases = Counter(['A', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T',
                         'T', 'T', 'T', 'T', 'T', 'G', 'G', 'G', 'G', 'G', 'G',
                         'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'C', 'C',
                         'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
                         'C', 'C'])
        self.assertTrue(isMinorVariantPosition(bases, 10, 0.2))
