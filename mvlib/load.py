import glob

from mvlib.common import DATADIR
from mvlib.minorVariants import MinorVariantInfo


def load(minMeanCoverage=None):
    """
    Return a generator of MinorVariantInfo instances of all available json
    files in data, provided they pass the filtering.

    @param minMeanCoverage: If not C{None} a C{float} minimum mean coverage
        that an alignment needs to have to be returned.
    """
    files = glob.glob(DATADIR + '*.json')

    for file in files:
        mvi = MinorVariantInfo(jsonFile=file)
        if minMeanCoverage:
            if mvi.meanCoverage > minMeanCoverage:
                yield mvi
        else:
            yield mvi
