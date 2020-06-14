#!/usr/bin/env python

import argparse

from mvlib.minorVariants import MinorVariantInfo


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=('Take a bamFile and write MinorVariantInfo.countsPerBase '
                     'json to stdout.'))

    parser.add_argument('bamFile', help='A bam file to be analysed.')

    parser.add_argument(
        '--sequencingTech', default=None,
        help='The sequencing technology used to create the reads in the bam '
             'file.')

    parser.add_argument(
        '--minBaseQuality', default=None, type=int,
        help='Minimum base quality. Bases below the minimum quality will not '
             'be output.')

    parser.add_argument(
        '--minMappingQuality', default=None, type=int,
        help='Only use reads above a minimum mapping quality.')

    args = parser.parse_args()

    mvi = MinorVariantInfo(bamFile=args.bamFile,
                           minBaseQuality=args.minBaseQuality,
                           minMappingQuality=args.minMappingQuality,
                           sequencingTech=args.sequencingTech)

    mvi.save()
