#!/usr/bin/env python

import argparse

from mvlib.minorVariants import MinorVariantInfo


parser = argparse.ArgumentParser(
    description=('Take a bamFile and write MinorVariantInfo.countsPerBase '
                 'json to stdout.'))

parser.add_argument('bamFile', help='A bam file to be analysed.')

args = parser.parse_args()

mvi = MinorVariantInfo(bamFile=args.bamFile)

mvi.save()
