# Minor variants

This repo contains code for minor variant analysis.


The code can take a `BAM` file of reads aligned against a reference and compute various statistics related to minor variant detection.


## Overview

To open a `BAM` file, and return a `MinorVariantInfo` object:

```
mvi = MinorVariantInfo(bamFile='data/complete-coverage-insertion-sorted.bam')
```

The `MinorVariantInfo` object has information about coverage per base (`coveragePerBase`), the frequency of the most common nucleotide at each position (`maxFreqPerBase`), as well as methods to compute the mean coverage, the richness, complexity, distance, and nucleotide diversity. It also has a save function to save the computed information as `json`, which can also be parsed using `MinorVariantInfo`, which will be much faster than parsing the `BAM` file.

`MinorVariantInfo` also takes arguments to specify the minimum base quality and the minimum mapping quality.

A `BAM` file can also be parsed and saved as json using the script `bin/generate-data.py`:

```
$ python bin/generate-data.py --h
usage: generate-data.py [-h] [--sequencingTech SEQUENCINGTECH] [--minBaseQuality MINBASEQUALITY] [--minMappingQuality MINMAPPINGQUALITY] bamFile

Take a bamFile and write MinorVariantInfo.countsPerBase json to stdout.

positional arguments:
  bamFile               A bam file to be analysed.

options:
  -h, --help            show this help message and exit
  --sequencingTech SEQUENCINGTECH
                        The sequencing technology used to create the reads in the bam file.
  --minBaseQuality MINBASEQUALITY
                        Minimum base quality. Bases below the minimum quality will not be output.
  --minMappingQuality MINMAPPINGQUALITY
                        Only use reads above a minimum mapping quality.
```