from os.path import dirname, join

import mvlib

TOPDIR = dirname(dirname(mvlib.__file__))
DATADIR = join(TOPDIR, 'data')

NTCOLORS = {
    'A': '#1f76b4',  # blue
    'T': '#ff7e0e',  # orange
    'G': 'gold',
    'C': 'firebrick'  # '#2ca02c', green
}
