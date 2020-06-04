from os.path import dirname, join

import mvlib

TOPDIR = dirname(dirname(mvlib.__file__))
DATADIR = join(TOPDIR, 'data')
