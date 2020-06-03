#!/usr/bin/env python

###############################################.
#          __main__ file for the code         #
###############################################.

from __future__ import absolute_import

import sys
from pyconfort import pyconfort

# If we are running from a wheel, add the wheel to sys.path
# This allows the usage python pip-*.whl/pip install pip-*.whl

if __package__ != 'pyconfort':
    print('pyCONFORT is not installed! Use: pip install pyconfort (anywhere, using a terminal) or python setup.py install (from the downloaded /pyCONFORT/pyconfort folder).')

if __name__ == '__main__':
    pyconfort.main()
    sys.exit()
