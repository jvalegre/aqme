#!/usr/bin/env python

###############################################.
#          __main__ file for the code         #
###############################################.

from __future__ import absolute_import

import sys
from aqme import aqme

# If we are running from a wheel, add the wheel to sys.path
# This allows the usage python pip-*.whl/pip install pip-*.whl

if __package__ != 'aqme':
    print('aqme is not installed! Use: pip install aqme (anywhere, using a terminal) or python setup.py install (from the downloaded /aqme/aqme folder).')

if __name__ == '__main__':
    aqme.main()
    sys.exit()
