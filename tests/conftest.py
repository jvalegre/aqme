#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest

try:
    import DBGEN
    BASEPATH = os.path.join(DBGEN.__path__[0])
except ImportError:
    here = os.path.dirname(os.path.abspath(__file__))
    BASEPATH = os.path.normpath(os.path.join(here, '..', 'pyCONFORT'))

def datapath(path):
    return os.path.join(BASEPATH, 'examples', path)
