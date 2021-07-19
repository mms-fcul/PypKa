__author__ = "Pedro B. P. S. Reis"
__email__ = "pdreis@fc.ul.pt"

import os, sys

#sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from ._version import __version__
from .main import Titration, getTitrableSites
