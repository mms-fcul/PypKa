import sys
import os
import pathlib

# sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "../")
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))

from pypka.cli import CLI

if __name__ == "__main__":
    CLI()
