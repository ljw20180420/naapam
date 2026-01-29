import subprocess
import sys
from importlib import resources


def parse():
    with resources.as_file(resources.files(".") / "parse.sh") as pf:
        subprocess.run(args=["bash", pf.as_posix(), sys.argv[1]])


def align():
    with resources.as_file(resources.files(".") / "align.sh") as pf:
        subprocess.run(args=["bash", pf.as_posix(), sys.argv[1]])
