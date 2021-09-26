import os
import subprocess
from pypka.config import Config

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def rinse_pdb(pdb_to_clean, pdb_cleaned, ff, ffout, logfile="LOG_pdb2pqr", hopt=True):
    try:
        # TODO: Port pdb2pqr to py3 and import it as a module
        cmd = (
            "python2 {0}/pdb2pqr/pdb2pqr.py {1} {2} "
            "--ff {3} --ffout {4} --drop-water -v --chain {6} > {5} 2>&1 ".format(
                SCRIPT_DIR,
                pdb_to_clean,
                pdb_cleaned,
                ff,
                ffout,
                logfile,
                "" if hopt else "--noopt",
            )
        )
        subprocess.run(
            cmd,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        raise Exception(
            "pdb2pqr did not run successfully\nMessage: {}".format(
                e.stderr.decode("ascii")
            )
        )


def add_tautomers(sites_addHtaut, ff_family, outputpqr):

    logfile = "LOG_addHtaut"
    try:
        # TODO rewrite addHtaut as python module
        cmd = "{}/addHtaut cleaned.pqr {} {} > {}".format(
            SCRIPT_DIR,
            ff_family,
            sites_addHtaut,
            outputpqr,
            logfile,
        )
        subprocess.run(
            cmd,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        raise Exception(
            "addHtaut did not run successfully\nMessage: {}".format(
                e.stderr.decode("ascii")
            )
        )