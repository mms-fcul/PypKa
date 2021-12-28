import os
from pypka.config import Config


def checkDelPhiErrors(filename, mode=None):
    exceptions = []
    if mode == "readFiles":
        exceptions = ["part of system outside the box!", "has a net charge of"]
    elif mode == "runDelPhi":
        exceptions = ["part of system outside the box!"]

    exit_trigger = False
    errors = ""
    with open(filename) as f:
        for line in f:
            ignore_error = False
            if "WARNING".lower() in line.lower() or "ERROR".lower() in line.lower():
                for exception in exceptions:
                    if exception in line:
                        ignore_error = True
                if not ignore_error:
                    exit_trigger = True
                    errors += line
    if exit_trigger:
        raise Exception(
            "The following errors have been found on {0}: \n{1}".format(
                filename, errors
            )
        )
    if not Config.debug:
        if os.path.isfile(filename):
            os.remove(filename)
        focusing_log = "{}_focusing".format(filename)
        if os.path.isfile(focusing_log):
            os.remove(focusing_log)
