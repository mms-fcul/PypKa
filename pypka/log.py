import os
import sys

from pypka.config import Config


class Log:
    def __init__(self):
        self.f_log = "LOG"
        self.stdout = None
        self.stderr = None
        self.stdout_file = None
        self.stderr_file = None

    def report2log(self, info, stdout=False):
        with open(self.f_log, "a") as logfile:
            logfile.write(info + "\n")
        if stdout:
            print(info)

    @staticmethod
    def raise_required_param_error(parameter):
        raise IOError(
            'Required input parameter "{0}" ' "is not defined.".format(parameter)
        )

    @staticmethod
    def raise_input_param_error(parameter, complaint, explanation):
        raise ValueError(
            'Input parameter "{0}" is not {1}\n '
            "{2}".format(parameter, complaint, explanation)
        )

    def report_warning(self, info, stdout=False):
        warning = "warning: {0}".format(info)
        self.report2log(warning, stdout=stdout)


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
