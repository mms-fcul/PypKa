import os
import sys
from config import Config

class Log:
    def __init__(self):
        self.f_log       = "LOG"
        self.stdout      = None
        self.stderr      = None
        self.stdout_file = None
        self.stderr_file = None

    def report2log(self, info, stdout=False):
        with open(self.f_log, 'a') as logfile:
            logfile.write(info + '\n')
        if stdout:
            print(info)

    def raise_required_param_error(self, parameter):
        raise IOError('Required input parameter "{0}" '
                      'is not defined.'.format(parameter))


    def raise_input_param_error(self, parameter, complaint, explanation):
        raise ValueError('Input parameter "{0}" is not {1}\n '
                      '{2}'.format(parameter, complaint, explanation))

    def report_warning(self, info, stdout=False):
        warning = 'warning: {0}'.format(info)
        self.report2log(warning, stdout=stdout)


def redirectOutput(mode, outputname):
    if mode == 'start':
        config.stdout = os.dup(sys.stdout.fileno())
        config.stdout_file = open(outputname, 'w')
        os.dup2(config.stdout_file.fileno(), sys.stdout.fileno())
    elif mode == 'stop':
        if not config.stdout:
            raise Exception('Output redirection has not been started, '
                            'thus it can not be ended.')
        sys.stdout.flush()
        os.dup2(config.stdout, sys.stdout.fileno())
        config.stdout_file.close()


def redirectErr(mode, outputname):
    if mode == 'start':
        config.stderr = os.dup(sys.stderr.fileno())
        config.stderr_file = open(outputname, 'w')
        os.dup2(config.stderr_file.fileno(), sys.stderr.fileno())
    elif mode == 'stop':
        if not config.stderr:
            raise Exception('Error redirection has not been started, '
                            'thus it can not be ended.')
        sys.stderr.flush()
        os.dup2(config.stderr, sys.stderr.fileno())
        config.stderr_file.close()


def checkDelPhiErrors(filename, mode=None):
    exceptions = []
    if mode == 'readFiles':
        exceptions = ['part of system outside the box!',
                      'has a net charge of']
    elif mode == 'runDelPhi':
        exceptions = ['part of system outside the box!']

    exit_trigger = False
    errors = ''
    with open(filename) as f:
        for line in f:
            ignore_error = False
            if 'WARNING' in line or 'Warning' in line:
                for exception in exceptions:
                    if exception in line:
                        ignore_error = True
                if not ignore_error:
                    exit_trigger = True
                    errors += line
    if exit_trigger:
        raise Exception('The following errors have been found on {0}: \n{1}'.format(filename, errors))
    if not Config.debug:
        if os.path.isfile(filename):
            os.remove(filename)
        focusing_log = '{}_focusing'.format(filename)
        if os.path.isfile(focusing_log):
            os.remove(focusing_log)
