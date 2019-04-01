import config
import os
import sys


def reportToLOG(info):
    with open(config.f_log, 'a') as logfile:
        logfile.write(info + '\n')


def requiredParameterError(parameter):
    raise IOError('Required input parameter {0} '
                  'is not defined.'.format(parameter))


def inputVariableError(parameter, complaint, explanation):
    raise IOError('Input parameter {0} is not {1}\n '
                  '{2}'.format(parameter, complaint, explanation))


def reportWarning(info):
    warning = 'warning: {0}'.format(info)
    reportToLOG(warning)

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

def checkDelPhiErrors(filename, remove=False):
    exit_trigger = False
    errors = ''
    with open(filename) as f:
        for line in f:
            if ('WARNING' in line or 'Warning' in line) and \
               'part of system outside the box!' not in line:
                exit_trigger = True
                errors += line
    if exit_trigger:
        raise Exception('The following errors have been found on {0}: \n{1}'.format(filename, errors))
    if remove:
        os.remove(filename)
