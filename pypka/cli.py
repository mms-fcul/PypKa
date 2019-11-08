import argparse
import os
import config
import log
from copy import copy

def drange(dmin, dmax, step):
    """Decimal Range

    Requires:
      dmin (float or int) is the range minimum value
      dmax (float or int) is the range maximum value
      step (float or int) is the range step value

    Ensures:
      lrange (list)

    Example:
      drange(5, 7.5, 0.5)
        -> [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]

    """
    lrange = []
    dmin, dmax = float(dmin), float(dmax)
    step = float(step)
    inv_step = step ** -1
    dmin_int = int(int(dmin) * inv_step)
    dmax_int = int(int(dmax) * inv_step) + 1
    for i in range(dmin_int, dmax_int):
        lrange.append(i / inv_step)
    return lrange


def setParameter(param, value):
    if param in list(config.input_conversion.keys()):
        param = config.input_conversion[param]
    config.params[param] = value


def getParameter(param):
    if param in list(config.input_conversion.keys()):
        param = config.input_conversion[param]
    return config.params[param]


def inputParametersFilter(settings):
    """
    Check if input parameters are valid

    Mandatory Parameters:
      - Dielectric Constant
      - Ionic Strength

    """
    config.params = copy(config.default_params)
    config.pid = os.getpid()
    config.script_dir = os.path.dirname(__file__)
    param_names = list(settings.keys())
    # MANDATORY #
    mandatory_params = ('structure', 'epsin', 'ionicstr',
                        'pbc_dimensions', 'ncpus', 'sites_A')
    for param in mandatory_params:
        if param not in param_names:
            log.requiredParameterError(param)
        if type(settings[param]) in (float, int):
            settings[param] = str(settings[param])
        if len(settings[param]) == 0:
            log.requiredParameterError(param)

    # Check parameter conditions: parameter is integer
    integer_params = ('gsize', 'seed', 'ncpus')
    for param in integer_params:
        if param in param_names:
            try:
                param_value = int(settings[param])
            except ValueError:
                log.inputVariableError(param,
                                       'an integer.', '')
            setParameter(param, param_value)

    # Check parameter conditions: parameter is float
    float_params = ('scaleP', 'scaleM', 'ionicstr', 'convergence', 'nlit',
                    'nonit', 'relfac', 'relpar', 'pHstep', 'epssol',
                    'temp', 'epsin', 'slice', 'cutoff')
    for param in float_params:
        if param in param_names:
            try:
                param_value = float(settings[param])
            except ValueError:
                log.inputVariableError(param,
                                       'a float.', '')
            setParameter(param, param_value)

    # Check parameter conditions: parameter > 0
    # These parameters have already been checked for type int or float
    great_params = ('scaleP', 'scaleM', 'convergence', 'pHstep', 'gsize',
                    'ncpus', 'temp', 'grid_fill', 'pH_step')
    for param in great_params:
        if param in param_names:
            param_value = getParameter(param)
            if not param_value > 0:                    
                log.inputVariableError(param,
                                       'greater than zero.', '')
            setParameter(param, param_value)

    # Check parameters conditions: parameter is boolean
    bool_params = ('pbx', 'pby', 'clean_pdb', 'keep_ions')
    for param in bool_params:
        if param in param_names:
            param_value = settings[param]
            if param_value == 'yes':
                param_value = True
            elif param_value == 'no':
                param_value = False
            elif param_value not in (True, False):
                log.inputVariableError(param,
                                       'either "yes" or "no".', '')

            setParameter(param, param_value)

    # Check particular parameter conditions
    if 'bndcon' in param_names:
        if str(settings['bndcon']) not in ('1', '2', '3', '4'):

            log.inputVariableError('bndcon',
                                   '1 (zero), 2(dipolar),'
                                   ' 3(focusing) or 4 (coulombic).', '')
        
        else:
            setParameter('bndcon', int(settings['bndcon']))

    if 'precision' in param_names:
       if settings['precision'] not in ('single', 'double'):

           log.inputVariableError('precision',
                                  'either "single" or "double".', '')
       else:
           setParameter('precision', settings['precision'])

    if 'ffID' in param_names:
        setParameter('ffID', settings['ffID'])
        if settings['ffID'] not in ('G54A7'):  # for now only GROMOS FF
            log.inputVariableError('ffID',
                                   'equal to "G54A7".', '')
    
    if 'ffinput' in param_names:
        setParameter('ffinput', settings['ffinput'])
        if settings['ffinput'] not in ('GROMOS', 'AMBER', 'CHARMM'):
            log.inputVariableError('ffinput',
                                   'either "GROMOS", "AMBER" or "CHARMM".', '')
    
    file_path = os.path.join(config.script_dir, config.params['ffID'])
    config.f_crg = '{0}/DataBaseT.crg'.format(file_path)
    config.f_siz = '{0}/DataBaseT.siz'.format(file_path)

    if settings['pbc_dimensions'] not in ('0', '2'):
        log.inputVariableError('pbc_dimensions',
                               'either "0" or "2".', '')
    else:
        param_value = int(settings['pbc_dimensions'])
        setParameter('pbc_dimensions', param_value)

    if 'nanoshaper' in settings:
        if settings['nanoshaper'] in ('0', '1', '2', '3', '4', '-1'):
            log.inputVariableError('nanoshaper',
                                   'either "-1" to turn of nanoshaper or "0" to turn on nanoshaper.', '')
        else:
            param_value = int(settings['nanoshaper'])
            config.nanoshaper = param_value
        
    # Needs to accept both a single value and a range
    if 'pH' in param_names:
        pH_parts = settings['pH'].split(',')
        if len(pH_parts) > 1:
            try:
                pHmin = float(pH_parts[0])
                pHmax = float(pH_parts[1])
                setParameter('pHmin', pHmin)
                setParameter('pHmax', pHmax)
            except ValueError:
                log.inputVariableError('pH',
                                       'a float.', '')
        else:
            try:
                pHmin = float(pH_parts[0])
                pHmax = float(pH_parts[0])
                setParameter('pHmin', pHmin)
                setParameter('pHmax', pHmax)
            except ValueError:
                log.inputVariableError('pH',
                                       'a float.', '')
            setParameter('pH', [param_value])
        if pHmin >= pHmax:
            log.inputVariableError('pHmax',
                                   'a float greater than pHmin.', '')

    # Declare IO Files
    # Input .pdb File
    config.f_in = settings['structure']
    f_in_parts = settings['structure'].split('.')
    if len(f_in_parts) <= 1:
        log.inputVariableError('structure',
                               'a string containing a file extension.',
                               'Ex: structure.pdb or structure.gro')

    extension = f_in_parts[-1].lower().replace('pqr', 'pdb')
    if extension not in ('gro', 'pdb'):
        log.inputVariableError('structure',
                               'a string containing a valid file extension.',
                               'Ex: structure.pdb or structure.gro or structure.pqr')

    config.f_in_extension = extension

    # Output pKs File
    if 'output' in param_names:
        config.f_out = settings['output']
    else:
        outputname = settings['structure'].split('.')[0]
        config.f_out = outputname

    # Output Titration File
    if 'titration_output' in param_names:
        config.f_prot_out = settings['titration_output']

    if 'structure_output' in param_names:
        error_raise = False
        if len(settings['structure_output']) == 1:
            settings['structure_output'] = settings['structure_output'].split(',')
        if len(settings['structure_output']) != 2:
            error_raise = True

        outfilename = settings['structure_output'][0]
        try:
            pH = float(settings['structure_output'][1])
        except:
            error_raise = True
        if type(outfilename) != str:
            error_raise = True

        if error_raise:
            log.inputVariableError('structure_output',
                               'a tuple containing a filename and the desired pH value.',
                               'Ex: ("structure.pdb", 7)')

        if pH < getParameter('pHmin') or pH > getParameter('pHmax'):
            message = 'pH value for output structure not in range [pHmin, pHmax].'
            log.inputVariableError('structure_output', message, '')
        config.f_structure_out = outfilename
        config.f_structure_out_pH = pH

        if config.sites != {}:
            config.sites = {}
            warning = 'When using f_structure_out all titratable '\
                      'residues are included in the calculation'
            log.reportWarning(warning)

    # Output log File
    if 'logfile' in param_names:
        config.f_log = settings['logfile']  # default: "LOG"

    # Check coherence between variables
    if getParameter('pbc_dim') == 2:
        if getParameter('relfac') != 0.2 and 'relfac' not in settings:
            setParameter('relfac', 0.2)
        if getParameter('nonit') != 5 and 'nonit' not in settings:
            setParameter('nonit', 5)

    if 'lipid_definition' in settings:
        for i in settings['lipid_definition']:
            resname = settings['lipid_definition'][i]
            config.lipids[i] = resname
            if resname in config.lipid_residues:
                resname_i = config.lipid_residues.index(resname)
                del config.lipid_residues[resname_i]
    return


def readSettings(filename):
    """Reads the settings file.

    This file should have the following format:
     - commented lines should being with a '#'
     - every parameter should be declared as such: name = value

    All parameter values are interpreted as strings, however, a type
    check is later performed for each declared input value.

    All parameter names not recognizable are reported as a warning.
    """
    parameters = {}
    parameters['lipid_definition'] = {}
    with open(filename) as f:
        nline = 0
        for line in f:
            nline += 1
            if len(line.strip()) > 0 and line[0] != '#':
                parts = line.split('=')
                param_name = parts[0].strip()
                param_value = '='.join(parts[1:]).strip()
                if 'lipid_definition' in param_name:
                    parts = param_value.split(':')
                    old_name = parts[0]
                    new_name = parts[1]
                    parameters['lipid_definition'][old_name] = new_name
                elif (len(parts) != 2 or
                      not len(param_name) > 0 or
                      not len(param_value) > 0):
                    raise IOError('Incorrect format in line {0} of file {1}: '
                                  '\n{1}#{0}: {2}'.format(nline, filename, line))
                else:
                    parameters[param_name] = param_value

    # Search for all titrable sites in different chains
    sites = {}
    for param_name in parameters:
        if 'site' in param_name:
            chain = param_name.split('_')[1]
            chain_sites = parameters[param_name].split(', ')
            if chain_sites != ['all']:
                sites[chain] = chain_sites
    config.sites = sites

    return parameters


def checkParsedInput():
    """Gets the CLI arguments and interprets them"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description="""
Object-Oriented Script to Calculate the pKint of each site \
as well as the pairwise energies
Requires:
DelPhi2Py module installation
Nanoshaper and cppSolver are optional libraries

Objects:
DelPhiParams stores the DelPhi input parameters like a .prm file

TitratingMolecule is the molecule which has more than one Site

Site

Tautomer

Example:
python pypka.py test.pdb test.dat -o pKas.out --debug

    """)

    # Mandatory Arguments
    parser.add_argument('settings', help=' settings file name',
                        default="settings.dat", action='store')

    # Optional Arguments
    parser.add_argument('--debug', help='activation of the debug mode '
                        'to print extra information', action='store_true')

    args = parser.parse_args()

    # Apply some criteria to input arguments
    if not os.path.isfile(args.settings):
        raise IOError('File {0} does not exist.'.format(args.settings))

    # Read Settings File
    config.f_dat = args.settings
    parameters = readSettings(args.settings)

    return parameters, args.debug
