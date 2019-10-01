import config
from multiprocessing import Pool, Manager, Value
from time import time
from datetime import timedelta, datetime
from sys import stdout


def configRefresh(config, pb_time, njobs):
    config.pb_time = pb_time
    config.njobs = njobs

def printTimeLeft(time1, time2, message):
    config.njobs.value -= 1
    njobs = config.njobs.value
    total_jobs = config.total_jobs
    step = total_jobs - njobs
    config.pb_time.append(time2 - time1)
    left = sum(config.pb_time) / len(config.pb_time) * njobs / config.params['ncpus']
    if left > 3600:
        left_time = f'~{int(left / 3600.0)}h'
    elif left > 60:
        left_time = f'~{int(left / 60.0)}m'
    else:
        left_time = f'{int(left)}s'
    
    end = datetime.now() + timedelta(seconds=left)
    end_time = end.strftime('%H:%M:%S %d/%m/%Y')

    stdout.write(f'\r{message} '
                 f'Run {step:5} of {total_jobs:<10} '
                 f'Ends in {left_time:5} at {end_time:5}')
    stdout.flush()

def startPoolProcesses(targetFunction, iterable_job_arguments_list,
                       ncpus, assign='distributed', merged_results=False):
    """Start a pool of ncpus processes that execute the targetFunction
    with the argument iterable_job_arguments_list

    Args:
      targetFunction (function): function to be execute in each process
      iterable_job_arguments_list (iterable): object enumerating
        the arguments to be passed to the targetFunction.
        when using assing='ordered', it can not be a generator
      ncpus  (int): number of processes to launch
      assign (str): mode to split the jobs among processes
        Example: iterable_job_arguments_list = [0, 1, 2, 3, 4], ncpus = 2
        'distributed' p0 = [0, 2, 4] p1 = [1, 3]
        'ordered' p0 = [0, 1, 2] p1 = [3, 4]

    Ensures:
      unpacked_results (list): list of lists containing the returned
        results from the executed jobs. the order is dependent
        only on $assign and $ncpus but not on the execution time.
    """
    jobs = [[] for i in range(ncpus)]
    results = []
    i = -1
    if assign == 'distributed':
        config.njobs = 0
        for tautomer in iterable_job_arguments_list:
            i += 1
            jobs[i % ncpus].append(i)
            config.njobs += 1
    elif assign == 'ordered':
        max_njobs = int(len(iterable_job_arguments_list) / ncpus)
        ncores_extra_load = len(iterable_job_arguments_list) % ncpus
        core = 0
        core_jobs = 0
        extra_load = False
        config.njobs = 0
        for tautomer in iterable_job_arguments_list:
            i += 1
            if core_jobs == max_njobs:
                if ncores_extra_load > core and not extra_load:
                    core_jobs -= 1
                    extra_load = True
                else:
                    core_jobs = 0
                    core += 1
                    extra_load = False
            jobs[core].append(i)
            core_jobs += 1
            config.njobs += 1
        if config.debug:
            print(f'max_njobs = {max_njobs}, ncores_extra_load = {ncores_extra_load}')

    config.total_jobs = config.njobs
    pb_time = Manager().list()
    njobs = Value('i', config.njobs)
    if config.debug:
        print(f'ncpus = {ncpus}, njobs = {config.njobs}')
        for i, job in enumerate(jobs):
            print(f'ncore {i}: njobs = {len(job)}')
    pool = Pool(processes=ncpus, initializer=configRefresh,
                initargs=(config, pb_time, njobs))
    for job in jobs:
        result = pool.apply_async(targetFunction, args=(job, ))
        results.append(result)
        # Easier debug of the loop but fails afterwards
        #targetFunction(job)

    pool.close()
    pool.join()
    #print 'exit'
    #exit()
    
    unpacked_results = []
    for results_percore in results:
        result = results_percore.get()
        if merged_results:
            unpacked_results += result
        else:
            unpacked_results.append(result)
    return unpacked_results


def runDelPhiSims(job_list):
    """Run sets of simulations of tautomers included in job_list

    Args:
      job_list (list): list of tautomer numbers to be analysed.
        The use of numbers is needed since class objects
        can not be pickled and therefore can not be parsed
        in multiprocessing.

    Ensures:
      results (list)
        tauname (str):   the tautomer name
        sitenum (int):   site id number
        esolvM  (float): the Molecule solvation energy
        sitpotM (list):  the site potential in the Molecule
          the index corresponds to the atoms in the site
        esolvS  (float): the Site solvation energy
        sitpotS (list):  the site potential in the Site
          the index corresponds to the atoms in the site
    """
    results = []
    for tau_number in job_list:
        time1 = time()
        tauname, sitenum, esolvM, sitpotM, esolvS, sitpotS = calcPotential(tau_number)
        results.append([tauname, sitenum, esolvM, sitpotM, esolvS, sitpotS])
        time2 = time()

        message = f'PB Runs: {tauname:3} {sitenum:<10}'

        printTimeLeft(time1, time2, message)

    stdout.write(f'\rPB Runs Ended{"":>80}')
    stdout.flush()

    return results


def calcPotential(taut):
    """Run two DelPhi simulations: one for the tautomer and
    other for the same tautomer when all other sites are neutral.

    Args:
      taut (int): number of the tautomer

    Ensures:
      tauname (str):   the tautomer name
      sitenum (int):   site id number
      esolvM  (float): the Molecule solvation energy
      sitpotM (list):  the site potential in the Molecule
        the index corresponds to the atoms in the site
      esolvS  (float): the Site solvation energy
      sitpotS (list):  the site potential in the Site
        the index corresponds to the atoms in the site
    """
    # get the tautomer object from the tautomer number
    # transformation needed due to multiprocessing
    taut = config.tit_mole.getTautomerNumber(taut)

    # Whole Molecule with zero charge except for the tautomer being evaluated
    e_solvationM, sitpotM = taut.CalcPotentialTitratingMolecule()

    if config.debug:
        print(('finished Whole Molecule', taut._name, e_solvationM))

    # Single Site Tautomer
    e_solvationS, sitpotS = taut.CalcPotentialTautomer()

    if config.debug:
        print(('finished Tautomer', taut._name, e_solvationS))

    return (taut._name, taut.getSiteResNumber(), e_solvationM,
            sitpotM, e_solvationS, sitpotS)


def runInteractionCalcs(job_list):
    """Run interaction calculations between
    sites in job_list

    Args:
      job_list (list)
        interaction_number (int): interaction number to be analysed.
          The use of numbers is needed since class objects
          can not be pickled and therefore can not be parsed
          in multiprocessing.

    Ensures:
      results (list)
        interaction_energies (str): .dat formatted energy
          interactions between all tautomers of two sites
    """
    results = []
    for interaction_number in job_list:
        interaction_energies = config.tit_mole.calcInteractionNumber(interaction_number)
        results += interaction_energies

    return results


def runMCCalcs(job_list):
    """
    """
    results = []
        
    for pH_index in job_list:
        time1 = time()
        mc_output = config.tit_mole.parallelMCrun(pH_index)
        results.append(mc_output)
        time2 = time()

        message = f'MC'

        printTimeLeft(time1, time2, message)

    return results
