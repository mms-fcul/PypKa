from datetime import datetime, timedelta
from multiprocessing import Manager, Pool, Value
from sys import stdout
from time import time

from pypka.config import Config
from pypka.mc.mc import MCrun


def configRefresh(configs, pb_time, njobs):
    configs.pb_time = pb_time
    configs.njobs = njobs


def printTimeLeft(time1, time2, message):
    Config.parallel_params.njobs.value -= 1
    njobs = Config.parallel_params.njobs.value
    total_jobs = Config.parallel_params.total_jobs
    step = total_jobs - njobs
    Config.parallel_params.pb_time.append(time2 - time1)
    pb_time = Config.parallel_params.pb_time
    left = sum(pb_time) / len(pb_time) * njobs / Config.pypka_params["ncpus"]
    if left > 3600:
        left_time = "~{0}h".format(int(left / 3600.0))
    elif left > 60:
        left_time = "~{0}m".format(int(left / 60.0))
    else:
        left_time = "{0}s".format(int(left))

    end = datetime.now() + timedelta(seconds=left)
    end_time = end.strftime("%H:%M:%S %d/%m/%Y")

    stdout.write(
        "\r{0} "
        "Run {1:5} of {2:<10} "
        "Ends in {3:5} at {4:5}".format(message, step, total_jobs, left_time, end_time)
    )
    stdout.flush()


def startPoolProcesses(
    targetFunction,
    iterable_job_arguments_list,
    ncpus,
    assign="distributed",
    merged_results=False,
):
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
    if assign == "distributed":
        Config.parallel_params.njobs = 0
        for tautomer in iterable_job_arguments_list:
            i += 1
            jobs[i % ncpus].append(i)
            Config.parallel_params.njobs += 1
    elif assign == "ordered":
        max_njobs = int(len(iterable_job_arguments_list) / ncpus)
        ncores_extra_load = len(iterable_job_arguments_list) % ncpus
        core = 0
        core_jobs = 0
        extra_load = False
        Config.parallel_params.njobs = 0
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
            Config.parallel_params.njobs += 1
        if Config.debug:
            print(
                "max_njobs = {}, ncores_extra_load = {}".format(
                    max_njobs, ncores_extra_load
                )
            )

    Config.parallel_params.total_jobs = Config.parallel_params.njobs
    pb_time = Manager().list()
    njobs = Value("i", Config.parallel_params.njobs, lock=False)
    Config.parallel_params.njobs = njobs

    if Config.debug:
        print("ncpus = {}, njobs = {}".format(ncpus, Config.parallel_params.njobs))
        for i, job in enumerate(jobs):
            print("ncore {}: njobs = {}".format(i, len(job)))

    with Pool(
        processes=ncpus,
        initializer=configRefresh,
        initargs=(Config.parallel_params, pb_time, njobs),
    ) as pool:
        for job in jobs:
            result = pool.apply_async(targetFunction, args=(job,))
            results.append(result)
            # Easier debug of the loop but fails afterwards
            # targetFunction(job)

        pool.close()
        pool.join()
    # print('exit')
    # exit()

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
        tauname, sitenum, chain, esolvM, sitpotM, esolvS, sitpotS = calcPotential(
            tau_number
        )
        results.append([tauname, sitenum, esolvM, sitpotM, esolvS, sitpotS])
        time2 = time()

        message = "PB Runs: {:1} {:3} {:<10} ".format(chain, tauname, sitenum)

        printTimeLeft(time1, time2, message)

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
    taut = Config.parallel_params.all_tautomers_order[taut]

    # Whole Molecule with zero charge except for the tautomer being evaluated
    e_solvationM, sitpotM = taut.CalcPotentialTitratingMolecule()

    if Config.debug:
        print(("finished Whole Molecule", taut.name, e_solvationM))

    # Single Site Tautomer
    e_solvationS, sitpotS = taut.CalcPotentialTautomer()

    if Config.debug:
        print(("finished Tautomer", taut.name, e_solvationS))

    return (
        taut.name,
        taut.getSiteResNumber(),
        taut.molecule.chain,
        e_solvationM,
        sitpotM,
        e_solvationS,
        sitpotS,
    )


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
        site1, site2 = Config.parallel_params.site_interactions[interaction_number]
        interaction_energies = site1.calc_interaction_between(site2)
        results += interaction_energies
    return results


def parallelMCrun(i):
    params = Config.mc_params
    parallel_params = Config.parallel_params
    pH = params.pH_values[i]

    sites = Config.parallel_params.all_sites
    nsites = len(sites)

    avgs, pmean, count, cur_states = MCrun(
        nsites,
        parallel_params.npossible_states,
        parallel_params.possible_states_g,
        parallel_params.possible_states_occ,
        parallel_params.interactions,
        parallel_params.interactions_look,
        params["mcsteps"],
        params["eqsteps"],
        params["seed"],
        params["couple_min"],
        pH,
    )
    return (avgs, count, pmean, cur_states)


def runMCCalcs(job_list):
    results = []

    for pH_index in job_list:
        time1 = time()
        mc_output = parallelMCrun(pH_index)
        results.append(mc_output)
        time2 = time()

        message = "MC"

        printTimeLeft(time1, time2, message)

    return results
