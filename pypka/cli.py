import argparse
import os

from pypka.main import Titration
from pypka.__init__ import __version__


def read_settings(filename):
    """
    Reads the settings file.

    This file should have the following format:
      - commented lines should being with a '#'
      - every parameter should be declared as such: name = value

    All parameter values are interpreted as strings, however, a type
        check is later performed for each declared input value.
    All parameter names not recognizable are reported as a warning.
    """
    parameters = {}
    parameters["lipid_definition"] = {}
    with open(filename) as f:
        nline = 0
        for line in f:
            nline += 1
            line = line.strip()
            if line and line[0] != "#":
                no_comments = line.split("#")[0]
                parts = no_comments.split("=")
                param_name = parts[0].strip()
                param_value = "=".join(parts[1:]).strip()
                if "lipid_definition" in param_name:
                    parts = param_value.split(":")
                    old_name = parts[0]
                    new_name = parts[1]
                    parameters["lipid_definition"][old_name] = new_name
                elif len(parts) != 2 or not param_name or not param_value:
                    raise IOError(
                        "Incorrect format in line {0} of file {1}: "
                        "\n{1}#{0}: {2}".format(nline, filename, line)
                    )
                else:
                    parameters[param_name] = param_value
    # Search for all titrable sites in different chains
    sites = {}
    fixed_sites = {}
    for param_name in parameters:
        if "_" in param_name:
            chain = param_name.split("_")[-1]
        else:
            chain = " "
        if param_name.startswith("site"):
            sites_str = parameters[param_name]
            chain_sites = [i.strip() for i in sites_str.split(",")]
            sites[chain] = chain_sites
        elif param_name.startswith("fixed_site"):
            fixed_sites_str = parameters[param_name]
            fixed_sites[chain] = eval(fixed_sites_str)

    if " " in sites and sites[" "] == ["all"]:
        if len(sites) != 1:
            raise Exception(
                '"sites" parameter is incorrectly defined.\n '
                'Incompatible parameters: "sites" and the multi-chain nomenclature "sites_X" where X is the chain'
                "Please choose one of the two ways to define titratble sites.\n"
                "Automatic Example: sites = all\n"
                "Residues NTR, 1, 5 and 8 in chain A \nExample: sites_A = 1N, 1, 5, 8"
            )
        sites = "all"

    if not sites:
        raise Exception(
            'Missing mandatory "sites" parameter\nTitratate all sites: sites = all\nTitrate all sites only from chain A: sites_A = all'
        )

    return sites, fixed_sites, parameters


def check_cli_args():
    """Gets the CLI arguments and interprets them"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
PypKa
A python module for flexible Poisson-Boltzmann based pKa calculations with proton tautomerism

Documentation can be found at https://pypka.readthedocs.io/en/latest/

PypKa is distributed under a LGPL-3.0, however delphi4py depends on DelPhi which is proprietary.
To use DelPhi the user is required to download the DelPhi license https://honiglab.c2b2.columbia.edu/software/cgi-bin/software.pl?input=DelPhi

Example:
python3 pypka.py params.dat --debug""",
    )

    # Mandatory Arguments
    parser.add_argument(
        "settings", help=" settings file name", default="settings.dat", action="store"
    )

    # Optional Arguments
    parser.add_argument(
        "--debug",
        help="activation of the debug mode " "to print extra information",
        action="store_true",
    )

    parser.add_argument("--version", action="version", version=__version__)

    args = parser.parse_args()

    # Apply some criteria to input arguments
    if not os.path.isfile(args.settings):
        raise IOError("File {0} does not exist.".format(args.settings))

    # Read Settings File
    sites, fixed_sites, parameters = read_settings(args.settings)

    return sites, fixed_sites, parameters, args.debug


def CLI():
    # Read command line arguments
    sites, fixed_sites, parameters, debug = check_cli_args()


    Titration(parameters, sites=sites, fixed_sites=fixed_sites, debug=debug)
    print("CLI exited successfully")
