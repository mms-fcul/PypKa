import json
import numpy as np
from pypka.config import Config
from pypka.concurrency import runMCCalcs, startPoolProcesses
from pypka.constants import KBOLTZ, PKAPLACEHOLDER, MAXNPKHALFS


def calcpKhalfs(pH, nsites, avgs, pmean, pKs, dpH):
    for site in range(nsites):
        mean = avgs[site]
        p = pmean[site]

        if p > 0.5 and mean <= 0.5 or p < 0.5 and mean >= 0.5:
            pKhalf = pH - dpH * (mean - 0.5) / (mean - p)
            for i in range(MAXNPKHALFS):
                if pKs[site][i] == PKAPLACEHOLDER:
                    pKs[site][i] = pKhalf
                    break

        elif p > 1.5 and mean <= 1.5 or p < 1.5 and mean >= 1.5:
            pKhalf = pH - dpH * (mean - 1.5) / (mean - p)
            for i in range(MAXNPKHALFS):
                if pKs[site][i] == PKAPLACEHOLDER:
                    pKs[site][i] = pKhalf
                    break
    return pKs


class MonteCarlo:
    """Monte Carlo Simulation"""

    def __init__(self, sites):
        self.sites = sites
        self.nsites = len(sites)

        self.init_mc_vars()
        self.run()
        self.get_pkas()
        self.get_tit_states()
        self.get_coupled_sites()

    def init_mc_vars(self):
        """Initialize need variables"""

        def resize_list_of_lists(listn, maxsize, filler=None):
            for i in listn:
                diff = maxsize - len(i)
                for _ in range(diff):
                    i.append(filler)

        states_ddG = [[] for _ in self.sites]
        possible_states_g = [[] for _ in self.sites]
        possible_states_occ = [[] for _ in self.sites]

        temperature = float(Config.pypka_params["temp"])
        isite = -1
        for site in self.sites:
            isite += 1
            itaut = 0
            for tautomer in site.iterOrderedTautomersWithoutRef():
                dg = tautomer.dg / (KBOLTZ * temperature)
                possible_states_g[isite].append(dg)

                ddg = tautomer.dG_solvationM - tautomer.dG_solvationS + tautomer.dG_back
                states_ddG[isite].append(ddg)

                if site.type == "c":
                    prot_state = 0
                elif site.type == "a":
                    prot_state = 1

                possible_states_occ[isite].append(prot_state)
                itaut += 1

            if site.type == "c":
                prot_state = 1
            elif site.type == "a":
                prot_state = 0
            possible_states_occ[isite].append(prot_state)
            possible_states_g[isite].append(0.0)

        maxstates = max(Config.parallel_params.npossible_states)
        interactions_look = Config.parallel_params.interactions_look

        resize_list_of_lists(possible_states_g, maxstates)
        resize_list_of_lists(possible_states_occ, maxstates, filler=-500)
        resize_list_of_lists(interactions_look, maxstates, filler=-500)

        Config.parallel_params.possible_states_g = possible_states_g
        Config.parallel_params.possible_states_occ = possible_states_occ
        Config.parallel_params.states_ddG = states_ddG

    def run(self):
        """Run a MC sim"""
        params = Config.mc_params
        self.dpH = params["pHstep"]
        self.pH_values = params.pH_values

        if Config.pypka_params.save_mc_energies:
            all_sites = []
            for site in Config.parallel_params.all_sites:
                site_id = f"{site.molecule.chain}_{site.res_name}_{site.res_number}"
                all_sites.append(site_id)

            to_json_dict = {
                "all_sites": all_sites,
                "npossible_states": Config.parallel_params.npossible_states,
                "possible_states_g": Config.parallel_params.possible_states_g,
                "states_ddG": Config.parallel_params.states_ddG,
                "possible_states_occ": Config.parallel_params.possible_states_occ,
                "interactions": Config.parallel_params.interactions,
                "interactions_look": Config.parallel_params.interactions_look,
            }

            json_object = json.dumps(to_json_dict, indent=4)

            with open(Config.pypka_params.save_mc_energies, "w") as f_out:
                f_out.write(json_object)

        ncpus = min(Config.pypka_params["ncpus"], self.nsites)
        results = startPoolProcesses(
            runMCCalcs,
            self.pH_values,
            ncpus,
            assign="ordered",
            merged_results=True,
        )

        self.counts_all = []
        self.avgs_all = []
        self.final_states = []
        self.coupled_sites_pairs = []

        for i in results:
            self.avgs_all.append(i[0])
            self.counts_all.append(i[1])
            self.final_states.append(i[3])
            self.coupled_sites_pairs.append(i[4])
        self.coupled_sites_pairs = self.coupled_sites_pairs[0]

    def get_pkas(self):
        mcsteps = Config.mc_params["mcsteps"]

        pKs = np.array(
            [[PKAPLACEHOLDER for ii in range(MAXNPKHALFS)] for i in range(self.nsites)]
        )

        self.total_tit_curve = {}
        for pHstep, pH in enumerate(self.pH_values):
            ph_avgs_prots = self.avgs_all[pHstep] / mcsteps
            self.total_tit_curve[pH] = sum(ph_avgs_prots) / self.nsites

            if pHstep > 0:
                pKs = calcpKhalfs(
                    pH, self.nsites, ph_avgs_prots, prev_prots, pKs, self.dpH
                )

            prev_prots = ph_avgs_prots[:]

        pKas = pKs

        text_pks = ""
        text_prots = "#pH      "
        c = -1
        for i in pKas:
            c += 1
            pK = i[0]
            site = self.sites[c]
            site.setpK(pK)
            chain = site.molecule.chain
            sitename = site.getName()
            resnumb = site.getResNumber(correct_icode=True)
            text_prots += "{0}_{2}_{1} ".format(resnumb, sitename, chain)
            text_pks += "{0:>5} {1:3} {2:20} {3:3}\n".format(
                resnumb, sitename, str(pK), chain
            )
        self.text_pks = text_pks
        self.text_prots = text_prots

    def get_tit_states(self):
        mcsteps = Config.mc_params["mcsteps"]

        for pHstep, pH in enumerate(self.pH_values):
            self.text_prots += "\n{pH:5.2f}".format(pH=pH)

            ph_avgs_prots = self.avgs_all[pHstep] / mcsteps
            ph_taut_probs = self.counts_all[pHstep] / mcsteps

            for c, site in enumerate(self.sites):
                mean = ph_avgs_prots[c]
                res_state_dist = list(ph_taut_probs[c])

                ntauts = site.getNTautomers()
                ref_i = ntauts
                prot_state = site.getRefProtState()
                if (mean > 0.5 and prot_state == 1) or (
                    mean <= 0.5 and prot_state == -1
                ):
                    state_i = ref_i
                else:
                    max_prob = max(res_state_dist[:ref_i])
                    state_i = res_state_dist.index(max_prob)

                most_prob_state = state_i + 1

                if mean != PKAPLACEHOLDER:
                    self.text_prots += "\t{mean:7.4f}".format(mean=mean)
                else:
                    self.text_prots += "\t-"

                site.most_prob_states[pH] = most_prob_state
                site.final_states[pH] = self.final_states[pHstep][c]
                site.tit_curve[pH] = mean
                site.states_prob[pH] = res_state_dist

    def get_coupled_sites(self):
        """Get coupled sites"""
        self.coupled_sites_dict = {}
        for sitei1, sitei2 in self.coupled_sites_pairs:
            site1 = Config.parallel_params.all_sites[sitei1]
            site2 = Config.parallel_params.all_sites[sitei2]

            if site1 not in self.coupled_sites_dict:
                self.coupled_sites_dict[site1] = []
            if site2 not in self.coupled_sites_dict:
                self.coupled_sites_dict[site2] = []
            self.coupled_sites_dict[site1].append(site2)
            self.coupled_sites_dict[site2].append(site1)

        self.text_coupled_sites = ""
        for site1, sites in self.coupled_sites_dict.items():
            to_print = "{:3} {:3} {:6} | {:6} -> ".format(
                site1.molecule.chain, site1.res_name, site1.res_number, len(sites)
            )

            for site2 in sites:
                to_print += "{:3} {:3} {:6}, ".format(
                    site2.molecule.chain, site2.res_name, site2.res_number
                )
            self.text_coupled_sites += to_print[:-2] + "\n"
