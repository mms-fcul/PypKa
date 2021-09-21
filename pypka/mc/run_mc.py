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
    def __init__(self, sites):
        self.sites = sites
        self.nsites = len(sites)

        self.init_mc_vars()
        self.run()
        self.get_pkas()
        self.get_tit_states()

    def init_mc_vars(self):
        def resize_list_of_lists(listn, maxsize, filler=None):
            for i in listn:
                diff = maxsize - len(i)
                for _ in range(diff):
                    i.append(filler)

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

    def run(self):
        params = Config.mc_params
        self.dpH = params["pHstep"]
        self.pH_values = params.pH_values

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

        for i in results:
            self.avgs_all.append(i[0])
            self.counts_all.append(i[1])
            self.final_states.append(i[3])

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
        text_prots = "#pH       "
        c = -1
        for i in pKas:
            c += 1
            pK = i[0]
            site = self.sites[c]
            site.setpK(pK)
            chain = site.molecule.chain
            sitename = site.getName()
            resnumb = site.getResNumber()
            if sitename in ("NTR", "CTR"):
                text_prots += "     {0:3}".format(sitename)
            else:
                text_prots += "{0:5d}{1:3s}".format(resnumb, sitename)
            text_pks += "{0:5} {1:3} {2:20} {3:3}\n".format(
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
