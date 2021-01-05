# cython: language_level=3
# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False

import numpy as np
from libc.stdio cimport printf, fflush, stdout

cdef int mcsteps, eqsteps, MAXNPKHALFS, MAXDELTA, seed
cdef float LN10, couple_min


LN10        = 2.302585092994
MAXDELTA    = 50
MAXNPKHALFS = 5


cdef extern from "stdlib.h":
    double drand48()
    void srand48(long int seedval)

cdef void mc_step(int nsites, int[:] cur_states, int[:] npossible_states, 
                  float[:, ::1] interactions, int[:, ::1] interactions_lookup,
                  float[:, ::1] possible_states_u, int npairs, int[:] pair1, int[:] pair2):
    cdef int site, site2, state, newstate, site1i, site1newi, \
             state1, state2, site2i, site1, isite1, isite1new, \
             isite2, isite2new, site3, state3, site3i, newstate1, newstate2

    cdef float dU, newG, G, dG, new_interaction, \
               old_interaction, U1_new, U2_new, U1, U2

    for site in range(nsites):
        state  = cur_states[site]
        newstate = int(drand48() * npossible_states[site])

        site1i    = interactions_lookup[site][state]
        site1newi = interactions_lookup[site][newstate]

        dU = possible_states_u[site][newstate] - possible_states_u[site][state]

        for site2 in range(nsites):
            if site != site2:
                state2 = cur_states[site2]
                site2i = interactions_lookup[site2][state2]

                dU += interactions[site1newi][site2i] - interactions[site1i][site2i]

        if accept_move(dU):
            cur_states[site] = newstate

    for npair in range(npairs):
        site1 = pair1[npair]
        site2 = pair2[npair]

        state1 = cur_states[site1]
        state2 = cur_states[site2]

        newstate1 = int(drand48() * npossible_states[site1])
        newstate2 = int(drand48() * npossible_states[site2])
        
        U1     = possible_states_u[site1][state1]
        U2     = possible_states_u[site2][state2]
        U1_new = possible_states_u[site1][newstate1]
        U2_new = possible_states_u[site2][newstate2]

        isite1    = interactions_lookup[site1][state1]
        isite1new = interactions_lookup[site1][newstate1]
        isite2    = interactions_lookup[site2][state2]
        isite2new = interactions_lookup[site2][newstate2]

        old_interaction = interactions[isite1][isite2]
        new_interaction = interactions[isite1new][isite2new]

        dU = U1_new - U1 + U2_new - U2 + new_interaction - old_interaction

        for site3 in range(nsites):
            if site3 != site1 and site3 != site2:
                state3 = cur_states[site3]
                site3i = interactions_lookup[site3][state3]

                dU += interactions[isite1new][site3i] - interactions[isite1][site3i] \
                    + interactions[isite2new][site3i] - interactions[isite2][site3i]

        if accept_move(dU):
            cur_states[site1] = newstate1
            cur_states[site2] = newstate2


cdef int accept_move(float e):
    if e > MAXDELTA:
        return 0
    if e <= 0:
        return 1
    if invexp(e) > drand48():
        return 1
    else:
        return 0


cdef void compute_statistics(float t, int nsites, int [:] cur_states,
                             int [:] avgs, int [:, ::1] possible_states_occ,
                             int [:, ::1] count):
    cdef int site, state
    for site in range(nsites):
        state = cur_states[site]
        avgs[site] += possible_states_occ[site][state]
        count[site][state] += 1


cdef void initialize_u(float pH, int nsites, float [:, :]
                       possible_states_u, int [:] npossible_states,
                       int [:, :] possible_states_occ,
                       float [:, :] possible_states_g, int [:] avgs):
    cdef int site, state
    cdef float u
    for site in range(nsites):
        for state in range(npossible_states[site]):
            u = LN10 * pH * possible_states_occ[site][state] + possible_states_g[site][state]
            possible_states_u[site][state] = u
            #print site, state, LN10, pH, possible_states_occ[site][state], possible_states_g[site][state], possible_states_u[site][state]
        avgs[site] = 0


cdef float invexp(float x):
    """
    This function is a piecewise rational polynomial approximation of
    exp(-x) in the intervals [0,6], [6,13] and [13,70]. Thus, it may
    give significant (even drastic) errors outside the range [0,70].
    """
    cdef float x2, x3, x4
    if x > 13:
        return 8.194236147130614e-10 - 1.3290994520804703e-11 * x
    else:
        x2 = x * x
        x3 = x2 * x
        if x > 6:
            return (-0.0013245823657199278 + 0.00027464252539452071 * x - 0.000019314947607346905 * x2 + 4.598224667374957e-7 * x3) / (1.0 - 0.5165170691890946 * x + 0.09211442135429947 * x2 - 0.006143102546214945 * x3) 
        else:
            x4 = x2 * x2
            return (0.9999965470613797 - 0.3960827416191208 * x + 0.06303500815508939 * x2 - 0.00476617578304489 * x3 + 0.00014392025197088043 * x4) / (1.0 + 0.6038220689877429 * x + 0.16732494517488303 * x2 + 0.026354026827091058 * x3 + 0.00289071552898347 * x4)


cdef select_pairs(int nsites, npossible_states,
                  interactions_lookup,
                  interactions, float couple_min):
    cdef float ggmax, ggmax_tmp
    cdef int npairs, site1, site2, state1, state2, s1, s2, coup
    pair1 = []
    pair2 = []

    #print "### Coupled site pairs, with max|gg| >= {0} pK units:\n".format(couple_min)
    npairs = 0
    for site1 in range(nsites):
        for site2 in range(site1 + 1, nsites):
            coup = 0
            ggmax = 0.0
            for state1 in range(npossible_states[site1]):
                for state2 in range(npossible_states[site2]):
                    s1 = interactions_lookup[site1][state1]
                    s2 = interactions_lookup[site2][state2]
                    ggmax_tmp = interactions[s1][s2]
                    if abs(ggmax_tmp) >= LN10 * couple_min:
                        coup = 1
                        if abs(ggmax_tmp) > ggmax:
                            ggmax = ggmax_tmp
            if coup != 0:
                pair1.append(site1)
                pair2.append(site2)
                npairs += 1
                #print "###   {} {}  {}\n".format(site1, site2, ggmax / LN10)
            
    #print "### Total number of coupled site pairs = {}\n#\n".format(npairs)

    return npairs, pair1, pair2


cpdef MCrun(int nsites, npossible_states_aux,
            possible_states_g_aux, possible_states_occ_aux,
            interactions_aux, interactions_lookup_aux,
            int mcsteps_aux, int eqsteps, int seed, 
            float couple_min, float pH):

    global mcsteps
    mcsteps = mcsteps_aux

    npairs_aux, pair1_aux, pair2_aux = select_pairs(nsites, npossible_states_aux,
                                                    interactions_lookup_aux,
                                                    interactions_aux, couple_min)

    cdef int npairs = npairs_aux

    pair1_arr = np.array(pair1_aux, 'intc')
    cdef int [:] pair1 = pair1_arr

    pair2_arr = np.array(pair2_aux, 'intc')
    cdef int [:] pair2 = pair2_arr

    cdef int i, t, pHstep

    maxstates = max(npossible_states_aux)

    avgs_arr = np.zeros(nsites, dtype='intc')
    cdef int [:] avgs = avgs_arr

    npossible_states_arr = np.array(npossible_states_aux, 'intc')
    cdef int [:] npossible_states = npossible_states_arr

    possible_states_u_arr = np.array([np.zeros(maxstates) for i in range(nsites)], 'single')
    cdef float [:, ::1] possible_states_u = possible_states_u_arr

    possible_states_occ_arr = np.array(possible_states_occ_aux, 'intc')
    cdef int [:, ::1] possible_states_occ = possible_states_occ_arr

    possible_states_g_arr = np.array(possible_states_g_aux, 'single')
    cdef float [:, ::1] possible_states_g = possible_states_g_arr

    cur_states_arr = np.zeros(nsites, dtype='intc')
    cdef int [:] cur_states = cur_states_arr

    interactions_arr = np.array(interactions_aux, 'single')
    cdef float [:, ::1] interactions = interactions_arr

    interactions_lookup_arr = np.array(interactions_lookup_aux, 'intc')
    cdef int [:, ::1] interactions_lookup = interactions_lookup_arr

    pKs_arr = np.array([[100.0 for ii in range(MAXNPKHALFS)] for i in range(nsites)], 'single')
    cdef float [:, ::1] pKs = pKs_arr

    pmean_arr = np.array([0.0 for i in range(nsites)], 'single')
    cdef float [:] pmean = pmean_arr

    count_arr = np.array([np.zeros(maxstates) for i in range(nsites)], 'intc')
    cdef int [:, ::1] count = count_arr

    srand48(seed)

    initialize_u(pH, nsites, possible_states_u, npossible_states,
                 possible_states_occ, possible_states_g, avgs)

    for i in range(eqsteps):
        mc_step(nsites, cur_states, npossible_states,
                interactions, interactions_lookup,
                possible_states_u, npairs, pair1, pair2)
    for t in range(mcsteps):
        mc_step(nsites, cur_states, npossible_states,
                interactions, interactions_lookup,
                possible_states_u, npairs, pair1, pair2)
        compute_statistics(t, nsites, cur_states, avgs,
                           possible_states_occ, count)

    return np.asarray(avgs), np.asarray(pmean), np.asarray(count), np.asarray(cur_states)
