import pytest
import os
import subprocess as sb
import sys
from pypka.tests.builder.check_diffs import compareFiles

sys.path.insert(1, '../')

ncpus = 4

# run coverage.sh to generate coverage

def runTest(path, ncpus, results):
    os.system("rm -f {0}/*out {0}/clean*pqr {0}/TMP.gro {0}/delphi_in_*pdb".format(path))
    results_lines = results.split('\n')[1:-1]
    #"python3 -m coverage erase; "
    #"python3 -m coverage run ../../../pypka.py parameters.dat;"
    sb.Popen("cd {0}; "
             "sed -i 's/ncpus .*/ncpus           = {1}/' parameters.dat; "
             "python3 ../../../pypka.py parameters.dat;"
             "rm -f *-*.pdb *.profl *.prm *.xvg *.frc *.crg cent "
             "contributions interactions.dat pkint tmp.sites".format(path, ncpus), shell=True).wait()
    checkOutput('{0}/pKas.out'.format(path), results_lines)

def checkOutput(filename, results_lines):
    with open(filename) as f:
        c = -1
        for line in f:
            c += 1
            line = line.strip()
            assert line == results_lines[c]
    assert c + 1 == len(results_lines)

def checkStructureOutput(filename):
    problems = compareFiles(f'builder/{filename}', filename)
    if problems:
        raise Exception(f'Problems found with {filename}')
    else:
        os.remove(filename)


def checkAPIResult(pKa, results):
    i = -1
    for site in pKa:
        result = '{0} {1} {2} {3}'.format(site.res_name,
                                          site.res_number,
                                          site.pK,
                                          site.getProtState(7))
        i += 1
        assert result == results[i]
    assert i + 1 == len(results)

class TestCLI(object):
    def test_cli_ktp_gro(self):
        path = "ktp/ktp_gro"
        results = """
1 NTR 7.9006117525417885
1 TYR 9.88518370784027
2 CTR 3.123787446504994
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_allsites(self):
        path = "ktp/ktp_pdb_allsites"
        results = """
1 NTR 7.903609060282417
1 TYR 9.783802064791741
2 CTR 3.3021160843607587
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_allsites_noclean(self):
        path = "ktp/ktp_pdb_allsites_noclean"
        results = """
1 NTR 7.9006117525417885
1 TYR 9.88518370784027
2 CTR 3.123787446504994
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_onlytermini(self):
        path = 'ktp/ktp_pdb_onlytermini'
        results = """
1 NTR 7.900161128244222
2 CTR 3.2974246565617125
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_onlytermini_noclean(self):
        path = 'ktp/ktp_pdb_onlytermini_noclean'
        results = """
1 NTR 8.167721738527177
2 CTR 3.390519247038917
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_gro(self):
        path = "lyso/lyso_gro"
        results = """
1 NTR 7.636678712935952
18 ASP 3.140817760292259
35 GLU 4.663003222341569
48 ASP 2.964461211978307
66 ASP 3.255524771073022
129 CTR 2.3530479084733646
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_sites(self):
        path = "lyso/lyso_pdb_sites"
        results = """
1 NTR 7.506180977656024
18 ASP 3.1878485151212437
35 GLU 4.608229404531882
48 ASP 2.2579192691339616
66 ASP 1.820122723288205
129 CTR 2.3052342318764727
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_sites_noclean(self):
        path = "lyso/lyso_pdb_sites_noclean"
        results = """
1 NTR 7.393383742911153
18 ASP 2.7812887935898885
35 GLU 4.271359971006614
48 ASP 2.6547679911654636
66 ASP 3.010625039606261
129 CTR 1.9407914909927175
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_all(self):
        path = "lyso/lyso_pdb_all"
        results = """
1 NTR 7.385486289829849
1 LYS 10.385584365325077
7 GLU 2.8826946555866253
13 LYS None
15 HIS 5.867371279679988
18 ASP 2.699499030381383
20 TYR 9.911190805007351
23 TYR 9.392016503352243
24 SER None
33 LYS 10.512639059944679
35 GLU 4.175352330803779
36 SER None
40 THR None
43 THR None
47 THR None
48 ASP 1.859099446300949
50 SER None
51 THR None
52 ASP 2.1516355521689436
53 TYR 10.738445987683207
60 SER None
66 ASP 2.337319641232289
69 THR None
72 SER None
81 SER None
85 SER None
86 SER None
87 ASP 2.0567544019613804
89 THR None
91 SER None
96 LYS None
97 LYS None
100 SER None
101 ASP 3.31453580803786
116 LYS 10.253852282222178
118 THR None
119 ASP 2.3967126151939904
129 CTR 1.7030741123202797
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_all_noclean(self):
        path = "lyso/lyso_pdb_all_noclean"
        results = """
1 NTR 7.384356794782545
1 LYS 10.389604559472225
7 GLU 2.883256437556686
13 LYS None
15 HIS 5.857142857142857
18 ASP 2.704063319217584
20 TYR 9.912713556689594
23 TYR 9.394213400114952
24 SER None
33 LYS 10.516970339406788
35 GLU 4.177328004243067
36 SER None
40 THR None
43 THR None
47 THR None
48 ASP 1.8664671985815602
50 SER None
51 THR None
52 ASP 2.1504670878557155
53 TYR 10.748483872190425
60 SER None
66 ASP 2.3538356252272292
69 THR None
72 SER None
81 SER None
85 SER None
86 SER None
87 ASP 2.047742733457019
89 THR None
91 SER None
96 LYS None
97 LYS None
100 SER None
101 ASP 3.3119866025485094
116 LYS 10.254489298936823
118 THR None
119 ASP 2.397993604492454
129 CTR 1.702982198203926
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_multichain_all(self):
        path = "lyso/lyso_multichain_all"
        results = """
1 NTR 7.183462163467873
1 LYS 10.396268714909544
7 GLU 2.7881048585788903
13 NTR 2.4753604406285437
13 LYS None
15 HIS 5.766867861454174
18 ASP 2.409258450943305
20 CTR None
20 TYR 9.978007803338805
21 NTR None
23 TYR 8.51584948914053
24 SER None
33 LYS 10.516964051164953
35 GLU 4.053515554030321
36 SER None
40 THR None
41 NTR 6.30871506642378
43 THR None
47 THR None
48 ASP 1.8451897981859995
50 SER None
51 THR None
52 ASP 2.002075936318864
53 TYR 10.705104182281517
60 SER None
66 ASP 2.2240888633956937
69 THR None
72 SER None
81 SER None
85 SER None
86 SER None
87 ASP 1.891905748004584
89 THR None
91 SER None
96 LYS None
97 LYS None
100 SER None
101 ASP 3.338910637657156
116 LYS 10.240017754608644
118 THR None
119 ASP 2.355701271398538
129 CTR 1.8967363230805352
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_multichain_sites(self):
        path = "lyso/lyso_multichain_sites"
        results = """
1 NTR 7.2083941115933605
1 ASP 3.0432687937460243
35 GLU 4.612318525602017
35 ASP 2.2543511178409195
66 ASP 1.8262090921859637
1 CTR 2.702574220217964
        """
        runTest(path, ncpus, results)

    def test_cli_pHLIP_gro(self):
        path = "pHLIP/pHLIP_gro"
        results = """
769 NTR 10.7649336974175
770 CYS None
771 GLU 3.997796390277407
782 ASP None
793 ASP None
799 ASP 1.9435900082576383
801 ASP 0.9864795312376277
802 GLU 4.057797348110034
804 CTR 3.837918510075153
        """
        runTest(path, ncpus, results)

    def test_cli_pHLIP_pdb_all(self):
        path = "pHLIP/pHLIP_pdb_all"
        results = """
769 NTR 11.20572255781535
770 CYS 24.212240582368022
771 GLU 4.006244338659009
776 TYR 18.59753191489362
780 TYR 24.40723419370116
782 ASP None
786 THR None
787 THR None
793 ASP 22.721932396697188
799 ASP 0.4165757481296758
801 ASP 1.4362119419712347
802 GLU 3.945716473628793
804 CTR 4.773451582666378
804 THR None
        """
        runTest(path, ncpus, results)

    def test_cli_nucleosome_pdb_all(self):
        path = "nucleosome/nucleosome_pdb_all"
        results = """
1 NTR 7.882393397524072
3 THR None
4 LYS 10.290327706758951
6 THR None
9 LYS 8.516685074155886
10 SER None
11 THR None
14 LYS 10.06376595512564
18 LYS 7.63923256980438
22 THR None
23 LYS 9.672642410098653
27 LYS 10.839709485474273
28 SER None
32 THR None
36 LYS 10.414026770775237
37 LYS 11.254282087233317
39 HIS 8.760644085644085
41 TYR 13.72589107854297
45 THR None
50 GLU 3.0977474008471315
54 TYR 10.836061671427526
56 LYS 10.461906877554647
57 SER None
58 THR None
59 GLU 2.972853139247919
64 LYS 10.749247888651274
73 GLU 3.4225222172013856
77 ASP 4.501544770250534
79 LYS 10.883174614330393
80 THR None
81 ASP 2.58587883648028
86 SER None
87 SER None
94 GLU 2.336822582673205
96 SER None
97 GLU 0.2736797240504321
99 TYR 12.748583252553873
105 GLU 4.3537843607220434
106 ASP 2.758237453090849
107 THR None
110 CYS 13.30342405281947
113 HIS 3.462506496399139
115 LYS 10.421432084021216
118 THR None
122 LYS 8.99167708946439
123 ASP None
133 GLU 0.7110435350307528
135 CTR 3.5171617220137366
"""
        runTest(path, ncpus, results)


class TestAPI(object):
    def test_api_ktp_gro(self):
        from ..pypka import Titration
        os.system('rm -f LOG* *out *gro *pdb *pqr *crg *sites cent contributions interactions.dat pkint')
        results = """
NTR 2001 7.9006117525417885 ('undefined', 0.8883240472631488)
TYR 1 9.88518370784027 ('protonated', 0.99869907894847)
CTR 2002 3.123787446504994 ('deprotonated', 0.0001329626608755893)
        """
        results = results.split('\n')[1:-1]
        parameters = {
            'structure'     : 'ktp/ktp_gro/ktp.gro',     # MANDATORY
            'epsin'         : 2,
            'ionicstr'      : 0.1,
            'pbc_dimensions': 0,
            'temp'          : 310,
            'grid_fill'     : 0.8,         # FUTURE VERSION
            'ncpus'         : ncpus,
            'pH'            : '-5,15',
            'pHstep'        : 0.2,
            'logfile'       : 'LOGFILE',
            'scaleM'        : 4,
            'scaleP'        : 1,
            'gsize'         : 81,
            'convergence'   : 0.01,
            'nlit'          : 500,
            'cutoff'        : -1,
            'relfac'        : 0.0,
            'output'        : 'pKas.out',
            'clean_pdb'     : False
        }
        sites = {'A': ('1N', '1', '2C')}
        pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_allsites_clean(self):
        from ..pypka import Titration
        os.system('rm -f *out *gro *pdb *pqr *crg *sites cent contributions interactions.dat pkint')
        results = """
NTR 2001 7.903967864316001 ('undefined', 0.8890883751899552)
TYR 1 9.781315673289184 ('protonated', 0.9983481661893414)
CTR 2002 3.3021115745568306 ('deprotonated', 0.0002004585145523701)
        """
        results = results.split('\n')[1:-1]
        parameters = {
            'structure'     : 'ktp/ktp_pdb_allsites/ktp.pdb',
            'clean_pdb'     : 'yes',
            'epsin'         : 2,
            'ionicstr'      : 0.1,
            'pbc_dimensions': 0,
            'temp'          : 310,
            'grid_fill'     : 0.8,
            'ncpus'         : ncpus,
            'pH'            : '0,15',
            'pHstep'        : 0.2,
            'logfile'       : 'LOGFILE',
            'scaleM'        : 4,
            'scaleP'        : 1,
            'gsize'        : 81,
            'convergence'   : 0.01,
            'nlit'          : 300,
            'cutoff'        : -1,
            'relfac'        : 0.0,
            'output'        : 'pKas.out',
            'clean_pdb'     : True
        }
        sites = 'all'
        pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_allsites_noclean(self):
        from ..pypka import Titration
        os.system('rm -rf LOG* *out *gro *pdb *pqr *crg *sites '
                  'cent contributions interactions.dat pkint LOG* __pycache__ *pyc')
        results = """
NTR 2001 7.9006117525417885 ('undefined', 0.8883240472631488)
TYR 1 9.88518370784027 ('protonated', 0.99869907894847)
CTR 2002 3.123787446504994 ('deprotonated', 0.0001329626608755893)
        """
        results = results.split('\n')[1:-1]
        parameters = {
            'structure'     : 'ktp/ktp_pdb_allsites_noclean/ktp_noclean.pdb',     # MANDATORY
            'clean_pdb'     : 'no',
            'epsin'         : 2,
            'ionicstr'      : 0.1,
            'pbc_dimensions': 0,
            'temp'          : 310,
            'grid_fill'     : 0.8,         # FUTURE VERSION
            'ncpus'         : ncpus,
            'pH'            : '-5,15',
            'pHstep'        : 0.2,
            'logfile'       : 'LOGFILE',
            'scaleM'        : 4,
            'scaleP'        : 1,
            'gsize'         : 81,
            'convergence'   : 0.01,
            'nlit'          : 500,
            'cutoff'        : -1,
            'relfac'        : 0.0,
            'output'        : 'pKas.out',
            'clean_pdb'     : False
        }
        sites = 'all'
        pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_onlytermini(self):
        from ..pypka import Titration
        os.system('rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent '
                  'contributions interactions.dat pkint LOG* __pycache__ *pyc')
        results = """
NTR 2001 7.900161128244222 ('undefined', 0.8882210711812845)
CTR 2002 3.2974246565617125 ('deprotonated', 0.000198307219057625)
        """
        results = results.split('\n')[1:-1]
        parameters = {
            'structure'     : 'ktp/ktp_pdb_onlytermini/ktp.pdb',     # MANDATORY
            'clean_pdb'     : 'yes',
            'epsin'         : 2,
            'ionicstr'      : 0.1,
            'pbc_dimensions': 0,
            'temp'          : 310,
            'grid_fill'     : 0.8,         # FUTURE VERSION
            'ncpus'         : ncpus,
            'pH'            : '0,15',
            'pHstep'        : 0.25,
            'logfile'       : 'LOGFILE',
            'scaleM'        : 4,
            'scaleP'        : 1,
            'gsize'         : 81,
            'convergence'   : 0.01,
            'nlit'          : 300,
            'cutoff'        : -1,
            'relfac'        : 0.0,
            'output'        : 'pKas.out',
            'clean_pdb'     : True
        }
        sites = {' ': ('1N', '2C')}
        pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_onlytermini_noclean(self):
        from ..pypka import Titration
        os.system('rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent '
                  'contributions interactions.dat pkint LOG* __pycache__ *pyc')
        results = """
NTR 2001 8.167721738527177 ('protonated', 0.9363612426447844)
CTR 2002 3.390519247038917 ('deprotonated', 0.0002457041692844058)
        """
        results = results.split('\n')[1:-1]
        parameters = {
            'structure'     : 'ktp/ktp_pdb_onlytermini_noclean/ktp_noclean.pdb',     # MANDATORY
            'clean_pdb'     : 'no',
            'epsin'         : 2,
            'ionicstr'      : 0.1,
            'pbc_dimensions': 0,
            'temp'          : 310,
            'grid_fill'     : 0.8,         # FUTURE VERSION
            'ncpus'         : ncpus,
            'pH'            : '0,15',
            'pHstep'        : 0.25,
            'logfile'       : 'LOGFILE',
            'scaleM'        : 4,
            'scaleP'        : 1,
            'gsize'         : 81,
            'convergence'   : 0.01,
            'nlit'          : 300,
            'cutoff'        : -1,
            'relfac'        : 0.0,
            'output'        : 'pKas.out',
            'clean_pdb'     : False
        }
        sites = {' ': ('1N', '2C')}
        pKa = Titration(parameters, sites=sites)

        checkAPIResult(pKa, results)


class TestBuilder(object):
    def test_lyso_pH0(self):
        outfile = 'amber1.pdb'
        from ..pypka import Titration

        parameters = {'structure'     : 'builder/4lzt.pdb',
                      'epsin'         : 15,
                      'ionicstr'      : 0.1,
                      'pbc_dimensions': 0,
                      'ncpus'         : ncpus,
                      'output'        : 'pKas.out',
                      'titration_output': 'titration.out',
                      'structure_output': (outfile, 1, "amber")
                      }
        pKa = Titration(parameters)
        checkStructureOutput(outfile)

    def test_lyso_pH7(self):
        outfile = 'amber7.pdb'
        from ..pypka import Titration

        parameters = {'structure'     : 'builder/4lzt.pdb',
                      'epsin'         : 15,
                      'ionicstr'      : 0.1,
                      'pbc_dimensions': 0,
                      'ncpus'         : ncpus,
                      'output'        : 'pKas.out',
                      'titration_output': 'titration.out',
                      'structure_output': (outfile, 7, "amber")
                      }
        pKa = Titration(parameters)
        checkStructureOutput(outfile)

    def test_lyso_pH12(self):
        outfile = 'amber12.pdb'
        from ..pypka import Titration

        parameters = {'structure'     : 'builder/4lzt.pdb',
                      'epsin'         : 15,
                      'ionicstr'      : 0.1,
                      'pbc_dimensions': 0,
                      'ncpus'         : ncpus,
                      'output'        : 'pKas.out',
                      'titration_output': 'titration.out',
                      'structure_output': (outfile, 12, "amber")
                      }
        pKa = Titration(parameters)
        checkStructureOutput(outfile)

    def test_lyso_pH12_gromos(self):
        outfile = 'gromos12.pdb'
        from ..pypka import Titration

        parameters = {'structure'     : 'builder/4lzt.pdb',
                      'epsin'         : 15,
                      'ionicstr'      : 0.1,
                      'pbc_dimensions': 0,
                      'ncpus'         : ncpus,
                      'output'        : 'pKas.out',
                      'titration_output': 'titration.out',
                      'structure_output': (outfile, 12, "gromos_CpH")
                      }
        pKa = Titration(parameters)
        checkStructureOutput(outfile)
