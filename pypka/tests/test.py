import pytest
import os
import subprocess as sb
import sys
from pypka.tests.builder.check_diffs import compareFiles

sys.path.insert(1, '../')

ncpus = 2

# run coverage.sh to generate coverage

def runTest(path, ncpus, results):
    os.system("rm -f {0}/*out {0}/clean*pqr {0}/TMP.gro {0}/delphi_in_*pdb".format(path))
    results_lines = results.split('\n')[1:-1]
    sb.Popen("cd {0}; "
             "sed -i 's/ncpus .*/ncpus           = {1}/' parameters.dat; "
             "python3 -m coverage erase; "
             "python3 -m coverage run ../../../pypka.py parameters.dat;"
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
        result = '{0} {1} {2}'.format(site, pKa[site], pKa.getProtState(site, 7))
        i += 1
        assert result == results[i]
    assert i + 1 == len(results)

class TestCLI(object):
    def test_cli_ktp_gro(self):
        path = "ktp/ktp_gro"
        results = """
1 NTR 7.900646106133701
1 TYR 9.88522556057786
2 CTR 3.123850715655238
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_allsites(self):
        path = "ktp/ktp_pdb_allsites"
        results = """
1 NTR 7.9037275985663085
1 TYR 9.783097263247802
2 CTR 3.3040036750415207
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_allsites_noclean(self):
        path = "ktp/ktp_pdb_allsites_noclean"
        results = """
1 NTR 7.900646106133701
1 TYR 9.88522556057786
2 CTR 3.123850715655238
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_onlytermini(self):
        path = 'ktp/ktp_pdb_onlytermini'
        results = """
1 NTR 7.9001039320976965
2 CTR 3.2990234444483972
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_onlytermini_noclean(self):
        path = 'ktp/ktp_pdb_onlytermini_noclean'
        results = """
1 NTR 8.167718884417926
2 CTR 3.390519247038917
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_gro(self):
        path = "lyso/lyso_gro"
        results = """
1 NTR 7.63661126098686
18 ASP 3.141601921175308
35 GLU 4.663011344104504
48 ASP 2.9639603215645716
66 ASP 3.2554706185808753
129 CTR 2.353209368853925
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_sites(self):
        path = "lyso/lyso_pdb_sites"
        results = """
1 NTR 7.505477696120323
18 ASP 3.1651915544951166
35 GLU 4.610836750719833
48 ASP 2.258441855151931
66 ASP 1.871335236896869
129 CTR 2.3040367135982893
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_sites_noclean(self):
        path = "lyso/lyso_pdb_sites_noclean"
        results = """
1 NTR 7.39335708373817
18 ASP 2.7816106018108915
35 GLU 4.271223894654857
48 ASP 2.6544965096521933
66 ASP 3.0112861335754557
129 CTR 1.9409245098602885
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_all(self):
        path = "lyso/lyso_pdb_all"
        results = """
1 NTR 7.384144013092099
1 LYS 10.391515199306076
7 GLU 2.880921222215613
13 LYS 100.0
15 HIS 5.861356475494928
18 ASP 2.6945728178749397
20 TYR 9.910664161785833
23 TYR 9.392332137365454
24 SER 100.0
33 LYS 10.510393195704687
35 GLU 4.179904580152671
36 SER 100.0
40 THR 100.0
43 THR 100.0
47 THR 100.0
48 ASP 1.8672505767012688
50 SER 100.0
51 THR 100.0
52 ASP 2.1472304416360166
53 TYR 10.743180477033484
60 SER 100.0
66 ASP 2.3464474331287817
69 THR 100.0
72 SER 100.0
81 SER 100.0
85 SER 100.0
86 SER 100.0
87 ASP 2.0476639876603344
89 THR 100.0
91 SER 100.0
96 LYS 100.0
97 LYS 100.0
100 SER 100.0
101 ASP 3.31081873557565
116 LYS 10.255162473280647
118 THR 100.0
119 ASP 2.3955112075754372
129 CTR 1.7035346555851842
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_all_noclean(self):
        path = "lyso/lyso_pdb_all_noclean"
        results = """
1 NTR 7.384144013092099
1 LYS 10.391515199306076
7 GLU 2.880921222215613
13 LYS 100.0
15 HIS 5.861356475494928
18 ASP 2.6945728178749397
20 TYR 9.910664161785833
23 TYR 9.392332137365454
24 SER 100.0
33 LYS 10.510393195704687
35 GLU 4.179904580152671
36 SER 100.0
40 THR 100.0
43 THR 100.0
47 THR 100.0
48 ASP 1.8672505767012688
50 SER 100.0
51 THR 100.0
52 ASP 2.1472304416360166
53 TYR 10.743180477033484
60 SER 100.0
66 ASP 2.3464474331287817
69 THR 100.0
72 SER 100.0
81 SER 100.0
85 SER 100.0
86 SER 100.0
87 ASP 2.0476639876603344
89 THR 100.0
91 SER 100.0
96 LYS 100.0
97 LYS 100.0
100 SER 100.0
101 ASP 3.31081873557565
116 LYS 10.255162473280647
118 THR 100.0
119 ASP 2.3955112075754372
129 CTR 1.7035346555851842
        """
        runTest(path, ncpus, results)

    def test_cli_pHLIP_gro(self):
        path = "pHLIP/pHLIP_gro"
        results = """
769 NTR 10.767904983032684
770 CYS 100.0
771 GLU 3.997527399091152
782 ASP 100.0
793 ASP 100.0
799 ASP 1.9413177479539034
801 ASP 0.9902989821882952
802 GLU 4.059942648076733
804 CTR 3.8417877441756705
        """
        runTest(path, ncpus, results)

    def test_cli_pHLIP_pdb_all(self):
        path = "pHLIP/pHLIP_pdb_all"
        results = """
769 NTR 11.20670245619075
770 CYS 24.20614026887048
771 GLU 4.006259451388756
776 TYR 18.773742141263668
780 TYR 24.34882079633006
782 ASP 100.0
786 THR 100.0
787 THR 100.0
793 ASP 22.671456552380604
799 ASP 0.4467670504871568
801 ASP 1.4713225156432284
802 GLU 3.954407177042764
804 CTR 4.7834128769177315
804 THR 100.0
        """
        runTest(path, ncpus, results)

    def test_cli_nucleosome_pdb_all(self):
        path = "nucleosome/nucleosome_pdb_all"
        results = """
1 NTR 7.889916847060441
3 THR 100.0
4 LYS 10.28997470844457
6 THR 100.0
9 LYS 8.464272616136919
10 SER 100.0
11 THR 100.0
14 LYS 10.058901584961298
18 LYS 7.657196150519031
22 THR 100.0
23 LYS 9.677646644354413
27 LYS 10.839003677114341
28 SER 100.0
32 THR 100.0
36 LYS 10.421760695545597
37 LYS 11.257852058413444
39 HIS 8.75823395797842
41 TYR 13.695628162189125
45 THR 100.0
50 GLU 3.08273635080864
54 TYR 10.820982492349298
56 LYS 10.450701578551742
57 SER 100.0
58 THR 100.0
59 GLU 2.948171586494172
64 LYS 10.751162708759072
73 GLU 3.4186971463876126
77 ASP 4.4852618344376145
79 LYS 10.886947935368044
80 THR 100.0
81 ASP 2.576220362622036
86 SER 100.0
87 SER 100.0
94 GLU 2.330335757057314
96 SER 100.0
97 GLU 0.2605436274114134
99 TYR 12.73404896055955
105 GLU 4.157654966392831
106 ASP 2.3647680694777
107 THR 100.0
110 CYS 13.245581710160186
113 HIS 3.458425170817721
115 LYS 10.41972280689237
118 THR 100.0
122 LYS 8.980250289771082
123 ASP 100.0
133 GLU 0.6664895285221452
135 CTR 3.5103777044371105
"""
        runTest(path, ncpus, results)

                
class TestAPI(object):
    def test_api_ktp_gro(self):
        from ..pypka import Titration
        os.system('rm -f LOG* *out *gro *pdb *pqr *crg *sites cent contributions interactions.dat pkint')
        results = """
NTR 7.900646106133701 ('undefined', 0.8883318942980771)
1 9.88522556057786 ('protonated', 0.9986992041484555)
CTR 3.123850715655238 ('deprotonated', 0.00013298203005695505)
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
            'output'        : 'pKas.out'
        }
        sites = {'A': ('1N', '1', '2C')}
        pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_allsites_clean(self):
        from ..pypka import Titration
        os.system('rm -f *out *gro *pdb *pqr *crg *sites cent contributions interactions.dat pkint')
        results = """
NTR 7.904092836420442 ('undefined', 0.8891167479866279)
1 9.780713528788414 ('protonated', 0.9983458781425558)
CTR 3.3033136197599586 ('deprotonated', 0.00020101400276712804)
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
            'output'        : 'pKas.out'
        }
        sites = 'all'           
        pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_allsites_noclean(self):
        from ..pypka import Titration
        os.system('rm -rf LOG* *out *gro *pdb *pqr *crg *sites '
                  'cent contributions interactions.dat pkint LOG* __pycache__ *pyc')
        results = """
NTR 7.900646106133701 ('undefined', 0.8883318942980771)
1 9.88522556057786 ('protonated', 0.9986992041484555)
CTR 3.123850715655238 ('deprotonated', 0.00013298203005695505)
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
            'output'        : 'pKas.out'
        }
        sites = 'all'           
        pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_onlytermini(self):
        from ..pypka import Titration
        os.system('rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent '
                  'contributions interactions.dat pkint LOG* __pycache__ *pyc')
        results = """
NTR 7.9001039320976965 ('undefined', 0.8882079948571108)
CTR 3.2990234444483972 ('deprotonated', 0.0001990384562145051)
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
            'output'        : 'pKas.out'
        }
        sites = {'A': ('1N', '2C')}
        pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_onlytermini_noclean(self):
        from ..pypka import Titration
        os.system('rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent '
                  'contributions interactions.dat pkint LOG* __pycache__ *pyc')
        results = """
NTR 8.167718884417926 ('protonated', 0.9363608510357995)
CTR 3.390519247038917 ('deprotonated', 0.0002457041692844058)
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
            'output'        : 'pKas.out'
        }
        sites = {'A': ('1N', '2C')}
        pKa = Titration(parameters, sites=sites)

        checkAPIResult(pKa, results)

class TestBuilder(object):
    def test_lyso_pH1(self):
        outfile = 'gromos1.pdb'
        from ..pypka import Titration

        parameters = {'structure'     : 'builder/4lzt.pdb',
                      'epsin'         : 15,
                      'ionicstr'      : 0.1,
                      'pbc_dimensions': 0,
                      'ncpus'         : ncpus,
                      'output'        : 'pKas.out',
                      'titration_output': 'titration.out',
                      'structure_output': (outfile, 1)
                      }
        sites = {'A': ('1N', '2C')} # overwritten
        pKa = Titration(parameters, sites=sites)
        checkStructureOutput(outfile)

    def test_lyso_pH7(self):
        outfile = 'gromos7.pdb'
        from ..pypka import Titration

        parameters = {'structure'     : 'builder/4lzt.pdb',
                      'epsin'         : 15,
                      'ionicstr'      : 0.1,
                      'pbc_dimensions': 0,
                      'ncpus'         : ncpus,
                      'output'        : 'pKas.out',
                      'titration_output': 'titration.out',
                      'structure_output': (outfile, 7)
                      }
        pKa = Titration(parameters)
        checkStructureOutput(outfile)

    def test_lyso_pH12(self):
        outfile = 'gromos12.pdb'
        from ..pypka import Titration
        
        parameters = {'structure'     : 'builder/4lzt.pdb',
                      'epsin'         : 15,
                      'ionicstr'      : 0.1,
                      'pbc_dimensions': 0,
                      'ncpus'         : ncpus,
                      'output'        : 'pKas.out',
                      'titration_output': 'titration.out',
                      'structure_output': (outfile, 12)
                      }
        pKa = Titration(parameters)
        checkStructureOutput(outfile)
