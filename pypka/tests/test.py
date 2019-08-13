import pytest
import os
import subprocess as sb
import sys
sys.path.insert(1, '../')

ncpus = 1

# run coverage.sh to generate coverage

def runTest(path, ncpus, results):
    os.system("rm -f {0}/*out {0}/clean*pqr {0}/TMP.gro {0}/delphi_in_*pdb".format(path))
    results_lines = results.split('\n')[1:-1]
    sb.Popen("cd {0}; "
             "sed -i 's/ncpus .*/ncpus           = {1}/' parameters.dat; "
             "python3 -m coverage erase; "
             "python3 -m coverage run ../../../pypka.py parameters.dat".format(path, ncpus), shell=True).wait()
    checkOutput('{0}/pKas.out'.format(path), results_lines)    

def checkOutput(filename, results_lines):
    with open(filename) as f:
        c = -1
        for line in f:
            c += 1
            line = line.strip()
            assert line == results_lines[c]
    assert c + 1 == len(results_lines)

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
769 NTR 11.202571597654359
770 CYS 24.21238816620012
771 GLU 4.010786085385137
776 TYR 18.77569436939371
780 TYR 24.348420389190892
782 ASP 100.0
786 THR 100.0
787 THR 100.0
793 ASP 22.658592170556883
799 ASP 0.37064577017674316
801 ASP 1.4856797431942301
802 GLU 3.965613455945925
804 CTR 4.289631904096492
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
