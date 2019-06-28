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
             "python -m coverage erase; "
             "python -m coverage run ../../../pypka.py parameters.dat".format(path, ncpus), shell=True).wait()
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
1 NTR 7.90064610613
1 TYR 9.88522556058
2 CTR 3.12385071566
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_allsites(self):
        path = "ktp/ktp_pdb_allsites"
        results = """
1 NTR 7.90372759857
1 TYR 9.78309726325
2 CTR 3.30400367504
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_allsites_noclean(self):
        path = "ktp/ktp_pdb_allsites_noclean"
        results = """
1 NTR 7.90064610613
1 TYR 9.88522556058
2 CTR 3.12385071566
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_onlytermini(self):
        path = 'ktp/ktp_pdb_onlytermini'
        results = """
1 NTR 7.9001039321
2 CTR 3.29902344445
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_onlytermini_noclean(self):
        path = 'ktp/ktp_pdb_onlytermini_noclean'
        results = """
1 NTR 8.16771888442
2 CTR 3.39051924704
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_gro(self):
        path = "lyso/lyso_gro"
        results = """
1 NTR 7.63661126099
18 ASP 3.14160192118
35 GLU 4.6630113441
48 ASP 2.96396032156
66 ASP 3.25547061858
129 CTR 2.35320936885
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_sites(self):
        path = "lyso/lyso_pdb_sites"
        results = """
1 NTR 7.50547769612
18 ASP 3.1651915545
35 GLU 4.61083675072
48 ASP 2.25844185515
66 ASP 1.8713352369
129 CTR 2.3040367136
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_sites_noclean(self):
        path = "lyso/lyso_pdb_sites_noclean"
        results = """
1 NTR 7.39335708374
18 ASP 2.78161060181
35 GLU 4.27122389465
48 ASP 2.65449650965
66 ASP 3.01128613358
129 CTR 1.94092450986
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_all(self):
        path = "lyso/lyso_pdb_all"
        results = """
1 NTR 7.38414401309
1 LYS 10.3915151993
7 GLU 2.88092122222
13 LYS 100.0
15 HIS 5.86135647549
18 ASP 2.69457281787
20 TYR 9.91066416179
23 TYR 9.39233213737
24 SER 100.0
33 LYS 10.5103931957
35 GLU 4.17990458015
36 SER 100.0
40 THR 100.0
43 THR 100.0
47 THR 100.0
48 ASP 1.8672505767
50 SER 100.0
51 THR 100.0
52 ASP 2.14723044164
53 TYR 10.743180477
60 SER 100.0
66 ASP 2.34644743313
69 THR 100.0
72 SER 100.0
81 SER 100.0
85 SER 100.0
86 SER 100.0
87 ASP 2.04766398766
89 THR 100.0
91 SER 100.0
96 LYS 100.0
97 LYS 100.0
100 SER 100.0
101 ASP 3.31081873558
116 LYS 10.2551624733
118 THR 100.0
119 ASP 2.39551120758
129 CTR 1.70353465559
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_all_noclean(self):
        path = "lyso/lyso_pdb_all_noclean"
        results = """
1 NTR 7.38414401309
1 LYS 10.3915151993
7 GLU 2.88092122222
13 LYS 100.0
15 HIS 5.86135647549
18 ASP 2.69457281787
20 TYR 9.91066416179
23 TYR 9.39233213737
24 SER 100.0
33 LYS 10.5103931957
35 GLU 4.17990458015
36 SER 100.0
40 THR 100.0
43 THR 100.0
47 THR 100.0
48 ASP 1.8672505767
50 SER 100.0
51 THR 100.0
52 ASP 2.14723044164
53 TYR 10.743180477
60 SER 100.0
66 ASP 2.34644743313
69 THR 100.0
72 SER 100.0
81 SER 100.0
85 SER 100.0
86 SER 100.0
87 ASP 2.04766398766
89 THR 100.0
91 SER 100.0
96 LYS 100.0
97 LYS 100.0
100 SER 100.0
101 ASP 3.31081873558
116 LYS 10.2551624733
118 THR 100.0
119 ASP 2.39551120758
129 CTR 1.70353465559
        """
        runTest(path, ncpus, results)

    def test_cli_pHLIP_gro(self):
        path = "pHLIP/pHLIP_gro"
        results = """
769 NTR 10.767904983
770 CYS 100.0
771 GLU 3.99752739909
782 ASP 100.0
793 ASP 100.0
799 ASP 1.94131774795
801 ASP 0.990298982188
802 GLU 4.05994264808
804 CTR 3.84178774418
        """
        runTest(path, ncpus, results)

    def test_cli_pHLIP_pdb_all(self):
        path = "pHLIP/pHLIP_pdb_all"
        results = """
769 NTR 11.2025715977
770 CYS 24.2123881662
771 GLU 4.01078608539
776 TYR 18.7756943694
780 TYR 24.3484203892
782 ASP 100.0
786 THR 100.0
787 THR 100.0
793 ASP 22.6585921706
799 ASP 0.370645770177
801 ASP 1.48567974319
802 GLU 3.96561345595
804 CTR 4.2896319041
        """
        runTest(path, ncpus, results)


class TestAPI(object):
    def test_api_ktp_gro(self):
        from pypka import Titration
        os.system('rm -f LOG* *out *gro *pdb *pqr *crg *sites cent contributions interactions.dat pkint')
        results = """
NTR 7.90064610613 ('undefined', 0.8883318942980771)
1 9.88522556058 ('protonated', 0.9986992041484555)
CTR 3.12385071566 ('deprotonated', 0.00013298203005695505)
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
        from pypka import Titration
        os.system('rm -f *out *gro *pdb *pqr *crg *sites cent contributions interactions.dat pkint')
        results = """
NTR 7.90409283642 ('undefined', 0.8891167479866279)
1 9.78071352879 ('protonated', 0.9983458781425558)
CTR 3.30331361976 ('deprotonated', 0.00020101400276712804)
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
        from pypka import Titration
        os.system('rm -rf LOG* *out *gro *pdb *pqr *crg *sites '
                  'cent contributions interactions.dat pkint LOG* __pycache__ *pyc')
        results = """
NTR 7.90064610613 ('undefined', 0.8883318942980771)
1 9.88522556058 ('protonated', 0.9986992041484555)
CTR 3.12385071566 ('deprotonated', 0.00013298203005695505)
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
        from pypka import Titration
        os.system('rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent '
                  'contributions interactions.dat pkint LOG* __pycache__ *pyc')
        results = """
NTR 7.9001039321 ('undefined', 0.8882079948571108)
CTR 3.29902344445 ('deprotonated', 0.0001990384562145051)
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
        from pypka import Titration
        os.system('rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent '
                  'contributions interactions.dat pkint LOG* __pycache__ *pyc')
        results = """
NTR 8.16771888442 ('protonated', 0.9363608510357995)
CTR 3.39051924704 ('deprotonated', 0.0002457041692844058)
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
