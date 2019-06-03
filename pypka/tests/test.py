import pytest
import os
import subprocess as sb
import sys
sys.path.insert(1, '../')

ncpus = 2

# run coverage.sh to generate coverage

def runTest(path, ncpus, results):
    os.system("rm -f {0}/*out".format(path))
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
1 NTR 7.89621829987
1 TYR 9.88696575165
2 CTR 3.12329292297
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_allsites(self):
        path = "ktp/ktp_pdb_allsites"
        results = """
1 NTR 7.900203228
1 TYR 9.77924633026
2 CTR 3.30362272263
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_allsites_noclean(self):
        path = "ktp/ktp_pdb_allsites_noclean"
        results = """
1 NTR 7.89621829987
1 TYR 9.88696575165
2 CTR 3.12329292297
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_onlytermini(self):
        path = 'ktp/ktp_pdb_onlytermini'
        results = """
1 NTR 7.904296875
2 CTR 3.29493141174
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_onlytermini_noclean(self):
        path = 'ktp/ktp_pdb_onlytermini_noclean'
        results = """
1 NTR 8.17041015625
2 CTR 3.38864183426
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_gro(self):
        path = "lyso/lyso_gro"
        results = """
1 NTR 7.63068819046
18 ASP 3.14008665085
35 GLU 4.66894435883
48 ASP 2.96325278282
66 ASP 3.25400090218
129 CTR 2.35684490204
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_sites(self):
        path = "lyso/lyso_pdb_sites"
        results = """
1 NTR 7.50067424774
18 ASP 3.16751217842
35 GLU 4.61393737793
48 ASP 2.26052093506
66 ASP 1.87332427502
129 CTR 2.30435156822
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_sites_noclean(self):
        path = "lyso/lyso_pdb_sites_noclean"
        results = """
1 NTR 7.39441871643
18 ASP 2.78099513054
35 GLU 4.27633476257
48 ASP 2.65356683731
66 ASP 3.01925420761
129 CTR 1.94424641132
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_all(self):
        path = "lyso/lyso_pdb_all"
        results = """
1 NTR 7.38174772263
1 LYS 10.3868246078
7 GLU 2.88459134102
13 LYS 100.0
15 HIS 5.86267566681
18 ASP 2.69885849953
20 TYR 9.90721607208
23 TYR 9.38904476166
24 SER 100.0
33 LYS 10.5100259781
35 GLU 4.17546081543
36 SER 100.0
40 THR 100.0
43 THR 100.0
47 THR 100.0
48 ASP 1.86072945595
50 SER 100.0
51 THR 100.0
52 ASP 2.15181279182
53 TYR 10.7421417236
60 SER 100.0
66 ASP 2.34874129295
69 THR 100.0
72 SER 100.0
81 SER 100.0
85 SER 100.0
86 SER 100.0
87 ASP 2.05946564674
89 THR 100.0
91 SER 100.0
96 LYS 100.0
97 LYS 100.0
100 SER 100.0
101 ASP 3.31438612938
116 LYS 10.2509069443
118 THR 100.0
119 ASP 2.39909434319
129 CTR 1.69917798042
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_all_noclean(self):
        path = "lyso/lyso_pdb_all_noclean"
        results = """
1 NTR 7.38174772263
1 LYS 10.3868246078
7 GLU 2.88459134102
13 LYS 100.0
15 HIS 5.86267566681
18 ASP 2.69885849953
20 TYR 9.90721607208
23 TYR 9.38904476166
24 SER 100.0
33 LYS 10.5100259781
35 GLU 4.17546081543
36 SER 100.0
40 THR 100.0
43 THR 100.0
47 THR 100.0
48 ASP 1.86072945595
50 SER 100.0
51 THR 100.0
52 ASP 2.15181279182
53 TYR 10.7421417236
60 SER 100.0
66 ASP 2.34874129295
69 THR 100.0
72 SER 100.0
81 SER 100.0
85 SER 100.0
86 SER 100.0
87 ASP 2.05946564674
89 THR 100.0
91 SER 100.0
96 LYS 100.0
97 LYS 100.0
100 SER 100.0
101 ASP 3.31438612938
116 LYS 10.2509069443
118 THR 100.0
119 ASP 2.39909434319
129 CTR 1.69917798042
        """
        runTest(path, ncpus, results)

    def test_cli_pHLIP_gro(self):
        path = "pHLIP/pHLIP_gro"
        results = """
769 NTR 11.0565090179
770 CYS 100.0
771 GLU 4.03883981705
782 ASP 100.0
793 ASP 100.0
799 ASP 1.80626475811
801 ASP 0.836304783821
802 GLU 3.98396849632
804 CTR 3.63486742973
"""
        runTest(path, ncpus, results)

    def test_cli_pHLIP_pdb_all(self):
        path = "pHLIP/pHLIP_pdb_all"
        results = """
769 NTR 11.2083435059
770 CYS 24.2029476166
771 GLU 4.00404262543
776 TYR 18.7731571198
780 TYR 24.3460330963
782 ASP 100.0
786 THR 100.0
787 THR 100.0
793 ASP 22.6629772186
799 ASP 0.369246691465
801 ASP 1.48556101322
802 GLU 3.96952009201
804 CTR 4.29328298569
"""
        runTest(path, ncpus, results)


class TestAPI(object):
    def test_api_ktp_gro(self):
        from pypka import Titration
        os.system('rm -f LOG* *out *gro *pdb *pqr *crg *sites cent contributions interactions.dat pkint')
        results = """
NTR 7.89621829987 ('undefined', 0.8873165161846842)
1 9.88696575165 ('protonated', 0.9987043991888108)
CTR 3.12329292297 ('deprotonated', 0.00013281136488231348)
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
NTR 7.90208864212 ('undefined', 0.888660962962448)
1 9.78845977783 ('protonated', 0.9983750726401587)
CTR 3.3024392128 ('deprotonated', 0.00020060977015866723)
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
NTR 7.89621829987 ('undefined', 0.8873165161846842)
1 9.88696575165 ('protonated', 0.9987043991888108)
CTR 3.12329292297 ('deprotonated', 0.00013281136488231348)
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
NTR 7.904296875 ('undefined', 0.8891630578319204)
CTR 3.29493141174 ('deprotonated', 0.0001971722409805938)
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
NTR 8.17041015625 ('protonated', 0.9367291213385088)
CTR 3.38864183426 ('deprotonated', 0.00024464456585429595)
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
