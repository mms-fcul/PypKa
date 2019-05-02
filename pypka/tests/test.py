import pytest
import os
import subprocess as sb
import sys
sys.path.insert(1, '../')

ncpus = 8

def erase_old(directory):
    os.system("rm -f {0}/*out".format(directory))

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
        erase_old("ktp/ktp_gro")
        results = """
1 NTR 7.89621829987
1 TYR 9.88696575165
2 CTR 3.12329292297
        """
        results_lines = results.split('\n')[1:-1]
        sb.Popen('cd ktp/ktp_gro/; bash run.sh', shell=True).wait()
        checkOutput('ktp/ktp_gro/pKas.out', results_lines)

    def test_cli_ktp_pdb_allsites(self):
        erase_old("ktp/ktp_pdb_allsites")
        results = """
1 NTR 7.900203228
1 TYR 9.77924633026
2 CTR 3.30362272263
        """
        results_lines = results.split('\n')[1:-1]
        sb.Popen('cd ktp/ktp_pdb_allsites/; bash run.sh', shell=True).wait()
        checkOutput("ktp/ktp_pdb_allsites/pKas.out", results_lines)

    def test_cli_ktp_pdb_allsites_noclean(self):
        erase_old("ktp/ktp_pdb_allsites_noclean")
        results = """
1 NTR 7.89621829987
1 TYR 9.88696575165
2 CTR 3.12329292297
        """
        results_lines = results.split('\n')[1:-1]
        sb.Popen('cd ktp/ktp_pdb_allsites_noclean/; bash run_noclean.sh', shell=True).wait()
        checkOutput("ktp/ktp_pdb_allsites_noclean/pKas.out", results_lines)

    def test_cli_ktp_pdb_onlytermini(self):
        erase_old('ktp_pdb_onlytermini')
        results = """
1 NTR 7.904296875
2 CTR 3.29493141174
        """
        results_lines = results.split('\n')[1:-1]
        sb.Popen('cd ktp/ktp_pdb_onlytermini; bash run.sh', shell=True).wait()
        checkOutput("ktp/ktp_pdb_onlytermini/pKas.out", results_lines)

    def test_cli_ktp_pdb_onlytermini_noclean(self):
        erase_old('ktp_pdb_onlytermini_noclean')
        results = """
1 NTR 8.17041015625
2 CTR 3.38864183426
        """
        results_lines = results.split('\n')[1:-1]
        sb.Popen('cd ktp/ktp_pdb_onlytermini_noclean; bash run.sh', shell=True).wait()
        checkOutput("ktp/ktp_pdb_onlytermini_noclean/pKas.out", results_lines)

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
	    'structure'     : 'ktp/ktp_gro/TMP.gro',     # MANDATORY
	    'epsin'         : 2,
	    'ionicstr'      : 0.1,
	    'pbc_dimensions': 0,
	    'temp'          : 310,
	    'grid_fill'     : 0.8,         # FUTURE VERSION
	    'ncpus'         : 1,
	    'pH'            : '-5,15',
            'pHstep'        : 0.2,
	    'logfile'       : 'LOGFILE',
	    'scaleM'        : 4,
	    'scaleP'        : 1,
	    'gsizeM'        : 81,
	    'convergence'   : 0.01,
	    'nlit'          : 500,
	    'cutoff'        : -1,            
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
	    'ncpus'         : 1,
	    'pH'            : '0,15',
            'pHstep'        : 0.2,
	    'logfile'       : 'LOGFILE',
	    'scaleM'        : 4,
	    'scaleP'        : 1,
	    'gsizeM'        : 81,
	    'convergence'   : 0.01,
	    'nlit'          : 300,
	    'cutoff'        : -1,            
	    'output'        : 'pKas.out'
	}
	sites = 'all'		
	pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_allsites_noclean(self):
        from pypka import Titration
        os.system('rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent contributions interactions.dat pkint LOG* __pycache__ *pyc')
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
	    'ncpus'         : 1,
	    'pH'            : '-5,15',
            'pHstep'        : 0.2,
	    'logfile'       : 'LOGFILE',
	    'scaleM'        : 4,
	    'scaleP'        : 1,
	    'gsizeM'        : 81,
	    'convergence'   : 0.01,
	    'nlit'          : 500,
	    'cutoff'        : -1,            
	    'output'        : 'pKas.out'
	}
	sites = 'all'		
	pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_onlytermini(self):
        from pypka import Titration
        os.system('rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent contributions interactions.dat pkint LOG* __pycache__ *pyc')
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
	    'ncpus'         : 1,
	    'pH'            : '0,15',
            'pHstep'        : 0.25,
	    'logfile'       : 'LOGFILE',
	    'scaleM'        : 4,
	    'scaleP'        : 1,
	    'gsizeM'        : 81,
	    'convergence'   : 0.01,
	    'nlit'          : 300,
	    'cutoff'        : -1,            
	    'output'        : 'pKas.out'
	}
	sites = {'A': ('1N', '2C')}
	pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)
                    
    def test_api_ktp_pdb_onlytermini_noclean(self):
        from pypka import Titration
        os.system('rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent contributions interactions.dat pkint LOG* __pycache__ *pyc')
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
	    'ncpus'         : 1,
	    'pH'            : '0,15',
            'pHstep'        : 0.25,
	    'logfile'       : 'LOGFILE',
	    'scaleM'        : 4,
	    'scaleP'        : 1,
	    'gsizeM'        : 81,
	    'convergence'   : 0.01,
	    'nlit'          : 300,
	    'cutoff'        : -1,            
	    'output'        : 'pKas.out'
	}
	sites = {'A': ('1N', '2C')}
	pKa = Titration(parameters, sites=sites)

        checkAPIResult(pKa, results)
                    


