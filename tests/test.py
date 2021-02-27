import pytest
import os
import subprocess as sb
import sys
from .builder.check_diffs import compareFiles

sys.path.insert(1, "../")
from pypka.pypka import Titration

ncpus = -1


def runTest(path, ncpus, results, delete_extra=""):
    os.system(
        "rm -f {0}/*out {0}/clean*pqr {0}/TMP.gro {0}/delphi_in_*pdb {1}".format(
            path, delete_extra
        )
    )
    results_lines = results.split("\n")[1:-1]
    # "python3 -m coverage erase; "
    # "python3 -m coverage run ../../../pypka.py parameters.dat;"
    sb.Popen(
        "cd {0}; "
        "sed -i 's/ncpus .*/ncpus           = {1}/' parameters.dat; "
        "python3 ../../../pypka/pypka.py parameters.dat;"
        "rm -f *-*.pdb *.profl *.prm *.xvg *.frc *.crg cent "
        "contributions interactions.dat pkint tmp.sites".format(path, ncpus),
        shell=True,
    ).wait()
    checkOutput("{0}/pKas.out".format(path), results_lines)


def checkOutput(filename, results_lines):
    with open(filename) as f:
        c = -1
        for line in f:
            c += 1
            line = line.strip()
            print(line)
            assert line == results_lines[c]
    assert c + 1 == len(results_lines)


def checkStructureOutput(filename):
    problems = compareFiles("builder/{0}".format(filename), filename)
    if problems:
        raise Exception("Problems found with {0}".format(filename))
    else:
        os.remove(filename)


def check_file_diff(f1, f2):

    with open(f1) as f:
        content1 = f.read()

    with open(f2) as f:
        content2 = f.read()

    assert content1 == content2


def checkAPIResult(pKa, results):
    i = -1
    for site in pKa:
        result = "{0} {1} {2} {3}".format(
            site.res_name, site.res_number, site.pK, site.getProtState(7)
        )
        i += 1
        assert result == results[i]
    assert i + 1 == len(results)


class TestCLI(object):
    def test_cli_ktp_gro(self):
        path = "ktp/ktp_gro"
        results = """
1 NTR 7.932515230635335    A
1 TYR 9.884779718110096    A
2 CTR 2.8373655675051106   A
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_allsites(self):
        path = "ktp/ktp_pdb_allsites"
        results = """
1 NTR 7.934083369392759
1 TYR 9.783344629543752
2 CTR 3.010420419338569
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_allsites_noclean(self):
        path = "ktp/ktp_pdb_allsites_noclean"
        results = """
1 NTR 7.932515230635335
1 TYR 9.884779718110096
2 CTR 2.8373655675051106
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_onlytermini(self):
        path = "ktp/ktp_pdb_onlytermini"
        results = """
1 NTR 7.931056701030927
2 CTR 3.0073796144818266
        """
        runTest(path, ncpus, results)

    def test_cli_ktp_pdb_onlytermini_noclean(self):
        path = "ktp/ktp_pdb_onlytermini_noclean"
        results = """
1 NTR 8.199054937290231
2 CTR 3.099830508474576
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_gro(self):
        path = "lyso/lyso_gro"
        results = """
1 NTR 7.668774975293259    A
18 ASP 3.2083443047169387   A
35 GLU 4.6773087182961435   A
48 ASP 3.017856790123458    A
66 ASP 3.312043638451601    A
129 CTR 2.056                A
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_sites(self):
        path = "lyso/lyso_pdb_sites"
        results = """
1 NTR 7.5351439903431645
18 ASP 3.255433112350178
35 GLU 4.619474578486473
48 ASP 2.3134230246724226
66 ASP 1.891644424338385
129 CTR 2.00831923042665
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_sites_noclean(self):
        path = "lyso/lyso_pdb_sites_noclean"
        results = """
1 NTR 7.423739542051959
18 ASP 2.8497022976277133
35 GLU 4.278395263521122
48 ASP 2.7135698073794363
66 ASP 3.0788890257420864
129 CTR 1.6449628904042497
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_all(self):
        path = "lyso/lyso_pdb_all"
        results = """
1 NTR 7.414526765648602
1 LYS 10.36729217344299
7 GLU 2.894038631238357
13 LYS None
15 HIS 5.726972226653197
18 ASP 2.7716275732876303
20 TYR 9.911490578356108
23 TYR 9.396414012870652
24 SER None
33 LYS 10.487372183063636
35 GLU 4.185039523383046
36 SER None
40 THR None
43 THR None
47 THR None
48 ASP 1.9185078450208133
50 SER None
51 THR None
52 ASP 2.2085225009956195
53 TYR 10.740697784064121
60 SER None
66 ASP 2.403900927381685
69 THR None
72 SER None
81 SER None
85 SER None
86 SER None
87 ASP 2.1193747587803937
89 THR None
91 SER None
96 LYS None
97 LYS None
100 SER None
101 ASP 3.3738081246564176
116 LYS 10.226543135843928
118 THR None
119 ASP 2.45892517060784
129 CTR 1.3758302160470264
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_all_noclean(self):
        path = "lyso/lyso_pdb_all_noclean"
        results = """
1 NTR 7.410777679165742
1 LYS 10.370347660204002
7 GLU 2.8965110938030474
13 LYS None
15 HIS 5.726922173858075
18 ASP 2.7683778888070547
20 TYR 9.915309193926785
23 TYR 9.386319106388996
24 SER None
33 LYS 10.489262258753374
35 GLU 4.181233181233181
36 SER None
40 THR None
43 THR None
47 THR None
48 ASP 1.9239332323182952
50 SER None
51 THR None
52 ASP 2.2139391755149567
53 TYR 10.743145797340945
60 SER None
66 ASP 2.4136309313031163
69 THR None
72 SER None
81 SER None
85 SER None
86 SER None
87 ASP 2.1255277147994684
89 THR None
91 SER None
96 LYS None
97 LYS None
100 SER None
101 ASP 3.3811188811188813
116 LYS 10.241441458807202
118 THR None
119 ASP 2.4695061553700297
129 CTR 1.373601454300118
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_multichain_all(self):
        path = "lyso/lyso_multichain_all"
        results = """
1 NTR 7.258148916332084    A
1 LYS 10.649226123820188   A
7 GLU 2.894240130337168    A
12 CTR None                 A
13 NTR 5.320197570158755    B
13 LYS None                 B
15 HIS 6.919149690360516    B
18 ASP 2.5703887911687997   B
20 CTR None                 B
20 TYR 10.15103009061378    B
21 NTR None                 C
23 TYR 8.597807220550798    C
24 SER None                 C
33 LYS 10.619279181732052   C
35 GLU 4.248735012108444    C
36 SER None                 C
40 CTR 1.112605534954853    C
40 THR None                 C
41 NTR 9.446173994761361    D
43 THR None                 D
47 THR None                 D
48 ASP 1.9209061407888166   D
50 SER None                 D
51 THR None                 D
52 ASP 2.2624183934147033   D
53 TYR 10.886784944775512   D
60 SER None                 D
66 ASP 2.384563181468579    D
69 THR None                 D
72 SER None                 D
81 SER None                 D
85 SER None                 D
86 SER None                 D
87 ASP 2.2879831342234715   D
89 THR None                 D
91 SER None                 D
96 LYS None                 D
97 LYS None                 D
100 SER None                 D
101 ASP 3.445209202990001    D
116 LYS 10.246746860341958   D
118 THR None                 D
119 ASP 2.4371008876850992   D
129 CTR 1.6374364742984553   D
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_multichain_sites(self):
        path = "lyso/lyso_multichain_sites"
        results = """
1 NTR 7.242283991665559    A
1 ASP 3.115283494887798    B
35 GLU 4.622119696447052    C
35 ASP 2.3149351529666955   D
66 ASP 1.8855224015144685   D
1 CTR 2.395771670190275    D
        """
        runTest(path, ncpus, results)

    def test_cli_pHLIP_gro(self):
        path = "pHLIP/pHLIP_gro"
        results = """
769 NTR 10.792369232976968   A
770 CYS None                 A
771 GLU 4.010003260890415    A
782 ASP None                 A
793 ASP None                 A
799 ASP 2.0326332157953524   A
801 ASP 1.0540321754278101   A
802 GLU 4.079350871731009    A
804 CTR 3.5223596867999754   A
        """
        runTest(path, ncpus, results)

    def test_cli_pHLIP_pdb_all(self):
        path = "pHLIP/pHLIP_pdb_all"
        results = """
769 NTR 11.23426435877262
770 CYS 24.30125978205222
771 GLU 4.017264339389327
776 TYR 18.59752245862884
780 TYR 24.40040678638754
782 ASP None
786 THR None
787 THR None
793 ASP 22.780211657381866
799 ASP 0.48126431060228975
801 ASP 1.5060378425314513
802 GLU 3.9599154489122252
804 CTR 4.475875871985698
804 THR None
        """
        runTest(path, ncpus, results)

    def test_cli_nucleosome_pdb_all(self):
        path = "nucleic_acids/nucleosome"
        results = """
1 NTR 9.189724382638616    A
4 LYS 10.318058252427184   A
9 LYS 8.413024233705876    A
14 LYS 10.187375962864543   A
18 LYS 7.836933797909408    A
23 LYS 9.685264765422453    A
27 LYS 10.48775419193721    A
36 LYS 10.308033923034834   A
37 LYS 10.201180692209746   A
39 HIS 6.135722229945781    A
41 TYR 9.245076801890509    A
50 GLU 2.607285466040596    A
54 TYR 10.42479042961928    A
56 LYS 10.321385265478597   A
59 GLU 2.4567940354147253   A
64 LYS 9.999453522539955    A
73 GLU 3.560929696504034    A
77 ASP 4.62020624303233     A
79 LYS 10.673193776123073   A
81 ASP 2.242447157106334    A
94 GLU 4.89061127700685     A
97 GLU 1.5794541347257711   A
99 TYR 13.09163800968489    A
105 GLU 4.697366553094833    A
106 ASP 3.2938242668701623   A
110 CYS 10.885418159667479   A
113 HIS 6.689373344486268    A
115 LYS 10.516852808498253   A
122 LYS 10.044322818009617   A
123 ASP 0.1261840913839939   A
133 GLU 1.8462956429705415   A
135 CTR 3.471648476661097    A
"""
        runTest(path, ncpus, results, delete_extra="nucleosome_4.pdb")
        check_file_diff("builder/nucleosome_4.pdb", f"{path}/nucleosome_4.pdb")

    def test_cli_crispr_pdb_all(self):
        path = "nucleic_acids/crispr"
        results = """
3 NTR 8.005239278383925    B
3 LYS 11.142589338019919   B
4 LYS 10.640888208269525   B
5 TYR 9.671807437301815    B
10 ASP 3.879042759106106    B
23 ASP 3.3729027636391056   B
24 GLU 3.8071811109808062   B
25 TYR 9.92113539726605     B
26 LYS 11.204305457779094   B
30 LYS 11.640715174833089   B
31 CTR 3.8620041827655878   B
31 LYS 12.08578953064685    B
"""
        runTest(path, ncpus, results, delete_extra="crispr_10.pdb")
        check_file_diff("builder/crispr_10.pdb", f"{path}/crispr_10.pdb")


class TestAPI(object):
    def test_api_ktp_gro(self):
        os.system(
            "rm -f LOG* *out *gro *pdb *pqr *crg *sites cent contributions interactions.dat pkint"
        )
        results = """
NTR 2001 7.932515230635335 ('undefined', 0.8954064318344567)
TYR 1 9.884779718110096 ('protonated', 0.9986978698182121)
CTR 2002 2.8373655675051106 ('deprotonated', 6.875997409570579e-05)
        """
        results = results.split("\n")[1:-1]
        parameters = {
            "structure": "ktp/ktp_gro/ktp.gro",  # MANDATORY
            "epsin": 2,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "temp": 310,
            "grid_fill": 0.8,  # FUTURE VERSION
            "ncpus": ncpus,
            "pH": "-5,15",
            "pHstep": 0.2,
            "logfile": "LOGFILE",
            "scaleM": 4,
            "scaleP": 1,
            "gsize": 81,
            "convergence": 0.01,
            "nlit": 500,
            "cutoff": -1,
            "relfac": 0.0,
            "output": "pKas.out",
            "clean_pdb": False,
        }
        sites = {"A": ("1N", "1", "2C")}
        pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_allsites_clean(self):

        os.system(
            "rm -f *out *gro *pdb *pqr *crg *sites cent contributions interactions.dat pkint"
        )
        results = """
NTR 2001 7.934819040827242 ('undefined', 0.8959021976401332)
TYR 1 9.78114442380826 ('protonated', 0.9983475157921265)
CTR 2002 3.010253456221198 ('deprotonated', 0.00010237855405755827)
        """
        results = results.split("\n")[1:-1]
        parameters = {
            "structure": "ktp/ktp_pdb_allsites/ktp.pdb",
            "clean_pdb": "yes",
            "epsin": 2,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "temp": 310,
            "grid_fill": 0.8,
            "ncpus": ncpus,
            "pH": "0,15",
            "pHstep": 0.2,
            "logfile": "LOGFILE",
            "scaleM": 4,
            "scaleP": 1,
            "gsize": 81,
            "convergence": 0.01,
            "nlit": 300,
            "cutoff": -1,
            "relfac": 0.0,
            "output": "pKas.out",
            "clean_pdb": True,
        }
        sites = "all"
        pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_allsites_noclean(self):

        os.system(
            "rm -rf LOG* *out *gro *pdb *pqr *crg *sites "
            "cent contributions interactions.dat pkint LOG* __pycache__ *pyc"
        )
        results = """
NTR 2001 7.932515230635335 ('undefined', 0.8954064318344567)
TYR 1 9.884779718110096 ('protonated', 0.9986978698182121)
CTR 2002 2.8373655675051106 ('deprotonated', 6.875997409570579e-05)
        """
        results = results.split("\n")[1:-1]
        parameters = {
            "structure": "ktp/ktp_pdb_allsites_noclean/ktp_noclean.pdb",  # MANDATORY
            "clean_pdb": "no",
            "epsin": 2,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "temp": 310,
            "grid_fill": 0.8,  # FUTURE VERSION
            "ncpus": ncpus,
            "pH": "-5,15",
            "pHstep": 0.2,
            "logfile": "LOGFILE",
            "scaleM": 4,
            "scaleP": 1,
            "gsize": 81,
            "convergence": 0.01,
            "nlit": 500,
            "cutoff": -1,
            "relfac": 0.0,
            "output": "pKas.out",
            "clean_pdb": False,
        }
        sites = "all"
        pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_onlytermini(self):

        os.system(
            "rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent "
            "contributions interactions.dat pkint LOG* __pycache__ *pyc"
        )
        results = """
NTR 2001 7.931056701030927 ('undefined', 0.8950914882162796)
CTR 2002 3.0073796144818266 ('deprotonated', 0.00010170339324314826)
        """
        results = results.split("\n")[1:-1]
        parameters = {
            "structure": "ktp/ktp_pdb_onlytermini/ktp.pdb",  # MANDATORY
            "clean_pdb": "yes",
            "epsin": 2,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "temp": 310,
            "grid_fill": 0.8,  # FUTURE VERSION
            "ncpus": ncpus,
            "pH": "0,15",
            "pHstep": 0.25,
            "logfile": "LOGFILE",
            "scaleM": 4,
            "scaleP": 1,
            "gsize": 81,
            "convergence": 0.01,
            "nlit": 300,
            "cutoff": -1,
            "relfac": 0.0,
            "output": "pKas.out",
            "clean_pdb": True,
        }
        sites = {" ": ("1N", "2C")}
        pKa = Titration(parameters, sites=sites)
        checkAPIResult(pKa, results)

    def test_api_ktp_pdb_onlytermini_noclean(self):

        os.system(
            "rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent "
            "contributions interactions.dat pkint LOG* __pycache__ *pyc"
        )
        results = """
NTR 2001 8.199054937290231 ('protonated', 0.9405274528487423)
CTR 2002 3.099830508474576 ('deprotonated', 0.00012582758427877238)
        """
        results = results.split("\n")[1:-1]
        parameters = {
            "structure": "ktp/ktp_pdb_onlytermini_noclean/ktp_noclean.pdb",  # MANDATORY
            "clean_pdb": "no",
            "epsin": 2,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "temp": 310,
            "grid_fill": 0.8,  # FUTURE VERSION
            "ncpus": ncpus,
            "pH": "0,15",
            "pHstep": 0.25,
            "logfile": "LOGFILE",
            "scaleM": 4,
            "scaleP": 1,
            "gsize": 81,
            "convergence": 0.01,
            "nlit": 300,
            "cutoff": -1,
            "relfac": 0.0,
            "output": "pKas.out",
            "clean_pdb": False,
        }
        sites = {" ": ("1N", "2C")}
        pKa = Titration(parameters, sites=sites)

        checkAPIResult(pKa, results)


class TestBuilder(object):
    def test_lyso_pH0(self):
        outfile = "amber1.pdb"

        parameters = {
            "structure": "builder/4lzt.pdb",
            "epsin": 15,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "ncpus": ncpus,
            "output": "pKas.out",
            "titration_output": "titration.out",
            "structure_output": (outfile, 1, "amber"),
        }
        pKa = Titration(parameters)
        checkStructureOutput(outfile)

    def test_lyso_pH7(self):
        outfile = "amber7.pdb"

        parameters = {
            "structure": "builder/4lzt.pdb",
            "epsin": 15,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "ncpus": ncpus,
            "output": "pKas.out",
            "titration_output": "titration.out",
            "structure_output": (outfile, 7, "amber"),
        }
        pKa = Titration(parameters)
        checkStructureOutput(outfile)

    def test_lyso_pH12(self):
        outfile = "amber12.pdb"

        parameters = {
            "structure": "builder/4lzt.pdb",
            "epsin": 15,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "ncpus": ncpus,
            "output": "pKas.out",
            "titration_output": "titration.out",
            "structure_output": (outfile, 12, "amber"),
        }
        pKa = Titration(parameters)
        checkStructureOutput(outfile)

    def test_lyso_pH12_gromos(self):
        outfile = "gromos12.pdb"

        parameters = {
            "structure": "builder/4lzt.pdb",
            "epsin": 15,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "ncpus": ncpus,
            "output": "pKas.out",
            "titration_output": "titration.out",
            "structure_output": (outfile, 12, "gromos_CpH"),
        }
        pKa = Titration(parameters)
        checkStructureOutput(outfile)
