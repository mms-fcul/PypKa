import pytest
import os
import subprocess as sb
import sys
from pypka.tests.builder.check_diffs import compareFiles

sys.path.insert(1, "../")

ncpus = -1

# run coverage.sh to generate coverage


def runTest(path, ncpus, results):
    os.system(
        "rm -f {0}/*out {0}/clean*pqr {0}/TMP.gro {0}/delphi_in_*pdb".format(path)
    )
    results_lines = results.split("\n")[1:-1]
    # "python3 -m coverage erase; "
    # "python3 -m coverage run ../../../pypka.py parameters.dat;"
    sb.Popen(
        "cd {0}; "
        "sed -i 's/ncpus .*/ncpus           = {1}/' parameters.dat; "
        "python3 ../../../pypka.py parameters.dat;"
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
    problems = compareFiles(f"builder/{filename}", filename)
    if problems:
        raise Exception(f"Problems found with {filename}")
    else:
        os.remove(filename)


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
1 NTR 7.22153540243797     A
1 LYS 10.373619702838642   A
7 GLU 2.80311227400602     A
13 NTR 2.5473232026418158   B
13 LYS None                 B
15 HIS 5.629283755997259    B
18 ASP 2.498573271797926    B
20 CTR None                 B
20 TYR 9.980641127902214    B
21 NTR None                 C
23 TYR 8.515707586535786    C
24 SER None                 C
33 LYS 10.494265428727472   C
35 GLU 4.0597973173189      C
36 SER None                 C
40 THR None                 C
41 NTR 6.339947119585919    D
43 THR None                 D
47 THR None                 D
48 ASP 1.8867926817431817   D
50 SER None                 D
51 THR None                 D
52 ASP 2.077072679795452    D
53 TYR 10.712174156778474   D
60 SER None                 D
66 ASP 2.286443132931232    D
69 THR None                 D
72 SER None                 D
81 SER None                 D
85 SER None                 D
86 SER None                 D
87 ASP 1.9621765854063706   D
89 THR None                 D
91 SER None                 D
96 LYS None                 D
97 LYS None                 D
100 SER None                 D
101 ASP 3.4122321499169765   D
116 LYS 10.22102559312535    D
118 THR None                 D
119 ASP 2.4155240252772145   D
129 CTR 1.5108601216333621   D
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
        path = "nucleosome/nucleosome_pdb_all"
        results = """
1 NTR 7.911460961415548    A
3 THR None                 A
4 LYS 10.267237072195853   A
6 THR None                 A
9 LYS 8.489974126778785    A
10 SER None                 A
11 THR None                 A
14 LYS 10.039472714543646   A
18 LYS 7.614718902169101    A
22 THR None                 A
23 LYS 9.658094527182026    A
27 LYS 10.817635575614299   A
28 SER None                 A
32 THR None                 A
36 LYS 10.398197647777607   A
37 LYS 11.234171216257204   A
39 HIS 8.619424255379792    A
41 TYR 13.723361232240821   A
45 THR None                 A
50 GLU 3.102139215251607    A
54 TYR 10.83160066303917    A
56 LYS 10.443918445922296   A
57 SER None                 A
58 THR None                 A
59 GLU 2.9811499919185387   A
64 LYS 10.72633048559488    A
73 GLU 3.425082973139858    A
77 ASP 4.563932625411937    A
79 LYS 10.861795321034732   A
80 THR None                 A
81 ASP 2.6472602739726026   A
86 SER None                 A
87 SER None                 A
94 GLU 2.3503348789283875   A
96 SER None                 A
97 GLU 0.2868671554139638   A
99 TYR 12.745682322485207   A
105 GLU 4.3560364988689955   A
106 ASP 2.835617991424802    A
107 THR None                 A
110 CYS 13.400614520901298   A
113 HIS 3.320579413270905    A
115 LYS 10.39996134111197    A
118 THR None                 A
122 LYS 8.964280662045551    A
123 ASP None                 A
133 GLU 0.7206369197553076   A
135 CTR 3.22337196972467     A
"""
        runTest(path, ncpus, results)


class TestAPI(object):
    def test_api_ktp_gro(self):
        from ..pypka import Titration

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
        from ..pypka import Titration

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
        from ..pypka import Titration

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
        from ..pypka import Titration

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
        from ..pypka import Titration

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
        from ..pypka import Titration

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
        from ..pypka import Titration

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
        from ..pypka import Titration

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
        from ..pypka import Titration

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
