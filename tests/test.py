import pytest
import os
import subprocess as sb
import sys
from .builder.check_diffs import compareFiles

sys.path.insert(1, "../")
from pypka import Titration

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
        "python3 ../../../pypka/__main__.py parameters.dat;"
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
            assert line == results_lines[c].strip()
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


def checkAPIResult(pKa, results, pH=7, extended=False):
    results = results.split("\n")[1:-1]
    i = -1
    for site in pKa:
        result = "{0} {1} {2} {3}".format(
            site.res_name, site.res_number, site.pK, site.getProtState(pH)
        )
        if extended:
            result += " {0} {1} {2}".format(
                site.getMostProbTaut(pH),
                site.getFinalState(pH),
                site.getTitrationCurve(),
            )
        i += 1
        print(result)
        assert result == results[i]

    # raise Exception
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
18 ASP 3.208344304716938    A
35 GLU 4.677308718296143    A
48 ASP 3.017856790123457    A
66 ASP 3.3120436384516005   A
129 CTR 2.056                A
        """
        runTest(path, ncpus, results)

    def test_cli_lyso_pdb_sites(self):
        path = "lyso/lyso_pdb_sites"
        results = """
1 NTR 7.535143990343163
18 ASP 3.2554331123501776
35 GLU 4.619474578486472
48 ASP 2.313423024672422
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
    1 ASP 3.1153294367693944   B  
   35 GLU 4.622119696447051    C  
   35 ASP 2.3149249683601516   D  
   66 ASP 1.8855327203893997   D  
    1 CTR 2.395745255014842    D
        """
        runTest(path, ncpus, results)

    def test_cli_pHLIP_gro(self):
        path = "pHLIP/pHLIP_gro"
        results = """
769 NTR 10.809909308908642   A
770 CYS None                 A
771 GLU 4.010693080721711    A
782 ASP None                 A
793 ASP None                 A
799 ASP 2.0325970005050396   A
801 ASP 1.038801379054068    A
802 GLU 4.124791182626394    A
804 CTR 3.543877365540363    A
        """
        runTest(path, ncpus, results)

    def test_cli_pHLIP_pdb_all(self):
        path = "pHLIP/pHLIP_pdb_all"
        results = """
769 NTR 11.234137974163346
770 CYS 24.296543354078075
771 GLU 4.015179252479023
776 TYR 18.595372225793035
780 TYR 24.399890927119483
782 ASP None
786 THR None
787 THR None
793 ASP 22.778483296053107
799 ASP 0.48049158811597664
801 ASP 1.5057450374588557
802 GLU 3.961777777777778
804 CTR 4.477023346303502
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

    def test_cli_isoelectric_point_P00698(self):
        path = "isoletric_point/P00698"
        results = """
1 NTR 7.436243954862977    A
1 LYS 10.38343658298011    A
7 GLU 3.270669540958661    A
13 LYS 11.185447421177328   A
15 HIS 5.721977625405991    A
18 ASP 3.1467348396556263   A
20 TYR 10.234760852612425   A
23 TYR 9.62216181031512     A
24 SER None                 A
33 LYS 9.066692907995618    A
35 GLU 4.760520208857088    A
36 SER None                 A
40 THR None                 A
43 THR None                 A
47 THR None                 A
48 ASP 2.0749896184192314   A
50 SER None                 A
51 THR None                 A
52 ASP 2.837970294910756    A
53 TYR 11.216992599478031   A
60 SER None                 A
66 ASP 2.621881674928465    A
69 THR None                 A
72 SER None                 A
81 SER None                 A
85 SER None                 A
86 SER None                 A
87 ASP 2.5711121711871336   A
89 THR None                 A
91 SER None                 A
96 LYS 11.164805773721078   A
97 LYS 11.096740732123468   A
100 SER None                 A
101 ASP 3.749534008581276    A
116 LYS 10.104388635210553   A
118 THR None                 A
119 ASP 2.70851766088339     A
129 CTR 1.9159387844979447   A

Predicted Isoelectric Point: 11.337484705471072
        """
        runTest(path, ncpus, results)

    def test_cli_isoelectric_point_P00921(self):
        path = "isoletric_point/P00921"
        results = """    
2 NTR None                 A
3 HIS 6.9526728209531985   A
4 HIS None                 A
7 TYR None                 A
9 LYS None                 A
10 HIS 6.36920164389817     A
14 GLU None                 A
15 HIS 6.070748761775274    A
17 HIS None                 A
18 LYS None                 A
19 ASP None                 A
26 GLU None                 A
32 ASP None                 A
34 ASP None                 A
36 LYS None                 A
41 ASP None                 A
45 LYS None                 A
51 TYR None                 A
53 GLU None                 A
64 HIS None                 A
69 GLU None                 A
70 TYR None                 A
71 ASP None                 A
72 ASP None                 A
75 ASP None                 A
76 LYS None                 A
80 LYS None                 A
81 ASP None                 A
88 TYR None                 A
94 HIS None                 A
96 HIS None                 A
101 ASP None                 A
102 ASP None                 A
106 GLU None                 A
107 HIS None                 A
110 ASP None                 A
112 LYS None                 A
113 LYS None                 A
114 TYR None                 A
117 GLU None                 A
119 HIS None                 A
122 HIS None                 A
126 LYS None                 A
127 TYR None                 A
129 ASP None                 A
138 ASP None                 A
148 LYS None                 A
151 ASP None                 A
158 LYS None                 A
161 ASP None                 A
164 ASP None                 A
167 LYS None                 A
169 LYS None                 A
171 LYS None                 A
174 ASP None                 A
179 ASP None                 A
189 ASP None                 A
190 TYR None                 A
193 TYR None                 A
204 GLU None                 A
212 LYS None                 A
213 GLU None                 A
224 LYS None                 A
233 GLU None                 A
235 GLU None                 A
237 GLU None                 A
251 LYS None                 A
260 CTR None                 A
260 LYS None                 A

Predicted Isoelectric Point: 6.528696803645706
"""
        runTest(path, ncpus, results)

    def test_cli_proteins_2A25(self):
        path = "proteins/2A25"
        results = """
126 TYR 10.098011676580745   A
127 SER None                 A
128 CYS 12.482018561484919   A
130 CYS 11.494419051404346   A
134 SER None                 A
135 CYS 8.748730472535566    A
136 LYS 10.48507993867718    A
140 SER None                 A
142 ASP 3.4184905980010165   A
147 HIS 5.303316661200394    A
150 HIS 6.001450138291333    A
152 HIS 4.514518909600259    A
153 LYS 10.731265301777668   A
154 SER None                 A
156 THR None                 A
157 THR None                 A
161 GLU 4.111811478157238    A
162 ASP 3.107573956380595    A
168 THR None                 A
169 ASP 3.8422890048587366   A
177 ASP 5.102473498233215    A
183 SER None                 A
184 CYS 11.578526853252647   A
188 HIS 6.605217900398844    A
194 GLU 5.5565424748735674   A
195 LYS 11.20974523427757    A
214 THR None                 A
216 LYS 10.35628958773197    A
219 GLU 3.8965649304737213   A
223 TYR 11.754048723897911   A
226 GLU 3.6358084181272132   A
230 HIS 6.274154674360508    A
235 THR None                 A
237 GLU 4.242028277853261    A
239 THR None                 A
242 SER None                 A
244 HIS 6.086237822248795    A
245 GLU 4.073873949934569    A
249 THR None                 A
254 SER None                 A
255 ASP 3.3012508394895903   A
256 CYS 13.259905381431105   A
260 ASP 3.3551277624309392   A
261 THR None                 A
262 SER None                 A
269 GLU 4.332567084078712    A
278 THR None                 A
280 SER None                 A
282 CTR 3.137268188302425    A
282 CYS 10.441979429808733   A
59 NTR 7.640076928505893    B
59 LYS 11.248239751242004   B
67 CTR 3.124341781263022    B
"""
        runTest(path, ncpus, results)

    def test_cli_proteins_4A60(self):
        path = "proteins/4A60"
        results = """
-9 NTR 7.469081336238199    A
-7 THR None                 A
-6 GLU 3.30028880866426     A
-3 TYR 9.663368821020182    A
0 SER None                 A
3 GLU 4.391259703067299    A
8 THR None                 A
10 LYS 10.668427252346392   A
13 SER None                 A
14 SER None                 A
15 GLU 4.113518927621972    A
18 GLU 3.393629148019927    A
19 ASP 3.3952060310500647   A
20 TYR 13.400726190916146   A
22 LYS 10.587554806681887   A
23 GLU 4.166726112258306    A
38 LYS 10.604896474538332   A
40 THR None                 A
42 THR None                 A
44 SER None                 A
46 ASP 3.5814576918993524   A
48 LYS 10.551811942647316   A
51 THR None                 A
54 THR None                 A
55 GLU 3.9517168554054547   A
56 SER None                 A
57 SER None                 A
60 ASP 1.9627320811744386   A
61 THR None                 A
62 LYS 10.474231284632936   A
64 SER None                 A
66 LYS 10.698242405317393   A
69 GLU 3.790949851826194    A
70 GLU 3.965377855887522    A
72 ASP 2.3844213594409776   A
73 GLU 1.7597915806851012   A
74 THR None                 A
75 THR None                 A
77 ASP 1.0748239786785472   A
80 LYS 10.584316060021923   A
82 LYS 11.756319579674468   A
83 SER None                 A
84 THR None                 A
86 THR None                 A
88 GLU 4.062438981868898    A
91 SER None                 A
94 HIS 5.3611310915303      A
97 LYS 11.144764397905758   A
101 LYS 10.389757328030639   A
102 GLU 3.695859810911345    A
103 THR None                 A
104 THR None                 A
106 LYS 10.85503817491885    A
108 LYS 10.54367650587953    A
111 ASP 3.599820110870443    A
112 GLU 4.283998646820027    A
113 LYS 11.49138842312459    A
117 GLU 3.1609913635365587   A
118 CYS 10.448528923424321   A
119 LYS 10.210060404082483   A
125 SER None                 A
126 THR None                 A
129 TYR 9.756175027262813    A
130 GLU 4.113634640175931    A
131 CTR 2.6923277480358756   A
131 LYS 11.507075386341429   A
"""
        runTest(path, ncpus, results)

    def test_cli_proteins_4ZKS(self):
        path = "proteins/4ZKS"
        results = """
16 NTR 9.215638435208982    U  
20 GLU 4.37962882945167     U  
25 GLU 3.7921010035610228   U  
34 TYR None                 U  
37 HIS 4.864467149951481    U  
40 TYR None                 U  
57 HIS 8.524140613604317    U  
60A ASP None                 U  
60B TYR 9.86359298649037     U  
61 LYS 10.986353128340252   U  
62 LYS None                 U  
63A GLU 3.2147847829134926   U  
63 ASP None                 U  
64 TYR None                 U  
67 TYR 9.412326846647545    U  
80 GLU None                 U  
82 LYS 10.636570288480408   U  
84 GLU 3.1707694436408005   U  
86 GLU 3.6547056272935903   U  
91 HIS None                 U  
92 LYS 10.517627816176859   U  
93 ASP None                 U  
94 TYR None                 U  
97 ASP None                 U  
99 HIS 5.970971131673429    U  
100 HIS 6.093368749833151    U  
102 ASP None                 U  
107 LYS None                 U  
110A LYS 10.50891425562974    U  
110B GLU 3.5389448387523803   U  
127 TYR 9.718214825049       U  
129 ASP None                 U  
137 GLU 4.776213266745278    U  
143 LYS 10.517736890729177   U  
144 GLU 3.3573313905576323   U  
148 ASP None                 U  
149 TYR 9.668698742197352    U  
151 TYR None                 U  
153 GLU None                 U  
156 LYS None                 U  
161 LYS 10.555811406875236   U  
165 HIS 5.37000161108426     U  
167 GLU 3.9046268136412268   U  
170B HIS 6.50534220000193     U  
171 TYR None                 U  
172 TYR None                 U  
175 GLU 5.945245596130807    U  
179 LYS 10.768538703344614   U  
185 ASP None                 U  
187 LYS None                 U  
189 ASP None                 U  
194 ASP None                 U  
223 LYS 10.997764043814692   U  
223A ASP None                 U  
224 LYS None                 U  
228 TYR None                 U  
233 HIS 5.734327813928046    U  
241 HIS 5.936690031059966    U  
243 CTR None                 U  
243 LYS 10.666443863806688   U  
1 NTR 6.887950067428067    P  
7 GLU None                 P  
9 HIS 7.387228722464566    P  
12 CTR 3.1105660862151527   P
"""
        runTest(path, ncpus, results)

    def test_cli_proteins_5J6W(self):
        path = "proteins/5J6W"
        results = """
1 NTR 6.970007685858021    A
1 LYS 10.40120113865167    A
15 LYS 10.682987974098058   A
19 CTR 3.2598815614268695   A
19 LYS 10.667911153119093   A
"""
        runTest(path, ncpus, results)

    def test_cli_proteins_1S5M(self):
        path = "proteins/1S5M"
        results = """
1 NTR 7.385571353343865    A
1 SER None                 A
"""
        runTest(path, ncpus, results)

    def test_cli_proteins_4NU1(self):
        path = "proteins/4NU1"
        results = """
383 CTR 2.963029176201373    A
383 NTR 7.371109691614323    B
"""
        runTest(path, ncpus, results)


class TestAPI(object):
    def test_api_ktp_gro(self):
        os.system(
            "rm -f LOG* *out *gro *pdb *pqr *crg *sites cent contributions interactions.dat pkint"
        )
        results = """
NTR 5001 7.932515230635335 ('undefined', 0.8954064318344567)
TYR 1 9.884779718110096 ('protonated', 0.9986978698182121)
CTR 5002 2.8373655675051106 ('deprotonated', 6.875997409570579e-05)
        """
        parameters = {
            "structure": "ktp/ktp_gro/ktp.gro",
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
NTR 5001 7.934819040827242 ('undefined', 0.8959021976401332)
TYR 1 9.78114442380826 ('protonated', 0.9983475157921265)
CTR 5002 3.010253456221198 ('deprotonated', 0.00010237855405755827)
        """
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
            "remove_hs": False,
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
NTR 5001 7.932515230635335 ('undefined', 0.8954064318344567)
TYR 1 9.884779718110096 ('protonated', 0.9986978698182121)
CTR 5002 2.8373655675051106 ('deprotonated', 6.875997409570579e-05)
        """
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
NTR 5001 7.931056701030927 ('undefined', 0.8950914882162796)
CTR 5002 3.0073796144818266 ('deprotonated', 0.00010170339324314826)
        """
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
            "remove_hs": False,
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
NTR 5001 8.199054937290231 ('protonated', 0.9405274528487423)
CTR 5002 3.099830508474576 ('deprotonated', 0.00012582758427877238)
        """
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

    def test_api_antibody(self):
        os.system(
            "rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent "
            "contributions interactions.dat pkint LOG* __pycache__ *pyc"
        )
        results = """
NTR 5001 None ('undefined', 'pk Not In Range') (4, (0.883505, 0.883505)) 3 {7.0: 0.883505}
GLU 6 None ('undefined', 'pk Not In Range') (5, (0.998505, 0.998505)) 4 {7.0: 0.001495}
SER 7 None ('undefined', 'pk Not In Range') (3, (1.0, 0.45602)) 2 {7.0: 1.0}
SER 17 None ('undefined', 'pk Not In Range') (1, (1.0, 0.33859)) 0 {7.0: 1.0}
SER 21 None ('undefined', 'pk Not In Range') (3, (1.0, 0.46141)) 1 {7.0: 1.0}
CYS 22 None ('undefined', 'pk Not In Range') (2, (0.99653, 0.360205)) 1 {7.0: 0.99653}
SER 25 None ('undefined', 'pk Not In Range') (2, (1.0, 0.460985)) 1 {7.0: 1.0}
LYS 218 None ('undefined', 'pk Not In Range') (4, (0.999715, 0.999715)) 3 {7.0: 0.999715}
SER 219 None ('undefined', 'pk Not In Range') (1, (0.9999999999999999, 0.44376)) 0 {7.0: 1.0}
CTR 5220 None ('undefined', 'pk Not In Range') (5, (0.999825, 0.999825)) 4 {7.0: 0.000175}
CYS 220 None ('undefined', 'pk Not In Range') (3, (0.9984500000000001, 0.376525)) 2 {7.0: 0.99845}
NTR 5001 None ('undefined', 'pk Not In Range') (4, (0.89596, 0.89596)) 3 {7.0: 0.89596}
THR 5 None ('undefined', 'pk Not In Range') (3, (1.0, 0.50905)) 0 {7.0: 1.0}
SER 7 None ('undefined', 'pk Not In Range') (1, (1.0, 0.44146)) 1 {7.0: 1.0}
SER 9 None ('undefined', 'pk Not In Range') (3, (1.0, 0.41821)) 2 {7.0: 1.0}
SER 12 None ('undefined', 'pk Not In Range') (3, (1.0, 0.36406)) 0 {7.0: 1.0}
SER 14 None ('undefined', 'pk Not In Range') (2, (1.0, 0.382485)) 2 {7.0: 1.0}
ASP 17 None ('undefined', 'pk Not In Range') (5, (0.99974, 0.99974)) 4 {7.0: 0.00026}
THR 20 None ('undefined', 'pk Not In Range') (2, (1.0, 0.39265)) 0 {7.0: 1.0}
THR 22 None ('undefined', 'pk Not In Range') (3, (1.0, 0.365115)) 2 {7.0: 1.0}
CYS 23 None ('undefined', 'pk Not In Range') (2, (0.99302, 0.358615)) 2 {7.0: 0.99302}
GLU 214 None ('undefined', 'pk Not In Range') (5, (0.99813, 0.99813)) 4 {7.0: 0.00187}
CTR 5215 None ('undefined', 'pk Not In Range') (5, (0.999775, 0.999775)) 4 {7.0: 0.000225}
CYS 215 None ('undefined', 'pk Not In Range') (2, (0.99854, 0.37496)) 1 {7.0: 0.99854}
NTR 5332 None ('undefined', 'pk Not In Range') (4, (0.854645, 0.854645)) 3 {7.0: 0.854645}
HIS 332 None ('undefined', 'pk Not In Range') (2, (0.875545, 0.70737)) 1 {7.0: 0.124455}
THR 333 None ('undefined', 'pk Not In Range') (2, (1.0, 0.46669)) 2 {7.0: 1.0}
CYS 336 None ('undefined', 'pk Not In Range') (3, (0.985125, 0.336865)) 0 {7.0: 0.985125}
GLU 340 None ('undefined', 'pk Not In Range') (5, (0.99765, 0.99765)) 4 {7.0: 0.00235}
CYS 525 None ('undefined', 'pk Not In Range') (1, (0.9649099999999999, 0.335965)) 0 {7.0: 0.96491}
CTR 5527 None ('undefined', 'pk Not In Range') (5, (0.99998, 0.99998)) 4 {7.0: 2e-05}
LYS 527 None ('undefined', 'pk Not In Range') (4, (0.999805, 0.999805)) 3 {7.0: 0.999805}
        """
        params = {
            "structure": "proteins/Ab269/Ab269-RBD.pdb",
            "ncpus": ncpus,
            "epsin": 15,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "pH": "7",
            "structure_output": ("proteins/Ab269/Ab269-RBD_7.pdb", 7, "AMBER"),
        }
        pKa = Titration(params)

        checkAPIResult(pKa, results, extended=True)
        check_file_diff("builder/Ab269-RBD_7.pdb", "proteins/Ab269/Ab269-RBD_7.pdb")

    def test_api_penta_asp(self):
        os.system(
            "rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent "
            "contributions interactions.dat pkint LOG* __pycache__ *pyc"
        )
        results = """
ASP 4 None ('undefined', 'pk Not In Range') (5, (0.879815, 0.879815)) 4 {4.2: 0.120185}
"""
        params = {
            "structure": "penta/ASP/asp1.gro",
            "ncpus": ncpus,
            "epsin": 2,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "pH": "4.2",
            "clean_pdb": False,
            "sts": "sts_cphmd",
        }
        sites = {"A": [4]}
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=4.2, extended=True)

        results = """
ASP 4 None ('undefined', 'pk Not In Range') (1, (0.5596, 0.33636)) 4 {4.2: 0.5596}
"""
        params["structure"] = "penta/ASP/asp2.gro"
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=4.2, extended=True)

    def test_api_penta_ctr(self):
        os.system(
            "rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent "
            "contributions interactions.dat pkint LOG* __pycache__ *pyc"
        )
        results = """
CTR 5006 None ('undefined', 'pk Not In Range') (5, (0.886665, 0.886665)) 4 {4.0: 0.113335}
"""
        params = {
            "structure": "penta/CTR/ctr1.gro",
            "ncpus": ncpus,
            "epsin": 2,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "pH": "4.0",
            "clean_pdb": False,
            "sts": "sts_cphmd",
        }
        sites = {"A": ["6C"]}
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=4.0, extended=True)

        results = """
CTR 5006 None ('undefined', 'pk Not In Range') (5, (0.57106, 0.57106)) 4 {4.0: 0.42894}
"""
        params["structure"] = "penta/CTR/ctr2.gro"
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=4.0, extended=True)

    def test_api_penta_cys(self):
        os.system(
            "rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent "
            "contributions interactions.dat pkint LOG* __pycache__ *pyc"
        )
        results = """
CYS 4 None ('undefined', 'pk Not In Range') (3, (0.925685, 0.479055)) 2 {8.2: 0.925685}
"""
        params = {
            "structure": "penta/CYS/protein_000.gro",
            "ncpus": ncpus,
            "epsin": 2,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "pH": "8.2",
            "clean_pdb": False,
            "sts": "sts_cphmd",
        }
        sites = {"A": ["4"]}
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=8.2, extended=True)

        results = """
CYS 4 None ('undefined', 'pk Not In Range') (2, (0.61114, 0.28937)) 3 {8.2: 0.61114}
"""
        params["structure"] = "penta/CYS/protein_001.gro"
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=8.2, extended=True)

    def test_api_penta_glu(self):
        os.system(
            "rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent "
            "contributions interactions.dat pkint LOG* __pycache__ *pyc"
        )
        results = """
GLU 4 None ('undefined', 'pk Not In Range') (5, (0.73255, 0.73255)) 4 {4.5: 0.26745}
"""
        params = {
            "structure": "penta/GLU/protein_000.gro",
            "ncpus": ncpus,
            "epsin": 2,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "pH": "4.5",
            "clean_pdb": False,
            "sts": "sts_cphmd",
        }
        sites = {"A": ["4"]}
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=4.5, extended=True)

        results = """
GLU 4 None ('undefined', 'pk Not In Range') (5, (0.55045, 0.55045)) 4 {4.5: 0.44955}
"""
        params["structure"] = "penta/GLU/protein_003.gro"
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=4.5, extended=True)

    def test_api_penta_his(self):
        os.system(
            "rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent "
            "contributions interactions.dat pkint LOG* __pycache__ *pyc"
        )
        results = """
HIS 4 None ('undefined', 'pk Not In Range') (3, (0.859865, 0.859865)) 0 {6.2: 0.859865}
"""
        params = {
            "structure": "penta/HIS/protein_000.gro",
            "ncpus": ncpus,
            "epsin": 2,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "pH": "6.2",
            "clean_pdb": False,
            "sts": "sts_cphmd",
        }
        sites = {"A": ["4"]}
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=6.2, extended=True)

        results = """
HIS 4 None ('undefined', 'pk Not In Range') (2, (0.933485, 0.900285)) 1 {6.2: 0.066515}
"""
        params["structure"] = "penta/HIS/protein_003.gro"
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=6.2, extended=True)

    def test_api_penta_lys(self):
        os.system(
            "rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent "
            "contributions interactions.dat pkint LOG* __pycache__ *pyc"
        )
        results = """
LYS 4 None ('undefined', 'pk Not In Range') (2, (0.651625, 0.219385)) 3 {10.7: 0.348375}
"""
        params = {
            "structure": "penta/LYS/protein_000.gro",
            "ncpus": ncpus,
            "epsin": 2,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "pH": "10.7",
            "clean_pdb": False,
            "sts": "sts_cphmd",
        }
        sites = {"A": ["4"]}
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=10.7, extended=True)

        results = """
LYS 4 None ('undefined', 'pk Not In Range') (2, (0.898065, 0.71377)) 3 {10.7: 0.101935}
"""
        params["structure"] = "penta/LYS/protein_006.gro"
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=10.7, extended=True)

    def test_api_penta_ntr(self):
        os.system(
            "rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent "
            "contributions interactions.dat pkint LOG* __pycache__ *pyc"
        )
        results = """
NTR 5001 None ('undefined', 'pk Not In Range') (4, (0.608115, 0.608115)) 3 {7.7: 0.608115}
"""
        params = {
            "structure": "penta/NTR/protein_000.gro",
            "ncpus": ncpus,
            "epsin": 2,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "pH": "7.7",
            "clean_pdb": False,
            "sts": "sts_cphmd",
        }
        sites = {"A": ["1N"]}
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=7.7, extended=True)

        results = """
NTR 5001 None ('undefined', 'pk Not In Range') (2, (0.736045, 0.447095)) 2 {7.7: 0.263955}
"""
        params["structure"] = "penta/NTR/protein_007.gro"
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=7.7, extended=True)

    def test_api_penta_tyr(self):
        os.system(
            "rm -rf LOG* *out *gro *pdb *pqr *crg *sites cent "
            "contributions interactions.dat pkint LOG* __pycache__ *pyc"
        )
        results = """
TYR 4 None ('undefined', 'pk Not In Range') (3, (0.52272, 0.52272)) 0 {9.5: 0.47728}
"""
        params = {
            "structure": "penta/TYR/protein_001.gro",
            "ncpus": ncpus,
            "epsin": 2,
            "ionicstr": 0.1,
            "pbc_dimensions": 0,
            "pH": "9.5",
            "clean_pdb": False,
            "sts": "sts_cphmd",
        }
        sites = {"A": ["4"]}
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=9.5, extended=True)

        results = """
TYR 4 None ('undefined', 'pk Not In Range') (2, (0.7652749999999999, 0.390635)) 0 {9.5: 0.765275}
"""
        params["structure"] = "penta/TYR/protein_005.gro"
        pKa = Titration(params, sites=sites)
        checkAPIResult(pKa, results, pH=9.5, extended=True)


# class TestBuilder(object):
#    def test_lyso_pH0(self):
#        outfile = "amber1.pdb"
#
#        parameters = {
#            "structure": "builder/4lzt.pdb",
#            "epsin": 15,
#            "ionicstr": 0.1,
#            "pbc_dimensions": 0,
#            "ncpus": ncpus,
#            "output": "pKas.out",
#            "titration_output": "titration.out",
#            "structure_output": (outfile, 1, "amber"),
#        }
#        pKa = Titration(parameters)
#        checkStructureOutput(outfile)
#
#    def test_lyso_pH7(self):
#        outfile = "amber7.pdb"
#
#        parameters = {
#            "structure": "builder/4lzt.pdb",
#            "epsin": 15,
#            "ionicstr": 0.1,
#            "pbc_dimensions": 0,
#            "ncpus": ncpus,
#            "output": "pKas.out",
#            "titration_output": "titration.out",
#            "structure_output": (outfile, 7, "amber"),
#        }
#        pKa = Titration(parameters)
#        checkStructureOutput(outfile)
#
#    def test_lyso_pH12(self):
#        outfile = "amber12.pdb"
#
#        parameters = {
#            "structure": "builder/4lzt.pdb",
#            "epsin": 15,
#            "ionicstr": 0.1,
#            "pbc_dimensions": 0,
#            "ncpus": ncpus,
#            "output": "pKas.out",
#            "titration_output": "titration.out",
#            "structure_output": (outfile, 12, "amber"),
#        }
#        pKa = Titration(parameters)
#        checkStructureOutput(outfile)
#
#    def test_lyso_pH12_gromos(self):
#        outfile = "gromos12.pdb"
#
#        parameters = {
#            "structure": "builder/4lzt.pdb",
#            "epsin": 15,
#            "ionicstr": 0.1,
#            "pbc_dimensions": 0,
#            "ncpus": ncpus,
#            "output": "pKas.out",
#            "titration_output": "titration.out",
#            "structure_output": (outfile, 12, "gromos_CpH"),
#        }
#        pKa = Titration(parameters)
#        checkStructureOutput(outfile)
