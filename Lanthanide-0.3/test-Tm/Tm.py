import Numeric
from Spectrum import Spectrum
from Lanthanide import *
from Utilities import *
from Ion import Ion

from Matrix import *
from Hamiltonian import Hamiltonian
import TermNames
from parms import *

DB = {
    "element":    "Tm",
    "config":     "f12",
    "levels" :    ["3H6", "3F4", "3H5", "3H4", "3F3", "3F2", "1G4", "1D2",
		   "3P0", "1I6", "3P1", "3P2", "1S0"],
    "first":      0,
    "last":       50000,
    "resolution": 2.0,
    "glass": {
     881:  [ 2.99,  29.e-6,  9.88, 0.10, 13.33, 0.10, 22.23, 0.10, "ZBLANP" ],
    1030:  [ 1.07,  14.e-6,  8.08, 0.10, 11.98, 0.10, 12.83, 0.10, "ZBLANP" ],
    1142:  [ 0.300, 3.9e-6,  8.70, 0.10, 15.35, 0.10, 17.25, 0.10, "ZBLANP" ],
    1104:  [ 0.107, 1.4e-6,  8.45, 0.10, 14.70, 0.10, 17.80, 0.10, "ZBLANP" ],
    1011:  [ 0.00,  0.0e-6,  8.30, 0.10, 15.33, 0.10, 18.98, 0.10, "ZBLANP" ],
    },
    "cauchy" : {
    "ZBLAN":  [1.35123e-5, 2.94780e-3, 1.48965, -1.30933e-3, -3.23335e-6],
    "ZBLANP": [1.35123e-5, 2.94780e-3, 1.49985, -1.30933e-3, -3.23335e-6],
    },
    "absorption": {
    "meas" : [
    [ "20000928-030.csv", 1104, 1, 293 ],  #  0
    [ "20000928-032.csv",  881, 1, 293 ],  #  1
    [ "20000928-033.csv", 1030, 1, 293 ],  #  2
    [ "20000929-002.csv", 1142, 1, 293 ],  #  3
    [ "20001004-004.csv",  881, 2, 293 ],  #  4
    [ "20001004-006.csv",  881, 1, 293 ],  #  5
    [ "20001004-007.csv",  881, 2, 293 ],  #  6
    [ "20001113-002.csv",  881, 1, 293 ],  #  7
    [ "20001113-003.csv",  881, 3, 293 ],  #  8
    [ "20001113-004.csv", 1142, 3, 293 ],  #  9
    [ "20001113-005.csv", 1030, 3, 293 ],  # 10
    ],
    "band" : [ # first, last, corr(meas, span, c/p/b, c/n, data)
    [ 4500, 7000, [        # 00
    [ 0, (), "c", "c", [] ],
    [ 1, (), "c", "c", [] ],
    [ 2, (), "c", "c", [] ],
    [ 5, (), "c", "c", [] ],
    ]],
    [ 7000, 9500, [        # 01
    [ 0, (), "c", "c", [] ],
    [ 1, (), "c", "c", [] ],
    [ 2, (), "c", "c", [] ],
    [ 5, (), "c", "c", [] ],
    ]],
    [ 11800, 13700, [      # 02
    [ 5, (), "c", "c", [] ],
    [ 7, (), "c", "c", [] ],
    [ 8, (), "c", "c", [] ],
    [ 9, (), "c", "c", [] ],
    ]],
    [ 13700, 16000, [      # 03
    [ 0, (), "c", "c", [] ],
    [ 1, (), "c", "c", [] ],
    [ 2, (), "c", "c", [] ],
    [ 5, (), "c", "c", [] ],
    ]],
    [ 20000, 22500, [      # 04
    [ 0, (), "c", "c", [] ],
    [ 1, (), "c", "c", [] ],
    [ 2, (), "c", "c", [] ],
    [ 5, (), "c", "c", [] ],
    [ 6, (), "c", "c", [] ],
    ]],
    [ 26500, 29500, [      # 05
    [ 0, (), "c", "c", [] ],
    [ 1, (), "c", "c", [] ],
    [ 2, (), "c", "c", [] ],
    [ 4, (), "c", "c", [] ],
    [ 5, (), "c", "c", [] ],
    ]],
    [ 33000, 37300, [      # 06
    [ 4, (), "c", "c", [] ],
    [ 5, (), "c", "c", [] ],
    ]],
    [ 37300, 42000, [      # 07
    [ 4, (), "c", "c", [] ],
    [ 5, (), "c", "c", [] ],
    ]],
    ],
    "line" : [ # initial, final, band, split(pos, width), peak(pos, width)
    [ 0, 1, 0,    # 00
      [],
      [(6000, 100)]],
    [ 0, 2, 1,    # 01
      [],
      [(8270, 50)]],
    [ 0, 3, 2,    # 02
      [],
      [(12640, 50)]],
    [ 0, 4, 3,    # 03
      [(14580, 260)],
      [(14620, 200)]],
    [ 0, 5, 3,    # 04
      [(15180, 150)],
      [(15180, 50)]],
    [ 0, 6, 4,    # 05
      [],
      [(21250, 70), (21580, 100)]],
    [ 0, 7, 5,    # 06
      [],
      [(28000, 200)]],
    [ 0, 8, 6,    # 07
      [(34630, 270)],
      [(34700, 200)]],
    [ 0, 9, 6,    # 08
      [(35170, 430)],
      [(35200, 100)]],
    [ 0, 10, 6,   # 09
      [(36570, 340)],
      [(36600, 200)]],
    [ 0, 11, 7,   # 10
      [],
      [(38350, 300)]],
    ],
    },
    "emission": {
    "meas" : [
    [ "20010605.001", 1030, 0, 293 ],  #  0  2070-2670 nm  (Ge-)
    [ "20010606.001", 1030, 0, 293 ],  #  1  1330-2670 nm  (-)
    [ "20010606.002", 1030, 0, 293 ],  #  2  1320-2680 nm  (-)
    [ "20010612.001", 1030, 0, 293 ],  #  3   650-950 nm   (-)
    [ "20010612.002", 1030, 0, 293 ],  #  4   650-950 nm   (-)
    [ "20010612.003", 1030, 0, 293 ],  #  5   650-950 nm   (-)
    [ "20010612.004", 1030, 0, 293 ],  #  6   650-950 nm   (-)
    [ "20010612.005", 1030, 0, 293 ],  #  7   650-950 nm   (-)
    [ "20010612.006", 1030, 0, 293 ],  #  8  1300-1600 nm  (Si)
    [ "20010612.007", 1030, 0, 293 ],  #  9  1300-1600 nm  (Si)
    ],
    "background" : [
    [ "19991221.003",  5600, 6500, 0.2, (1293,1673) ], # 00  850-2000 nm (GaAs)
    [ "19991221.004",  2500, 5600, 0,   (1293,1673) ], # 01 1700-4000 nm (Ge-)
    [ "19991221.008",  2200, 2500, 0.2, (1293,1673) ], # 02 3500-5500 nm (InAs)
    [ "19991221.008",  2000, 2200, 0.2, (1293,1673) ], # 03 3500-5500 nm (InAs)
    ],
    "band": [ # first, last, corr(meas, span, c/p/b, c/n, data)
    [ 3900, 4600, [     # 00:  0, 1, 2
    [ 1, (), "c", "n", [1, 2.9, 5] ],
    [ 2, (), "c", "n", [1, 2.9, 5] ],
    ]],
    [ 4600, 6400, [     # 01:  1, 2
    [ 1, (), "c", "n", [1, 2.9, 5] ],
    [ 2, (), "c", "n", [1, 2.9, 5] ],
    ]],
    [ 6400, 7500, [     # 02:  1, 2, 8, 9
    [ 8, (), "c", "n", [0, 0, 5] ],
    [ 9, (), "c", "n", [0, 0, 5] ],
    ]],
    [ 11000, 14000, [   # 03:  3, 4, 5, 6, 7
    [ 3, (), "l", "n", [(12040, 12130)] ],
    [ 4, (), "l", "n", [(12030, 12140)] ],
    [ 6, (), "l", "n", [(13100, 13400)] ],
    [ 7, (), "l", "n", [(13070, 13170)] ],
    ]],
    ],
    "line": [ # initial, final, band, split(pos, width), peak(pos, width)
    [ 3, 2, 0,      # 00
      [],
      [(4340, 50)]],
    [ 1, 0, 1,      # 01
      [],
      [(5500, 100)]],
    [ 3, 1, 2,      # 02
      [],
      [(6850, 100)]],
    [ 3, 0, 3,       # 03
      [(12450,100)],
      [(12450,100)]],
    ],
    },
    }


def absorption(setup):
    element  = DB["element"]
    setup.absorption.LaTeXOsc(element+"-osc.tab")
    setup.absorption.LaTeXPeak(element+"-peak.tab")
    gnu = {
	"ground": "right",
	"xrange": (0, 40000),
	"yrange": (0, 0.7),
	"mxtics": 1,
	"mytics": 1,
	"shift": [(4, 400), (7, -400)],
	}
    setup.absorption.GnuplotSpect(element+".gnu", element+".data", gnu)


def levelfit(setup, ion):
    element  = DB["element"]
    cauchy   = DB["cauchy"]["ZBLANP"]
    meas     = setup.absorption.getcmlist()
    fmeas    = setup.absorption.getosclist()
    bandspan = setup.absorption.getbandlist()

    #for num,name,f,df in fmeas:
    #    print "%5s: %g" % (name, f*1e8)

    parmFit(ion, element, "",  "1", FIT_SINGLE, meas, fmeas, cauchy)
    parmFit(ion, element, "1", "2", FIT_SINGLE, meas, fmeas, cauchy)
    #parmFit(ion, element, "2", "3", FIT_SINGLE, meas, fmeas, cauchy)
    
    xmax = 40000
    ymax = 1
    #ion.plotOsc(element + "-fit.gnu", xmax, ymax, meas, fmeas)
    #ion.doubleLaTeX(element + "-best.tab", meas, fmeas)
    #ion.vectorLaTeX(element + "-intermediate.tab", meas, fmeas)
    #ion.emissionLaTeX(element + "-emission.tab", meas, fmeas, 0)
    #ion.absorptionLaTeX(element + "-absorption.tab", meas, fmeas, 0)

    #setup.regrate(ion.Aed, ion.Amd)


def emission(setup):
    element  = DB["element"]
    setup.regrate()
    setup.emission.LaTeXOsc(element+"Em-osc.tab")
    setup.emission.LaTeXPeak(element+"Em-peak.tab")
    gnu = {
	"xrange": (2000, 15000),
	"yrange": (0, 0.9),
	"mxtics": 2,
	"mytics": 1,
	}
    setup.emission.GnuplotSpect(element+"Em.gnu", element+"Em.data", gnu)


def sigma(setup):
    element  = DB["element"]
    path = "calc/%s/level/sigma" % element
    sigma = setup.sigmaarray(4)
    sigma = setup.fillMcCumber(sigma)
    setup.WriteBinData(path, sigma)

    
def test(setup):
    setup.doTest(0,1)
    setup.doTest(0,3)


#setup = Spectrum(DB)
ion = Ion(DB["config"], 1)
config = ion.config
term = ion.term
key = "slj"


def f02():
    ion = Ion("f02", 1)
    config = ion.config
    term = ion.term
    
    key = "slj"
    m = Matrix(config, term, key, "H2")
    m = SparseToTri(m)
    w,v,info = diagonalize(DIAG_FAST, DIAG_VEC, m)
    print w
    # n = TermNames.short(term, key)
    # showMatrix(m, n)

def e():
    ion = Ion("f12", 1)
    config = ion.config
    term = ion.term
    
    key = "slj"
    m = Matrix(config, term, key, "H2")
    m = SparseToTri(m)
    w,v,info = diagonalize(DIAG_FAST, DIAG_VEC, m)
    print w
    # n = TermNames.short(term, key)
    # showMatrix(m, n)


#a = ["", "", 225, "", 1089, "", 184041/25.0]
#l = 3
#for k in [0]:
    #cmat = (2*l+1) * wign3j(l,k,l,0,0,0)
    #m = Matrix(config, term, key, "H3:%i" % k)
    #m = Matrix(config, term, key, "H2")
    #m["data"] *= a[k]
    #m["data"] += 990
    #m = SparseToTri(m)
    #w,v,info = diagonalize(DIAG_FAST, DIAG_VEC, m)
    #w += 300
    #print w
    #showMatrix(m, n)
#LLL = "SPDFGHIK"
#k = 2
#for L in range(7):
#    cmat = (2*l+1) * wign3j(l,k,l,0,0,0)
#    print "    %s:  %6.2f" % (LLL[L], SIGNEXP(L) * cmat * cmat * \
#                              wign6j(3,3,k,3,3,L)*a[k])

#absorption(setup)
#levelfit(setup, ion)
#emission(setup)
#sigma(setup)

#parms = BaseParms["Tm"]
h = Hamiltonian(config, term, "product", parms)
e,v = h.diagonalize()
#n,v = h.intermediate()
#for i in range(len(e)):
#    print "%s:  %g" % (n[i], e[i])
for i in range(len(e)):
    print "%8.0f" % e[i]
