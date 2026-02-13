# (c) 01/2001 by Reinhard Caspary

from Lanthanide import *
from Ion import *
#from ZBLANcauchy import *
import Matrix

PrMeas = [
    [  0, "3H4",   200,  1 ],
    [  2, "3H6",  4487,  1 ],
    [  3, "3F2",  5215,  1 ],
    [  4, "3F3",  6568,  1 ],
    [  5, "3F4",  7031,  1 ],
    [  6, "1G4", 10001,  1 ],
    [  7, "1D2", 17047,  1 ],
    [  8, "3P0", 20927,  1 ],
    [  9, "3P1", 21514,  1 ],
    [ 11, "3P2", 22746,  1 ],
    [ 12, "1S0", 46986,  1 ],
    ]

fmeas = [
    [1, '3H5', 1.9618909300147944e-06, 2.8242516174098058e-08],
    [2, '3H6', 1.3158832643586664e-06, 5.8432411589751691e-08],
    [3, '3F2', 2.29451295506235e-06, 6.3319377294992039e-08],
    [4, '3F3', 6.3702382496830319e-06, 1.6602196614642903e-07],
    [5, '3F4', 2.4415770852447936e-06, 1.4683014926241082e-07],
    [6, '1G4', 2.640652006750471e-07, 4.7583407603229392e-09],
    [7, '1D2', 1.9916417404027569e-06, 1.0654613781849282e-07],
    [8, '3P0', 1.7998263011564285e-06, 2.3152403492628923e-07],
    [9, '3P1', 5.4274727760643657e-06, 3.4335984984657339e-07],
    [-10, '3P2', 9.2672765333618606e-06, 3.122763244585255e-07],
    ]

afmeas = [
    [  2, '3H6', 102.2e-8, 1e-8 ],
    [  3, '3F2', 282.0e-8, 1e-8 ],
    [  4, '3F3', 829.0e-8, 1e-8 ],
    [  5, '3F4', 520.2e-8, 1e-8 ],
    [  6, '1G4',  29.7e-8, 1e-8 ],
    [  7, '1D2', 212.2e-8, 1e-8 ],
    [  8, '3P0', 280.5e-8, 1e-8 ],
    [  9, '3P1', 283.4e-8, 1e-8 ],
    [ 10, '1I6', 144.1e-8, 1e-8 ],
    ]

afmeas = [
    [  2, '3H6',  97.1e-8, 1e-8 ],
    [  3, '3F2', 281.4e-8, 1e-8 ],
    [  4, '3F3', 839.6e-8, 1e-8 ],
    [  5, '3F4', 522.6e-8, 1e-8 ],
    [  6, '1G4',  20.7e-8, 1e-8 ],
    [  7, '1D2', 157.1e-8, 1e-8 ],
    [  8, '3P0', 280.1e-8, 1e-8 ],
    [  9, '3P1', 283.0e-8, 1e-8 ],
    [ 10, '1I6', 141.0e-8, 1e-8 ],
    ]


#PrParms[0][0] = 11393.2
#TmParms[0][0] = 271695.9
#ion = Ion("Pr3+")
#ion.setParms(PrParms)
#ion.setCauchy(cauchy)

#ion.FitLevels(PrMeas)
#ion.ShowLevels(PrMeas)

#ion.ShowJOdata(afmeas)
#ion.FitOsc(afmeas)

#ion.FitDouble(PrMeas, afmeas)

#ion = Ion("Sm3+", PrParms, cauchy)
#termSplit(ion.config)
#sys.exit()

#key = "product"
#a = Matrix.Matrix(ion.config, ion.term, key, "U1:2:-2")
#print ion.termNames(key)
#n = ion.termShortNames(key)
#showMatrix(a[:20,:20], n[:20])
#sys.exit()

#names = ion.termNames()
#short = ion.termShortNames()
#for i in range(len(names)):
#    print "  %-12s:  %s" % (short[i], names[i])

#ion.Fit(PrMeas)
#ion.saveParms("testParms.new")
#on.showLevels(PrMeas)
#ion.showLevels(meas)

#ion.showJOdata()
#ion.JOfit(fmeas)

import Config
import Term
from Utilities import *

l = 3
config = Config.Config(3, 13)
l = config["l"]
term = {}
term["product"] = Term.termProduct(config)
term["sljm"],term["vectors"] = Term.termSLJM(config, term)

#ion = Ion("f13")
#term = ion.term
key = "sljm"
key = "product"
syms   = term[key]["symmetry"]
#jvals  = term[key]["values"][:,syms.index("J2")]
#mvals  = term[key]["values"][:,syms.index("Jz")]
matrix = Matrix.Matrix(config, term, key, "S2", 0)

#print jvals
#print mvals
print term[key]["values"]
showMatrix(matrix)

