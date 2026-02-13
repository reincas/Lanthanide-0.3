# (c) 01/2001 by Reinhard Caspary

from Lanthanide import *
from Ion import *


parms = [
    [ 76400,       "H1:2",     1 ],
    [ 54900,       "H1:4",     1 ],
    [ 37700,       "H1:6",     1 ],
    [  1025,       "H2",       1 ],
    [    20.5,     "H3:0",     1 ],
    [  -560,       "H3:1",     1 ],
    [  1475,       "H3:2",     1 ],
    [   300,       "H4:2",     1 ],
    [    35,       "H4:3",     1 ],
    [    58,       "H4:4",     1 ],
    [  -310,       "H4:6",     1 ],
    [   350,       "H4:7",     1 ],
    [   320,       "H4:8",     1 ],
    [   152,       "cf:2:0",   1 ],
    [  -294,       "cf:4:0",   1 ],
    [  -592,       "cf:6:0",   1 ],
    [   384,       "cf:6:6",   1 ],
    ]

Pm = Ion("Pm3+", parms)
Pm.Term(TERM_SLJ)
#Pm.Fit(energies, 30000)
Pm.hamiltonDiag()
for i in range(len(Pm.energies)):
    print "   %12.6f" % Pm.energies[i]
