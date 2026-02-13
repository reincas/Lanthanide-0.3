# (c) 01/2001 by Reinhard Caspary

from Lanthanide import *

# V. V. Ravi Kanth Kumar, 1994
parms = [
#    [  10765.18, "offset",   0 ],
    [  10765.18, "offset",   0 ],
    [  67982.00, "H1:2",     0 ],
    [  49260.00, "H1:4",     0 ],
    [  31189.00, "H1:6",     0 ],
    [    751.00, "H2",       0 ],
    [     14.09, "H3:0",     0 ],
    [   -806.00, "H3:1",     0 ],
    [   2080.00, "H3:2",     0 ],
    [      1.72, "H5fix",    0 ],
    [     63.18, "H6fix",    1 ],
#    [     63.18, "H6fix",    1 ],
#    [    674.00, "H6fix",    1 ],
]
energies = [
    [  0, "3H4",     5,  1 ],
    [  3, "3F2",  5105,  1 ],
    [  4, "3F3",  6502,  1 ],
    [  5, "3F4",  6935,  1 ],
    [  7, "1D2", 16930,  1 ],
    [  8, "3P0", 20821,  1 ],
    [  9, "3P1", 21335,  1 ],
    [ 11, "3P2", 22471,  1 ],
    ]
aenergies = [
    [  0, "3H4",     0,  1 ],
    [  3, "3F2",  5118,  1 ],
    [  4, "3F3",  6502,  1 ],
    [  5, "3F4",  6930,  1 ],
    [  7, "1D2", 16926,  1 ],
    [  8, "3P0", 20816,  1 ],
    [  9, "3P1", 21340,  1 ],
    [ 11, "3P2", 22472,  1 ],
    ]


l = 3
electrons = 2
Pr = Ion(l, electrons, parms)
Pr.buildTerm(TERM_SLJ)

diff = 1
while diff > 0.001:
    Pr.parms[0][2] = 0
    Pr.parms[-1][2] = 1
    p,chi2 = Pr.termFit(energies)
    Pr.parms[0][2] = 1
    Pr.parms[-1][2] = 0
    p,nchi2 = Pr.termFit(energies)
    diff = chi2-nchi2
    
