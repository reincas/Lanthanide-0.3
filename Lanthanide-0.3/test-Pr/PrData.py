# Copyright 2000,2003 Reinhard Caspary <r.caspary@tu-bs.de>
#
# This file is part of the Python package Lanthanide.
# 
# Lanthanide is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# Lanthanide is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Lanthanide; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

DB = {
  "Element": "Pr",
  "Config" : "f02",
  "Ground" : "3H4",
  "Glass": {
    1109:  [ 3.00,  47.e-6,  7.98, 0.10, 14.08, 0.10, 18.80, 0.10 ],
     958:  [ 1.05,  13.e-6,  8.30, 0.10, 11.70, 0.10, 19.28, 0.10 ],
    1159:  [ 0.220, 3.2e-6,  8.80, 0.10, 14.30, 0.10, 16.75, 0.10 ],
    1107:  [ 0.191, 2.8e-6,  7.93, 0.10, 14.18, 0.10, 18.78, 0.10 ],
    1118:  [ 0.112, 1.7e-6,  8.20, 0.10, 15.15, 0.10, 18.63, 0.10 ],
    1011:  [ 0.00,  0.0e-6,  8.30, 0.10, 15.33, 0.10, 18.98, 0.10 ],
  },
  "Data" : [
    [ "data/meas-00.csv", 1109, 1 ],  #  0
    [ "data/meas-01.csv",  958, 1,],  #  1
    [ "data/meas-02.ifs", 1109, 1,],  #  2
    [ "data/meas-03.ifs", 1109, 1,],  #  3
    [ "data/meas-04.ifs", 1011, 1,],  #  4 (undoped)
  ],
  "Spectrum" : [
    { "Span"  : (1400,3500),  # 00
      "Meas"  : [2, 3, 4],
      "Corr"  : [("p",0,2), ("p",1,2)],
      "Sigma" : [("c",0), ("c",1)],
      "Label" : ["3H5"],
      "Assign": [ 1 ],
      "Peak"  : [(2140,50,0), (2350,100,0)],
    },
    { "Span"  : (3900,5780),  # 01
      "Rubber": (3900,8000),
      "Meas"  : [0, 1],
      "Corr"  : [("c",0), ("c",1)],
      "Sigma" : [("c",0), ("c",1)],
      "Fit"   : [[(4520, 370)],
		 [(4870, 190), (5140, 160)]],
      "Label" : ["3H6", "3F2"],
      "Assign": [ 2, 3 ],
      "Peak"  : [(5150,200,1)],
    },
    { "Span"  : (5780,8000),  # 02
      "Rubber": (3900,8000),
      "Meas"  : [0, 1],
      "Corr"  : [("c",0), ("c",1)],
      "Sigma" : [("c",0), ("c",1)],
      "Fit"   : [[(6470, 260)],
		 [(6940, 170)]],
      "Label" : ["3F3", "3F4"],
      "Assign": [ 4, 5 ],
      "Peak"  : [(6500,150,0), (6950,100,1)],
    },
    { "Span"  : (8500,11500),  # 03
      "Meas"  : [0, 1],
      "Corr"  : [("c",0), ("c",1)],
      "Sigma" : [("c",0), ("c",1)],
      "Label" : ["1G4"],
      "Assign": [ 6 ],
      "Peak"  : [(9850,300,0)],
    },
    { "Span"  : (15500,18500),  # 04
      "Meas"  : [0, 1],
      "Corr"  : [("c",0), ("c",1)],
      "Sigma" : [("c",0), ("c",1)],
      "Label" : ["1D2"],
      "Assign": [ 7 ],
      "Peak"  : [(17000,300,0)],
      },
    { "Span"  : (20000,24000),  # 05
      "Meas"  : [0, 1],
      "Corr"  : [("c",0), ("c",1)],
      "Sigma" : [("c",0), ("c",1)],
      "Fit"   : [[(20860, 130)],
		 [(21520, 400)],
		 [(22650, 300)]],
      "Label" : [ "3P0", ("3P1", "1I6"), "3P2" ],
      "Assign": [ 8, (9, 10), 11 ],
      "Peak"  : [(20900,100,0), (21450,200,1), (22600,200,2)],
      },
    ],
  }
