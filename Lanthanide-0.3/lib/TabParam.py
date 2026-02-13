# -*- coding: latin-1 -*-
import sys
import os

KEYS = {
    "base":    "base",
    "H1:2":    "$F^2$",
    "H1:4":    "$F^4$",
    "H1:6":    "$F^6$",
    "H2":      "$\\zeta$",
    "H3:0":    "$\\alpha$",
    "H3:1":    "$\\beta$",
    "H3:2":    "$\\gamma$",
    "H4:2":    "$T^2$",
    "H4:3":    "$T^3$",
    "H4:4":    "$T^4$",
    "H4:6":    "$T^6$",
    "H4:7":    "$T^7$",
    "H4:8":    "$T^8$",
    "H5fix":   "$M^0$",
    "H6fix":   "$P^2$",
    "dk":      "$\\dkbar$",
    "Omega:1": "$\\Omega_1$",
    "Omega:2": "$\\Omega_2$",
    "Omega:3": "$\\Omega_3$",
    "Omega:4": "$\\Omega_4$",
    "Omega:5": "$\\Omega_5$",
    "Omega:6": "$\\Omega_6$",
    "df":      "$\\dfbar$",
    }

def TabParam(db):
    element = db["Element"]
    dir = "energy/%s" % element
    
    files = []
    for fn in os.listdir(dir):
        if fn[-5:] == ".parm":
            files.append("%s/%s" % (dir, fn))
    files.sort()
    
    pkeys = []
    okeys = []
    list = []
    for fn in files:
        fp = open(fn, "r")
        exec(fp.read())
        fp.close()
        
        if not list:
            akeys = []
            for val,key,flag in oldparms:
                akeys.append(key)
            list.append([akeys, oldparms, -1, [], [], -1])
        else:
            akeys = list[0][0]
            aparm = list[0][1]
            for val,key,flag in oldparms:
                if not key in akeys:
                    akeys.append(key)
                    aparm.append([val, key, flag])
            list[0][0] = akeys
            list[0][1] = aparm
            
        keys1 = []
        for val,key,flag in parms:
            keys1.append(key)
            if not key in pkeys:
                pkeys.append(key)
        keys2 = []
        for val,key in omega:
            keys2.append(key)
            if not key in okeys:
                okeys.append(key)
    
        list.append([keys1, parms, dk, keys2, omega, df])
    
    cols = (len(list)-1) * "r"
    s = "\\begin{ParmsTab}{%s}\n" % cols
    
    s += "  parameter\n"
    for i in range(len(list)):
        if i < 1:
            val = "initial"
            col = "c|"
        else:
            val = "step~%i" % i
            if i < len(list)-1:
                col = "c"
            else:
                col = "c|"
        s += "  & \\multicolumn{1}{%s}{%s}\n" % (col, val)
    s += "  & \\multicolumn{1}{c|}{unit} \\\\\n"
        
    s += "  \\hline\n"
    for key in pkeys:
        if not key in KEYS.keys():
            raise RuntimeError, "Unknown parameter %s!" % key
        line = [ "%-10s" % KEYS[key] ]
        for keys1,parms,dk,keys2,omega,df in list:
            if key in keys1:
                i = keys1.index(key)
                val = parms[i][0]
                flag = parms[i][2]
                if flag:
                    val = "$%.2f$" % val
                else:
                    val = "$(%.2f)$" % val
            else:
                val = ""
            line.append("%14s" % val)
        line.append("$\\cm$")
        line = " & ".join(line)
        line = "  %s \\\\\n" % line
        s += line
    s += "  \\hline\n"
    
    line = [ "%-10s" % KEYS["dk"] ]
    for keys1,parms,dk,keys2,omega,df in list:
        if dk >= 0:
            val = "$%.1f$" % dk
        else:
            val = ""
        line.append("%14s" % val)
    line.append("$\\cm$")
    line = " & ".join(line)
    line = "  %s \\\\\n" % line
    s += line
    s += "  \\hline\n"
    
    for key in okeys:
        if not key in KEYS.keys():
            raise RuntimeError, "Unknown parameter %s!" % key
        line = [ "%-10s" % KEYS[key] ]
        for keys1,parms,dk,keys2,omega,df in list:
            if key in keys2:
                i = keys2.index(key)
                val = omega[i][0]
                val = "$%.3f$" % (val*1e24)
            else:
                val = ""
            line.append("%14s" % val)
        line.append("$\\spm$")
        line = " & ".join(line)
        line = "  %s \\\\\n" % line
        s += line
    s += "  \\hline\n"
    
    line = [ "%-10s" % KEYS["df"] ]
    for keys1,parms,dk,keys2,omega,df in list:
        if df >= 0:
            val = "$%.1f$" % df
        else:
            val = ""
        line.append("%14s" % val)
    line.append("$10^{-8}$")
    line = " & ".join(line)
    line = "  %s \\\\\n" % line
    s += line
    
    s += "\\end{ParmsTab}\n"
    
    fp = open("%s-parms.tab" % element, "w")
    fp.write(s)
    fp.close()
