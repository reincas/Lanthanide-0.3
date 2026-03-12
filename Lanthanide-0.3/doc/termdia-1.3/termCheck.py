import string
import re
import types

SyntaxError = "SyntaxError"


class Check:

  # Global variables.

  section = ""
  command = ""

  # special match function.

  def fullmatch(self, expr, str, num):
    if len(str) != num:
      return 0
    m = re.match(expr, str)
    if not m:
      return 0
    if m.end()-m.start() != num:
      return 0
    return 1

  # Convert a value or a list of values from mm to pt.

  def mm2pt(self, v):
    if type(v) != types.ListType:
      v = v * 72.0 / 25.4
    else:
      for i in range(len(v)):
	v[i] = v[i] * 72.0 / 25.4
    return v


  # Check wether c is a key in dictionary d and raise an error, if it
  # is not and flag == 1.

  def check_cmd(self, c, d, flag):
    self.command = c
    if d.has_key(c):
      v = d[c]
    else:
      v = []
      if flag:
        s = self.section
	raise SyntaxError, "Command %s is missing in section %s!" % (c, s)
    return v


  # Convert the list of strings v, containing the values of command c
  # in section s, to a list of integers. The argument n specifies the
  # number of integers required, n=0 allows a list of unlimited
  # length, in the case n=1, the integer value is returned instead of
  # a list.

  def check_int(self, v, n):
    if n > 0 and len(v) != n:
      s = self.section
      c = self.command
      raise SyntaxError, \
	    "Command %s in section %s has wrong number of args!" % (c, s)
    try:
      for p in range(len(v)):
        v[p] = string.atoi(v[p])
    except:
      s = self.section
      c = self.command
      raise SyntaxError, \
	    "Command %s in section %s requires integer values!" % (c, s)
    if n == 1:
      v = v[0]
    return v


  # Convert the list of strings v, containing the values of command c
  # in section s, to a list of floating point numbers. The argument n
  # specifies the number of numbers required, n=0 allows a list of
  # unlimited length, in the case n=1, the value is returned instead
  # of a list.

  def check_float(self, v, n):
    if n > 0 and len(v) != n:
      s = self.section
      c = self.command
      raise SyntaxError, \
	    "Command %s in section %s has wrong number of args!" % (c, s)
    try:
      for p in range(len(v)):
        v[p] = string.atof(v[p])
    except:
      s = self.section
      c = self.command
      raise SyntaxError, \
	    "Command %s in section %s requires float values!" % (c, s)
    if n == 1:
      v = v[0]
    return v


  # Check the list of strings v, containing the values of command c in
  # section s. The argument n specifies the number of strings
  # required, n=0 allows a list of unlimited length, in the case n=1,
  # the string itself is returned instead of a list.

  def check_str(self, v, n):
    if n > 0 and len(v) != n:
      s = self.section
      c = self.command
      raise SyntaxError, \
	    "Command %s in section %s has wrong number of args!" % (c, s)
    if n == 1:
      v = v[0]
    return v


  # Checks, wether r is a number or a range "<num1>-<num2>", the
  # numbers are indices for the list of numbers e. The returned value
  # is the average (e[num1]+e[num2])/2. The current section and
  # command are specified in the arguments s and c, respectively.

  def check_range(self, r, e):
    if re.match("^[0-9]+$", r):
      a = b = string.atoi(r)
    else:
      if re.match("^[0-9]+-[0-9]+$", r):
        p = re.search("-", r).start()
        a = string.atoi(r[:p])
        b = string.atoi(r[p+1:])
      else:
        s = self.section
        c = self.command
        raise SyntaxError, \
	      "Command %s in section %s contains wrong range!" % (c, s)
    if a > len(e) or b > len(e):
      s = self.section
      c = self.command
      raise SyntaxError, \
	    "Command %s in section %s contains wrong range!" % (c, s)
    return int((e[a]+e[b])/2)


  # Checks, if str is a string in correct horz, vert, align
  # syntax. The current section and command are specified in the
  # arguments s and c, respectively.

  def check_hva(self, str):
    if not self.fullmatch("[lcr][tcb][lcr]", str, 3):
      s = self.section
      c = self.command
      raise SyntaxError, \
	    "Command %s in section %s specifies wrong position!" % (c, s)
    return str


  # Checks, if str is "l" or "r". The current section and command are
  # specified in the arguments s and c, respectively.

  def check_lr(self, str):
    if not self.fullmatch("[lr]", str, 1):
      s = self.section
      c = self.command
      raise SyntaxError, \
	    "Command %s in section %s specifies wrong position!" % (c, s)
    return str


  # Checks, if str is matching "[lr][ud]". The current section and
  # command are specified in the arguments s and c, respectively.

  def check_sd(self, str):
    if not self.fullmatch("[lr][ud]", str, 2):
      s = self.section
      c = self.command
      raise SyntaxError, \
	    "Command %s in section %s specifies wrong position!" % (c, s)
    return str


  # Checks, if str is a valid arrow type. The current section and
  # command are specified in the arguments s and c, respectively.

  def check_arrow(self, str):
    if not self.fullmatch("[nbcda]", str, 1):
      s = self.section
      c = self.command
      raise SyntaxError, \
	    "Command %s in section %s specifies wrong arrow type!" % (c, s)
    return str


  def __init__(self, section):
    self.section = section
    self.command = ""
