import string
import re

ReadFileFail = "ReadFileFail"
SyntaxError = "SyntaxError"
ParseError = "ParseError"


class Config:

  # Global variables.

  confname = ""
  config   = {}
  sections = []


  # Parse the string given as the argument and split it into a list of
  # strings. One or more spaces are interpreted as delimiters. Strings
  # containing spaces must be enclosed by doublequotes.
  
  def split(self, str, lineno):
    s = ""
    v = []
    quote = 0
    for c in str:
      if c == "\"":
        if quote:
	  quote = 0
        else:
	  quote = 1
        if s or not quote:
	  v.append(s)
	  s = ""
        continue

      if quote:
        s = s + c
        continue
	
      if c == " ":
        if s:
	  v.append(s)
	  s = ""
        continue

      s = s + c

    if quote:
      raise SyntaxError, "Closing quotes missing in line %d!" % lineno
    if s:
      v.append(s)

    return v


  # Read the configuration file 'self.confname' and build the
  # dictionary 'self.config' containing an other dictionary per
  # section containing the commands and values. No interpretation and
  # syntax checking of the values is done here.
  
  def parse(self):
    # open config file:
    try:
      fp = open(self.confname, "r")
    except:
      raise ReadFileFail, "Can't open file '%s'!" % self.confname
  
    # scan config file:
    c = {}
    line = ""
    section = ""
    lineno = 0
    sections = []
    while 1:
  
      # read a complete line including continuation lines
      # and strip comments:
      if not line or line[-1] != "\\":
        line = ""
      else:
        line = line[0:-1] + " "
      line = line + fp.readline()
      lineno = lineno + 1
      if not line:
        break
      line = re.sub("\\\\?#.*", "", line, 1)
      line = string.expandtabs(line, 1)
      line = string.strip(line)
      if not line or line[-1] == "\\":
        continue
  
      # section labels are included in square brackets:
      if re.match("^\[.*\]$", line):
        section = line[1:-1]
        if not c.has_key(section):
          c[section] = {}
	  if not section == "font" and \
	     not section == "axis" and \
	     not section == "global":
	    sections.append(section)
        continue

      # the first nonempty line must contain a section identifier:
      if not section:
        raise ParseError, "Missing section on line %d!" % lineno
  
      # split command name and values:
      p = string.find(line, "=")
      if p == -1:
        raise ParseError, "Missing '=' on line %d!" % lineno
      if p == 0:
        raise ParseError, "Missing command on line %d!" % lineno
      if p == len(line):
        raise ParseError, "Missing values on line %d!" % lineno
      cmd = string.strip(line[:p])
      val = string.strip(line[p+1:])
      val = self.split(val, lineno)
      if not c[section].has_key(cmd):
        c[section][cmd] = val
      else:
        c[section][cmd] = c[section][cmd] + val
  
    # close the file and exit:
    fp.close()
    self.config = c
    self.sections = sections

  

  def __init__(self, confname):
    self.confname = confname
    self.parse()
