import os
import copy
import string
import re
import math
import stdenc
import isoenc
from termCheck import Check

ReadFileFail = "ReadFileFail"
ParseError = "ParseError"

small = 1.0 / math.sqrt(2)


class Font(Check):

  # Global variables.

  libpath  = ""
  fontpath = ""
  fontname = ""
  fontsize = 0
  encoding = ""

  fontlist = {}
  afmname  = ""
  metric   = {}


  # Parse the Fontmap of ghostscript and create the dictionary
  # 'self.fontlist' mapping a fontname (always starting with '/') to
  # another fontname or the name of the file containing this font.

  def get_fontlist(self):
    fn = self.libpath + "/Fontmap.GS"
    if not os.path.isfile(fn):
      fn = self.basepath + "/Fontmap"
      if not os.path.isfile(fn):
        raise ReadFileFail, "Didn't find fontmap!"

    # open file fn:
    try:
      fp = open(fn, "r")
    except:
      raise ReadFileFail, "Can't open file '%s'!" % fn

    fontlist = {}
    lineno = 0
    while 1:
      line = fp.readline()
      lineno = lineno + 1
      if not line:
        break
      line = re.sub("%.*", "", line)
      line = string.strip(line)
      if not line:
        continue

      line = string.split(line)
      if len(line) != 3:
	if (len(line) == 1) and (line[0] == ";"):
	  continue
        raise ParseError, \
  	    "Wrong number of items on line %d in file %s!" % (lineno, fn)

      if len(line[0]) > 0 and line[0][0] == "/":
        key = line[0]
      else:
        if len(line[0]) > 1 and line[0][0] == "(" and line[0][-1] == ")":
	  key = "/" + line[0][1:-1]
        else:
	  raise ParseError, \
  	      "First font name wrong on line %d in file %s!" % (lineno, fn)

      if len(line[1]) > 0 and line[1][0] == "/":
        value = line[1]
      else:
        if len(line[1]) > 1 and line[1][0] == "(" and line[1][-1] == ")":
	  value = line[1][1:-1]
        else:
	  raise ParseError, \
  	      "Second font name wrong on line %d in file %s!" % (lineno, fn)

      if line[2] != ";":
        raise ParseError, \
	      "Semicolon missing on line %d in file %s!" % (lineno, fn)

      fontlist[key] = value

    fp.close()
    self.fontlist = fontlist


  # Get the filename 'self.afmname' of the AFM file of the font
  # 'self.fontname' and test the existence of this file.

  def get_afmname(self):
    file = ""
    name = "/" + self.fontname
    while not file:
      if not self.fontlist.has_key(name):
        raise FontError, "Font %s is missing!" % self.fontname
      name = self.fontlist[name]
      if name[0] != "/":
        file = name

    (file, ext) = os.path.splitext(file)
    fn = self.fontpath + "/" + file + ".afm"
    if not os.path.isfile(fn):
      raise ReadFileFail, "File '%s' is not a valid filename!" % fn
    self.afmname = fn


  # Parse the AFM file and build the dictionary 'self.metric'
  # containing the character metrics contained in the AFM file. Key of
  # the dictionary is the character name, value is a array containing
  # five items: llx, lly, urx, ury and width.

  def get_metric(self):
    # open file AFM file:
    fn = self.afmname
    try:
      fp = open(fn, "r")
    except:
      raise ReadFileFail, "Can't open file '%s'!" % self.afmname

    line = ""
    lineno = 0
    while not re.match("^StartCharMetrics +[0-9]+$", line):
      line = fp.readline()
      lineno = lineno + 1
      if not line:
        break
      line = string.strip(line)

    if not line:
      raise ParseError, "'StartCharMetrics' not found in file %s!" % fn

    metric = {}
    while 1:
      line = fp.readline()
      lineno = lineno + 1
      if not line:
        break
      line = string.strip(line)
      if line == "EndCharMetrics":
        break
      line = string.split(line)

      if len(line) < 15:
        raise ParseError, \
  	    "Wrong number of items on line %d in file %s!" % (lineno, fn)

      if line[0]  != "C" or \
         not self.isint(line[1]) or \
         line[2]  != ";" or \
         line[3]  != "WX" or \
         not self.isint(line[4]) or \
         line[5]  != ";" or \
         line[6]  != "N" or \
         line[8]  != ";" or \
         line[9]  != "B" or \
         not self.isint(line[10]) or \
         not self.isint(line[11]) or \
         not self.isint(line[12]) or \
         not self.isint(line[13]) or \
         line[14] != ";":
        raise ParseError, "Syntax error on line %d in file %s!" % (lineno, fn)

      llx = string.atoi(line[10])/1000.0
      lly = string.atoi(line[11])/1000.0
      urx = string.atoi(line[12])/1000.0
      ury = string.atoi(line[13])/1000.0
      width = string.atoi(line[4])/1000.0
      metric[line[7]] = [ llx, lly, urx, ury, width ]

    if not line:
      raise ParseError, "'EndCharMetrics' not found in file %s!" % fn

    fp.close()
    self.metric = metric


  # Return true, if string str represents an integer number, otherwise
  # false.

  def isint(self, str):
    return re.match("^-?[0-9]+$", str)


  # Define the font in normal size (/F0) and small size (/F1).

  def header(self):
    name = self.fontname
    size = self.fontsize
    font = "%s-%s" % (name, self.encoding)
    enc  = "StandardEncoding"
    if self.encoding == "iso":
      enc = "ISOLatin1Encoding"

    s=""
    s=s+"%%%%%% BeginFont %s\n" % name
    s=s+"  /%s findfont\n" % name
    s=s+"  dup length dict begin\n"
    s=s+"    { 1 index /FID ne { def } { pop pop } ifelse } forall\n"
    s=s+"    /Encoding %s def\n" % enc
    s=s+"    currentdict\n"
    s=s+"  end\n"
    s=s+"  /%s exch definefont pop\n" % font
    s=s+"  /F0 /%s findfont %f scalefont def\n" % (font, size)
    s=s+"  /F1 /%s findfont %f scalefont def\n" % (font, size*small)
    s=s+"%%%%%% EndFont %s\n" % name
    return s


  # Return the character metrics of character char in encoding enc as
  # contained in the metrics dictionary metric, scaled by the fontsize
  # fs.

  def char_bbox(self, char, size):
    code = ord(char)
    if self.encoding == "iso":
      name = isoenc.vector[code]
    else:
      name = stdenc.vector[code]
    m = copy.copy(self.metric[name])
    fs = self.fontsize
    if size == "small":
      fs = fs * small
    for i in range(5):
      m[i] = m[i] * fs
    return m


  # Initialize the global variables with data contained in the config
  # dictionary 'config'.

  def __init__(self, config, basepath, libpath, fontpath):
    Check.__init__(self, "font")

    self.basepath = basepath
    self.libpath = libpath
    self.fontpath = fontpath

    list = self.check_cmd("fontname", config, 1)
    self.fontname = self.check_str(list, 1)

    list = self.check_cmd("fontsize", config, 1)
    self.fontsize = self.check_float(list, 1)

    list = self.check_cmd("encoding", config, 1)
    self.encoding = self.check_str(list, 1)

    if not self.metric:
      self.get_fontlist()
      self.get_afmname()
      self.get_metric()
