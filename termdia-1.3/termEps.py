import time
from termConfig import Config
from termFont import Font
from termAxis import Axis
from termDiagram import Diagram
from termGlobal import Global
from termUtil import BBoxSet, Line, Arrow, Text

class Eps:

  confname = ""
  config = ""
  font = ""
  axis = ""
  diagram = []
  glob = ""
  box = []


  def bbox(self):
    if self.box:
      return self.box

    x = 0
    y = 0
    set = BBoxSet()
    bbox = self.axis.bbox()
    set.append(x, y, bbox)
    x = x + bbox[4]
    for i in range(len(self.diagram)):
      bbox = self.diagram[i].bbox()
      set.append(x, y, bbox)
      x = x + bbox[4]
    self.box = set.bbox(0)
    return self.box


  def header(self):
    bb = self.bbox()
    bbox = "%d %d %d %d" % (int(bb[0]-1.5), int(bb[1]-1.5), 
			    int(bb[2]+1.5), int(bb[3]+1.5))
    date = time.asctime(time.localtime(time.time()))

    s=""
    s=s+"%!PS-Adobe-3.0 EPSF-3.0\n"
    s=s+"%%%%BoundingBox: %s\n" % bbox
    s=s+"%%Creator: Reinhard Caspary\n"
    s=s+"%%%%Title: %s\n" % self.confname
    s=s+"%%%%CreationDate: %s\n" % date
    s=s+"%%Pages: 1\n"
    s=s+"%%%%DocumentFonts: %s\n" % self.font.fontname
    s=s+"%%%%DocumentNeededFonts: %s\n" % self.font.fontname
    s=s+"%%LanguageLevel: 2\n"
    s=s+"%%EndComments\n"
    return s


  def show(self):
    line = Line(0, 0, 0, 0, "", self.font)
    arrow = Arrow(0, 0, 0, 0, "", "", self.font)
    text = Text(0, 0, "", "", 0, 0, 0, self.font)

    s=""
    s=s+self.header()
    s=s+"\n"
    s=s+"100 dict begin\n"
    s=s+"%%BeginProlog\n"
    s=s+line.header()
    s=s+arrow.header()
    s=s+text.header()
    s=s+"/black     { 0.0 0.0 0.0 setrgbcolor } bind def\n"
    s=s+"/white     { 1.0 1.0 1.0 setrgbcolor } bind def\n"
    s=s+"/red       { 1.0 0.0 0.0 setrgbcolor } bind def\n"
    s=s+"/orange    { 0.9 0.7 0.0 setrgbcolor } bind def\n"
    s=s+"/darkred   { 0.7 0.0 0.0 setrgbcolor } bind def\n"
    s=s+"/green     { 0.0 1.0 0.0 setrgbcolor } bind def\n"
    s=s+"/darkgreen { 0.0 0.7 0.0 setrgbcolor } bind def\n"
    s=s+"/blue      { 0.0 0.0 1.0 setrgbcolor } bind def\n"
    s=s+"/darkblue  { 0.0 0.0 0.7 setrgbcolor } bind def\n"
    s=s+"/lightblue { 0.0 1.0 1.0 setrgbcolor } bind def\n"
    s=s+"%%EndProlog\n"
    s=s+"\n"
    s=s+"%%BeginSetup\n"
    s=s+self.font.header()
    s=s+"%%EndSetup\n"
    s=s+"\n"
    s=s+"%%Page: 1\n"
    s=s+"gsave\n"

    x = 0
    y = 0
    bbox = self.axis.bbox()
    s=s+"%f %f translate\n" % (x, y)
    s=s+self.axis.show()
    x = bbox[4]
    for i in range(len(self.diagram)):
      bbox = self.diagram[i].bbox()
      s=s+"%f %f translate\n" % (x, y)
      s=s+self.diagram[i].show()
      x = bbox[4]
    s=s+"grestore\n"

    if self.glob:
      s=s+"gsave\n"
      s=s+self.glob.show()
      s=s+"grestore\n"

    s=s+"%%EndPage: 1\n"
    s=s+"end\n"
    s=s+"\n"
    s=s+"showpage\n"
    s=s+"%%EOF\n"
    return s


  def __init__(self, confname, basepath, libpath, fontpath):
    self.confname = confname
    self.config = Config(confname)
    self.font = Font(self.config.config["font"], basepath, libpath, fontpath)
    self.axis = Axis(self.config.config["axis"], self.font)

    scale  = self.axis.height[1] / (1000.0*self.axis.tics[0])
    depth  = self.axis.height[0]
    height = self.axis.height[1] + self.axis.height[2]

    self.diagram = []
    for i in range(len(self.config.sections)):
      section = self.config.sections[i]
      self.diagram.append(Diagram(self.config.config[section], section, 
				  scale, depth, height, self.font))

    if self.config.config.has_key("global"):
      self.glob = Global(self.config.config["global"], self.font)

    self.box = []
