import time
import math
import string
import re
import copy

linewidth   = 1.0 / 15
superscript = 0.4
subscript   = 0.25
boxskip     = 0.5


##### Class: PointSet ########################################################

class PointSet:

  x = []
  y = []

  def append(self, x, y):
    self.x.append(x)
    self.y.append(y)

  def bbox(self, rot):
    xx = copy.copy(self.x)
    yy = copy.copy(self.y)
    rot = rot * math.pi / 180
    a11 = math.cos(rot)
    a21 = math.sin(rot)
    a12 = -a21
    a22 = a11
    for i in range(len(self.x)):
      x = a11*xx[i] + a12*yy[i]
      y = a21*xx[i] + a22*yy[i]
      xx[i] = x
      yy[i] = y
    return [ min(xx), min(yy), max(xx), max(yy), 0 ]

  def clean(self):
    self.x = []
    self.y = []

  def __init__(self):
    self.clean()



##### Class: BBoxSet #########################################################

class BBoxSet(PointSet):

  def append(self, x, y, bbox):
    llx = bbox[0] + x
    lly = bbox[1] + y
    urx = bbox[2] + x
    ury = bbox[3] + y
    PointSet.append(self, llx, lly)
    PointSet.append(self, llx, ury)
    PointSet.append(self, urx, lly)
    PointSet.append(self, urx, ury)

  def __init__(self):
    PointSet.__init__(self)



##### Class: Line ############################################################

class Line:

  x = 0
  y = 0
  len = 0
  rot = 0
  lw = 0
  color = "blue"

  def header(self):
    s=""
    s=s+"%%% x y len rot lw color  <Line>  -\n"
    s=s+"  /Line { gsave cvx exec\n"
    s=s+"    5 3 roll translate setlinewidth 0 setlinecap rotate\n"
    s=s+"    0 0 moveto 0 lineto stroke\n"
    s=s+"  grestore } def \n"
    s=s+"%%% <Line>\n"
    return s

  def bbox(self):
    dy = 0.5 * self.lw
    set = PointSet()
    set.append(0, -dy)
    set.append(self.len, -dy)
    set.append(self.len, dy)
    set.append(0, dy)
    bb = set.bbox(self.rot)
    return [self.x+bb[0], self.y+bb[1], self.x+bb[2], self.y+bb[3], 0]

  def show(self):
    return "%f %f %f %f %f /%s Line\n" % \
	   (self.x, self.y, self.len, self.rot, self.lw, self.color)

  def __init__(self, x, y, len, rot, lw, font, color="black"):
    self.x = x
    self.y = y
    self.len = len
    self.rot = rot
    self.lw = font.fontsize * linewidth
    self.color = color
    if lw == "thick":
      self.lw = self.lw * 2.0



##### Class: Arrow ###########################################################

class Arrow:

  x = 0
  y = 0
  len = 0
  rot = 0
  lw = 0
  color = 0
  scale = 0
  p1 = 0
  p2 = 0
  h = 0
  size = ""

  def header(self):
    s=""
    s=s+"%%% name  <edef> -\n"
    s=s+"  /edef { exch def } def\n"
    s=s+"%%% <edef>\n"
    s=s+"%%% scale  <arrowhead> -\n"
    s=s+"  /arrowhead { 1 dict begin\n"
    s=s+"    0 0 moveto -13 s mul 4 s mul lineto\n"
    s=s+"    -11 s mul 0 lineto -13 s mul -4 s mul lineto\n"
    s=s+"    closepath fill\n"
    s=s+"  end } def\n"
    s=s+"%%% <arrowhead>\n"
    s=s+"%%% arc prog1 prog2 prog3 prog4  <quadrant> -\n"
    s=s+"  /quadrant { 5 dict begin\n"
    s=s+"    /p4 exch def /p3 exch def /p2 exch def /p1 exch def\n"
    s=s+"    /a exch def\n"
    s=s+"    a 90 lt {p1} { a 180 lt {p2} { a 270 lt {p3} {p4}\n"
    s=s+"      ifelse } ifelse } ifelse\n"
    s=s+"  end } def\n"
    s=s+"%%% <quadrant>\n"
    s=s+"%%% x y len rot lw scale color  <Arrow>  -\n"
    s=s+"  /Arrow { gsave cvx exec 1 dict begin\n"
    s=s+"    /s exch def\n"
    s=s+"    5 3 roll translate setlinewidth 0 setlinecap rotate\n"
    s=s+"    dup 0 0 moveto s 11 mul sub 0 lineto stroke\n"
    s=s+"    0 translate\n"
    s=s+"    arrowhead\n"
    s=s+"  end grestore } def\n"
    s=s+"%%% <Arrow>\n"
    s=s+"%%% x y len rot lw scale color  <Darrow>  -\n"
    s=s+"  /Darrow { gsave cvx exec 1 dict begin\n"
    s=s+"    /s exch def\n"
    s=s+"    5 3 roll translate setlinewidth 0 setlinecap\n"
    s=s+"    [10 5] 0 setdash rotate\n"
    s=s+"    dup 0 0 moveto s 11 mul sub 0 lineto stroke\n"
    s=s+"    0 translate\n"
    s=s+"    arrowhead\n"
    s=s+"  end grestore } def\n"
    s=s+"%%% <Darrow>\n"
    s=s+"%%% x y len1 num r len2 rot lw scale color  <Carrow>  -\n"
    s=s+"  /Carrow { gsave cvx exec 8 dict begin\n"
    s=s+"    /s edef /w edef /d edef /l2 edef /r edef /n edef /l1 edef\n"
    s=s+"    translate w setlinewidth 0 setlinecap d rotate\n"
    s=s+"    0 0 moveto l1 w 2 div add 0 lineto stroke\n"
    s=s+"    l1 0 translate\n"
    s=s+"    1 1 n { /i edef 0 0 moveto r 0 r\n"
    s=s+"      i dup 2 idiv 2 mul eq { 180 0 arc } { 180 0 arcn } ifelse\n"
    s=s+"      r 2 mul 0 translate\n"
    s=s+"    } for\n"
    s=s+"    w -2 div 0 l2 0 w s Arrow\n"
    s=s+"  end grestore } def\n"
    s=s+"%%% <Carrow>\n"
    s=s+"%%% x y len rot lw scale color  <Marrow>  -\n"
    s=s+"  /Marrow { gsave cvx exec 6 dict begin\n"
    s=s+"    /s exch def\n"
    s=s+"    setlinewidth 0 setlinecap\n"
    s=s+"    /rot exch dup sin exch cos atan def\n"
    s=s+"    /d exch 2 div def translate\n"
    s=s+"    rot dup dup 180 sub dup rot quadrant rotate\n"
    s=s+"    rot d d neg dup d quadrant 0 translate\n"
    s=s+"    rot 90 rot sub rot 90 sub 270 rot sub rot 270 sub quadrant\n"
    s=s+"    /beta exch 0.5 mul def\n"
    s=s+"    /r d beta sin div def\n"
    s=s+"    0 r beta cos mul neg r 90 beta sub 90 beta add arc stroke\n"
    s=s+"    rot d d neg dup d quadrant 0 translate\n"
    s=s+"    /betas beta 5.5 s mul r atan sub def\n"
    s=s+"    rot betas neg 180 betas add dup betas neg quadrant rotate\n"
    s=s+"    arrowhead\n"
    s=s+"  end grestore } def\n"
    s=s+"%%% <Marrow>\n"
    return s

  def bbox(self):
    dy = 0.5 * self.lw
    set = PointSet()
    if self.f:
      h = 2*self.f + self.lw
      set.append(0, -h)
      set.append(0, h)
      set.append(self.len, -h)
      set.append(self.len, h)
    else:
      set.append(0, -dy)
      set.append(0, dy)
      set.append(self.len, 0)
      set.append(self.p1, -self.h)
      set.append(self.p1, self.h)
    bb = set.bbox(self.rot)
    return [self.x+bb[0], self.y+bb[1], self.x+bb[2], self.y+bb[3], 0]

  def show(self):
    if self.size == "c":
      num = int((self.len - 19.0*self.scale) / (2.0*self.f))
      len1 = (self.len - 11.0*self.scale - num*2.0*self.f) / 2.0
      len2 = 11.0*self.scale + len1
      return "%f %f %f %d %f %f %f %f %f /%s Carrow\n" % \
             (self.x, self.y, len1, num, self.f, len2, \
              self.rot, self.lw, self.scale, self.color)
    if self.size == "d":
      return "%f %f %f %f %f %f /%s Darrow\n" % \
             (self.x, self.y, self.len, self.rot, self.lw, self.scale,
              self.color)
    if self.size == "a":
      return "%f %f %f %f %f %f /%s Marrow\n" % \
             (self.x, self.y, self.len, self.rot, self.lw, self.scale,
              self.color)
    return "%f %f %f %f %f %f /%s Arrow\n" % \
           (self.x, self.y, self.len, self.rot, self.lw, self.scale,
            self.color)


  def __init__(self, x, y, len, rot, size, color, font):
    self.x = x
    self.y = y
    self.len = len
    self.rot = rot
    if size:
      self.size = size
    else:
      self.size = "n"
    
    if color:
      self.color = color
    else:
      self.color = "black"

    self.lw = font.fontsize * linewidth
    self.scale = self.lw * 0.7
    if size == "b":
      self.lw = self.lw * 2.0
      self.scale = self.scale * math.sqrt(2)

    self.p1 = len - 13*self.lw
    self.p2 = len - 11*self.lw
    self.h = 4 * self.lw

    self.f = 0
    if size == "c":
      self.f = 2*self.lw


##### Class: Text ############################################################

class Text:

  x = 0
  y = 0
  color = "black"
  text = ""
  hva = ""
  lskip = 0
  clip = 0
  rot = 0
  font = ""
  lw = 0
  x0 = 0
  y0 = 0
  w = 0
  offset = []
  width = []
  box = []
  rbox = []

  def header(self):
    s=""
    s=s+"%%% string x y font  <Show>  -\n"
    s=s+"  /Show { setfont moveto show } def\n"
    s=s+"%%% <Show>\n"
    return s

  def bbox(self):
    if not self.text:
      return [0, 0, 0, 0, 0]

    if self.rbox:
      return self.rbox

    self.text = re.sub("^[|~]*", "", self.text, 1)
    self.text = re.sub("\|*$", "", self.text)
    self.text = re.sub("\|+", "|", self.text)
    self.text = re.sub("~+", "~", self.text)
    self.text = re.sub("\|*~\|*", "~", self.text)

    x = 0
    y = 0
    u = superscript * self.font.fontsize
    d = subscript * self.font.fontsize
    
    self.offset = []
    self.width = []
    cset = BBoxSet()
    lset = BBoxSet()
    pos = "norm"
    line = 0
    for i in range(len(self.text)):
      c = self.text[i]

      # Superscript:
      if c == "^":
	if pos == "super":
	  pos = "norm"
	  continue
	if pos == "norm":
	  pos = "super"
	  continue
	raise SyntaxError, "Syntax error in string '%s'!" % self.text

      # Subscript:
      if c == "_":
	if pos == "sub":
	  pos = "norm"
	  continue
	if pos == "norm":
	  pos = "sub"
	  continue
	raise SyntaxError, "Syntax error in string '%s'!" % str

      # Newline or horizontal rule:
      if c == "|" or c == "~":
	if pos == "norm":
	  bbox = cset.bbox(0)
	  self.offset.append(bbox[0])
	  bbox[2] = bbox[2] - bbox[0]
	  bbox[0] = 0
	  self.width.append(bbox[2])
	  lset.append(0, -line*self.lskip*self.font.fontsize, bbox)
	  cset.clean()
	  x = 0
	  y = 0
	  line = line + 1
	  continue
	raise SyntaxError, "Syntax error in string '%s'!" % str

      # Normal character:
      size = "norm"
      if pos != "norm":
	size = "small"

      dy = 0
      if pos == "super":
	dy = superscript*self.font.fontsize
      if pos == "sub":
	dy = -subscript*self.font.fontsize

      bbox = copy.copy(self.font.char_bbox(c, size))
      cset.append(x, y+dy, bbox)
      x = x + bbox[4]

    # Garbage collection ;-)
    if pos != "norm":
      raise SyntaxError, "Syntax error in string '%s'!" % str

    bbox = cset.bbox(0)
    self.offset.append(bbox[0])
    bbox[2] = bbox[2] - bbox[0]
    bbox[0] = 0
    self.width.append(bbox[2]-bbox[0])
    lset.append(0, -line*self.lskip*self.font.fontsize, bbox)
    cset.clean()

    bbox = lset.bbox(0)
    self.box = bbox

    horz = string.find("lcr", self.hva[0])
    vert = string.find("bct", self.hva[1])
    align = string.find("lcr", self.hva[2])
    self.x0 = -bbox[0] - (bbox[2]-bbox[0])*horz/2.0
    self.y0 = -bbox[1] - (bbox[3]-bbox[1])*vert/2.0
    self.w  = bbox[2]-bbox[0]

    lset.clean()
    lset.append(self.x0, self.y0, self.box)
    self.rbox = lset.bbox(self.rot)

    self.rbox[0] = self.rbox[0] + self.x
    self.rbox[1] = self.rbox[1] + self.y
    self.rbox[2] = self.rbox[2] + self.x
    self.rbox[3] = self.rbox[3] + self.y
    return self.rbox

  def show(self):
    if not self.text:
      return ""

    if not self.box:
      self.bbox()

    s=""
    s=s+"%%%%%% BeginText (%s)\n" % self.text

    if self.clip >= 0:
      x = self.rbox[0]-self.clip
      y = self.rbox[1]-self.clip
      dx = (self.rbox[2]-self.rbox[0]) + 2*self.clip
      dy = (self.rbox[3]-self.rbox[1]) + 2*self.clip
      s=s+"gsave 1 setgray %f %f %f %f rectfill grestore\n" % (x, y, dx, dy)

    s=s+"gsave %s %f %f translate %f rotate\n" % \
       (self.color, self.x, self.y, self.rot)

    bbox = self.font.char_bbox("0", "normal")
    rskip = bbox[3] + (self.lskip*self.font.fontsize-(bbox[3]-bbox[1]))/2.0

    line = 0
    align = string.find("lcr", self.hva[2])
    x = self.x0 - self.offset[line] + (self.w-self.width[line])*align/2.0
    y = self.y0 - line*self.lskip*self.font.fontsize
    u = superscript * self.font.fontsize
    d = subscript * self.font.fontsize

    pos = "norm"
    for i in range(len(self.text)):
      c = self.text[i]

      # Superscript:
      if c == "^":
	if pos == "super":
	  pos = "norm"
	  continue
	if pos == "norm":
	  pos = "super"
	  continue
	raise SyntaxError, "Syntax error in string '%s'!" % self.text

      # Subscript:
      if c == "_":
	if pos == "sub":
	  pos = "norm"
	  continue
	if pos == "norm":
	  pos = "sub"
	  continue
	raise SyntaxError, "Syntax error in string '%s'!" % self.text

      # Newline or horizontal rule:
      if c == "|" or c == "~":
	if pos == "norm":
	  line = line + 1
	  x = self.x0 - self.offset[line] + (self.w-self.width[line])*align/2.0
	  y = self.y0 - line*self.lskip*self.font.fontsize
	  if c == "~":
	    rule = Line(self.x0, y+rskip, self.w, 0, "thin", self.font)
	    s=s+rule.show()
	  continue
	raise SyntaxError, "Syntax error in string '%s'!" % self.text

      # Normal character:
      font = "F0"
      size = "norm"
      if pos != "norm":
	font = "F1"
	size = "small"

      dy = 0
      if pos == "super":
	dy = superscript*self.font.fontsize
      if pos == "sub":
	dy = -subscript*self.font.fontsize

      bbox = copy.copy(self.font.char_bbox(c, size))
      if c == "(": c = "\\("
      if c == ")": c = "\\)"
      s=s+"  (%s) %f %f %s Show\n" % (c, x, y+dy, font)
      x = x + bbox[4]

    # Garbage collection ;-)
    # llx = self.box[0] + self.x0
    # lly = self.box[1] + self.y0
    # urx = self.box[2] + self.x0
    # ury = self.box[3] + self.y0
    # s=s+"0.3 setlinewidth\n"
    # s=s+"%f 0 moveto %f 0 lineto stroke\n" % (llx, urx)
    # s=s+"%f %f moveto %f %f lineto\n" % (llx, lly, urx, lly)
    # s=s+"%f %f lineto %f %f lineto\n" % (urx, ury, llx, ury)
    # s=s+"closepath stroke\n"
    s=s+"grestore\n"
    s=s+"%%%%%% EndText (%s)\n" % self.text
    return s

  def __init__(self, x, y, text, hva, lskip, clip, rot, font, color=""):
    self.x = x
    self.y = y
    self.text = text
    self.hva = hva
    self.lskip = lskip
    self.rot = rot
    self.clip = font.fontsize * clip
    if color:
      self.color = color
    else:
      self.color = "black"
    self.font = font
    self.lw = font.fontsize * linewidth
    self.x0 = 0
    self.y0 = 0
    self.w = 0
    self.offset = []
    self.width = []
    self.box = []
    self.rbox = []
