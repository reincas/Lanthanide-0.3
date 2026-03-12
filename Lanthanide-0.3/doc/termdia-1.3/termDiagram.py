from termCheck import Check
from termUtil import *

SyntaxError = "SyntaxError"


class Diagram(Check):

  scale  = 0
  depth  = 0
  height = 0
  width  = 0
  energy = []
  num    = 0
  lskip  = 0
  rskip  = 0
  font   = ""
  box    = []
  offset = 0
  Energy = []
  Arrow  = []
  Atext  = []

  def bbox(self):
    if self.box:
      return self.box

    set = BBoxSet()
    for i in range(len(self.Energy)):
      set.append(0, 0, self.Energy[i].bbox())
    for i in range(len(self.Label)):
      set.append(0, 0, self.Label[i].bbox())
    for i in range(len(self.Arrow)):
      set.append(0, 0, self.Arrow[i].bbox())
    for i in range(len(self.Atext)):
      set.append(0, 0, self.Atext[i].bbox())
    for i in range(len(self.Text)):
      set.append(0, 0, self.Text[i].bbox())
    self.box = set.bbox(0)

    skip = boxskip * self.font.fontsize
    self.offset = skip - self.box[0]
    self.box[2] = skip + self.box[2] - self.box[0]
    self.box[4] = self.box[2] + skip
    self.box[0] = skip

    return self.box


  def show(self):
    s=""
    s=s+"%%%%%% BeginDiagram %s\n" % self.section
    s=s+"gsave\n"
    s=s+"%f 0 translate\n" % self.offset

    for i in range(len(self.Energy)):
      s=s+self.Energy[i].show()
    for i in range(len(self.Label)):
      s=s+self.Label[i].show()
    for i in range(len(self.Atext)):
      s=s+self.Atext[i].show()
    for i in range(len(self.Arrow)):
      s=s+self.Arrow[i].show()
    for i in range(len(self.Text)):
      s=s+self.Text[i].show()

    s=s+"grestore\n"
    s=s+"%%%%%% EndDiagram %s\n" % self.section
    return s


  def __init__(self, config, section, scale, depth, height, font):
    Check.__init__(self, section)

    fs    = font.fontsize
    space = 0.4*fs

    self.section = section
    self.scale = scale
    self.depth = depth
    self.height = height
    self.font = font
    self.color = "black"

    ## color ##
    list = self.check_cmd("color", config, 0)
    if len(list) > 1:
        raise SyntaxError, "Wrong color command in section %s!" % section
    if list:
      self.color = list[0]
  
    ## width ##
    list = self.check_cmd("width", config, 1)
    self.width = self.check_float(list, 1)
    self.width = self.mm2pt(self.width)
  
    ## energy ##
    list = self.check_cmd("energy", config, 1)
    self.energy = self.check_int(list, 0)
    self.Energy = []
    for e in self.energy:
      y = scale * e
      self.Energy.append(Line(0, y, self.width, 0, "thin", self.font, self.color))
  
    ## label ##
    list = self.check_cmd("label", config, 0)
    self.Label = []
    if list:
      if len(list) % 4 != 0:
        raise SyntaxError, "Wrong label command in section %s!" % section
      for i in range(len(list)/4):
        pos   = list[4*i]
        level = list[4*i+1]
        skip  = list[4*i+2]
        text  = list[4*i+3]
        pos   = self.check_lr(pos)
        level = self.check_range(level, self.energy)
        skip  = self.check_float([skip], 1)
	y = level * scale
	if pos == "l":
	  x = -space
	  self.Label.append(Text(x, y, text, "rcr", skip, -1, 0, font))
	else:
	  x = self.width + space
	  self.Label.append(Text(x, y, text, "lcl", skip, -1, 0, font))
  
    ## text ##
    list = self.check_cmd("text", config, 0)
    self.Text = []
    if list:
      if len(list) % 5 != 0:
        raise SyntaxError, "Wrong text command in section %s!" % section
      for i in range(len(list)/5):
        pos   = list[5*i]
        dx    = list[5*i+1]
        dy    = list[5*i+2]
        skip  = list[5*i+3]
        text  = list[5*i+4]
        pos   = self.check_hva(pos)
        dx    = self.check_int([dx], 1) # unbenutzt !
        dy    = self.check_int([dy], 1) # unbenutzt !
        skip  = self.check_float([skip], 1)

	horz = string.find("lcr", pos[0])
	vert = string.find("bct", pos[1])

	x = self.width * horz / 2.0
	y = -self.depth + (self.depth+self.height) * vert / 2.0
	self.Text.append(Text(x, y, text, pos, skip, -1, 0, font))

    ## num ##
    list = self.check_cmd("num", config, 1)
    self.num = self.check_int(list, 1)
  
    ## lskip ##
    list = self.check_cmd("lskip", config, 1)
    self.lskip = self.check_float(list, 1)
    self.lskip = self.lskip*font.fontsize
  
    ## rskip ##
    list = self.check_cmd("rskip", config, 1)
    self.rskip = self.check_float(list, 1)
    self.rskip = self.rskip*font.fontsize
  
    ## arrow ##
    list = self.check_cmd("arrow", config, 0)
    if self.num > 1:
      xskip = float(self.width - self.lskip - self.rskip) / (self.num - 1)
    else:
      xskip = 0.0
    if list:
      if len(list) % 7 != 0:
        raise SyntaxError, "Wrong arrow command in section %s!" % section
      self.Arrow = []
      self.Atext = []
      for i in range(len(list)/7):
        size  = list[7*i]
        x     = list[7*i+1]
        y1    = list[7*i+2]
        y2    = list[7*i+3]
        pos   = list[7*i+4]
        text  = list[7*i+5]
        color = list[7*i+6]
        size  = self.check_arrow(size)
        x     = self.lskip+self.check_float([x], 1)*xskip
        y1    = self.check_range(y1, self.energy)*scale
        y2    = self.check_range(y2, self.energy)*scale
        pos   = self.check_sd(pos)

	dir   = 90.0*(y2-y1)/abs(y2-y1)
	dx = 0.2*font.fontsize
	if pos[0] == "l":
	  dx = -dx
	rot = 90
	if pos[1] == "d":
	  rot = -rot
	if pos == "lu": hva = "cbc"
	if pos == "ld": hva = "ctc"
	if pos == "ru": hva = "ctc"
	if pos == "rd": hva = "cbc"
	y = (y1 + y2) / 2.0

	self.Arrow.append(Arrow(x, y1, abs(y2-y1), dir, size, color, font))
	self.Atext.append(Text(x+dx, y, text, hva, 1.3, 0.1, rot, font, color))

    ## farrow ##
    list = self.check_cmd("farrow", config, 0)
    if self.num > 1:
      xskip = float(self.width - self.lskip - self.rskip) / (self.num - 1)
    else:
      xskip = 0.0
    if list:
      if len(list) % 7 != 0:
        raise SyntaxError, "Wrong farrow command in section %s!" % section
      for i in range(len(list)/7):
        size  = list[7*i]
        x     = list[7*i+1]
        y1    = list[7*i+2]
        dy    = list[7*i+3]
        pos   = list[7*i+4]
        text  = list[7*i+5]
        color = list[7*i+6]
        size  = self.check_arrow(size)
        x     = self.lskip+self.check_float([x], 1)*xskip
        y1    = self.check_range(y1, self.energy)*scale
        dy    = self.check_float([dy], 1)*scale
        pos   = self.check_sd(pos)

	dir   = 90.0*(dy)/abs(dy)
	dx = 0.2*font.fontsize
	if pos[0] == "l":
	  dx = -dx
	rot = 90
	if pos[1] == "d":
	  rot = -rot
	if pos == "lu": hva = "cbc"
	if pos == "ld": hva = "ctc"
	if pos == "ru": hva = "ctc"
	if pos == "rd": hva = "cbc"
	y = y1 + dy/2.0

	self.Arrow.append(Arrow(x, y1, dy, dir, size, color, font))
	self.Atext.append(Text(x+dx, y, text, hva, 1.3, 0.1, rot, font, color))

    self.box = []
    self.offset = 0
