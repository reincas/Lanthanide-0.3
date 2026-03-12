import math

from termCheck import Check
from termUtil import *

SyntaxError = "SyntaxError"


class Global(Check):

  font   = ""
  offset = 0
  Items = []

  def bbox(self):
    if self.box:
      return self.box

    set = BBoxSet()
    for i in range(len(self.Items)):
      set.append(0, 0, self.Items[i].bbox())
    self.box = set.bbox(0)

    skip = boxskip * self.font.fontsize
    self.offset = skip - self.box[0]
    self.box[2] = skip + self.box[2] - self.box[0]
    self.box[4] = self.box[2] + skip
    self.box[0] = skip

    return self.box


  def show(self):
    s=""
    s=s+"%%%%%% BeginGlobal\n"
    s=s+"gsave\n"
    s=s+"%f 0 translate\n" % self.offset

    for i in range(len(self.Items)):
      s=s+self.Items[i].show()

    s=s+"grestore\n"
    s=s+"%%%%%% EndGlobal\n"
    return s


  def __init__(self, config, font):
    Check.__init__(self, "global")

    fs    = font.fontsize
    space = 0.4*fs

    self.font = font
    self.color = "black"
    self.Items = []

    ## text ##
    list = self.check_cmd("text", config, 0)
    if list:
      if len(list) % 5 != 0:
        raise SyntaxError, "Wrong text command in section global!"
      for i in range(len(list)/5):
        pos   = list[5*i]
        x     = list[5*i+1]
        y     = list[5*i+2]
        skip  = list[5*i+3]
        text  = list[5*i+4]
        pos   = self.check_hva(pos)
        x     = self.check_int([x], 1)
        y     = self.check_int([y], 1)
        skip  = self.check_float([skip], 1)

	horz = string.find("lcr", pos[0])
	vert = string.find("bct", pos[1])

	self.Items.append(Text(x, y, text, pos, skip, -1, 0, font))

    ## arrow ##
    list = self.check_cmd("arrow", config, 0)
    if list:
      if len(list) % 6 != 0:
        raise SyntaxError, "Wrong arrow command in section global!"
      self.arrow = []
      self.Arrow = []
      self.Atext = []
      for i in range(len(list)/6):
        size  = list[6*i]
        x0    = list[6*i+1]
        y0    = list[6*i+2]
        x1    = list[6*i+3]
        y1    = list[6*i+4]
        color = list[6*i+5]
        size  = self.check_arrow(size)
        x0    = self.check_float([x0], 1)
        y0    = self.check_float([y0], 1)
        x1    = self.check_float([x1], 1)
        y1    = self.check_float([y1], 1)

        l     = math.sqrt(pow((x1-x0),2)+pow((y1-y0),2))
        rot   = math.atan2(y1-y0, x1-x0)*180/math.pi

	self.Items.append(Arrow(x0, y0, l, rot, size, color, font))

