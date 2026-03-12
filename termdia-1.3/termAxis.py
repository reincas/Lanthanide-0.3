from termCheck import Check
from termUtil import *


class Axis(Check):

  height = []
  tics   = []
  scale  = 0
  label  = ""
  font   = ""
  box    = []
  offset = 0
  Arrow  = ""
  Label  = ""
  Tics   = []
  Nums   = []


  def bbox(self):
    if self.box:
      return self.box

    set = BBoxSet()
    set.append(0, 0, self.Arrow.bbox())
    set.append(0, 0, self.Label.bbox())
    for i in range(len(self.Tics)):
      set.append(0, 0, self.Tics[i].bbox())
      set.append(0, 0, self.Nums[i].bbox())
    self.box = set.bbox(0)

    skip = 1.5 * boxskip * self.font.fontsize
    self.offset = skip - self.box[0]
    self.box[4] = skip - self.box[0] + skip
    self.box[2] = skip + self.box[2] - self.box[0]
    self.box[0] = skip
    
    return self.box


  def show(self):
    s=""
    s=s+"%%% BeginAxis\n"
    s=s+"gsave\n"
    s=s+"%f 0 translate\n" % self.offset

    s=s+self.Arrow.show()
    s=s+self.Label.show()
    for i in range(len(self.Tics)):
      s=s+self.Tics[i].show()
      s=s+self.Nums[i].show()

    s=s+"grestore\n"
    s=s+"%%% EndAxis\n"
    return s


  def __init__(self, config, font):
    Check.__init__(self, "axis")

    list = self.check_cmd("height", config, 1)
    self.height = self.check_float(list, 3)
    self.height = self.mm2pt(self.height)

    list = self.check_cmd("tics", config, 1)
    self.tics = self.check_int(list, 2)
 
    list = self.check_cmd("scale", config, 1)
    self.scale = self.check_int(list, 1)
    self.scale = float(self.height[1]) / (self.scale * self.tics[0])

    list = self.check_cmd("label", config, 1)
    self.label = self.check_str(list, 1)

    self.font = font
    self.box = []
    self.offset = 0

    fs    = font.fontsize
    x0    = 0.5*fs
    x1    = 0.5*fs
    x2    = -0.9*fs
    y1    = self.height[0]
    y2    = self.height[1]
    y3    = self.height[2]
    max   = self.tics[0]
    delta = self.tics[1]
    bb0   = self.font.char_bbox("0", "normal")
    sy    = 0.5 * (bb0[3]-bb0[1])

    self.Arrow = Arrow(0, -y1, y1+y2+y3, 90, "", "", self.font)
    self.Label = Text(x0, y2+y3, self.label, "ltc", 1.3, -1, 0, self.font)
    self.Tics = []
    self.Nums = []
    for yi in range(0, (int(max/delta)+1)*delta, delta):
      y = y2 * yi / max
      s = "%d" % yi
      self.Tics.append(Line(0, y, x1, 180, "thin", self.font))
      self.Nums.append(Text(x2, y-sy, s, "rbl", 1, -1, 0, self.font))
