import ROOT,os

class makeEnv_TCanvas(object):

  def __init__(self,plotName, logy=False):
    self.plotName = plotName
    self.logy = logy

  def __enter__(self):
    self.canvas = ROOT.TCanvas( "canvas" , "canvas" , 10 , 10 , 1000, 750 )
    if self.logy: self.canvas.SetLogy(0)
    return self

  def __exit__(self,*exc):
    # If directory for plotName doesn't exist yet, make it
    plotNameComponents = self.plotName.split('/')
    plotDir = '/'.join(plotNameComponents[:-1]) 
    if not os.path.isdir(plotDir):
      print "Making plot directory {0}".format(plotDir)
      os.system( "mkdir %s" % plotDir )
    self.canvas.Print(self.plotName)
    del self.canvas

def localDrawErrorSummary( plotter , hist , xaxis_label ):
  
  box = ROOT.TBox(1,0,2,0.16)
  box.SetFillColor(ROOT.kGray)
  box.SetFillStyle(3001)
  box.Draw()

  hist.GetXaxis().SetTitle(xaxis_label)
  hist.GetXaxis().SetTitleSize(0.05)

  plotter.DrawErrorSummary(hist,"TL",True,True,0.00001,False,"",True,"",True)  

def writeHist(hist,outFile):
  #exec('localHist = {0}'.format(histName))
  #localHist.SetName(histName)
  outFile.cd()
  print 'Writing {0} to output file'.format(hist.GetName())
  hist.Write()

## Function to rebin a TH1D according to a reference histogram.
## If hist_to_rebin has more bins than referenceHist, the first extra bin will become the overflow bin
## (this is actually the main purpose of this function).
def Rebin(hist_to_rebin, referenceHist, newname = ""):
  newHist = referenceHist.Clone(newname)
  for i in range(newHist.GetNbinsX() + 2):
    newHist.SetBinContent(i, hist_to_rebin.GetBinContent(i))
    newHist.SetBinError(i, hist_to_rebin.GetBinError(i))
  return newHist

## Function to draw a histogram with its overflow bin.
def DrawWithOverflow(hist, canvas, draw_options):
    # Get number of bins
    nbins = hist.GetNbinsX()

    # Divide canvas if it has not already been divided.
    pads = canvas.GetListOfPrimitives()
    npads = len(pads)
    if npads < 2:
      canvas.Divide(2, 1)
    
    # cd to first pad and set pad size and margin.
    canvas.cd(1)
    ROOT.gPad.SetPad(0, 0, 10.0/11, 1.0)
    ROOT.gPad.SetRightMargin(0.02)

    ROOT.gStyle.SetTitleSize(0.05, "t")
    ROOT.gStyle.SetTitleY(1.0)
    ROOT.gStyle.SetTitleX(0.55)

    # Draw the original histogram and obtain view range for use with overflow bin.
    hist.Draw(draw_options)
    ROOT.gPad.Update()
    miny = ROOT.gPad.GetUymin()
    maxy = ROOT.gPad.GetUymax()

    # cd to second pad and set pad size and margin.
    canvas.cd(2)
    ROOT.gPad.SetPad(10.0/11, 0, 1.0, 1.0)
    ROOT.gPad.SetLeftMargin(0.0)
  
    # Create a new histogram for overflow
    overflow = ROOT.TH1D(hist.GetName() + "_overflow", "Overflow", 1, hist.GetBinLowEdge(nbins), hist.GetBinLowEdge(nbins + 1))
    
    # Fill the overflow bin
    overflow.SetBinContent(1, hist.GetBinContent(nbins + 1))
    overflow.SetBinError(1, hist.GetBinError(nbins + 1))

    # Set overflow histogram to same style and y range as the main one.
    # Remove axes ticks and numbering since they will be the same as for the main histogram.
    # Set title to appropriate size and position.
    # Draw and update.
    overflow.SetFillColor(hist.GetFillColor())
    overflow.SetFillStyle(hist.GetFillStyle())
    overflow.SetLineColor(hist.GetLineColor())
    overflow.SetLineStyle(hist.GetLineStyle())
    overflow.SetLineWidth(hist.GetLineWidth())
    overflow.SetMarkerColor(hist.GetMarkerColor())
    overflow.SetMarkerSize(hist.GetMarkerSize())
    overflow.SetMarkerStyle(hist.GetMarkerStyle())
    overflow.SetNdivisions(0, "xy")
    overflow.GetYaxis().SetRangeUser(miny, maxy)
    ROOT.gStyle.SetTitleSize(0.3, "t")
    ROOT.gStyle.SetTitleY(0.295)
    ROOT.gStyle.SetTitleX(0.4)
    overflow.Draw(draw_options)
    ROOT.gPad.Update()

    # Update the canvas. Must return the overflow to prevent it from being deleted by Python prematurely.
    canvas.cd(0)
    canvas.Modified()
    canvas.Update()
    return overflow
