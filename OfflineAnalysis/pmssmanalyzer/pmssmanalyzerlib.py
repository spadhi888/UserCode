# -----------------------------------------------------------------------------
#  File:        pmssmanalyzerlib.py
#  Description: Analyzer for ntuples created by Mkntuple
#  Created:     Fri Dec  3 11:22:14 2010 by mkntanalyzer.py
#  Author:      Sezen Sekmen
#  $Revision: 1.19 $
# -----------------------------------------------------------------------------
from ROOT import *
from time import sleep
from string import *
from PhysicsTools.Mkntuple.AutoLoader import *
import os, sys, re
# -----------------------------------------------------------------------------
# -- Classes, procedures and functions
# -----------------------------------------------------------------------------
class skimFile:
	def __init__(self, stream, filename):
		print "events will be skimmed to file", filename
		self.filename = filename
		self.file = TFile(filename, "recreate")
		self.tree = stream.tree().CloneTree(0)
		self.hist = TH1F("counts", "", 1, 0, 1)
		self.hist.SetBit(TH1.kCanRebin)
		self.hist.SetStats(0)
		self.b_weight = 0

	def addEvent(self, weight=-1.0):
		self.file = self.tree.GetCurrentFile()
		self.file.cd()
		if weight > -1:
			if self.b_weight == 0:
				self.weight = Double()
				self.b_weight = self.tree.Branch("eventWeight", self.weight,
												 "eventWeight/D")
		self.tree.Fill()

	def count(self, cond):
		self.hist.Fill(cond, 1)
		
	def close(self):
		print "==> events skimmed to file", self.filename
		self.file = self.tree.GetCurrentFile()
		self.file.cd()
		self.file.Write("", TObject.kOverwrite)
		self.hist.Write("", TObject.kOverwrite)
		self.file.Close()
# -----------------------------------------------------------------------------
class histogramFile:
	def __init__(self, histfilename):
		self.filename = histfilename
		self.file = TFile(self.filename, "recreate")

	def close(self):
		print "==> histograms saved to file", self.filename
		self.file.cd()
		self.file.Write("", TObject.kOverwrite)
		self.file.ls()
		self.file.Close()
# -----------------------------------------------------------------------------
class commandLine:
	def __init__(self):
		pass

def decodeCommandLine():
	argv = sys.argv
	argc = len(argv)

	cl = commandLine()
	cl.progname = split(os.path.basename(argv[0]),'.')[0]

	if argc > 1:
		cl.filelist = argv[1]
	else:
		cl.filelist = "filelist.txt"

	if argc > 2: 
		cl.histfilename = argv[2] # 2nd (optional) command line argument
	else:
		cl.histfilename = cl.progname + "_histograms"

	if argc > 3:
		cl.skimfilename = argv[3]
	else:
		cl.skimfilename = cl.progname + "_skim"

	# Make sure extension is ".root"
	cl.histfilename = os.path.basename(cl.histfilename)
	cl.histfilename = split(cl.histfilename, ".")[0] + ".root"

	cl.skimfilename = os.path.basename(cl.skimfilename)
	cl.skimfilename = split(cl.skimfilename, ".")[0] + ".root"
	return cl
# -----------------------------------------------------------------------------
def error(message):
	print "** error ** " + message
	sys.exit(0)
# -----------------------------------------------------------------------------
#  Read ntuple filenames from file list

def getFilenames(filelist):
	if not os.path.exists(filelist):
		error("unable to open file: " + filelist)

	# Get ntuple file names
	filenames = filter(lambda x: x != "",
					   map(strip, open(filelist).readlines()))
	v = vector("string")()
	for filename in filenames: v.push_back(filename)
	return v
# -----------------------------------------------------------------------------
TEXTFONT = 42
TEXTSIZE = 0.031
#------------------------------------------------------------------------------
def setStyle():
	gROOT.SetStyle("Pub")
	style = gROOT.GetStyle("Pub")

	# For the canvas
	style.SetCanvasBorderMode(0)
	style.SetCanvasColor(kWhite)
	style.SetCanvasDefH(500)
	style.SetCanvasDefW(500)
	style.SetCanvasDefX(0)
	style.SetCanvasDefY(0)

	# For the pad
	style.SetPadBorderMode(0)
	style.SetPadColor(kWhite)
	style.SetPadGridX(kFALSE)
	style.SetPadGridY(kTRUE)
	style.SetGridColor(kGreen)
	style.SetGridStyle(3)
	style.SetGridWidth(1)

	# For the frame
	style.SetFrameBorderMode(0)
	style.SetFrameBorderSize(1)
	style.SetFrameFillColor(0)
	style.SetFrameFillStyle(0)
	style.SetFrameLineColor(1)
	style.SetFrameLineStyle(1)
	style.SetFrameLineWidth(1)

	# For the histogram
	style.SetHistLineColor(1)
	style.SetHistLineStyle(0)
	style.SetHistLineWidth(1)
	style.SetEndErrorSize(2)
	style.SetErrorX(0.)
	style.SetMarkerSize(0.1)
	style.SetMarkerStyle(20)

	# For the fit/function
	style.SetOptFit(1)
	style.SetFitFormat("5.4g")
	style.SetFuncColor(2)
	style.SetFuncStyle(1)
	style.SetFuncWidth(1)

	#For the date
	style.SetOptDate(0)

	# For the statistics box
	style.SetOptFile(0)
	style.SetOptStat("")
	# To display the mean and RMS
	#style.SetOptStat("mr") 
	style.SetStatColor(kWhite)
	style.SetStatFont(TEXTFONT)
	style.SetStatFontSize(TEXTSIZE)
	style.SetStatTextColor(1)
	style.SetStatFormat("6.4g")
	style.SetStatBorderSize(1)
	style.SetStatH(0.2)
	style.SetStatW(0.3)

	# Margins
	style.SetPadTopMargin(0.05)
	style.SetPadBottomMargin(0.16)
	style.SetPadLeftMargin(0.16)
	style.SetPadRightMargin(0.16)

	# For the global title
	style.SetOptTitle(0)
	style.SetTitleFont(TEXTFONT)
	style.SetTitleColor(1)
	style.SetTitleTextColor(1)
	style.SetTitleFillColor(10)
	style.SetTitleFontSize(TEXTSIZE*1.1)

	# For the axis titles
	style.SetTitleColor(1, "XYZ")
	style.SetTitleFont(TEXTFONT, "XYZ")
	style.SetTitleSize(TEXTSIZE*1.2, "XYZ") # 0,05
	style.SetTitleXOffset(1.25) # 0.9
	style.SetTitleYOffset(1.25) # 1.25

	# For the axis labels
	style.SetLabelColor(1, "XYZ")
	style.SetLabelFont(TEXTFONT, "XYZ")
	style.SetLabelOffset(0.006, "XYZ")
	style.SetLabelSize(TEXTSIZE*1.2, "XYZ")

	# For the axis
	style.SetAxisColor(1, "XYZ")
	style.SetStripDecimals(kTRUE)
	style.SetTickLength(0.03, "XYZ")
	style.SetNdivisions(505, "XYZ")

	# To get tick marks on the opposite side of the frame
	style.SetPadTickX(1)  
	style.SetPadTickY(1)

	# Change for log plots
	style.SetOptLogx(0)
	style.SetOptLogy(0)
	style.SetOptLogz(0)

	# Postscript options
	style.SetPaperSize(20.,20.)

	style.cd()
# -----------------------------------------------------------------------------
#  Define variables to be read
# -----------------------------------------------------------------------------
cmdline = decodeCommandLine()

#  Get names of ntuple files to be processed and open chain of ntuples

filenames = getFilenames(cmdline.filelist)
stream = itreestream(filenames, "Analysis GEN")
if not stream.good(): error("unable to open ntuple file(s)")

CaloTower_E	= vector("float")(989)
CaloTower_ET	= vector("float")(989)
CaloTower_E_em	= vector("float")(989)
CaloTower_E_had	= vector("float")(989)
CaloTower_Eta	= vector("float")(989)
CaloTower_Phi	= vector("float")(989)
CaloTower_size	= Long()
ETmis_ET	= vector("float")(3)
ETmis_Phi	= vector("float")(3)
ETmis_Px	= vector("float")(3)
ETmis_Py	= vector("float")(3)
ETmis_size	= Long()
Electron_Charge	= vector("int")(9)
Electron_E	= vector("float")(9)
Electron_EHoverEE	= vector("float")(9)
Electron_EtRatio	= vector("float")(9)
Electron_Eta	= vector("float")(9)
Electron_EtaCalo	= vector("float")(9)
Electron_IsolFlag	= vector("int")(9)
Electron_PT	= vector("float")(9)
Electron_Phi	= vector("float")(9)
Electron_PhiCalo	= vector("float")(9)
Electron_Px	= vector("float")(9)
Electron_Py	= vector("float")(9)
Electron_Pz	= vector("float")(9)
Electron_SumEt	= vector("float")(9)
Electron_SumPt	= vector("float")(9)
Electron_size	= Long()
Event_CouplingQCD	= vector("double")(3)
Event_CouplingQED	= vector("double")(3)
Event_Nparticles	= vector("int")(3)
Event_Number	= vector("long")(3)
Event_ProcessID	= vector("int")(3)
Event_ScalePDF	= vector("double")(3)
Event_Weight	= vector("double")(3)
Event_size	= Long()
FP420hits	= Long()
FP420hits_E	= Double()
FP420hits_S	= Double()
FP420hits_T	= Double()
FP420hits_Tx	= Double()
FP420hits_Ty	= Double()
FP420hits_X	= Double()
FP420hits_Y	= Double()
FP420hits_genE	= Double()
FP420hits_genEta	= Double()
FP420hits_genPT	= Double()
FP420hits_genPhi	= Double()
FP420hits_genPx	= Double()
FP420hits_genPy	= Double()
FP420hits_genPz	= Double()
FP420hits_pid	= Long()
FP420hits_q2	= Double()
FP420hits_side	= Long()
FP420hits_size	= Long()
Jet_Btag	= vector("int")(27)
Jet_E	= vector("float")(27)
Jet_EHoverEE	= vector("float")(27)
Jet_Eta	= vector("float")(27)
Jet_NCalo	= vector("int")(27)
Jet_NTracks	= vector("int")(27)
Jet_PT	= vector("float")(27)
Jet_Phi	= vector("float")(27)
Jet_Px	= vector("float")(27)
Jet_Py	= vector("float")(27)
Jet_Pz	= vector("float")(27)
Jet_size	= Long()
Muon_Charge	= vector("int")(9)
Muon_E	= vector("float")(9)
Muon_EHoverEE	= vector("float")(9)
Muon_EtRatio	= vector("float")(9)
Muon_Eta	= vector("float")(9)
Muon_EtaCalo	= vector("float")(9)
Muon_IsolFlag	= vector("int")(9)
Muon_PT	= vector("float")(9)
Muon_Phi	= vector("float")(9)
Muon_PhiCalo	= vector("float")(9)
Muon_Px	= vector("float")(9)
Muon_Py	= vector("float")(9)
Muon_Pz	= vector("float")(9)
Muon_SumEt	= vector("float")(9)
Muon_SumPt	= vector("float")(9)
Muon_size	= Long()
Particle_Charge	= vector("float")(2537)
Particle_D1	= vector("int")(2537)
Particle_D2	= vector("int")(2537)
Particle_E	= vector("float")(2537)
Particle_Eta	= vector("float")(2537)
Particle_M	= vector("float")(2537)
Particle_M1	= vector("int")(2537)
Particle_M2	= vector("int")(2537)
Particle_PID	= vector("int")(2537)
Particle_PT	= vector("float")(2537)
Particle_Phi	= vector("float")(2537)
Particle_Px	= vector("float")(2537)
Particle_Py	= vector("float")(2537)
Particle_Pz	= vector("float")(2537)
Particle_Status	= vector("int")(2537)
Particle_T	= vector("float")(2537)
Particle_X	= vector("float")(2537)
Particle_Y	= vector("float")(2537)
Particle_Z	= vector("float")(2537)
Particle_size	= Long()
Photon_E	= vector("float")(49)
Photon_EHoverEE	= vector("float")(49)
Photon_Eta	= vector("float")(49)
Photon_PT	= vector("float")(49)
Photon_Phi	= vector("float")(49)
Photon_Px	= vector("float")(49)
Photon_Py	= vector("float")(49)
Photon_Pz	= vector("float")(49)
Photon_size	= Long()
RP220hits	= Long()
RP220hits_E	= Double()
RP220hits_S	= Double()
RP220hits_T	= Double()
RP220hits_Tx	= Double()
RP220hits_Ty	= Double()
RP220hits_X	= Double()
RP220hits_Y	= Double()
RP220hits_genE	= Double()
RP220hits_genEta	= Double()
RP220hits_genPT	= Double()
RP220hits_genPhi	= Double()
RP220hits_genPx	= Double()
RP220hits_genPy	= Double()
RP220hits_genPz	= Double()
RP220hits_pid	= Long()
RP220hits_q2	= Double()
RP220hits_side	= Long()
RP220hits_size	= Long()
TauJet_Charge	= vector("float")(7)
TauJet_E	= vector("float")(7)
TauJet_EHoverEE	= vector("float")(7)
TauJet_Eta	= vector("float")(7)
TauJet_NCalo	= vector("int")(7)
TauJet_NTracks	= vector("int")(7)
TauJet_PT	= vector("float")(7)
TauJet_Phi	= vector("float")(7)
TauJet_Px	= vector("float")(7)
TauJet_Py	= vector("float")(7)
TauJet_Pz	= vector("float")(7)
TauJet_size	= Long()
Tracks_Charge	= vector("float")(357)
Tracks_E	= vector("float")(357)
Tracks_Eta	= vector("float")(357)
Tracks_EtaOuter	= vector("float")(357)
Tracks_PT	= vector("float")(357)
Tracks_Phi	= vector("float")(357)
Tracks_PhiOuter	= vector("float")(357)
Tracks_Px	= vector("float")(357)
Tracks_Py	= vector("float")(357)
Tracks_Pz	= vector("float")(357)
Tracks_Vx	= vector("float")(357)
Tracks_Vy	= vector("float")(357)
Tracks_Vz	= vector("float")(357)
Tracks_size	= Long()
ZDChits_E	= vector("float")(9)
ZDChits_T	= vector("float")(9)
ZDChits_genE	= vector("float")(9)
ZDChits_genEta	= vector("float")(9)
ZDChits_genPT	= vector("float")(9)
ZDChits_genPhi	= vector("float")(9)
ZDChits_genPx	= vector("float")(9)
ZDChits_genPy	= vector("float")(9)
ZDChits_genPz	= vector("float")(9)
ZDChits_hadronic_hit	= vector("int")(9)
ZDChits_pid	= vector("int")(9)
ZDChits_side	= vector("int")(9)
ZDChits_size	= Long()

stream.select("CaloTower.E", CaloTower_E)
stream.select("CaloTower.ET", CaloTower_ET)
stream.select("CaloTower.E_em", CaloTower_E_em)
stream.select("CaloTower.E_had", CaloTower_E_had)
stream.select("CaloTower.Eta", CaloTower_Eta)
stream.select("CaloTower.Phi", CaloTower_Phi)
stream.select("CaloTower_size", CaloTower_size)
stream.select("ETmis.ET", ETmis_ET)
stream.select("ETmis.Phi", ETmis_Phi)
stream.select("ETmis.Px", ETmis_Px)
stream.select("ETmis.Py", ETmis_Py)
stream.select("ETmis_size", ETmis_size)
stream.select("Electron.Charge", Electron_Charge)
stream.select("Electron.E", Electron_E)
stream.select("Electron.EHoverEE", Electron_EHoverEE)
stream.select("Electron.EtRatio", Electron_EtRatio)
stream.select("Electron.Eta", Electron_Eta)
stream.select("Electron.EtaCalo", Electron_EtaCalo)
stream.select("Electron.IsolFlag", Electron_IsolFlag)
stream.select("Electron.PT", Electron_PT)
stream.select("Electron.Phi", Electron_Phi)
stream.select("Electron.PhiCalo", Electron_PhiCalo)
stream.select("Electron.Px", Electron_Px)
stream.select("Electron.Py", Electron_Py)
stream.select("Electron.Pz", Electron_Pz)
stream.select("Electron.SumEt", Electron_SumEt)
stream.select("Electron.SumPt", Electron_SumPt)
stream.select("Electron_size", Electron_size)
stream.select("Event.CouplingQCD", Event_CouplingQCD)
stream.select("Event.CouplingQED", Event_CouplingQED)
stream.select("Event.Nparticles", Event_Nparticles)
stream.select("Event.Number", Event_Number)
stream.select("Event.ProcessID", Event_ProcessID)
stream.select("Event.ScalePDF", Event_ScalePDF)
stream.select("Event.Weight", Event_Weight)
stream.select("Event_size", Event_size)
stream.select("FP420hits", FP420hits)
stream.select("FP420hits.E", FP420hits_E)
stream.select("FP420hits.S", FP420hits_S)
stream.select("FP420hits.T", FP420hits_T)
stream.select("FP420hits.Tx", FP420hits_Tx)
stream.select("FP420hits.Ty", FP420hits_Ty)
stream.select("FP420hits.X", FP420hits_X)
stream.select("FP420hits.Y", FP420hits_Y)
stream.select("FP420hits.genE", FP420hits_genE)
stream.select("FP420hits.genEta", FP420hits_genEta)
stream.select("FP420hits.genPT", FP420hits_genPT)
stream.select("FP420hits.genPhi", FP420hits_genPhi)
stream.select("FP420hits.genPx", FP420hits_genPx)
stream.select("FP420hits.genPy", FP420hits_genPy)
stream.select("FP420hits.genPz", FP420hits_genPz)
stream.select("FP420hits.pid", FP420hits_pid)
stream.select("FP420hits.q2", FP420hits_q2)
stream.select("FP420hits.side", FP420hits_side)
stream.select("FP420hits_size", FP420hits_size)
stream.select("Jet.Btag", Jet_Btag)
stream.select("Jet.E", Jet_E)
stream.select("Jet.EHoverEE", Jet_EHoverEE)
stream.select("Jet.Eta", Jet_Eta)
stream.select("Jet.NCalo", Jet_NCalo)
stream.select("Jet.NTracks", Jet_NTracks)
stream.select("Jet.PT", Jet_PT)
stream.select("Jet.Phi", Jet_Phi)
stream.select("Jet.Px", Jet_Px)
stream.select("Jet.Py", Jet_Py)
stream.select("Jet.Pz", Jet_Pz)
stream.select("Jet_size", Jet_size)
stream.select("Muon.Charge", Muon_Charge)
stream.select("Muon.E", Muon_E)
stream.select("Muon.EHoverEE", Muon_EHoverEE)
stream.select("Muon.EtRatio", Muon_EtRatio)
stream.select("Muon.Eta", Muon_Eta)
stream.select("Muon.EtaCalo", Muon_EtaCalo)
stream.select("Muon.IsolFlag", Muon_IsolFlag)
stream.select("Muon.PT", Muon_PT)
stream.select("Muon.Phi", Muon_Phi)
stream.select("Muon.PhiCalo", Muon_PhiCalo)
stream.select("Muon.Px", Muon_Px)
stream.select("Muon.Py", Muon_Py)
stream.select("Muon.Pz", Muon_Pz)
stream.select("Muon.SumEt", Muon_SumEt)
stream.select("Muon.SumPt", Muon_SumPt)
stream.select("Muon_size", Muon_size)
stream.select("Particle.Charge", Particle_Charge)
stream.select("Particle.D1", Particle_D1)
stream.select("Particle.D2", Particle_D2)
stream.select("Particle.E", Particle_E)
stream.select("Particle.Eta", Particle_Eta)
stream.select("Particle.M", Particle_M)
stream.select("Particle.M1", Particle_M1)
stream.select("Particle.M2", Particle_M2)
stream.select("Particle.PID", Particle_PID)
stream.select("Particle.PT", Particle_PT)
stream.select("Particle.Phi", Particle_Phi)
stream.select("Particle.Px", Particle_Px)
stream.select("Particle.Py", Particle_Py)
stream.select("Particle.Pz", Particle_Pz)
stream.select("Particle.Status", Particle_Status)
stream.select("Particle.T", Particle_T)
stream.select("Particle.X", Particle_X)
stream.select("Particle.Y", Particle_Y)
stream.select("Particle.Z", Particle_Z)
stream.select("Particle_size", Particle_size)
stream.select("Photon.E", Photon_E)
stream.select("Photon.EHoverEE", Photon_EHoverEE)
stream.select("Photon.Eta", Photon_Eta)
stream.select("Photon.PT", Photon_PT)
stream.select("Photon.Phi", Photon_Phi)
stream.select("Photon.Px", Photon_Px)
stream.select("Photon.Py", Photon_Py)
stream.select("Photon.Pz", Photon_Pz)
stream.select("Photon_size", Photon_size)
stream.select("RP220hits", RP220hits)
stream.select("RP220hits.E", RP220hits_E)
stream.select("RP220hits.S", RP220hits_S)
stream.select("RP220hits.T", RP220hits_T)
stream.select("RP220hits.Tx", RP220hits_Tx)
stream.select("RP220hits.Ty", RP220hits_Ty)
stream.select("RP220hits.X", RP220hits_X)
stream.select("RP220hits.Y", RP220hits_Y)
stream.select("RP220hits.genE", RP220hits_genE)
stream.select("RP220hits.genEta", RP220hits_genEta)
stream.select("RP220hits.genPT", RP220hits_genPT)
stream.select("RP220hits.genPhi", RP220hits_genPhi)
stream.select("RP220hits.genPx", RP220hits_genPx)
stream.select("RP220hits.genPy", RP220hits_genPy)
stream.select("RP220hits.genPz", RP220hits_genPz)
stream.select("RP220hits.pid", RP220hits_pid)
stream.select("RP220hits.q2", RP220hits_q2)
stream.select("RP220hits.side", RP220hits_side)
stream.select("RP220hits_size", RP220hits_size)
stream.select("TauJet.Charge", TauJet_Charge)
stream.select("TauJet.E", TauJet_E)
stream.select("TauJet.EHoverEE", TauJet_EHoverEE)
stream.select("TauJet.Eta", TauJet_Eta)
stream.select("TauJet.NCalo", TauJet_NCalo)
stream.select("TauJet.NTracks", TauJet_NTracks)
stream.select("TauJet.PT", TauJet_PT)
stream.select("TauJet.Phi", TauJet_Phi)
stream.select("TauJet.Px", TauJet_Px)
stream.select("TauJet.Py", TauJet_Py)
stream.select("TauJet.Pz", TauJet_Pz)
stream.select("TauJet_size", TauJet_size)
stream.select("Tracks.Charge", Tracks_Charge)
stream.select("Tracks.E", Tracks_E)
stream.select("Tracks.Eta", Tracks_Eta)
stream.select("Tracks.EtaOuter", Tracks_EtaOuter)
stream.select("Tracks.PT", Tracks_PT)
stream.select("Tracks.Phi", Tracks_Phi)
stream.select("Tracks.PhiOuter", Tracks_PhiOuter)
stream.select("Tracks.Px", Tracks_Px)
stream.select("Tracks.Py", Tracks_Py)
stream.select("Tracks.Pz", Tracks_Pz)
stream.select("Tracks.Vx", Tracks_Vx)
stream.select("Tracks.Vy", Tracks_Vy)
stream.select("Tracks.Vz", Tracks_Vz)
stream.select("Tracks_size", Tracks_size)
stream.select("ZDChits.E", ZDChits_E)
stream.select("ZDChits.T", ZDChits_T)
stream.select("ZDChits.genE", ZDChits_genE)
stream.select("ZDChits.genEta", ZDChits_genEta)
stream.select("ZDChits.genPT", ZDChits_genPT)
stream.select("ZDChits.genPhi", ZDChits_genPhi)
stream.select("ZDChits.genPx", ZDChits_genPx)
stream.select("ZDChits.genPy", ZDChits_genPy)
stream.select("ZDChits.genPz", ZDChits_genPz)
stream.select("ZDChits.hadronic_hit", ZDChits_hadronic_hit)
stream.select("ZDChits.pid", ZDChits_pid)
stream.select("ZDChits.side", ZDChits_side)
stream.select("ZDChits_size", ZDChits_size)

