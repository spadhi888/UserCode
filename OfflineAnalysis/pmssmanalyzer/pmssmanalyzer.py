#!/usr/bin/env python
# -----------------------------------------------------------------------------
#  File:        pmssmanalyzer.py
#  Description: Analyzer for ntuples created by Mkntuple
#  Created:     Fri Dec  3 11:22:14 2010 by mkntanalyzer.py
#  Author:      Sezen Sekmen
#  $Revision: 1.19 $
# -----------------------------------------------------------------------------
from ROOT import *
from string import *
from pmssmanalyzerlib import *
import os, sys, re
# -----------------------------------------------------------------------------
# -- Procedures and functions
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
def main():

	# Get number of events
	nevents = stream.size()
	print "Number of events:", nevents

	# -------------------------------------------------------------------------
	# Book histograms etc.
	# -------------------------------------------------------------------------
	setStyle()

	hfile = histogramFile(cmdline.histfilename)

	# -------------------------------------------------------------------------
	# Loop over events
	# -------------------------------------------------------------------------
	for entry in xrange(nevents):
		stream.read(entry)

		# -- Event selection
		#if not SUSY: continue

	stream.close()
	hfile.close()
# -----------------------------------------------------------------------------
main()
