########################################################################
#       Makefile for ISAJET on a generic Unix system.
#
#       Makes ISAJET library, executables, and auxiliary files from 
# isajet.car in the working directory.
#
#       Usage: Edit the Installation Parameters section. Uncomment the
# appropriate definition of NL for your version of make and comment out
# .IGNORE if possible. Set MAKETMP, CERN, [CERNBIN], [CERNLIB], PDFLIB,
# RANLUX, MACHINE, F77, FFLAGS, LINK, LFLAGS, RANLIB, MAKEFILE, YPATCHY,
# and YCRADLE correctly for your system. Then type 
#       make [>& make.log &]            # make everything
#       -or-
#       make isazeb [>& make.log &]     # Zebra version
#       make clean                      # clean up directory
#       make superclean                 # delete all created files
# It should not be necessary to change anything else. Use s5make
# rather than make on DEC Ultrix.
#
#       The following files should be produced and should remain after
# the make clean step:
#       isajet.tex   documentation
#       libisajet.a  library 
#       aldata.o     block data 
#       isadecay.dat decay table
#       isajet.x     executable to generate events
#       isasusy.x    executable to calculate SUSY decays (MSSM input)
#       isasugra.x   executable to calculate SUSY decays (SUGRA input)
#
# Thanks to Arthur Kreymer (FNAL) for suggesting many improvements to
# the original version.
#
#       ISATOOLS is a package of programs interfaced to ISASUGRA to 
# calculate SUSY contributions to various quantities:
#       DAMU: Anomalous muon magnetic moment
#       ISABSG: Branching ratio for b -> s gamma
#       ISADET: Cross section for direct dark matter detection
#       ISARED: Cold dark matter relic density.
# Set the parameters ISATOOLS and LIBTOOLS. First type
#       make isatools
# This will take some time -- there are 1.4M lines of code. Then proceed 
# as above to make the rest of Isajet.
#
#       22 Nov 1996: nypatchy does not seem to work with the cradle
# piped to standard input. Use a temporary file YCRADLE instead.
#
#       16 Nov 1999: Multiline echo fails on RedHat Linux 6.1. Use 
# multiple single-line echos instead -- this is not pretty but works
# everywhere. Hence NL is no longer needed.
#
#       6 April 2005: Added Isatools.
#
########################################################################

#       Do not change this -- you will break everything.
SHELL = /bin/csh

########## Installation parameters #####################################

#       Some versions of make (e.g. Linux) treat inconsequential warnings
# as errors and stop. Ignore all errors; comment out if possible:
.IGNORE:

#       Name for temporary subdirectory. This will be deleted and 
# recreated:
MAKETMP = ./MakeTmpDir

#       Link with CERN libraries. This is recommended for HP, IBMRT, or 
# LINUX to get date and time information but is not required:
#CERN = NOCERN
CERN = CERN

#       Directories for Cernlib, terminated with a /. These must be set
# if you select CERN or PDFLIB:
#CERNDIR = /cosmic/hd1/hires/cernlib/2006_binary/
CERNDIR = /tas07/disk00/jribnik/cernlib/2006/
CERNBIN = $(CERNDIR)bin/
CERNLIB = $(CERNDIR)lib/

#       PDFLIB support. PDFLIB (part of the CERN library) contains many
# partion distributions but produces a larger executable image:
#PDFLIB =
PDFLIB = PDFLIB

#       RANLUX support. RANLUXprovides better random numbers than the
# standard RANF, and each seed gives an independent sequence. It is
# more difficult to restart and slightly slower.
RANLUX = RANLUX
#RANLUX = NORANLUX

#       Computer type. Most only differ in minor ways, e.g., the system
# clock interface. Uncomment the right one:
# ANSI          ANSI Fortran
# DECS          DECStation with Ultrix
# HPUX          HP/700 with HPUX
# IBMRT         IBM RS/6000 with AIX
# LINUX         PC with LINUX and f2c/gcc or g77 or Absoft f77
# OSF           Digital OSF on Alpha processor
# SGI           Silicon Graphics with IRIX
# SUN           Sparc or PC with Solaris
#MACHINE = ANSI
#MACHINE = DECS
#MACHINE = HPUX
#MACHINE = IBMRT
MACHINE = LINUX
#MACHINE = OSF
#MACHINE = SGI
#MACHINE = SUN

#       Fortran 77 compiler.
#F77 = f77
#       For HPUX; note that f77 on HPUX is brain dead. fort77 is also an
# alternate f2c/gcc interface on Linux:
#F77 = fort77
#       GNU Fortran compiler
#F77 = g77
#       GNU 4.x Fortran compiler
F77 = gfortran

#       Fortran compile flags: Want -c, highest safe optimization, all
# variables static, and postpended underscore for Cernlib compatibility.
# Uncomment the right one or make your own:
#       DECS
#FFLAGS = -c -O -static -w
#       HPUX, -O is broken
#FFLAGS = -c -O2 -K +ppu -w
#       IBMRT
#FFLAGS = -c -O -qextname -w
#       LINUX with f2c/gcc or g77. The -fno-silent displays progress...
#FFLAGS = -c -O -w -fno-automatic -fno-silent
FFLAGS = -c -O -w -fno-automatic 
#       LINUX with gfortran
#FFLAGS = -c -O -w -fno-automatic
#       LINUX with Absoft Fortran
#FFLAGS = -c -O -w -f -s
#       OSF
#FFLAGS= -c -O -w
#       SGI
#FFLAGS = -c -O -static -w
#       SUN - the second option traps floating errors
#FFLAGS = -c -O -w
#FFLAGS = -c -g -w -fnonstd -ftrap=common

#       Use fsplit or CERN's fcasplit to split the whole Fortran source into
# individual subroutines:
FSPLIT = $(CERNBIN)fcasplit
#FSPLIT = fsplit

#       Link command for Fortran, c.f. $(F77):
LINK = $(F77)

#       Extra libraries (e.g. Cernlib) and flags for link. If you want 
# PDFLIB support, you must also define $(PDFLIB) above. The Zebra
# version requires libpacklib.a.
#LFLAGS =
#LFLAGS = -L$(CERNLIB) -lmathlib -lkernlib
LFLAGS = -L$(CERNLIB) -lpdflib804 -lmathlib -lkernlib
#LFLAGS = -L$(CERNLIB) -lpdflib -lpacklib -lmathlib -lkernlib

#       Ranlib command if required (it is for Sun 4.x and Linux):
RANLIB = ranlib
#RANLIB = $(ECHO)

#       Name for this Makefile:
MAKEFILE = Makefile

#       Patchy commands. YPATCHY can be either ypatchy or nypatchy,
# perhaps with a full path. A temporary cradle file YCRADLE is required;
# it is deleted and recreated:
YPATCHY = $(CERNBIN)ypatchy
#YPATCHY = $(CERNBIN)nypatchy
YCRADLE = ./tmpcradle.cra

#       Echo command. /bin/echo works on most (all?) systems; plain echo
# fails on AIX 4.1:
ECHO = /bin/echo

#       Remove command:
RM = /bin/rm

#       Parameter to use Isatools

#ISATOOLS = ISATOOLS
ISATOOLS =

#LIBTOOLS = -lisared
LIBTOOLS =

########## End of installation parameters ##############################

all:    isajet.tex isadecay.dat libisajet.a aldata.o isajet.x \
        isasusy.x isasugra.x

isazeb: isajet.tex isadecay.dat libisazeb.a aldata.o isazeb.x

########## Extract files ###############################################

isajet.f:       isajet.car $(MAKEFILE)                                ;\
  $(RM) isajet.f                                                      ;\
  $(RM) isajet.lis                                                    ;\
  $(ECHO) "+USE,*ISAJET,$(PDFLIB),$(MACHINE)"            > $(YCRADLE) ;\
  $(ECHO) "+USE,$(CERN),$(RANLUX),$(ISATOOLS)"          >> $(YCRADLE) ;\
  $(ECHO) "+EXE."                                       >> $(YCRADLE) ;\
  $(ECHO) "+PAM,T=C."                                   >> $(YCRADLE) ;\
  $(ECHO) "+QUIT"                                       >> $(YCRADLE) ;\
  $(YPATCHY) 'isajet.car isajet.f $(YCRADLE) isajet.lis .GO'          ;\
  $(RM) $(YCRADLE)

isazeb.f:       isajet.car $(MAKEFILE)                                ;\
  $(RM) isazeb.f                                                      ;\
  $(RM) isazeb.lis                                                    ;\
  $(ECHO) "+USE,*ISAZEB,*ISZRUN,$(PDFLIB),$(MACHINE)"    > $(YCRADLE) ;\
  $(ECHO) "+USE,P=ISARUN,D=ISARUN,T=INHIBIT"            >> $(YCRADLE) ;\
  $(ECHO) "+USE,P=ISZRUN,D=ISZRUN,T=INHIBIT"            >> $(YCRADLE) ;\
  $(ECHO) "+EXE."                                       >> $(YCRADLE) ;\
  $(ECHO) "+PAM,T=C."                                   >> $(YCRADLE) ;\
  $(ECHO) "+QUIT."                                      >> $(YCRADLE) ;\
  $(YPATCHY) 'isajet.car isazeb.f $(YCRADLE) isazeb.lis .GO'          ;\
  $(RM) $(YCRADLE)

aldata.f:       isajet.car $(MAKEFILE)                                ;\
  $(RM) aldata.f                                                      ;\
  $(RM) aldata.lis                                                    ;\
  $(ECHO) "+USE,*ISAJET,$(MACHINE)"                      > $(YCRADLE) ;\
  $(ECHO) "+USE,$(CERN)"                                >> $(YCRADLE) ;\
  $(ECHO) "+USE,$(RANLUX)"                              >> $(YCRADLE) ;\
  $(ECHO) "+EXE,ISADATA"                                >> $(YCRADLE) ;\
  $(ECHO) "+PAM,T=C."                                   >> $(YCRADLE) ;\
  $(ECHO) "+QUIT."                                      >> $(YCRADLE) ;\
  $(YPATCHY) 'isajet.car aldata.f $(YCRADLE) aldata.lis .GO'          ;\
  $(RM) $(YCRADLE)

isajet.tex:     isajet.car $(MAKEFILE)                                ;\
  $(RM) isajet.tex                                                    ;\
  $(RM) isatex.lis                                                    ;\
  $(ECHO) "+USE,*ISATEXT,$(MACHINE)"                     > $(YCRADLE) ;\
  $(ECHO) "+USE,$(CERN)"                                >> $(YCRADLE) ;\
  $(ECHO) "+USE,$(RANLUX)"                              >> $(YCRADLE) ;\
  $(ECHO) "+EXE"                                        >> $(YCRADLE) ;\
  $(ECHO) "+PAM,T=C."                                   >> $(YCRADLE) ;\
  $(ECHO) "+QUIT"                                       >> $(YCRADLE) ;\
  $(YPATCHY) 'isajet.car isajet.tex $(YCRADLE) isatex.lis .GO'        ;\
  $(RM) $(YCRADLE)

isadecay.dat:   isajet.car $(MAKEFILE)                                ;\
  $(RM) isadecay.dat                                                  ;\
  $(RM) isadecay.lis                                                  ;\
  $(ECHO) "+USE,*ISADECAY"                               > $(YCRADLE) ;\
  $(ECHO) "+EXE"                                        >> $(YCRADLE) ;\
  $(ECHO) "+PAM,T=C."                                   >> $(YCRADLE) ;\
  $(ECHO) "+QUIT"                                       >> $(YCRADLE) ;\
  $(YPATCHY) 'isajet.car isadecay.dat $(YCRADLE) isadecay.lis .GO'    ;\
  $(RM) $(YCRADLE) 

ssrun.f:  isajet.car $(MAKEFILE)                                      ;\
  $(RM) ssrun.f                                                       ;\
  $(RM) ssrun.lis                                                     ;\
  $(ECHO) "+USE,*ISASUSY,$(MACHINE)"                     > $(YCRADLE) ;\
  $(ECHO) "+EXE,ISASSRUN"                               >> $(YCRADLE) ;\
  $(ECHO) "+PAM,T=C."                                   >> $(YCRADLE) ;\
  $(ECHO) "+QUIT"                                       >> $(YCRADLE) ;\
  $(YPATCHY) 'isajet.car ssrun.f $(YCRADLE) ssrun.lis .GO'            ;\
  $(RM) $(YCRADLE)

sugrun.f:       isajet.car $(MAKEFILE)                                ;\
  $(RM) sugrun.f                                                      ;\
  $(RM) sugrun.lis                                                    ;\
  $(ECHO) "+USE,*ISASUGRA,$(MACHINE),$(ISATOOLS)"        > $(YCRADLE) ;\
  $(ECHO) "+EXE,ISASSRUN"                               >> $(YCRADLE) ;\
  $(ECHO) "+PAM,T=C."                                   >> $(YCRADLE) ;\
  $(ECHO) "+QUIT"                                       >> $(YCRADLE) ;\
  $(YPATCHY) 'isajet.car sugrun.f $(YCRADLE) sugrun.lis .GO'          ;\
  $(RM) $(YCRADLE)

isarun.f:       isajet.car $(MAKEFILE)                                ;\
$(RM) isarun.f                                                        ;\
$(RM) isarun.lis                                                      ;\
  $(ECHO) "+USE,*ISAZEB,*ISZRUN,$(PDFLIB),$(MACHINE)"    > $(YCRADLE) ;\
  $(ECHO) "+EXE,P=ISARUN,D=ISARUN"                      >> $(YCRADLE) ;\
  $(ECHO) "+PAM,T=C."                                   >> $(YCRADLE) ;\
  $(ECHO) "+QUIT"                                       >> $(YCRADLE) ;\
  $(YPATCHY) 'isajet.car isarun.f $(YCRADLE) isarun.lis .GO'          ;\
  $(RM) $(YCRADLE)

runjet.f:       $(MAKEFILE)                                           ;\
  $(RM) runjet.f                                                      ;\
  $(ECHO) "      PROGRAM RUNJET"                            >runjet.f ;\
  $(ECHO) "      CHARACTER*60 FN"                          >>runjet.f ;\
  $(ECHO) "      READ 1000, FN"                            >>runjet.f ;\
  $(ECHO) "1000  FORMAT(A)"                                >>runjet.f ;\
  $(ECHO) "      PRINT 1020, FN"                           >>runjet.f ;\
  $(ECHO) "1020  FORMAT(1X,'Data file      = ',A)"         >>runjet.f ;\
  $(ECHO) "      OPEN(2,FILE=FN,STATUS='NEW',FORM='UNFORMATTED')" >>runjet.f ;\
  $(ECHO) "      READ 1000, FN"                            >>runjet.f ;\
  $(ECHO) "      PRINT 1030, FN"                           >>runjet.f ;\
  $(ECHO) "1030  FORMAT(1X,'Parameter file = ',A)"         >>runjet.f ;\
  $(ECHO) "      OPEN(3,FILE=FN,STATUS='OLD',FORM='FORMATTED')" >>runjet.f ;\
  $(ECHO) "      READ 1000, FN"                            >>runjet.f ;\
  $(ECHO) "      PRINT 1040, FN"                           >>runjet.f ;\
  $(ECHO) "1040  FORMAT(1X,'Listing file   = ',A)"         >>runjet.f ;\
  $(ECHO) "      OPEN(4,FILE=FN,STATUS='NEW',FORM='FORMATTED')" >>runjet.f ;\
  $(ECHO) "      READ 1000, FN"                            >>runjet.f ;\
  $(ECHO) "      OPEN(1,FILE=FN,STATUS='OLD',FORM='FORMATTED')" >>runjet.f ;\
  $(ECHO) "      CALL ISAJET(-1,2,3,4)"                    >>runjet.f ;\
  $(ECHO) "      STOP"                                     >>runjet.f ;\
  $(ECHO) "      END"                                      >>runjet.f ;\
  cd susylha                                                          ;\
  make                                                                ;\
  cd ..

########## Compile #####################################################

# Note HP compiler cannot compile arbitrarily many files!
libisajet.a:    isajet.f                                              ;\
  $(RM) -r $(MAKETMP)                                                 ;\
  mkdir  $(MAKETMP)                                                   ;\
  cd     $(MAKETMP)                                                   ;\
  $(FSPLIT) ../isajet.f                                               ;\
  $(F77) $(FFLAGS) [a-i,A-I]*.f                                       ;\
  $(F77) $(FFLAGS) [j-r,J-R]*.f                                       ;\
  $(F77) $(FFLAGS) [s,S]*.f                                           ;\
  $(F77) $(FFLAGS) [t-z,T-Z]*.f                                       ;\
  ar r libisajet.a *.o                                                ;\
  $(RANLIB) libisajet.a                                               ;\
  mv libisajet.a ..                                                   ;\
  cd ..

aldata.o:       aldata.f                                              ;\
  $(F77) $(FFLAGS) aldata.f

runjet.o:       runjet.f                                              ;\
  $(F77) $(FFLAGS) runjet.f

isajet.x:       runjet.o aldata.o libisajet.a                         ;\
  $(LINK) -o isajet.x runjet.o aldata.o -L. -lisajet $(LFLAGS)

isasusy.x:      ssrun.f aldata.o libisajet.a                          ;\
  $(F77) $(FFLAGS) ssrun.f ;\
  $(LINK) -o isasusy.x ssrun.o aldata.o -L. -lisajet $(LFLAGS)

isasugra.x:     sugrun.f aldata.o libisajet.a                         ;\
  $(F77) $(FFLAGS) sugrun.f ;\
  $(LINK) -o isasugra.x sugrun.o aldata.o -L. -lisajet $(LIBTOOLS) $(LFLAGS)

libisazeb.a:    isazeb.f                                              ;\
  $(RM) -r $(MAKETMP)                                                 ;\
  mkdir  $(MAKETMP)                                                   ;\
  cd     $(MAKETMP)                                                   ;\
  $(FSPLIT) ../isazeb.f                                               ;\
  $(F77) $(FFLAGS) [a-i,A-I]*.f                                       ;\
  $(F77) $(FFLAGS) [j-r,J-R]*.f                                       ;\
  $(F77) $(FFLAGS) [s,S]*.f                                           ;\
  $(F77) $(FFLAGS) [t-z,T-Z]*.f                                       ;\
  ar r libisazeb.a *.o                                                ;\
  $(RANLIB) libisazeb.a                                               ;\
  mv libisazeb.a ..                                                   ;\
  cd ..

isarun.o:       isarun.f                                              ;\
  $(F77) $(FFLAGS) isarun.f

isazeb.x:       isarun.o aldata.o libisazeb.a                         ;\
  $(LINK) -o isazeb.x isarun.o aldata.o -L. -lisazeb $(LFLAGS)

########## Install/clean ###############################################

clean: ;\
  set nonomatch                                                       ;\
  if (-f aldata.o) mv aldata.o aldata.obj                             ;\
  $(RM) *.f *.o                                                       ;\
  $(RM) *.lis                                                         ;\
  if(-f $(YCRADLE)) $(RM) $(YCRADLE)                                  ;\
  if (-f aldata.obj) mv aldata.obj aldata.o                           ;\
  $(RM) -r $(MAKETMP)

superclean:                                                           ;\
  set nonomatch                                                       ;\
  $(RM) *.f *.o                                                       ;\
  $(RM) *.lis                                                         ;\
  $(RM) $(YCRADLE)                                                    ;\
  $(RM) *.obj                                                         ;\
  $(RM) -r $(MAKETMP)                                                 ;\
  $(RM) isajet.f libisajet.a isajet.x aldata.f \
        isadecay.dat isajet.tex isasusy.x isasugra.x \
        isazeb.x libisazeb.a isarun.o isarun.f


########################################################################
# Isatools: Tools interfaced to Isajet to calculate Delta(a_mu),
# BF(b->s gamma), sigma(Cdm), and Omega h^2 for CDM.
# Requires isared.tar, tar file of CDM code.
########################################################################

isatools:       libisared.a

MAKETOOLS = ./MakeRedDir

libisared.a:  isajet.car isared.tar $(MAKEFILE)                       ;\
  set nonomatch                                                       ;\
  $(RM) -r $(MAKETOOLS)                                               ;\
  mkdir $(MAKETOOLS)                                                  ;\
  cd $(MAKETOOLS)                                                     ;\
  tar xvf ../isared.tar                                               ;\
  $(F77) $(FFLAGS) *.f                                                ;\
  ar r libisared.a *.o                                                ;\
  $(RANLIB) libisared.a                                               ;\
  mv libisared.a ..                                                   ;\
  cd ..

cleantools:                                                           ;\
  set nonomatch                                                       ;\
  $(RM) -r $(MAKETOOLS)
