CDECK  ID>, LOGERR.
      SUBROUTINE LOGERR(IMSG,I,IERR)
C
C          ERROR MESSAGES
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
C          Jet limits
      INTEGER MXLIM
      PARAMETER (MXLIM=8)
      INTEGER MXLX12
      PARAMETER (MXLX12=12*MXLIM)
      COMMON/JETLIM/PMIN(MXLIM),PMAX(MXLIM),PTMIN(MXLIM),PTMAX(MXLIM),
     $YJMIN(MXLIM),YJMAX(MXLIM),PHIMIN(MXLIM),PHIMAX(MXLIM),
     $XJMIN(MXLIM),XJMAX(MXLIM),THMIN(MXLIM),THMAX(MXLIM),
     $SETLMJ(12*MXLIM)
      SAVE /JETLIM/
      COMMON/FIXPAR/FIXP(MXLIM),FIXPT(MXLIM),FIXYJ(MXLIM),
     $FIXPHI(MXLIM),FIXXJ(MXLIM),FIXQM,FIXQT,FIXYW,FIXXW,FIXPHW
      SAVE /FIXPAR/
      COMMON/SGNPAR/CTHS(2,MXLIM),THS(2,MXLIM),YJS(2,MXLIM),XJS(2,MXLIM)
      SAVE /SGNPAR/
      REAL      PMIN,PMAX,PTMIN,PTMAX,YJMIN,YJMAX,PHIMIN,PHIMAX,XJMIN,
     +          XJMAX,THMIN,THMAX,BLIMS(12*MXLIM),CTHS,THS,YJS,XJS
      LOGICAL SETLMJ
      LOGICAL FIXQM,FIXQT,FIXYW,FIXXW,FIXPHW
      LOGICAL FIXP,FIXPT,FIXYJ,FIXPHI,FIXXJ
      EQUIVALENCE(BLIMS(1),PMIN(1))
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
      COMMON/JETPAR/P(3),PT(3),YJ(3),PHI(3),XJ(3),TH(3),CTH(3),STH(3)
     1 ,JETTYP(3),SHAT,THAT,UHAT,QSQ,X1,X2,PBEAM(2)
     2 ,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,JWTYP
     3 ,ALFQSQ,CTHW,STHW,Q0W
     4 ,INITYP(2),ISIGS,PBEAMS(5)
      SAVE /JETPAR/
      INTEGER   JETTYP,JWTYP,INITYP,ISIGS
      REAL      P,PT,YJ,PHI,XJ,TH,CTH,STH,SHAT,THAT,UHAT,QSQ,X1,X2,
     +          PBEAM,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,
     +          ALFQSQ,CTHW,STHW,Q0W,PBEAMS
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
      COMMON/DYLIM/QMIN,QMAX,QTMIN,QTMAX,YWMIN,YWMAX,XWMIN,XWMAX,THWMIN,
     2  THWMAX,PHWMIN,PHWMAX
     3  ,SETLMQ(12)
      SAVE /DYLIM/
      LOGICAL SETLMQ
      EQUIVALENCE(BLIM1(1),QMIN)
      REAL      QMIN,QMAX,QTMIN,QTMAX,YWMIN,YWMAX,XWMIN,XWMAX,THWMIN,
     +          THWMAX,PHWMIN,PHWMAX,BLIM1(12)
      INTEGER MXKEYS
      PARAMETER (MXKEYS=20)
      COMMON/KEYS/IKEYS,KEYON,KEYS(MXKEYS)
      COMMON/XKEYS/REAC
      SAVE /KEYS/,/XKEYS/
      LOGICAL KEYS
      LOGICAL KEYON
      CHARACTER*8 REAC
      INTEGER   IKEYS
      INTEGER MXGOQ,MXGOJ
      PARAMETER (MXGOQ=85,MXGOJ=8)
      COMMON/Q1Q2/GOQ(MXGOQ,MXGOJ),GOALL(MXGOJ),GODY(4),STDDY,
     $GOWW(25,2),ALLWW(2),GOWMOD(25,MXGOJ)
      SAVE /Q1Q2/
      LOGICAL GOQ,GOALL,GODY,STDDY,GOWW,ALLWW,GOWMOD

C
C        ERRORS IN JET PARAMETERS
C
      IERR=IERR+1
      IF(IMSG.EQ.0) WRITE(ITLIS,81)
81    FORMAT(//5X,'DEFAULT LIMITS HAVE BEEN SET')
      IF(IMSG.EQ.1) WRITE(ITLIS,1001) I,PMIN(I),PMAX(I)
 1001 FORMAT(//10X,'BAD LIMITS FOR P(',I2,')=',2E12.4)
      IF(IMSG.EQ.2) WRITE(ITLIS,1002) I,PTMIN(I),PTMAX(I)
 1002 FORMAT(//10X,'BAD LIMITS FOR PT(',I2,')=',2E12.4)
      IF(IMSG.EQ.3) WRITE(ITLIS,1003) I,THMIN(I),THMAX(I)
 1003 FORMAT(//10X,'BAD LIMITS FOR THETA(',I2,')=',2E12.4)
      IF(IMSG.EQ.4) WRITE(ITLIS,1004) I,XJMIN(I),XJMAX(I)
 1004 FORMAT(//10X,'BAD LIMITS FOR X(',I2,')=',2E12.4)
      IF(IMSG.EQ.5) WRITE(ITLIS,1005) I,XJ(I),P(I)
 1005 FORMAT(//5X,'X AND P FOR JET',I2,' ARE INCOMPATIBLE',2E12.4)
      IF(IMSG.EQ.6) WRITE(ITLIS,1006) I,THMIN(I),THMAX(I)
 1006 FORMAT(//10X,'LIMITS FOR THETA MUST BE .GT.0 AND .LT.PI. PRESENT'
     C  ,' LIMITS FOR JET NO.',I3,' ARE',2E12.4)
      IF(IMSG.EQ.7) WRITE(ITLIS,1007) I,XJ(I),X1,X2
 1007 FORMAT(//5X,'FIXED X VALUE FOR JET NO.',I3,' IS',E12.4,2X,
     C  'THIS IS INCOMPATIBLE WITH ALLOWED X LIMITS',2E12.4)
C
C           ERRORS IN W(Z0) PARAMETERS
C
      IF(IMSG.EQ.101) WRITE(ITLIS,901) XW,XWMIN,XWMAX
  901 FORMAT(//5X,'CHOICE OF PARAMETERS GIVES A FIXED XW',E12.4,
     C  ' ,THIS VALUE IS INCOMPATIBLE WITH THE LIMITS',2E12.4)
      IF(IMSG.EQ.102) WRITE(ITLIS,902) YW,YWMIN,YWMAX
  902 FORMAT(//5X,'CHOICE OF PARAMETERS GIVES A FIXED YW',
     C  E12.4,' ,THIS VALUE IS INCOMPATIBLE WITH THE LIMITS ')
      IF(IMSG.EQ.103) WRITE(ITLIS,903) QMW,QMIN,QMAX
  903 FORMAT(//5X,'CHOICE OF PARAMETERS GIVES A FIXED QMW',
     C  E12.4,' ,THIS VALUE IS INCOMPATIBLE WITH THE LIMITS',
     C  E12.4)
      IF(IMSG.EQ.104) WRITE(ITLIS,904) XW,YW,QTW
  904 FORMAT(//5X,'FIXED VALUES FOR XW,YW,AND QTW',3E12.4,
     C  ' ARE UNPHYSICAL')
      IF(IMSG.EQ.105) WRITE(ITLIS,905) QTW,QTMIN,QTMAX
  905 FORMAT(//5X,'CHOICE OF PARAMETERS GIVES A FIXED QTW',E12.4
     C  ,' ,THIS VALUE IS INCOMPATIBLE WITH THE LIMITS',2E12.4)
      IF(IMSG.EQ.106) WRITE(ITLIS,906) XW,YW,QMW
  906 FORMAT(//5X,'FIXED VALUS FOR XW,YW,AND QMW',3E12.4,
     C  ' ARE UNPHYSICAL')
      IF(IMSG.EQ.107) WRITE(ITLIS,907) QTMIN,QTMAX
  907 FORMAT(//5X,'BAD LIMITS FOR QTW',2E12.4)
      IF(IMSG.EQ.108) WRITE(ITLIS,908) QMIN,QMAX
  908 FORMAT(//5X,'BAD LIMITS FOR QMW',2E12.4)
      IF(IMSG.EQ.109) WRITE(ITLIS,909) THWMIN,THWMAX
  909 FORMAT(//5X,'BAD LIMITS FOR THW',2E12.4,2X,' REMEMBER TH MUST',
     C  ' BE IN RADIANS AND LIE BETWEEN 0 AND PI')
      IF(IMSG.EQ.110) WRITE(ITLIS,910) PHWMIN,PHWMAX
  910 FORMAT(//5X,'BAD LIMITS FOR PHW',2E12.4,' ,REMEMBER PHW MUST',
     C  ' BE IN RADIANS AND PHMAX-PHMIN MUST BE LESS THAN 2PI')
      IF(IMSG.EQ.111) WRITE(ITLIS,911) XWMIN,XWMAX
  911 FORMAT(//5X,'BAD LIMITS FOR XW',2E12.4)
      IF(IMSG.EQ.112) WRITE(ITLIS,912) YWMIN,YWMAX
  912 FORMAT(//5X,'BAD LIMITS FOR YW',2E12.4)
      IF(IMSG.EQ.113) WRITE(ITLIS,913)
  913 FORMAT(//5X,'SORRY, BUT YOU CANNOT FIX THETA FOR DRELLYAN EVENTS.'
     C,'  THINK OF SOMETHING ELSE.')
      IF(IMSG.EQ.114) WRITE(ITLIS,914)
  914 FORMAT(//5X,'YOU CANNOT FIX PARAMETERS FOR THE DECAY OF A',
     C  ' DRELL YAN JET')
      IF(IMSG.EQ.115) WRITE(ITLIS,915)
  915 FORMAT(//5X,'YOU CANNOT FIX QTW,QMW,YW AND XW SIMULTANEUOSLY')
C
C       ERRORS IN E+E- PARAMETERS
C
      IF(IMSG.EQ.116)
     1WRITE(ITLIS,631) THMIN(1),THMAX(1),THMIN(2),THMAX(2)
631   FORMAT(//10X,'THETA LIMITS',2E12.4,' FOR JET 1 AND',2E12.4
     C  ,' FOR JET 2 ARE INCOMPATIBLE')
C
      RETURN
      END
