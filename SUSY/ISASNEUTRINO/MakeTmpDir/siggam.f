CDECK  ID>, SIGGAM.
      SUBROUTINE SIGGAM
C
C          Compute D(SIGMA)/D(PT**2)D(Y1)D(Y2) for gamma + jet and
C          gamma + gamma.
C
C          SIGMA    = cross section summed over quark types allowed by
C                     JETTYPE card.
C          SIGS(I)  = partial cross section for I1 + I2 --> I3 + I4.
C          INOUT(I) = IOPAK**3*I4 + IOPAK**2*I3 + IOPAK*I2 + I1
C                     using JETTYPE code.
C
C          Cross sections from Berger, Bratten, and Field, Nucl. Phys.
C          B239, 52 (1984), Table 2. Masses are neglected.
C
C
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/QCDPAR/ALAM,ALAM2,CUTJET,ISTRUC
      SAVE /QCDPAR/
      INTEGER   ISTRUC
      REAL      ALAM,ALAM2,CUTJET
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
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
      INTEGER MXGOQ,MXGOJ
      PARAMETER (MXGOQ=85,MXGOJ=8)
      COMMON/Q1Q2/GOQ(MXGOQ,MXGOJ),GOALL(MXGOJ),GODY(4),STDDY,
     $GOWW(25,2),ALLWW(2),GOWMOD(25,MXGOJ)
      SAVE /Q1Q2/
      LOGICAL GOQ,GOALL,GODY,STDDY,GOWW,ALLWW,GOWMOD
      INTEGER   MXSIGS,IOPAK
      PARAMETER (MXSIGS=3000,IOPAK=100)
      COMMON/JETSIG/SIGMA,SIGS(MXSIGS),NSIGS,INOUT(MXSIGS),SIGEVT
      SAVE /JETSIG/
      INTEGER   NSIGS,INOUT
      REAL      SIGMA,SIGS,SIGEVT
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
      COMMON/WCON/SIN2W,WMASS(4),WGAM(4),AQ(12,4),BQ(12,4),COUT(4),
     1MATCH(25,4),WCBR(25,4),CUTOFF,CUTPOW,TBRWW(4,2),RBRWW(12,4,2),EZ,
     2AQDP(12,4),BQDP(12,4),EZDP,WFUDGE
      SAVE /WCON/
      DOUBLE PRECISION AQDP,BQDP,EZDP
      INTEGER   MATCH
      REAL      SIN2W,WMASS,WGAM,AQ,BQ,COUT,WCBR,CUTOFF,CUTPOW,TBRWW,
     +          RBRWW,EZ,WFUDGE
      COMMON/WCON2/CUMWBR(25,3)
      REAL CUMWBR
C
      REAL BBF1,BBF2,BBF3,S,T,U,FJAC,STRUC,SIG0,SIG,BBF3TU,BBF3UT
      INTEGER I,IH,IQ,IFL
      REAL X(2),QSAVE(13,2)
      INTEGER LISTJ(13)
      EQUIVALENCE (X(1),X1),(S,SHAT),(T,THAT),(U,UHAT)
      DATA LISTJ/9,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6/
C
C          Cross sections with couplings and Jacobean removed.
      BBF1(S,T,U)=8./9.*(U/T+T/U)
      BBF2(S,T,U)=2./3.*(U/T+T/U)
      BBF3(S,T,U)=-1./3.*(U/S+S/U)
C
C          Initialize cross sections.
C
      SIGMA=0.
      NSIGS=0
      DO 100 I=1,MXSIGS
      SIGS(I)=0.
100   CONTINUE
C
C          Kinematics and structure functions for CH and lighter quarks
C
      CALL TWOKIN(0.,0.,0.,0.)
      FJAC=SHAT/SCM*UNITS*PI/SHAT**2
      IF(X1.GE.1.0.OR.X2.GE.1.0) RETURN
      DO 110 IH=1,2
      DO 110 IQ=1,9
        QSAVE(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
110   CONTINUE
C
C          Compute cross sections summed over all types allowed by
C          JETTYPE card.
C
      IF(.NOT.(GOQ(26,1).OR.GOQ(26,2))) RETURN
C
C          Gluon-photon
C
      IF((GOQ(1,1).AND.GOQ(26,2)).OR.(GOQ(26,1).AND.GOQ(1,2))) THEN
        SIG0=.5*FJAC*ALFQSQ*ALFA*BBF1(S,T,U)
        DO 210 I=1,4
          IFL=LISTJ(2*I)
          SIG=SIG0*AQ(IFL,1)**2*QSAVE(2*I,1)*QSAVE(2*I+1,2)
          IF(GOQ(26,1).AND.GOQ(1,2)) CALL SIGFIL(SIG,2*I,2*I+1,26,1)
          IF(GOQ(1,1).AND.GOQ(26,2)) CALL SIGFIL(SIG,2*I,2*I+1,1,26)
          SIG=SIG0*AQ(IFL,1)**2*QSAVE(2*I+1,1)*QSAVE(2*I,2)
          IF(GOQ(26,1).AND.GOQ(1,2)) CALL SIGFIL(SIG,2*I+1,2*I,26,1)
          IF(GOQ(1,1).AND.GOQ(26,2)) CALL SIGFIL(SIG,2*I+1,2*I,1,26)
210     CONTINUE
      ENDIF
C
C          Photon-photon
C
      IF(GOQ(26,1).AND.GOQ(26,2)) THEN
        SIG0=.5*FJAC*ALFA**2*BBF2(S,T,U)
        DO 220 I=1,4
          IFL=LISTJ(2*I)
          SIG=SIG0*AQ(IFL,1)**4*QSAVE(2*I,1)*QSAVE(2*I+1,2)
          CALL SIGFIL(SIG,2*I,2*I+1,26,26)
          SIG=SIG0*AQ(IFL,1)**4*QSAVE(2*I+1,1)*QSAVE(2*I,2)
          CALL SIGFIL(SIG,2*I+1,2*I,26,26)
220     CONTINUE
      ENDIF
C
C          Quark-photon
C
      BBF3TU=.5*FJAC*ALFA*ALFQSQ*BBF3(S,T,U)
      BBF3UT=.5*FJAC*ALFA*ALFQSQ*BBF3(S,U,T)
      DO 230 I=2,9
        IFL=IABS(LISTJ(I))
        IF(GOQ(26,1).AND.GOQ(I,2)) THEN
          SIG=BBF3TU*AQ(IFL,1)**2*QSAVE(I,1)*QSAVE(1,2)
          CALL SIGFIL(SIG,I,1,26,I)
          SIG=BBF3UT*AQ(IFL,1)**2*QSAVE(1,1)*QSAVE(I,2)
          CALL SIGFIL(SIG,1,I,26,I)
        ENDIF
        IF(GOQ(I,1).AND.GOQ(26,2)) THEN
          SIG=BBF3UT*AQ(IFL,1)**2*QSAVE(I,1)*QSAVE(1,2)
          CALL SIGFIL(SIG,I,1,I,26)
          SIG=BBF3TU*AQ(IFL,1)**2*QSAVE(1,1)*QSAVE(I,2)
          CALL SIGFIL(SIG,1,I,I,26)
        ENDIF
230   CONTINUE
C
      RETURN
      END
