CDECK  ID>, SIGHSS.
      SUBROUTINE SIGHSS
C
C          Compute the integrated MSSM Higgs cross section
C          d(sigma)/d(QMW**2)d(YW)
C          Since SUSY Higgs are always narrow, can use the widths to
C          determine couplings and ignore interference with continuum.
C
C          SIGMA    = cross section summed over quark types allowed by
C                     JETTYPE and WTYPE cards.
C          SIGS(I)  = partial cross section for I1 + I2 --> I3 + I4.
C          INOUT(I) = IOPAK**3*I4 + IOPAK**2*I3 + IOPAK*I2 + I1
C                     using JETTYPE code from LISTSS.
C
C          Ver 7.18: Correct GOQ's and include TBRWW for W/Z modes.
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
      COMMON/QSAVE/QSAVE(29,2)
      SAVE /QSAVE/
      REAL      QSAVE
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
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
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
      COMMON/HCON/ANWWWW(4,4,4),ADWWWW(2,4),AIWWWW(4)
     $,HMASS,HGAM,HGAMS(29),ETAHGG,MATCHH(29),ZSTARS(4,2)
     $,IHTYPE,HGAMSS(85,85)
      SAVE /HCON/
      DOUBLE PRECISION ANWWWW,ADWWWW,AIWWWW
      INTEGER   MATCHH,IHTYPE
      REAL      HMASS,HGAM,HGAMS,ETAHGG,ZSTARS,HGAMSS
C
      REAL X(2)
      REAL AMASS,STRUC
      REAL AM1,AM2,S,T,U,Q2SAVE,YHAT,EY,ANEFF,QMW2,QZW,EHAT,SIG0,SIG,
     $AMW
      INTEGER JT1,JT2,I,J,IH,IQ,I1,I2,JTGL,JTOFF
      EQUIVALENCE (S,SHAT),(T,THAT),(U,UHAT),(X(1),X1)
C
C          Kinematics (identical to Drell-Yan)
C
      QMW2=QMW**2
      QTMW=SQRT(QMW2+QTW**2)
      Q0W=QTMW*COSH(YW)
      QZW=QTMW*SINH(YW)
      QW=SQRT(QZW**2+QTW**2)
      IF(QW.NE.0.) THEN
        CTHW=QZW/QW
        STHW=QTW/QW
        IF(ABS(CTHW).LT.1.) THEN
          THW=ACOS(CTHW)
        ELSE
          CTHW=0.
          STHW=1.
          THW=.5*PI
        ENDIF
      ELSE
        CTHW=0.
        STHW=1.
        THW=.5*PI
      ENDIF
      EHAT=QMW
      SHAT=QMW**2
      QSQ=SHAT
      ANEFF=4.+QSQ/(QSQ+AMASS(5)**2)+QSQ/(QSQ+AMASS(6)**2)
      ALFQSQ=12.*PI/((33.-ANEFF)*ALOG(QSQ/ALAM2))
      Q2SAVE=QSQ
      YHAT=YW
      EY=EXP(YHAT)
      X1=EHAT/ECM*EY
      X2=EHAT/(ECM*EY)
C
C          Initialize
C
      SIGMA=0.
      NSIGS=0
      DO 100 I=1,MXSIGS
100   SIGS(I)=0
      IF(X1.GE.1..OR.X2.GE.1.) RETURN
C
C          Compute structure functions
C
      DO 110 IH=1,2
        DO 120 IQ=1,13
120     QSAVE(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
        DO 130 IQ=14,26
130     QSAVE(IQ,IH)=0.
110   CONTINUE
C
C          gl + gl -> Higgs
C
      JTGL=52
      SIG0=PI*HMASS**2/(8*S**2)*HGAMSS(JTGL,JTGL)*X1*X2*UNITS
     $/((S-HMASS**2)**2+(HMASS*HGAM)**2)
      SIG0=SIG0*QSAVE(1,1)*QSAVE(1,2)
      DO 200 I=1,85
        DO 210 J=1,85
          IF(HGAMSS(I,J).EQ.0) GO TO 210
          IF(.NOT.(GOQ(I,1).AND.GOQ(J,2))) GO TO 210
          SIG=SIG0*HGAMSS(I,J)
C          Include W/Z branching ratios
          IF((I.GE.78.AND.I.LE.80).AND.(J.GE.78.AND.J.LE.80)) THEN
            SIG=SIG*TBRWW(I-76,1)*TBRWW(J-76,2)
          ENDIF
          CALL SIGFIL(SIG,JTGL,JTGL,I,J)
210     CONTINUE
200   CONTINUE
C
C          qk + qb -> Higgs
C
      JTOFF=51
C          Note I1,I2 run over quarks; JT1,JT2,I,J over LISTSS
      DO 300 I1=2,13
        AM1=AMASS(I1/2)
        JT1=I1+JTOFF
        DO 310 I2=2,13
          AM2=AMASS(I2/2)
          JT2=I2+JTOFF
          IF(HGAMSS(JT1,JT2).LE.0) GO TO 310
          SIG0=4*PI*HMASS**2/(9*S**2)*HGAMSS(JT1,JT2)*X1*X2*UNITS
     $    /((S-HMASS**2)**2+(HMASS*HGAM)**2)
          SIG0=SIG0*QSAVE(I1,1)*QSAVE(I2,2)
C          Decay partial cross sections
          DO 320 I=1,85
            DO 330 J=1,85
              IF(HGAMSS(I,J).EQ.0) GO TO 330
              IF(.NOT.(GOQ(I,1).AND.GOQ(J,2))) GO TO 330
              SIG=SIG0*HGAMSS(I,J)
C          Include W/Z branching ratios
              IF((I.GE.78.AND.I.LE.80).AND.(J.GE.78.AND.J.LE.80)) THEN
                SIG=SIG*TBRWW(I-76,1)*TBRWW(J-76,2)
              ENDIF
              CALL SIGFIL(SIG,JT1,JT2,I,J)
330         CONTINUE
320       CONTINUE
310     CONTINUE
300   CONTINUE
C
      RETURN
      END
