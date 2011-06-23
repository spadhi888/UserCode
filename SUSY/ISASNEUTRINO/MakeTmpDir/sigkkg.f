CDECK  ID>, SIGKKG.
      SUBROUTINE SIGKKG
C
C          Compute the KK graviton direct production cross-section
C          d(sigma)/d(m**2)d(pT**2)d(y3)d(y4)
C          X-sections: G.F.Giudice et al.  hep-ph/9811291
C          Kinematics: sigdy.car (Drell-Yan + jet)
C
      IMPLICIT NONE
C
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
      COMMON/NODCAY/NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOGRAV,
     $NOB,NOTAU
      LOGICAL NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOGRAV,
     $NOB,NOTAU
      SAVE /NODCAY/
C          KKGravity common
      COMMON/KKGRAV/NEXTRAD,MASSD,KKGSD,SURFD,UVCUT
      INTEGER NEXTRAD
      REAL    MASSD,KKGSD,SURFD
      LOGICAL UVCUT
      SAVE /KKGRAV/
C
      REAL X(2)
      REAL Z,S,T,U,QMW2,QZW,EHAT,Q2SAVE,YHAT,EY,P3Z,P1,P2,AMASS,ANEFF,
     $SIG0,DENOM,QT2CUT,SIGT,SIGU,FAC,PROP,FACTOR,SIG,AMT,AMT2,SWT,
     $P1WT,P2WT,X1WT,X2WT,TWT,UWT,Q2,QFCN,STRUC,XX,ACOSH,ATANH,P2M,P1M
      REAL AMI2,AMF2,EFWT
      REAL AJLWT,AJLZT1,AJLZT2,A2,A2B2,QQ,TM2
      INTEGER I,IQ,IH,IQ1,IFL,IQ2,IW,IQ3
      INTEGER NZERO(4)
      REAL AMFAC(13)
      INTEGER NUTYP(25)
      INTEGER IFL1,IFL2
      REAL TERM
      REAL KKGF1,KKGF2,KKGF3,SIG1,SIG2,XG,YG,F1,F2T,F2U,F3
      EQUIVALENCE (S,SHAT),(T,THAT),(U,UHAT),(X(1),X1)
C          Electric charge:
      REAL CHARGE
      EXTERNAL CHARGE
C
C          Kinematics: (Drell-Yan plus jet)
C
      QMW2=QMW**2
      QTMW=SQRT(QMW2+QTW**2)
      Q0W=QTMW*COSH(YW)
      QZW=QTMW*SINH(YW)
      QW=SQRT(QZW**2+QTW**2)
C          Protect against errors
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
c          Drell-Yan plus jet
      P3Z=P(3)*CTH(3)
      SHAT=QMW2+2.*Q0W*P(3)-2.*QZW*P3Z+2.*PT(3)**2
      P1=.5*(P(3)+P3Z+Q0W+QZW)
      P2=.5*(P(3)-P3Z+Q0W-QZW)
      X1=P1/HALFE
      X2=P2/HALFE
      THAT=-2.*P1*(P(3)-P3Z)
      UHAT=-2.*P2*(P(3)+P3Z)
      QSQ=QTW**2
      QSQ=AMAX1(QSQ,4.)
      ANEFF=4.+QSQ/(QSQ+AMASS(5)**2)+QSQ/(QSQ+AMASS(6)**2)
      ALFQSQ=12.*PI/((33.-2.*ANEFF)*ALOG(QSQ/ALAM2))
      Q2SAVE=QSQ
      QSQ=SHAT
C
C          Initialize
C
      SIGMA=0.
      NSIGS=0
      DO 100 I=1,MXSIGS
        SIGS(I)=0.
100   CONTINUE
      IF(X1.GE.1..OR.X2.GE.1.) RETURN
C
C          Structure functions
C
      DO 110 IH=1,2
        DO 120 IQ=1,11
          QSAVE(IQ,IH)=STRUC(X(IH),QSQ,IQ,IDIN(IH))/X(IH)
120     CONTINUE
        QSAVE(12,IH)=0
        QSAVE(13,IH)=0
110   CONTINUE
      QSQ=Q2SAVE
C
      IF((THAT/SHAT).EQ.0.) RETURN
      IF(ABS(THAT/SHAT+1).LT.1.E-06) RETURN
      F1=KKGF1(SHAT,THAT,QMW2)
      F2T=KKGF2(SHAT,THAT,QMW2)
      F2U=KKGF2(SHAT,UHAT,QMW2)
      F3=KKGF3(SHAT,THAT,QMW2)
      IF(F1.LE.0.OR.F2T.LE.0.OR.F2U.LE.0.OR.F3.LE.0) RETURN
C
      SIG0=UNITS*0.5*KKGSD*ALFQSQ*QMW**(NEXTRAD-2)/SCM
C
C          Jet 3 = gamma:
C
      IF(GOQ(26,3)) THEN
        SIG1=UNITS*0.5*KKGSD*ALFA*QMW**(NEXTRAD-2)/SCM
        SIG1=SIG1*F1/48.0
C          qk + qb --> gamma + KKG
        DO 410 IFL=1,5
          IQ1=2*IFL
          IQ2=IQ1+1
          SIG2=SIG1*ABS(CHARGE(IFL))
          SIG=SIG2*QSAVE(IQ1,1)*QSAVE(IQ2,2)
          IF(UVCUT.AND.SHAT.GE.(MASSD**2)) SIG=SIG*(MASSD**2/SHAT)**2
          CALL SIGFIL(SIG,IQ1,IQ2,5,26)
          SIG=SIG2*QSAVE(IQ2,1)*QSAVE(IQ1,2)
          IF(UVCUT.AND.SHAT.GE.(MASSD**2)) SIG=SIG*(MASSD**2/SHAT)**2
          CALL SIGFIL(SIG,IQ2,IQ1,5,26)
410     CONTINUE
      ENDIF
C
C          Jet 3 = gluon:
C
      IF(GOQ(1,3)) THEN
        SIG1=SIG0*F1/36.0
C          qk + qb --> gl + KKG
        DO 210 IFL=1,5
          IQ1=2*IFL
          IQ2=IQ1+1
          SIG=SIG1*QSAVE(IQ1,1)*QSAVE(IQ2,2)
          IF(UVCUT.AND.SHAT.GE.(MASSD**2)) SIG=SIG*(MASSD**2/SHAT)**2
          CALL SIGFIL(SIG,IQ1,IQ2,5,1)
          SIG=SIG1*QSAVE(IQ2,1)*QSAVE(IQ1,2)
          IF(UVCUT.AND.SHAT.GE.(MASSD**2)) SIG=SIG*(MASSD**2/SHAT)**2
          CALL SIGFIL(SIG,IQ1,IQ1,5,1)
210     CONTINUE
C          gl + gl --> gl + KKG
        SIG1=SIG0*F3*3.0/16.0
        SIG=SIG1*QSAVE(1,1)*QSAVE(1,2)
          IF(UVCUT.AND.SHAT.GE.(MASSD**2)) SIG=SIG*(MASSD**2/SHAT)**2
        CALL SIGFIL(SIG,1,1,5,1)
      ENDIF
C
C          Jet 3 = quark:
C
      SIGT=SIG0*F2T/96.0
      SIGU=SIG0*F2U/96.0
C          qk + gl --> qk + KKG
      DO 310 IQ1=2,11
        IQ3=IQ1
        IF(GOQ(IQ3,3)) THEN
          SIG=SIGU*QSAVE(IQ1,1)*QSAVE(1,2)
          IF(UVCUT.AND.SHAT.GE.(MASSD**2)) SIG=SIG*(MASSD**2/SHAT)**2
          CALL SIGFIL(SIG,IQ1,1,5,IQ3)
          SIG=SIGT*QSAVE(IQ1,2)*QSAVE(1,1)
          IF(UVCUT.AND.SHAT.GE.(MASSD**2)) SIG=SIG*(MASSD**2/SHAT)**2
          CALL SIGFIL(SIG,1,IQ1,5,IQ3)
        ENDIF
310   CONTINUE
C
      RETURN
      END
