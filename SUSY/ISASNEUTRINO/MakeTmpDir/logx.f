CDECK  ID>, LOGX.
      LOGICAL FUNCTION LOGX(IERR)
C
C         SET AND CHECK LIMITS FOR JET FEYNMAN X
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
      DATA UNDEF/-.9E9/
C
      HALFPI=PI/2.
      LOGX=.TRUE.
C
      DO 40 I=1,NJET
      FIXXJ(I)=.FALSE.
      IF(FIXYJ(I).AND.(FIXP(I).OR.FIXPT(I)))FIXXJ(I)=.TRUE.
      IF(FIXXJ(I)) GOTO 40
C
      IF(XJMIN(I).LT.UNDEF.AND.XJMAX(I).LT.UNDEF) THEN
        XJMAX(I)=1.0
        XJMIN(I)=-1.0
      ENDIF
C
      IF(XJMAX(I).LT.UNDEF) FIXXJ(I)=.TRUE.
      IF(FIXXJ(I)) XJMAX(I)=XJMIN(I)
C
      IF(.NOT.FIXXJ(I)) THEN
        IF(THMIN(I).LT.HALFPI) X1=PMAX(I)*COS(THMIN(I))/HALFE
        IF(THMIN(I).GE.HALFPI) X1=PMIN(I)*COS(THMIN(I))/HALFE
        IF(THMAX(I).GT.HALFPI) X2=PMAX(I)*COS(THMAX(I))/HALFE
        IF(THMAX(I).LT.HALFPI) X2=PMIN(I)*COS(THMAX(I))/HALFE
        IF(X1.LT.XJMAX(I)) XJMAX(I)=X1
        IF(X2.GT.XJMIN(I)) XJMIN(I)=X2
      ELSE
C
        XJ(I)=XJMIN(I)
C
        IF(FIXP(I)) THEN
          CTH(I)=XJ(I)*HALFE/P(I)
          IF(ABS(CTH(I)).LE.1.0) THEN
            STH(I)=SQRT(1.-CTH(I)**2)
            TH(I)=ATAN2(STH(I),CTH(I))
            YJ(I)=-ALOG(TAN(TH(I)/2.))
            FIXYJ(I)=.TRUE.
            PT(I)=P(I)*STH(I)
            FIXPT(I)=.TRUE.
            YJMIN(I)=YJ(I)
            YJMAX(I)=YJ(I)
            PTMIN(I)=PT(I)
            PTMAX(I)=PT(I)
          ELSE
            LOGX=.FALSE.
            CALL LOGERR(5,I,IERR)
          ENDIF
        ENDIF
C
        IF(FIXPT(I)) THEN
          TH(I)=ATAN(PT(I)/XJ(I)/HALFE)
          FIXYJ(I)=.TRUE.
          YJ(I)=-ALOG(TAN(TH(I)/2.))
          CTH(I)=COS(TH(I))
          STH(I)=SIN(TH(I))
          P(I)=PT(I)/STH(I)
          FIXP(I)=.TRUE.
          YJMIN(I)=YJ(I)
          YJMAX(I)=YJ(I)
          PMAX(I)=P(I)
          PMIN(I)=P(I)
        ENDIF
C
        IF(FIXYJ(I)) THEN
          FIXPT(I)=.TRUE.
          P(I)=XJ(I)*HALFE/CTH(I)
          PT(I)=P(I)*STH(I)
          FIXP(I)=.TRUE.
          PTMIN(I)=PT(I)
          PTMAX(I)=PT(I)
          PMAX(I)=P(I)
          PMIN(I)=P(I)
        ENDIF
C
      ENDIF
C
   40 CONTINUE
C
      RETURN
      END
