CDECK  ID>, NOGOOD.
      LOGICAL FUNCTION NOGOOD(KK)
C
C          Insure proper distribution and check kinematics.
C          Select jet types.
C
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      INTEGER MXKEYS
      PARAMETER (MXKEYS=20)
      COMMON/KEYS/IKEYS,KEYON,KEYS(MXKEYS)
      COMMON/XKEYS/REAC
      SAVE /KEYS/,/XKEYS/
      LOGICAL KEYS
      LOGICAL KEYON
      CHARACTER*8 REAC
      INTEGER   IKEYS
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
      COMMON/WSIG/SIGLLQ
      SAVE /WSIG/
      REAL      SIGLLQ
      COMMON/WGEN/PTGN(3,3),QGEN(3,3),PTSEL(3),QSEL(3),SIGSL(3),NKL,NKH
     1,EMSQ,EMGAM,KSEL,QSELWT(3)
      SAVE /WGEN/
      INTEGER   NKL,NKH,KSEL
      REAL      PTGN,QGEN,PTSEL,QSEL,SIGSL,EMSQ,EMGAM,QSELWT
      COMMON/DYLIM/QMIN,QMAX,QTMIN,QTMAX,YWMIN,YWMAX,XWMIN,XWMAX,THWMIN,
     2  THWMAX,PHWMIN,PHWMAX
     3  ,SETLMQ(12)
      SAVE /DYLIM/
      LOGICAL SETLMQ
      EQUIVALENCE(BLIM1(1),QMIN)
      REAL      QMIN,QMAX,QTMIN,QTMAX,YWMIN,YWMAX,XWMIN,XWMAX,THWMIN,
     +          THWMAX,PHWMIN,PHWMAX,BLIM1(12)
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
      INTEGER   MXSIGS,IOPAK
      PARAMETER (MXSIGS=3000,IOPAK=100)
      COMMON/JETSIG/SIGMA,SIGS(MXSIGS),NSIGS,INOUT(MXSIGS),SIGEVT
      SAVE /JETSIG/
      INTEGER   NSIGS,INOUT
      REAL      SIGMA,SIGS,SIGEVT
      COMMON/PTPAR/PTFUN1,PTFUN2,PTGEN1,PTGEN2,PTGEN3,SIGMAX
      SAVE /PTPAR/
      REAL      PTFUN1,PTFUN2,PTGEN1,PTGEN2,PTGEN3,SIGMAX
      COMMON/HCON/ANWWWW(4,4,4),ADWWWW(2,4),AIWWWW(4)
     $,HMASS,HGAM,HGAMS(29),ETAHGG,MATCHH(29),ZSTARS(4,2)
     $,IHTYPE,HGAMSS(85,85)
      SAVE /HCON/
      DOUBLE PRECISION ANWWWW,ADWWWW,AIWWWW
      INTEGER   MATCHH,IHTYPE
      REAL      HMASS,HGAM,HGAMS,ETAHGG,ZSTARS,HGAMSS
      COMMON/XMSSM/GOMSSM,GOSUG,GOGMSB,GOAMSB,AL3UNI,GOMMAM,GOHCAM
     $,XGLSS,XMUSS,XHASS,XTBSS
     $,XQ1SS,XDRSS,XURSS,XL1SS,XERSS
     $,XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
     $,XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
     $,XM1SS,XM2SS,XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
     $,XLAMGM,XMESGM,XN5GM,XCMGV,XMGVTO
     $,XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
     $,XMN3NR,XMAJNR,XANSS,XNRSS,XSBCS,
     $XCQAM,XCDAM,XCUAM,XCLAM,XCEAM,XCHDAM,XCHUAM,
     $XL1AM,XL2AM,XL3AM
      SAVE /XMSSM/
      REAL XGLSS,XMUSS,XHASS,XTBSS
     $,XQ1SS,XDRSS,XURSS,XL1SS,XERSS
     $,XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS
     $,XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS
     $,XM1SS,XM2SS
     $,XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU
     $,XLAMGM,XMESGM,XN5GM,XCMGV,XMGVTO
     $,XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM
     $,XMN3NR,XMAJNR,XANSS,XNRSS,XSBCS,
     $XCQAM,XCDAM,XCUAM,XCLAM,XCEAM,XCHDAM,XCHUAM,
     $XL1AM,XL2AM,XL3AM
      LOGICAL GOMSSM,GOSUG,GOGMSB,GOAMSB,AL3UNI,GOMMAM,GOHCAM
C
      REAL RANF,SIGINV,SUM,TRY,BRANCH
      INTEGER KK,I,II,K,IFL
C
      NOGOOD=.TRUE.
      GO TO (1,2,3,4,5,6),KK
C
C          TWOJET, SUPERSYM, WPAIR or PHOTON events
C
1     CONTINUE
      IF(KEYS(1)) THEN
        CALL SIGQCD
      ELSEIF(KEYS(5)) THEN
        CALL SIGSSY
      ELSEIF(KEYS(6)) THEN
        CALL SIGWW
      ELSEIF(KEYS(8)) THEN
        CALL SIGGAM
      ELSEIF(KEYS(10)) THEN
        CALL SIGWH
      ENDIF
      IF(SIGMA.LE.0) RETURN
      IF(SIGMAX*RANF().GT.SIGMA) RETURN
      NOGOOD=.FALSE.
      SIGINV=1./SIGMA
      SUM=0.
      TRY=RANF()
      DO 100 I=1,NSIGS
        SUM=SUM+SIGS(I)*SIGINV
        IF(SUM.LT.TRY) GO TO 100
C          Find reaction
        ISIGS=I
        SIGEVT=SIGS(ISIGS)
        II=INOUT(I)
        DO 110 K=1,2
        INITYP(K)=MOD(II,IOPAK)
110     II=II/IOPAK
        DO 120 K=1,2
        JETTYP(K)=MOD(II,IOPAK)
120     II=II/IOPAK
        RETURN
100   CONTINUE
      RETURN
C
C          DRELLYAN events--test of SIGDY
C
2     CONTINUE
      IF(KEYS(3)) THEN
        CALL SIGDY
      ELSEIF(KEYS(7).AND..NOT.GOMSSM) THEN
        CALL SIGH
      ELSEIF(KEYS(7).AND.GOMSSM) THEN
        CALL SIGHSS
      ELSEIF(KEYS(9)) THEN
        CALL SIGTC
      ELSEIF(KEYS(11)) THEN
        CALL SIGKKG
      ENDIF
      IF(SIGMA.LE.0.) RETURN
      IF(SIGSL(KSEL)*RANF().GT.SIGMA) RETURN
      NOGOOD=.FALSE.
      SIGINV=1./SIGMA
      SUM=0.
      TRY=RANF()
C          Find reaction.
      DO 200 I=1,NSIGS
        SUM=SUM+SIGS(I)*SIGINV
        IF(SUM.LT.TRY) GO TO 200
        ISIGS=I
        SIGEVT=SIGS(ISIGS)
        GO TO 210
200   CONTINUE
C          Unpack INOUT to find JETTYP and INITYP
210   IF(KEYS(3).OR.KEYS(11)) THEN
        II=INOUT(I)
        DO 220 K=1,2
        INITYP(K)=MOD(II,IOPAK)
220     II=II/IOPAK
        JWTYP=MOD(II,IOPAK)
        II=II/IOPAK
        JETTYP(3)=MOD(II,IOPAK)
      ELSEIF(KEYS(7).OR.KEYS(9)) THEN
        II=INOUT(ISIGS)
        DO 230 I=1,2
        INITYP(I)=MOD(II,IOPAK)
230     II=II/IOPAK
        DO 240 I=1,2
        JETTYP(I)=MOD(II,IOPAK)
240     II=II/IOPAK
      ENDIF
      RETURN
C
C          DRELLYAN events--test of SIGDY2
C
3     CONTINUE
      IF(KEYS(3)) THEN
        CALL SIGDY2
        IFL=JETTYP(1)/2
        BRANCH=(AQ(IFL,JWTYP)**2+BQ(IFL,JWTYP)**2)/COUT(JWTYP)
      ELSEIF(KEYS(7).AND..NOT.GOMSSM) THEN
        CALL SIGH2
        BRANCH=1.
      ELSEIF(KEYS(7).AND.GOMSSM) THEN
        SIGLLQ=SIGMA/(4*PI)
        NOGOOD=.FALSE.
        RETURN
      ELSEIF(KEYS(9)) THEN
        CALL SIGTC2
        BRANCH=1.
      ENDIF
      IF(SIGLLQ.GT.SIGS(ISIGS)*BRANCH*3.*RANF()/(4.*PI))
     1NOGOOD=.FALSE.
      RETURN
C
C          DRELLYAN events--test of kinematics
C
4     CONTINUE
      DO 400 I=1,2
        IF(P(I).LT.PMIN(I).OR.P(I).GT.PMAX(I)) GO TO 410
        IF(PT(I).LT.PTMIN(I).OR.PT(I).GT.PTMAX(I)) GO TO 410
        IF(YJ(I).LT.YJMIN(I).OR.YJ(I).GT.YJMAX(I)) GO TO 410
        IF(PHI(I).LT.PHIMIN(I).OR.PHI(I).GT.PHIMAX(I)) GO TO 410
400   CONTINUE
      NOGOOD=.FALSE.
410   RETURN
C
5     CONTINUE
6     CONTINUE
      RETURN
C
      END
