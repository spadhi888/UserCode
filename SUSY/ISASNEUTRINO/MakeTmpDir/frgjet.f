CDECK  ID>, FRGJET.
      SUBROUTINE FRGJET(JET)
C
C          Hadronize all partons in /JETSET/ corresponding to jet JET.
C
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
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
      INTEGER MXJETS
      PARAMETER (MXJETS=10)
      COMMON/PJETS/PJETS(5,MXJETS),IDJETS(MXJETS),QWJET(5),IDENTW
     $,PPAIR(5,4),IDPAIR(4),JPAIR(4),NPAIR,IFRAME(MXJETS)
      SAVE /PJETS/
      INTEGER   IDJETS,IDENTW,IDPAIR,JPAIR,NPAIR,IFRAME
      REAL      PJETS,QWJET,PPAIR
      COMMON/PINITS/PINITS(5,2),IDINIT(2)
      SAVE /PINITS/
      INTEGER   IDINIT
      REAL      PINITS
      INTEGER   MXPTCL,IPACK
      PARAMETER (MXPTCL=4000,IPACK=10000)
      COMMON/PARTCL/NPTCL,PPTCL(5,MXPTCL),IORIG(MXPTCL),IDENT(MXPTCL)
     1,IDCAY(MXPTCL)
      SAVE /PARTCL/
      INTEGER   NPTCL,IORIG,IDENT,IDCAY
      REAL      PPTCL
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
      INTEGER   MXJSET,JPACK
      PARAMETER (MXJSET=400,JPACK=1000)
      COMMON/JETSET/NJSET,PJSET(5,MXJSET),JORIG(MXJSET),JTYPE(MXJSET),
     $JDCAY(MXJSET)
      SAVE /JETSET/
      INTEGER   NJSET,JORIG,JTYPE,JDCAY
      REAL      PJSET
      COMMON/JWORK/ZZC(MXJSET),JMATCH(MXJSET),TNEW,P1CM(4),
     1J1,J2,J3,J4,J5,E1CM,E2CM,E3CM,E4CM,E5CM
      SAVE /JWORK/
      LOGICAL TNEW
      EQUIVALENCE (J1,JJ(1)),(E1CM,EE(1))
      INTEGER   JMATCH,J1,J2,J3,J4,J5,JJ(5)
      REAL      ZZC,P1CM,E1CM,E2CM,E3CM,E4CM,E5CM,EE(5)
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
      COMMON/FRAME/FRAME(5,3),N0JETS,N0W,N0PAIR
      SAVE /FRAME/
      INTEGER   N0JETS,N0W,N0PAIR
      REAL      FRAME
C
      REAL ROT(3,3),POLD(5),PNEW(5),PSUM(5)
      REAL CPHI,SPHI,AMSUM,ESUM,PJ,CTHJ,STHJ,PTJ
      INTEGER K,K1,K2,IP,NPLV1,IFAIL,NBEGIN,JET,NFRAG,NFRGMX,JETJ,
     $JTABS,NFIRST,J
C
      DATA PSUM/5*0./
C
C          NFRAG counter protects against possible infinite loop.
C
      NFRAG=0
      NFRGMX=10*MXJSET
201   NBEGIN=NPTCL+1
      NFRAG=NFRAG+1
C
C          Loop over partons
C
      ESUM=0.
      DO 220 J=1,NJSET
        IF(JDCAY(J).NE.0) GO TO 220
        JETJ=JORIG(J)/JPACK
        IF(JETJ.NE.JET) GO TO 220
        ESUM=ESUM+PJSET(4,J)
C
C          Generate Field-Feynman jet for each quark or gluon, or...
C
        JTABS = IABS(JTYPE(J))
        IF(JTABS.LT.10) THEN
          NFIRST=NPTCL+1
          CALL JETGEN(J)
          IF(NPTCL.LT.NFIRST) GO TO 220
C
C          Rotate hadrons to parton direction
C
          PTJ=PJSET(1,J)**2+PJSET(2,J)**2
          PJ=SQRT(PTJ+PJSET(3,J)**2)
          PTJ=SQRT(PTJ)
C          Following is to fix occasional bug on 32-bit machines
          IF(PJ.GT.0.) THEN
            CTHJ=PJSET(3,J)/PJ
            STHJ=PTJ/PJ
          ELSE
            CTHJ=1.
            STHJ=0.
          ENDIF
          IF(PTJ.GT.0.) THEN
            CPHI=PJSET(1,J)/PTJ
            SPHI=PJSET(2,J)/PTJ
          ELSE
            CPHI=SIGN(1.,PJSET(3,J))
            SPHI=0.
          ENDIF
          ROT(1,1)=CPHI*CTHJ
          ROT(2,1)=SPHI*CTHJ
          ROT(3,1)=-STHJ
          ROT(1,2)=-SPHI
          ROT(2,2)=CPHI
          ROT(3,2)=0.
          ROT(1,3)=CPHI*STHJ
          ROT(2,3)=SPHI*STHJ
          ROT(3,3)=CTHJ
          DO 230 IP=NFIRST,NPTCL
            DO 235 K=1,3
              POLD(K)=PPTCL(K,IP)
              PPTCL(K,IP)=0
235         CONTINUE
            DO 240 K1=1,3
            DO 240 K2=1,3
240         PPTCL(K1,IP)=PPTCL(K1,IP)+ROT(K1,K2)*POLD(K2)
230       CONTINUE
C
C          ... hadronize all other partons with delta function.
C
        ELSE
          IF((IABS(JTYPE(J)).EQ.80.OR.IABS(JTYPE(J)).EQ.90).AND.
     $    .NOT.KEYS(2).AND..NOT.KEYS(12)) GO TO 210
          IF(NPTCL.GE.MXPTCL) GO TO 9999
          NPTCL=NPTCL+1
          DO 255 K=1,5
            PPTCL(K,NPTCL)=PJSET(K,J)
255       CONTINUE
          IORIG(NPTCL)=-J
          IDENT(NPTCL)=JTYPE(J)
          IDCAY(NPTCL)=0
        ENDIF
220   CONTINUE
C
C          Sum masses and insert jet label
C
      AMSUM=0.
      DO 260 IP=NBEGIN,NPTCL
        AMSUM=AMSUM+PPTCL(5,IP)
        IORIG(IP)=ISIGN(IABS(IORIG(IP))+IPACK*JET,IORIG(IP))
260   CONTINUE
C
C          Require sum of masses less than jet energy.
C
      IF(AMSUM.GT.ESUM.AND.NBEGIN.NE.NPTCL.AND.NFRAG.LT.NFRGMX) THEN
        NPTCL=NBEGIN-1
        GO TO 201
      ENDIF
C
C          For WPAIR events rescale jet to W mass.
C
      IF((KEYS(6).OR.KEYS(7).OR.KEYS(9).OR.KEYS(10)).AND.JET.LT.10)
     $ THEN
        IF(IABS(JTYPE(JET+N0JETS-1)).LT.80) RETURN
        IF(AMSUM.GE.PJSET(5,JET+N0JETS-1)) THEN
          IF(NFRAG.GT.NFRGMX) RETURN
          NPTCL=NBEGIN-1
          GO TO 201
        ENDIF
        PSUM(4)=PJSET(5,JET+N0JETS-1)
        PSUM(5)=PSUM(4)
        NPLV1=NPTCL
        CALL RESCAL(NBEGIN,NPLV1,PSUM,IFAIL)
      ENDIF
C
210   RETURN
C
C          Error
C
9999  CALL PRTEVT(0)
      WRITE(ITLIS,9998) NPTCL
9998  FORMAT(//' ERROR IN FRGJET ... NPTCL > ',I6)
      RETURN
      END
