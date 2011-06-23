CDECK  ID>, SETDKY.
      SUBROUTINE SETDKY(LPRINT)
C
C          Read in decay table from tape ITDKY and set up /DKYTAB/.
C          Then append forced decay modes and set LOOK to negative
C          number pointing to LOOK2, which points to table.
C          Forced decays for antiparticles are stored in conjugated
C          form so that DECAY can always conjugate them.
C
C          Logical flag LPRINT controls printing of table.
C
C          Ver 7.41: Check version of decay table. Also read matrix
C          element flags and save in MELEM:
C          MELEM=0: Phase space
C          MELEM=1: Dalitz decay
C          MELEM=2: omega/phi decay
C          MELEM=3: V-A
C          MELEM=4: V-A plus W propagator (for top)
C          MELEM=5: tau -> ell nu nu
C          MELEM=6: tau -> nu pi/K
C          MELEM=7: tau -> nu rho/a1
C          MELEM=8: tau -> tau (for NOTAU)
C          MELEM=9: H -> W f fbar
C
C          Ver 7.52: add NOB and NOTAU flags
C
      IMPLICIT NONE
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      INTEGER   MXFORC
      PARAMETER (MXFORC=40)
      COMMON/FORCE/NFORCE,IFORCE(MXFORC),MFORCE(5,MXFORC)
     $,LOOK2(2,MXFORC),LOOKST(MXFORC),MEFORC(MXFORC)
      SAVE /FORCE/
      INTEGER   NFORCE,IFORCE,MFORCE,LOOK2,LOOKST,MEFORC
C          LOOK must be dimensioned to the maximum value of INDEX.
      INTEGER   MXLOOK
      PARAMETER (MXLOOK=500)
      INTEGER   MXDKY
      PARAMETER (MXDKY=3000)
      COMMON/DKYTAB/LOOK(MXLOOK),CBR(MXDKY),MODE(5,MXDKY),MELEM(MXDKY)
      SAVE /DKYTAB/
      INTEGER   LOOK,MODE,MELEM
      REAL      CBR
      COMMON/NODCAY/NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOGRAV,
     $NOB,NOTAU
      LOGICAL NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOGRAV,
     $NOB,NOTAU
      SAVE /NODCAY/
C          MXSS         =  maximum number of modes
C          NSSMOD       = number of modes
C          ISSMOD       = initial particle
C          JSSMOD       = final particles
C          GSSMOD       = width
C          BSSMOD       = branching ratio
C          MSSMOD       = decay matrix element pointer
C          LSSMOD       = logical flag used internally by SSME3
      INTEGER MXSS
      PARAMETER (MXSS=1000)
      COMMON/SSMODE/NSSMOD,ISSMOD(MXSS),JSSMOD(5,MXSS),GSSMOD(MXSS)
     $,BSSMOD(MXSS),MSSMOD(MXSS),LSSMOD
      INTEGER NSSMOD,ISSMOD,JSSMOD,MSSMOD
      REAL GSSMOD,BSSMOD
      LOGICAL LSSMOD
      SAVE /SSMODE/
C          SM ident code definitions. These are standard ISAJET but
C          can be changed.
      INTEGER IDUP,IDDN,IDST,IDCH,IDBT,IDTP
      INTEGER IDNE,IDE,IDNM,IDMU,IDNT,IDTAU
      INTEGER IDGL,IDGM,IDW,IDZ,IDH
      PARAMETER (IDUP=1,IDDN=2,IDST=3,IDCH=4,IDBT=5,IDTP=6)
      PARAMETER (IDNE=11,IDE=12,IDNM=13,IDMU=14,IDNT=15,IDTAU=16)
      PARAMETER (IDGL=9,IDGM=10,IDW=80,IDZ=90,IDH=81)
C          SUSY ident code definitions. They are chosen to be similar
C          to those in versions < 6.50 but may be changed.
      INTEGER ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1
      INTEGER ISNEL,ISEL,ISNML,ISMUL,ISNTL,ISTAU1
      INTEGER ISUPR,ISDNR,ISSTR,ISCHR,ISBT2,ISTP2
      INTEGER ISNER,ISER,ISNMR,ISMUR,ISNTR,ISTAU2
      INTEGER ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,ISGL
      INTEGER ISHL,ISHH,ISHA,ISHC
      INTEGER ISGRAV
      INTEGER IDTAUL,IDTAUR
      PARAMETER (ISUPL=21,ISDNL=22,ISSTL=23,ISCHL=24,ISBT1=25,ISTP1=26)
      PARAMETER (ISNEL=31,ISEL=32,ISNML=33,ISMUL=34,ISNTL=35,ISTAU1=36)
      PARAMETER (ISUPR=41,ISDNR=42,ISSTR=43,ISCHR=44,ISBT2=45,ISTP2=46)
      PARAMETER (ISNER=51,ISER=52,ISNMR=53,ISMUR=54,ISNTR=55,ISTAU2=56)
      PARAMETER (ISGL=29)
      PARAMETER (ISZ1=30,ISZ2=40,ISZ3=50,ISZ4=60,ISW1=39,ISW2=49)
      PARAMETER (ISHL=82,ISHH=83,ISHA=84,ISHC=86)
      PARAMETER (ISGRAV=91)
      PARAMETER (IDTAUL=10016,IDTAUR=20016)
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
      INTEGER MXKEYS
      PARAMETER (MXKEYS=20)
      COMMON/KEYS/IKEYS,KEYON,KEYS(MXKEYS)
      COMMON/XKEYS/REAC
      SAVE /KEYS/,/XKEYS/
      LOGICAL KEYS
      LOGICAL KEYON
      CHARACTER*8 REAC
      INTEGER   IKEYS
C
      INTEGER IMODE(6),LOOP,IOLD,I,IRES,ITYPE,K,J,IPOINT
      INTEGER IFL1,IFL2,IFL3,JSPIN,INDEX,ID1,IDANTI,KTYPE,IRES2
      REAL    BR
      CHARACTER*8 LABEL,LMODE(6),LRES
      CHARACTER*8 IBLANK,LREAD(10),IQUIT
      LOGICAL LPRINT
      INTEGER NOUT,NTHAD
      PARAMETER (NOUT=33)
      PARAMETER (NTHAD=12)
      INTEGER IDOUT(NOUT),ITHAD(NTHAD),IDUMMY(5),MEOUT
      REAL SUMBR,SUMBR2,SUMGAM
      CHARACTER*40,V,VOLD,VISAJE
C
      DATA IDOUT/
     $IDTP,ISGL,ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1,ISUPR,ISDNR,
     $ISSTR,ISCHR,ISBT2,ISTP2,ISEL,ISMUL,ISTAU1,ISNEL,ISNML,ISNTL,
     $ISER,ISMUR,ISTAU2,ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,
     $ISHL,ISHH,ISHA,ISHC/
      DATA IQUIT/'////'/,IBLANK/' '/
      DATA ITHAD/-160,-260,-360,
     $  1160,1260,2260,2160,1360,2360,3160,3260,3360/
C
C          Print header for table.
C
      IF(LPRINT) WRITE(ITLIS,10)
10    FORMAT('1',30('*')/' *',28X,'*'/
     1' *',5X,'ISAJET DECAY TABLE',5X,'*'/
     2' *',28X,'*'/' ',30('*')//
     33X,'PART',16X,'DECAY MODE',16X,'CUM BR',10X,'IDENT',18X,
     4'DECAY IDENT'/)
C
C          Initialize. LOOP is the decay mode counter.
C
      LOOP=0
      IOLD=0
      DO 100 I=1,MXLOOK
        LOOK(I)=0
100   CONTINUE
      DO 110 I=1,MXFORC
        LOOK2(1,I)=0
        LOOK2(2,I)=0
110   CONTINUE
C
C          Read in table, checking for valid version.
C
      IF(NODCAY.OR.ITDKY.EQ.0) RETURN
      REWIND ITDKY
C
      VOLD=VISAJE()
      READ(ITDKY,*) V
      IF(V.NE.VOLD) THEN
        WRITE(ITLIS,2000) V,VOLD
2000    FORMAT(//
     $  '  ***WARNING: DECAY TABLE DOES NOT MATCH ISAJET VERSION'/
     $  '  ***DECAY VERSION  : ',A40/
     $  '  ***PROGRAM VERSION: ',A40)
      ENDIF
C
200   LOOP=LOOP+1
      IF(LOOP.GT.MXDKY) GO TO 9999
220   DO 210 I=1,5
        IMODE(I)=0
        LMODE(I)=IBLANK
210   CONTINUE
      READ(ITDKY,*) IRES,ITYPE,BR,IMODE
C
      IF(IRES.NE.0) THEN
        IF(NOPI0.AND.IRES.EQ.110) GO TO 220
        IF(NOETA.AND.IRES.EQ.220) GO TO 220
        IF(NOB.AND.IRES.GT.100.AND.IRES.LT.1000.AND.
     $  MOD(IRES/10,10).EQ.5) GO TO 220
        IF(NOTAU.AND.IRES.EQ.16) GO TO 220
        IF(IRES.NE.IOLD) THEN
          CALL FLAVOR(IRES,IFL1,IFL2,IFL3,JSPIN,INDEX)
          LOOK(INDEX)=LOOP
        ENDIF
        IOLD=IRES
        CBR(LOOP)=BR
        MELEM(LOOP)=ITYPE
        DO 240 I=1,5
          MODE(I,LOOP)=IMODE(I)
          IF(IMODE(I).NE.0) LMODE(I)=LABEL(IMODE(I))
240     CONTINUE
        LRES=LABEL(IRES)
        IF(LPRINT) WRITE(ITLIS,20) LRES,(LMODE(K),K=1,5),
     1  BR,IRES,(IMODE(K),K=1,5)
20      FORMAT(3X,A5,4X,5(A5,2X),F8.5,10X,I5,4X,5(I7,2X))
        GO TO 200
      ENDIF
C
C          Add TAU -> TAUL,TAUR if NOTAU
C
      IF(NOTAU) THEN
        IRES=16
        LRES=LABEL(IRES)
        LOOP=LOOP+1
        CALL FLAVOR(16,IFL1,IFL2,IFL3,JSPIN,INDEX)
        LOOK(INDEX)=LOOP
        BR=.5
        CBR(LOOP)=BR
        MELEM(LOOP)=8
        MODE(1,LOOP)=IDTAUL
        LMODE(1)=LABEL(IDTAUL)
        DO 241 I=2,5
          MODE(I,LOOP)=0
          LMODE(I)=LABEL(0)
241     CONTINUE
        IF(LPRINT) WRITE(ITLIS,20) LRES,(LMODE(K),K=1,5),
     1  BR,IRES,(IMODE(K),K=1,5)
        LOOP=LOOP+1
        BR=1.0
        CBR(LOOP)=BR
        MELEM(LOOP)=8
        MODE(1,LOOP)=IDTAUR
        LMODE(1)=LABEL(IDTAUR)
        DO 242 I=2,5
          MODE(I,LOOP)=0
242     CONTINUE
        IF(LPRINT) WRITE(ITLIS,20) LRES,(LMODE(K),K=1,5),
     1  BR,IRES,(IMODE(K),K=1,5)
      ENDIF
C
C          Add HIGGS for E+E- or WHIGGS
C
      IF((KEYS(2).OR.KEYS(10)).AND..NOT.GOMSSM) THEN
          SUMGAM=0
          SUMBR=0
          DO 244 J=1,NSSMOD
            IF(ISSMOD(J).EQ.81.AND.GSSMOD(J).GT.0) THEN
              SUMGAM=SUMGAM+GSSMOD(J)
            ENDIF
244       CONTINUE
          DO 245 J=1,NSSMOD
            IF(ISSMOD(J).EQ.81.AND.GSSMOD(J).GT.0) THEN
              BSSMOD(J)=GSSMOD(J)/SUMGAM
            ENDIF
245       CONTINUE
          DO 246 J=1,NSSMOD
            IF(ISSMOD(J).EQ.81.AND.BSSMOD(J).GT.0) THEN
              SUMBR=SUMBR+BSSMOD(J)
            ENDIF
246       CONTINUE
C          If modes exist, add them
          IF(SUMBR.LE.0) GO TO 249
          IRES=81
          LRES=LABEL(IRES)
          CALL FLAVOR(IRES,IFL1,IFL2,IFL3,JSPIN,INDEX)
          LOOK(INDEX)=LOOP+1
          SUMBR2=0
          DO 247 J=1,NSSMOD
            IF(ISSMOD(J).EQ.81.AND.BSSMOD(J).GT.0) THEN
              LOOP=LOOP+1
              SUMBR2=SUMBR2+BSSMOD(J)
              BR=SUMBR2/SUMBR
              CBR(LOOP)=BR
              MELEM(LOOP)=MSSMOD(J)
              DO 248 K=1,5
                MODE(K,LOOP)=JSSMOD(K,J)
                LMODE(K)=LABEL(MODE(K,LOOP))
248           CONTINUE
              IF(LPRINT) WRITE(ITLIS,20) LRES,(LMODE(K),K=1,5),
     $        BR,IRES,(MODE(K,LOOP),K=1,5)
            ENDIF
247       CONTINUE
249     CONTINUE
      END IF
C
C          Add MSSM decay modes if applicable, OR H_SM FOR WHIGGS
C
      IF(GOMSSM) THEN
        DO 250 I=1,NOUT
C          Check for modes
          SUMBR=0
          DO 251 J=1,NSSMOD
            IF(ISSMOD(J).EQ.IDOUT(I).AND.BSSMOD(J).GT.0) THEN
              SUMBR=SUMBR+BSSMOD(J)
            ENDIF
251       CONTINUE
C          If modes exist, add them
          IF(SUMBR.LE.0) GO TO 250
          IRES=IDOUT(I)
          LRES=LABEL(IRES)
          CALL FLAVOR(IRES,IFL1,IFL2,IFL3,JSPIN,INDEX)
          LOOK(INDEX)=LOOP+1
          SUMBR2=0
          DO 252 J=1,NSSMOD
            IF(ISSMOD(J).EQ.IDOUT(I).AND.BSSMOD(J).GT.0) THEN
              LOOP=LOOP+1
              SUMBR2=SUMBR2+BSSMOD(J)
              BR=SUMBR2/SUMBR
              CBR(LOOP)=BR
              MELEM(LOOP)=MSSMOD(J)
              DO 253 K=1,5
                MODE(K,LOOP)=JSSMOD(K,J)
                LMODE(K)=LABEL(MODE(K,LOOP))
253           CONTINUE
              IF(LPRINT) WRITE(ITLIS,20) LRES,(LMODE(K),K=1,5),
     $        BR,IRES,(MODE(K,LOOP),K=1,5)
            ENDIF
252       CONTINUE
250     CONTINUE
C
C          Top hadron decays
C
        DO 260 I=1,NTHAD
C          Check for modes
          SUMBR=0
          DO 261 J=1,NSSMOD
            IF(ISSMOD(J).EQ.6.AND.BSSMOD(J).GT.0) THEN
              SUMBR=SUMBR+BSSMOD(J)
            ENDIF
261       CONTINUE
C          If modes exist, add them -- conjugate for antimesons
          IF(SUMBR.LE.0) GO TO 260
          IRES=IABS(ITHAD(I))
          LRES=LABEL(IRES)
          CALL FLAVOR(IRES,IFL1,IFL2,IFL3,JSPIN,INDEX)
          LOOK(INDEX)=LOOP+1
          SUMBR2=0
          DO 262 J=1,NSSMOD
            IF(ISSMOD(J).EQ.6.AND.BSSMOD(J).GT.0) THEN
              LOOP=LOOP+1
              SUMBR2=SUMBR2+BSSMOD(J)
              BR=SUMBR2/SUMBR
              CBR(LOOP)=BR
              IF(IABS(JSSMOD(1,J)).LT.20.AND.IABS(JSSMOD(2,J)).LT.20
     $        .AND.IABS(JSSMOD(3,J)).LT.20.AND.IABS(JSSMOD(4,J)).LT.20
     $        .AND.IABS(JSSMOD(5,J)).LT.20) THEN
                MELEM(LOOP)=4
              ELSE
                MELEM(LOOP)=0
              ENDIF
              DO 263 K=1,5
                IF(ITHAD(I).GT.0) THEN
                  MODE(K,LOOP)=JSSMOD(K,J)
                ELSE
                  MODE(K,LOOP)=IDANTI(JSSMOD(K,J))
                ENDIF
                LMODE(K)=LABEL(MODE(K,LOOP))
263           CONTINUE
              IF(LPRINT) WRITE(ITLIS,20) LRES,(LMODE(K),K=1,5),
     $        BR,IRES,(MODE(K,LOOP),K=1,5)
            ENDIF
262       CONTINUE
260     CONTINUE
      ENDIF
C
C          Set forced decay modes.
C          LOOK(INDEX) = -IRES, where LOOK2(K,IRES) points to entries in
C          decay table for IDENT>0 and IDENT<0.
C          LOOKST(IRES) = standard LOOK value.
C
      IF(NFORCE.EQ.0) GO TO 400
C          Append each forced decay to table
      IRES=0
      DO 310 I=1,NFORCE
        IF(IFORCE(I).EQ.0) GO TO 310
        LOOP=LOOP+1
        IF(LOOP.GT.MXDKY) GO TO 9999
        CALL FLAVOR(IFORCE(I),IFL1,IFL2,IFL3,JSPIN,INDEX)
        IF(IFORCE(I).GT.0) THEN
          KTYPE=1
        ELSE
          KTYPE=2
        ENDIF
C
        IF(LOOK(INDEX).GE.0) THEN
          IRES=IRES+1
          IF(IRES.GT.MXFORC) GO TO 9998
          LOOKST(IRES)=LOOK(INDEX)
          LOOK2(KTYPE,IRES)=LOOP
          LOOK2(3-KTYPE,IRES)=LOOKST(IRES)
          LOOK(INDEX)=-IRES
        ELSE
          IRES2=-LOOK(INDEX)
          IF(IRES2.GT.MXFORC) GO TO 9998
          LOOK2(KTYPE,IRES2)=LOOP
        ENDIF
C          Set forced decay mode - conjugate if necessary
        IF(KTYPE.EQ.1) THEN
          DO 320 K=1,5
320       MODE(K,LOOP)=MFORCE(K,I)
        ELSE
          DO 330 K=1,5
330       MODE(K,LOOP)=IDANTI(MFORCE(K,I))
        ENDIF
        CBR(LOOP)=1.
C          Set matrix element flag
        CALL ORDER(IFORCE(I),MFORCE(1,I),IDUMMY,MEOUT,.FALSE.)
        MELEM(LOOP)=MEOUT
        MEFORC(I)=MEOUT
310   CONTINUE
C
400   RETURN
C
C          Errors
C
9999  WRITE(ITLIS,3001) LOOP
3001  FORMAT(//' ***** ERROR IN SETDKY ... DECAY COUNTER LOOP = ',
     $I6,' *****')
      STOP 99
9998  WRITE(ITLIS,3002) IRES
3002  FORMAT(//' ***** ERROR IN SETDKY ... FORCE COUNTER IRES = ',
     $I6,' *****')
      STOP 99
      END
