CDECK  ID>, HDCYAD
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :  Peter Richardson
C-----------------------------------------------------------------------
      SUBROUTINE HDCYAD(HIGGS,WIDTH,BRFRAC,DECAY)
C-----------------------------------------------------------------------
C--Subroutine to add Higgs modes to the ISAJET decay tables
C-----------------------------------------------------------------------
      IMPLICIT NONE
C--Common block containing ISAJET decay modes
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
C
      DOUBLE PRECISION WIDTH,BRFRAC
      INTEGER HIGGS,DECAY(3),J
      NSSMOD = NSSMOD+1
      IF(NSSMOD.GT.MXSS) THEN
        WRITE(*,*) 'TOO MANY MODES'
        WRITE(*,*) 'RERUN WITH INCREASED MXSS'
        WRITE(*,*) 'STOPPING'
        STOP
      ENDIF
      ISSMOD(NSSMOD) = HIGGS
      DO J=1,3
        JSSMOD(J,NSSMOD) = DECAY(J)
      ENDDO
      DO J=4,5
        JSSMOD(J,NSSMOD) = 0
      ENDDO
      GSSMOD(NSSMOD) = WIDTH*BRFRAC
      BSSMOD(NSSMOD) = BRFRAC
      MSSMOD(NSSMOD) = 0
      END
