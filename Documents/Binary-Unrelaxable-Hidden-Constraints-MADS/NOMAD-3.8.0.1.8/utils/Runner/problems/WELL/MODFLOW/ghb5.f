      SUBROUTINE GHB5AL(ISUM,LENX,LCBNDS,NBOUND,MXBND,IN,IOUT,IGHBCB,
     1        NGHBVL,IGHBAL,IFREFM)
C
C-----VERSION 0943 21FEB1996 GHB5AL
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR HEAD-DEPENDENT BOUNDARIES
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      COMMON /GHBCOM/GHBAUX(5)
      CHARACTER*16 GHBAUX
      CHARACTER*80 LINE
C     ------------------------------------------------------------------
C
C1------IDENTIFY PACKAGE AND INITIALIZE # OF GENERAL HEAD BOUNDS.
      WRITE(IOUT,1)IN
1     FORMAT(1X,/1X,'GHB5 -- GHB PACKAGE, VERSION 5, 9/1/93',
     1' INPUT READ FROM UNIT',I3)
      NBOUND=0
C
C2------READ MAXIMUM NUMBER OF BOUNDS AND UNIT OR FLAG FOR
C2------CELL-BY-CELL FLOW TERMS.
      READ(IN,'(A)') LINE
      IF(IFREFM.EQ.0) THEN
         READ(LINE,'(2I10)') MXBND,IGHBCB
         LLOC=21
      ELSE
         LLOC=1
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXBND,R,IOUT,IN)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IGHBCB,R,IOUT,IN)
      END IF
      WRITE(IOUT,3) MXBND
    3 FORMAT(1X,'MAXIMUM OF',I5,' HEAD-DEPENDENT BOUNDARY NODES')
      IF(IGHBCB.LT.0) WRITE(IOUT,7)
    7 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE PRINTED WHEN ICBCFL NOT 0')
      IF(IGHBCB.GT.0) WRITE(IOUT,8) IGHBCB
    8 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT',I3)
C
C3------READ AUXILIARY PARAMETERS AND CBC ALLOCATION OPTION.
      IGHBAL=0
      NAUX=0
   10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'CBCALLOCATE' .OR.
     1   LINE(ISTART:ISTOP).EQ.'CBC') THEN
         IGHBAL=1
         WRITE(IOUT,11)
   11    FORMAT(1X,'MEMORY IS ALLOCATED FOR CELL-BY-CELL BUDGET TERMS')
         GO TO 10
      ELSE IF(LINE(ISTART:ISTOP).EQ.'AUXILIARY' .OR.
     1        LINE(ISTART:ISTOP).EQ.'AUX') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
         IF(NAUX.LT.5) THEN
            NAUX=NAUX+1
            GHBAUX(NAUX)=LINE(ISTART:ISTOP)
            WRITE(IOUT,12) GHBAUX(NAUX)
   12       FORMAT(1X,'AUXILIARY BOUNDARY PARAMETER: ',A)
         END IF
         GO TO 10
      END IF
      NGHBVL=5+NAUX+IGHBAL
C
C4------ALLOCATE SPACE IN THE X ARRAY FOR THE BNDS ARRAY.
      LCBNDS=ISUM
      ISP=NGHBVL*MXBND
      ISUM=ISUM+ISP
C
C5------PRINT AMOUNT OF SPACE USED BY THE GHB PACKAGE.
      WRITE(IOUT,14) ISP
   14 FORMAT(1X,I10,' ELEMENTS IN X ARRAY ARE USED BY GHB')
      ISUM1=ISUM-1
      WRITE(IOUT,15) ISUM1,LENX
   15 FORMAT(1X,I10,' ELEMENTS OF X ARRAY USED OUT OF ',I10)
      IF(ISUM1.GT.LENX) WRITE(IOUT,16)
   16 FORMAT(1X,'   ***X ARRAY MUST BE DIMENSIONED LARGER***')
C
C6------RETURN.
      RETURN
      END
      SUBROUTINE GHB5RP(BNDS,NBOUND,MXBND,IN,IOUT,NGHBVL,IGHBAL,IFREFM)
C
C-----VERSION 0946 21FEB1996 GHB5RP
C     ******************************************************************
C     READ DATA FOR GHB
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      DIMENSION BNDS(NGHBVL,MXBND)
      COMMON /GHBCOM/GHBAUX(5)
      CHARACTER*16 GHBAUX
      CHARACTER*151 LINE
C     ------------------------------------------------------------------
C
C1------READ ITMP (# OF GENERAL HEAD BOUNDS OR FLAG TO REUSE DATA).
      IF(IFREFM.EQ.0) THEN
         READ(IN,'(I10)') ITMP
      ELSE
         READ(IN,*) ITMP
      END IF
C
C2------TEST ITMP
      IF(ITMP.GE.0) GO TO 50
C
C2A-----IF ITMP<0 THEN REUSE DATA FROM LAST STRESS PERIOD.
      WRITE(IOUT,7)
    7 FORMAT(1X,/1X,'REUSING HEAD-DEPENDENT BOUNDS FROM LAST STRESS',
     1      ' PERIOD')
      GO TO 260
C
C3------IF ITMP=>0 THEN IT IS THE # OF GENERAL HEAD BOUNDS.
   50 NBOUND=ITMP
C
C4------IF MAX NUMBER OF BOUNDS IS EXCEEDED THEN STOP.
      IF(NBOUND.LE.MXBND) GO TO 100
      WRITE(IOUT,99) NBOUND,MXBND
   99 FORMAT(1X,/1X,'NBOUND(',I4,') IS GREATER THAN MXBND(',I4,')')
C
C4A-----ABNORMAL STOP.
      STOP
C
C5------PRINT # OF GENERAL HEAD BOUNDS THIS STRESS PERIOD.
  100 WRITE(IOUT,101) NBOUND
  101 FORMAT(1X,//1X,I5,' HEAD-DEPENDENT BOUNDARY NODES')
C
C6------IF THERE ARE NO GENERAL HEAD BOUNDS THEN RETURN.
      IF(NBOUND.EQ.0) GO TO 260
C
C7------READ & PRINT DATA FOR EACH GENERAL HEAD BOUNDARY.
      NAUX=NGHBVL-5-IGHBAL
      MAXAUX=NGHBVL-IGHBAL
      IF(NAUX.GT.0) THEN
         WRITE(IOUT,103) (GHBAUX(JJ),JJ=1,NAUX)
         WRITE(IOUT,104) ('------------------',JJ=1,NAUX)
      ELSE
         WRITE(IOUT,103)
         WRITE(IOUT,104)
      END IF
  103 FORMAT(1X,/1X,'LAYER   ROW   COL   ELEVATION   CONDUCTANCE   ',
     1           'BOUND NO.',:5(2X,A))
  104 FORMAT(1X,55('-'),5A)
      DO 250 II=1,NBOUND
C7A-----READ THE REQUIRED DATA WITH FIXED OR FREE FORMAT.
      READ(IN,'(A)') LINE
      IF(IFREFM.EQ.0) THEN
         READ(LINE,'(3I10,2F10.0)') K,I,J,(BNDS(JJ,II),JJ=4,5)
         LLOC=51
      ELSE
         LLOC=1
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,K,R,IOUT,IN)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I,R,IOUT,IN)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,J,R,IOUT,IN)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,N,BNDS(4,II),IOUT,IN)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,N,BNDS(5,II),IOUT,IN)
      END IF
C7B-----READ ANY AUXILIARY DATA WITH FREE FORMAT, AND PRINT ALL VALUES.
      IF(NAUX.GT.0) THEN
         DO 110 JJ=1,NAUX
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,N,BNDS(JJ+5,II),IOUT,IN)
  110    CONTINUE
         WRITE (IOUT,115) K,I,J,BNDS(4,II),BNDS(5,II),II,
     1         (BNDS(JJ,II),JJ=6,MAXAUX)
      ELSE
         WRITE (IOUT,115) K,I,J,BNDS(4,II),BNDS(5,II),II
      END IF
  115 FORMAT(1X,I4,I7,I6,G13.4,G14.4,I8,:5(2X,G16.5))
      BNDS(1,II)=K
      BNDS(2,II)=I
      BNDS(3,II)=J
  250 CONTINUE
C
C8------RETURN.
  260 RETURN
      END
      SUBROUTINE GHB5FM(NBOUND,MXBND,BNDS,HCOF,RHS,IBOUND,
     1              NCOL,NROW,NLAY,NGHBVL)
C
C-----VERSION 1352 28AUG1992 GHB5FM
C     ******************************************************************
C     ADD GHB TERMS TO RHS AND HCOF
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      DIMENSION BNDS(NGHBVL,MXBND),HCOF(NCOL,NROW,NLAY),
     1         RHS(NCOL,NROW,NLAY),IBOUND(NCOL,NROW,NLAY)
C     ------------------------------------------------------------------
C
C1------IF NBOUND<=0 THEN THERE ARE NO GENERAL HEAD BOUNDS. RETURN.
      IF(NBOUND.LE.0) RETURN
C
C2------PROCESS EACH ENTRY IN THE GENERAL HEAD BOUND LIST (BNDS).
      DO 100 L=1,NBOUND
C
C3------GET COLUMN, ROW AND LAYER OF CELL CONTAINING BOUNDARY.
      IL=BNDS(1,L)
      IR=BNDS(2,L)
      IC=BNDS(3,L)
C
C4------IF THE CELL IS EXTERNAL THEN SKIP IT.
      IF(IBOUND(IC,IR,IL).LE.0) GO TO 100
C
C5------SINCE THE CELL IS INTERNAL GET THE BOUNDARY DATA.
      HB=BNDS(4,L)
      C=BNDS(5,L)
C
C6------ADD TERMS TO RHS AND HCOF.
      HCOF(IC,IR,IL)=HCOF(IC,IR,IL)-C
      RHS(IC,IR,IL)=RHS(IC,IR,IL)-C*HB
  100 CONTINUE
C
C7------RETURN.
      RETURN
      END
      SUBROUTINE GHB5BD(NBOUND,MXBND,VBNM,VBVL,MSUM,BNDS,DELT,HNEW,
     1   NCOL,NROW,NLAY,IBOUND,KSTP,KPER,IGHBCB,ICBCFL,BUFF,IOUT,
     2   PERTIM,TOTIM,NGHBVL,IGHBAL)
C-----VERSION 1410 07APRIL1993 GHB5BD
C     ******************************************************************
C     CALCULATE VOLUMETRIC BUDGET FOR GHB
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 VBNM(MSUM),TEXT
      DOUBLE PRECISION HNEW,CC,CHB,RATIN,RATOUT,RRATE
      DIMENSION VBVL(4,MSUM),BNDS(NGHBVL,MXBND),HNEW(NCOL,NROW,NLAY),
     1           IBOUND(NCOL,NROW,NLAY),BUFF(NCOL,NROW,NLAY)
C
      DATA TEXT /' HEAD DEP BOUNDS'/
C     ------------------------------------------------------------------
C
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD) AND
C1------ACCUMULATORS (RATIN AND RATOUT).
      ZERO=0.
      RATOUT=ZERO
      RATIN=ZERO
      IBD=0
      IF(IGHBCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IGHBCB.GT.0) IBD=ICBCFL
      IBDLBL=0
C
C2------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IF(IBD.EQ.2) CALL UBDSV2(KSTP,KPER,TEXT,IGHBCB,NCOL,NROW,NLAY,
     1          NBOUND,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C
C3------CLEAR THE BUFFER.
      DO 50 IL=1,NLAY
      DO 50 IR=1,NROW
      DO 50 IC=1,NCOL
      BUFF(IC,IR,IL)=ZERO
50    CONTINUE
C
C4------IF NO BOUNDARIES, SKIP FLOW CALCULATIONS.
      IF(NBOUND.EQ.0) GO TO 200
C
C5------LOOP THROUGH EACH BOUNDARY CALCULATING FLOW.
      DO 100 L=1,NBOUND
C
C5A-----GET LAYER, ROW AND COLUMN OF EACH GENERAL HEAD BOUNDARY.
      IL=BNDS(1,L)
      IR=BNDS(2,L)
      IC=BNDS(3,L)
      RATE=ZERO
C
C5B-----IF CELL IS NO-FLOW OR CONSTANT-HEAD, THEN IGNORE IT.
      IF(IBOUND(IC,IR,IL).LE.0) GO TO 99
C
C5C-----GET PARAMETERS FROM BOUNDARY LIST.
      HB=BNDS(4,L)
      C=BNDS(5,L)
      CC=C
C
C5D-----CALCULATE THE FOW RATE INTO THE CELL.
      CHB=C*HB
      RRATE=CHB - CC*HNEW(IC,IR,IL)
      RATE=RRATE
C
C5E-----PRINT THE INDIVIDUAL RATES IF REQUESTED(IGHBCB<0).
      IF(IBD.LT.0) THEN
         IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
   61    FORMAT(1X,/1X,A,'   PERIOD',I3,'   STEP',I3)
         WRITE(IOUT,62) L,IL,IR,IC,RATE
   62    FORMAT(1X,'BOUNDARY',I4,'   LAYER',I3,'   ROW',I4,'   COL',I4,
     1       '   RATE',1PG15.6)
         IBDLBL=1
      END IF
C
C5F-----ADD RATE TO BUFFER.
      BUFF(IC,IR,IL)=BUFF(IC,IR,IL)+RATE
C
C5G-----SEE IF FLOW IS INTO AQUIFER OR OUT OF AQUIFER.
      IF(RATE)94,99,96
C
C5H------FLOW IS OUT OF AQUIFER SUBTRACT RATE FROM RATOUT.
94    RATOUT=RATOUT-RRATE
      GO TO 99
C
C5I-----FLOW IS INTO AQIFER; ADD RATE TO RATIN.
96    RATIN=RATIN+RRATE
C
C5J-----IF SAVING CELL-BY-CELL FLOWS IN LIST, WRITE FLOW.  OR IF
C5J-----RETURNING THE FLOW IN THE BNDS ARRAY, COPY FLOW TO BNDS.
99    IF(IBD.EQ.2) CALL UBDSVA(IGHBCB,NCOL,NROW,IC,IR,IL,RATE,IBOUND,
     1                        NLAY)
      IF(IGHBAL.NE.0) BNDS(NGHBVL,L)=RATE
100   CONTINUE
C
C6------IF CELL-BY-CELL TERMS WILL BE SAVED AS A 3-D ARRAY, THEN CALL
C6------UTILITY MODULE UBUDSV TO SAVE THEM.
      IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IGHBCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
C
C7------MOVE RATES, VOLUMES AND LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVL(3,MSUM)=RIN
      VBVL(1,MSUM)=VBVL(1,MSUM)+RIN*DELT
      VBVL(4,MSUM)=ROUT
      VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
      VBNM(MSUM)=TEXT
C
C8------INCREMENT THE BUDGET TERM COUNTER.
      MSUM=MSUM+1
C
C9------RETURN.
      RETURN
      END
