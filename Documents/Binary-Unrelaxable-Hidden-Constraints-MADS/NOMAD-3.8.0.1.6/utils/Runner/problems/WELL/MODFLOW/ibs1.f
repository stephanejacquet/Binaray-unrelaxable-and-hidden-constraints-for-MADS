C FILE IBS1.F77 -- June, 1996 -- 3 statements in the version documented
C   in TWRI 6-A2 have been modified in order to correct a problem.
C   Although subsidence is only meant to be active for layers in which
C   IBQ>0, some of the subroutines performed subsidence calculations when
C   IBQ<0.  Note that this was a problem only if negative IBQ values
C   were specified.  That is, the code has always worked correctly for
C   IBQ=0 and IBQ>0.
      SUBROUTINE IBS1AL(ISUM,LENX,LCHC,LCSCE,LCSCV,LCSUB,
     1                  NCOL,NROW,NLAY,IIBSCB,IIBSOC,ISS,IN,IOUT)
C
C-----VERSION 07JUN1996 IBS1AL
C-----VERSION 01AUG1996 -- modified to allow 200 layers instead of 80
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR INTERBED STORAGE PACKAGE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      DIMENSION IBQ1(200)
      COMMON /IBSCOM/ IBQ(200)
C     ------------------------------------------------------------------
C
C1------IDENTIFY PACKAGE.
      WRITE(IOUT,1)IN
    1 FORMAT(1H0,'IBS1 -- INTERBED STORAGE PACKAGE, VERSION 1,',
     1     ' 06/07/96',' INPUT READ FROM UNIT',I3)
C
C2------CHECK TO SEE THAT INTERBED STORAGE OPTION IS APPROPRIATE
      IF(ISS.EQ.0) GO TO 100
C
C3------IF INAPPROPRIATE PRINT A MESSAGE & CANCEL OPTION.
      WRITE(IOUT,8)
    8 FORMAT(1X,'INTERBED STORAGE INAPPROPRIATE FOR STEADY-STATE',
     1 ' PROBLEM.',/,1X,'OPTION CANCELLED, SIMULATION CONTINUING.')
      IN=0
      RETURN
C
C4------READ FLAG FOR STORING CELL-BY-CELL STORAGE CHANGES AND
C4------FLAG FOR PRINTING AND STORING COMPACTION, SUBSIDENCE, AND
C4------CRITICAL HEAD ARRAYS.
  100 READ(IN,3) IIBSCB,IIBSOC
    3 FORMAT(2I10)
C
C5------IF CELL-BY-CELL TERMS TO BE SAVED THEN PRINT UNIT NUMBER.
      IF(IIBSCB.GT.0) WRITE(IOUT,105) IIBSCB
  105 FORMAT(1X,'CELL-BY-CELL FLOW TERMS WILL BE SAVED ON UNIT',I3)
C
C5A-----IF OUTPUT CONTROL FOR PRINTING ARRAYS IS SELECTED PRINT MESSAGE.
      IF(IIBSOC.GT.0) WRITE(IOUT,106)
  106 FORMAT(1X,'OUTPUT CONTROL RECORDS FOR IBS1 PACKAGE WILL BE ',
     1 'READ EACH TIME STEP.')
C
C6------READ INDICATOR AND FIND OUT HOW MANY LAYERS HAVE INTERBED STORAGE.
      READ(IN,110) (IBQ(K),K=1,NLAY)
  110 FORMAT(40I2)
      NAQL=0
      DO 120 K=1,NLAY
      IF(IBQ(K).LE.0) GO TO 120
      NAQL=NAQL+1
      IBQ1(NAQL)=K
  120 CONTINUE
C
C7------IDENTIFY WHICH LAYERS HAVE INTERBED STORAGE.
      WRITE(IOUT,130) (IBQ1(K),K=1,NAQL)
  130 FORMAT(1X,'INTERBED STORAGE IN LAYER(S) ',80I2)
C
C8------ALLOCATE SPACE FOR THE ARRAYS HC, SCE, SCV, AND SUB.
      IRK=ISUM
      NA=NROW*NCOL*NAQL
      LCHC=ISUM
      ISUM=ISUM+NA
      LCSCE=ISUM
      ISUM=ISUM+NA
      LCSCV=ISUM
      ISUM=ISUM+NA
      LCSUB=ISUM
      ISUM=ISUM+NA
C
C9------CALCULATE & PRINT AMOUNT OF SPACE USED BY PACKAGE.
  300 IRK=ISUM-IRK
      WRITE(IOUT,4)IRK
    4 FORMAT(1X,I8,' ELEMENTS OF X ARRAY USED FOR INTERBED STORAGE')
      ISUM1=ISUM-1
      WRITE(IOUT,5)ISUM1,LENX
    5 FORMAT(1X,I8,' ELEMENTS OF X ARRAY USED OUT OF',I8)
      IF(ISUM1.GT.LENX)WRITE(IOUT,6)
    6 FORMAT(1X,'   ***X ARRAY MUST BE MADE LARGER***')
C
C10-----RETURN.
      RETURN
      END
      SUBROUTINE IBS1RP(DELR,DELC,HNEW,HC,SCE,SCV,SUB,NCOL,NROW,
     1                  NLAY,NODES,IIBSOC,ISUBFM,ICOMFM,IHCFM,
     2                  ISUBUN,ICOMUN,IHCUN,IN,IOUT)
C
C-----VERSION 1117 02JUN1988 IBS1RP
C-----VERSION 01AUG1996 -- modified to allow 200 layers instead of 80
C     ******************************************************************
C     READ INTERBED STORAGE DATA
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*4 ANAME
      DOUBLE PRECISION HNEW
      DIMENSION HNEW(NODES),HC(NODES),SCE(NODES),
     1       SCV(NODES),SUB(NODES),ANAME(6,4),
     2       DELR(NCOL),DELC(NROW)
C
      COMMON /IBSCOM/ IBQ(200)
C
      DATA ANAME(1,1),ANAME(2,1),ANAME(3,1),ANAME(4,1),ANAME(5,1),
     1  ANAME(6,1) /'   P','RECO','NSOL','IDAT','ION ','HEAD'/
      DATA ANAME(1,2),ANAME(2,2),ANAME(3,2),ANAME(4,2),ANAME(5,2),
     1  ANAME(6,2) /'ELAS','TIC ','INTE','RBED',' STO','RAGE'/
      DATA ANAME(1,3),ANAME(2,3),ANAME(3,3),ANAME(4,3),ANAME(5,3),
     1  ANAME(6,3) /' VIR','GIN ','INTE','RBED',' STO','RAGE'/
      DATA ANAME(1,4),ANAME(2,4),ANAME(3,4),ANAME(4,4),ANAME(5,4),
     1  ANAME(6,4) /'    ',' STA','RTIN','G CO','MPAC','TION'/
C     ------------------------------------------------------------------
C
C1------READ IN STORAGE AND CRITICAL HEAD ARRAYS
      NIJ=NROW*NCOL
      KQ=0
      DO 60 K=1,NLAY
      IF(IBQ(K).LE.0) GO TO 60
      KQ=KQ+1
      LOC=1+(KQ-1)*NIJ
      CALL U2DREL(HC(LOC),ANAME(1,1),NROW,NCOL,K,IN,IOUT)
      CALL U2DREL(SCE(LOC),ANAME(1,2),NROW,NCOL,K,IN,IOUT)
      CALL U2DREL(SCV(LOC),ANAME(1,3),NROW,NCOL,K,IN,IOUT)
      CALL U2DREL(SUB(LOC),ANAME(1,4),NROW,NCOL,K,IN,IOUT)
   60 CONTINUE
C
C2------LOOP THROUGH ALL CELLS WITH INTERBED STORAGE.
      KQ=0
      DO 80 K=1,NLAY
      IF(IBQ(K).LE.0) GO TO 80
      KQ=KQ+1
      NQ=(KQ-1)*NIJ
      NK=(K-1)*NIJ
      DO 70 IR=1,NROW
      NQR=NQ+(IR-1)*NCOL
      NKR=NK+(IR-1)*NCOL
      DO 70 IC=1,NCOL
      LOC=NQR+IC
      LOCH=NKR+IC
C
C3------MULTIPLY STORAGE BY AREA TO GET STORAGE CAPACITY.
      AREA=DELR(IC)*DELC(IR)
      SCE(LOC)=SCE(LOC)*AREA
      SCV(LOC)=SCV(LOC)*AREA
C
C4------MAKE SURE THAT PRECONSOLIDATION HEAD VALUES
C4------ARE CONSISTANT WITH STARTING HEADS.
      IF(HC(LOC).GT.HNEW(LOCH)) HC(LOC)=HNEW(LOCH)
   70 CONTINUE
   80 CONTINUE
C
C5------INITIALIZE AND READ OUTPUT FLAGS.
      ICOMFM=0
      ISUBFM=0
      IHCFM=0
      ICOMUN=0
      ISUBUN=0
      IHCUN=0
      IF(IIBSOC.LE.0) GO TO 200
      READ(IN,100) ISUBFM,ICOMFM,IHCFM,ISUBUN,ICOMUN,IHCUN
  100 FORMAT(6I10)
      WRITE(IOUT,110) ISUBFM,ICOMFM,IHCFM
  110 FORMAT(1H0,'    SUBSIDENCE PRINT FORMAT IS NUMBER',I4/
     1          '     COMPACTION PRINT FORMAT IS NUMBER',I4/
     2          '  CRITICAL HEAD PRINT FORMAT IS NUMBER',I4)
      IF(ISUBUN.GT.0) WRITE(IOUT,120) ISUBUN
  120 FORMAT(1H0,'    UNIT FOR SAVING SUBSIDENCE IS',I4)
      IF(ICOMUN.GT.0) WRITE(IOUT,130) ICOMUN
  130 FORMAT(1H ,'    UNIT FOR SAVING COMPACTION IS',I4)
      IF(IHCUN.GT.0)  WRITE(IOUT,140) IHCUN
  140 FORMAT(1H ,' UNIT FOR SAVING CRITICAL HEAD IS',I4)
C
C6------RETURN
  200 RETURN
      END
      SUBROUTINE IBS1FM(RHS,HCOF,HNEW,HOLD,HC,SCE,SCV,
     1                  IBOUND,NCOL,NROW,NLAY,DELT)
C
C-----VERSION 07JUN1996 IBS1FM
C-----VERSION 01AUG1996 -- modified to allow 200 layers instead of 80
C     ******************************************************************
C        ADD INTERBED STORAGE TO RHS AND HCOF
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      DOUBLE PRECISION HNEW
      DIMENSION RHS(NCOL,NROW,NLAY),HCOF(NCOL,NROW,NLAY),
     1          IBOUND(NCOL,NROW,NLAY),HNEW(NCOL,NROW,NLAY),
     2          HOLD(NCOL,NROW,NLAY),HC(NCOL,NROW,NLAY),
     3          SCE(NCOL,NROW,NLAY),SCV(NCOL,NROW,NLAY)
C
      COMMON /IBSCOM/ IBQ(200)
C     ------------------------------------------------------------------
C
C1------INITIALIZE
       TLED=1./DELT
      KQ=0
C
C2------FIND LAYERS WITH INTERBED STORAGE
      DO 110 K=1,NLAY
      IF(IBQ(K).LE.0) GO TO 110
      KQ=KQ+1
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
      IF(IBOUND(J,I,K).LE.0) GO TO 100
C
C3------DETERMINE STORAGE CAPACITIES FOR CELL AT START AND END OF STEP
      RHO1=SCE(J,I,KQ)*TLED
      RHO2=RHO1
      HCTMP=HC(J,I,KQ)
      IF(HNEW(J,I,K).LT.HCTMP) RHO2=SCV(J,I,KQ)*TLED
C
C4------ADD APPROPRIATE TERMS TO RHS AND HCOF
      RHS(J,I,K)=RHS(J,I,K)-HCTMP*(RHO2-RHO1)-RHO1*HOLD(J,I,K)
      HCOF(J,I,K)=HCOF(J,I,K)-RHO2
  100 CONTINUE
  110 CONTINUE
C
C5------RETURN
      RETURN
      END
      SUBROUTINE IBS1BD(IBOUND,HNEW,HOLD,HC,SCE,SCV,SUB,DELR,DELC,
     1      NCOL,NROW,NLAY,DELT,VBVL,VBNM,MSUM,KSTP,KPER,IIBSCB,
     2      ICBCFL,BUFF,IOUT)
C-----VERSION 07JUN1996 IBS1BD
C-----VERSION 01AUG1996 -- modified to allow 200 layers instead of 80
C     ******************************************************************
C     CALCULATE VOLUMETRIC BUDGET FOR INTERBED STORAGE
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*4 TEXT,VBNM
      DOUBLE PRECISION HNEW
      DIMENSION IBOUND(NCOL,NROW,NLAY),HOLD(NCOL,NROW,NLAY),
     1          HNEW(NCOL,NROW,NLAY),HC(NCOL,NROW,NLAY),
     2          SCE(NCOL,NROW,NLAY),SCV(NCOL,NROW,NLAY),
     3          SUB(NCOL,NROW,NLAY),VBVL(4,20),VBNM(4,20),
     4          BUFF(NCOL,NROW,NLAY),DELR(NCOL),DELC(NROW)
      DIMENSION TEXT(4)
C
      COMMON /IBSCOM/ IBQ(200)
      DATA TEXT(1),TEXT(2),TEXT(3),TEXT(4) /'INTE','RBED',' STO','RAGE'/
C     ------------------------------------------------------------------
C
C1------INITIALIZE CELL-BY-CELL FLOW TERM FLAG (IBD) AND
C1------ACCUMULATORS (STOIN AND STOUT).
      IBD=0
      STOIN=0.
      STOUT=0.
C
C2------TEST TO SEE IF CELL-BY-CELL FLOW TERMS ARE NEEDED.
      IF(ICBCFL.EQ.0  .OR. IIBSCB.LE.0 ) GO TO 10
C
C3------CELL-BY-CELL FLOW TERMS ARE NEEDED SET IBD AND CLEAR BUFFER.
      IBD=1
      DO 5 IL=1,NLAY
      DO 5 IR=1,NROW
      DO 5 IC=1,NCOL
      BUFF(IC,IR,IL)=0.
    5 CONTINUE
C
C4------RUN THROUGH EVERY CELL IN THE GRID WITH INTERBED STORAGE.
   10 KQ=0
      TLED=1./DELT
      DO 110 K=1,NLAY
      IF(IBQ(K).LE.0) GO TO 110
      KQ=KQ+1
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
C
C5------CALCULATE FLOW FROM STORAGE (VARIABLE HEAD CELLS ONLY)
      IF(IBOUND(J,I,K).LE.0) GO TO 100
      HHOLD=HOLD(J,I,K)
      HHNEW=HNEW(J,I,K)
      HHC=HC(J,I,KQ)
C
C6------GET STORAGE CAPACITIES AT BEGINNING AND END OF TIME STEP.
      SBGN=SCE(J,I,KQ)
      SEND=SBGN
      IF(HHNEW.LT.HHC) SEND=SCV(J,I,KQ)
C
C7------CALCULATE VOLUME CHANGE IN INTERBED STORAGE FOR TIME STEP.
      STRG=HHC*(SEND-SBGN)+SBGN*HHOLD-SEND*HHNEW
C
C8------ACCUMULATE SUBSIDENCE ASSOCIATED WITH CHANGE IN STORAGE
      SUB(J,I,KQ)=SUB(J,I,KQ)+STRG/(DELR(J)*DELC(I))
C
C9------IF C-B-C FLOW TERMS ARE TO BE SAVED THEN ADD RATE TO BUFFER.
      IF(IBD.EQ.1) BUFF(J,I,K)=BUFF(J,I,K)+STRG*TLED
C
C10-----SEE IF FLOW IS INTO OR OUT OF STORAGE.
      IF(STRG)94,100,96
   94 STOUT=STOUT-STRG
      GO TO 100
   96 STOIN=STOIN+STRG
  100 CONTINUE
  110 CONTINUE
C
C11-----IF C-B-C FLOW TERMS WILL BE SAVED CALL UBUDSV TO RECORD THEM.
      IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IIBSCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
C
C12-----MOVE RATES,VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 VBVL(3,MSUM)=STOIN*TLED
      VBVL(4,MSUM)=STOUT*TLED
      VBVL(1,MSUM)=VBVL(1,MSUM)+STOIN
      VBVL(2,MSUM)=VBVL(2,MSUM)+STOUT
      VBNM(1,MSUM)=TEXT(1)
      VBNM(2,MSUM)=TEXT(2)
      VBNM(3,MSUM)=TEXT(3)
      VBNM(4,MSUM)=TEXT(4)
C
C13-----INCREMENT BUDGET TERM COUNTER
      MSUM=MSUM+1
C
C14-----UPDATE PRECONSOLIDATION HEAD ARRAY
      KQ=0
      DO 310 K=1,NLAY
      IF(IBQ(K).LE.0) GO TO 310
      KQ=KQ+1
      DO 300 I=1,NROW
      DO 300 J=1,NCOL
      IF(IBOUND(J,I,K).LE.0) GO TO 300
      HHNEW=HNEW(J,I,K)
      IF(HHNEW.LT.HC(J,I,KQ)) HC(J,I,KQ)=HHNEW
  300 CONTINUE
  310 CONTINUE
C
C15-----RETURN
      RETURN
      END
      SUBROUTINE IBS1OT(NCOL,NROW,NLAY,PERTIM,TOTIM,KSTP,KPER,NSTP,
     1           BUFF,SUB,HC,IIBSOC,ISUBFM,ICOMFM,IHCFM,ISUBUN,
     2           ICOMUN,IHCUN,IN,IOUT)
C-----VERSION 07JUN1996 IBS1OT
C-----VERSION 01AUG1996 -- modified to allow 200 layers instead of 80
C     ******************************************************************
C     PRINT AND STORE SUBSIDENCE, COMPACTION AND CRITICAL HEAD.
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*4 TEXT
      DIMENSION HC(NCOL,NROW,NLAY),SUB(NCOL,NROW,NLAY),
     1          BUFF(NCOL,NROW,NLAY),TEXT(4,3)
      COMMON /IBSCOM/ IBQ(200)
      DATA TEXT(1,1),TEXT(2,1),TEXT(3,1),TEXT(4,1) /'    ','  SU',
     1     'BSID','ENCE'/,TEXT(1,2),TEXT(2,2),TEXT(3,2),TEXT(4,2)
     2     /'    ','  CO','MPAC','TION'/,TEXT(1,3),TEXT(2,3),
     3     TEXT(3,3),TEXT(4,3) /'   C','RITI','CAL ','HEAD'/
C     ------------------------------------------------------------------
C
C1------INITIALIZE FLAGS FOR PRINTING AND SAVING SUBSIDENCE, COMPACTION,
C1------AND CRITICAL HEAD
      ISUBPR=0
      ICOMPR=0
      IHCPR=0
      ISUBSV=0
      ICOMSV=0
      IHCSV=0
      IF(KSTP.EQ.NSTP) ISUBPR=1
C2------READ FLAGS FOR PRINTING AND SAVING.
      IF(IIBSOC.LE.0) GO TO 28
      READ(IN,10) ISUBPR,ICOMPR,IHCPR,ISUBSV,ICOMSV,IHCSV
   10 FORMAT(6I10)
      WRITE(IOUT,15) ISUBPR,ICOMPR,IHCPR,ISUBSV,ICOMSV,IHCSV
   15 FORMAT(1H0,'FLAGS FOR PRINTING AND STORING SUBSIDENCE, ',
     1 'COMPACTION, AND CRITICAL HEAD:'/
     2 '   ISUBPR    ICOMPR    IHCPR     ISUBSV    ICOMSV    IHCSV   '/
     3 ' ------------------------------------------------------------'/
     4 I6,5I10)
C
C3------PRINT AND STORE SUBSIDENCE, FIRST, CLEAR OUT BUFF.
   28 IF(ISUBPR.LE.0.AND.ISUBSV.LE.0) GO TO 100
      DO 30 IR=1,NROW
      DO 30 IC=1,NCOL
      BUFF(IC,IR,1)=0.
   30 CONTINUE
C
C4------SUM COMPACTION IN ALL LAYERS TO GET SUBSIDENCE.
      KQ=0
      DO 50 K=1,NLAY
      IF(IBQ(K).LE.0) GO TO 50
      KQ=KQ+1
      DO 40 I=1,NROW
      DO 40 J=1,NCOL
      BUFF(J,I,1)=BUFF(J,I,1)+SUB(J,I,KQ)
   40 CONTINUE
   50 CONTINUE
C
C5------PRINT SUBSIDENCE.
      IF(ISUBPR.LE.0) GO TO 60
      IF(ISUBFM.LT.0) CALL ULAPRS(BUFF,TEXT(1,1),KSTP,KPER,NCOL,NROW,1,
     1          -ISUBFM,IOUT)
      IF(ISUBFM.GE.0) CALL ULAPRW(BUFF,TEXT(1,1),KSTP,KPER,NCOL,NROW,1,
     1           ISUBFM,IOUT)
C
C6------STORE SUBSIDENCE.
   60 IF(ISUBSV.LE.0) GO TO 100
      CALL ULASAV(BUFF,TEXT(1,1),KSTP,KPER,PERTIM,TOTIM,NCOL,NROW,1,
     1             ISUBUN)
C
C7------PRINT COMPACTION FOR ALL LAYERS WITH INTERBED STORAGE.
  100 IF(ICOMPR.LE.0) GO TO 140
      KQ=0
      DO 130 K=1,NLAY
      IF(IBQ(K).LE.0) GO TO 130
      KQ=KQ+1
      IF(ICOMFM.LT.0) CALL ULAPRS(SUB(1,1,KQ),TEXT(1,2),KSTP,KPER,NCOL,
     1          NROW,K,-ICOMFM,IOUT)
      IF(ICOMFM.GE.0) CALL ULAPRW(SUB(1,1,KQ),TEXT(1,2),KSTP,KPER,NCOL,
     1           NROW,K,ICOMFM,IOUT)
  130 CONTINUE
C
C8------SAVE COMPACTION FOR ALL LAYERS WITH INTERBED STORAGE.
  140 IF(ICOMSV.LE.0) GO TO 200
      KQ=0
      DO 160 K=1,NLAY
      IF(IBQ(K).LE.0) GO TO 160
      KQ=KQ+1
      CALL ULASAV(SUB(1,1,KQ),TEXT(1,2),KSTP,KPER,PERTIM,TOTIM,NCOL,
     1            NROW,K,ICOMUN)
  160 CONTINUE
C
C9------PRINT CRITICAL HEAD FOR ALL LAYERS WITH INTERBED STORAGE.
  200 IF(IHCPR.LE.0) GO TO 240
      KQ=0
      DO 230 K=1,NLAY
      IF(IBQ(K).LE.0) GO TO 230
      KQ=KQ+1
      IF(IHCFM.LT.0) CALL ULAPRS(HC(1,1,KQ),TEXT(1,3),KSTP,KPER,NCOL,
     1          NROW,K,-IHCFM,IOUT)
      IF(IHCFM.GE.0) CALL ULAPRW(HC(1,1,KQ),TEXT(1,3),KSTP,KPER,NCOL,
     1           NROW,K,IHCFM,IOUT)
  230 CONTINUE
C
C10-----SAVE CRITICAL HEAD FOR ALL LAYERS WITH INTERBED STORAGE.
  240 IF(IHCSV.LE.0) GO TO 300
      KQ=0
      DO 260 K=1,NLAY
      IF(IBQ(K).LE.0) GO TO 260
      KQ=KQ+1
      CALL ULASAV(HC(1,1,KQ),TEXT(1,3),KSTP,KPER,PERTIM,TOTIM,NCOL,
     1            NROW,K,IHCUN)
  260 CONTINUE
C
C11-----RETURN
  300 RETURN
      END
