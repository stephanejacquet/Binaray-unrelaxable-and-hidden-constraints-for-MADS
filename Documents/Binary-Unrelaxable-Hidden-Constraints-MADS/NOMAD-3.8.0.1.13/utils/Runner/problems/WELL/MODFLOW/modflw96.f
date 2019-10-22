C     ******************************************************************
C     MAIN CODE FOR U.S. GEOLOGICAL SURVEY MODULAR MODEL -- MODFLOW-96
C           BY MICHAEL G. MCDONALD AND ARLEN W. HARBAUGH
C     MODFLOW-88 documented in:
C           McDonald, M.G. and Harbaugh, A.W., 1988, A modular
C           three-dimensional finite-difference ground-water flow
C           model: U.S. Geological Survey Techniques of Water
C           Resources Investigations, Book 6, Chapter A1, 586 p.
C     MODFLOW-96 documented in:
C           Harbaugh, A.W. and McDonald, M.G., 1996, User's
C           documentation for the U.S. Geological Survey modular
C           finite-difference ground-water flow model: U.S. Geological
C           Survey Open-File Report 96-485
C-----VERSION 0950 23MAY1996 MAIN
C-----VERSION 1401 03DEC1996 -- added PCG2, STR1, IBS1, CHD1, GFD1,
C                               HFB1, TLK1, DE45, and RES1 as documented
C                               in USGS reports
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
C1------SPECIFY THE SIZE OF THE X ARRAY.  TO CHANGE THE SIZE OF THE
C1------X ARRAY, CHANGE VALUE OF LENX IN THE NEXT STATEMENT.
      PARAMETER (LENX=1500000)
      COMMON X(LENX)
      COMMON /FLWCOM/LAYCON(200)
      CHARACTER*16 VBNM(40)
      CHARACTER*80 HEADNG(2)
      DIMENSION VBVL(4,40),IUNIT(40)
      DOUBLE PRECISION DUMMY
      EQUIVALENCE (DUMMY,X(1))
      CHARACTER*20 CHEDFM,CDDNFM
      CHARACTER*80 FNAME
      LOGICAL EXISTS
      CHARACTER*4 CUNIT(40)
      DATA CUNIT/'BCF ','WEL ','DRN ','RIV ','EVT ','TLK ','GHB ',
     1           'RCH ','SIP ','DE4 ','SOR ','OC  ','PCG ','GFD ',
     2           '    ','HFB ','RES ','STR ','IBS ','CHD ','FHB ',
     3           '    ','    ','    ','    ','    ','    ','    ',
     4           '    ','    ','    ','    ','    ','    ','    ',
     5           '    ','    ','    ','    ','    '/
C     ------------------------------------------------------------------
C     set string for use if RCS ident command
      FNAME =
     &'$Id: modflw96.f,v 3.2 1998/01/09 19:19:39 rsregan Exp rsregan $'
      FNAME = 
     &    '@(#)MODFLOW-96 - Modular 3-D Finite-Difference GW Flow Model'
      FNAME = '@(#)MODFLOW-96 - USGS TWRI, Book 6, Chap. A1, McDonald an
     &d Harbaugh'
      FNAME = '@(#)MODFLOW-96 - USGS OFR 96-485, Harbaugh and McDonald'
      FNAME = 
     &     '@(#)MODFLOW-96 - Contact: h2osoft@usgs.gov'
      FNAME = '@(#)MODFLOW-96 - Version: 3.0 1996/12/03 includes MOC'
      FNAME =
     &     '@(#)MODFLOW-96 - Version: 3.1 1997/03/11 fixed HFB1FM call'
      FNAME = '@(#)MODFLOW-96 - Version: 3.2x 1998/01/09 includes FHB'
C     ------------------------------------------------------------------
      INUNIT=99
      IBUNIT=98
      IBOUTS=97
      IBATCH=0
      INQUIRE(FILE='modflow.bf',EXIST=EXISTS)
      IF(EXISTS) THEN
         IBATCH=1
         OPEN(UNIT=IBUNIT,FILE='modflow.bf',STATUS='OLD')
         OPEN(UNIT=IBOUTS,FILE='modbatch.rpt')
         WRITE(IBOUTS,*) ' USGS MODFLOW MODEL BATCH-MODE REPORT'
      END IF
C
C2------OPEN FILE OF FILE NAMES.
50    IF(IBATCH.GT.0) THEN
         READ(IBUNIT,'(A)',END=500) FNAME
         IF(FNAME.EQ.' ') GO TO 50
         WRITE(IBOUTS,'(1X,/1X,A)') FNAME
      ELSE
         WRITE(*,*) ' Enter the name of the NAME FILE:'
         READ(*,'(A)') FNAME
      END IF
      INQUIRE(FILE=FNAME,EXIST=EXISTS)
      IF(.NOT.EXISTS) THEN
         IF(IBATCH.GT.0) THEN
            WRITE(IBOUTS,*) ' Specified name file does not exist.'
            WRITE(IBOUTS,*) ' Processing will continue with the next ',
     1                      'name file in modflow.bf.'
         ELSE
            WRITE(*,*) ' File does not exist'
         END IF
         GO TO 50
      END IF
      OPEN(UNIT=INUNIT,FILE=FNAME,STATUS='OLD')
C
C3------DEFINE PROBLEM--ROWS,COLUMNS,LAYERS,STRESS PERIODS,PACKAGES.
      CALL BAS5DF(ISUM,HEADNG,NPER,ITMUNI,TOTIM,NCOL,NROW,NLAY,
     1        NODES,INBAS,IOUT,IUNIT,CUNIT,INUNIT,IXSEC,ICHFLG,IFREFM)
C
C4------ALLOCATE SPACE IN "X" ARRAY.
      CALL BAS5AL(ISUM,LENX,LCHNEW,LCHOLD,LCIBOU,LCCR,LCCC,LCCV,
     1              LCHCOF,LCRHS,LCDELR,LCDELC,LCSTRT,LCBUFF,LCIOFL,
     2              INBAS,ISTRT,NCOL,NROW,NLAY,IOUT,IAPART,IFREFM)
      IF(IUNIT(1).GT.0) CALL BCF5AL(ISUM,LENX,LCSC1,LCHY,
     1     LCBOT,LCTOP,LCSC2,LCTRPY,IUNIT(1),ISS,
     2     NCOL,NROW,NLAY,IOUT,IBCFCB,LCWETD,IWDFLG,LCCVWD,
     3     WETFCT,IWETIT,IHDWET,HDRY,IAPART,IFREFM)
      IF(IUNIT(2).GT.0) CALL WEL5AL(ISUM,LENX,LCWELL,MXWELL,NWELLS,
     1                 IUNIT(2),IOUT,IWELCB,NWELVL,IWELAL,IFREFM)
      IF(IUNIT(3).GT.0) CALL DRN5AL(ISUM,LENX,LCDRAI,NDRAIN,MXDRN,
     1                 IUNIT(3),IOUT,IDRNCB,NDRNVL,IDRNAL,IFREFM)
      IF(IUNIT(4).GT.0) CALL RIV5AL(ISUM,LENX,LCRIVR,MXRIVR,NRIVER,
     1            IUNIT(4),IOUT,IRIVCB,NRIVVL,IRIVAL,IFREFM)
      IF(IUNIT(5).GT.0) CALL EVT5AL(ISUM,LENX,LCIEVT,LCEVTR,LCEXDP,
     1            LCSURF,NCOL,NROW,NEVTOP,IUNIT(5),IOUT,IEVTCB,IFREFM)
      IF(IUNIT(6).GT.0) CALL TLK1AL(ISUM,LENX,NCOL,NROW,NLAY,
     1          LCRAT,LCZCB,LCA1,LCB1,LCALPH,LCBET,LCRM1,LCRM2,LCRM3,
     2          LCRM4,LCTL,LCTLK,LCSLU,LCSLD,NODES1,NM1,NM2,NUMC,
     3          NTM1,ITLKSV,ITLKRS,ITLKCB,ISS,IUNIT(6),IOUT)
      IF(IUNIT(7).GT.0) CALL GHB5AL(ISUM,LENX,LCBNDS,NBOUND,MXBND,
     1            IUNIT(7),IOUT,IGHBCB,NGHBVL,IGHBAL,IFREFM)
      IF(IUNIT(8).GT.0) CALL RCH5AL(ISUM,LENX,LCIRCH,LCRECH,NRCHOP,
     1            NCOL,NROW,IUNIT(8),IOUT,IRCHCB,IFREFM)
      IF(IUNIT(9).GT.0) CALL SIP5AL(ISUM,LENX,LCEL,LCFL,LCGL,LCV,
     1          LCHDCG,LCLRCH,LCW,MXITER,NPARM,NCOL,NROW,NLAY,
     2          IUNIT(9),IOUT,IFREFM)
      IF(IUNIT(10).GT.0) CALL DE45AL(ISUM,LENX,LCAU,LCAL,LCIUPP,
     1           LCIEQP,LCD4B,LCLRCH,LCHDCG,
     2           MXUP,MXLOW,MXEQ,MXBW,IUNIT(10),ITMX,ID4DIR,
     3           NCOL,NROW,NLAY,IOUT,ID4DIM)
      IF(IUNIT(11).GT.0) CALL SOR5AL(ISUM,LENX,LCA,LCRES,LCHDCG,LCLRCH,
     1       LCIEQP,MXITER,NCOL,NLAY,NSLICE,MBW,IUNIT(11),IOUT,IFREFM)
      IF(IUNIT(13).GT.0) CALL PCG2AL(ISUM,LENX,LCV,LCSS,LCP,LCCD,
     1       LCHCHG,LCLHCH,LCRCHG,LCLRCH,MXITER,ITER1,NCOL,NROW,NLAY,
     2       IUNIT(13),IOUT,NPCOND,LCIT1)
      IF(IUNIT(14).GT.0) CALL GFD1AL(ISUM,LENX,LCSC1,LCCDTR,LCCDTC,
     1     LCBOT,LCTOP,LCSC2,IUNIT(14),ISS,NCOL,NROW,NLAY,IOUT,IGFDCB)
      IF(IUNIT(16).GT.0) CALL HFB1AL(ISUM,LENX,LCHFBR,NHFB,IUNIT(16),      *HFB*
     1           IOUT)                                                     *HFB*
      IF(IUNIT(17).GT.0) CALL RES1AL(ISUM,LENX,LCIRES,LCIRSL,LCBRES,
     1  LCCRES,LCBBRE,LCHRES,LCHRSE,IUNIT(17),IOUT,NRES,IRESCB,
     2  NRESOP,IRESPT,NPTS,NCOL,NROW)
      IF(IUNIT(18).GT.0) CALL STR1AL(ISUM,LENX,LCSTRM,ICSTRM,MXSTRM,    STR1
     1                 NSTREM,IUNIT(18),IOUT,ISTCB1,ISTCB2,NSS,NTRIB,   STR1
     2                  NDIV,ICALC,CONST,LCTBAR,LCTRIB,LCIVAR,LCFGAR)   STR1
      IF (IUNIT(19).GT.0) CALL IBS1AL(ISUM,LENX,LCHC,LCSCE,LCSCV,       IBS
     1           LCSUB,NCOL,NROW,NLAY,IIBSCB,IIBSOC,ISS,IUNIT(19),IOUT) IBS
      IF(IUNIT(20).GT.0) CALL CHD1AL(ISUM,LENX,LCCHDS,NCHDS,MXCHD,      CHD
     1           IUNIT(20),IOUT)                                        CHD
      IF(IUNIT(21).GT.0) CALL FHB1AL(ISUM,LENX,LCFLLC,LCBDTM,LCFLRT,
     1          LCBDFV,LCBDHV,LCHDLC,LCSBHD,NBDTIM,NFLW,NHED,IUNIT(21),
     2          IOUT,IFHBCB,NFHBX1,NFHBX2,IFHBD3,IFHBD4,IFHBD5,
     3          IFHBSS,ISS)
C
C5------IF THE "X" ARRAY IS NOT BIG ENOUGH THEN STOP.
      IF(ISUM-1.GT.LENX) STOP
C
C6------READ AND PREPARE INFORMATION FOR ENTIRE SIMULATION.
      CALL BAS5RP(X(LCIBOU),X(LCHNEW),X(LCSTRT),X(LCHOLD),
     1       ISTRT,INBAS,HEADNG,NCOL,NROW,NLAY,VBVL,X(LCIOFL),
     2       IUNIT(12),IHEDFM,IDDNFM,IHEDUN,IDDNUN,IOUT,IPEROC,ITSOC,
     3       CHEDFM,CDDNFM,IBDOPT,IXSEC,LBHDSV,LBDDSV,IFREFM)
      IF(IUNIT(1).GT.0) CALL BCF5RP(X(LCIBOU),X(LCHNEW),X(LCSC1),
     1          X(LCHY),X(LCCR),X(LCCC),X(LCCV),X(LCDELR),
     2     X(LCDELC),X(LCBOT),X(LCTOP),X(LCSC2),X(LCTRPY),IUNIT(1),
     3     ISS,NCOL,NROW,NLAY,IOUT,X(LCWETD),IWDFLG,X(LCCVWD))
      IF(IUNIT(6).GT.0) CALL TLK1RP(X(LCRAT),X(LCZCB),X(LCA1),X(LCB1),
     1          X(LCALPH),X(LCBET),X(LCRM1),X(LCRM2),X(LCRM3),X(LCRM4),
     2          NODES1,NM1,NM2,NUMC,NTM1,ITLKRS,DELTM1,X(LCBUFF),
     3          X(LCDELC),X(LCDELR),TLKTIM,NROW,NCOL,IUNIT(6),IOUT)
      IF(IUNIT(9).GT.0) CALL SIP5RP(NPARM,MXITER,ACCL,HCLOSE,X(LCW),
     1          IUNIT(9),IPCALC,IPRSIP,IOUT,IFREFM)
      IF(IUNIT(10).GT.0) CALL DE45RP(IUNIT(10),MXITER,NITER,ITMX,
     1            ACCL,HCLOSE,IFREQ,IPRD4,IOUT,MUTD4)
      IF(IUNIT(11).GT.0) CALL SOR5RP(MXITER,ACCL,HCLOSE,IUNIT(11),
     1         IPRSOR,IOUT,IFREFM)
      IF(IUNIT(13).GT.0) CALL PCG2RP(MXITER,ITER1,HCLOSE,RCLOSE,
     1         NPCOND,NBPOL,RELAX,IPRPCG,IUNIT(13),IOUT,MUTPCG,
     2         NITER,X(LCIT1),DAMP)
      IF(IUNIT(14).GT.0) CALL GFD1RP(X(LCIBOU),X(LCHNEW),X(LCSC1),
     1          X(LCCDTR),X(LCCDTC),X(LCCR),X(LCCC),X(LCCV),X(LCDELR),
     2          X(LCDELC),X(LCBOT),X(LCTOP),X(LCSC2),
     3          IUNIT(14),ISS,NCOL,NROW,NLAY,NODES,IOUT)
      IF(IUNIT(16).GT.0) CALL HFB1RP(X(LCCR),X(LCCC),X(LCDELR),            *HFB*
     1         X(LCDELC),X(LCHFBR),IUNIT(16),NCOL,NROW,NLAY,NODES,         *HFB*
     1         NHFB,IOUT)                                                  *HFB*
      IF(IUNIT(19).GT.0) CALL IBS1RP(X(LCDELR),X(LCDELC),X(LCHNEW),     IBS
     1      X(LCHC),X(LCSCE),X(LCSCV),X(LCSUB),NCOL,NROW,NLAY,          IBS
     2      NODES,IIBSOC,ISUBFM,ICOMFM,IHCFM,ISUBUN,ICOMUN,IHCUN,       IBS
     3      IUNIT(19),IOUT)                                             IBS
      IF(IUNIT(21).GT.0) CALL FHB1RP(X(LCIBOU),NROW,NCOL,NLAY,
     &          X(LCFLLC),X(LCBDTM),NBDTIM,X(LCFLRT),NFLW,NHED,
     &          X(LCHDLC),X(LCSBHD),IUNIT(21),IOUT,
     &          NFHBX1,NFHBX2,IFHBD3,IFHBD5)
C
C7------SIMULATE EACH STRESS PERIOD.
      DO 300 KPER=1,NPER
      KKPER=KPER
C
C7A-----READ STRESS PERIOD TIMING INFORMATION.
      CALL BAS5ST(NSTP,DELT,TSMULT,PERTIM,KKPER,INBAS,IOUT,IFREFM)
C
C7B-----READ AND PREPARE INFORMATION FOR STRESS PERIOD.
      IF(IUNIT(2).GT.0) CALL WEL5RP(X(LCWELL),NWELLS,MXWELL,IUNIT(2),
     1             IOUT,NWELVL,IWELAL,IFREFM)
      IF(IUNIT(3).GT.0) CALL DRN5RP(X(LCDRAI),NDRAIN,MXDRN,IUNIT(3),
     1                 IOUT,NDRNVL,IDRNAL,IFREFM)
      IF(IUNIT(4).GT.0) CALL RIV5RP(X(LCRIVR),NRIVER,MXRIVR,IUNIT(4),
     1            IOUT,NRIVVL,IRIVAL,IFREFM)
      IF(IUNIT(5).GT.0) CALL EVT5RP(NEVTOP,X(LCIEVT),X(LCEVTR),
     1            X(LCEXDP),X(LCSURF),X(LCDELR),X(LCDELC),NCOL,NROW,
     1            IUNIT(5),IOUT,IFREFM)
      IF(IUNIT(7).GT.0) CALL GHB5RP(X(LCBNDS),NBOUND,MXBND,IUNIT(7),
     1              IOUT,NGHBVL,IGHBAL,IFREFM)
      IF(IUNIT(8).GT.0) CALL RCH5RP(NRCHOP,X(LCIRCH),X(LCRECH),
     1            X(LCDELR),X(LCDELC),NROW,NCOL,IUNIT(8),IOUT,IFREFM)
      IF(IUNIT(17).GT.0) CALL RES1RP(X(LCIRES),X(LCIRSL),X(LCBRES),
     1   X(LCCRES),X(LCBBRE),X(LCHRSE),X(LCIBOU),X(LCDELR),X(LCDELC),
     2   NRES,NRESOP,NPTS,NCOL,NROW,NLAY,PERLEN,DELT,NSTP,TSMULT,
     3   IUNIT(17),IOUT)
      IF(IUNIT(18).GT.0) CALL STR1RP(X(LCSTRM),X(ICSTRM),NSTREM,        STR1
     1                     MXSTRM,IUNIT(18),IOUT,X(LCTBAR),NDIV,NSS,    STR1
     2                     NTRIB,X(LCIVAR),ICALC,IPTFLG)                STR1
      IF(IUNIT(20).GT.0) CALL CHD1RP(X(LCCHDS),NCHDS,MXCHD,X(LCIBOU),   CHD
     1            NCOL,NROW,NLAY,PERLEN,DELT,NSTP,TSMULT,IUNIT(20),IOUT)CHD
C
C7C-----SIMULATE EACH TIME STEP.
      DO 200 KSTP=1,NSTP
      KKSTP=KSTP
C
C7C1----CALCULATE TIME STEP LENGTH. SET HOLD=HNEW..
      CALL BAS5AD(DELT,TSMULT,TOTIM,PERTIM,X(LCHNEW),X(LCHOLD),KKSTP,
     1             NCOL,NROW,NLAY)
      IF(IUNIT(6).GT.0) CALL TLK1AD(X(LCRAT),X(LCZCB),X(LCA1),X(LCB1),
     1          X(LCALPH),X(LCBET),X(LCRM1),X(LCRM2),X(LCRM3),X(LCRM4),
     2          X(LCTL),X(LCTLK),X(LCSLU),X(LCSLD),NM1,NM2,NUMC,NTM1,
     3          DELTM1,X(LCHNEW),X(LCIBOU),X(LCTOP),
     4          NROW,NCOL,NLAY,DELT,TLKTIM,IUNIT(6),IOUT)
      IF(IUNIT(20).GT.0) CALL CHD1FM(NCHDS,MXCHD,X(LCCHDS),X(LCIBOU),   CHD
     1          X(LCHNEW),X(LCHOLD),PERLEN,PERTIM,DELT,NCOL,NROW,NLAY)  CHD
      IF(IUNIT(1).GT.0) CALL BCF5AD(X(LCIBOU),X(LCHOLD),X(LCBOT),
     1             X(LCWETD),IWDFLG,ISS,NCOL,NROW,NLAY)
      IF(IUNIT(17).GT.0) CALL RES1AD(X(LCHRES),X(LCHRSE),X(LCIRES),
     1 X(LCBRES),X(LCDELR),X(LCDELC),NRES,IRESPT,NCOL,NROW,
     1      PERLEN,PERTIM,TOTIM,KKSTP,KKPER,IOUT)
      IF(IUNIT(21).GT.0) CALL FHB1AD(X(LCHNEW),X(LCHOLD),NCOL,NROW,NLAY,
     &          ISS,TOTIM,DELT,X(LCBDTM),NBDTIM,X(LCFLRT),
     &          X(LCBDFV),X(LCBDHV),NFLW,X(LCSBHD),X(LCHDLC),NHED,
     &          NFHBX1,NFHBX2,IFHBD3,IFHBD4,IFHBD5,IFHBSS)
C
C7C2----ITERATIVELY FORMULATE AND SOLVE THE EQUATIONS.
      DO 100 KITER=1,MXITER
      KKITER=KITER
C
C7C2A---FORMULATE THE FINITE DIFFERENCE EQUATIONS.
      CALL BAS5FM(X(LCHCOF),X(LCRHS),NODES)
      IF(IUNIT(1).GT.0) CALL BCF5FM(X(LCHCOF),X(LCRHS),X(LCHOLD),
     1          X(LCSC1),X(LCHNEW),X(LCIBOU),X(LCCR),X(LCCC),X(LCCV),
     2          X(LCHY),X(LCTRPY),X(LCBOT),X(LCTOP),X(LCSC2),
     3          X(LCDELR),X(LCDELC),DELT,ISS,KKITER,KKSTP,KKPER,NCOL,
     4          NROW,NLAY,IOUT,X(LCWETD),IWDFLG,X(LCCVWD),WETFCT,
     5          IWETIT,IHDWET,HDRY,X(LCBUFF))
      IF(IUNIT(14).GT.0) CALL GFD1FM(X(LCHCOF),X(LCRHS),X(LCHOLD),
     1          X(LCSC1),X(LCHNEW),X(LCIBOU),X(LCCR),X(LCCC),X(LCCV),
     2          X(LCCDTR),X(LCCDTC),X(LCBOT),X(LCTOP),X(LCSC2),
     3          DELT,ISS,KKITER,KKSTP,KKPER,NCOL,NROW,NLAY,IOUT)
      IF(IUNIT(16).GT.0) CALL HFB1FM(X(LCHNEW),X(LCCR),X(LCCC),            *HFB*
     1          X(LCBOT),X(LCTOP),X(LCDELR),X(LCDELC),X(LCHFBR),           *HFB*
     2          NCOL,NROW,NLAY,NHFB)                                       *HFB*
      IF(IUNIT(6).GT.0) CALL TLK1FM(X(LCRAT),X(LCTL),X(LCTLK),X(LCSLU),
     1          X(LCSLD),NUMC,X(LCHNEW),X(LCIBOU),X(LCTOP),X(LCCV),
     2          X(LCHCOF),X(LCRHS),NROW,NCOL,NLAY)
      IF(IUNIT(2).GT.0) CALL WEL5FM(NWELLS,MXWELL,X(LCRHS),X(LCWELL),
     1           X(LCIBOU),NCOL,NROW,NLAY,NWELVL)
      IF(IUNIT(3).GT.0) CALL DRN5FM(NDRAIN,MXDRN,X(LCDRAI),X(LCHNEW),
     1         X(LCHCOF),X(LCRHS),X(LCIBOU),NCOL,NROW,NLAY,NDRNVL)
      IF(IUNIT(4).GT.0) CALL RIV5FM(NRIVER,MXRIVR,X(LCRIVR),X(LCHNEW),
     1            X(LCHCOF),X(LCRHS),X(LCIBOU),NCOL,NROW,NLAY,NRIVVL)
      IF(IUNIT(5).GT.0) CALL EVT5FM(NEVTOP,X(LCIEVT),X(LCEVTR),
     1            X(LCEXDP),X(LCSURF),X(LCRHS),X(LCHCOF),X(LCIBOU),
     1            X(LCHNEW),NCOL,NROW,NLAY)
      IF(IUNIT(7).GT.0) CALL GHB5FM(NBOUND,MXBND,X(LCBNDS),X(LCHCOF),
     1            X(LCRHS),X(LCIBOU),NCOL,NROW,NLAY,NGHBVL)
      IF(IUNIT(8).GT.0) CALL RCH5FM(NRCHOP,X(LCIRCH),X(LCRECH),
     1            X(LCRHS),X(LCIBOU),NCOL,NROW,NLAY)
      IF(IUNIT(17).GT.0) CALL RES1FM(X(LCIRES),X(LCIRSL),X(LCBRES),
     1   X(LCCRES),X(LCBBRE),X(LCHRES),X(LCIBOU),X(LCHNEW),X(LCHCOF),
     2   X(LCRHS),NRES,NRESOP,NCOL,NROW,NLAY)
      IF(IUNIT(18).GT.0) CALL STR1FM(NSTREM,X(LCSTRM),X(ICSTRM),        STR1
     1                     X(LCHNEW),X(LCHCOF),X(LCRHS),X(LCIBOU),      STR1
     2              MXSTRM,NCOL,NROW,NLAY,IOUT,NSS,X(LCTBAR),           STR1
     3              NTRIB,X(LCTRIB),X(LCIVAR),X(LCFGAR),ICALC,CONST)    STR1
      IF(IUNIT(19).GT.0) CALL IBS1FM(X(LCRHS),X(LCHCOF),X(LCHNEW),      IBS
     1       X(LCHOLD),X(LCHC),X(LCSCE),X(LCSCV),X(LCIBOU),             IBS
     2       NCOL,NROW,NLAY,DELT)                                       IBS
      IF(IUNIT(21).GT.0) CALL FHB1FM(X(LCRHS),X(LCIBOU),X(LCFLLC),
     1 X(LCBDFV),NFLW,NCOL,NROW,NLAY,IFHBD4)
C
C7C2B---MAKE ONE CUT AT AN APPROXIMATE SOLUTION.
      IF(IUNIT(9).GT.0) CALL SIP5AP(X(LCHNEW),X(LCIBOU),X(LCCR),X(LCCC),
     1     X(LCCV),X(LCHCOF),X(LCRHS),X(LCEL),X(LCFL),X(LCGL),X(LCV),
     2     X(LCW),X(LCHDCG),X(LCLRCH),NPARM,KKITER,HCLOSE,ACCL,ICNVG,
     3     KKSTP,KKPER,IPCALC,IPRSIP,MXITER,NSTP,NCOL,NROW,NLAY,NODES,
     4     IOUT)
      IF(IUNIT(10).GT.0) CALL DE45AP(X(LCHNEW),X(LCIBOU),X(LCAU),
     1  X(LCAL),X(LCIUPP),X(LCIEQP),X(LCD4B),MXUP,MXLOW,MXEQ,MXBW,
     2  X(LCCR),X(LCCC),X(LCCV),X(LCHCOF),X(LCRHS),ACCL,KKITER,ITMX,
     3  MXITER,NITER,HCLOSE,IPRD4,ICNVG,NCOL,NROW,NLAY,IOUT,X(LCLRCH),
     4  X(LCHDCG),IFREQ,KKSTP,KKPER,DELT,NSTP,ID4DIR,ID4DIM,MUTD4)
      IF(IUNIT(11).GT.0) CALL SOR5AP(X(LCHNEW),X(LCIBOU),X(LCCR),
     1     X(LCCC),X(LCCV),X(LCHCOF),X(LCRHS),X(LCA),X(LCRES),X(LCIEQP),
     2     X(LCHDCG),X(LCLRCH),KKITER,HCLOSE,ACCL,ICNVG,KKSTP,KKPER,
     3     IPRSOR,MXITER,NSTP,NCOL,NROW,NLAY,NSLICE,MBW,IOUT)
      IF(IUNIT(13).GT.0) CALL PCG2AP(X(LCHNEW),X(LCIBOU),X(LCCR),
     1      X(LCCC),X(LCCV),X(LCHCOF),X(LCRHS),X(LCV),X(LCSS),X(LCP),
     2      X(LCCD),X(LCHCHG),X(LCLHCH),X(LCRCHG),X(LCLRCH),KKITER,
     3      NITER,HCLOSE,RCLOSE,ICNVG,KKSTP,KKPER,IPRPCG,MXITER,ITER1,
     4      NPCOND,NBPOL,NSTP,NCOL,NROW,NLAY,NODES,RELAX,IOUT,MUTPCG,
     5      0,0,SN,SP,SR,X(LCIT1),DAMP)
C
C7C2C---IF CONVERGENCE CRITERION HAS BEEN MET STOP ITERATING.
      IF(ICNVG.EQ.1) GO TO 110
  100 CONTINUE
      KITER=MXITER
  110 CONTINUE
C
C7C3----DETERMINE WHICH OUTPUT IS NEEDED.
      CALL BAS5OC(NSTP,KKSTP,ICNVG,X(LCIOFL),NLAY,IBUDFL,ICBCFL,
     1   IHDDFL,IUNIT(12),IOUT,KKPER,IPEROC,ITSOC,IBDOPT,IXSEC,IFREFM)
C
C7C4----CALCULATE BUDGET TERMS. SAVE CELL-BY-CELL FLOW TERMS.
      MSUM=1
      IF(IUNIT(6).GT.0) CALL TLK1BD(X(LCRAT),X(LCTL),X(LCTLK),
     1          X(LCSLU),X(LCSLD),NUMC,ITLKCB,X(LCHNEW),X(LCBUFF),
     2          X(LCIBOU),X(LCTOP),X(LCCV),VBNM,VBVL,MSUM,NCOL,NROW,
     3          NLAY,DELT,KSTP,KPER,ICBCFL,IOUT)
C7C4A---THE ORIGINAL BCF BUDGET MODULE HAS BEEN REPLACED BY THREE
C7C4A---SUBMODULES: SBCF5S, SBCF5F, AND SBCF5B .
      IF(IUNIT(1).GT.0) THEN
         CALL SBCF5S(VBNM,VBVL,MSUM,X(LCHNEW),X(LCIBOU),X(LCHOLD),
     1     X(LCSC1),X(LCTOP),X(LCSC2),DELT,ISS,NCOL,NROW,NLAY,KKSTP,
     2     KKPER,IBCFCB,ICBCFL,X(LCBUFF),IOUT,PERTIM,TOTIM)
         CALL SBCF5F(VBNM,VBVL,MSUM,X(LCHNEW),X(LCIBOU),X(LCCR),
     1     X(LCCC),X(LCCV),X(LCTOP),DELT,NCOL,NROW,NLAY,KKSTP,KKPER,
     2     IBCFCB,X(LCBUFF),IOUT,ICBCFL,PERTIM,TOTIM,ICHFLG)
         IBDRET=0
         IC1=1
         IC2=NCOL
         IR1=1
         IR2=NROW
         IL1=1
         IL2=NLAY
         DO 155 IDIR=1,3
         CALL SBCF5B(X(LCHNEW),X(LCIBOU),X(LCCR),X(LCCC),X(LCCV),
     1      X(LCTOP),NCOL,NROW,NLAY,KKSTP,KKPER,IBCFCB,X(LCBUFF),
     2      IOUT,ICBCFL,DELT,PERTIM,TOTIM,IDIR,IBDRET,ICHFLG,
     3      IC1,IC2,IR1,IR2,IL1,IL2)
155      CONTINUE
      END IF
      IF(IUNIT(14).GT.0) CALL GFD1BD(VBNM,VBVL,MSUM,X(LCHNEW),
     1     X(LCIBOU),X(LCHOLD),X(LCSC1),X(LCCR),X(LCCC),X(LCCV),
     2     X(LCTOP),X(LCSC2),DELT,ISS,NCOL,NROW,NLAY,KKSTP,KKPER,
     3     IGFDCB,ICBCFL,X(LCBUFF),IOUT)
      IF(IUNIT(2).GT.0) CALL WEL5BD(NWELLS,MXWELL,VBNM,VBVL,MSUM,
     1     X(LCWELL),X(LCIBOU),DELT,NCOL,NROW,NLAY,KKSTP,KKPER,IWELCB,
     1     ICBCFL,X(LCBUFF),IOUT,PERTIM,TOTIM,NWELVL,IWELAL)
      IF(IUNIT(3).GT.0) CALL DRN5BD(NDRAIN,MXDRN,VBNM,VBVL,MSUM,
     1     X(LCDRAI),DELT,X(LCHNEW),NCOL,NROW,NLAY,X(LCIBOU),KKSTP,
     2     KKPER,IDRNCB,ICBCFL,X(LCBUFF),IOUT,PERTIM,TOTIM,NDRNVL,
     3     IDRNAL)
      IF(IUNIT(4).GT.0) CALL RIV5BD(NRIVER,MXRIVR,X(LCRIVR),X(LCIBOU),
     1     X(LCHNEW),NCOL,NROW,NLAY,DELT,VBVL,VBNM,MSUM,KKSTP,KKPER,
     2     IRIVCB,ICBCFL,X(LCBUFF),IOUT,PERTIM,TOTIM,NRIVVL,IRIVAL)
      IF(IUNIT(5).GT.0) CALL EVT5BD(NEVTOP,X(LCIEVT),X(LCEVTR),
     1     X(LCEXDP),X(LCSURF),X(LCIBOU),X(LCHNEW),NCOL,NROW,NLAY,
     2     DELT,VBVL,VBNM,MSUM,KKSTP,KKPER,IEVTCB,ICBCFL,X(LCBUFF),IOUT,
     3     PERTIM,TOTIM)
      IF(IUNIT(7).GT.0) CALL GHB5BD(NBOUND,MXBND,VBNM,VBVL,MSUM,
     1     X(LCBNDS),DELT,X(LCHNEW),NCOL,NROW,NLAY,X(LCIBOU),KKSTP,
     2     KKPER,IGHBCB,ICBCFL,X(LCBUFF),IOUT,PERTIM,TOTIM,NGHBVL,
     3     IGHBAL)
      IF(IUNIT(8).GT.0) CALL RCH5BD(NRCHOP,X(LCIRCH),X(LCRECH),
     1     X(LCIBOU),NROW,NCOL,NLAY,DELT,VBVL,VBNM,MSUM,KKSTP,KKPER,
     2     IRCHCB,ICBCFL,X(LCBUFF),IOUT,PERTIM,TOTIM)
      IF(IUNIT(17).GT.0) CALL RES1BD(X(LCIRES),X(LCIRSL),X(LCBRES),
     1      X(LCCRES),X(LCBBRE),X(LCHRES),X(LCIBOU),X(LCHNEW),
     2      X(LCBUFF),VBVL,VBNM,MSUM,KSTP,KPER,NRES,NRESOP,
     3      NCOL,NROW,NLAY,DELT,IRESCB,ICBCFL,IOUT)
      IF(IUNIT(18).GT.0) CALL STR1BD(NSTREM,X(LCSTRM),X(ICSTRM),        STR1
     1   X(LCIBOU),MXSTRM,X(LCHNEW),NCOL,NROW,NLAY,DELT,VBVL,VBNM,MSUM, STR1
     2   KKSTP,KKPER,ISTCB1,ISTCB2,ICBCFL,X(LCBUFF),IOUT,NTRIB,NSS,     STR1
     3   X(LCTRIB),X(LCTBAR),X(LCIVAR),X(LCFGAR),ICALC,CONST,IPTFLG)    STR1
      IF(IUNIT(19).GT.0) CALL IBS1BD(X(LCIBOU),X(LCHNEW),X(LCHOLD),     IBS
     1      X(LCHC),X(LCSCE),X(LCSCV),X(LCSUB),X(LCDELR),X(LCDELC),     IBS
     2      NCOL,NROW,NLAY,DELT,VBVL,VBNM,MSUM,KSTP,KPER,IIBSCB,        IBS
     3      ICBCFL,X(LCBUFF),IOUT)                                      IBS
      IF(IUNIT(21).GT.0) CALL FHB1BD(X(LCFLLC),X(LCBDFV),NFLW,
     1     VBNM,VBVL,MSUM,X(LCIBOU),DELT,NCOL,NROW,NLAY,KKSTP,KKPER,
     2     IFHBCB,ICBCFL,X(LCBUFF),IOUT,IFHBD4)
C
C7C5---PRINT AND OR SAVE HEADS AND DRAWDOWNS. PRINT OVERALL BUDGET.
      CALL BAS5OT(X(LCHNEW),X(LCSTRT),ISTRT,X(LCBUFF),X(LCIOFL),
     1     MSUM,X(LCIBOU),VBNM,VBVL,KKSTP,KKPER,DELT,PERTIM,TOTIM,
     2     ITMUNI,NCOL,NROW,NLAY,ICNVG,IHDDFL,IBUDFL,IHEDFM,IHEDUN,
     3     IDDNFM,IDDNUN,IOUT,CHEDFM,CDDNFM,IXSEC,LBHDSV,LBDDSV)
C
C7C5A--PRINT AND OR SAVE SUBSIDENCE, COMPACTION, AND CRITICAL HEAD.
      IF(IUNIT(19).GT.0) CALL IBS1OT(NCOL,NROW,NLAY,PERTIM,TOTIM,KSTP,  IBS
     1      KPER,NSTP,X(LCBUFF),X(LCSUB),X(LCHC),IIBSOC,ISUBFM,ICOMFM,  IBS
     2      IHCFM,ISUBUN,ICOMUN,IHCUN,IUNIT(19),IOUT)                   IBS
C
C7C6----IF ITERATION FAILED TO CONVERGE THEN STOP.
      IF(ICNVG.EQ.0) STOP
  200 CONTINUE
  300 CONTINUE
C
C7C7----WRITE RESTART RECORDS
C7C7A---WRITE RESTART RECORDS FOR TRANSIENT-LEAKAGE PACKAGE
      IF(IUNIT(6).GT.0) CALL TLK1OT(X(LCRM1),X(LCRM2),
     1     X(LCRM3),X(LCRM4),NM1,NM2,ITLKSV,DELTM1,TLKTIM,IOUT)
C
C8------END OF SIMULATION
      IF(IBATCH.GT.0) THEN
         WRITE(IBOUTS,*) ' Normal termination of simulation.'
         DO 400 I=1,IBOUTS-1
            INQUIRE(UNIT=I,OPENED=EXISTS)
            IF(EXISTS) CLOSE(I)
  400    CONTINUE
         GO TO 50
      END IF
  500 STOP
C
      END
