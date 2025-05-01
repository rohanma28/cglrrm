      Subroutine InitOntMain(StartYear, StartMonth, StartDay)
C
c This routine is called by "MISCUTIL.FOR" and is ONLY called once.  
c Therefore all initializations should be performed here using the 
c variable "UStLawrSolution" to determine which Lake Ontario 
c module is called:
c
C     1=outflow eqn, 2=1958D, 3=IS Plan, 4=Plan 1998, 5=Preproject
c
c
C This subroutine exits to kluge Paul's initialization
C data in COMMON blocks with the stuff that Matt set
C up for reading them in.  Someday it should be better
C integrated.  It gets called from CheckPar().
c
c the Verbosity outputs: V5 = verbosity level 5
c                      : V9 = verbosity level 9
C  
      implicit integer*4 (i-n)
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      integer   StartYear, StartMonth, StartDay
C     Apparently unused variables removed by mmm 13Mar09 to hush compiler
C     integer*4 fmon,fday,fyr,totday,elim
      integer*4 fmon,fday,fyr,       elim
      integer*4 nlimMon(20),nlimDay(20),nlimYear(20),nlim
C     Apparently unused variables removed by mmm 13Mar09 to hush compiler
C     integer*4 day,yr,daye,yre,md(8),jd(8)
      integer*4 day,yr,         md(8),jd(8)
      integer*4 isea(12,31),iwt(12,31),ksi(8),jsi(13)
C     Apparently unused variables removed by mmm 13Mar09 to hush compiler
C     integer*4 q_stlaw_ont_diff,q_stlouis,q_desprairies
C     integer*4 q_milleiles,ice_retard
C     integer*4 isup,iflo,qdev,erieout,Locsup
      integer*4 isup,iflo,qdev

      character*2 sympx
      integer   msi(21)
      logical   ElimDev,ForceFlow,SetFlow,UseDatabase,SiDatabase
      logical   DCompute
c
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
      CHARACTER    CTMP*255, CWRK*80, LLT*1
      CHARACTER*2  P58DAppLimitSym
      INTEGER      P58DForecastYear, P58DForecastMon, P58DForecastDay 
      INTEGER      CritDepths,P58DTotalSup,MaxILimit,MinPLimit
      INTEGER      P58DRunType,P58DWeightedSupply,P58DRuleCurve
      INTEGER      P58DAppLimit,SupplyIndicator(21)
      Integer      P58DAdjustment,P58DrcAdjustment,P58DActualFlow
      integer      P58DTotalDev,P58DLSLflow,prevstlou
C     Some compilers squawk about common blocks that don't total a multiple of
C     8 bytes, so use ipad to satisfy them when necessary.
      integer*4 ipd, ipe
      DOUBLE PRECISION P58DActualLevel, P58DPlanLevel, Onstlv1
      DOUBLE PRECISION PreprojectLevel
      INTEGER      MidLakeSolution, BOMRound , FlowRound, LevelRound,
     &             RoutingInterval(3), RoutingIncrement(3),
     &             StMarysSolution, UStLawrSolution, P77ACompat
      CHARACTER*256 P58DSimulateFile, P58DSIDataFile
      CHARACTER*256 OutFileName, filevar, CritKfile
      LOGICAL      ElimDev1,ForceFlow1,SetFlow1,PreProjConditions
      logical      P58DUseDatabase, P58DSiDatabase, ComputeMontreal
      COMMON/MISINT/MidLakeSolution, BOMRound, FlowRound, LevelRound,
     &              RoutingInterval, RoutingIncrement,
     &              StMarysSolution, UStLawrSolution, P77ACompat
      common/DComp/DCompute
      common/initpass/Elimdev,ForceFlow,SetFlow,UseDatabase,SiDatabase
      common/ae/nlimMon,nlimDay,nlimYear,nlim
      common/az/fmon,fday,fyr
      common/d1/mon,day,yr,md,jd
      common/d2/xdislev,AveOntLev
      common/p1/ipd,nm1,nd1,ny1,plevc,dlevc,itotdev,nwt1785,nwts,limplan
      common/p2/iasi,isi,icsi,elim
      common/q1/isea,iwt,isup,iflo,qdev
      common/r1/iprnt,kcrt,iforce,isetflow
      common/t1/ipe,nm,nd,ny,plev,dlev,iws,iws1,nrcp,nrcd,lasi
C     MMM 10Oct2007 - Block t2 is different size here than in Ontario.for - lmax and lmin seem unused
C      common/t2/iadj,lmax,lmin,n13,ksi,jsi
      common/t2/iadj,n13,ksi,jsi
      common/precond1/beglev
      common/OntBeglv/onstlv1
      common/r2/prevstlou
      common/maxmin/maxilimit0,minplimit0
      COMMON /CritK/  CritKfile
C
      COMMON /INI58D/ CritDepths, P58DTotalSup, MaxILimit, MinPLimit,
     &                P58DWeightedSupply, P58DRuleCurve, P58DAppLimit,
     &                P58DRunType, P58DForecastYear, P58DForecastMon,
     &                P58DForecastDay, SupplyIndicator, 
     &                P58DPlanLevel, P58DActualLevel,
     &                P58DSimulateFile, P58DSIDataFile, ElimDev1, 
     &                ForceFlow1, SetFlow1, PreProjConditions,
     &                P58DUseDatabase, P58DSiDatabase
      COMMON /Ini58dc/P58DAppLimitSym
      COMMON /INI58D1/P58DAdjustment,P58DrcAdjustment,P58DActualFlow, 
     &                P58DTotalDev,P58DLSLflow
      COMMON /INI58D2/ComputeMontreal
      common /PreCond/PreprojectLevel
      COMMON /RCPARI/ GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /LOGCOD/ LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/ CTMP,LLT,CWRK
C
      mon=StartMonth
      day=StartDay
      yr =StartYear
      fmon=mon
      fday=day
      fyr=yr
      DCompute=ComputeMontreal
      if(ontverbosity.ge.5)then
          filevar=outfilename('VOntario',8,ilen)
          open(753,file=filevar)
          filevar=outfilename('VStlaw_Up',9,ilen)
          open(759,file=filevar)
          write(759,106)
          if (ComputeMontreal) then
              filevar=outfilename('VStlaw_Dn',9,ilen)
              open(760,file=filevar)
              write(760,107)
          endif
          if(ontverbosity.ge.9)then
              filevar=outfilename('Vdetails',8,ilen)
              open(812,file=filevar)
          endif
      endif
      if(UStLawrSolution .EQ. 2)then
c
c Choosing Plan 1958-D
c
          if(ontverbosity.ge.5)then
              write(753,104)
          endif
          elim=0
          if(ElimDev)then
c
c     The total number of dates is set at 20 (see dimension statements)
c
              elim=1
              i=1
              open(450,file=CritKfile)
              read(450,*)
              read(450,*)
              do while (i.ge.1)
                  read(450,*,end=10)nlimMon(i),nlimDay(i),nlimYear(i)
                  i=i+1
              end do
   10         nlim=i-1
              close(450)
          endif
          iforce=0
          if(ForceFlow)then
              iforce=1
          endif
          isetflow=0
          if(SetFlow)then
              isetflow=1
          endif
          kcrt=CritDepths
          isupcx=P58DTotalSup * .1
          iwsx=P58DWeightedSupply * .1
          iadjx=P58DAdjustment * .1
          iadjrc=P58DrcAdjustment * .1
          nrcx=0
          plevx=P58DPlanLevel
          dlevx=onstlv1
          limplnx=P58DAppLimit *.1
          sympx=P58DAppLimitSym
          limdisx=P58DActualFlow *.1
          itotdev=P58DTotalDev * .1
          prevstlou=P58DLSLflow
          maxilimit0=maxilimit*.1
          minplimit0=minplimit*.1
          do i=2,21
              msi(i)=SupplyIndicator(i)
          end do
c
          isup=0
          iflo=0
          qdev=0.
c
c     Getting the First line of Data
c
          call firstline(UseDatabase,SiDatabase,
     +        P58DSIDataFile(1:LastKC(P58DSIDataFile,256)),
     +        isupcx,iwsx,nrcx,iadjx,iadjrc,plevx,dlevx,
     +        limplnx,limdisx,sympx,msi)
      Elseif(UStLawrSolution .EQ. 5)then
c
c Choosing Preproject
c
          if(ontverbosity.ge.5)then
              write(753,201)
              write(753,202)
          endif
c
c     Preproject Initial Conditions
c
          beglev=PreprojectLevel
c
c     variable "iforce" must be set to 1 in order force the preproject
c     flow thru the St. Lawrence River.
          iforce=1
c
c     Must get the initial Preproject level
          adj=.0017*(ny-1903)
          xdislev=beglev+adj
c
c     Getting the next week's date, since the start date (fmon,fday,fyr)
c     is a Saturday, you only need to get the following Friday's date
c     which is the 7th day (1st day is current day).
          call date10(fmon,fday,fyr,jd,md,7,'f')
          nm=md(7)
          nd=jd(7)
          ny=fyr
          if(nm.eq.12.and.md(7).eq.1)ny=ny+1
      endif
      return
  104 format(37x,'Lake Ontario Regulation - Plan 1958-D',//,
     + 12x,'lake',21x,'chan',6x,'adj. basc',
     + 30x,'eow',10x,'dev.',14x,'eow',/,
     + 12x,'ont. 1785',6x,'norm sup.  in',7x,
     + 'sup. rule  sea. sea. app.    delta    plan    dis.  dis. ',
     + 'tot. delta   dis.',/,
     + 4x,'week    sply   x  wt.s wt.s ind. ind. adjm ind. curv  adjm',
     + ' adjm lim.   storage   level   flow  flow dev. stage  level',/,
     + '   ending   cms1 wt.s',7(' cms1'),2x,'cms1',2(' cms1'),
     + '  cms1 metre metre   cms1  cms1 cms1 metre  metre  ',/)
  106 format(4x,'Week',7x,'Plan 1958-D, Upstream of Montreal (metres)',
     + /,1x,'  Ending    Out  L.Ont  Ogden   Card   Iroq Morris LSault',
     + /,1x,'----------  --- ------ ------ ------ ------ ------ ------')
  107 format(4x,'Week',7x,'Plan 1958-D, Downstream Montreal (metres)',
     + /,1x,'  Ending    Out Jetty1 Varen  Sorel  Pierre 3River Batis',
     + /,1x,'----------  --- ------ ------ ------ ------ ------ ------')
  201 format(22x,'Preproject Computations')
  202 format(9x,'Date   Beglev Supply Outflo   Diff   Difm Endlev',/,
     +    6x,'--------- ------ ------ ------ ------ ------ ------')
      end
      SUBROUTINE OntarioPlans(SDayOfMon, EDayOfMon, QMon, Monn, Year, 
     &    DaysElapsed, NewWeek, NewQMon, NewMon, SecsInRtPeriod, 
     &    NiagaraConverted, FlowFromWelland, SupplyON, ONNTS, 
     &    OntRetard, OntDeviation, StLouisRetard, StLouisFlow,  
     &    LSRoughness, DesPrFlow, MilleFlow, NumOntForceGates, 
     &    OntForceEOPLevel, OntForceFlow, LouisOntFlowDiff, 
     &    Qdpmi, Qstfran, Qstmau, QRich, RJetty1, RVarennes, 
     &    RSorel, RPierre, R3Rivers, RBatiscan, Tidal, 
     &    Other05Q, Other06n, Other06Q,
     &    LakeOntarioFlow, OntBOP, OntEOP, OntMeanLev)
C
c This routine is called by "CONTROL.FOR"
C
      IMPLICIT NONE
C
C
C     Declarations regarding function arguments:
C
      INTEGER       SDayOfMon, EDayOfMon, QMon, Monn, Year, DaysElapsed
C                   (Input) Starting and ending date for the Ontario/St.
C                   Lawrence routng interval, the quarter-month number
C                   containing the end day of the routing period, the
C                   month and year containing the end day, and the 
C                   of days elapsed since the start of the run 
C                   (potentially rolled over if very long simulation).
C
      LOGICAL       NewWeek, NewQMon, NewMon
C                   (Input) .TRUE. if the current Lake Ontario/St. 
C                   Lawrence routing period begins a new week, new 
C                   quarter-month, or new month (compared to the 
C                   previous Lake Ontario routing period).
C
      INTEGER       SecsInRtPeriod
C                   (Input) The lenght of the Lake Ontario routing 
C                   interval, in seconds.
C
      DOUBLE PRECISION NiagaraConverted, FlowFromWelland
C                   (Input) The flow coming from the Niagara River and 
C                   the Welland Canal, adjusted to the Lake Ontario 
C                   routing interval.
C
      DOUBLE PRECISION SupplyON
C                   (Input) The Lake Ontario NBS and miscellanous 
C                   supplies and diversions, adjusted to the Lake 
C                   Ontario routing interval.
C
      DOUBLE PRECISION ONNTS 
C                   (Input) The sum of NiagaraConverted, SupplyON, and
C                   FlowFromWelland.  This includes all Lake Ontario
C                   water balance terms, except outflow.  
C
      DOUBLE PRECISION OntRetard
C                   (Input) Total Lake Ontario outflow retardation,
C                   in cubic meters per second, adjusted to the 
C                   Lake Ontario routing interval.
C
      DOUBLE PRECISION OntDeviation, NumOntForceGates, OntForceEOPLevel,
     &              StLouisRetard, StLouisFlow, LSRoughness,
     &              OntForceFlow, LouisOntFlowDiff,
     &              DesPrFlow, MilleFlow,
     &              Qdpmi, Qstfran, Qstmau, QRich,
     &              RJetty1, RVarennes, RSorel, RPierre,
     &              R3Rivers, RBatiscan, Tidal, Other05Q,
     &              Other06n, Other06Q
C                   (Input) more inputs; likely to change
C
      DOUBLE PRECISION OntBOP
C                   (Input) Beginning-Of-Period level for Lake Ontario
C
      DOUBLE PRECISION OntEOP
C                   (Input/Output) End-Of-Period level for Lake Ontario
C
      DOUBLE PRECISION OntMeanLev
C                   (Output) Average Lake Ontario level for period
C
      DOUBLE PRECISION LakeOntarioFlow
C                   (Output) The Lake Ontario outflow for the routing
C                   period.
C
C     Declarations regarding local variables or functions:
C
C      INTEGER       I, J
      INTEGER       I
C                   Scratch
C
C                   
      CHARACTER*(*) Prefix
C                   Used at start of lines written to detailed.log, in
C                   order to help identify Ontario/St. Lawrence messages.
C
      DOUBLE PRECISION GetArea
C                   Function that returns surface area for a lake
C
      CHARACTER*16  MakeIntervalName
C                   Returns a string with the interval name.  For output.
C 
C     Apparently unused line commented by mmm 13Mar09 to hush compiler
C      DOUBLE PRECISION Round
C                   Returns value rounded to a certain accuracy,
C                   according to the old Lake Survey method.
C
C
C     Declarations regarding COMMON block variables:
C
      INTEGER       RoutingInterval(3), RoutingIncrement(3)
C                   Routing interval to use on Superior, middle lakes,
C                   and Ontario (1=daily, 2=weekly, 3=qmonthly, 
C                   4=monthly), and the number of increments within each 
C                   interval.  The Lake Ontario element for each of
C                   these arrays is the third one.
C
      INTEGER       UStLawrSolution
C                   Indicates how to determine Lake Ontario outflows. 
C                   1 means use some outflow equation, 2 means use
C                   plan 1958D. 
C
      INTEGER       BOMRound, FlowRound, LevelRound
C                   Will these be used? 
C
      INTEGER       MidLakeSolution, StMarysSolution, P77ACompat
C                   Not used; just in same COMMON as UStLawrSolution
C
      DOUBLE PRECISION ONK,ONYM,ONA,ONB,ONWT,ONC
C                   Discharge equation coefficients for Lake Ontario.  May
C                   not apply at all.  We can eliminate them if never used.
C                   See descriptions in CHECKPAR and OUTFLOW
C
      CHARACTER*80  OutputHeader
C                   Used to identify columns in ontario output file.
C
      INTEGER       TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &              MissingVal, SecsInUniformMon

      DOUBLE PRECISION MPFT, HALF, ONE, ZERO
C                   These are some handy constants.  
C                   See descriptions in Block Data General.  
C
      CHARACTER*3   MonthAbbrev(12)
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      INTEGER       GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &              SupVerbosity, MLVerbosity, OntVerbosity
      CHARACTER     CTMP*255, CWRK*80, LLT*1
C                   See descriptions in Block Data General.  
C
C
      COMMON /MISINT/   MidLakeSolution, BOMRound, FlowRound, 
     &                  LevelRound, RoutingInterval, RoutingIncrement,
     &                  StMarysSolution, UStLawrSolution, P77ACompat
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /ONDAT/    ONK,ONYM,ONA,ONB,ONWT,ONC
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
      COMMON /MONAME/   MonthAbbrev
      COMMON /ONTOUT/   OutputHeader
C
C
      PARAMETER (Prefix=' O/SL ==> ')
C
C
C
C     Initialize BOP levels as the EOP levels from previous period
      OntBOP = OntEOP
C
      IF( OntVerbosity .GE. 5 ) THEN  
C
         WRITE(CTMP,4105) Prefix, LLT
         CALL LogIt(LOGD,0,0)
C
C        Put in I the month number for month that period started 
         I = Monn
         IF( EDayOfMon .LT. SDayOfMon ) THEN 
            I = I - 1
            IF( I .EQ. 0 ) I =12
         ENDIF
         WRITE(CTMP,4110)Prefix, CWRK(1:I), EDayOfMon,MonthAbbrev(Monn),
     &      Year, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4115) Prefix, DaysElapsed, QMon, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4120) Prefix, SDayOfMon, MonthAbbrev(I), 
     &      SecsInRtPeriod, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4125) Prefix, NewWeek, NewQMon, NewMon, LLT
         CALL LogIt(LOGD,0,0)
         CWRK = MakeIntervalName( RoutingInterval(3), .FALSE., I)
C         WRITE(CTMP,4127) Prefix, CWRK(1:I), LLT
C         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4130) Prefix, NiagaraConverted, FlowFromWelland, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4135) Prefix, SupplyON, ONNTS, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4140) Prefix, OntRetard, StLouisRetard, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4145) Prefix, OntDeviation, NumOntForceGates, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4146) Prefix, OntForceEOPLevel, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4147) Prefix, OntForceFlow, LouisOntFlowDiff, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4150) Prefix, DesPrFlow, MilleFlow, LSRoughness, LLT
         CALL LogIt(LOGD,0,0)

         WRITE(CTMP,4151) Prefix, 'DPMI        ', Qdpmi, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4151) Prefix, 'St. Francois', Qstfran, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4151) Prefix, 'St. Maurice ', Qstmau, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4151) Prefix, 'Richelieu   ', Qrich, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4152) Prefix, 'Jetty 1     ', RJetty1, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4152) Prefix, 'Varennes    ', RVarennes, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4152) Prefix, 'Sorel       ', RSorel, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4152) Prefix, 'LacSt.Pierre', RPierre, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4152) Prefix, '3 Rivers    ', R3Rivers, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4152) Prefix, 'Batiscan    ', RBatiscan, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,4153) Prefix, 'Tidal Wave  ', Tidal, LLT
         CALL LogIt(LOGD,0,0)
 4151 FORMAT(A, A, '  flow:', F12.6,' m3s',A1)
 4152 FORMAT(A, A, ' Roughness:', F12.6,A1)
 4153 FORMAT(A, A, '     Value:', F12.6,A1)
         WRITE(CTMP,4155) Prefix, OntBOP, GetArea(5,OntBOP), LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
c
c**********************************************************
c
c     Calling the various Lake Ontario Plans
c
      if(UStLawrSolution .EQ. 2)then
c
c Calling Plan 1958-D
c
          call Plan58d(ONNTS, OntDeviation, LouisOntFlowDiff, 
     &        StLouisRetard, DesPrFlow, MilleFlow, OntForceFlow,
     &        OntMeanLev, LakeOntarioFlow, Qdpmi, Qstfran, Qstmau,
     &        QRich, RJetty1, RVarennes, RSorel, RPierre, R3Rivers, 
     &        RBatiscan, Tidal)
      elseif(UStLawrSolution .EQ. 5)then
c
c Calling Preproject Conditions
c
          call Pre_project(Onnts, LouisOntFlowDiff, StLouisRetard,
     &        DesPrFlow, MilleFlow, OntMeanLev, LakeOntarioFlow,
     &        Qdpmi, Qstfran, Qstmau, QRich, RJetty1, RVarennes, 
     &        RSorel, RPierre, R3Rivers, RBatiscan, Tidal)
      endif
c
c**********************************************************
C
C -------------- End of Week averaging and logging:
C      WRITE(CWRK, 4210) EDayOfMon, Monn, Year, Round(OntBOP, -3), 
C     &   Round(OntMeanLev, -3), Round(OntEOP, -3), 
C     &   IDNINT(Round(LakeOntarioFlow,0)), LLT  
C     Copy line to CTMP and send it to Ontario output file
C      CTMP = CWRK
C      CALL LogIt(LOGO,0,0)
C      IF( OntVerbosity .GE. 4 ) THEN  
C        Write header line first
C         WRITE(CTMP,4000) Prefix, OutputHeader, LLT
C         CALL LogIt(LOGD,0,0)
C        Send line to detailed.log
C         CTMP = Prefix//CWRK
C         CALL LogIt(LOGD,0,0)
C      ENDIF
C
C
C
C     Send the answers off to the accumulator for derived averages
      CALL TallyMeans(5, 3, LakeOntarioFlow, OntMeanLev)
C
C
C 4000 FORMAT(A,A, A1)
      return
 4105 FORMAT(A, A1)
 4110 FORMAT(A, 'Lake Ontario/St. Lawrence ', A, 
     &   ' routing for period ending ', I2, A3, I5, '.', A1)
 4115 FORMAT(A, 'Period end date occurs on day ', I9, 
     &   ' of simulation, within quarter-month ', I1, '.', A1)
 4120 FORMAT(A, 'Period started on ', I2, A3, ', with duration of ',
     &   I7, ' seconds.', A1)
 4125 FORMAT(A, 'Since the previous routing period, do we ',
     &   'start a new week (', L1, '), qmonth (', L1, '), month (', 
     &   L1, ')?', A1)
c4127 FORMAT(A, 'Using ', A, ' routing interval.', A1)
 4130 FORMAT(A, 'Niagara River:', F12.6, ' m3s,   Welland Canal:', 
     &   F12.6, ' m3s', A1)
 4135 FORMAT(A, 'NBS, et al:', F13.6, ' m3s,   Net Total Supply:', 
     &   F13.6, ' m3s', A1)
 4140 FORMAT(A, 'Ontario outflow retardation:', F12.6, 
     &   ' m3s,   St. Louis retardation:', F12.6, ' m3s', A1)
 4145 FORMAT(A, 'Ontario Deviation:', F12.6, 
     &   ' m3s,   Forced gate setting:', F11.4, A1)
 4146 FORMAT(A, 'Ontario EOP level to force (BOP for next period):', 
     &   F12.6, ' m', A1)
 4147 FORMAT(A, 'Ontario forced outflow:', F12.6, 
     &   ' m3s,   St Louis/Ontario Flow Difference:', F12.6, A1)
 4150 FORMAT(A, 'Des Pr flow:', F12.6, ' m3s  Mille flow:', F12.6, 
     &   ' m3s,   LS Roughness:', F13.6, A1)
C     Apparently unused line commented by mmm 13Mar09 to hush compiler
C4154 FORMAT(A, A, ' Roughness:', F12.6, '  flow:', F12.6,' m3s',A1)
 4155 FORMAT(A, 'Lake Ontario Beginning-Of-Period level:', F12.6, 
     &   '  corresponding area:', F13.1, ' m2', A1)
c4205 FORMAT(A, '  Period End     BOP    Mean     EOP Outflow', A1)
c4210 FORMAT( 1X, 2(I3), I5, 3F8.3, I8, A1)
      end
