C-----------------------------------------------------------------------
      SUBROUTINE CHECKPAR()
C
      IMPLICIT NONE
C
C     Declarations regarding local variables or functions:
C
      CHARACTER    TmpStr*256
      CHARACTER*9  DayOfWeekName(7)
      CHARACTER    Unspecified*12
C     CHARACTER    SolveString*(*)
      DOUBLE PRECISION TMPA, TMPZ
      DOUBLE PRECISION Areas(8)
      INTEGER      I, J, K, L, M, F(5), U(5), K1, idd
      LOGICAL      LOGUSE(8)
C                  These are scratch variables.
C
      CHARACTER    CHNKST*16
C                  This is the calling subroutine name passed to ErrMess() in
C                  case of error messages.
C
      INTEGER      ENDSTR
C                  Function that returns the end position of the current word.
C 
      INTEGER      LastKC
C                  Function that returns the position of last 'keyboard char'
C 
      INTEGER      FirstKC
C                  Function that returns the position of first 'keyboard char'
C 
      INTEGER      TSOutputFileName
C                  Formulates filename and path for output files.
C
      INTEGER      DateSequence
C                  Returns the number of days since 1Jan0000
C
      DOUBLE PRECISION GETAREA
C                  Function that returns surface area for a lake
C
      CHARACTER*13 NameOfLoop(3)
C 
C
C
C     Declarations regarding COMMON block variables:
C
      INTEGER      P58DStartDay, P58DStartMonth, P58DStartYear
C                  Indicate the date of the first Saturday of the simulation
C
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                  See description in Block Data Genral
C
      INTEGER      NLN(16), LkToSlv2LN(8), NRN(16), LkSlv2RtLoopI(8)
C                  See descriptions in Block Data Genral
C
      CHARACTER*21 LAKENAME(16), RiverName(16)
C                  See description in Block Data Genral
C
      CHARACTER*256 FileNames(100)
C                  See description in Block Data Genral
C
      CHARACTER*256 RepeatLines(80)
C                  See description in Block Data Genral
C
      CHARACTER    CTMP*255, CWRK*80
C                  These are just buffers used throughout program.
C
      CHARACTER    LLT*1
C                  This is the character that marks the end of a line sent
C                  to the LogIt subroutine.
C
      CHARACTER*3  MonthAbbrev(12)
C                  See description in Block Data Genral
C
      CHARACTER*80 RCTITLE(3), MLHEADER
C                  These are descriptive lines for use in output.
C
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C                  See descriptions in Block Data Genral
C
      INTEGER      MidLakeSolution, StMarysSolution, UStLawrSolution 
      INTEGER      BOMRound, FlowRound, LevelRound, P77ACompat
C                  See descriptions in Block Data Genral.
C
      INTEGER      RoutingInterval(3), RoutingIncrement(3)
C                  Routing interval to use on Superior, middle lakes,
C                  and Ontario (1=daily, 2=weekly, 3=qmonthly, 
C                  4=monthly), and the number of increments within each 
C                  interval.
C
      INTEGER      PeriodEnd(6,3), Day(3), Q(4,3)
      DOUBLE PRECISION dT(3)
      LOGICAL      NumIncsVary(3) 
C                  See descriptions in RUNIT()
C
      INTEGER      MonthLength
      CHARACTER*8  MonthLengthName(2)
C                   See description in Block Data Genral
C
      CHARACTER*8  RCRUNTYP, IntervalName(4)
C                  See description in Block Data Genral
C
      CHARACTER*9  MLSolutionName(2)
C                  Strings set to 'interative' and 'matrix', used for output.
C
      CHARACTER*16 SMSolutionName(24), USLSolutionName(5)
      LOGICAL      WriteTSOutput(4,2,8), WriteErieOutflow
C                  See description in Block Data Genral
C
      LOGICAL      LakesToSolve(8)
C                  See description in Block Data Genral
C
      INTEGER      StartYear, StartMonth, StartDay, 
     &             EndYear, EndMonth, EndDay
C                  See descriptions in Block Data Genral
C
      INTEGER      DaysSinceSaturday
C                  The start date's offset from the previous Saturday 
C
      INTEGER      IntervalNameLen(4)
C                  See description in Block Data Genral
C
      DOUBLE PRECISION SGX(17), SGY(17), SGZ(17)
      INTEGER       MaxSGSetting
      INTEGER      NUMAVSZ(8)
C                  See descriptions in Block Data Genral
C
      DOUBLE PRECISION AVSZ(8,5,2)
C                  See descriptions in Block Data Genral
C
      DOUBLE PRECISION SPSTLV, MHSTLV, SCSTLV, ERSTLV, CGSTLV, ONSTLV
C                  See descriptions in Block Data Genral
C
      DOUBLE PRECISION SPBOMLev, MHBOMLev, SCBOMLev, ERBOMLev
C                  The beginning-of-month lake levels used by Plan 1977A
C
      REAL         PrevSupGates 
C                  Initial gates open at Comp Works
C 
      INTEGER      PrevStMarysFlow
C                  Initial Lake Superior outflow (average for month 0)
C 
      REAL         SupGatesOpen
C                  The number of Superior gates open in current routing 
C                  period.
C
      DOUBLE PRECISION SPK,SPYM,SPA,SPB,SPWT,SPC
      DOUBLE PRECISION MHK,MHYM,MHA,MHB,MHWT,MHC
      DOUBLE PRECISION SCK,SCYM,SCA,SCB,SCWT,SCC
      DOUBLE PRECISION ERK,ERYM,ERA,ERB,ERWT,ERC
      DOUBLE PRECISION ONK,ONYM,ONA,ONB,ONWT,ONC
C                  These are the stage-discharge equation K, ym, a, b, zc,
C                  and wt, as described by GLERL's Tech Memo ERL-GLERL-xx,
C                  and in the function OUTFLOW.  Used for when the lake
C                  outflow is determined from an outflow equation, instead
C                  of a regulation plan.
C
      DOUBLE PRECISION FWRK
      DIMENSION FWRK(75,36)
C                  This is a common block used for 64-bit calculations.  It is just
C                  used as scratch until execution reaches RUNIT() in CONTROL.FOR
C
      INTEGER      TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &             MissingVal, SecsInUniformMon
      DOUBLE PRECISION MPFT, HALF, ONE, ZERO
C                  See descriptions in Block Data General.  
C
      logical   EElimDev,FForceFlow,SSetFlow,Usedatabase,SiDatabase
      INTEGER   CritDepths, P58DTotalSup, MaxILimit, MinPLimit,
     &          P58DForecastDay, P58DRunType,
     &          P58DForecastYear, P58DForecastMon,  
     &          P58DWeightedSupply, P58DRuleCurve, P58DAppLimit,
     &          SupplyIndicator(21)
      Integer   P58DAdjustment, P58DrcAdjustment, P58DActualFlow, 
     &          P58DTotalDev, P58DLSLflow
      DOUBLE PRECISION P58DActualLevel, P58DPlanLevel
      DOUBLE PRECISION PreprojectLevel, onstlv1
      CHARACTER*2  P58DAppLimitSym
      CHARACTER*256 P58DSimulateFile, P58DSIDataFile
      LOGICAL   ElimDev, SetFlow, PreProjConditions, ForceFlow
      logical   P58DUseDatabase, P58DSiDatabase, ComputeMontreal
      INTEGER   CritBOMRound, CritFlowRound, CritLevelRound, 
     &                  Q1887Method, CritBMethod
      DOUBLE PRECISION CritCThresholdEl
C               See descriptions in Block Data General.  
C
C               These variables included just beacuse they are
C               in same COMMON as the rounding variables that we need.
      INTEGER      P77AWinterMax, P77AMaxChange, P77AMinimumFlow(12),
     &             P77ABaseFlow(12), P77AMonthsToAve(12),
     &             P77ASup05NBS(12), P77ASup50NBS(12), P77ASup95NBS(12), 
     &             P77ALLOgokiDiv(12),P77AStMarysRetrd(12), 
     &             P77AMH05NBS(12), P77AMH50NBS(12), P77AMH95NBS(12), 
     &             P77AChicagoDiv(12), P77AStClairRetrd(12), 
     &             P77AStC50NBS(12), P77ADetroitRetrd(12),
     &             P77AErie50NBS(12), P77AWellandFlow(12), 
     &             P77ANiagaraRetrd(12),
     &             P77AMidLakesInc, P77AMidLakesTol, 
     &             P77AMonthRound, P77AFlowRound, P77ALevelRound
      DOUBLE PRECISION P77ASupTargetLev(12), P77ASupRangeFact(12),
     &             P77AMHTargetLev(12), P77AMHRangeFact(12), P77AAValue,
     &             P77ACGIP, P77ALowTrigger, P77AHighTrigger, 
     &             P77ASPK, P77ASPYM, P77ASPA, P77ASPB,P77ASPWT,P77ASPC, 
     &             P77AMHK, P77AMHYM, P77AMHA, P77AMHB,P77AMHWT,P77AMHC, 
     &             P77ASCK, P77ASCYM, P77ASCA, P77ASCB,P77ASCWT,P77ASCC, 
     &             P77AERK, P77AERYM, P77AERA, P77AERB,P77AERWT,P77AERC
      LOGICAL      P77AUniformMonth
      INTEGER      COUNTER_JUNEFLOW, Last_JuneFlow
      INTEGER      IndexGatesSide
      DOUBLE PRECISION   ForcedSupGates(12)
      INTEGER      ForcedSideChan(12)
C
C               These variables included just beacuse they are
C               in same COMMON as LStoSWP and SupCritRnd
      INTEGER   ForecastMonths, LStoSWP, SupCritRnd, MonthsChecked
      INTEGER   SPOPNBS06(12,6), SPOPNBS07(12,6), OPLLO(12), 
     &          MHOPNBS06(12,6), MHOPNBS07(12,6), OPChicago, 
     &          SCRRet(12), SCOPNBS06(12,6), SCOPNBS07(12,6), 
     &          DetRet(12), EROPNBS06(12,6), EROPNBS07(12,6),
     &          OPWelland(12), NiaRet(12),
     &          SCRRet_P2012(12), SCOPNBS08(12,6), SPOPNBS08(12,6), 
     &          DetRet_P2012(12), EROPNBS08(12,6), MHOPNBS08(12,6),
     &          OPWelland_P2012(12), NiaRet_P2012(12)  
      LOGICAL   BackCompat
C
      CHARACTER*256 CritKfile
C
C
      COMMON /CritK/    CritKfile
      common/initpass/  EElimdev,FForceFlow,SSetFlow,Usedatabase,
     &                  SiDatabase
      COMMON /SUPCOM/   SPOPNBS06, OPLLO, MHOPNBS06, OPChicago, SCRRet,
     &                  SCOPNBS06, DetRet, EROPNBS06, OPWelland, NiaRet,
     &                  SPOPNBS07, MHOPNBS07, SCOPNBS07, EROPNBS07,
     &                  SPOPNBS08, MHOPNBS08, SCOPNBS08, EROPNBS08,
     &                  OPWelland_P2012, SCRRet_P2012, DetRet_P2012,
     &                  NiaRet_P2012,
     &                  ForecastMonths,
     &                  LStoSWP, SupCritRnd, BackCompat, MonthsChecked
      COMMON /OUTSTF/   WriteTSOutput, WriteErieOutflow
      COMMON /FSPACE/   FWRK
      COMMON /MONAME/   MonthAbbrev
      COMMON /FSAT/     P58DStartDay, P58DStartMonth, P58DStartYear
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /LNAMEC/   LAKENAME, RiverName
      COMMON /LNAMEI/   NLN, LkToSlv2LN, NRN, LkSlv2RtLoopI
      COMMON /RCAVSZ/   NUMAVSZ, AVSZ
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /STAGE/    MHK,SCK,ERK,MHYM,SCYM,ERYM,MHA,SCA,ERA,
     &                  MHB,SCB,ERB,MHWT,SCWT,ERWT,MHC,SCC,ERC
      COMMON /SPDAT/    SPK,SPYM,SPA,SPB,SPWT,SPC,
     &                  SGX, SGY, SGZ, MaxSGSetting,
     &                  SupGatesOpen
      COMMON /ONDAT/    ONK,ONYM,ONA,ONB,ONWT,ONC
      COMMON /MISINT/   MidLakeSolution, BOMRound, FlowRound, 
     &                  LevelRound, RoutingInterval, RoutingIncrement,
     &                  StMarysSolution, UStLawrSolution, P77ACompat
      COMMON /INITLV/   SPSTLV, MHSTLV, SCSTLV, ERSTLV, CGSTLV, ONSTLV,
     &                  PrevStMarysFlow, PrevSupGates, 
     &                  SPBOMLev, MHBOMLev, SCBOMLev, ERBOMLev  
      COMMON /FNAMEC/   FileNames, RepeatLines
      COMMON /RCRUNC/   RCRUNTYP, IntervalName, MonthLengthName,
     &                  MLSolutionName, SMSolutionName, USLSolutionName
      COMMON /RCRUNI/   LakesToSolve, 
     &                  IntervalNameLen, MonthLength,
     &                  StartYear, StartMonth, StartDay,
     &                  EndYear, EndMonth, EndDay
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
      COMMON /CLDATE/   dT, Day, Q, PeriodEnd, NumIncsVary,
     &                  DaysSinceSaturday
      COMMON /INI58D/   CritDepths, P58DTotalSup, MaxILimit, MinPLimit,
     &                  P58DWeightedSupply, P58DRuleCurve, P58DAppLimit,
     &                  P58DRunType, P58DForecastYear, P58DForecastMon,
     &                  P58DForecastDay, SupplyIndicator, 
     &                  P58DPlanLevel, P58DActualLevel,
     &                  P58DSimulateFile, P58DSIDataFile, ElimDev, 
     &                  ForceFlow, SetFlow, PreProjConditions,
     &                  P58DUseDatabase, P58DSiDatabase
      COMMON /Ini58dc/  P58DAppLimitSym
      COMMON /INI58D1/  P58DAdjustment,P58DrcAdjustment,P58DActualFlow, 
     &                  P58DTotalDev,P58DLSLflow
      COMMON /INI58D2/  ComputeMontreal
      common /PreCond/  PreprojectLevel
      common /OntBeglv/ onstlv1
      COMMON /TITLES/   RCTITLE, MLHEADER
      COMMON /SUPCRT/   CritBOMRound, CritFlowRound, CritLevelRound,
     &                  Q1887Method, CritBMethod, CritCThresholdEl
      COMMON /P1977A/   P77ASPK, P77ASPYM, P77ASPA, P77ASPB, P77ASPWT, 
     &                  P77ASPC, P77AMHK, P77AMHYM, P77AMHA, P77AMHB, 
     &                  P77AMHWT, P77AMHC, P77ASCK, P77ASCYM, P77ASCA,
     &                  P77ASCB, P77ASCWT, P77ASCC, P77AERK, P77AERYM,
     &                  P77AERA, P77AERB, P77AERWT, P77AERC, P77ACGIP,
     &                  P77ASupTargetLev, P77AMHTargetLev, P77AAValue, 
     &                  P77ASupRangeFact, P77AMHRangeFact, 
     &                  P77ALowTrigger, P77AHighTrigger, P77AMaxChange,
     &                  P77AMonthsToAve, P77AMinimumFlow, P77AWinterMax,
     &                  P77ASup05NBS, P77ASup50NBS, P77ASup95NBS,
     &                  P77AMH05NBS, P77AMH50NBS, P77AMH95NBS,
     &                  P77AStC50NBS, P77AErie50NBS, P77AWellandFlow,
     &                  P77AChicagoDiv, P77ALLOgokiDiv, P77ABaseFlow,
     &                  P77AStMarysRetrd, P77AStClairRetrd, 
     &                  P77ADetroitRetrd, P77ANiagaraRetrd, 
     &                  P77AMidLakesInc, P77AUniformMonth, 
     &                  P77AMonthRound, P77AFlowRound, P77ALevelRound,
     &                  P77AMidLakesTol
      COMMON /STURGEON/ COUNTER_JUNEFLOW, Last_JuneFlow 
      COMMON /FORGATE/ ForcedSupGates     
      COMMON /FORSIDE/ ForcedSideChan 
      COMMON /INFORCE/ IndexGatesSide     
C
C
      DATA NameOfLoop /'Lake Superior','Middle Lakes ','Lake Ontario '/
C
      DATA DayOfWeekName / 'Tuesday', 'Wednesday', 'Thursday', 'Friday',
     &                     'Saturday', 'Sunday', 'Monday' /
C
      DATA AREAS /8.21D10, 1.174D11, 5.78D10, 5.96D10, 1.114D9, 2.57D10,
     &            1.896D10, 0.0D0/
C    These values used simply to test reasonability of user input
C
      CHNKST='CheckPar'//LLT
C
      Unspecified='Unspecified'//LLT
C
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvv   Set up Verbosity variables    vvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C     Use General Verbosity unless overridden with input parameter
      IF ( TSDVerbosity .EQ. MissingVal) TSDVerbosity = GenVerbosity
      IF ( ChkVerbosity .EQ. MissingVal) ChkVerbosity = GenVerbosity
      IF ( SupVerbosity .EQ. MissingVal) SupVerbosity = GenVerbosity
      IF ( MLVerbosity  .EQ. MissingVal) MLVerbosity  = GenVerbosity
      IF ( OntVerbosity .EQ. MissingVal) OntVerbosity = GenVerbosity
C
C
C                     <END OF Verbosity SECTION>  
C
C
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvv   Set up where CheckPar() sends messages    vvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C     Send messages to console if ChkVerbosity > 3, else just send to
C     detailed log.
      DO 6255, I=1,8
 6255    LOGUSE(I) = LOGCD(I)
      IF(ChkVerbosity .LE. 3)THEN
         DO 6256, I=1,8
 6256      LOGUSE(I) = LOGD(I)
      ENDIF
C
      WRITE(CTMP, '(A,A1)')
     &   ' Checking input for completeness and consistency . . .',LLT
      CALL LogIt(LOGUSE,2,0)
C
C                     <END OF CheckPar LOG SECTION>  
C
C
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvv  Run Titles  vvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C
      J = 0
      DO 2104, I = 1,3
 2104    IF( FirstKC(RCTitle(I), 80) .GT. 0) J = 1
C
      IF ( J .EQ. 0 ) THEN
         WRITE(CTMP,5324) LLT
         CALL LogIt(LOGCDE,1,0)
         WRITE(CTMP,5445) 
     &      ' No titles given for the run (Title 1, etc)', LLT
         CALL LogIt(LOGCDE,0,0)
         CALL SAVETOD
         CALL ADIOS
      ENDIF
C
      DO 2105, I = 1,3
         WRITE(CTMP,'(A,A1)') RCTitle(I)(1:79), LLT
 2105    CALL LogIt(LOGUSE,0,I/3)
C
C
C                     <END OF RUN TITLES SECTION>
C
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvv   figure out which lakes to solve   vvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C      LOGICAL      LakesToSolve(8)
C                  an 8 element logical array listing which lakes to solve for.
C                     1          Lake Superior
C                     2   Lakes Michigan-Huron
C                     3         Lake St. Clair
C                     4              Lake Erie
C                     5           Lake Ontario
C                     6      Lake St. Lawrence
C
C     Keep track of total lakes to solve in M
      M = 0
C
      IF( LakesToSolve(1) ) THEN
         M = M + 1
         WRITE(CTMP, 5322)LAKENAME(1)(1:NLN(1)),LLT
         CALL LogIt(LOGUSE,0,0)
         IF(StMarysSolution .EQ. MissingVal) THEN
            CALL ErrMess(1120,CHNKST,LOGCDE)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
         WRITE(CTMP, '(1X,A,A,A,A)') LAKENAME(1)(1:NLN(1)),
     &      ' outflows computed using ', 
     &      SMSolutionName(StMarysSolution),LLT
         CALL LogIt(LOGUSE,0,0)
      ENDIF
C
      IF( LakesToSolve(2) .OR. LakesToSolve(3) .OR. LakesToSolve(4))THEN
C        Use J as temporary flag inside 2024 loop
         J = 0
C        If solving for any of the middle lakes, must solve for all 3
         DO 2024, I=2,4
            IF( LakesToSolve(I) )THEN
               WRITE(CTMP, 5322)
     &            LAKENAME(LkToSlv2LN(I))(1:NLN(LkToSlv2LN(I))),LLT
               CALL LogIt(LOGUSE,0,0)
            ELSE
C              Hmmm.  Trying to *not* solve for one?
               IF( J .EQ. 0) THEN
                  CALL ErrMess(1113,CHNKST,LOGCDE)
                  J = 1
               ENDIF
               WRITE(CTMP,'(3X,A,A1)')
     &            LAKENAME(LkToSlv2LN(I))(1:NLN(LkToSlv2LN(I))),LLT
               CALL LogIt(LOGCDE,0,0)
            ENDIF
 2024       M = M + 1
         IF(J .EQ. 1)THEN
            CALL SAVETOD
            CALL ADIOS
         ENDIF
C
         IF( FirstKC(MLHEADER, 80) .EQ. 0) THEN
            WRITE(CTMP,5323) 'Routing Title', LLT
            CALL LogIt(LOGCDE,0,0)
         ENDIF
C
C        See if a routing technique is specified:
         IF( MidLakeSolution .EQ. MissingVal) THEN
            WRITE(CTMP,5324) LLT
            CALL LogIt(LOGCDE,1,0)
            WRITE(CTMP,5445) ' No solution technique specified '//
     &         'for middle lakes routing (Mid-Lake Solution)', LLT
            CALL LogIt(LOGCDE,0,0)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
C
C        Use F array for scratch, with length of MLSolutionName for each 
C        value of MidlakesRouting:
         F(1) = 9
         F(2) = 6
         WRITE(CTMP,5339)
     &      MLSolutionName(MidLakeSolution)(1:F(MidLakeSolution)), LLT
         CALL LogIt(LOGD,1,1)
      ENDIF
 5339 FORMAT(' Using ',A,' procedure for middle lake routing.',A1)
C  
      IF( LakesToSolve(5) ) THEN
         M = M + 1
         WRITE(CTMP, 5322)LAKENAME(7)(1:NLN(7)),LLT
         CALL LogIt(LOGUSE,0,0)
         IF(UStLawrSolution .EQ. MissingVal) THEN
            CALL ErrMess(1131,CHNKST,LOGCDE)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
         WRITE(CTMP, '(1X,A,A,A,A)') LAKENAME(7)(1:NLN(7)),
     &      ' outflows computed using ', 
     &      USLSolutionName(UStLawrSolution),LLT
         CALL LogIt(LOGUSE,0,0)
      ENDIF
 5322 FORMAT(' Will attempt to solve for ',A,'.',A1)
C
      IF (M .EQ. 0 ) THEN
         WRITE(CTMP,5324) LLT
         CALL LogIt(LOGCDE,1,0)
         WRITE(CTMP,5445) ' Failed to specify which lakes to route'//
     &      ' through (Lakes To Solve For)', LLT
         CALL LogIt(LOGCDE,0,0)
         CALL SAVETOD
         CALL ADIOS
      ENDIF
C
C     Write an empty line
      CTMP(1:1) = LLT
      CALL LogIt(LOGUSE,0,0)
C
C                     <END OF LAKES TO SOLVE SECTION>  
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvv   surface area checks   vvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C                  Published areas are:
C                     1          Lake Superior:    82100 km2 
C                     2   Lakes Michigan-Huron:   117400 km2 
C                     3          Lake Michigan:    57800 km2 
C                     4             Lake Huron:    59600 km2 
C                     5         Lake St. Clair:     1114 km2 
C                     6              Lake Erie:    25700 km2 
C                     7           Lake Ontario:    18960 km2 
C                     8         not used (yet):        0 km2
C                  These figures are from the 1977 "Coordinated Great
C                  Lakes Physical Data" book, page 8.  This order must
C                  be the same as the first eight lakes in LAKENAME and NLN.
C
C 
C     Check the surface area information: 
C
C     Loop through each lake and see if more than A-Z point 
      DO 2099, M=1,8
C        
         L = LkToSlv2LN(M)
C        Do nothing for the not-implemented lakes
         IF ( L .LE. 0 ) GOTO 2099
C 
         IF ( .NOT. LakesToSolve(M) ) THEN
C           Skip to next lake if no area and not solving
            IF ( NUMAVSZ(L) .EQ. MissingVal ) GOTO 2099
C           If fell through, then gave area for lake we are NOT modeling
            CALL ErrMess(1158,CHNKST,LOGCDE)
            WRITE(CTMP,5331) NUMAVSZ(L), LAKENAME(L)(1:NLN(L)),LLT
            CALL LogIt(LOGCE,1,0)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
C
C        At this point, LakesToSolve(M) is .TRUE.
         IF ( NUMAVSZ(L)  .EQ. MissingVal) THEN
C           Failed to specify area for lake we are modeling
            CALL ErrMess(1159,CHNKST,LOGCDE)
            WRITE(CTMP,5328) LAKENAME(L)(1:NLN(L)),LLT
            CALL LogIt(LOGCDE,1,0)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
C  
         IF ( NUMAVSZ(L) .EQ. 1 ) THEN
C           Constant area
            WRITE(CTMP,5334)AVSZ(L,1,2), LAKENAME(L)(1:NLN(L)),LLT
            CALL LogIt(LOGUSE,0,0)
C    
         ELSE
C           At least two points specified:
            WRITE(CTMP,5331) NUMAVSZ(L), LAKENAME(L)(1:NLN(L)),LLT
            CALL LogIt(LOGUSE,1,0)
            WRITE(CTMP, 5332) LLT
            CALL LogIt(LOGUSE,0,0)
C
            IF(ChkVerbosity .GE. 9)THEN
C           Print out values in memory before sorting
               WRITE(CTMP, 5330) '   Before',LLT
               CALL LogIt(LOGD,0,0)
               DO 2091, J=1,5
                  WRITE(CTMP,5333)J, AVSZ(L,J,1), AVSZ(L,J,2),LLT
 2091             CALL LogIt(LOGD,0,0)
               WRITE(CTMP, 5330) '   After',LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
C           Sort array   
C           Start at highest index and work to bottom (lowest Z at index 1)
            DO 2094, J = NUMAVSZ(L), 2, -1
C              Choose arbitrary initial index of great elevation
               K = 1
C              Loop through other elevations (after arbitrary initial one)
               DO 2093, I = 2, J
C                 If higher than one indexed by K, assign this index to K
 2093             IF ( AVSZ(L,I,1) .GT. AVSZ(L,K,1) ) K = I
C              Save current values indexed by J into TMPA and TMPZ
               TMPZ = AVSZ(L,J,1)
               TMPA = AVSZ(L,J,2)
C              Put values for highest remaining elev into J's data pair
               AVSZ(L,J,1) = AVSZ(L,K,1)
               AVSZ(L,J,2) = AVSZ(L,K,2)
C              Put values that were in J into K
               AVSZ(L,K,1) = TMPZ
               AVSZ(L,K,2) = TMPA
C
 2094          CONTINUE
C
C           OK, now should be sorted in increasing elevation
C           See if area are also increasing, and at least 1mm between elevs:
            J = 0
            DO 2096, I=2,NUMAVSZ(L)
               IF(AVSZ(L,I-1,2) .GT. AVSZ(L,I,2))THEN
                   J = 1049
                   GOTO 2097
               ENDIF
               IF( (AVSZ(L,I,1) - AVSZ(L,I-1,1))*1000 .LT. .501) THEN
                   J = 1050
                   GOTO 2097
               ENDIF
 2096          CONTINUE
 2097       IF( J .NE. 0 ) THEN
               CALL ErrMess(J,CHNKST,LOGCDE)
               WRITE(CTMP,
     &            '(1X,A,A,F7.3,A1,F17.3,A,F7.3,A1,F17.3,A2,A1)')
     &            LAKENAME(L)(1:NLN(L)), ' points (',AVSZ(L,I-1,1),',',
     &            AVSZ(L,I-1,2), ') and (', AVSZ(L,I,1), ',', 
     &            AVSZ(L,I,2), ').', LLT
                  CALL LogIt(LOGCDE,1,0)
                  CALL SAVETOD
                  CALL ADIOS
            ENDIF
            IF(ChkVerbosity .GE. 9)THEN
               I = 5
            ELSE
               I = NUMAVSZ(L)
            ENDIF
            DO 2098, J=1,I
               WRITE(CTMP,5333)J, AVSZ(L,J,1), AVSZ(L,J,2), LLT
 2098          CALL LogIt(LOGUSE,0,0)
         ENDIF
C
C
C        Test to see if in ballpark of official area
         I = NUMAVSZ(L)
         IF ( AVSZ(L,1,2) .LT. 0.49D0*Areas(L) .OR. 
     &        AVSZ(L,I,2) .GT. 2.1D0*Areas(L) ) THEN
            CALL IssueWarning(1185,LogCDE,L,CHNKST)
            IF ( AVSZ(L,1,2) .LT. 0.49D0*Areas(L)) THEN  
               TMPA = AVSZ(L,1,2)
            ELSE 
               TMPA = AVSZ(L,I,2)
            ENDIF 
            WRITE(CTMP, 5325)  TMPA, Areas(L), LLT
            CALL LogIt(LOGCDE,0,0)
         ENDIF
C
 2099    CONTINUE
C
 5328 FORMAT(' No area read for ',A,'!',A1)
 5330 FORMAT(1X,A,' sorting . . .',A1)
 5331 FORMAT(' Read ',I1,' surface area vs. elevation data points for ',
     &   A,' . . .',A1)
 5332 FORMAT('  Point #',5X,'Elev m IGLD-85',23X,'Area in m2',A1)
 5333 FORMAT(I5,16X,F7.3,1X,F32.3,A1)
 5334 FORMAT(' Using constant area of ',F18.3,' sq. meters for ',A,A1)
 5325 FORMAT('   - specified lake area of ', D18.12, 
     &   ' differs greatly from official value of ', D12.6, ' m2.',A1)
C
C                     <END OF SURFACE AREA SECTION>  
C

C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvv   Global Rounding Parameters      vvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C 
      CWRK(1:1) = LLT
C
      IF( FlowRound .EQ. MissingVal) 
     &   WRITE(CWRK, 5329) 'intermediate flows', 'Flow Round', LLT
      IF( LevelRound .EQ. MissingVal) 
     &   WRITE(CWRK, 5329) 'intermediate levels', 'Level Round', LLT
      IF( BOMRound .EQ. MissingVal) 
     &   WRITE(CWRK, 5329) 'beginning-of-month levels', 'BOM Round', LLT
       
      IF (CWRK(1:1) .NE. LLT ) THEN
         WRITE(CTMP,5324) LLT
         CALL LogIt(LOGCDE,1,0)
         CTMP = CWRK
         CALL LogIt(LOGCDE,0,0)
         CALL SAVETOD
         CALL ADIOS
      ENDIF
C
      WRITE(CTMP,5445) ' Global rounding settings', LLT
      CALL LogIt(LOGUSE,1,0)
      WRITE(CTMP,5327) '           Intermediate flows: ', FlowRound, LLT
      CALL LogIt(LOGUSE,0,0)
      WRITE(CTMP,5327) '          Intermediate levels: ', LevelRound,LLT
      CALL LogIt(LOGUSE,0,0)
      WRITE(CTMP,5327) '    Beginning-of-month levels: ', BOMRound, LLT
      CALL LogIt(LOGUSE,0,1)

 5329 FORMAT(' Failed to specify how to round ',A,'! (',A,')',A1)
 5327 FORMAT(A,I3, A1)
 5445 FORMAT(A,A1)
C
C                     <END OF GLOBAL ROUNDING SECTION>
C
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvv   Uniform or Actual Month Length  vvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C 
      J = 0
      DO 2175, I = 1,3
 2175    IF( RoutingInterval(I) .GE. 3 ) J = I
C
      IF( J .EQ. 0 ) THEN
C        Only daily or weekly routing; month length irrelevant
         IF( MonthLength .EQ. MissingVal) THEN
            WRITE(CTMP,2176) 'Month Length', LLT
            CALL LogIt(LOGCDE,1,0)
            WRITE(CTMP,2177) 'Specified ', MonthLengthName(MonthLength)
     &          (1:ENDSTR(MonthLengthName(MonthLength),8)), ' month '//
     &          'length, but there is no monthly or qmonthly routing!', 
     &           LLT
            CALL LogIt(LOGCDE,1,0)
            CALL SAVETOD
            CALL ADIOS
C        ELSE
C           could log a note that month length not needed, and not given
         ENDIF
      ELSE
C        Have some monthly or qmonthly routing; need month length 
         IF( MonthLength .EQ. MissingVal) THEN
            WRITE(CTMP,5445) ' Month length not specified! '//
     &         'Must choose actual or uniform.', LLT
            CALL LogIt(LOGCDE,1,0)
            CALL SAVETOD
            CALL ADIOS
         ELSE
            WRITE(CTMP,2177)
     &         ' Change in storage calculations based on ',
     &         MonthLengthName(MonthLength)
     &         (1:ENDSTR(MonthLengthName(MonthLength),8)),
     &         ' month/quarter-month length.', LLT
            CALL LogIt(LOGUSE,0,1)
         ENDIF
      ENDIF
C
C   
C
C                     <END OF MONTH LENGTH SECTION>
C
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvv   routing intervals and increments   vvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C                          LakesToSolve:
C                     1          Lake Superior
C                     2   Lakes Michigan-Huron
C                     3         Lake St. Clair
C                     4              Lake Erie
C                     5           Lake Ontario
C                     6      Lake St. Lawrence
C
C     Put index to LakesToSolve for into soemthing corresponding to loop
      F(1) = 1
      F(2) = 2
      F(3) = 5
C
C     Pre-set the Ontario routing increments for Plan 1958D to 1 (because
C     outflow is constant) in case not specified explicitly in input.  
      IF( LakesToSolve(5) .AND. RoutingIncrement(3) .EQ. MissingVal)
     &   RoutingIncrement(3) = 1
C     Look at the increments used within each time step:
C
      DO 2268, I = 1,3
C
C
         IF( RoutingInterval(I) .EQ. MissingVal) THEN
C           No routing interval specified
            IF ( LakesToSolve(F(I)) ) THEN
C              Routing through it, but didn't specify  :-(
               WRITE(CTMP,5324) LLT
               CALL LogIt(LOGCDE,1,0)
               WRITE(CTMP,2177) ' No routing interval (i.e. TIMESTEP)',
     &             ' specified for ', NameOfLoop(I),LLT
               CALL LogIt(LOGCDE,0,0)
               CALL SAVETOD
               CALL ADIOS
C           ELSE
C              Not routing through it, didn't specify  :-)
            ENDIF
         ELSE
C           Routing interval was specified
            IF ( LakesToSolve(F(I)) ) THEN
C              Routing through it, did specify  :-)
               WRITE(CTMP,5338) IntervalName(RoutingInterval(I))
     &            (1:IntervalNameLen(RoutingInterval(I))), 
     &            NameOfLoop(I), LLT
               CALL LogIt(LOGUSE,0,0)
            ELSE
C              Not routing through it, did specify    :-(
               WRITE(CTMP,2176) 'Routing Time Step', LLT
               CALL LogIt(LOGCDE,1,0)
               WRITE(CTMP,2177) ' Specified a timestep for ',
     &             NameOfLoop(I),', but not routing through it.', LLT
               CALL LogIt(LOGCDE,0,0)
               CALL SAVETOD
               CALL ADIOS
            ENDIF
         ENDIF
C
C
         IF(RoutingIncrement(I) .EQ. MissingVal ) THEN
C           Increments for routing not specified. 
            IF ( LakesToSolve(F(I)) ) THEN
C              Routing through it, but didn't specify  :-(
               WRITE(CTMP,5324) LLT
               CALL LogIt(LOGCDE,1,0)
               WRITE(CTMP,2177) ' Increments per ',
     &            IntervalName(RoutingInterval(I))
     &            (1:IntervalNameLen(RoutingInterval(I))),
     &            ' routing interval not specified for '//NameOfLoop(I),
     &            LLT 
               CALL LogIt(LOGCDE,0,0)
               CALL SAVETOD
               CALL ADIOS
C           ELSE
C              Not routing through it, didn't specify  :-)
            ENDIF
         ELSE
C           Routing increments were specified
            IF ( LakesToSolve(F(I)) ) THEN
C              Routing through it, did specify  :-)
               IF(RoutingIncrement(I) .EQ. 0 ) THEN
C                 Increments for routing specified as hourly. 
C                 If monthly or qmonthly routing, and using actual month
C                 length, then number of increments changes between 
C                 intervals.  Otherwise stays constant and remains as 
C                 set here.  NumIncsVary tells Runit() which ones that 
C                 need to be recomputed at the start of every interval. 
C           
                  IF( RoutingInterval(I) .EQ. 4) THEN
                     IF( MonthLength .EQ. 1) THEN
C                       Uniform month length, dT equals 59 mins 52.5 s
                        RoutingIncrement(I) = 732
                        NumIncsVary(I) = .FALSE.
                     ELSE 
                        NumIncsVary(I) = .TRUE.
                     ENDIF
C          
                  ELSE IF( RoutingInterval(I) .EQ. 3) THEN
                     IF( MonthLength .EQ. 1) THEN
C                       Uniform month length, dT equals 59 mins 52.5 secs
                        RoutingIncrement(I) = 183
                        NumIncsVary(I) = .FALSE.
                     ELSE 
                        NumIncsVary(I) = .TRUE.
                     ENDIF
C
                  ELSE IF( RoutingInterval(I) .EQ. 2) THEN
                     RoutingIncrement(I) = 168
                     NumIncsVary(I) = .FALSE.
C
                  ELSE IF( RoutingInterval(I) .EQ. 1) THEN
                     RoutingIncrement(I) = 24
                     NumIncsVary(I) = .FALSE.
C
                  ELSE
                     STOP 'Sideshow Mel'
                  ENDIF
                  WRITE(CTMP,2263) NameOfLoop(I), LLT
               ELSE
C                 Number of increments given explicitly
                  NumIncsVary(I) = .FALSE.
                  WRITE(CTMP,5335) RoutingIncrement(I),NameOfLoop(I),LLT
               ENDIF
               CALL LogIt(LOGUSE,0,I/3)
C
            ELSE
C              Not routing through it, but did specify    :-(
               WRITE(CTMP,2176) 'Increments per Time Step', LLT
               CALL LogIt(LOGCDE,1,0)
               WRITE(CTMP,2177) 
     &             ' Specified number of increments per timestep for ',
     &             NameOfLoop(I), ', but not routing through it.', LLT
               CALL LogIt(LOGCDE,0,0)
               CALL SAVETOD
               CALL ADIOS
            ENDIF
         ENDIF
C
 2268    CONTINUE
C
 2263 FORMAT(' Setting increments per ', A, ' routing interval',
     &   ' so that dT is one hour.',A1)
 5335 FORMAT(' Using ',I3,' increments per ', A,' routing interval.',A1)
 5338 FORMAT(' Using ',A,' routing interval for ', A, A1)
C
C
C             <END OF ROUTING INTERVALS AND INCREMENTS SECTION>
C
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvv   viability of dates  vvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C     Set K as error flag
      K = 0
C
C     Look for -9999 in dates     
C
      IF( StartYear .EQ. MissingVal .OR. StartMonth .EQ. MissingVal .OR.
     &   StartDay .EQ. MissingVal ) THEN
         K = 1166
         WRITE(CWRK,6344)' Bad or missing Start Date', LLT
C
      ELSE IF(EndYear .EQ. MissingVal .OR. EndMonth .EQ. MissingVal .OR.
     &   EndDay .EQ. MissingVal ) THEN
         K = 1166
         WRITE(CWRK,6344)' Bad or missing End Date', LLT
C
C     Look at whether start date is greater than ending date:
C
      ELSE IF( EndYear .LT. StartYear .OR. (EndYear .EQ. StartYear
     &   .AND. EndMonth .EQ. StartMonth .AND. EndDay .LE. StartDay) .OR.
     &   (EndYear .EQ. StartYear .AND. EndMonth .LT. StartMonth)) THEN
         K = 1095
         WRITE(CWRK,'(A,I2,A3,I5,A,I2,A3,I5,A1)')' Start date: ', 
     &      StartDay, MonthAbbrev(StartMonth), StartYear, 
     &      '      End date: ', EndDay,
     &      MonthAbbrev(EndMonth), EndYear, LLT
C     ELSE
C        Dates are apparently legitimate; report this condition?
      ENDIF
C   
      IF( K .NE. 0) THEN
         CALL ErrMess(K,CHNKST,LOGCDE)
         CTMP = CWRK
         CALL LogIt(LOGCDE,1,0)
         CALL SAVETOD
         CALL ADIOS
      ENDIF
C
C     Get day of week, relative to Jan 1 of year 0.  Saturday=4
      I = MOD( DateSequence( StartDay, StartMonth, StartYear ), 7)
      DaysSinceSaturday = I - 4
      IF( DaysSinceSaturday .LT. 0) 
     &   DaysSinceSaturday = DaysSinceSaturday + 7
      WRITE(CTMP,3976) DayOfWeekName(I + 1), DaysSinceSaturday, LLT 
      CALL LogIt(LOGD,1,1)
 3976 FORMAT(' Start Date falls on a ', A, ' (', I1, 
     &   ' days since Saturday)', A1)
C
C
C                       <END OF DATE CHECKING SECTION>
C
C
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvv   starting levels and supplies  vvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C     These assignments get cleared in RUNIT
      FWRK(1,1) = SPSTLV
      FWRK(1,2) = MHSTLV
      FWRK(1,3) = MHSTLV
      FWRK(1,4) = MHSTLV
      FWRK(1,5) = SCSTLV
      FWRK(1,6) = ERSTLV
      FWRK(1,7) = ONSTLV
      FWRK(1,8) = MissingVal
      FWRK(1,9) = CGSTLV
      FWRK(1,10) = MissingVal
      FWRK(1,11) = MissingVal
      FWRK(1,12) = MissingVal
      FWRK(1,13) = MissingVal
      FWRK(1,14) = MissingVal
      FWRK(1,15) = MissingVal
      FWRK(1,16) = MissingVal
      IF(ChkVerbosity .GE. 9)THEN
C        Print out all starting levels in memory.
C
         DO 2288, I = 1,16
           WRITE(CTMP, 5337) FWRK(1,I), LAKENAME(I)(1:NLN(I)), LLT
 2288      CALL LogIt(LOGD,0,0)
      ENDIF
C
 5337 FORMAT(' Using starting elevation of ',F7.3,
     &   ' meters IGLD-85 for ',A,'.',A1)
C 
C     K .GT. 0 means to jump to 5320 with fatal error message
      K = 0
C     Don't actually *need* any particular supply or retardation
C     component, but would be very unusual not to have NBS and 
C     certain diversions.  Just spit out a warning if not specified.
C
C     I refers to lake that caused the message
      I = 1
C
C     Check lake superior 
      IF(LakesToSolve(1))THEN
C
C        Check initial level
         IF(IDNINT(SPSTLV) .EQ. MissingVal ) K = 1063
C        Check SUPNBSDATA
         IF(FileNames(42)(1:12) .EQ. Unspecified .AND.  
     &      FileNames(41)(1:12) .EQ. Unspecified) 
     &      CALL IssueWarning(1122,LogCDE,1,CHNKST)
C        Check LONGLACDIVDATA, OGOKIDIVDATA, COMBLLODIVDATA
         IF(FileNames(44)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(45)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(46)(1:12) .EQ. Unspecified .AND.
     &      FileNames(41)(1:12) .EQ. Unspecified) 
     &      CALL IssueWarning(1123,LogCDE,1,CHNKST)
C        Check STMRIVRETRDATA, STMRIVICEDATA, STMRIVWEEDDATA
         IF(FileNames(50)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(51)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(52)(1:12) .EQ. Unspecified) 
     &      CALL IssueWarning(1161,LogCDE,1,CHNKST)
C
      ELSE
         IF(IDNINT(SPSTLV) .NE. MissingVal )  K = 1114
      ENDIF          
      IF(K .GT. 0) GOTO 5320
C
C     Check lakes mh
      I = 2
      IF(LakesToSolve(2))THEN
C        Check initial level
         IF(IDNINT(MHSTLV) .EQ. MissingVal ) K = 1063
C        Check MHNBSDATA
         IF(FileNames( 1)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(31)(1:12) .EQ. Unspecified) 
     &      CALL IssueWarning(1081,LogCDE,2,CHNKST)
C        Check CHICAGODIVDATA
         IF(FileNames(22)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(31)(1:12) .EQ. Unspecified) 
     &      CALL IssueWarning(1084,LogCDE,2,CHNKST)
C        Check STCRIVRETRDATA, STCRIVICEDATA, STCRIVWEEDDATA
         IF(FileNames(16)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(34)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(35)(1:12) .EQ. Unspecified) 
     &      CALL IssueWarning(1161,LogCDE,2,CHNKST)

C        See if Superior is modeled.  If not, then need St. Marys flow
         IF( .NOT. LakesToSolve(1) .AND. 
     &       FileNames(19)(1:12) .EQ. Unspecified .AND. 
     &       FileNames(31)(1:12) .EQ. Unspecified) 
     &       CALL IssueWarning(1184,LogCDE,2,CHNKST)
C
      ELSE
         IF(MHSTLV .NE. MissingVal )  K = 1114
      ENDIF
      IF(K .GT. 0) GOTO 5320
C
C
C     Check lake st clair
      I = 5
      IF(LakesToSolve(3))THEN
C        Check initial level
         IF(IDNINT(SCSTLV) .EQ. MissingVal) K = 1063
C        Check STCNBSDATA
         IF(FileNames( 2)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(32)(1:12) .EQ. Unspecified) 
     &      CALL IssueWarning(1082,LogCDE,5,CHNKST)
C        Check DETRIVRETRDATA, DETRIVICEDATA, DETRIVWEEDDATA 
         IF(FileNames(17)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(36)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(37)(1:12) .EQ. Unspecified) 
     &      CALL IssueWarning(1161,LogCDE,5,CHNKST)
      ELSE
         IF(IDNINT(SCSTLV) .NE. MissingVal) K = 1114
C        should also squawk if gave STCNBSDATA!
      ENDIF
      IF(K .GT. 0) GOTO 5320
C
C     Check lake erie
      I = 6
      IF(LakesToSolve(4))THEN
C        Check initial level
         IF(IDNINT(ERSTLV) .EQ. MissingVal) K = 1063
C        Check ERINBSDATA
         IF(FileNames( 3)(1:12) .EQ. Unspecified .AND.  
     &      FileNames(33)(1:12) .EQ. Unspecified) 
     &      CALL IssueWarning(1083,LogCDE,6,CHNKST)
C        Check WELLANDCANALDATA
         IF(FileNames(21)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(33)(1:12) .EQ. Unspecified) 
     &      CALL IssueWarning(1085,LogCDE,6,CHNKST)
C        Check NIARIVRETRDATA, NIARIVICEDATA, NIARIVWEEDDATA
         IF(FileNames(18)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(38)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(39)(1:12) .EQ. Unspecified) 
     &      CALL IssueWarning(1161,LogCDE,6,CHNKST)
      ELSE
         IF(IDNINT(ERSTLV) .NE. MissingVal) K = 1114
      ENDIF
      IF(K .GT. 0) GOTO 5320
C
C     Check lake ontario
      I = 7
      IF(LakesToSolve(5))THEN
C
C        Check initial level
         IF(IDNINT(ONSTLV) .EQ. MissingVal) K = 1063
C        Check ONTNBSDATA
         IF(FileNames(53)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(57)(1:12) .EQ. Unspecified) 
     &      CALL IssueWarning(1132,LogCDE,7,CHNKST)
C        Check WELLANDCANALDATA
         IF(FileNames(21)(1:12) .EQ. Unspecified .AND. 
     &      FileNames(57)(1:12) .EQ. Unspecified) 
     &      CALL IssueWarning(1133,LogCDE,7,CHNKST)
C
      ELSE
         IF(IDNINT(ONSTLV) .NE. MissingVal) K = 1114
      ENDIF          
C

C
 5320 IF ( K .GT. 0) THEN
         CALL ErrMess(K,CHNKST,LOGCDE)
         WRITE(CTMP,5326) LAKENAME(I)(1:NLN(I)), LLT
         CALL LogIt(LOGCDE,1,0)
         IF(K .EQ. 1063) THEN
            WRITE(CTMP,5336) LAKENAME(I)(1:NLN(I)), FWRK(1,I), LLT
            CALL LogIt(LOGCDE,0,0)
         ENDIF
         CALL SAVETOD
         CALL ADIOS
 5326    FORMAT(' Inconsistent parameters for ',A,A1)
 5336    FORMAT(' Missing or bogus level for ',A,': ',F9.3,A1)
      ENDIF
C
C
C
C                     <END OF WHAT'S NEEDED LEVELS SECTION>
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvv   Middle Lakes Routing   vvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      IF( LakesToSolve(2) .OR. LakesToSolve(3) .OR. LakesToSolve(4))THEN
C
C        Use K and CWRK(1:1) as flags that error was found
         CWRK(1:1) = LLT
         K = 0
C
         IF(LakesToSolve(2) .AND. (
     &      IDNINT(MHK)  .EQ. MissingVal .OR. 
     &      IDNINT(MHYM) .EQ. MissingVal .OR. 
     &      IDNINT(MHA)  .EQ. MissingVal .OR. 
     &      IDNINT(MHB)  .EQ. MissingVal .OR. 
     &      IDNINT(MHWT) .EQ. MissingVal .OR. 
     &      IDNINT(MHC)  .EQ. MissingVal ) )
     &      WRITE(CWRK,2176) 'St. Clair River', LLT
C
         IF(LakesToSolve(3) .AND. (
     &      IDNINT(SCK)  .EQ. MissingVal .OR. 
     &      IDNINT(SCYM) .EQ. MissingVal .OR. 
     &      IDNINT(SCA)  .EQ. MissingVal .OR. 
     &      IDNINT(SCB)  .EQ. MissingVal .OR. 
     &      IDNINT(SCWT) .EQ. MissingVal .OR. 
     &      IDNINT(SCC)  .EQ. MissingVal ) )
     &      WRITE(CWRK,2176) 'Detroit River', LLT
C
         IF(LakesToSolve(4) .AND. (
     &      IDNINT(ERK)  .EQ. MissingVal .OR. 
     &      IDNINT(ERYM) .EQ. MissingVal .OR. 
     &      IDNINT(ERA)  .EQ. MissingVal .OR. 
     &      IDNINT(ERB)  .EQ. MissingVal .OR. 
     &      IDNINT(ERWT) .EQ. MissingVal .OR. 
     &      IDNINT(ERC)  .EQ. MissingVal ) )
     &      WRITE(CWRK,2176) 'Niagara River', LLT
C
         IF(CWRK(1:1) .NE. LLT) K = 1162
 
         IF( IDNINT(ERB * 1000) .EQ. 0) THEN
C           Not Double-gage, don't specify CGIPSTLV and CGIPAVELEVELDATA
            IF( IDNINT(CGSTLV) .NE. MissingVal) K = 1101
            IF(FileNames(20)(1:12) .NE. Unspecified) K = 1102
         ELSE
C           Double-gage!  
C           Must specify CGIPAVELEVELDATA
            IF(FileNames(20)(1:12) .EQ. Unspecified) K = 1100
C           For the matrix solution only, must specify CGIP start level
            IF(MidLakeSolution .EQ. 1 .AND. 
     &         IDNINT(CGSTLV) .NE. MissingVal) K = 1112
            IF(MidLakeSolution .EQ. 2 .AND. 
     &         IDNINT(CGSTLV) .EQ. MissingVal) K = 1099
         ENDIF
         IF(K .NE. 0 .AND. K .NE. 1162)
     &      WRITE(CWRK,2176)'Niagara River',LLT
C
         IF(K .NE. 0) THEN
            CALL ErrMess(K,CHNKST,LOGCDE)
            CTMP = CWRK
            CALL LogIt(LOGCDE,0,0)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
      ENDIF
C
C                     <END OF Middle Lakes Routing SECTION>
C
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvv  Superior and  Plan 1977A parameters   vvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C     See if using outflow equation to get St Marys Flows
      IF( LakesToSolve(1) .AND. StMarysSolution .EQ. 1 ) THEN
         IF(IDNINT(SPK) .EQ. MissingVal .OR.
     &      IDNINT(SPYM) .EQ. MissingVal .OR. 
     &      IDNINT(SPA) .EQ. MissingVal .OR. 
     &      IDNINT(SPB) .EQ. MissingVal .OR. 
     &      IDNINT(SPWT) .EQ. MissingVal .OR. 
     &      IDNINT(SPC) .EQ. MissingVal ) THEN
            CALL ErrMess(1160,CHNKST,LOGCDE)
            CALL SAVETOD
            CALL ADIOS
C        ELSE
C           make sure all plan 77a stuff is Missing (except stuff with defaults)
C
         ENDIF
      ENDIF
C
      IF( LakesToSolve(1) .AND. StMarysSolution .NE. 1) THEN
C        check sw pier gate equation
         IF( LStoSWP .EQ. MissingVal) THEN
            CALL ErrMess(1169,CHNKST,LOGCDE)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
      ENDIF
C
C     See if using a Superior Regulation Plan
C     XX TODO: this clause used to just refer to Plan 1977A, but 
C     XX has been extended to cover any plan.  The checks performed 
C     XX are probably not appropriate to these other plans, and this 
C     XX logic should be cleaned up when he opportunity presents.
      IF( LakesToSolve(1) .AND. StMarysSolution .GE. 2) THEN
C
C        must have some gate equations defined
         IF(IDNINT(SPK) .NE. MissingVal .OR.
     &      IDNINT(SPYM) .NE. MissingVal .OR. 
     &      IDNINT(SPA) .NE. MissingVal .OR. 
     &      IDNINT(SPB) .NE. MissingVal .OR. 
     &      IDNINT(SPWT) .NE. MissingVal .OR. 
     &      IDNINT(SPC) .NE. MissingVal ) THEN
            CALL ErrMess(1163,CHNKST,LOGCDE)
            WRITE(CTMP,2176) 'Lake Superior', LLT
            CALL LogIt(LOGCDE,1,0)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
C
C        XX If a regulation plan is added that doesn't need to solve
C        XX for the middle lakes, this requirement should be ignored
         IF(.NOT. LakesToSolve(2) .OR. .NOT. LakesToSolve(3) .OR. 
     &      .NOT. LakesToSolve(4) ) THEN
            CALL ErrMess(1167,CHNKST,LOGCDE)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
C
C        XX This requirement is appropriate for any plan that involves 
C        XX operational supplies (which are also used to check Criterion
C        XX B, etc.).  Currently, this would include Plan 77A, Plan 2012,
C        XX and the '55 Modified Rule of '49 (as implemented).
         IF( P77ACompat .EQ. MissingVal) THEN
            WRITE(CTMP,5324) LLT
            CALL LogIt(LOGCDE,1,0)
            WRITE(CTMP,2177) ' P77A Compat Mode not specified!', LLT
            CALL LogIt(LOGCDE,0,0)
            CALL SAVETOD
            CALL ADIOS
         ELSE IF( P77ACompat .EQ. -1) THEN
            WRITE(CTMP,5445) ' Using VirtCarl-style operational '//
     &         'supplies (go month-by-month and choose first value)',LLT
         ELSE IF ( P77ACompat .EQ. 0) THEN
            WRITE(CTMP,5445) ' Using actual Superior supplies as '//
     &         'operational supplies (same as previous models)',LLT
         ELSE IF ( P77ACompat .EQ. 1 ) THEN
            WRITE(CTMP,5445) ' Using Plan 1977A 50% supplies',LLT
         ELSE IF ( P77ACompat .EQ. 2 ) THEN
C            WRITE(CTMP,5445) ' Using 6-month sequence of 50% '//
C     &         'operational Sup. supplies',LLT
            WRITE(CTMP,'(A,A1)')' P77A Compat Mode 2 not supported.',LLT
            CALL LogIt(LOGCDE,1,0)
            CALL SAVETOD
            CALL ADIOS
         ELSE IF( P77ACompat .EQ. 3) THEN
            WRITE(CTMP,5445) ' Using VirtCarl-style operational '//
     &         'supplies from matrices used before 2007)',LLT
         ELSE 
            STOP 'Tod Flanders'
         ENDIF
         CALL LogIt(LOGUSE,1,0)
C
         IF( RoutingIncrement(1) .EQ. 1) THEN
            WRITE(CTMP,5445) ' Assuming constant St. Marys River '//
     &         'flows during each routing period',LLT
         ELSE IF( RoutingIncrement(1) .GT. 1 .OR. NumIncsVary(1)) THEN
            WRITE(CTMP,5445) ' Allowing St. Marys River flows to'//
     &         'vary with Superior level within each routing period',LLT
         ELSE 
            STOP 'Ralph Wiggins'
         ENDIF
         CALL LogIt(LOGUSE,0,0)
C
C        See if non-monthly Plan 1977A routing
         IF( RoutingInterval(1) .NE. 4 )CALL ErrMess(1118,CHNKST,LogUse)
C
         IF( SupCritRnd .EQ. MissingVal) THEN
            CALL ErrMess(1168,CHNKST,LOGCDE)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
         IF( SupCritRnd .EQ. 1) THEN
C            Round Criterion B and C calculations to P77A standards
             CritBOMRound = P77AMONTHROUND
             CritFlowRound = P77AFLOWROUND
             CritLevelRound = P77ALEVELROUND
         ELSE
C            Round Criterion B and C calculations to general settings
             CritBOMRound = BOMRound
             CritFlowRound = FlowRound
             CritLevelRound = LevelRound
         ENDIF
C        Make sure last month's flow is provided
         IF( PrevStMarysFlow .EQ. MissingVal) THEN
            CALL ErrMess(1130,CHNKST,LOGCDE)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
C
         IF( Q1887Method .EQ. MissingVal) THEN
            CALL ErrMess(1175,CHNKST,LOGCDE)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
C
         IF( CritBMethod .EQ. MissingVal) THEN
            CALL ErrMess(1176,CHNKST,LOGCDE)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
C
         CWRK = ' CGLRRM will not'
         IF(RoutingIncrement(1) .EQ. 1 .AND. Q1887Method .EQ. 1
     &      .AND. CritBMethod .EQ. 1 .AND. BOMRound .EQ. -2
     &      .AND. RoutingInterval(1) .EQ. 4 
     &      .AND. LStoSWP .EQ. 1 .AND. P77ACompat .LE. 0) THEN 
            BackCompat = .TRUE.
            K = 12
         ELSE
            K = 16
         ENDIF
         WRITE(CTMP,5445) CWRK(1:K)//' replicate the answers '//
     &         'produced with the old upper lakes models.',LLT
         CALL LogIt(LOGUSE,0,1)
C
C        Check for valid input gate setting
         IF( NINT(PrevSupGates ) .NE. MissingVal ) THEN
            IF ((PrevSupGates .GT. 16) .OR. (PrevSupGates .LE. 0)) THEN
               WRITE(CWRK,5663) PrevSupGates, LLT
               K = 1185
               GOTO 5650
            ELSE IF ((PrevSupGates .NE. NINT(PrevSupGates)) .AND.
     &         (PrevSupGates .NE. 0.5)) THEN
               WRITE(CWRK,5663) PrevSupGates, LLT
               K = 1185
               GOTO 5650
            END IF
         ENDIF
C
C        See if beginning-of-month levels are known (needed for P77A)
C        Use K as a flag (i.e., K .NE. 0, abort run)
         K = 0
         CWRK(1:1) = LLT
         IF(StartDay .EQ. 1)THEN
C           Run starts at beginning-of-month, don't need SPBOMLev
            IF( IDNINT(SPBOMLev) .EQ. MissingVal ) THEN
C              BOM not specified ;-), set to start level
               SPBOMLev = SPSTLV
            ELSE 
C              Unnecessary BOM level exists; spit out a warning 
               WRITE(CTMP,5658) LakeName(1)(1:NLN(1)), LLT
               CALL LogIt(LOGUSE,1,0)
C              If unnecessary BOM level exists, it's OK as long as it's
C              the same as start level.
C              ?? Should be able to do direct floating point compares in 
C              ?? this context; failures due to tiny differences of 
C              ?? precision will be obvious because of spurious error   
C              ?? messages. If happens, do something like this line:
C              IF( IDNINT(SPBOMLev * SecsInDay) .NE. IDNINT(SPSTLV * SecsInDay ) THEN
               IF( SPBOMLev .NE. SPSTLV )  THEN
C                 Unnecessary BOM level contradicts start level; abort
                  WRITE(CWRK,5659) LakeName(1)(1:NLN(1)), SPSTLV, 
     &               SPBOMLev, LLT
                  K = 1164
                  GOTO 5650
               ENDIF
            ENDIF
            IF( IDNINT(MHBOMLev) .EQ. MissingVal ) THEN
C              BOM not specified ;-), set to start level
               MHBOMLev = MHSTLV
            ELSE 
C              Unnecessary BOM level exists; spit out a warning 
               WRITE(CTMP,5658) LakeName(2)(1:NLN(2)), LLT
               CALL LogIt(LOGUSE,1,0)
               IF( MHBOMLev .NE. MHSTLV ) THEN
                  WRITE(CWRK,5659) LakeName(2)(1:NLN(2)), MHSTLV, 
     &               MHBOMLev, LLT
                  K = 1164
                  GOTO 5650
               ENDIF
            ENDIF
            IF( IDNINT(SCBOMLev) .EQ. MissingVal ) THEN
C              BOM not specified ;-), set to start level
               SCBOMLev = SCSTLV
            ELSE 
C              Unnecessary BOM level exists; spit out a warning 
               WRITE(CTMP,5658) LakeName(5)(1:NLN(5)), LLT
               CALL LogIt(LOGUSE,1,0)
               IF( SCBOMLev .NE. SCSTLV ) THEN
                  WRITE(CWRK,5659) LakeName(5)(1:NLN(5)), SCSTLV, 
     &               SCBOMLev, LLT
                  K = 1164
                  GOTO 5650
               ENDIF
            ENDIF
            IF( IDNINT(ERBOMLev) .EQ. MissingVal ) THEN
C              BOM not specified ;-), set to start level
               ERBOMLev = ERSTLV
            ELSE 
C              Unnecessary BOM level exists; spit out a warning 
               WRITE(CTMP,5658) LakeName(6)(1:NLN(6)), LLT
               CALL LogIt(LOGUSE,1,0)
               IF( ERBOMLev .NE. ERSTLV ) THEN
                  WRITE(CWRK,5659) LakeName(6)(1:NLN(6)), ERSTLV, 
     &               ERBOMLev, LLT
                  K = 1164
                  GOTO 5650
               ENDIF
            ENDIF
C
C           At this spot, have nothing extra.  Could report that
C           BOM levels are same as initial levels.
C
         ELSE
            CTMP(1:1) = LLT
C           Run doesn't start at BOM; must give BOM levels explicitly
            IF( IDNINT(SPBOMLev) .EQ. MissingVal ) THEN
C              Necessary BOM level not given
               WRITE(CWRK,5657) LakeName(1)(1:NLN(1)), LLT
               K = 1093
            ELSE IF( IDNINT(MHBOMLev) .EQ. MissingVal ) THEN
C              Necessary BOM level not given
               WRITE(CWRK,5657) LakeName(2)(1:NLN(2)), LLT
               K = 1093
            ELSE IF( IDNINT(SCBOMLev) .EQ. MissingVal ) THEN
C              Necessary BOM level not given
               WRITE(CWRK,5657) LakeName(5)(1:NLN(5)), LLT
               K = 1093
            ELSE IF( IDNINT(ERBOMLev) .EQ. MissingVal ) THEN
C              Necessary BOM level not given
               WRITE(CWRK,5657) LakeName(6)(1:NLN(6)), LLT
               K = 1093
            ELSE
C              At this spot, could report that BOM levels have
C              been specified apart from initial levels, as required.
            ENDIF
C
         ENDIF
C
C        If running Plan 2012, need to check that the year of the last
C        June Sturgeon flow is specified.  Then need to calculate 
C        counter value based on that year and the starting year/month 
C        of the simulation.
         IF ( StMarysSolution .EQ. 4 .AND. 
     &      Last_JuneFlow .EQ. MissingVal) THEN
C           Write some error message and exit the program
            WRITE(CWRK,5660) LLT
            K = 1183
            GOTO 5650
         ELSE
C           Check that JuneFlow occurred in the past
            IF ( Last_JuneFlow .GT. StartYear ) THEN
C              Write some error message and die
               WRITE(CWRK,5661) Last_JuneFlow, StartYear, LLT
               K = 1184
               GOTO 5650
            ELSE IF (Last_JuneFlow .EQ. StartYear .AND. 
     &         StartMonth .LE. 6) THEN
C              Write some error message and die
               WRITE(CWRK,5662) Last_JuneFlow, MonthAbbrev(StartMonth),
     &            StartYear, LLT
               K = 1184
               GOTO 5650
            ELSE
C              JuneFlow is a valid number, so calculate the counter
               COUNTER_JUNEFLOW = StartYear - Last_JuneFlow
               IF (StartMonth .LE. 6) THEN
C                 It hasn't been a full year yet, so decrement counter
                  COUNTER_JUNEFLOW = COUNTER_JUNEFLOW - 1
               END IF
            END IF
         END IF
C
 5650    IF ( K .NE. 0 ) THEN
            CALL ErrMess(K,CHNKST,LOGCDE)
            IF(CWRK(1:1) .NE. LLT) THEN
                CTMP = CWRK
                CALL LogIt(LOGUSE,0,0)
            ENDIF
            CALL SAVETOD
            CALL ADIOS
         ENDIF
C
      ENDIF
C
 5657 FORMAT(' No Plan 1977A beginning-of-month level specified for ', 
     &       A, A1) 
 5658 FORMAT(' Warning: Unnecessary BOM level given for ',A, A1)
 5659 FORMAT(1X, A,' BOM level ', F10.4,
     &   ' contradicts initial level ', F10.4, A1)
 5660 FORMAT(1X, 'No value specified for year of last June sturgeon ',
     &   'flow.', A1)
 5661 FORMAT(1X, I5, ' specified for last sturgeon flow, but ',
     &   'simulation starts in ', I5, A1)
 5662 FORMAT(1X, 'Last sturgeon flow given as June, ', I5,
     &   ', but simulation starts in ', A3, ', ', I5, A1)
 5663 FORMAT(1X, 'Previous gate setting of ', A, ' is invalid.', A1)
C8206 FORMAT('  . . . but run starts at beginning of the month,',
C    &       ' so use initial level.', A1) 
C
C
C
C                         <END OF PLAN 1977A SECTION>
C
C
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvv   Lake Ontario parameters   vvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
      IF( LakesToSolve(5) ) THEN
C
C
          IF(UStLawrSolution .EQ. 2)THEN
C
C              P58DStartDay, P58DStartMonth, P58DStartYear are the 
C              variables corresponding to the beginning of Plan 1958D 
C              regulation.  Plan 1958D works off weeks defined from 
C              Saturday at 0000 to  Friday at 2400.  If StartDay is a 
C              Saturday, then simply set P58DStartDay equal to StartDay. 
C              If not, then set P58DStartDay to the **previous** Saturday,
C              per Paul.  
               P58DStartYear = StartYear
               P58DStartMonth = StartMonth
C              Get days until Saturday in K.
C              K = 7 - DaysSinceSaturday
C              IF( K .GT. 6) K = 0
C               P58DStartDay = StartDay + K
C              CALL GetQtrMonBounds (StartYear, StartMonth, U)
C              IF( P58DStartDay .GT. U(4) ) THEN
C                  P58DStartDay = P58DStartDay -  U(4)
C                  P58DStartMonth = P58DStartMonth + 1
C                  IF( P58DStartMonth .GT. 12 ) THEN
C                      P58DStartMonth = 1
C                      P58DStartYear = P58DStartYear + 1
C                  ENDIF
C              ENDIF
C              Commented stuff above was for case when date was shifting 
C              to NEXT Saturday.
               P58DStartDay = StartDay - DaysSinceSaturday
               IF( P58DStartDay .LT. 1 ) THEN
                   P58DStartMonth = P58DStartMonth - 1
                   IF( P58DStartMonth .LT. 1 ) THEN
                       P58DStartMonth = 12
                       P58DStartYear = P58DStartYear - 1
                   ENDIF
                   CALL GetQtrMonBounds(P58DStartYear,P58DStartMonth, U)
                   P58DStartDay = U(4) + P58DStartDay
               ENDIF
               WRITE(CTMP,3377)  P58DStartDay,
     &            MonthAbbrev(P58DStartMonth), P58DStartYear, LLT
              CALL LogIt(LOGCD,1,1)
 3377 FORMAT(' First Ontario regulation week begins on Saturday, ',
     &          I2, A3, I5, A1)
C
C     Using Plan 1958D for Lake Ontario Regulation
C
              WRITE(CTMP,6344) 'Plan 58D initial conditions: '//LLT
              CALL LogIt(LOGD,1,0)
              k1=RoutingInterval(3)
              IF (RoutingInterval(3) .NE. 2) THEN
                  WRITE(CTMP,5445)' If using Plan 58D, Ontario '//
     &           'routing interval must be weekly, but you specified: ',
     &            IntervalName(k1)(1:IntervalNameLen(k1)),LLT
                  CALL LogIt(LOGCDE,0,0)
                  CALL ErrMess(1147,CHNKST,LOGCDE)
                  CALL SAVETOD
                  CALL ADIOS
              ENDIF
C
C
C
c--- Write statements for the detailed log file
              WRITE(CTMP,6347) 'Using weekly routing through Lake '//
     &        'Ontario, with: ', RoutingIncrement(3), ' increments', LLT
              CALL LogIt(LOGD,0,0)
              IF (RoutingIncrement(3) .NE. 1 ) THEN
                  CALL LogIt(LOGCE,1,0)
                  CALL ErrMess(1181,CHNKST,LOGCDE)
                  CALL SAVETOD
                  CALL ADIOS
              ENDIF
              WRITE(CTMP,6347) 'P58D Run Type: ', P58DRunType, ' ', LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6346) 'Start level: ', ONSTLV, ' m', LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6346)'Corresponding area: ',GetArea(5,ONSTLV),
     &            ' m2', LLT
              CALL LogIt(LOGD,0,0)
c             WRITE(CTMP,6277) ' Forecast date: ', P58DForecastDay, 
c    &            MonthAbbrev(P58DForecastMon), P58DForecastYear, LLT
c             CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6348) 'Use DATABASE to obtain '//
     &            'Starting Conditions: ', P58DUseDatabase, LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6348) 'Use SI DATABASE to obtain '//
     &            'SI Starting Conditions: ', P58DSiDatabase, LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6347) 'Critical depths: ', CritDepths, ' ', LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6347) 'App Limit: ', P58DAppLimit, ' m3s',  LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6347) 'Total Supply: ', P58DTotalSup,' m3s',LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6347) 'Maximum I limit: ', MaxILimit,' m3s',LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6347) 'Minimum P limit: ', MinPLimit,' m3s',LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6347) 'Weighted supply: ',P58DWeightedSupply,
     &            ' m3s',LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6347) 'Total Deviations: ', P58DTotalDev,
     &            ' m3s',LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6347)'Supply Adjustment: ',P58DAdjustment,
     &            ' m3s',LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6347)'RC Adjustment: ',P58DrcAdjustment,
     &            ' m3s',LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6347) 'Rule curve: ', P58DRuleCurve, ' m3s',LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6346) 'Plan Level: ', P58DPlanLevel, ' m',LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6346) 'Actual Level: ',P58DActualLevel,' m',LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6347)'Actual flow: ',P58DActualflow,' m3s',LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6347)'Lake St. Louis Flow: ',P58DLSLflow,
     &            ' m3s',LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6348)'Set plan level = '// 
     &            'actual level (elim dev): ', ElimDev, LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6348)'Set flow = plan flow if total dev < 0: ', 
     &            SetFlow, LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6348) 'Ignore above levels for '//
     &            'forced flows below: ', ForceFlow, LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6348) 'Pre-project conditions if '//
     &           '74.15 > level > 75.37: ', PreProjConditions, LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6343) 'Input file for Simulation: ', 
     &            P58DSimulateFile(1:LastKC(P58DSimulateFile,256)), LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6343)'P58D App Limit Sym: ',
     &            P58DAppLimitSym(1:2),LLT
              CALL LogIt(LOGD,0,0)
              WRITE(CTMP,6348) 'Compute River Profiles '//
     &            'Pass Montreal: ', ComputeMontreal, LLT
              CALL LogIt(LOGD,0,0)
C
C          
              IF( ElimDev ) THEN
                 IF( CritKFile(1:1) .EQ. ' ' ) THEN
                    WRITE(CTMP,5324) LLT
                    CALL LogIt(LOGCDE,1,0)
                    WRITE(CTMP,5445)' No data specified for '//
     &                 'CritKFile',LLT
                    CALL LogIt(LOGCDE,0,0)
                    CALL SAVETOD
                    CALL ADIOS
                ELSE
                    WRITE(CTMP,6343) 'Input file for Criterion K: ', 
     &                CritKFile(1:LastKC(CritKFile,256)), LLT
                    CALL LogIt(LOGD,0,0)
                ENDIF
              ELSE
                 IF( CritKFile(1:1) .NE. ' ' ) THEN
                    WRITE(CTMP,2176) 
     &                 CritKFile(1:LastKC(CritKFile,256)),LLT
                    CALL LogIt(LOGCDE,1,0)
                    WRITE(CTMP,5445)' Parameter P58D Eliminate Dev'//
     &                 ' must be .TRUE. in order to specify '//
     &                 'parameter CritKFile',LLT
                    CALL LogIt(LOGCDE,0,0)
                    CALL SAVETOD
                    CALL ADIOS
                ENDIF
              ENDIF

c
c Checking to for Lake Ontario Supplies data
c    Filenames(89) = "Onttotalsupdata", if exists, will provide the total
c                    supplies to Lake Ontario.  This will skip Niagara River
c                    and the Welland values
c    Filenames(53) = "Ontnbsdata", if exists, will be the LOCAL Lake Ontario
c                    supplies and MUST be added along with Niagara River and
c                    the Welland values
c
              if( (filenames(89)(1:12).ne.Unspecified) .and.
     &            (filenames(53)(1:12).ne.Unspecified) )then
c----------------
C                 open(881,file='zfile.txt')
C                 write(881,1881)filenames(89)
C                 write(881,1881)filenames(53)
C1881 format(1x,a)                  
C                 close(881)
c----------------
c     
c    Both filenames cannot exists together!!!
                  WRITE(CTMP,5324) LLT
                  CALL LogIt(LOGCDE,1,0)
                  WRITE(CTMP,5445)' No data specified for '//
     &                '[ONLY 1 Lake Ontario Supplies data can exists]',
     &                 LLT
                  CALL LogIt(LOGCDE,0,0)
                  CALL SAVETOD
                  CALL ADIOS
              endif
              if( (filenames(89)(1:12).eq.Unspecified) .and.
     &            (filenames(53)(1:12).eq.Unspecified) ) then
c     
c    One of these file MUST exists!!!
                  WRITE(CTMP,5324) LLT
                  CALL LogIt(LOGCDE,1,0)
                  WRITE(CTMP,5445)' No data specified for '//
     &                '[MUST have a Lake Ontario SUPPLIES data file]',
     &                LLT
                  CALL LogIt(LOGCDE,0,0)
                  CALL SAVETOD
                  CALL ADIOS
              endif

C             If "UseDatabase" or "SiDatabase" variable is true, then it will read 
C             the supply indicators (ie. SI) from the database (ie. "wk5mSI.txt")
              IF(P58DUseDatabase .OR. P58DSiDatabase)THEN
                 WRITE(CTMP,6343) 'Supply indicator database file: ', 
     &               P58DSIDataFile(1:LastKC(P58DSIDataFile,256)), LLT
                 CALL LogIt(LOGD,0,0)
                 IF( P58DSIDataFile(1:1) .EQ. ' ' ) THEN 
                    WRITE(CTMP,5324) LLT
                    CALL LogIt(LOGCDE,1,0)
                    WRITE(CTMP,5445)' No data specified for '//
     &                  'P58DSIDataFile',LLT
                    CALL LogIt(LOGCDE,0,0)
                    CALL SAVETOD
                    CALL ADIOS
                 ENDIF
C                TODO:  check to see if SupplyIndicator were read in, and quit if so.
              ELSE
C                SI must be specified by parameter in input file
                 WRITE(CTMP,6349) (SupplyIndicator(I), I=1,21), LLT
                 CALL LogIt(LOGD,0,0)
C                TODO:  check to see if P58DSIDataFile was given, and quit if so.
              ENDIF
c
c These statements is no longer needed due to "COMPUTEMONTREAL" statements!
c Must specify StLouisRetrData or -9999 causes math error
c
c             IF(FileNames(62)(1:12) .EQ. Unspecified)THEN
c                 WRITE(CTMP,5324) LLT
c                 CALL LogIt(LOGCDE,1,0)
c                 WRITE(CTMP,5445)' No data specified for '//
c    &                'StLouisRetard',LLT
c                 CALL LogIt(LOGCDE,0,0)
c                 CALL SAVETOD
c                 CALL ADIOS
c             ENDIF
c
c Need to check if the other datafiles are included if "ComputeMontreal"
c is set to TRUE.
c   
              if (computemontreal) then
                  IF(FileNames(76)(1:12) .EQ. Unspecified)THEN
                      WRITE(CTMP,5324) LLT
                      CALL LogIt(LOGCDE,1,0)
                      WRITE(CTMP,5445)' No data specified for '//
     &                    'LSL-Ont. Diff. Data',LLT
                      CALL LogIt(LOGCDE,0,0)
                      CALL SAVETOD
                      CALL ADIOS
                  ENDIF
                  do 919, idd=78,88
                      IF(FileNames(idd)(1:12) .EQ. Unspecified)THEN
                          WRITE(CTMP,5324) LLT
                          CALL LogIt(LOGCDE,1,0)
                          if(idd.eq.78)then
                              WRITE(CTMP,5445)' No data specified for '
     &                            //'[Flow DPMI Data]',LLT
                          elseif(idd.eq.79)then
                              WRITE(CTMP,5445)' No data specified for '
     &                            //'[Flow StFran Data]',LLT
                          elseif(idd.eq.80)then
                              WRITE(CTMP,5445)' No data specified for '
     &                            //'[Flow StMau Data]',LLT
                          elseif(idd.eq.81)then
                              WRITE(CTMP,5445)' No data specified for '
     &                            //'[Flow Rich Data',LLT
                          elseif(idd.eq.82)then
                              WRITE(CTMP,5445)' No data specified for '
     &                            //'[Rough Jetty1 Data]',LLT
                          elseif(idd.eq.83)then
                              WRITE(CTMP,5445)' No data specified for '
     &                            //'[Rough Varen Data]',LLT
                          elseif(idd.eq.84)then
                              WRITE(CTMP,5445)' No data specified for '
     &                            //'[Rough Sorel Data]',LLT
                          elseif(idd.eq.85)then
                              WRITE(CTMP,5445)' No data specified for '
     &                            //'[Rough Pierre Data]',LLT
                          elseif(idd.eq.86)then
                              WRITE(CTMP,5445)' No data specified for '
     &                            //'[Rough 3Rivers Data]',LLT
                          elseif(idd.eq.87)then
                              WRITE(CTMP,5445)' No data specified for '
     &                            //'[Rough Batis Data]',LLT
                          elseif(idd.eq.88)then
                              WRITE(CTMP,5445)' No data specified for '
     &                            //'[tidal Data]',LLT
                          endif
                          CALL LogIt(LOGCDE,0,0)
                          CALL SAVETOD
                          CALL ADIOS
                      endif
  919             continue
              ENDIF
c
              WRITE(CTMP, 3457) P58DStartDay, 
     &             MonthAbbrev(P58DStartMonth), P58DStartYear, LLT
              CALL LOGIT(LOGD,0,1)
c
c        Run Paul's initialization routine
c
              EElimDev=Elimdev
              FForceFlow=ForceFlow
              SSetFlow=SetFlow
              UseDatabase=P58DUseDatabase
              SiDatabase=P58DSiDatabase
              onstlv1=onstlv
          endif
      elseif(UStLawrSolution .EQ. 5)then
          WRITE(CTMP,6344) 'Preproject initial conditions: '//LLT
          CALL LogIt(LOGD,1,0)
          WRITE(CTMP,6346) 'Preproject Level: ',PreprojectLevel,' m',LLT
C         XXX if PreprojectLevel not defined exit with message
          CALL LogIt(LOGD,0,0)
          WRITE(CTMP, 3457) I, MonthAbbrev(J), K, LLT
          CALL LOGIT(LOGD,0,1)
      ENDIF
      I = StartDay - DaysSinceSaturday
      J = StartMonth
      K = StartYear
      IF(I .LE. 0) THEN
          J = J - 1
          IF( J .LT. 1) THEN
              J = 12
              K = K - 1
          ENDIF
          CALL GetQtrMonBounds( K, J, U)
          I = U(4) + I
      ENDIF
      call InitOntMain (K, J, I)
 3457 FORMAT('Ontario initialization data applies to week beginning ',
     &   I2,A3,I5,A1)
 6343 FORMAT(A51,2X,A,A1)
 6344 FORMAT(A,A1)
 6346 FORMAT(A51,F17.4,A,A1)
 6347 FORMAT(A51,I12,A,A1)
 6348 FORMAT(A51,11X,L1,A1)
 6349 FORMAT(' Supply Indicators: ',21I4,A1)
C
C
C                         <END OF PLAN 1958D SECTION>
C
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvv   output file checks   vvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
      DO 2027, L=1,8 
         DO 2027, I=1,4
            DO 2027, J=1,2
               IF( WriteTSOutput(I,J,L)) THEN
C                 So we asked for this output.
C                     K = TSOutputFileName(I, J, L, TmpStr)
C                     WRITE(CTMP, 
C     &          '(A,I1,A,I1,A,I1,1X,A,1X,A,1X,A,1X,I1,1X,A,1X,A,A1)') 
C     &                  ' Evaluating I=', I, ' J=', J, ' L=', L,
C     &                  LakeName(LkToSlv2LN(L))(1:NLN(LkToSlv2LN(L))),
C     &              IntervalName(I)(1:IntervalNameLen(I)),TmpStr(1:K),
C     &                 RoutingInterval(LkSlv2RtLoopI(L)),
C     &                 IntervalName(RoutingInterval(LkSlv2RtLoopI(L)))
C     &          (1:IntervalNameLen(RoutingInterval(LkSlv2RtLoopI(L)))),
C     &                   LLT
C                     CALL LogIt(LOGD,0,0)
C
C                 Are we solving the lake? 
                  IF( .NOT. LakesToSolve(L)) THEN
                     CALL ErrMess(1128,CHNKST,LOGCDE)
                     K = TSOutputFileName(I, J, L, TmpStr)
                     WRITE(CTMP, '(A,A,A,A,A)') ' Not solving ',
     &                  LakeName(LkToSlv2LN(L))(1:NLN(LkToSlv2LN(L))),
     &                  ', so can not make output file ',TmpStr(1:K),LLT
                     CALL LogIt(LOGCDE,0,0)
                     CALL SAVETOD
                     CALL ADIOS
                  ENDIF 
C                 Are we asking for output with a shorter interval than
C                 the routing computations?  It's OK to write quarter-
C                 monthly output from a daily run, but not vice versa.
                  IF(RoutingInterval(LkSlv2RtLoopI(L)) .GT. I) THEN
                     CALL ErrMess(1129,CHNKST,LOGCDE)
                     WRITE(CTMP, '(A,A,A,A,1X,A,A,A1)') 
     &                  ' Unable to generate ',
     &                  IntervalName(I)(1:IntervalNameLen(I)),
     &                  ' output from ', 
     &                  LakeName(LkToSlv2LN(L))(1:NLN(LkToSlv2LN(L))),
     &                  IntervalName(RoutingInterval(LkSlv2RtLoopI(L)))
     &                  (1:IntervalNameLen(RoutingInterval(
     &                  LkSlv2RtLoopI(L)))), ' routing.', LLT
                     CALL LogIt(LOGCDE,0,0)
                     K = TSOutputFileName(I, J, L, TmpStr)
                     WRITE(CTMP, '(A,A,A1)') 
     &                  ' Improper output file:', TmpStr(1:K), LLT
                     CALL LogIt(LOGCDE,0,0)
                     CALL SAVETOD
                     CALL ADIOS
                  ENDIF 
               ENDIF 
 2027    CONTINUE    
C
      L = 0
C     See if we are writing Lake Erie output files.  Set L=1 if we are.
      DO 2029, I=1,4
         DO 2029, J=1,2
 2029      IF( WriteTSOutput(I,J,4)) L = 1
      IF(WriteErieOutflow)THEN
         IF( L .GT. 0)THEN
            WRITE(CTMP,5445) ' Lake Erie Outflow files will be '// 
     &          'written instead of Niagara River flow files.', LLT
            CALL LogIt(LOGCD,1,0)
         ELSE
            WRITE(CTMP,2176) 'WriteErieOutflow', LLT
            CALL LogIt(LOGCDE,1,0)
            WRITE(CTMP,2177) ' Specified ', 'WriteErieOutflow: True', 
     &          ', but did not call for Erie output files!', LLT
            CALL LogIt(LOGCDE,1,0)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
      ELSE
          IF( L .GT. 0)THEN
            WRITE(CTMP,5445) ' Niagara River flow files will be '// 
     &          'written instead of Lake Erie Outflow files.', LLT
            CALL LogIt(LOGCD,1,0)
          ENDIF
      ENDIF 
C
C                     <END OF OUTPUT FILE SECTION>  
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvv   establish summary input filenames   vvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C     FileNames(31 -33) are summary input filenames for MH, St Clair, 
C     and Erie.  If they have been specified, issue message.
C     If not, then fill them in, and issue message.
C     Use U array for scratch, with FileNames index to Summary Input 
C     files for each lake.
C
C     Lake Superior
      U(1) = 41
C     Lakes Michigan-Huron
C      F(2) = 2
      U(2) = 31
C     Lake St. Clair
C      F(3) = 5
      U(3) = 32
C     Lake Erie
C      F(4) = 6
      U(4) = 33
C     Lake Ontario
C      F(5) = 7
      U(5) = 57
C
C     Use J to hold the desired routing interval.  Init to superior.
      J = RoutingInterval(1)
      DO 5345, I=1,5
         IF( LakesToSolve(I) ) THEN
C           Very kludgy logic with these subscripts.  Will need more
C           duct tape if looping changes.
            IF( I .GT. 1) J = RoutingInterval(2)
            IF( I .GT. 4) J = RoutingInterval(3)
            IF( FileNames(U(I))(1:12) .EQ. Unspecified )THEN
C              No intermediate supply file given; use default like 
C              'monthly.m-h'.  Terminate the name specially to 
C              indicate not a specified filename
C              FileNames(U(I))= IntervalName(J)(1:IntervalNameLen(J))//
C    &            '.'//FileExtension(LkToSlv2LN(I))//LLT//'default'
C
C              No intermediate supply file given; use default like 
C              'output\mhqtsd.s11'.  Terminate the name specially to 
C              indicate not a specified filename
               K = TSOutputFileName(J, 3, I, FileNames(U(I)))
               FileNames(U(I)) = FileNames(U(I))(1:K)//LLT//'default'
               WRITE(CTMP,5350) 
     &            LAKENAME(LkToSlv2LN(I))(1:NLN(LkToSlv2LN(I))), LLT
               CALL LogIt(LOGD,1,0)
               WRITE(CTMP, '(A)')' so program will build one.'//LLT
               CALL LogIt(LOGD,0,0)
               WRITE(CTMP,5351) 'default', 
     &            FileNames(U(I))(1:K), LLT
            ELSE
C              Existing intermediate supply file specified.
               WRITE(CTMP,5353) 
     &            LAKENAME(LkToSlv2LN(I))(1:NLN(LkToSlv2LN(I))), LLT
               CALL LogIt(LOGD,1,0)
               WRITE(CTMP,5351) 'input-specified', 
     &            FileNames(U(I))(1:LastKC(FileNames(U(I)),256)), LLT
            ENDIF
            CALL LogIt(LOGD,0,0)
         ENDIF
 5345    CONTINUE
C 
C
 5350 FORMAT(' No existing summary inflow/outflow file specified for ',
     &   A, ',', A1)
 5351 FORMAT(' Will use ',A,' filename for time-series data summary: {', 
     &   A, '}', A1)
C 5352 FORMAT(' filename: {', A, '}', A1)
 5353 FORMAT(' Existing summary inflow/outflow file specified for ',
     &   A, '.', A1)
 5323 FORMAT(' Note: optional input not used (', A, ')', A1)
 5324 FORMAT(' Fatal condition: required input not received!', A1)
 2176 FORMAT(' Inconsistent parameters for ',A,A1)
 2177 FORMAT (A,A,A,A1)
C5422 FORMAT(' Inconsistent input error!  ', 
C    &           'Unnecessary parameter given: ', A, A1)
C
C                    <END OF ESTABLISH SUMMARY INPUT FILES>
C
C
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     vvvvvvvvvvv   Adjust LOGCOD Flags for Optional Files   vvvvvvvvvvv
C     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C                  col 1 - stdout (unit 6)
C                  col 2 - detailed log file (unit 90)
C                  col 3 - critical failure log (unit 91)
C                  col 4 - summary output file (unit 92)
C                  col 5 - Middle Lakes output file (unit 93)
C                  col 6 - Lake Superior output file (unit 94)
C                  col 7 - Lake Ontario output file (unit 95)
C                  col 8 - St. Lawrence output file (unit 96)
C
C     LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
C    &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
C    &             LOGM(8), LOGP(8), LOGO(8), 
C    &             LOGDM(8), LOGDP(8), LOGDO(8), 
C    &             LOGSM(8), LOGSP(8), LOGSO(8),
C    &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
C    &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C
      IF(.NOT. LakesToSolve(2) .AND. .NOT. LakesToSolve(3) .AND. 
     &   .NOT. LakesToSolve(4) ) THEN
C        Turn off all elements for middle lakes output
         LOGM(5) = .FALSE.
         LOGDM(5) = .FALSE.
         LOGSM(5) = .FALSE.
         LOGSDM(5) = .FALSE.
         LOGSMPO(5) = .FALSE.
         LOGDMPO(5) = .FALSE.
         LOGSDMPO(5) = .FALSE.
      ENDIF
C
      IF( .NOT. LakesToSolve(1) ) THEN
C        Turn off all elements for Superior output 
         LOGP(6) = .FALSE.
         LOGDP(6) = .FALSE.
         LOGSP(6) = .FALSE.
         LOGSDP(6) = .FALSE.
         LOGSMPO(6) = .FALSE.
         LOGDMPO(6) = .FALSE.
         LOGSDMPO(6) = .FALSE.
      ENDIF
C
      IF( .NOT. LakesToSolve(5) ) THEN
C        Turn off all elements for Ontario output 
         LOGO(7) = .FALSE.
         LOGDO(7) = .FALSE.
         LOGSO(7) = .FALSE.
         LOGSDO(7) = .FALSE.
         LOGSMPO(7) = .FALSE.
         LOGDMPO(7) = .FALSE.
         LOGSDMPO(7) = .FALSE.
      ENDIF
C
C     St. Lawrence output not enabled yet, so don't need to adjust
C
C
C
C                    <END OF ESTABLISH SUMMARY INPUT FILES>
C
C
      WRITE(CTMP, '(A)')' End of module CHECKPAR()'//LLT
      CALL LogIt(LOGD,1,1)
C
      END
C
C
C----&------------------------------------------------------------------
C
C
C
C
C
      SUBROUTINE IssueWarning (MessageNum, LogsToWrite, LakeNum,
     &   CallingFunction)
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C       
      INTEGER      MessageNum, LakeNum
C 
C
      LOGICAL      LogsToWrite(8)
C 
C
      CHARACTER*16 CallingFunction
C 
C
C      
C     Declarations regarding COMMON block variables:
C
      INTEGER      NLN(16), LkToSlv2LN(8), NRN(16), LkSlv2RtLoopI(8)
C                  See descriptions in Block Data Genral
C
      CHARACTER*21 LAKENAME(16), RiverName(16)
C                  See description in Block Data Genral
C
      CHARACTER    CTMP*255, CWRK*80, LLT*1
C                  These are just buffers used throughout program.
C
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /LNAMEC/   LAKENAME, RiverName
      COMMON /LNAMEI/   NLN, LkToSlv2LN, NRN, LkSlv2RtLoopI
C
C
      CALL ErrMess(MessageNum, CallingFunction, LogsToWrite)
      WRITE(CTMP,5325) LAKENAME(LakeNum)(1:NLN(LakeNum)), LLT
      CALL LogIt(LogsToWrite,0,1)
 5325 FORMAT(' Warning regarding parameters for ',A,A1)
      END
C
C
C
C
C-----------------------------------------------------------------------
C
C
      LOGICAL FUNCTION ReadWeekly(Line, Year, Month, Day, Float)
C
C     This routine reads a line from file FUNIT according to standard
C     format for weekly data (one line each week, 1 values on line),
C     i.e.:  yr, mon, day, value ==> I5,I3,I3,1(F8)
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER Line*(*)
C               (Input) The data line passed to this subroutine.
C
      INTEGER   Year, Month, Day
C               (Output) The year, month, and day read from the 
C               data line if parsed OK, or the IOSTAT value if READ() 
C               fails (the function returns in this case .FALSE.).
C
      DOUBLE PRECISION Float
C               (Output) The 64-bit floating point value read from
C               the data line.
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER   I, LDY, LDM
C 
C
      ReadWeekly = .FALSE.
C
      READ(Line,'(I5, I3, I3, D8.0)',ERR=233,IOSTAT=I) Year, Month,
     &    Day, Float
C
C     See if got bogus month/year
      IF( Year .GE. 1 .AND. Month .GE. 1 .AND. Month .LE. 12) THEN
C
C        Get last day of month in LDM
         CALL NUMDAYS(Year,Month,LDY,LDM)
C        See if got bogus day
         IF( Day .GE. 1 .AND. Day .LE. LDM ) THEN
            ReadWeekly = .TRUE.
            RETURN
         ENDIF
      ENDIF
C   
  233 Year = I
      END
C      
C
C
C-----------------------------------------------------------------------
C
C
      LOGICAL FUNCTION ReadDaily(Line, Year, Month, Qmonth, Floats)
C
C     This routine reads a line from file FUNIT according to standard
C     format for daily data (one line each qmonth, 7-8 values on line),
C     i.e.:  yr, mon, qmonth, values ==> I5,I3,I2,8(F8)
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER Line*(*)
C               (Input) The data line passed to this subroutine.
C
      INTEGER   Year, Month, Qmonth
C               (Output) The year, month, and quarter-month read from the 
C               data line if parsed OK, or the IOSTAT value if READ() 
C               fails (the function returns in this case .FALSE.).
C
      DOUBLE PRECISION Floats(8)
C               (Output) The 7 or 8 64-bit floating point values read from
C               the data line.  Compilers may squawk about re-dimensioning
C               here, since usually uses a column from the much bigger
C               array of scratch 64-bit floats (FWRK), but the usage seems
C               consistent with the ANSI F77 standard, so don't let the
C               warnings bother you too much.
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER   I 
C 
C
      ReadDaily = .FALSE.
C
      READ(Line,'(I5, I3, I2, 8D8.0)',ERR=233,IOSTAT=I) Year, Month,
     &    Qmonth, Floats
C
C     See if got bogus dates
      IF( Year .GE. 1 .AND. (Month .GE. 1 .AND. Month .LE. 12) .AND. 
     &   (Qmonth .GE. 1 .AND. Qmonth .LE. 4) ) THEN
         ReadDaily = .TRUE.
         RETURN
      ENDIF
   
  233 Year = I
      END
C      
C
C
C
C-----------------------------------------------------------------------
C
C
      LOGICAL FUNCTION READMONTHLY(LINE, YEAR, FLOATS)
C
C     This routine reads a line from file FUNIT according to the standard
C     format for quarter-monthly data, ie, (I5, 12D8.0) as of Feb 2008
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER LINE*(*)
C               (Input) The data line passed to this subroutine.
C
      INTEGER   YEAR
C               (Output) The year read from the data line if parsed OK,
C               or the IOSTAT value if READ() fails (the function returns
C               in this case .FALSE.).
C
      DOUBLE PRECISION FLOATS(12)
C               (Output) The 12 64-bit floating point values read from
C               the data line.  Compilers may squawk about re-dimensioning
C               here, since usually uses a column from the much bigger
C               array of scratch 64-bit floats (FWRK), but the usage seems
C               consistent with the ANSI F77 standard, so don't let the
C               warnings bother you too much.
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER   I 
C
C
      READMONTHLY = .FALSE.
C
      I = ICHAR( LINE(5:5) )
C             
      IF(I .GE. 48 .AND. I .LE. 57) THEN
C        It's a number - must be five-digit year format
         READ(LINE,'(I5, 12D8.0)',ERR=233,IOSTAT=I)YEAR,FLOATS

      ELSE IF( LINE(5:5) .EQ. ' ' ) THEN
C        Assume it's the old 4-digit year format
         READ(LINE,'(I4, 12D8.0)',ERR=233,IOSTAT=I)YEAR,FLOATS

      ELSE
C        Set a recognizable flag value
         I = -9999
      ENDIF
C
C
C     Test if got bogus YEAR
      IF( YEAR .GE. 1 .AND. YEAR .LE. 99999 ) THEN
         READMONTHLY = .TRUE.
         RETURN
      ENDIF
   
  233 YEAR = I
      READMONTHLY = .FALSE.
      END
C      
C
C
C
C-----------------------------------------------------------------------
C
      LOGICAL FUNCTION READQMONTHLY(LINE, YEAR, QTR, FLOATS)
C
C     This routine reads four lines from file FUNIT according to the 
C     standard format for quarter-monthly data, ie, (I5, I2, 12D8.0)
C     or (I4, I2, 12D8.0) as of feb 2008.
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER LINE*(*)
C               (Input) The data line passed to this subroutine.
C
      INTEGER   YEAR
C               (Output) The year read from the data line.
C
      INTEGER   QTR
C               (Output) The quarter of the month read from the data line.
C               Quarter-monthly data comes in four records per year, with
C               the first quarter of each month first record, second 
C               quarter-month values on the second, etc.  
C
      DOUBLE PRECISION FLOATS(48)
C               (Output) The 48 64-bit floating point values read from
C               the data lines.
C
      INTEGER I, J
      DOUBLE PRECISION TMPFLOATS(12)
C               Scratch variables. 
C
      READQMONTHLY = .FALSE.
C
      I = ICHAR( LINE(5:5) )
C             
      IF(I .GE. 48 .AND. I .LE. 57) THEN
C        It's a number - must be five-digit year format
         READ(LINE,'(I5, I2, 12D8.0)',ERR=233)YEAR,QTR,TMPFLOATS

      ELSE IF( LINE(5:5) .EQ. ' ' ) THEN
C        Assume it's the old 4-digit year format
         READ(LINE,'(I4, I2, 12D8.0)',ERR=233)YEAR,QTR,TMPFLOATS

      ELSE
C        Set a recognizable flag value
         I = -9999
      ENDIF
C
C
C     Re-order the values to be sequential, and store in FLOATS
      DO 455, J=1,12
  455    FLOATS( (J - 1) * 4 + QTR ) = TMPFLOATS(J)
C
C      
C      WRITE(90, '(A)')' READQMONTHLY():'
C      WRITE(90, '(A)')'    INLINE:', LINE
C      WRITE(90, '(A,12F9.1)')'    FLOATS:',(FLOATS(I), I=1,12)
C     Test if got bogus YEAR
      IF( YEAR .GE. 1 .AND. YEAR .LE. 99999 ) THEN
         READQMONTHLY = .TRUE.
         RETURN
      ENDIF
   
      YEAR = I
      READQMONTHLY = .FALSE.
  233 END
C
C
C-----------------------------------------------------------------------
C
      CHARACTER*256 FUNCTION OutFileName( RootName, RootNameLen, 
     &   OutNameLen)
C
C     This routine pre-pends RCOUTDIR (if applicable) and appends 
C     RCOUTEXT (if applicable) to the root filename given as argument 1.
C     The resulting filespec is returned in OutFileName, with the 
C     position of the final character in OutNameLen.
C
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER RootName*(*)
C               (Input)  The base name (less extension) to use. For 
C               example, if constructing an output file name of
C               "/wcds/simulate/plan9/superior.pl9", then RootName
C               would be "superior".  The "/wcds/simulate/plan9/"
C               comes from RCOUTDIR, while the ".pl9" comes from
C               RCOUTEXT.
C               
      INTEGER   RootNameLen
C               (Input)  The number of characters in RootName.
C
      INTEGER   OutNameLen
C               (Output)  The number of characters in OutFileName.
C
C
C     Declarations regarding local variables or functions:
C
C     (none)
C
C     Declarations regarding COMMON block variables:
C
C      CHARACTER CTMP*255, CWRK*80
C               These are just buffers used throughout program.
C
C      CHARACTER LLT*1
C               This is the character that marks the end of a line sent
C               to the LogIt subroutine.
C
      CHARACTER RCOUTDIR*128, RCOUTEXT*16
C               See descriptions in Block Data Genral
C
      INTEGER   NOPDIR, NOPEXT, TSOutputFileNum(4,2,8)
      LOGICAL   OnlyOneOpen
C               See descriptions in Block Data Genral
C
C
      COMMON /OUTNAM/   RCOUTDIR, RCOUTEXT
C      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /OUTNMI/   NOPDIR, NOPEXT, TSOutputFileNum, OnlyOneOpen
C
C      write(*, '(1X,A,I3)')' NOPDIR=',NOPDIR
C     IF( NOPDIR .GT. 0)  
C    &   write(*, '(1X,A,A,A)')' ==>',RCOUTDIR(1:NOPDIR),'<=='
      OutNameLen = NOPDIR
      IF( NOPDIR .GT. 0) OutFileName = RCOUTDIR(1:NOPDIR)
C
C     write(*, '(1X,A,A,A)')' ==>',RootName(1:RootNameLen),'<=='
C 
      OutFileName(OutNameLen + 1:) = RootName(1:RootNameLen)
      OutNameLen = OutNameLen + RootNameLen
C
      IF( NOPEXT .GT. 0) OutFileName(OutNameLen + 1:LEN(OutFileName)) = 
     &   RCOUTEXT(1:NOPEXT)
C
      OutNameLen = OutNameLen + NOPEXT
C
C     TODO:  dump more gracefully here
      IF(OutNameLen .GT. 256) STOP 'Name 2 long'
C
      END
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SAVETOD()
C
C     This routine writes the date and time to the console and detailed
C     and failure logs.  Usually called just before ADIOS.
C
      IMPLICIT NONE
C
C      
C     Declarations regarding local variables or functions:
C
      CHARACTER TOD*8,DOD*9
C               These receive the current system date and time.
C
C
C     Declarations regarding COMMON block variables:
C
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C                  See descriptions in Block Data Genral
C
      CHARACTER CTMP*255, CWRK*80
C               These are just buffers used throughout program.
C
      CHARACTER LLT*1
C               This is the character that marks the end of a line sent
C               to the LogIt subroutine.
C
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
C
C
      CALL MyTime(TOD)
      CALL MyDate(DOD)
C
      WRITE(CTMP,2016)TOD,DOD,LLT
      CALL LogIt(LOGCDE,1,1)
C
 2016 FORMAT(' Time of death is ',A,', on ',A,'.  Your fortune:',A1)
C
      END
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LogIt(UFLAGS,BEFORE,AFTER)
C
C
C     This is the module that writes messages to all the log files.
C     The message to be logged must be written to common variable CTMP, 
C     and terminated with character LLT, before calling this subroutine.
C     The message is passed this way instead of as a parameter in hopes
C     of avoiding long, unnecessary automatic variables in the calling
C     routines, and minimizing function overhead while the arguments are
C     passed.  Of course, I (Matt) know squat about how FORTRAN actually
C     passes function arguments and manipulates the stack, and will 
C     cheerfully accept correction from anyone who understands such
C     efficiency issues.  (It kinda appears that arguments aren't copied
C     at all, but are actually passed by reference ALL the time, in which
C     case this arrangement probably doesn't offer much advantage.)
C
C     The first argument is an 8 element logical array that
C     tells which logs to write the message to.  Each position in this
C     string corresponds to a unit number in integer array USPEC.  These
C     are the log files, and should already be opened.  If the value in the
C     string is '0' (ASCII char decimal 48), then nothing happens.  If
C     anything else, the contents of CTMP get written to the corresponding
C     unit in USPEC.  The second argument specifies the number of blank 
C     lines to write before writing the contents of CTMP.  The third
C     arguments specifies the number of blank lines to write after CTMP.
C
      IMPLICIT NONE
C
C
C
C     Declarations regarding subroutine arguments:
C
      LOGICAL UFLAGS(8)
C               Flag block indicating which units to write message to.
C               Each column indicates an output unit.  A .TRUE. in the 
C               column means write the line to that unit.  The columns
C               correspond to the following units:
C                  col 1 - stdout (unit 6)
C                  col 2 - detailed log file (unit 90)
C                  col 3 - critical failure log (unit 91)
C                  col 4 - summary output file (unit 92)
C                  col 5 - Middle Lakes output file (unit 93)
C                  col 6 - Lake Superior output file (unit 94)
C                  col 7 - Lake Ontario output file (unit 95)
C                  col 8 - St. Lawrence output file (unit 96)
C      
      INTEGER BEFORE,AFTER
C
C      
C     Declarations regarding local variables or functions:
C
      INTEGER I,J,K,USPEC(8)
C
C      
C
C     Declarations regarding COMMON block variables:
C
      CHARACTER    CTMP*255, CWRK*80
C                  These are just buffers used throughout program.
C
      CHARACTER    LLT*1
C                  See description in INITPARAMS
C
C
      COMMON /WKCHAR/   CTMP,LLT,CWRK
C
C
C
      DATA    USPEC     / 6, 90, 91, 92, 93, 94, 95, 96/
C
      DO 90, J=1,8
C
         IF(UFLAGS(J))THEN
C
            DO 35, I=1, BEFORE
   35          WRITE( USPEC(J), 95, ERR=80, IOSTAT=K)
C
C           Find position of line-terminating character.  I had trouble with
C           garbage already in the static variable CTMP, so need a delimiter.
C       
            I = INDEX(CTMP,LLT)
C
            IF (I .GT. 1) THEN
C
C              Line is at least 1 character long.  Write it out.
C
               WRITE(USPEC(J),96, ERR=80, IOSTAT=K) CTMP(1:I-1)
C
            ELSEIF (I .EQ. 1) THEN
C
C              Must be a blank line.
C
               WRITE(USPEC(J),95, IOSTAT=K)
C
            ELSE
C
C              If CTMP always ends with LLT before each call to LogIt(),
C              then this won't happen.  But if we bungle the call, this 
C              will spit everything in CTMP to the log(s).
C
               WRITE(USPEC(J),96, ERR=80, IOSTAT=K) 
     &            ' Oops!  CTMP badly formed:'
               WRITE(USPEC(J),96, ERR=80, IOSTAT=K) CTMP
C
            ENDIF
C      
            DO 65, I=1, AFTER
   65          WRITE(USPEC(J),95, ERR=80, IOSTAT=K)
C
         ENDIF  
   90 CONTINUE
C
   95 FORMAT(' ')
   96 FORMAT(1X,A)
C
      RETURN
   80       CALL BuildInfo
            I = INDEX(CTMP,LLT)
            WRITE(*,'(/,A,I2,A,I5,A,/,A,A,/,A,A,/,/,A,/,/,A,A,/,/)')
     &         ' Uh-oh!  Could not write to unit ',USPEC(J),
     &         '.  For what its worth, IOSTAT=',K,',',
     &         ' although that code is only meaningful if you',
     &         ' have the documentation for',
     &         ' the compiler that produced this binary.  Build info',
     &         ' should be as follows:',CTMP(1:I),
     &         ' These errors are generally VERY BAD things, so this',
     &         ' run will abort now.'
            CALL ADIOS()
C
      END
C      
C      
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION LASTKC (STRING, STRLEN)
C
      IMPLICIT NONE
C
      INTEGER STRLEN, ASEQNO, I
      CHARACTER STRING*(*)
C
C     Start at end of string, and work backwards, looking for a regular
C     "keyboard" character (ASCII char between decimal 33 and 126, not
C     including the space character).
C
C     THIS ROUTINE RETURNS A VALUE OF 1 IF NO "KEYBOARD" CHARACTERS FOUND!
C     This behaviour is so it can used fearlessly as a string limit, such
C     as in SAMPLE(1:LASTKC(SAMPLE, 127)).  If returned 0 or some other
C     spurious value, then could have unpredictable results cropping up 
C     all over code.  So the moral is to always make sure that there is
C     at least one "keyboard" character in the string when you call this.
C
C     No doubt this could be a lot more efficient.  Maybe remove call
C     to ICHAR and just assume that the value in each byte is actually
C     the ASCII code.  But I don't think the CGLRRM will spend too much
C     of its time processing strings, so no big savings to had anyway.
C
      DO 10 I=STRLEN,1,-1
C        Get ASCII sequence number using ICHAR:
         ASEQNO=ICHAR(STRING(I:I))
   10    IF(ASEQNO .GE. 33 .AND. ASEQNO .LE. 126)GOTO 20
C
C     In case goes through loop, make sure don't return 0
      I=1
   20 LastKC = I
      END
C      
C-----------------------------------------------------------------------
C     
C
      SUBROUTINE ErrMess(MESSNO, EMETIC, LUNITS)
C
C
C     This module looks up a message (usually an error or warning) 
C     from the database and displays it.  The database is a simple text 
C     file with messages appearing in any order.  The message text is 
C     preceded and identified by a line of the form ('#',I4).  Columns 2-5 
C     hold a message id number between 1000 and 9999.  The messages are
C     stored in an external file so that the model can provide detailed
C     feedback and just go on and on without the ensuing blather taking
C     up memory in the executable.  It also allows the messages to be
C     changed/corrected without recompiling.  Even better, while we're
C     writing the darned model we just have to put a stub in the code and 
C     crank out a half-baked sentence or two that we can promise to finish
C     later.  ;-)
C
C
C      
      IMPLICIT      NONE
C
C
C
C     Declarations regarding subroutine arguments:
C
      INTEGER       MESSNO
C                   The message id number to look for and write out.      
C
      CHARACTER     EMETIC*16
C                   Name of the calling routine that initiates the message.
C      
      LOGICAL       LUNITS(8)
C                   The 8-element code array indicating which logs to 
C                   write to.
C
C
C     Declarations regarding local variables or functions:
C
      CHARACTER     CSTR*80
      INTEGER       I, M, L
C                   These are some scratch variables.
C
      CHARACTER*(*) CALTXT
C                   This is a message fragment used repeatedly.
C
      INTEGER       LASTKC
C                   Function returning the position of the last useful
C                   character in a string.
C
      INTEGER       LDBF
C                   Length of the name of the message database file.
C
C
C     Declarations regarding COMMON block variables:
C
      CHARACTER     DirSep*5,OSName(5)*8,MessDB*256
      CHARACTER*64  ModelerName
C                   See descriptions in Block Data Genral
C
      CHARACTER     CTMP*255, CWRK*80, LLT*1
C                   See descriptions in Block Data Genral
C
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /GENRLC/   OSName,DirSep,MessDB, ModelerName
C                                     
      PARAMETER (CALTXT=' The message originated from subroutine ')
C
      L = INDEX( EMETIC, LLT ) - 1
      IF(L .LT. 1) THEN
         L = 16
      ENDIF
C
C     
      LDBF = LastKC(MessDB, 256)
      IF(LDBF .LE. 1)THEN
C        No message database specified in input!
         WRITE(CTMP(1:126),'(A,A,A,I5,A1,A1)')
     &      ' You did not specify the location of the message ',
     &      'database.',LLT
         CALL LogIt(LUNITS,1,0)
      ELSE 
         OPEN(UNIT=99, FILE=MessDB(1:LDBF), STATUS='OLD', IOSTAT=I,
     &      ERR=6532)
C
C        If can't open the message database, squawk and return.
 6532    IF(I .NE. 0) THEN
C            WRITE(CTMP(1:126),'(A,A,A,I5,A1,A1)')
            WRITE(CTMP,'(A,A,A,I5,A1,A1)')
     &         ' Could not open message database file ',MessDB(1:LDBF),
     &         ': IOSTAT=',I,'.',LLT
            CALL LogIt(LUNITS,1,0)
            WRITE(CTMP,906)
     &         ' This code is specific to the compiler that produced', 
     &         ' your version of the model.',LLT
            CALL LogIt(LUNITS,0,0)
C           Ensure printing message below as well
            LDBF = 0
         ENDIF
      ENDIF
      IF(LDBF .LE. 1) THEN
         WRITE(CTMP,'(A,I4,A1,A1)')
     &      ' Please find the message database file and read #',MESSNO,
     &      '.',LLT
         CALL LogIt(LUNITS,0,0)
         WRITE(CTMP,908) CALTXT, EMETIC(1:L), '.', LLT
         CALL LogIt(LUNITS,0,0)
         WRITE(CTMP,906)
     &      ' HINT: the message database file and the model binary may',
     &      ' be',LLT
         CALL LogIt(LUNITS,0,0)
         WRITE(CTMP,907)' together in the same directory.',LLT
         CALL LogIt(LUNITS,0,1)
      RETURN
      ENDIF
C
      M=0
C           
C     Start loop to read in each line of database and look for #MESSNO.
C     Hope that there are less than 65533 lines in message database.
      DO 110 I= -32766, 32766, 1
C
C     Read a line from database into CSTR.  Jump to 111 if EOF.     
         READ(99,'(A5)',ERR=111)CSTR
C
C     If first char is a '#', read next four columns into integer M
         IF(CSTR(1:1) .EQ. '#')THEN
            READ(CSTR(2:5),909,ERR=109) M
         ENDIF
  109    CONTINUE
C
C     This is just some debugging stuff that might be handy later
C         WRITE(*,'(A,I7,A,A4,A,I4)') ' I=',I,'     CSTR(2:5)=',
C     &      CSTR(2:5),'   M=',M
C
C     If line has message number that we want, jump out of read loop
  110    IF(M .EQ. MESSNO) GOTO 111
C
C
C     If can't find the message number, do the best you can and abort.
C
  111 IF(M .NE. MESSNO) THEN
         IF(I .EQ. 32767) THEN      
            WRITE(CTMP,907)
     &      ' Oops!  Message database too long?  Tell someone.'
            CALL LogIt(LUNITS,1,0)
         ENDIF
         WRITE(CTMP,'(A,I4,A,A,A1,A1)') 
     &      ' Oops!  Could not find message number ',
     &      MESSNO,' in file ',MessDB(1:LDBF),'.',LLT
         CALL LogIt(LUNITS,1,0)
         WRITE(CTMP,'(A,I4,A,A1)')
     &      ' The message should have been preceded by a line with #'
     &      ,MESSNO,' in first 5 columns.',LLT
         CALL LogIt(LUNITS,0,0)
         WRITE(CTMP,908)CALTXT,EMETIC(1:L),'.',LLT
         CALL LogIt(LUNITS,0,0)
         WRITE(CTMP,907) 
     &      ' Pity.  You might never know what hit you.',LLT
         CALL LogIt(LUNITS,0,1)
         RETURN
      ENDIF   
C
C     So now we're at the appropriate message block.  Read each line
C     and write it back to stdout until EOF or find another message id.
C  
  905 FORMAT(A,I4,A,A,A1,A1)
      WRITE(CTMP,905)' ATTENTION!  Message #',MESSNO,
     &   ' from subroutine ', EMETIC(1:L), ':', LLT
      CALL LogIt(LUNITS,1,0)
      DO 220 I= -32766, 32766, 1
C           
         READ(99,'(A79)',ERR=221)CSTR
C
C        If first char is a '#', read next four columns into integer M
         IF(CSTR(1:1) .EQ. '#')THEN
C
C        If valid message number after '#', jump out of displaying 
C        loop, else display message.
C
            READ(CSTR(2:5),909,ERR=219) M
            IF(M .GE. 1000 .AND. M .LE. 9999) GOTO 221
         ENDIF
C
  219    WRITE(CTMP,'(1X,A,A1)')CSTR(1:LASTKC(CSTR,79)),LLT
  220    CALL LogIt(LUNITS,0,0)
            
  221 CLOSE(99)
C        
  906 FORMAT(A,A,A1)
  907 FORMAT(A,A1)
  908 FORMAT(A,A,A1,A1)
  909 FORMAT(I4)
      RETURN
      END
C
C
C-----------------------------------------------------------------------
C         
      SUBROUTINE ADIOS()
C
C
C     This module exits the program.
C
C
      CHARACTER STIME*8,SDATE*9, MSG*2024, LT*1, AP*1
      CHARACTER*(*) C1, C2, C3, C4, C5, C6,C7,C8,C9,C10,C11,C12,C13,C14,
     &    C15,C16,C17,C18,C19,C20,C21
      INTEGER NUMMSG, I, J, ST
C
      COMMON /STIME/    STIME,SDATE
C
      PARAMETER (LT=CHAR(0))
      PARAMETER (AP=CHAR(39))
C
      PARAMETER (C1=' Some mornings it just doesn'//AP//
     &   't seem worth it to gnaw through the leather straps!'//LT)
      PARAMETER (C2=' Did you really expect it to be easy?'//LT)
      PARAMETER (C3=
     &   ' It is better to travel happily than to arrive.'//LT)
      PARAMETER (C4=' He'//AP//'s dead Jim. Grab his tricorder. I'
     &   //AP//'ll get his wallet.'//LT)
      PARAMETER (C5=' Game Over.'//LT)
      PARAMETER (C6=
     &   ' The floggings shall continue until morale improves.'//LT)
      PARAMETER (C7= ' Complex problems have simple,'//
     &   ' easy-to-understand wrong answers.'//LT)
      PARAMETER (C8= ' Inside every small problem are '//
     &   'large problems struggling to get out.'//LT)
      PARAMETER (C9=' Quailty is Job 1.1'//LT)
      PARAMETER (C10= ' If a project is not worth doing at all, it'
     &   //AP//'s not worth doing well.'//LT)
      PARAMETER (C11=' Politicians are like diapers; they need '//
     &   'changing often, for the same reasons.'//LT)
      PARAMETER (C12=' C Code.  C Code Run.  Run, Code, '//
     &   'RUN! PLEASE!!!!'//LT)
      PARAMETER (C13=' Luke: I don'//AP//'t... believe '//
     &   'it.            Yoda: That is why you fail.'//LT)
      PARAMETER (C14=' Badgers?  Badgers?!?  We don'//AP//'t need no '//
     &   'stinkin'//AP//' Badgers!!!'//LT)
      PARAMETER (C15=' One person'//AP//'s error is another person'//
     &   AP//'s data.'//LT) 
      PARAMETER (C16=' Hardware: The parts of a computer '//
     &   'system that can be kicked.'//LT) 
      PARAMETER (C17=
     &   ' Disinformation is worse than datinformation.'//LT)
      PARAMETER (C18=' As a computer, I find your faith in '//
     &   'technology amusing.'//LT)  
      PARAMETER (C19=' A young person knows the rules but the old '//
     &   'person knows the exceptions.'//LT)
      PARAMETER (C20=' Don'//AP//'t sweat the petty things, and don'//
     &   AP//'t pet the sweaty things.'//LT)  
      PARAMETER (C21=' If you eat a live toad each morning, nothing '//
     &   'worse will happen that day.'//LT)
C
      DATA NUMMSG /21/
C
      MSG = C1//C2//C3//C4//C5//C6//C7//C8//C9//C10//C11//C12//C13//
     &   C14//C15//C16//C17//C18//C19//C20//C21
C
C
   66 FORMAT(I2)
      READ(STIME(7:8), 66, ERR=88)J
      I=J
      READ(STIME(4:5), 66, ERR=88)J
      I=J+I
      READ(STIME(1:2), 66, ERR=88)J
      I=J+I
   88 NUMMSG=MOD(I,NUMMSG)+1
      
C     WRITE(*,'(A,I3,A,I3,A,I3)')' I=',I,'  J=',J,'  NUMMSG=',NUMMSG
      J=0
      DO 93 I=1,NUMMSG
         ST=J+1
         DO 91 J=ST,32766
   91       IF( MSG(J:J) .EQ. LT ) GOTO 93
   93    CONTINUE
C
      WRITE(*,'(1X,A,/)')MSG(ST:J-1)
      WRITE(*,'(1X,A)')' Run failed'
C
      STOP ' -1'
      END
C
C
C-------------------------------------------------------------------------
C
C
      SUBROUTINE GETQTRMONBOUNDS (YEAR, MONTH, QTR)
C
C     GetQtrMonthBounds determines the last day of each quarter-month based on    
C     the definitions received from USACE - Buffalo.                              
C
      IMPLICIT   NONE
C
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER      YEAR
C                  (Input) The year of interest.
C
      INTEGER      MONTH
C                  (Input) The month of interest.
C
      INTEGER      QTR(4)
C                  (Output) Array containing the ending day of each quarter-month
C                  in MONTH and YEAR. 
C
C                  
C     Don't need number of days in year, so just use QTR(1) as scratch for now:
      CALL NUMDAYS(YEAR, MONTH, QTR(1), QTR(4))
C
      IF (QTR(4) .LE. 29 ) THEN
         QTR(1) = 7
         QTR(2) = 14
         QTR(3) = 21
      ELSE
         QTR(1) = 8
         QTR(2) = 15
         QTR(3) = 23
      ENDIF
C
      END
C
C-------------------------------------------------------------------------
C
C
      INTEGER FUNCTION GetDaysInQMon(Year, Month, Qmonth)
C
C     Determines length of the specified quarter-month. 
C
      IMPLICIT    NONE
C
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER      Year, Month, Qmonth
C                  (Input) Identify the quarter-month of interest.
C
C
C     Declarations regarding local variables:
C
      INTEGER      QTRBNDS(4)
C                  Just scratch, holding the ending day of the month 
C                  for each quarter-month
C 
C
      CALL GetQtrMonBounds (Year, Month, QTRBNDS)
C
      IF ( Qmonth .EQ. 1 ) THEN
         GetDaysInQMon = QtrBnds(Qmonth)
      ELSE
         GetDaysInQMon = QtrBnds(Qmonth) - QtrBnds(Qmonth - 1)
      ENDIF
C
      END
C
C-------------------------------------------------------------------------
C
C
      INTEGER FUNCTION GETQMONNUM (YEAR, MONTH, DAY)
C
C     Determines which quarter-month the specified date is in. 
C
      IMPLICIT    NONE
C
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER      YEAR
C                  (Input) The year of interest.
C
      INTEGER      MONTH
C                  (Input) The month of interest.
C
      INTEGER      DAY
C                  (Input) The day of interest.
C
C
C     Declarations regarding local variables:
C
      INTEGER      QTRBNDS(4)
C                  Just scratch, holding the ending day of the month 
C                  for each quarter-month
C 
C
      CALL GETQTRMONBOUNDS (YEAR, MONTH, QTRBNDS)
C
      DO 33, GETQMONNUM = 1, 4
C
   33    IF( DAY .LE. QTRBNDS(GETQMONNUM) ) RETURN
C
C
      STOP 'Invalid DAY in GETQMONNUM'
C
      END
C
C-------------------------------------------------------------------------
C
C
C
      LOGICAL FUNCTION LEAPYEAR (YEAR)
C
C     Returns .TRUE. if the specified year is a leap year
C
      IMPLICIT   NONE
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER      YEAR
C                  (Input) The year of interest.
C
      LEAPYEAR = .FALSE.
C
      IF(  MOD(YEAR, 4) .EQ. 0 )THEN
C        Divisible by four, so might be a leap year
C
         IF( MOD(YEAR, 100) .NE. 0 .OR. MOD(YEAR, 400) .EQ. 0 ) THEN
C           It's either not a century, or divisible by 400, so must be a leap year
            LEAPYEAR = .TRUE.
         ENDIF
      ENDIF
C
      END
C
C
C
C-------------------------------------------------------------------------
C
C
      SUBROUTINE NUMDAYS(YEAR,MONTH,DAYSINYEAR,DAYSINMONTH)
C
C     Computes the number of days in the specified month and in the year.
C
      IMPLICIT   NONE
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER      YEAR
C                  (Input) The year of interest.
C
      INTEGER      MONTH
C                  (Input) The month of interest.
C
      INTEGER      DAYSINYEAR
C                  (Output) The number of days in YEAR.
C
      INTEGER      DAYSINMONTH
C                  (Output) The number of days in MONTH.
C
C     Declarations regarding local variables or functions:
C
      LOGICAL      LEAPYEAR
C                  Returns .TRUE. if a leap year
C
      INTEGER      DAYSEACHMONTH(12) 
C                  Contains the number of days in each month
C
C
      DATA DAYSEACHMONTH /31,28,31,30,31,30,31,31,30,31,30,31/
C
C     The following should never happen if called properly.
C     If it happens, logic is bad in calling routine.
      IF(MONTH .LT. 1 .OR. MONTH .GT. 12) STOP 'Ren - get help!'
C
C
      DAYSINYEAR = 365
      DAYSINMONTH = DAYSEACHMONTH(MONTH)
C
      IF( LEAPYEAR(YEAR) ) THEN
C         
         DAYSINYEAR = DAYSINYEAR + 1
C
         IF( MONTH .EQ. 2 ) THEN
            DAYSINMONTH = DAYSINMONTH + 1
         ENDIF
      ENDIF
C
      END
C
C
C-------------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
C
      INTEGER FUNCTION DateSequence (DAY, MONTH, YEAR)
C
C     This routine computes the number of the date, with days numbered
C     sequentially from 1 (corresponding to the base date).  the base
C     date must be 1 March on any LEAP YEAR or century boundary. Dates
C     must be between the base date and the base date plus 2,147,483,646.
C     This routine is set up to use only 4-byte integer variables and 
C     to avoid intermediate use of reals in the evaluation of expressions.
C     It is a modification of the previous version which used 2-byte integer
C     variables. For reference on the 2-byte version of this subroutine
C     see Dr. Dobb's Journal 80:66-70.  Adapted by Tim Hunter.
C
C     Note that this routine correctly handles the year 2000, as well as all
C     other leap year issues.
C
      IMPLICIT NONE
C
      INTEGER   BASEYEAR
      PARAMETER (BASEYEAR = 0)
C
      INTEGER   DAY, MONTH, YEAR
C
      INTEGER   NUMYEARS, TOTC, LEAPC, SKIPC, NONLEAPC
C
      IF (YEAR .LE. BASEYEAR) STOP 'YEAR .LE. BASEYEAR'
C
      NUMYEARS = YEAR - BASEYEAR - 1 + (MONTH + 9) / 12
      TOTC = (BASEYEAR + NUMYEARS)/100 - (BASEYEAR/100)
      SKIPC = 3 - MOD(BASEYEAR/100, 4)
      LEAPC = (TOTC + 3 - SKIPC) / 4
      NONLEAPC = TOTC - LEAPC
C
      DateSequence = (365 * NUMYEARS + (NUMYEARS/4) - NONLEAPC
     &                + (153 * MOD(MONTH + 9, 12) + 2) / 5 + DAY) 
      END
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION ENDSTR (STRING, STRLEN)
C
      IMPLICIT NONE
C
      INTEGER STRLEN, ASEQNO, I
      CHARACTER STRING*(*)
C      DIMENSION STRING (1,*)
C
C     Search from beginning of string looking for a regular "keyboard"  
C     character (ASCII char between decimal 33 and 126, not including 
C     the space character).  Return position of last "keyboard"
C     character before whitespace or control characters.
C
C     If the string contains no "keyboard" characters, 0 is returned.
C     DON'T USE THE RESULT OF THIS FUNCTION IN A STRING LIMIT UNLESS
C     YOU POSITIVELY KNOW THAT THE STRING CONTAINS A "KEYBOARD" CHARACTER!
C     For instance, SAMPLE(ENDSTR(SAMPLE, 127):) will cause unexpected
C     results if SAMPLE is all blanks or something.  So safe use of this
C     for parsing lines read in from data files would be two-step, as in:
C          I= ENDSTR(SAMPLE, 127)
C          IF(I .NE. 0)THEN
C              WRITE(90,2015)SAMPLE(1:I)
C          ENDIF
C
C     No doubt this could be a lot more efficient.  Maybe remove call
C     to ICHAR and just assume that the value in each byte is actually
C     the ASCII code.  But I don't think the CGLRRM will spend too much
C     of its time processing strings, so no big savings to had anyway.
C
      DO 10 I=1,STRLEN,1
C        Get ASCII sequence number using ICHAR:
         ASEQNO=ICHAR(STRING(I:I))
C        write(*,'(A,I2,A,A1)')' ENDSTR=',ENDSTR,'   ',STRING(1,ENDSTR)
   10    IF(ASEQNO .LT. 33 .OR. ASEQNO .GT. 126)GOTO 20
C
   20 IF(I.GT. STRLEN .OR. ASEQNO .LT. 33 .OR. ASEQNO .GT. 126)THEN
         I = I - 1
      ENDIF

      ENDSTR = I
      END
C
C-----------------------------------------------------------------------
C      
C      
      CHARACTER*16 FUNCTION MakeIntervalName( IntervalNumber, 
     &   Capitalize, Length)
C
C     Stupid little function that returns a string like "daily" or
C     "Daily" depending on the arguments.  The string is terminated by
C     a space character, for use with function ENDSTR().  The strings
C     returned by this function differ from the strings in array
C     IntervalName (in COMMON block RCRUNC) only by spelling out
C     "quarter-monthly" and optionally capitalizing.
C
      IMPLICIT NONE 
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER      IntervalNumber 
C                  (Input) A number from 1-4 designating the routing 
C                  interval of interest.  1=daily, 2=weekly, 
C                  3=qmonthly, 4=monthly.
C      
      LOGICAL      Capitalize 
C                  (Input) Tells whether to capitalize the first letter
C                  (and ninth letter for Quarter-Monthly). 
C      
      INTEGER      Length 
C                  (Output) The number of characters to use in the
C                  string being returned.
C      
C     Declarations regarding local variables:
C
      CHARACTER*26 IntervalNames
C                  Contains the text strings 'DAILY' 'MONTHLY' etc.  
C                  Quarter-monthly string runs from characters 12-26; 
C                  monthly string runs from characters 20-26; weekly
C                  runs from 6-11; daily string runs from 1-5.
C
      CHARACTER*16 TempString

      DATA IntervalNames /'dailyweeklyquarter-monthly'/
C
      IF ( IntervalNumber .EQ. 1 ) THEN
         TempString = IntervalNames(1:5)//' '
         Length = 5
C
      ELSE IF ( IntervalNumber .EQ. 2 ) THEN
         TempString = IntervalNames(6:11)//' '
         Length = 6
C
      ELSE IF ( IntervalNumber .EQ. 3 ) THEN
         TempString = IntervalNames(12:26)//' '
         IF(Capitalize) TempString(9:9) = 'M'
         Length = 15
C
      ELSE IF ( IntervalNumber .EQ. 4 ) THEN
         TempString = IntervalNames(20:26)//' '
         Length = 7
C
      ELSE
         WRITE(*,'(A, I8)') 
     &      ' Bad IntervalNumber in MakeIntervalName():', 
     &      IntervalNumber
         STOP 'The Brain'
      ENDIF
C
      IF(Capitalize) TempString(1:1) = 
     &   CHAR(ICHAR( TempString(1:1) ) - 32)
C
      MakeIntervalName = TempString

      END
C
C     
C----------------------------------------------------------------------
C
      INTEGER FUNCTION TSOutputFileName(TSNum, DataType, Lake, FilePath)
C
C     This routine determines the name to use for an output, and stores
C     it in FilePath.  Returns lenght of filespec.
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER       Lake, DataType
C                   (Input) Lake number corresponding to element in
C                   LakesToSolve, and type of data (1=level, 2=flow,
C                   3=time series data summary)
C
      INTEGER       TSNum
C                   (Input) Timeseries of output (1=daily, 2=weekly,
C                   3=qmonthly, 4=monthly)
C
      CHARACTER*256 FilePath
C                   (Output) Holds the filename.
C
C     Declarations regarding local variables or functions:
C
      CHARACTER*2   LakeID(2,8)
      CHARACTER*3   TypeID(3)
      CHARACTER     IntervalID(4)
C                   These hold the filename parts.
      INTEGER       FNameLen, DTI
C
      CHARACTER*256 OutFileName
C                   Function tha builds complete filespec.
C     
C                   COMMON variables
C
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      CHARACTER*21  LAKENAME(16), RiverName(16)
      INTEGER       NLN(16), LkToSlv2LN(8), NRN(16), LkSlv2RtLoopI(8)
      INTEGER       GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &              SupVerbosity, MLVerbosity, OntVerbosity
      LOGICAL       WriteTSOutput(4,2,8), WriteErieOutflow
      CHARACTER     CTMP*255, CWRK*80
      CHARACTER     LLT*1
C
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /LNAMEC/   LAKENAME, RiverName
      COMMON /LNAMEI/   NLN, LkToSlv2LN, NRN, LkSlv2RtLoopI
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /OUTSTF/   WriteTSOutput, WriteErieOutflow
C
      DATA LakeID / 'sp', 'sm', 'mh', 'cr', 'sc', 'dr', 'er', 'nr', 
     &   'on', 'sl', '??', '??', '??', '??', '??', '??' /
      DATA IntervalID / 'd', 'w', 'q', 'm' /
      DATA TypeID / 'mlv', 'flo', 'tsd' /
C
C     Use the lake abbreviations (sp, mh, etc) for tsd filenames.
C     DTI is 1nd index to LakeID
      IF( DataType .EQ. 2 ) THEN
         IF( Lake .EQ. 4 .AND. WriteErieOutflow ) THEN
C           Instead of Niagara River flow we want Erie outflow
C           Use lake name, not Niagara RIver name
            DTI = 1
         ELSE
            DTI = 2
         ENDIF
      ELSE
         DTI = 1
      ENDIF
C     
C
C     Build the complete filespec
C
      FilePath = OutFileName( 
     &   LakeID(DTI, Lake)//IntervalID(TSNum)//TypeID(DataType), 
     &   6, FNameLen)

      IF (GenVerbosity .GE. 5)THEN
         WRITE(CTMP, 5689) 
     &      LakeName(LkToSlv2LN(Lake))(1:NLN(LkToSlv2LN(Lake))), 
     &      FilePath(1:FNameLen), LLT
         CALL LogIt(LOGD,1,0)
      ENDIF

      TSOutputFileName = FNameLen

C
 5689 FORMAT(' Formulating output filename for ', A, ': ', A, A1)

      END
C-----------------------------------------------------------------------
C
      LOGICAL FUNCTION InitTSOutputFile( TSNum, DataType, Lake, 
     &   OriginalInterval)
C
C     This routine determines the name to use for an output file and
C     opens it, then writes the headers with information as specified
C     by the AdHoc people.  If it can't be opened, returns .FALSE.
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER       Lake, DataType
C                   (Input) Lake number corresponding to element in
C                   LakesToSolve, and type of data (1=level, 2=flow).
C
      INTEGER       TSNum, OriginalInterval
C                   (Input) Timeseries of output (1=daily, 2=weekly,
C                   3=qmonthly, 4=monthly), and original routing 
C                   interval that output is derived from (1-4, as
C                   above).   
C
C     Declarations regarding local variables or functions:
C
      INTEGER       I, J, F
C                   scratch
C
      INTEGER       Len
C                   Lenght of the output filename.
C
      CHARACTER*256 TmpStr
C                   Scratch
C
      INTEGER       TSOutputFileName
C                   Formulates filename and path for output files.
C
      CHARACTER*16  MakeIntervalName
C                   Returns a string with the interval name.  For output.
C 
      INTEGER       LASTKC
C                   Function returning the position of the last useful
C                   character in a string.
C
C
C
C     Declarations regarding COMMON block variables:
C
      INTEGER       TSOutputFileNum(4,2,8)
C                   The filespec numbers reserved for each type of 
C                   time-series output data.
C                   See descriptions in Block Data Genral
C
      CHARACTER*8   IntervalName(4)
      INTEGER       IntervalNameLen(4)
C                   See descriptions in Block Data Genral
C
      INTEGER       NOPDIR, NOPEXT
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      CHARACTER*21  LAKENAME(16), RiverName(16)
      INTEGER       NLN(16), LkToSlv2LN(8), NRN(16), LkSlv2RtLoopI(8)
      INTEGER       GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &              SupVerbosity, MLVerbosity, OntVerbosity
      CHARACTER     CTMP*255, CWRK*80
      CHARACTER     LLT*1
C
      LOGICAL       OnlyOneOpen
      INTEGER       StartYear, StartMonth, StartDay, 
     &              EndYear, EndMonth, EndDay
      CHARACTER*3   MonthAbbrev(12)
      CHARACTER     STIME*8,SDATE*9
      CHARACTER     DirSep*5,OSName(5)*8,MessDB*256
      CHARACTER*64  ModelerName
C                   See descriptions in Block Data Genral
C
      CHARACTER*80  RCTITLE(3), MLHEADER
C                   See descriptions in PARLINE()
C
C     These variables not used this function; declared because members
C     of same COMMON block as IntervalNameLen and IntervalName
      LOGICAL       LakesToSolve(8)
      CHARACTER*8   RCRUNTYP
      INTEGER       MonthLength
      CHARACTER*8   MonthLengthName(2)
      CHARACTER*9   MLSolutionName(2)
      CHARACTER*16  SMSolutionName(4), USLSolutionName(5)
      LOGICAL       WriteTSOutput(4,2,8), WriteErieOutflow
C
C
C
C
      COMMON /MONAME/   MonthAbbrev
      COMMON /OUTSTF/   WriteTSOutput, WriteErieOutflow
      COMMON /GENRLC/   OSName,DirSep,MessDB, ModelerName
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /LNAMEC/   LAKENAME, RiverName
      COMMON /LNAMEI/   NLN, LkToSlv2LN, NRN, LkSlv2RtLoopI
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /OUTNMI/   NOPDIR, NOPEXT, TSOutputFileNum, OnlyOneOpen
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /RCRUNI/   LakesToSolve, 
     &                  IntervalNameLen, MonthLength,
     &                  StartYear, StartMonth, StartDay,
     &                  EndYear, EndMonth, EndDay
      COMMON /RCRUNC/   RCRUNTYP, IntervalName, MonthLengthName,
     &                  MLSolutionName, SMSolutionName, USLSolutionName
      COMMON /TITLES/   RCTITLE, MLHEADER
      COMMON /STIME/    STIME,SDATE
C
C
      InitTSOutputFile = .FALSE.
C-----------------------------------------------------------------------
C     FYI, calling routine sends CTMP to LOGCDE on returning .FALSE.
C-----------------------------------------------------------------------
C
C     Build the complete filespec
      Len = TSOutputFileName(TSNum, DataType, Lake, TmpStr)
      IF(GenVerbosity .GE. 7) THEN
         WRITE(CTMP, 5689) 
     &      LakeName(LkToSlv2LN(Lake))(1:NLN(LkToSlv2LN(Lake))), 
     &      TSOutputFileNum(TSNum,DataType,Lake), TmpStr(1:Len), LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
C  Open the output file.
C
      OPEN( UNIT = TSOutputFileNum(TSNum,DataType,Lake),
     &   FILE = TmpStr(1:Len), IOSTAT=I) 
      IF ( I .NE. 0 ) THEN
C        There was an error opening file.  Squawk and punt.
C        If problem is too many files open, then won't be able to open
C        message database either, so close the open output files first.
         F = 0
         DO 723, J=10,89
C           If it doesn't backspace, not already open 
            BACKSPACE(J, ERR=723) 
            CLOSE(J, ERR=723)
C           Count the ones already open
            F = F + 1
C            WRITE(CTMP(3*F:3*F+2),'(I3)') J
  723       CONTINUE
C        Files 90-96 should still be open.
C         WRITE(CTMP(3*F+3:3*F+3),'(A1)') LLT
C         CALL LogIt(LOGCDE,0,1)
         F = F + 7
C
         IF(OnlyOneOpen)THEN 
C           Problem could be anything; issue general message.
            J = 1124
         ELSE  
C           *Very* good chance that simply too many files open.  
            J = 1125
         ENDIF
C
         CALL ErrMess(J,'InitTSOutputFile',LOGCDE)
         WRITE(CTMP,2015) I, CWRK(1:Len), LLT
         CALL LogIt(LOGCDE,1,1)
         WRITE(CTMP,2013) F, LLT
         CALL LogIt(LOGCDE,0,1)
         CALL BuildInfo
         CALL LogIt(LOGCE,0,0)
         CALL SAVETOD
         CALL ADIOS()
      ENDIF
C
C 
C     Write name, data interval, and data type as header comment,
C     then the "Units:" line
C
      CWRK = MakeIntervalName( TSNum, .TRUE., F)
      IF(DataType .EQ. 1) THEN
         WRITE( TSOutputFileNum(TSNum,DataType,Lake), 4064)
     &      LakeName(LkToSlv2LN(Lake))(1:NLN(LkToSlv2LN(Lake))), 
     &      CWRK(1:F), ' Levels'
         WRITE(TSOutputFileNum(TSNum,DataType,Lake),4069) 'meters'
      ELSE
         IF(Lake .EQ. 4 .AND. WriteErieOutflow)THEN
C           Clearly label the Erie outflow file if applicable
            WRITE( TSOutputFileNum(TSNum,DataType,Lake), 4064)
     &         'Lake Erie', CWRK(1:F), ' Outflows'
            WRITE(TSOutputFileNum(TSNum,DataType,Lake),4065)
     &      'The following values are sums of flows from the '//
     &      'Niagara River and Welland Canal'
         ELSE 
            WRITE( TSOutputFileNum(TSNum,DataType,Lake), 4064)
     &         RiverName(LkToSlv2LN(Lake))(1:NRN(LkToSlv2LN(Lake))), 
     &         CWRK(1:F), ' Flows'
         ENDIF
         WRITE(TSOutputFileNum(TSNum,DataType,Lake),4069) 'm3s'
      ENDIF
C 
C     Write the "Interval:" line
C
      WRITE(TSOutputFileNum(TSNum,DataType,Lake),4063) CWRK(1:F)
C
C     Write the title lines
C
      DO 733, F=1,3
  733     WRITE(TSOutputFileNum(TSNum,DataType,Lake),4065)
     &       RCTITLE(F)(1:LASTKC(RCTITLE(F),80))
C
C     If middle lake, write the MLHEADER line
C
      IF(Lake .GE. 2 .AND. Lake .LE. 6) 
     &   WRITE( TSOutputFileNum(TSNum,DataType,Lake),4065)
     &       MLHEADER(1:LASTKC(MLHEADER,80))
C
C     Write line with start/end simulation dates
C
      WRITE( TSOutputFileNum(TSNum,DataType,Lake), 4068)
     &   StartDay, MonthAbbrev(StartMonth), StartYear,
     &   EndDay, MonthAbbrev(EndMonth), EndYear
C 
C     Write line specifying which routing period used in simulation
C
      CWRK = MakeIntervalName( OriginalInterval, .FALSE., F)
      IF(OriginalInterval .EQ. TSNum) THEN
         WRITE( TSOutputFileNum(TSNum,DataType,Lake), 4067)
     &      'resulting directly from ', CWRK(1:F) 
      ELSE
         WRITE( TSOutputFileNum(TSNum,DataType,Lake), 4067)
     &      'averaged/converted from ', CWRK(1:F) 
      ENDIF
C 
C     Write line with modeler's name, date of run, and model version
C
      WRITE(TSOutputFileNum(TSNum,DataType,Lake),4066)
     &   ModelerName(1:LASTKC(ModelerName,64)), STIME, SDATE
C
      CALL BuildInfo
      WRITE(TSOutputFileNum(TSNum,DataType,Lake),4065)
     &   CTMP(2:INDEX(CTMP, LLT)-1)
      WRITE(TSOutputFileNum(TSNum,DataType,Lake),4065) ' '
C
C     If opening and closing for each write, close it up.
C    
      IF(OnlyOneOpen) CLOSE(TSOutputFileNum(TSNum,DataType,Lake))
C
C     Be positive when return to calling function 
C    
      InitTSOutputFile = .TRUE.
C
 2015 FORMAT(' Got IOSTAT=',I4,' when opening file ',A,' for write.',A1)
 2013 FORMAT(' Approximately ',I2,' files were already open.',A1)
 4063 FORMAT('Interval: ', A, /, '#')
 4064 FORMAT('#', /, '# ', A, ' Mean ', A, A)
 4065 FORMAT('# ', A)
 4066 FORMAT('#', /, '# Modeled by ', A, ' at ', A, ' on ', A, /, '#')
 4067 FORMAT('#', /, '# Values ', A, A, ' routing computations.')
 4068 FORMAT('#', /, '# Data generated from simulation starting on ', 
     &   I2, A3, I5, ', and ending on ', I2, A3, I5)
 4069 FORMAT('#', /, 'Units: ', A)
C
 5689 FORMAT(' Initializing ',A,' output file on unit ',I2, ': ', A, A1)
      END
C
C
C-----------------------------------------------------------------------
C
C
C

      SUBROUTINE WriteToDisk(Interval, Lake )
C
C     Writes level and flow in output buffers to disk.
C
      IMPLICIT NONE 
C
C     Declarations regarding function arguments:
C
      INTEGER       Lake, Interval
C                   (Input) Lake number corresponding to element in
C                   LakesToSolve, and routing interval number.
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER       I, J, K
C                   Scratch
C
      INTEGER       ErrNum
C                   Which message to use if failure
C
      INTEGER       TSOutputFileName
C                   Formulates filename and path for output files.
C
      CHARACTER*256 TmpStr
C                   Scratch
C
C     Declarations regarding COMMON block variables:
C
      CHARACTER*101 MonthOutputBuf(10,2,8)
      CHARACTER*103 QMonthOutputBuf(4,10,2,8)
      CHARACTER*19  WeekOutputBuf(52,2,8)
      CHARACTER*74  DailyOutputBuf(48,2,8)
C                   See descriptions in TallyMeans
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
      CHARACTER     CTMP*255, CWRK*80
      CHARACTER     LLT*1
      INTEGER       NOPDIR, NOPEXT, TSOutputFileNum(4,2,8)
      LOGICAL       OnlyOneOpen
      LOGICAL       WriteTSOutput(4,2,8), WriteErieOutflow
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /OUTTXT/   MonthOutputBuf, WeekOutputBuf
      COMMON /OUTTX2/   DailyOutputBuf
      COMMON /OUTTX3/   QMonthOutputBuf
      COMMON /OUTNMI/   NOPDIR, NOPEXT, TSOutputFileNum, OnlyOneOpen
      COMMON /OUTSTF/   WriteTSOutput, WriteErieOutflow
C
C 5659 FORMAT(' Writing output to unit ',I2, A1)
 5689 FORMAT(' Writing output to file {', A, '} on unit ',I2, A1)
 8653 FORMAT(A)
 9015 FORMAT(' Got IOSTAT=',I5,' when accessing file ==>',A,'<==.',A1)
C
      DO 596, J=1,2
C
C        See if we want the outputs . . .  
         IF(WriteTSOutput(Interval,J,Lake)) THEN
C
C           Get filename in TmpStr
            K = TSOutputFileName(Interval, J, Lake, TmpStr)
C
C           If opening and closing for each write, open it up.
            IF(OnlyOneOpen) THEN
               ErrNum = 1126
               OPEN( UNIT= TSOutputFileNum(Interval,J,Lake), ERR=9455,
     &           STATUS='OLD',FILE=TmpStr(1:K),IOSTAT=I,ACCESS='APPEND')
            ENDIF
C
            IF(GenVerbosity .GE. 3)THEN
               WRITE(CTMP, 5689) TmpStr(1:K), 
     &            TSOutputFileNum(Interval,J,Lake), LLT 
               CALL LogIt(LOGD,1,0)
            ENDIF
C
            ErrNum = 1127
C
            IF( Interval .EQ. 4) THEN
               DO 444, I = 1, 10
  444             IF( MonthOutputBuf(I,J,Lake)(1:1) .NE. '#') 
     &               WRITE( TSOutputFileNum(Interval,J,Lake), 8653,
     &                  ERR=9455) MonthOutputBuf(I, J, Lake)
C
            ELSEIF( Interval .EQ. 3) THEN
               DO 333, I = 1, 10
                  DO 333, K = 1, 4
  333             IF( QMonthOutputBuf(K,I,J,Lake)(1:1) .NE. '#') 
     &               WRITE( TSOutputFileNum(Interval,J,Lake), 8653,
     &                  ERR=9455) QMonthOutputBuf(K, I, J, Lake)
C
            ELSEIF( Interval .EQ. 2) THEN
C               DO 222, I = 1, 11
               DO 222, I = 1, 52
  222             IF( WeekOutputBuf(I,J,Lake)(1:1) .NE. '#') 
     &               WRITE( TSOutputFileNum(Interval,J,Lake), 8653,
     &                  ERR=9455) WeekOutputBuf(I, J, Lake)
C
            ELSE
               DO 111, I = 1, 48
  111             IF( DailyOutputBuf(I,J,Lake)(1:1) .NE. '#') 
     &               WRITE( TSOutputFileNum(Interval,J,Lake), 8653,
     &                  ERR=9455) DailyOutputBuf(I, J, Lake)
            ENDIF
C
C           If opening and closing for each write, close it now.
            IF(OnlyOneOpen) 
     &         CLOSE(TSOutputFileNum( Interval, J, Lake))
         ENDIF
  596    CONTINUE
C    
C     Exit normally
      RETURN
C    
C     Puke
 9455 CALL ErrMess(ErrNum,'WriteToDisk'//LLT//'    ',LOGCDE)
      WRITE(CTMP,9015) I, TmpStr(1:K), LLT
      CALL LogIt(LOGCDE,1,1)
      CALL BuildInfo
      CALL LogIt(LOGCE,0,0)
      CALL SAVETOD
      CALL ADIOS()
      END
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE WriteHighPrecision(LakeNum, MeanLev, EOPLev, MeanFlow)
C
C     This routine sends a specialy formatted string to the detailed log
C     containing the full double-precision values for mean and EOP lake
C     level, and mean outflow.  A tag identfies whether the period ends 
C     a month or a quarter-month, in order to easily distinguish the values 
C     using a filter program like grep.  The date information is read from
C     the PeriodEnd arrays in the COMMON block.
C 
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER          LakeNum
C               (Input)  Index for that lake into PeriodEnd
C
      DOUBLE PRECISION MeanLev, EOPLev, MeanFlow
C               (Input)  The results to print
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER       I, L
C                   Indicies to the routing interval and routing loop arrays
C
      CHARACTER*1   EOPType
C                	  Returns 'M' if period ends a month, or 'Q' if period
C                   ends a quarter-month, else returns 'W', or 'D'.
C
c     CHARACTER*1   PECODE
C               Code for whether period ends a day (or week), qm, or month
C 
c     INTEGER      QTRBNDS(4), D
C                  Holds the ending day of the month 
C                  for each quarter-month, and a scratch counter
C 
C
C     Declarations regarding COMMON block variables:
C
      CHARACTER*8   IntervalName(4)
      INTEGER       PeriodEnd(6,3)
      INTEGER       MonthLength, IntervalNameLen(4)
      INTEGER       NLN(16), LkToSlv2LN(8), NRN(16), LkSlv2RtLoopI(8)
      CHARACTER*21  LAKENAME(16), RiverName(16)
      CHARACTER     CTMP*255, CWRK*80, LLT*1
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C
C                   Declared because of other variables in COMMON block:
      DOUBLE PRECISION dT(3)
      CHARACTER*8   RCRUNTYP
      LOGICAL       NumIncsVary(3) 
      INTEGER       Day(3), Q(4,3), DaysSinceSaturday
      INTEGER       MidLakeSolution, StMarysSolution, UStLawrSolution 
      INTEGER       BOMRound, FlowRound, LevelRound, P77ACompat
      INTEGER       RoutingInterval(3), RoutingIncrement(3)
      LOGICAL       LakesToSolve(8)
      INTEGER       StartYear, StartMonth, StartDay, 
     &              EndYear, EndMonth, EndDay
      CHARACTER*8   MonthLengthName(2)
      CHARACTER*9   MLSolutionName(2)
      CHARACTER*16  SMSolutionName(4), USLSolutionName(5)
C                   See description in Block Data Genral
C
C
      COMMON /LNAMEC/   LAKENAME, RiverName
      COMMON /LNAMEI/   NLN, LkToSlv2LN, NRN, LkSlv2RtLoopI
      COMMON /RCRUNC/   RCRUNTYP, IntervalName, MonthLengthName,
     &                  MLSolutionName, SMSolutionName, USLSolutionName
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /CLDATE/   dT, Day, Q, PeriodEnd, NumIncsVary,
     &                  DaysSinceSaturday
      COMMON /MISINT/   MidLakeSolution, BOMRound, FlowRound, 
     &                  LevelRound, RoutingInterval, RoutingIncrement,
     &                  StMarysSolution, UStLawrSolution, P77ACompat
      COMMON /RCRUNI/   LakesToSolve, 
     &                  IntervalNameLen, MonthLength,
     &                  StartYear, StartMonth, StartDay,
     &                  EndYear, EndMonth, EndDay
C
C
      L = LkSlv2RtLoopI(LakeNum)
      I = RoutingInterval(L)
C
c     IF( I .EQ. 4) THEN
C        If monthly, then all periods get labeled 'EOM'
c        PECODE = 'M'
c     ELSEIF( I .EQ. 3) THEN
C        If qmonthly, then check if it the last qmonth (i.e., EOM)
c        PECODE = 'Q'
c        IF( PeriodEnd(1,L) .GT. 27) PECODE = 'M'
c     ELSEIF( I .EQ. 2) THEN
C        If weekly, then all periods get labeled 'EOW' regardless of when end
c        PECODE = 'W'
c     ELSE
C        If daily check coincidence with month/qmonth
c        PECODE = 'D'
c        CALL GETQTRMONBOUNDS (PeriodEnd(4,L), PeriodEnd(3,L), QTRBNDS)
c        DO 33, D = 1,3
c  33       IF( PeriodEnd(1,L) .EQ. QTRBNDS(D) ) PECODE = 'Q'
c        IF( PeriodEnd(1,L) .EQ. QTRBNDS(4) ) PECODE = 'M'
c     ENDIF
C
c     WRITE(CTMP, 5997) LakeName(LakeNum)(1:NLN(LakeNum)), 
c    &   IntervalName(I)(1:IntervalNameLen(I)), PeriodEnd(1,L), 
c    &   PeriodEnd(3,L), PeriodEnd(4,L), PECODE, MeanLev, EOPLev, 
c    &   MeanFlow, LLT
      WRITE(CTMP, 5997) LakeName(LakeNum)(1:NLN(LakeNum)), 
     &   IntervalName(I)(1:IntervalNameLen(I)), PeriodEnd(1,L), 
     &   PeriodEnd(3,L), PeriodEnd(4,L), EOPType(L), MeanLev, EOPLev, 
     &   MeanFlow, LLT
      CALL LogIt(LOGD,0,1)
C 
 5997 FORMAT(A,' high-precision mean level, EOP level, and ', 
     &   'mean outflow for ', A, ' period ending', 2I3, I5, ' (EO', A1,
     &    '): ', 2F17.10,F15.5,A1)
C
      END
C
C
C
C
C-----------------------------------------------------------------------
C
      CHARACTER*1 FUNCTION EOPType( RoutLoopIndex ) 
C
C     This routine returns a character indicating whether the current 
C     period for a given routing loop (i.e., 1=Superior, 2=middle lakes, 
C     3=Ontario) ends a quarter-month or month.  If not, returns a 'D'
C     for daily routing.  Weekly routing always returns a 'W'.  The date
C     information is read from the PeriodEnd arrays in the COMMON block.
C 
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER          RoutLoopIndex
C                      (Input)  Index into RoutLoopIndex
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER          I
C                      Index to the routing interval
C
      INTEGER          QTRBNDS(4), D
C                      Holds the ending day of the month for 
C                      each quarter-month, and a scratch counter
C 
C
C     Declarations regarding COMMON block variables:
C
      INTEGER       PeriodEnd(6,3), RoutingInterval(3)
C                   See description in Block Data Genral
C
C                   Declared because of other variables in COMMON blocks:
      DOUBLE PRECISION dT(3)
      LOGICAL       NumIncsVary(3) 
      INTEGER       Day(3), Q(4,3), DaysSinceSaturday
      INTEGER       MidLakeSolution, StMarysSolution, UStLawrSolution 
      INTEGER       BOMRound, FlowRound, LevelRound, P77ACompat
      INTEGER       RoutingIncrement(3)
C
C
C
      COMMON /CLDATE/   dT, Day, Q, PeriodEnd, NumIncsVary,
     &                  DaysSinceSaturday
      COMMON /MISINT/   MidLakeSolution, BOMRound, FlowRound, 
     &                  LevelRound, RoutingInterval, RoutingIncrement,
     &                  StMarysSolution, UStLawrSolution, P77ACompat
C
      I = RoutingInterval(RoutLoopIndex)
C
      IF( I .EQ. 4) THEN
C        If monthly, then all periods get labeled 'EOM'
         EOPType = 'M'
      ELSEIF( I .EQ. 3) THEN
C        If qmonthly, then check if it the last qmonth (i.e., EOM)
         EOPType = 'Q'
         IF( PeriodEnd(1,RoutLoopIndex) .GT. 27) EOPType = 'M'
      ELSEIF( I .EQ. 2) THEN
C        If weekly, then all periods get labeled 'EOW' regardless of when end
         EOPType = 'W'
      ELSE
C        If daily check coincidence with month/qmonth
         EOPType = 'D'
         CALL GETQTRMONBOUNDS (PeriodEnd(4,RoutLoopIndex), 
     &      PeriodEnd(3,RoutLoopIndex), QTRBNDS)
         DO 33, D = 1,3
   33       IF( PeriodEnd(1,RoutLoopIndex) .EQ. QTRBNDS(D)) EOPType= 'Q'
         IF( PeriodEnd(1,RoutLoopIndex) .EQ. QTRBNDS(4) ) EOPType = 'M'
      ENDIF
C
      END
C
C
C
C-----------------------------------------------------------------------
C
      CHARACTER*(*) FUNCTION RightJustify( String, FW ) 
C
C     This routine pads String with enough leading spaces (max=40) to
C     push right-justify String in field of length FW.  The calling
C     routine must provide the dimension of the return array, using the
C     same value as FW.  (i.e., if justifying in a 24-character field,
C     the function declaration must be):
C          CHARACTER*(24) FUNCTION RightJustify
C     This is sadly lacking in general-purpose utility, but is just 
C     about as flexible as FORTRAN 77 can be.
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER String*(*)
C               (Input)  The input characters to right-justify  
C               
      INTEGER   FW
C               (Input)  The intended field width
C
      CHARACTER*40 Spaces
      DATA Spaces /'                                        '/
C
      IF( FW .EQ. LEN(String)) THEN
         RightJustify = String
      ELSE
         RightJustify = Spaces(1:FW - LEN(String))//String
      ENDIF
C
      END
C-----------------------------------------------------------------------
      SUBROUTINE UnReadInputData( F, ParID, LogCode )
C
      IMPLICIT NONE
C
      INTEGER      F
C                  (Input) Index to FileNames being examined
C
      CHARACTER*24 ParID
C                  (Input) String to write out as Parameter ID
C
      LOGICAL      LogCode(8)
C                  (Input) Describes which files to write output to
C
      INTEGER      NumLines
C                  Number of input lines (for repeating data)
C                  
      INTEGER      L
C                  Length of filename
C                  
      INTEGER      I, M 
C                  Scratch
C
C     INTEGER      ENDSTR
C                  Function that returns the end position of the current word.
C
      INTEGER      LASTKC
C                  Function returning the position of the last useful
C                  Character in a string.
C
      CHARACTER    CTMP*255, CWRK*80, LLT*1
      CHARACTER*256 FileNames(100)
      CHARACTER*256 RepeatLines(80)
C                  See description in Block Data Genral
C
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /FNAMEC/   FileNames, RepeatLines

      L = LastKC(FileNames(F), 256 )
      IF(FileNames(F)(1:L) .EQ. 'Unspecified')THEN
C        No input was given.  Write 'unused' line
         Write(CTMP, 9002) '#'//ParID, '             (not used)', LLT
         CALL LogIt(LOGCode,0,0)
      ELSE IF(FileNames(F)(1:4) .EQ. 'NULL' )THEN
         Write(CTMP, 9002) ParID, 'Null', LLT
         CALL LogIt(LOGCode,0,0)
      ELSE IF(FileNames(F)(1:15) .EQ. 'Repeating Data'//LLT)THEN
C        Repeating data line given, must re-create it         
C        Get the number of lines to process.
         READ(FileNames(F)(16:18), 1333) NumLines
         DO 7222, I = 1, NumLines
            READ(FileNames(F)(16+I*3:18+I*3), 1333) M
            Write(CTMP, 9002) ParID, RepeatLines(M)(1:228), LLT
 7222       CALL LogIt(LOGCode,0,0)        
      ELSE
C        Real file had been read in
         Write(CTMP, 9002) ParID, FileNames(F)(1:L), LLT
         CALL LogIt(LOGCode,0,0)        
      ENDIF
 1333 FORMAT(I3)
 9002 FORMAT (A,': ',A,A1)
      END
C----&------------------------------------------------------------------

      SUBROUTINE Summarize()
C
      IMPLICIT NONE
C
C     Declarations regarding local variables or functions:
C
      INTEGER      I, J, K, F(2)
      CHARACTER*4  OutFileCodes
C                  These are scratch variables.
C
      CHARACTER*(24) RightJustify
C
      INTEGER      LASTKC
C                  Function returning the position of the last useful
C                  Character in a string.
C
      INTEGER      ENDSTR
C                  Function that returns the end position of the current word.
C 
C 
C
C
C     Declarations regarding COMMON block variables:
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                  See description in Block Data Genral
C
      CHARACTER    DirSep*5,OSName(5)*8,MessDB*256
      CHARACTER*64 ModelerName
C                  See descriptions in Block Data Genral
C
      INTEGER      NLN(16), LkToSlv2LN(8), NRN(16), LkSlv2RtLoopI(8)
C                  See descriptions in Block Data Genral
C
      CHARACTER*21 LAKENAME(16), RiverName(16)
C                  See description in Block Data Genral
C
      CHARACTER*256 FileNames(100)
C                  See description in Block Data Genral
C
      CHARACTER*256 RepeatLines(80)
C                  See description in Block Data Genral
C
      CHARACTER    CTMP*255, CWRK*80
C                  These are just buffers used throughout program.
C
      CHARACTER    LLT*1
C                  This is the character that marks the end of a line sent
C                  to the LogIt subroutine.
C
      CHARACTER*3  MonthAbbrev(12)
C                  See description in Block Data Genral
C
      CHARACTER*80 RCTITLE(3), MLHEADER
C                  These are descriptive lines for use in output.
C
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C                  See descriptions in Block Data Genral
C
      INTEGER      MidLakeSolution, StMarysSolution, UStLawrSolution 
      INTEGER      BOMRound, FlowRound, LevelRound, P77ACompat
C                  See descriptions in Block Data Genral.
C                  Declared because of other variables in COMMON block.
C
      INTEGER      RoutingInterval(3), RoutingIncrement(3)
C                  Routing interval to use on Superior, middle lakes,
C                  and Ontario (1=daily, 2=weekly, 3=qmonthly, 
C                  4=monthly), and the number of increments within each 
C                  interval.
C
      INTEGER      PeriodEnd(6,3), Day(3), Q(4,3), DaysSinceSaturday
      DOUBLE PRECISION dT(3)
      LOGICAL      NumIncsVary(3) 
C                  See descriptions in RUNIT()
C
      INTEGER      MonthLength
      CHARACTER*8  MonthLengthName(2)
C                   See description in Block Data Genral
C
      CHARACTER*8  RCRUNTYP, IntervalName(4)
C                  See description in Block Data Genral
C
      CHARACTER*9  MLSolutionName(2)
C                  Strings set to 'interative' and 'matrix', used for output.
C
      CHARACTER*16 SMSolutionName(24), USLSolutionName(5)
      LOGICAL      WriteTSOutput(4,2,8), WriteErieOutflow
      LOGICAL      LakesToSolve(8)
C                  See description in Block Data Genral
C
      INTEGER      StartYear, StartMonth, StartDay, 
     &             EndYear, EndMonth, EndDay
C                  See descriptions in Block Data Genral
C
      INTEGER      IntervalNameLen(4)
C                  See description in Block Data Genral
C
      DOUBLE PRECISION SGX(17), SGY(17), SGZ(17)
      INTEGER       MaxSGSetting
      INTEGER      NUMAVSZ(8)
C                  See descriptions in Block Data Genral
C
      DOUBLE PRECISION AVSZ(8,5,2)
C                  See descriptions in Block Data Genral
C
      DOUBLE PRECISION SPSTLV, MHSTLV, SCSTLV, ERSTLV, CGSTLV, ONSTLV
C                  See descriptions in Block Data Genral
C
      DOUBLE PRECISION SPBOMLev, MHBOMLev, SCBOMLev, ERBOMLev
C                  The beginning-of-month lake levels used by Plan 1977A
C
      REAL         PrevSupGates 
C                  Initial gates open at Comp Works
C 
      INTEGER      PrevStMarysFlow
C                  Initial Lake Superior outflow (average for month 0)
C 
      REAL         SupGatesOpen
C                  The number of Superior gates open in current routing 
C                  period.
C
      DOUBLE PRECISION SPK,SPYM,SPA,SPB,SPWT,SPC
      DOUBLE PRECISION MHK,MHYM,MHA,MHB,MHWT,MHC
      DOUBLE PRECISION SCK,SCYM,SCA,SCB,SCWT,SCC
      DOUBLE PRECISION ERK,ERYM,ERA,ERB,ERWT,ERC
      DOUBLE PRECISION ONK,ONYM,ONA,ONB,ONWT,ONC
C                  These are the stage-discharge equation K, ym, a, b, zc,
C                  and wt, as described by GLERL's Tech Memo ERL-GLERL-xx,
C                  and in the function OUTFLOW.  Used for when the lake
C                  outflow is determined from an outflow equation, instead
C                  of a regulation plan.
C
C
      INTEGER      TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &             MissingVal, SecsInUniformMon
      DOUBLE PRECISION MPFT, HALF, ONE, ZERO
C                  See descriptions in Block Data General.  
C
      INTEGER   CritDepths, P58DTotalSup, MaxILimit, MinPLimit,
     &          P58DForecastDay, P58DRunType,
     &          P58DForecastYear, P58DForecastMon,  
     &          P58DWeightedSupply, P58DRuleCurve, P58DAppLimit,
     &          SupplyIndicator(21)
      Integer   P58DAdjustment, P58DrcAdjustment, P58DActualFlow, 
     &          P58DTotalDev,P58DLSLflow
      DOUBLE PRECISION P58DActualLevel, P58DPlanLevel
      DOUBLE PRECISION PreprojectLevel
      CHARACTER*2  P58DAppLimitSym
      CHARACTER*256 P58DSimulateFile, P58DSIDataFile
      LOGICAL   ElimDev, SetFlow, PreProjConditions, ForceFlow
      LOGICAL   P58DUseDatabase, P58DSiDatabase, ComputeMontreal
      CHARACTER SDATE*9,STIME*8
      INTEGER   NOPDIR, NOPEXT, TSOutputFileNum(4,2,8)
      LOGICAL   OnlyOneOpen
      CHARACTER*6  UnitType(20)
      INTEGER      UnitTypeLen(20)
      INTEGER      BaseUnitType(20)
      DOUBLE PRECISION UnitConvFactor(20)
C               See descriptions in Block Data General.  
C
C               These variables included just beacuse they are
C               in same COMMON as LStoSWP and SupCritRnd
      INTEGER   ForecastMonths, LStoSWP, SupCritRnd, MonthsChecked
      INTEGER   SPOPNBS06(12,6), SPOPNBS07(12,6), OPLLO(12), 
     &          MHOPNBS06(12,6), MHOPNBS07(12,6), OPChicago, 
     &          SCRRet(12), SCOPNBS06(12,6), SCOPNBS07(12,6), 
     &          DetRet(12), EROPNBS06(12,6), EROPNBS07(12,6),
     &          OPWelland(12), NiaRet(12),
     &          SCRRet_P2012(12), SCOPNBS08(12,6), SPOPNBS08(12,6), 
     &          DetRet_P2012(12), EROPNBS08(12,6), MHOPNBS08(12,6),
     &          OPWelland_P2012(12), NiaRet_P2012(12)       
      INTEGER      P77AWinterMax, P77AMaxChange, P77AMinimumFlow(12),
     &             P77ABaseFlow(12), P77AMonthsToAve(12),
     &             P77ASup05NBS(12), P77ASup50NBS(12), P77ASup95NBS(12), 
     &             P77ALLOgokiDiv(12),P77AStMarysRetrd(12), 
     &             P77AMH05NBS(12), P77AMH50NBS(12), P77AMH95NBS(12), 
     &             P77AChicagoDiv(12), P77AStClairRetrd(12), 
     &             P77AStC50NBS(12), P77ADetroitRetrd(12),
     &             P77AErie50NBS(12), P77AWellandFlow(12), 
     &             P77ANiagaraRetrd(12),
     &             P77AMidLakesInc, P77AMidLakesTol, 
     &             P77AMonthRound, P77AFlowRound, P77ALevelRound
      DOUBLE PRECISION P77ASupTargetLev(12), P77ASupRangeFact(12),
     &             P77AMHTargetLev(12), P77AMHRangeFact(12), P77AAValue,
     &             P77ACGIP, P77ALowTrigger, P77AHighTrigger, 
     &             P77ASPK, P77ASPYM, P77ASPA, P77ASPB,P77ASPWT,P77ASPC, 
     &             P77AMHK, P77AMHYM, P77AMHA, P77AMHB,P77AMHWT,P77AMHC, 
     &             P77ASCK, P77ASCYM, P77ASCA, P77ASCB,P77ASCWT,P77ASCC, 
     &             P77AERK, P77AERYM, P77AERA, P77AERB,P77AERWT,P77AERC
      LOGICAL      P77AUniformMonth
      LOGICAL      AlreadySpecified(420)
      INTEGER      CritBOMRound, CritFlowRound, CritLevelRound, 
     &             Q1887Method, CritBMethod
      DOUBLE PRECISION CritCThresholdEl
      LOGICAL      BackCompat
C
      CHARACTER*256 CritKfile
C
C
C
C
      COMMON /CRITK/    CritKfile
      COMMON /SUPCRT/   CritBOMRound, CritFlowRound, CritLevelRound, 
     &                  Q1887Method, CritBMethod, CritCThresholdEl
      COMMON /PARMRD/   AlreadySpecified
      COMMON /SUPCOM/   SPOPNBS06, OPLLO, MHOPNBS06, OPChicago, SCRRet,
     &                  SCOPNBS06, DetRet, EROPNBS06, OPWelland, NiaRet,
     &                  SPOPNBS07, MHOPNBS07, SCOPNBS07, EROPNBS07,
     &                  SPOPNBS08, MHOPNBS08, SCOPNBS08, EROPNBS08,
     &                  OPWelland_P2012, SCRRet_P2012, DetRet_P2012,
     &                  NiaRet_P2012,
     &                  ForecastMonths,
     &                  LStoSWP, SupCritRnd, BackCompat, MonthsChecked
      COMMON /UNITSC/   UnitType
      COMMON /UNITSI/   UnitTypeLen, BaseUnitType, UnitConvFactor
      COMMON /OUTSTF/   WriteTSOutput, WriteErieOutflow
      COMMON /MONAME/   MonthAbbrev
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /LNAMEC/   LAKENAME, RiverName
      COMMON /LNAMEI/   NLN, LkToSlv2LN, NRN, LkSlv2RtLoopI
      COMMON /RCAVSZ/   NUMAVSZ, AVSZ
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /STAGE/    MHK,SCK,ERK,MHYM,SCYM,ERYM,MHA,SCA,ERA,
     &                  MHB,SCB,ERB,MHWT,SCWT,ERWT,MHC,SCC,ERC
      COMMON /SPDAT/    SPK,SPYM,SPA,SPB,SPWT,SPC,
     &                  SGX, SGY, SGZ, MaxSGSetting,
     &                  SupGatesOpen
      COMMON /ONDAT/    ONK,ONYM,ONA,ONB,ONWT,ONC
      COMMON /MISINT/   MidLakeSolution, BOMRound, FlowRound, 
     &                  LevelRound, RoutingInterval, RoutingIncrement,
     &                  StMarysSolution, UStLawrSolution, P77ACompat
      COMMON /INITLV/   SPSTLV, MHSTLV, SCSTLV, ERSTLV, CGSTLV, ONSTLV,
     &                  PrevStMarysFlow, PrevSupGates, 
     &                  SPBOMLev, MHBOMLev, SCBOMLev, ERBOMLev  
      COMMON /FNAMEC/   FileNames, RepeatLines
      COMMON /RCRUNC/   RCRUNTYP, IntervalName, MonthLengthName,
     &                  MLSolutionName, SMSolutionName,
     &                  USLSolutionName
      COMMON /RCRUNI/   LakesToSolve, 
     &                  IntervalNameLen, MonthLength,
     &                  StartYear, StartMonth, StartDay,
     &                  EndYear, EndMonth, EndDay
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
      COMMON /CLDATE/   dT, Day, Q, PeriodEnd, NumIncsVary,
     &                  DaysSinceSaturday
      COMMON /INI58D/   CritDepths, P58DTotalSup, MaxILimit, MinPLimit,
     &                  P58DWeightedSupply, P58DRuleCurve, P58DAppLimit,
     &                  P58DRunType, P58DForecastYear, P58DForecastMon,
     &                  P58DForecastDay, SupplyIndicator, 
     &                  P58DPlanLevel, P58DActualLevel,
     &                  P58DSimulateFile, P58DSIDataFile, ElimDev, 
     &                  ForceFlow, SetFlow, PreProjConditions,
     &                  P58DUseDatabase, P58DSiDatabase
      COMMON /Ini58dc/  P58DAppLimitSym
      COMMON /INI58D1/  P58DAdjustment,P58DrcAdjustment,P58DActualFlow, 
     &                  P58DTotalDev,P58DLSLflow
      COMMON /INI58D2/  ComputeMontreal
      common /precond/  PreprojectLevel

      COMMON /TITLES/   RCTITLE, MLHEADER
      COMMON /STIME/    STIME,SDATE
      COMMON /GENRLC/   OSName,DirSep,MessDB, ModelerName
      CHARACTER RCOUTDIR*128, RCOUTEXT*16
C               See descriptions in Block Data Genral
C
      COMMON /OUTNAM/   RCOUTDIR, RCOUTEXT
      COMMON /OUTNMI/   NOPDIR, NOPEXT, TSOutputFileNum, OnlyOneOpen
C
      COMMON /P1977A/   P77ASPK, P77ASPYM, P77ASPA, P77ASPB, P77ASPWT, 
     &                  P77ASPC, P77AMHK, P77AMHYM, P77AMHA, P77AMHB, 
     &                  P77AMHWT, P77AMHC, P77ASCK, P77ASCYM, P77ASCA,
     &                  P77ASCB, P77ASCWT, P77ASCC, P77AERK, P77AERYM,
     &                  P77AERA, P77AERB, P77AERWT, P77AERC, P77ACGIP,
     &                  P77ASupTargetLev, P77AMHTargetLev, P77AAValue, 
     &                  P77ASupRangeFact, P77AMHRangeFact, 
     &                  P77ALowTrigger, P77AHighTrigger, P77AMaxChange,
     &                  P77AMonthsToAve, P77AMinimumFlow, P77AWinterMax,
     &                  P77ASup05NBS, P77ASup50NBS, P77ASup95NBS,
     &                  P77AMH05NBS, P77AMH50NBS, P77AMH95NBS,
     &                  P77AStC50NBS, P77AErie50NBS, P77AWellandFlow,
     &                  P77AChicagoDiv, P77ALLOgokiDiv, P77ABaseFlow,
     &                  P77AStMarysRetrd, P77AStClairRetrd, 
     &                  P77ADetroitRetrd, P77ANiagaraRetrd, 
     &                  P77AMidLakesInc, P77AUniformMonth, 
     &                  P77AMonthRound, P77AFlowRound, P77ALevelRound,
     &                  P77AMidLakesTol
C
C
      DATA OutFileCodes /'DWQM'/
C     
C     Write a little separator block
      DO 711, K=1,79
  711    CWRK(K:K)='#'
      WRITE(CTMP,9000) CWRK(1:79), LLT
      DO 712, K=1,5
  712    CALL LogIt(LOGSD,0,0)
C
C     ==================     List system info
C
      WRITE(CTMP,9000) 'Complete Input Listing of CGLRRM Run at '//
     &   STIME//', '//SDATE, LLT
      CALL LogIt(LOGSD,1,1)
      WRITE(CTMP,9000) 'You should be able to reproduce the '//
     &   'results of this run simply by using', LLT
      CALL LogIt(LOGSD,0,0)
      WRITE(CTMP,9000) 'this section as an input parameter '//
     &   'file (provided the same external', LLT
      CALL LogIt(LOGSD,0,0)
      WRITE(CTMP,9000) 'time-series data files are available '//
     &   'in the same locations, and excepting', LLT
      CALL LogIt(LOGSD,0,0)
      WRITE(CTMP,9000) 'differences due to model version).', LLT
      CALL LogIt(LOGSD,0,0)
      CALL BuildInfo
      CTMP = '#   '//CTMP(1:251)
      CALL LogIt(LOGSD,0,0)
      WRITE(CTMP,9001) 'System Settings', LLT
      CALL LogIt(LOGSD,1,0)
      WRITE(CTMP,9002) RightJustify( 'Message Database', 24), 
     &   MessDB(1:LASTKC(MessDB,256)), LLT
      CALL LogIt(LOGSD,0,0)
      WRITE(CTMP,9003) RightJustify('Only One Open', 24),OnlyOneOpen,LLT
      CALL LogIt(LOGSD,0,0)
      WRITE(CTMP,9002) RightJustify( 'Modeler Name', 24), 
     &   ModelerName(1:LASTKC(ModelerName,64)), LLT
      CALL LogIt(LOGSD,0,0)
C
C
C
C     ==================     List general job info
C
      WRITE(CTMP,9001) 'Job Titles', LLT
      CALL LogIt(LOGSD,1,0)
      DO 733, K=1,3
         J = EndStr(RCTITLE(K),80)
         IF ( J .EQ. 0 ) THEN 
            WRITE(CTMP,9004) '#', K, '(not used)',LLT
         ELSE
            WRITE(CTMP,9004) ' ',K, 
     &         RCTITLE(K)(1:LastKC(RCTITLE(K),80)), LLT
         ENDIF
  733    CALL LogIt(LOGSD,0,0)
C
      CWRK(1:5) = 'SMCEO'
      DO 738, K=1,5
  738    IF(.NOT. LakesToSolve(K)) CWRK(K:K) = ' '
      WRITE(CTMP,9013) RightJustify( 'Lakes To Solve For', 24), 
     &   RightJustify( CWRK(1:5), 24), LLT
      CALL LogIt(LOGSD,1,0)
C
      WRITE(CTMP,9001) 'Simulation Period', LLT
      CALL LogIt(LOGSD,1,0)
      WRITE(CTMP,9006) RightJustify( 'Start Date', 24), StartYear, 
     &   StartMonth, StartDay, LLT
      CALL LogIt(LOGSD,0,0)
      WRITE(CTMP,9006) RightJustify( 'End Date', 24), EndYear, 
     &   EndMonth, EndDay, LLT
      CALL LogIt(LOGSD,0,0)
      WRITE(CTMP,9013) RightJustify( 'Month Length', 24), 
     &   RightJustify( MonthLengthName(MonthLength)
     &      (1:EndStr(MonthLengthName(MonthLength),8)), 24), LLT
      CALL LogIt(LOGSD,0,0)
C
      WRITE(CTMP,9001) 'Output File Characteristics', LLT
      CALL LogIt(LOGSD,1,0)
      IF( NOPDIR .GT. 0) THEN
         WRITE(CTMP,9002) RightJustify( 'Output Directory', 24), 
     &      RCOUTDIR(1:NOPDIR), LLT
      ELSE
         WRITE(CTMP,9001) '   Output Directory not specified', LLT
      ENDIF
      CALL LogIt(LOGSD,0,0)
      IF( NOPEXT .GT. 0) THEN
         WRITE(CTMP,9013) RightJustify( 'Output Extension', 24), 
     &      RightJustify( RCOUTEXT(1:NOPEXT), 24), LLT
      ELSE
         WRITE(CTMP,9001) '   Output File Extension not specified', LLT
      ENDIF
      CALL LogIt(LOGSD,0,0)
C
      WRITE(CTMP,9001) 'Misc Levels of Debugging Output', LLT
      CALL LogIt(LOGSD,1,0)
      WRITE(CTMP,9005) RightJustify( 'General Verbosity', 24), 
     &   GenVerbosity, LLT
      CALL LogIt(LOGSD,1,0)
      WRITE(CTMP,9005) RightJustify( 'TS Data Verbosity', 24), 
     &   TSDVerbosity, LLT
      CALL LogIt(LOGSD,0,0)
      WRITE(CTMP,9005) RightJustify( 'Check Verbosity', 24), 
     &   ChkVerbosity, LLT
      CALL LogIt(LOGSD,0,0)
C
      WRITE(CTMP,9001) 'General Rounding', LLT
      CALL LogIt(LOGSD,1,0)
      WRITE(CTMP,9005) RightJustify( 'Flow Round', 24), FlowRound, LLT
      CALL LogIt(LOGSD,0,0)
      WRITE(CTMP,9005) RightJustify( 'Level Round', 24), LevelRound, LLT
      CALL LogIt(LOGSD,0,0)
      WRITE(CTMP,9005) RightJustify( 'BOM Round', 24), BOMRound, LLT
      CALL LogIt(LOGSD,0,0)
C
C     Lake Superior
      IF(LakesToSolve(1)) THEN
         WRITE(CTMP,9001) LakeName(1)(1:NLN(1))//' data', LLT
         CALL LogIt(LOGSD,1,0)
         WRITE(CTMP,9005) RightJustify( 'Sup Verbosity', 24), 
     &      SupVerbosity, LLT
         CALL LogIt(LOGSD,0,0)
         WRITE(CTMP,9013) RightJustify( 'Sup. Time Step', 24), 
     &      RightJustify( IntervalName(RoutingInterval(1))
     &      (1:IntervalNameLen(RoutingInterval(1))), 24), LLT
         CALL LogIt(LOGSD,0,0)
         IF(NumIncsVary(1)) THEN
            WRITE(CTMP,9013) RightJustify( 'Sup. Inc.', 24), 
     &         RightJustify( 'Hourly', 24), LLT
         ELSE
            WRITE(CTMP,9005) RightJustify( 'Sup. Inc.', 24), 
     &         RoutingIncrement(1), LLT
         ENDIF
         CALL LogIt(LOGSD,0,0)
C        Figure out which output files were requested
         I = 0
         DO 164, K=1,4
            IF(WriteTSOutput(K, 1, 1)) THEN
               I = I + 1
               CWRK(I:I) = OutFileCodes(K:K)
            ENDIF
  164       CONTINUE
         IF( I .EQ. 0 ) THEN  
            WRITE(CTMP,9013) RightJustify( 'Sup Output Files', 24), 
     &         RightJustify( 'none', 24), LLT
         ELSE
            WRITE(CTMP,9013) RightJustify( 'Sup Output Files', 24), 
     &         RightJustify( CWRK(1:I), 24), LLT
         ENDIF
         CALL LogIt(LOGSD,1,0)
         WRITE(CTMP,9007) RightJustify( 'Sup. Start Level', 24), 
     &      SPSTLV, UnitType(7)(1:UnitTypeLen(7)), LLT
         CALL LogIt(LOGSD,0,0)
C        Print Area input
         IF(NUMAVSZ(1) .EQ. 1) THEN
            WRITE(CTMP,9009) RightJustify( 'Sup. Area/Elev', 24), 
     &         AVSZ(1,1,2), UnitType(9)(1:UnitTypeLen(9)), LLT
            CALL LogIt(LOGSD,0,0)
         ELSE
            DO 2098, J=1,NUMAVSZ(1)
               WRITE(CTMP,9008) RightJustify( 'Sup. Area/Elev', 24), 
     &            AVSZ(1,J,2), UnitType(9)(1:UnitTypeLen(9)), 
     &            AVSZ(1,J,1), UnitType(7)(1:UnitTypeLen(7)), LLT
 2098          CALL LogIt(LOGSD,0,0)
         ENDIF
C
         IF(StMarysSolution .EQ. 1) THEN
            WRITE(CTMP,9002) RightJustify( 'St. M. Riv Solution', 24), 
     &         RightJustify(SMSolutionName(1)(1:16), 24), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9010) RightJustify( 'St. M. Riv Eqn', 24), 
     &         SPK, SPYM, SPA, SPB, SPWT, SPC, LLT
            CALL LogIt(LOGSD,0,0)
         ELSE IF(StMarysSolution .EQ. 3) THEN
            WRITE(CTMP,9001) 'Superior Regulation Parameters', LLT
            CALL LogIt(LOGSD,1,0)
            WRITE(CTMP,9002) RightJustify( 'St. M. Riv Solution', 24), 
     &         RightJustify('1977A',24), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9005) RightJustify( 'P77A Compat. Mode', 24), 
     &         P77ACompat, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9011) RightJustify('Prev. St. M. Riv. Flow', 24), 
     &         PrevStMarysFlow, UnitType(2)(1:UnitTypeLen(2)), LLT
            CALL LogIt(LOGSD,0,0)
            IF( NINT(PrevSupGates) .NE. MissingVal ) THEN
               WRITE(CTMP,9012) RightJustify( 'Prev. Sup. Gates Open',
     &             24), PrevSupGates, LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
            WRITE(CTMP,9005) RightJustify( 'Sup. SW Pier Formula', 24), 
     &         LStoSWP, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9005) RightJustify( 'Sup. Criteria Round', 24), 
     &         SupCritRnd, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9005) RightJustify( 'Sup. Pre-Proj. Method', 24), 
     &         Q1887Method, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9005) RightJustify( 'Sup. Crit. B Method', 24), 
     &         CritBMethod, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9007) RightJustify( 'Sup. Crit. C Thresh', 24), 
     &         CritCThresholdEl, LLT
            CALL LogIt(LOGSD,0,0)
            IF( IDNINT(SPBOMLev) .NE. MissingVal ) THEN
               WRITE(CTMP,9007) RightJustify( 'Sup. BOM Level', 24), 
     &            SPBOMLev, UnitType(7)(1:UnitTypeLen(7)), LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
            IF( IDNINT(MHBOMLev) .NE. MissingVal ) THEN
               WRITE(CTMP,9007) RightJustify( 'MHu BOM Level', 24), 
     &            MHBOMLev, UnitType(7)(1:UnitTypeLen(7)), LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
            IF( IDNINT(SCBOMLev) .NE. MissingVal ) THEN
               WRITE(CTMP,9007) RightJustify( 'StC BOM Level', 24), 
     &            SCBOMLev, UnitType(7)(1:UnitTypeLen(7)), LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
            IF( IDNINT(ERBOMLev) .NE. MissingVal ) THEN
               WRITE(CTMP,9007) RightJustify( 'Eri. BOM Level', 24), 
     &            ERBOMLev, UnitType(7)(1:UnitTypeLen(7)), LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
C  
C           The following Plan 1977A values have defaults in model.
C           Must determine if override had been explicitly given.
C
C      DATA PARIDS(2785:2800) /'P77AMINIMUMFLOW '/                        175
            IF (AlreadySpecified(175) ) THEN
               WRITE(CTMP,9014) RightJustify('P77A Minimum Flow', 24), 
     &            (P77AMinimumFlow(J), J=1,12), LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
C      DATA PARIDS(3025:3040) /'P77ANIARIVEQN   '/                        190
            IF (AlreadySpecified(190) ) THEN
               WRITE(CTMP,9010) RightJustify( 'P77A Nia. Riv Eqn', 24), 
     &            P77AERK,P77AERYM,P77AERA,P77AERB,P77AERWT,P77AERC,LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
C      DATA PARIDS(3057:3072) /'P77ASTCRIVEQN   '/                        192
            IF (AlreadySpecified(192) ) THEN
               WRITE(CTMP,9010) RightJustify( 'P77A St. C. Riv Eqn',24), 
     &            P77AMHK,P77AMHYM,P77AMHA,P77AMHB,P77AMHWT,P77AMHC,LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
C      DATA PARIDS(3073:3088) /'P77ASTMRIVEQN   '/                        193
            IF (AlreadySpecified(193) ) THEN
               WRITE(CTMP,9010) RightJustify( 'P77A St. M. Riv Eqn',24), 
     &            P77ASPK,P77ASPYM,P77ASPA,P77ASPB,P77ASPWT,P77ASPC,LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
C      DATA PARIDS(3041:3056) /'P77ADETRIVEQN   '/                        191
            IF (AlreadySpecified(191) ) THEN
               WRITE(CTMP,9010) RightJustify( 'P77A Det. Riv Eqn', 24), 
     &            P77ASCK,P77ASCYM,P77ASCA,P77ASCB,P77ASCWT,P77ASCC,LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
C      DATA PARIDS(3089:3104) /'P77AMIDLAKEINC  '/                        194
            IF (AlreadySpecified(194) ) THEN
               WRITE(CTMP,9005) RightJustify( 'P77A Midlake Inc.', 24), 
     &            P77AMidLakesInc, LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
C      DATA PARIDS(3137:3152) /'P77AMONTHROUND  '/                        197
            IF (AlreadySpecified(197) ) THEN
               WRITE(CTMP,9005) RightJustify( 'P77A Month Round', 24), 
     &           P77AMonthRound, LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
C      DATA PARIDS(3153:3168) /'P77ALEVELROUND  '/                        198
            IF (AlreadySpecified(198) ) THEN
               WRITE(CTMP,9005) RightJustify( 'P77A Level Round', 24), 
     &           P77ALevelRound, LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
C      DATA PARIDS(3169:3184) /'P77AFLOWROUND   '/                        199
            IF (AlreadySpecified(199) ) THEN
               WRITE(CTMP,9005) RightJustify( 'P77A Flow Round', 24), 
     &           P77AFlowRound, LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
C      DATA PARIDS(3185:3200) /'P77AMONTHSTOAVG '/                        200
            IF (AlreadySpecified(200) ) THEN
               WRITE(CTMP,9014) RightJustify( 'P77A Months to Avg.',24), 
     &            (P77AMonthsToAve(J), J=1,12), LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
C      DATA PARIDS(3201:3216) /'P77AACONSTANT   '/                        201
            IF (AlreadySpecified(201) ) THEN
               WRITE(CTMP,9012) RightJustify( 'P77A A Constant', 24), 
     &            P77AAValue, LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
C      DATA PARIDS(3217:3232) /'P77AMAXCHANGE   '/                        202
            IF (AlreadySpecified(202) ) THEN
               WRITE(CTMP,9005) RightJustify( 'P77A Max. Change', 24), 
     &           P77AMaxChange, LLT
               CALL LogIt(LOGSD,0,0)
            ENDIF
C      DATA PARIDS(3233:3248) /'P77ALOWTRIGGER  '/                        203
            IF (AlreadySpecified(203) ) THEN
               WRITE(CTMP,9007) RightJustify( 'P77A Low Trigger', 24), 
     &            P77ALowTrigger, UnitType(7)(1:UnitTypeLen(7)), LLT
            ENDIF
C      DATA PARIDS(3249:3264) /'P77AHIGHTRIGGER '/                        204
            IF (AlreadySpecified(204) ) THEN
               WRITE(CTMP,9007) RightJustify( 'P77A High Trigger', 24), 
     &            P77AHighTrigger, UnitType(7)(1:UnitTypeLen(7)), LLT
            ENDIF
         ENDIF
C        Determine if Time-Series Data file given
         I = INDEX(FileNames(41),LLT) + 1
         CWRK = FileNames(41)(I:I+6)
         IF(FileNames(41)(I:I+6) .EQ. 'default') THEN
C           Pre-written TSD file not given; had to specify Time-Series Data 
C           individually.  Print time-series data inputs:
            CALL UnReadInputData( 42, 
     &         RightJustify('Sup. NBS Data', 24), LOGSD)
            CALL UnReadInputData( 43, 
     &         RightJustify('Sup. Con. Use Data', 24), LOGSD)
            CALL UnReadInputData( 44, 
     &         RightJustify('Long Lac Div. Data', 24), LOGSD)
            CALL UnReadInputData( 45, 
     &         RightJustify('Ogoki Div. Data', 24), LOGSD)
            CALL UnReadInputData( 46, 
     &         RightJustify('Comb. LLO Div. Data', 24), LOGSD)
            CALL UnReadInputData( 47, 
     &         RightJustify('Sup. Other Div. Data', 24), LOGSD)
            CALL UnReadInputData( 48, 
     &         RightJustify('Sup. Other Sup. Data', 24), LOGSD)
            CALL UnReadInputData( 49, 
     &         RightJustify('Sup. Deviation Data', 24), LOGSD)
            CALL UnReadInputData( 40, 
     &         RightJustify('Sup. Gates Data', 24), LOGSD)
            CALL UnReadInputData( 69, 
     &         RightJustify('Sup. Side Chan. Data', 24), LOGSD)
            CALL UnReadInputData( 68, 
     &         RightJustify('Sup. Outflow Data', 24), LOGSD)
            CALL UnReadInputData( 70, 
     &         RightJustify('Sup. EOP Level Data', 24), LOGSD)
            CALL UnReadInputData( 50, 
     &         RightJustify('St. M. Riv. Retr. Data', 24), LOGSD)
            CALL UnReadInputData( 51, 
     &         RightJustify('St. M. Riv. Ice Data', 24), LOGSD)
            CALL UnReadInputData( 52, 
     &         RightJustify('St. M. Riv. Weed Data', 24), LOGSD)
         ELSE
C           Pre-written TSD file was given; no other Time-Series Data given
            CALL UnReadInputData( 41, 
     &         RightJustify('Sup. Sum. Input File', 24), LOGSD)
         ENDIF
C
      ELSE
         WRITE(CTMP,9001) LakeName(1)(1:NLN(1))//' not modeled', LLT
         CALL LogIt(LOGSD,1,0)
      ENDIF
C
C
      IF(LakesToSolve(2) .OR. LakesToSolve(3) .OR. LakesToSolve(4)) THEN
C        Stuff pertaining to Middle lakes in general
         WRITE(CTMP,9001) 'Middle Lakes Settings', LLT
         CALL LogIt(LOGSD,1,0)
C
         WRITE(CTMP,9005) RightJustify( 'Midlake Verbosity', 24), 
     &      MLVerbosity, LLT
         CALL LogIt(LOGSD,1,0)
         WRITE(CTMP,9002) RightJustify( 'Routing Title', 24),
     &      MLHEADER(1:LASTKC(MLHEADER,80)),LLT
         CALL LogIt(LOGSD,0,0)
         WRITE(CTMP,9013) RightJustify( 'MidLake Time Step', 24), 
     &      RightJustify( IntervalName(RoutingInterval(2))
     &      (1:IntervalNameLen(RoutingInterval(2))), 24), LLT
         CALL LogIt(LOGSD,0,0)
         IF(NumIncsVary(2)) THEN
            WRITE(CTMP,9013) RightJustify( 'MidLake Inc.', 24), 
     &         RightJustify( 'Hourly', 24), LLT
         ELSE
            WRITE(CTMP,9005) RightJustify( 'MidLake Inc.', 24), 
     &         RoutingIncrement(2), LLT
         ENDIF
         CALL LogIt(LOGSD,0,0)
         F(1) = 9
         F(2) = 6
         WRITE(CTMP,9002) RightJustify( 'MidLake Solution', 24), 
     &         RightJustify(MLSolutionName(MidLakeSolution)
     &         (1:F(MidLakeSolution)), 24), LLT
         CALL LogIt(LOGSD,0,0)
      ENDIF
C
C     MHu stuff
      IF(LakesToSolve(2)) THEN
         WRITE(CTMP,9001) LakeName(2)(1:NLN(2))//' data', LLT
         CALL LogIt(LOGSD,1,0)
C        Figure out which output files were requested
         I = 0
         DO 264, K=1,4
            IF(WriteTSOutput(K, 1, 2)) THEN
               I = I + 1
               CWRK(I:I) = OutFileCodes(K:K)
            ENDIF
  264       CONTINUE
         IF( I .EQ. 0 ) THEN  
            WRITE(CTMP,9013) RightJustify( 'MHu Output Files', 24), 
     &         RightJustify( 'none', 24), LLT
         ELSE
            WRITE(CTMP,9013) RightJustify( 'MHu Output Files', 24), 
     &         RightJustify( CWRK(1:I), 24), LLT
         ENDIF
         CALL LogIt(LOGSD,1,0)
         WRITE(CTMP,9010) RightJustify( 'St. C. Riv Eqn', 24), 
     &      MHK, MHYM, MHA, MHB, MHWT, MHC, LLT
         CALL LogIt(LOGSD,0,0)
         WRITE(CTMP,9007) RightJustify( 'MHu Start Level', 24), 
     &      MHSTLV, UnitType(7)(1:UnitTypeLen(7)), LLT
         CALL LogIt(LOGSD,0,0)
C        Print Area input
         IF(NUMAVSZ(2) .EQ. 1) THEN
            WRITE(CTMP,9009) RightJustify( 'MHu Area/Elev', 24), 
     &         AVSZ(2,1,2), UnitType(9)(1:UnitTypeLen(9)), LLT
            CALL LogIt(LOGSD,0,0)
         ELSE
            DO 2198, J=1,NUMAVSZ(2)
               WRITE(CTMP,9008) RightJustify( 'MHu Area/Elev', 24), 
     &            AVSZ(2,J,2), UnitType(9)(1:UnitTypeLen(9)), 
     &            AVSZ(2,J,1), UnitType(7)(1:UnitTypeLen(7)), LLT
 2198          CALL LogIt(LOGSD,0,0)
         ENDIF
C        Determine if Time-Series Data file given
         I = INDEX(FileNames(31),LLT) + 1
         CWRK = FileNames(31)(I:I+6)
         IF(FileNames(31)(I:I+6) .EQ. 'default') THEN
C           Pre-written TSD file not given; had to specify Time-Series Data 
C           individually.  Print time-series data inputs:
C                         FileNames(19) = STMRIVFLOWDATA					 84
            CALL UnReadInputData( 19, 
     &         RightJustify('St. M. Riv. Flow Data', 24), LOGSD)
C                         FileNames(1)  = MHUNBSDATA                       38
            CALL UnReadInputData(  1, 
     &         RightJustify('MHu NBS Data', 24), LOGSD)
C                         FileNames(22) = CHICAGODIVDATA					 40
            CALL UnReadInputData( 22, 
     &         RightJustify('Chicago Div. Data', 24), LOGSD)
C                         FileNames(24) = MHUOTHERDIVDATA					 41
            CALL UnReadInputData( 24, 
     &         RightJustify('MHu Other Div. Data', 24), LOGSD)
C                         FileNames(23) = MHUCONUSEDATA					 39
            CALL UnReadInputData( 23, 
     &         RightJustify('MHu Con. Use Data', 24), LOGSD)
C                         FileNames(25) = MHUOTHERSUPDATA					 42
            CALL UnReadInputData( 25, 
     &         RightJustify('MHu. Other Sup. Data', 24), LOGSD)
C                         FileNames(4)  = MICCONUSEDATA					 49
            CALL UnReadInputData( 4, 
     &         RightJustify('Mic Con. Use Data', 24), LOGSD)
C                         FileNames(12) = MICOTHERDIVDATA					 50
            CALL UnReadInputData( 12, 
     &         RightJustify('Mic Other Div. Data', 24), LOGSD)
C                         FileNames(8)  = MICOTHERSUPDATA					 51
            CALL UnReadInputData( 8, 
     &         RightJustify('Mic. Other Sup. Data', 24), LOGSD)
C                         FileNames(5)  = HURCONUSEDATA					 56
            CALL UnReadInputData( 5, 
     &         RightJustify('Hur. Con. Use Data', 24), LOGSD)
C                         FileNames(13) = HUROTHERDIVDATA					 57
            CALL UnReadInputData( 13, 
     &         RightJustify('Hur. Other Div. Data', 24), LOGSD)
C                         FileNames(9)  = HUROTHERSUPDATA 				 58
            CALL UnReadInputData( 9, 
     &         RightJustify('Hur. Other Sup. Data', 24), LOGSD)
C                         FileNames(71) = MHUEOPLEVELDATA					 46
            CALL UnReadInputData( 71, 
     &         RightJustify('MHu EOP Level Data', 24), LOGSD)
C                         FileNames(16) = STCRIVRETRDATA					 92
            CALL UnReadInputData( 16, 
     &         RightJustify('St. C. Riv. Retr. Data', 24), LOGSD)
C                         FileNames(34) = STCRIVICEDATA					 90
            CALL UnReadInputData( 34,
     &         RightJustify('St. C. Riv. Ice Data', 24), LOGSD)
C                         FileNames(35) = STCRIVWEEDDATA					 91
            CALL UnReadInputData( 35,
     &         RightJustify('St. C. Riv. Ice Data', 24), LOGSD)
         ELSE
C           Pre-written TSD file was given; no other Time-Series Data given
            CALL UnReadInputData( 31, 
     &         RightJustify('MHu Sum. Input File', 24), LOGSD)
         ENDIF
      ELSE
         WRITE(CTMP,9001) LakeName(2)(1:NLN(2))//' not modeled', LLT
         CALL LogIt(LOGSD,1,0)
      ENDIF
C
C
C     Lake St Clair stuff
      IF(LakesToSolve(3)) THEN
         WRITE(CTMP,9001) LakeName(5)(1:NLN(5))//' data', LLT
         CALL LogIt(LOGSD,1,0)
C        Figure out which output files were requested
         I = 0
         DO 564, K=1,4
            IF(WriteTSOutput(K, 1, 3)) THEN
               I = I + 1
               CWRK(I:I) = OutFileCodes(K:K)
            ENDIF
  564       CONTINUE
         IF( I .EQ. 0 ) THEN  
            WRITE(CTMP,9013) RightJustify( 'St. C. Output Files', 24), 
     &         RightJustify( 'none', 24), LLT
         ELSE
            WRITE(CTMP,9013) RightJustify( 'St. C. Output Files', 24), 
     &         RightJustify( CWRK(1:I), 24), LLT
         ENDIF
         CALL LogIt(LOGSD,1,0)
         WRITE(CTMP,9010) RightJustify( 'Det. Riv Eqn', 24), 
     &      SCK, SCYM, SCA, SCB, SCWT, SCC, LLT
         CALL LogIt(LOGSD,0,0)
         WRITE(CTMP,9007) RightJustify( 'St. C. Start Level', 24), 
     &      SCSTLV, UnitType(7)(1:UnitTypeLen(7)), LLT
         CALL LogIt(LOGSD,0,0)
         IF(NUMAVSZ(5) .EQ. 1) THEN
            WRITE(CTMP,9009) RightJustify( 'St. C. Area/Elev', 24), 
     &         AVSZ(5,1,2), UnitType(9)(1:UnitTypeLen(9)), LLT
            CALL LogIt(LOGSD,0,0)
         ELSE
            DO 2298, J=1,NUMAVSZ(5)
               WRITE(CTMP,9008) RightJustify( 'St. C. Area/Elev', 24), 
     &            AVSZ(5,J,2), UnitType(9)(1:UnitTypeLen(9)), 
     &            AVSZ(5,J,1), UnitType(7)(1:UnitTypeLen(7)), LLT
 2298          CALL LogIt(LOGSD,0,0)
         ENDIF
C        Determine if Time-Series Data file given
         I = INDEX(FileNames(32),LLT) + 1
         CWRK = FileNames(32)(I:I+6)
         IF(FileNames(32)(I:I+6) .EQ. 'default') THEN
C           Pre-written TSD file not given; had to specify Time-Series Data 
C           individually.  Print time-series data inputs:
C      DATA PARIDS(1409:1424) /'STCRIVFLOWDATA  '/                         28
            CALL UnReadInputData( 28, 
     &         RightJustify('St. C. Riv. Flow Data', 24), LOGSD)
C      DATA PARIDS(0993:1008) /'STCNBSDATA      '/                         63
            CALL UnReadInputData( 2, 
     &         RightJustify('St. C. NBS Data', 24), LOGSD)
C      DATA PARIDS(1009:1024) /'STCCONUSEDATA   '/                         64
            CALL UnReadInputData( 6, 
     &         RightJustify('St. C. Con. Use Data', 24), LOGSD)
C      DATA PARIDS(1025:1040) /'STCOTHERDIVDATA '/                         65
            CALL UnReadInputData( 14, 
     &         RightJustify('St. C. Other Div. Data', 24), LOGSD)
C      DATA PARIDS(1041:1056) /'STCOTHERSUPDATA '/                         66
            CALL UnReadInputData( 10, 
     &         RightJustify('St. C. Other Sup. Data', 24), LOGSD)
C      DATA PARIDS(1105:1120) /'STCEOPLEVELDATA '/                         70
            CALL UnReadInputData( 72, 
     &         RightJustify('St. C. EOP Level Data', 24), LOGSD)
C      DATA PARIDS(1505:1520) /'DETRIVICEDATA   '/                         95
            CALL UnReadInputData( 36, 
     &         RightJustify('Det. Riv. Ice Data', 24), LOGSD)
C      DATA PARIDS(1521:1536) /'DETRIVWEEDDATA  '/                         96
            CALL UnReadInputData( 37, 
     &         RightJustify('Det. Riv. Weed Data', 24), LOGSD)
C      DATA PARIDS(1537:1552) /'DETRIVRETRDATA  '/                         97
            CALL UnReadInputData( 17, 
     &         RightJustify('Det. Riv. Retr. Data', 24), LOGSD)
         ELSE
C           Pre-written TSD file was given; no other Time-Series Data given
            CALL UnReadInputData( 32, 
     &         RightJustify('St. C. Sum. Input File', 24), LOGSD)
         ENDIF
      ELSE
         WRITE(CTMP,9001) LakeName(5)(1:NLN(5))//' not modeled', LLT
         CALL LogIt(LOGSD,1,0)
      ENDIF
C
C     Erie
      IF(LakesToSolve(4)) THEN
         WRITE(CTMP,9001) LakeName(6)(1:NLN(6))//' data', LLT
         CALL LogIt(LOGSD,1,0)
C        Figure out which output files were requested
         I = 0
         DO 664, K=1,4
            IF(WriteTSOutput(K, 1, 4)) THEN
               I = I + 1
               CWRK(I:I) = OutFileCodes(K:K)
            ENDIF
  664       CONTINUE
         IF( I .EQ. 0 ) THEN  
            WRITE(CTMP,9013) RightJustify( 'Eri. Output Files', 24), 
     &         RightJustify( 'none', 24), LLT
            CALL LogIt(LOGSD,1,0)
         ELSE
            WRITE(CTMP,9013) RightJustify( 'Eri. Output Files', 24), 
     &         RightJustify( CWRK(1:I), 24), LLT
            CALL LogIt(LOGSD,1,0)
            IF( WriteErieOutflow ) THEN
               WRITE(CTMP,9003) RightJustify( 'Write Erie Outflow', 24), 
     &            WriteErieOutflow, LLT
            ELSE
               WRITE(CTMP,9003) RightJustify( 'Write Erie Outflow', 24), 
     &            .FALSE., LLT
            ENDIF
            CALL LogIt(LOGSD,0,0)
         ENDIF
         WRITE(CTMP,9010) RightJustify( 'Nia. Riv Eqn', 24), 
     &      ERK, ERYM, ERA, ERB, ERWT, ERC, LLT
         CALL LogIt(LOGSD,0,0)
         IF(IDNINT(CGSTLV) .NE. MissingVal) THEN
            WRITE(CTMP,9007) RightJustify( 'CGIP Start Level', 24), 
     &         CGSTLV, UnitType(7)(1:UnitTypeLen(7)), LLT
            CALL LogIt(LOGSD,0,0)
         ENDIF
         WRITE(CTMP,9007) RightJustify( 'Eri. Start Level', 24), 
     &      ERSTLV, UnitType(7)(1:UnitTypeLen(7)), LLT
         CALL LogIt(LOGSD,0,0)
         IF(NUMAVSZ(6) .EQ. 1) THEN
            WRITE(CTMP,9009) RightJustify( 'Eri. Area/Elev', 24), 
     &         AVSZ(6,1,2), UnitType(9)(1:UnitTypeLen(9)), LLT
            CALL LogIt(LOGSD,0,0)
         ELSE
            DO 2398, J=1,NUMAVSZ(6)
               WRITE(CTMP,9008) RightJustify( 'Eri. Area/Elev', 24), 
     &            AVSZ(6,J,2), UnitType(9)(1:UnitTypeLen(9)), 
     &            AVSZ(6,J,1), UnitType(7)(1:UnitTypeLen(7)), LLT
 2398          CALL LogIt(LOGSD,0,0)
         ENDIF
C        Determine if Time-Series Data file given
         I = INDEX(FileNames(33),LLT) + 1
         CWRK = FileNames(33)(I:I+6)
            IF(FileNames(33)(I:I+6) .EQ. 'default') THEN
C           Pre-written TSD file not given; had to specify Time-Series Data 
C           individually.  Print time-series data inputs:
            CALL UnReadInputData( 29,
     &         RightJustify('Det. Riv. Flow Data', 24), LOGSD)
            CALL UnReadInputData(  3, 
     &         RightJustify('Eri. NBS Data', 24), LOGSD)
            CALL UnReadInputData( 11, 
     &         RightJustify('Eri. Other Sup. Data', 24), LOGSD)
            CALL UnReadInputData(  7, 
     &         RightJustify('Eri. Con. Use Data', 24), LOGSD)
            CALL UnReadInputData( 21, 
     &         RightJustify('Welland Canal Data', 24), LOGSD)
            CALL UnReadInputData( 15, 
     &         RightJustify('Eri. Other Div. Data', 24), LOGSD)
            CALL UnReadInputData( 73, 
     &         RightJustify('Eri. EOP Level Data', 24), LOGSD)
            CALL UnReadInputData( 20, 
     &         RightJustify('CGIP Avg. Level Data', 24), LOGSD)
            CALL UnReadInputData( 18, 
     &         RightJustify('Nia. Riv. Retr. Data', 24), LOGSD)
            CALL UnReadInputData( 38, 
     &         RightJustify('Nia. Riv. Ice Data', 24), LOGSD)
            CALL UnReadInputData( 39, 
     &         RightJustify('Nia. Riv. Weed Data', 24), LOGSD)
         ELSE
C           Pre-written TSD file was given; no other Time-Series Data given
            CALL UnReadInputData( 33, 
     &         RightJustify('Eri. Sum. Input File', 24), LOGSD)
         ENDIF
      ELSE
         WRITE(CTMP,9001) LakeName(6)(1:NLN(6))//' not modeled', LLT
         CALL LogIt(LOGSD,1,0)
      ENDIF


C
      IF(LakesToSolve(5)) THEN
         WRITE(CTMP,9001) LakeName(7)(1:NLN(7))//' data', LLT
         CALL LogIt(LOGSD,1,0)
         WRITE(CTMP,9005) RightJustify( 'Ont. Verbosity', 24), 
     &      OntVerbosity, LLT
         CALL LogIt(LOGSD,0,0)
C        Figure out which output files were requested
         I = 0
         DO 764, K=1,4
            IF(WriteTSOutput(K, 1, 5)) THEN
               I = I + 1
               CWRK(I:I) = OutFileCodes(K:K)
            ENDIF
  764       CONTINUE
         IF( I .EQ. 0 ) THEN  
            WRITE(CTMP,9013) RightJustify( 'Ont. Output Files', 24), 
     &         RightJustify( 'none', 24), LLT
         ELSE
            WRITE(CTMP,9013) RightJustify( 'Ont. Output Files', 24), 
     &         RightJustify( CWRK(1:I), 24), LLT
         ENDIF
         CALL LogIt(LOGSD,1,0)
         WRITE(CTMP,9013) RightJustify( 'Ont. Time Step', 24), 
     &      RightJustify( IntervalName(RoutingInterval(3))
     &      (1:IntervalNameLen(RoutingInterval(3))), 24), LLT
         CALL LogIt(LOGSD,0,0)
         IF(NumIncsVary(3)) THEN
            WRITE(CTMP,9013) RightJustify( 'Ont. Inc.', 24), 
     &         RightJustify( 'Hourly', 24), LLT
         ELSE
            WRITE(CTMP,9005) RightJustify( 'Ont. Inc.', 24), 
     &         RoutingIncrement(3), LLT
         ENDIF
         CALL LogIt(LOGSD,0,0)
         WRITE(CTMP,9007) RightJustify( 'Ont. Start Level', 24), 
     &      ONSTLV, UnitType(7)(1:UnitTypeLen(7)), LLT
         CALL LogIt(LOGSD,0,0)
         IF(NUMAVSZ(7) .EQ. 1) THEN
            WRITE(CTMP,9009) RightJustify( 'Ont. Area/Elev', 24), 
     &         AVSZ(7,1,2), UnitType(9)(1:UnitTypeLen(9)), LLT
            CALL LogIt(LOGSD,0,0)
         ELSE
            DO 2498, J=1,NUMAVSZ(7)
               WRITE(CTMP,9008) RightJustify( 'Ont. Area/Elev', 24), 
     &            AVSZ(7,J,2), UnitType(9)(1:UnitTypeLen(9)), 
     &            AVSZ(7,J,1), UnitType(7)(1:UnitTypeLen(7)), LLT
 2498          CALL LogIt(LOGSD,0,0)
         ENDIF
c
c--- Write statements for the SUMMARY log file
c
         IF(UStLawrSolution .EQ. 1) THEN
            WRITE(CTMP,9002)RightJustify( 'U. St. L. Riv Solution',24), 
     &         RightJustify(SMSolutionName(1)(1:16), 24), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9010) RightJustify( 'St. L. Riv Eqn', 24), 
     &         ONK, ONYM, ONA, ONB, ONWT, ONC, LLT
            CALL LogIt(LOGSD,0,0)
         ELSEIF(UStLawrSolution .EQ. 2) THEN
            WRITE(CTMP,9001) 'Plan 1958D Regulation Parameters', LLT
            CALL LogIt(LOGSD,1,0)
            WRITE(CTMP,9002)RightJustify( 'U. St. L. Riv Solution',24), 
     &         RightJustify('1958D',24), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9003) RightJustify( 'P58D Use Database', 24), 
     &         P58DUseDatabase, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9003) RightJustify( 'P58D SI Database', 24), 
     &         P58DSiDatabase, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9005) RightJustify( 'P58D Crit. Depths', 24), 
     &         CritDepths, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9003) RightJustify( 'P58D Eliminate Dev.', 24), 
     &         ElimDev, LLT
            IF(ElimDev) THEN
               WRITE(CTMP,9002)RightJustify( 'Criterion K File', 24), 
     &             CritKFile(1:LastKC(CritKFile,256)), LLT
            CALL LogIt(LOGSD,0,0)
            ENDIF
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9003) RightJustify( 'P58D Force Flow', 24), 
     &         ForceFlow, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9003) RightJustify( 'P58D Set Plan Flow', 24), 
     &         SetFlow, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9003) RightJustify( 'P58D Compute Montreal', 24), 
     &         ComputeMontreal, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9003) RightJustify( 'P58D Pre-Proj. Cond.', 24), 
     &         PreProjConditions, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9011) RightJustify( 'P58D Total Supply', 24), 
     &         P58DTotalSup, UnitType(2)(1:UnitTypeLen(2)), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9011) RightJustify( 'P58D Weighted Supply', 24), 
     &         P58DWeightedSupply, UnitType(2)(1:UnitTypeLen(2)), LLT
            CALL LogIt(LOGSD,0,0)
C
            WRITE(CTMP,9011) RightJustify('P58D Supply Adjustment',24), 
     &         P58DAdjustment, UnitType(2)(1:UnitTypeLen(2)), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9011) RightJustify( 'P58D RC Adjustment', 24), 
     &         P58DrcAdjustment, UnitType(2)(1:UnitTypeLen(2)), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9011) RightJustify( 'P58D Rule Curve', 24), 
     &         P58DRuleCurve, UnitType(2)(1:UnitTypeLen(2)), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9011) RightJustify( 'P58D App Limit', 24), 
     &         P58DAppLimit, UnitType(2)(1:UnitTypeLen(2)), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9011) RightJustify( 'P58D Actual Flow', 24), 
     &         P58DActualFlow, UnitType(2)(1:UnitTypeLen(2)), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9007) RightJustify( 'P58D Plan Level', 24), 
     &         P58DPlanLevel, UnitType(7)(1:UnitTypeLen(7)), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9007) RightJustify( 'P58D Actual Level', 24), 
     &         P58DActualLevel, UnitType(7)(1:UnitTypeLen(7)), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9011) RightJustify( 'P58D Total Deviations',24), 
     &         P58Dtotaldev, UnitType(2)(1:UnitTypeLen(2)), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9011) RightJustify('P58D LSL Flow',24), 
     &         P58DLSLflow, UnitType(2)(1:UnitTypeLen(2)), LLT
            CALL LogIt(LOGSD,0,0)
C
            WRITE(CTMP,9011) RightJustify( 'P58D Max. I Limit', 24), 
     &         MaxILimit, UnitType(2)(1:UnitTypeLen(2)), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9011) RightJustify( 'P58D Min. P Limit', 24), 
     &         MinPLimit, UnitType(2)(1:UnitTypeLen(2)), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9005) RightJustify( 'P58D Run Type', 24), 
     &         P58DRunType, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9015) RightJustify('P58D Sup. Indicator', 24), 
     &         (SupplyIndicator(J), J=2,21), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9006) RightJustify( 'P58D Forecast Date', 24), 
     &         P58DForecastYear, P58DForecastMon, P58DForecastDay, LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9002)RightJustify( 'P58D Simulate File', 24), 
     &         P58DSIMULATEFILE(1:LastKC(P58DSIMULATEFILE,79)), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9002)RightJustify( 'P58D App Limit Sym', 24), 
     &         P58DAppLimitSym(1:LASTKC(P58DAppLimitSym,2)), LLT
            CALL LogIt(LOGSD,0,0)
         ELSEIF(UStLawrSolution .EQ. 5) THEN
c
c     Preproject Conditions
c
            WRITE(CTMP,9001) 'Preproject Conditions Parameters', LLT
            CALL LogIt(LOGSD,1,0)
            WRITE(CTMP,9002)RightJustify( 'U. St. L. Riv Solution',24), 
     &         RightJustify('Preproject',24), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9007) RightJustify( 'Preproject Level', 24), 
     &         PreprojectLevel, UnitType(7)(1:UnitTypeLen(7)), LLT
            CALL LogIt(LOGSD,0,0)
            WRITE(CTMP,9003) RightJustify( 'P58D Compute Montreal', 24), 
     &         ComputeMontreal, LLT
            CALL LogIt(LOGSD,0,0)
         ENDIF
C
C        Determine if Time-Series Data file given
         I = INDEX(FileNames(57),LLT) + 1
         CWRK = FileNames(57)(I:I+6)
         IF(FileNames(57)(I:I+6) .EQ. 'default') THEN
C           Pre-written TSD file not given; had to specify Time-Series Data 
C           individually.  Print time-series data inputs:
            CALL UnReadInputData( 30,
     &         RightJustify('Nia. Riv. Flow Data', 24), LOGSD)
C           If Lake Erie solved, the Welland Canal already specified
            IF(.NOT. LakesToSolve(4)) CALL UnReadInputData( 21, 
     &         RightJustify('Welland Canal Data', 24), LOGSD)
            CALL UnReadInputData( 64, 
     &         RightJustify('Ont. Deviation Data', 24), LOGSD)
            CALL UnReadInputData( 53, 
     &         RightJustify('Ont. NBS Data', 24), LOGSD)
            CALL UnReadInputData( 56, 
     &         RightJustify('Ont. Other Sup. Data', 24), LOGSD)
            CALL UnReadInputData( 54, 
     &         RightJustify('Ont. Con. Use Data', 24), LOGSD)
            CALL UnReadInputData( 55, 
     &         RightJustify('Ont. Other Div. Data', 24), LOGSD)
            CALL UnReadInputData( 74, 
     &         RightJustify('Ont. EOP Level Data', 24), LOGSD)
            CALL UnReadInputData( 59, 
     &         RightJustify('St. L. Riv. Retr. Data', 24), LOGSD)
            CALL UnReadInputData( 60, 
     &         RightJustify('St. L. Riv. Ice Data', 24), LOGSD)
            CALL UnReadInputData( 61, 
     &         RightJustify('St. L. Riv. Weed Data', 24), LOGSD)
            CALL UnReadInputData( 67, 
     &         RightJustify('Ont. Gates Data', 24), LOGSD)
            CALL UnReadInputData( 75, 
     &         RightJustify('Ont. Force Flow Data', 24), LOGSD)
            CALL UnReadInputData( 58, 
     &         RightJustify('St. L. Riv. Flow Data', 24), LOGSD)
            CALL UnReadInputData( 66, 
     &         RightJustify('St. Louis Flow Data', 24), LOGSD)
            CALL UnReadInputData( 62, 
     &         RightJustify('St. Louis Retr. Data', 24), LOGSD)
            CALL UnReadInputData( 63, 
     &         RightJustify('L.S. Roughness Data', 24), LOGSD)
            CALL UnReadInputData( 65, 
     &         RightJustify('D. Pr. Riv. Flow Data', 24), LOGSD)
            CALL UnReadInputData( 77, 
     &         RightJustify('M. I. I. Riv. Flow Data', 24), LOGSD)
            CALL UnReadInputData( 76, 
     &         RightJustify('Louis-Ont. Diff. Data', 24), LOGSD)
         ELSE
C           Pre-written TSD file was given; no other Time-Series Data given
            CALL UnReadInputData( 57, 
     &         RightJustify('Ont. Sum. Input File', 24), LOGSD)
         ENDIF
      ELSE
         WRITE(CTMP,9001) LakeName(7)(1:NLN(7))//' not modeled', LLT
         CALL LogIt(LOGSD,1,0)
      ENDIF

C
 9000 FORMAT ('#', 4X, A, A1)
 9001 FORMAT ('#', '    ........................    ', A,
     &             '    ........................    ', A1)
 9002 FORMAT (A,': ',A,A1)
 9003 FORMAT (A,': ',23X,L1,A1)
 9004 FORMAT (A1,16X,'Title',I2,': ',A,A1)
 9005 FORMAT (A,': ',I24,A1)
 9006 FORMAT (A,': ',12X,I5,',',I3,',',I3,A1)
 9007 FORMAT (A,': ',F24.7,1X,A,A1)
 9008 FORMAT (A,': ',F24.3,1X,A,F18.3,1X,A,A1)
 9009 FORMAT (A,': ',F24.3,1X,A,A1)
 9010 FORMAT (A,': ',2(F10.5,','),3(F8.5,','),F7.3,A1)
 9011 FORMAT (A,': ',I24,1X,A,A1)
 9012 FORMAT (A,': ',F24.5,A1)
 9013 FORMAT (A,': ',A24,A1)
 9014 FORMAT (A,': ',11(I5,','),I5,A1)
 9015 FORMAT (A,': ',19(I4,','),I4,A1)
C
C     Write a little separator block
      DO 9711, K=1,79
 9711    CWRK(K:K)='#'
      WRITE(CTMP,9000) CWRK(1:79), LLT
      DO 9712, K=1,5
 9712    CALL LogIt(LOGSD,0,0)

      END
