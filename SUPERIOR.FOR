      SUBROUTINE RouteSuperior(NewMonth, SPSolveErrors)
C      
C     This is the main Lake Superior module.  This is where we decide
C     what regulation plan to use, how to handle forced flows or 
C     gates, how to solve for EOP level, etc.  This is a somewhat less
C     dismal as of Dec08 mess.  
C
      IMPLICIT NONE 
C
C
      LOGICAL       NewMonth
C                   .TRUE. if first time routing through the month
C
      INTEGER       SPSolveErrors       
C                   Counts how often the simultaneous solution of the
C                   routing equations failed to achieve a certain
C                   accuracy.
C
C     Declarations regarding local variables or functions:
C
      INTEGER       I, K
      DOUBLE PRECISION TmpA
C                   Scratch variables
C
      CHARACTER     CHNKST*16
C                   This is the calling subroutine name passed to
C                   ErrMess() in case of error messages.
C
      INTEGER       MaxIterations
C
      DOUBLE PRECISION Tolerance
C                   Threshold at which the difference in outflow between
C                   two iterations becomes close enough, in m3s.
C                    
C
      DOUBLE PRECISION CWFreeBoard
C                   and the freeboard available for a given water
C                   level (assumed to be SW Pier).
C
      DOUBLE PRECISION SideChanMaxFlow
C                   lesser of physcial max capacity from Jacob vs user-defined maximum
C
      DOUBLE PRECISION SupIndex1, SupIndex2, SupIndex3
C                   time-series of arbitrary indices used for Superior 
C                   regulation studies
C      
      REAL          PreviousGates
C
      INTEGER       SecsInMonth
C
      INTEGER       SPNWB, MHNWB, SCNWB, ERNWB, MHRet, SCRet, ERRet
C                   Values that hold routing data appropriate to setting
C                   of P77ACompat.  Not used after overtopping calls.
C
      DOUBLE PRECISION SPIMLV, SWPI, GateFlowI
C                   Mean levels of Lake Superior and SW Pier for an increment,
C                   mean gated flow for an increment
C
      DOUBLE PRECISION LastSPEOI, LevelToInteger, SPErrorFlow
      INTEGER       EOILevels2Ago
C                   Used for tracking iterative convergence
C
      DOUBLE PRECISION SWPier, GateFlow
C                   Level at SW Pier, and flow through gates.
C
      DOUBLE PRECISION OUTFLOW
C                   Returns connecting channel flow, given the
C                   parameters for the discharge equation.
C
      DOUBLE PRECISION CompWorksFreeboard
C                   Returns amount of clearance before overtopping
C                   the lowest point of the Compensating Works
C
      DOUBLE PRECISION GETAREA
C                   Function that returns surface area for a lake
C
      DOUBLE PRECISION ROUND
C                   Returns value rounded to a certain accuracy,
C                   according to the old Lake Survey method.
C
      CHARACTER*1   EOPType
C                	  Returns 'M' if period ends a month, or 'Q' if period
C                   ends a quarter-month, else returns 'W', or 'D'.
C
C
C     Declarations regarding COMMON block variables:
C
C
      REAL          SupGatesOpen
C                   The number of Superior gates open in current routing 
C                   period.
C
      INTEGER       StMarysSolution, P77ACompat
C                   1=outflow eqn, 2=1955 mod rule, 3=plan77a
C                   -1 = virtcarl, 0 = read-in, 1 = 6-month operational
C
      INTEGER       PeriodEnd(6,3)
C                   Tells the date that the current routing period ends
C                   on for each routing interval.  PeriodEnd(1,?) gives
C                   the day, PeriodEnd(2,?) gives the quartermonth number, 
C                   PeriodEnd(3,?) gives the month, PeriodEnd(4,?) gives 
C                   the year, PeriodEnd(5,?) gives the number of elapsed
C                   days at the final day of the interval.  PeriodEnd(6,?)
C                   gives the elapsed days for the end of the next period 
C                   (for loop control).  PeriodEnd(?,1) is for the Lake 
C                   Superior routing interval, PeriodEnd(?,2) is for the
C                   middle lakes, and PeriodEnd(?,3) is for Lake Ontario.
C
      INTEGER       Day(3), Q(4,3)
C                   Start date (of month) for the period being routed.
C                   Q is an array holding the end date of each quarter-
C                   month within the month (so Q(4,L) tells month length
C                   for lake L).  
C
      INTEGER       DaysSinceSaturday
C                   See description in subroutine RUNIT
C
      DOUBLE PRECISION dT(3)
C                   "delta time" associated with routing interval (or
C                   increment).  Used in storage-change calculations.
C                   First element is for Lake Superior, second for
C                   middle lakes, third for Lake Ontario.
C
      LOGICAL       NumIncsVary(3) 
C                   True if number of increments gets adjusted on the fly 
C                   during monthly or qmonthly routing, using actual
C                   month length, to make dT come out to 1 hour. 
C
      CHARACTER     SDATE*9,STIME*8
C                   See descriptions in INITPARAMS()
C
      DOUBLE PRECISION SPK,SPYM,SPA,SPB,SPWT,SPC
C                   See descriptions in CHECKPAR and OUTFLOW
C
      INTEGER       SideChanMaxCapacity(12), SideChanMinCapacity(12)
      DOUBLE PRECISION CWOvertopElev, CWMinimumFreeboard(12)
C                   Maximum amount of Superior outflow that does not pass
C                   through the Compensating Works gates specified for each 
C                   month; elevation at which flow begins going over the 
C                   lowest gates of the Compensating Works; and the minimum 
C                   safe height between the SW Pier level and CWOvertopElev 
C                   (specified for each month). 
C
      INTEGER       RoutingInterval(3), RoutingIncrement(3)
C                   Routing interval to use on Superior, middle lakes,
C                   and Ontario (1=daily, 2=weekly, 3=qmonthly, 
C                   4=monthly), and the number of increments within each 
C                   interval.
C
      CHARACTER*9   MLSolutionName(2)
C                   Strings set to 'interative' and 'matrix', used for 
C                   output.
C
      CHARACTER*8   IntervalName(4)
C                   an 8 element char array containign names of the
C                   routing  intervals supported by the model, for use
C                   in output.
C
      INTEGER       IntervalNameLen(4)
C                   Number of letters in each element of IntervalName
C
      CHARACTER*16  SMSolutionName(4), USLSolutionName(5)
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      INTEGER       GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &              SupVerbosity, MLVerbosity, OntVerbosity
      INTEGER       MonthLength
      CHARACTER*8   MonthLengthName(2)
      CHARACTER*8   RCRUNTYP
      CHARACTER*256  FileNames(100)
      CHARACTER*256 RepeatLines(80)
      CHARACTER*21  LAKENAME(16), RiverName(16)
      INTEGER       NLN(16), LkToSlv2LN(8), NRN(16), LkSlv2RtLoopI(8)
      INTEGER       TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &              MissingVal, SecsInUniformMon
      DOUBLE PRECISION MPFT, HALF, ONE, ZERO
      LOGICAL       LakesToSolve(8) 
      INTEGER       StartYear, StartMonth, StartDay, 
     &              EndYear, EndMonth, EndDay
      INTEGER       MidLakeSolution, UStLawrSolution
      INTEGER       BOMRound, FlowRound, LevelRound
      CHARACTER     CTMP*255, CWRK*80
      CHARACTER     LLT*1
      CHARACTER*3   MonthAbbrev(12)
      DOUBLE PRECISION SGX(17), SGY(17), SGZ(17)
      INTEGER       MaxSGSetting
      DOUBLE PRECISION SupPlanFlow_Dble
C                   This variable is used in rounding the Plan outflow 
C                   for the output file (see ~1380).
C
C                   See descriptions in Block Data Genral
C     
      INTEGER       SPOPNBS06(12,6), SPOPNBS07(12,6), OPLLO(12), 
     &              MHOPNBS06(12,6), MHOPNBS07(12,6), OPChicago, 
     &              SCRRet(12), SCOPNBS06(12,6), SCOPNBS07(12,6), 
     &              DetRet(12), EROPNBS06(12,6), EROPNBS07(12,6),
     &              OPWelland(12), NiaRet(12),
     &              SCRRet_P2012(12), SCOPNBS08(12,6), SPOPNBS08(12,6), 
     &              DetRet_P2012(12), EROPNBS08(12,6), MHOPNBS08(12,6),
     &              OPWelland_P2012(12), NiaRet_P2012(12)       
C                   See descriptions in Block Data Genral
C
      INTEGER       ForecastMonths, LStoSWP, SupCritRnd, MonthsChecked
C                   Which LS:SWP relationship to use, and rounding
C                   practice to use for Criterion B and C computations,
C                   and which month of OpSupply sequence we are in if
C                   compat=2. 
C
      LOGICAL       BackCompat
C                   .TRUE. if program is trying to reproduce results from
C                   older upper lakes models.
C
      CHARACTER*40  SPRegCodes
C                   Describes which Plan 1977A Requirements and 
C                   Criteria come into play, used for output.
C
      CHARACTER*120 P77AHeader
C                   Column descriptions for Superior output
C
      CHARACTER*120 P2012Header
C                   Column descriptions for Plan 2012 Superior output
C
      CHARACTER*8   BalancingFlow
      CHARACTER*8   SupBOMLevel 
C                   The flow called for by the balancing equation,
C                   and beginning-of-month level.  Used for output.
C
      CHARACTER     SupPreProject*6, SupGateFlow*6,
     &              SupSideChanFlow*6, SupPlanFlow*6, SWPierLevel*8
C                   Lake Superior regulation output strings
C
      DOUBLE PRECISION FWRK
C                   This is a common block used for 64-bit calculations.
      DIMENSION FWRK(225,12)
C                   Dimensioned it this way to provide a 225 element
C                   array for each pool.
      DOUBLE PRECISION SPBOP, SPEOP, SPBOI, SPEOI, SPNTS, SMICE, SPMLV, 
     &              MHBOP, SCBOP,ERBOP, STMRIVFLOW, SPOutI, SolutionTest
C                   Beginning of period level, end of period level,
C                   beginning of increment level, end of increment 
C                   level, net total supply, ice/weed retardation,
C                   mean level for period, mean flow in the St. Marys
C                   River for the period, mean flow in the St. Marys 
C                   River for the increment.
C                   SolutionTest is a threshold for flagging whether
C                   the accuracy of a solution is poor enough to flag
C                   with a warning message.  
C                   Note that these storages are actually aliased to
C                   memory in FWRK.
C
C
      DOUBLE PRECISION MHEndofMonthL, SCEndofMonthL, EREndofMonthL,
     &              SPEndofMonthL, SPMonthlyMeanF
C                   The end-of-month levels and mean monthly flow, 
C                   computed by TallyMeans, that get used as beginning
C                   of month levels and previous month's outflow by 
C                   Superior regulation plans.              
C                   Note that these storages are actually aliased to
C                   memory in FWRK.
C
      DOUBLE PRECISION NominalStMarys
C                   Nominal Lake Superior outflow for the regulation 
C                   period, as determined by a regulation plan (Plan
C                   1977A, 1955 Modified Rule of '49).  Note that a
C                   regulation plan may have a different period than
C                   the routing interval.  For instance, you could
C                   route daily or quarter-monthly NBS through Lake
C                   Superior in conjunction with Plan 1977A-determined
C                   monthly outflows.  The outflow remains constant 
C                   for each month, but the level will fluctuate a
C                   little each day or quarter-month.  Note that the
C                   model completely ignores the bigger question of 
C                   whether this type of calculation has any validity.
C                   We modellers are supposed to know if what we are 
C                   doing makes sense; after all, that's why they pay
C                   us the big bucks (ROTFL).  Note that these storages 
C                   are actually aliased to memory in FWRK.
C
      DOUBLE PRECISION SupForceFlow
C                   The St. Marys River flow to use if solving for
C                   Lake Superior levels, but forcing an outflow.
C
      DOUBLE PRECISION SupSideChannel, ComputedNonGate,
     &                 LastComputedNonGate
C                   The combined flow NOT going through the Compensating
C                   Works.  ComputedNonGate is the value determined by  
C                   the regulation plan, and the value persists until the 
C                   regulation plan cycles again (i.e., next month).  
C                   SupSideChannel starts out as the input read from the
C                   time series data.  If the input is -9999, the value
C                   of SupSideChannel is taken as ComputedNonGate.  If
C                   SupSideChannel comes in as anything else, it retains 
C                   that value as a forced amount. If solving for Lake 
C                   Superior levels using an outflow equation, neither is
C                   used.   LastComputedNonGate stores the value of 
C                   ComputedNonGate in case the latter is overriden by a
C                   forced sidechannel.
C
      DOUBLE PRECISION SupDevFlow, SupForceEOPLev
C                   Deviations from Lake Superior regulated discharge,
C                   and optional level to replace the modeled BOP level.
C                   
      DOUBLE PRECISION NumSupForceGates
C                   Variables used to store the number of gates open 
C                   at the Compensating Works if gate settings forced.
C
      DOUBLE PRECISION SPSolutionCheck, LastFlow, AdjOTFlow, CBUSSLevel
C                   Scratch variables.  SPSolutionCheck contains the 
C                   the accuracy of superior routing (under some 
C                   methods).  LastFlow is used for convergence
C                   testing during iterative solutions. 
C                   Note that these storages are actually aliased to
C                   memory in FWRK.
C
      INTEGER P77ASup05NBS(12), P77ASup50NBS(12), P77ASup95NBS(12), 
     &        P77ALLOgokiDiv(12),P77AStMarysRetrd(12), 
     &        P77AMH05NBS(12), P77AMH50NBS(12), P77AMH95NBS(12), 
     &        P77AChicagoDiv(12), P77AStClairRetrd(12), 
     &        P77AStC50NBS(12), P77ADetroitRetrd(12),
     &        P77AErie50NBS(12), P77AWellandFlow(12), 
     &        P77ANiagaraRetrd(12), P77ABaseFlow(12),
     &        P77AMidLakesInc, P77AMidLakesTol, P77AWinterMax,
     &        P77AMonthRound, P77AFlowRound, P77ALevelRound,
     &        P77AMonthsToAve(12), P77AMaxChange, P77AMinimumFlow(12)
      DOUBLE PRECISION  P77ACGIP, P77ALowTrigger, P77AHighTrigger, 
     &        P77ASPK, P77ASPYM, P77ASPA, P77ASPB, P77ASPWT, P77ASPC, 
     &        P77AMHK, P77AMHYM, P77AMHA, P77AMHB, P77AMHWT, P77AMHC, 
     &        P77ASCK, P77ASCYM, P77ASCA, P77ASCB, P77ASCWT, P77ASCC,
     &        P77AERK, P77AERYM, P77AERA, P77AERB, P77AERWT, P77AERC,
     &        P77ASupTargetLev(12), P77ASupRangeFact(12),
     &        P77AMHTargetLev(12), P77AMHRangeFact(12), P77AAValue
      LOGICAL P77AUniformMonth 
C             See descriptions in SUBROUTINE Plan1977A()
      DOUBLE PRECISION P2012SPTARGET(12), P2012MHTARGET(12),
     &                 P2012SPHIGHLEVEL, P2012SPLOWLEVEL,
     &                 P2012MHTARGETADJ
      DOUBLE PRECISION P2012PARA(6) 
      DOUBLE PRECISION ForcedSupGates(12)
      INTEGER          IndexGatesSide, ForcedSideChan(12)        
      DOUBLE PRECISION PrevSupMean
C                  Used to store SPMLV, primarily for P55MR49.
C       
C  Common block statements.
C
      COMMON /P2012/    P2012PARA, P2012SPTARGET, P2012MHTARGET, 
     &                  P2012SPHIGHLEVEL, P2012SPLOWLEVEL,
     &                  P2012MHTARGETADJ
      COMMON /FORGATE/  ForcedSupGates   
      COMMON /FORSIDE/  ForcedSideChan 
      COMMON /INFORCE/  IndexGatesSide   
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
      COMMON /CLDATE/   dT, Day, Q, PeriodEnd, NumIncsVary,
     &                  DaysSinceSaturday
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /STIME/    STIME,SDATE
      COMMON /SPDAT/    SPK, SPYM, SPA, SPB, SPWT, SPC,
     &                  SGX, SGY, SGZ, MaxSGSetting,
     &                  SupGatesOpen
      COMMON /RCRUNI/   LakesToSolve, 
     &                  IntervalNameLen, MonthLength,
     &                  StartYear, StartMonth, StartDay, 
     &                  EndYear, EndMonth, EndDay
      COMMON /RCRUNC/   RCRUNTYP, IntervalName, MonthLengthName,
     &                  MLSolutionName, SMSolutionName, USLSolutionName
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /MISINT/   MidLakeSolution, BOMRound, FlowRound, 
     &                  LevelRound, RoutingInterval, RoutingIncrement,
     &                  StMarysSolution, UStLawrSolution, P77ACompat
      COMMON /FNAMEC/   FileNames, RepeatLines
      COMMON /FSPACE/   FWRK
      COMMON /LNAMEC/   LAKENAME, RiverName
      COMMON /LNAMEI/   NLN, LkToSlv2LN, NRN, LkSlv2RtLoopI
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
      COMMON /SUPOUT/   SupPreProject, BalancingFlow, SupGateFlow,
     &                  SupPlanFlow, SupBOMLevel, SWPierLevel, 
     &                  SupSideChanFlow, SPRegCodes, P77AHeader,
     &                  P2012Header
      COMMON /MONAME/   MonthAbbrev
      COMMON /SUPCOM/   SPOPNBS06, OPLLO, MHOPNBS06, OPChicago, SCRRet,
     &                  SCOPNBS06, DetRet, EROPNBS06, OPWelland, NiaRet,
     &                  SPOPNBS07, MHOPNBS07, SCOPNBS07, EROPNBS07,
     &                  SPOPNBS08, MHOPNBS08, SCOPNBS08, EROPNBS08,
     &                  OPWelland_P2012, SCRRet_P2012, DetRet_P2012,
     &                  NiaRet_P2012,     
     &                  ForecastMonths,
     &                  LStoSWP, SupCritRnd, BackCompat, MonthsChecked
      COMMON /SUPIDX/   SupIndex1, SupIndex2, SupIndex3
      COMMON /SPOTOP/   SideChanMaxCapacity, CWOvertopElev, 
     &                  CWMinimumFreeboard, SideChanMinCapacity
      COMMON /SPROUT/   LastComputedNonGate, SPNWB
C                       Put these in a COMMON block in order to preserve
C                       the value between calls to this SUBROUTINE 
C                       (without using FWRK storage).
C
C
      EQUIVALENCE (SPNTS, FWRK(1,5))
      EQUIVALENCE (FWRK(2,5), SMICE)
      EQUIVALENCE (FWRK(3,5), SPBOI)
      EQUIVALENCE (FWRK(4,5), SPEOI)
      EQUIVALENCE (SPMLV, FWRK(5,5))
      EQUIVALENCE (STMRIVFLOW, FWRK(6,5))
      EQUIVALENCE (FWRK(7,5), SPBOP)
      EQUIVALENCE (FWRK(7,1), MHBOP)
      EQUIVALENCE (FWRK(7,2), SCBOP)
      EQUIVALENCE (FWRK(7,3), ERBOP)
      EQUIVALENCE (SPErrorFlow, FWRK(7,4))
      EQUIVALENCE (FWRK(8,5), SPEOP)
      EQUIVALENCE (FWRK(9,4), SolutionTest)
      EQUIVALENCE (FWRK(9,5), SPSolutionCheck)
C
C     FWRK(10,1-4), FWRK(11,1-4), FWRK(12,1-4), FWRK(13,1-4), and 
C     FWRK(14,1-4) get overwritten each increment in Middle Lakes
C     routing.  Use them for intermediate values if desired.
C     
      EQUIVALENCE (FWRK(10,5), LastFlow)
      EQUIVALENCE (FWRK(11,4), SupForceFlow)
      EQUIVALENCE (FWRK(11,5), NominalStMarys)
      EQUIVALENCE (FWRK(11,6), SupDevFlow)
      EQUIVALENCE (FWRK(12,1), SWPI)
      EQUIVALENCE (FWRK(12,2), GateFlowI)
      EQUIVALENCE (FWRK(12,5), ComputedNonGate)
      EQUIVALENCE (FWRK(12,6), SupSideChannel)
C
      EQUIVALENCE (FWRK(13,5), SPIMLV)
      EQUIVALENCE (FWRK(14,5), LastSPEOI)
      EQUIVALENCE (FWRK(15,4), LevelToInteger) 
      EQUIVALENCE (FWRK(15,5), SPOutI)
      EQUIVALENCE (FWRK(16,2), NumSupForceGates)
      EQUIVALENCE (FWRK(22,5), SPMonthlyMeanF)
      EQUIVALENCE (FWRK(30,1), MHEndofMonthL)
      EQUIVALENCE (FWRK(30,2), SCEndofMonthL)
      EQUIVALENCE (FWRK(30,3), EREndofMonthL)
      EQUIVALENCE (FWRK(30,5), SPEndofMonthL)
      EQUIVALENCE (FWRK(31,5), SupForceEOPLev)
C
C      PARAMETER (SumTitleStr = 'Summary Output File')
C
      DATA MaxIterations /9999/
C
C     Start with  value of SecsInUniformMon.  If using actual month
C     length, will get overridden each month.
      DATA SecsInMonth / 2629746 /
C
C
      CHNKST = 'RouteSuperior'//LLT
C
C
C     Consider equations solved if iterate to one FlowRound (ie, if
C     if FlowRound set to 1, must be 10 m3s or less difference between
C     iterations).  If FlowRound == -99, use .1 m3s.
      Tolerance = 10.0D0**FlowRound
      IF( Tolerance .LT. 0.1D0 ) Tolerance = 0.1D0
C
C     Initialize BOP levels as the EOP levels from previous period
      SPBOP = SPEOP
C
C     Clear the field for codes for this routing interval
      SPRegCodes = '                                        '
C
C     Let's save ourselves some grief by checking a little bit here:
      IF( IDNINT(SPBOP) .LE. ZERO .OR. SPBOP .GT. 1000) THEN
         WRITE(CTMP,4289) 'Superior BOP level', PeriodEnd(1,1),
     &       MonthAbbrev(PeriodEnd(3,1)), PeriodEnd(4,1), SPBOP, LLT
         CALL LogIt(LOGCDE,1,0)
         CALL ErrMess(1141,'SuperiorRoute'//LLT//'  ',LOGCDE)
         CALL SAVETOD
         CALL ADIOS()
      ENDIF
 4289 FORMAT(' Unbelievable ', A, ' for period ending ', I2, A3, I5,
     &   ': ', F17.5, A1)
C
C     Store Superior's mean level for use in the next month's regulation
C     (by P55MR49)
      PrevSupMean = SPMLV
C     Initialize mean lake level registers for period.  See 
C     function RouteMiddleLakes() for an explanation of the 
C     arithmetic. 
      SPMLV = SPBOP * HALF
C
C     Initialize EOI values with BOP values; they will be rolled 
C     into BOI values at beginning of loop through increments.
      SPEOI = SPBOP
C     write(*, '(A, 2F21.13, I9)') '  LevelToInteger=', LevelToInteger, 
C    &   SPEOI*LevelToInteger, IDNINT(SPEOI*LevelToInteger) 
C
C     Initialize PreviousGates with whatever left in SupGatesOpen, and
C     reset SupGatesOpen
      PreviousGates = SupGatesOpen
      SupGatesOpen = -1.0
C
C     Initialize period-average discharge values with 0.0D.  They 
C     will be computed as the mean of the outflows computed for 
C     each increment.  (The outflow for each increment is 
C     calculated using the mean level for the increment).
      STMRIVFLOW = ZERO
C
      IF(SupVerbosity .GE. 7 ) THEN
         WRITE(CTMP,'(2X,3(A,F15.8),A,F15.3,A,F9.3,A,F9.3,A1)') 
     &      ' SPBOP=', SPBOP, ' Superior NTS=', SPNTS,
     &      ' Ice=', SMICE, ' dT=', dT(1), ' SupGatesOpen=', 
     &      SupGatesOpen, ' ComputedNonGate=',ComputedNonGate, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,'(2X,5(A,F15.8),A1)')' SideChannel=',SupSideChannel,
     &      ' Forced Flow=', SupForceFlow, ' Forced Gates=', 
     &      NumSupForceGates, ' Deviation Flow=', SupDevFlow,
     &      ' SupForceEOPLev=', SupForceEOPLev, LLT
         CALL LogIt(LOGD,0,0)
         IF(SupVerbosity .GE. 9 ) THEN
            WRITE(CTMP,'(2X,3(A,F15.8),A1)') 
     &         ' SPEOP=', SPEOP, ' SPMLV=', SPMLV, ' StMRivFlow=', 
     &         StMRivFlow, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
      ENDIF
C
      IF (StMarysSolution .GE. 2 ) THEN
C        Numerous Lake Superior reg plans need to know BOM levels and 
C        the length of the month.  Just figure them out here whether 
C        really needed or not to avoid duplicate effort later.
C        
         IF(NewMonth) THEN
            IF(MonthLength .EQ. 2)  THEN
C              Actual Month length.  (Don't need number of days in year, 
C              so just use K as scratch)
               CALL NUMDAYS(PeriodEnd(4,1), PeriodEnd(3,1), K, I)
               SecsInMonth = I * SecsInDay
            ENDIF
C        
C
C           Write interesting data to character buffers for later output
C
            WRITE(SupBOMLevel,'(F8.3)') Round(SPEndofMonthL, -3)

C           TallyMeans populates SPEndofMonthL (i.e., FWRK(30,5)) and
C           similar values for the middle lakes.  For daily, qmonthly,
C           or monthly routing we can just use the value available there
C           when we hit a new month.  The weekly boundary will not quite
C           match up most of the time, and needs to be adjusted.
C           XX Monthly-Weekly Boundary Fudge:
C           XX If its weekly routing, just use level at end of last
C           XX week as EOM level.  Could estimate based on BOP, EOP
C           XX of last week. 
            IF(RoutingInterval(1) .EQ. 2) SPEndofMonthL = FWRK(8,5) 
            IF(RoutingInterval(2) .EQ. 2) THEN
               MHEndofMonthL = FWRK(8,1) 
               SCEndofMonthL = FWRK(8,2) 
               EREndofMonthL = FWRK(8,3) 
            ENDIF
C           Do preliminary magnitude checks for plans that use BOM levels
            CTMP(1:1) = LLT
            IF(EREndofMonthL .LE. 99 .OR. EREndofMonthL .GE. 999)
     &         WRITE(CTMP,4289) ' Erie BOM level', PeriodEnd(1,1),
     &            MonthAbbrev(PeriodEnd(3,1)), PeriodEnd(4,1),
     &            EREndofMonthL, LLT
            IF(SCEndofMonthL .LE. 99 .OR. SCEndofMonthL .GE. 999)
     &         WRITE(CTMP,4289) ' St. Clair BOM level', 
     &            PeriodEnd(1,1), MonthAbbrev(PeriodEnd(3,1)), 
     &            PeriodEnd(4,1), SCEndofMonthL, LLT
            IF(MHEndofMonthL .LE. 99 .OR. MHEndofMonthL .GE. 999)
     &         WRITE(CTMP,4289) ' Michigan-Huron BOM level', 
     &            PeriodEnd(1,1), MonthAbbrev(PeriodEnd(3,1)), 
     &            PeriodEnd(4,1), MHEndofMonthL, LLT
            IF(SPEndofMonthL .LE. 99 .OR. SPEndofMonthL .GE. 999)
     &         WRITE(CTMP,4289)' Superior BOM level', PeriodEnd(1,1),
     &            MonthAbbrev(PeriodEnd(3,1)), PeriodEnd(4,1),
     &            SPEndofMonthL, LLT
            IF( IDNINT(NominalStMarys) .LE. ZERO) 
     &         WRITE(CTMP,4289) 'monthly mean St. Marys River flow', 
     &            PeriodEnd(1,1), MonthAbbrev(PeriodEnd(3,1)),
     &            PeriodEnd(4,1), SPMonthlyMeanF, LLT
            IF( CTMP(1:1) .NE. LLT) THEN
               CALL LogIt(LOGCDE,1,0)
               CALL ErrMess(1141,'SuperiorRoute'//LLT//'  ',LOGCDE)
               CALL SAVETOD
               CALL ADIOS()
            ENDIF
         ENDIF
      ENDIF
C
C     Beginning of regulation plan calculations
C---------------------------------------------------------------------
C     Determine gates and total sidechannel
C
C
      IF (StMarysSolution .EQ. 1) CONTINUE
C        Don't need to do anything additional here, just route.
C
      IF (StMarysSolution .EQ. 2) THEN
C        Calculate outflow using 1955 Modified Rule of '49
         IF(NewMonth) THEN
C           Syntax is:
C           SUBROUTINE P55MR49(Month (input), Year (input), 
C                Previous Superior Mean Level (input), 
C                Superior Beginning of Month Level (input), 
C                Regulated Flow (output), Gate Setting (output),
C                Side Channel Flow (output))
            CALL P55MR49(PeriodEnd(3,1), PeriodEnd(4,1), PrevSupMean,
     &          SPBOP, NominalStMarys, SupGatesOpen, ComputedNonGate)
         ENDIF
C        XX The following code was snarfed from .EQ. 4 below.  
C        XX Should check to make sure we should even be checking for 
C        XX overtopping...
C        We need to identify the Superior water balance term for 
C        regulation in case the overtopping routine kicks in. 
C        For immediate convenience just use the operational supply
C        for whatever month we happen to be in.
         SPNWB=OPLLO(PeriodEnd(3,1))+ SPOPNBS07(PeriodEnd(3,1),1)
         IF( SupVerbosity .GE. 8 )THEN
            WRITE( CTMP, '(A,A,I5,A,F7.2,A,D18.10,A1)') 
     &            ' Final 55 Modified Rule of 49 results for month of ',
     &            MonthAbbrev(PeriodEnd(3,1)), PeriodEnd(4,1), 
     &            ': desired gates=', SupGatesOpen,
     &            ': desired regulated flow=', NominalStMarys, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF      
      ENDIF


      IF (StMarysSolution .EQ. 3) THEN
C        
         IF(NewMonth) THEN
C
C           Now sort out ways to remain compatible with existing
C           procedures.  
            CALL SetUpSupOperationalSupplies(P77ACompat, PeriodEnd(3,1),
     &         SPNTS, SPNWB, MHNWB, SCNWB, ERNWB, MHRet, SCRet, ERRet)
C
C           XX Modeling/Operational Compatibility Issue:
C           XX    Operational practice sets outflow from the previous
C           XX month to the 'nominal regulated flow', calculated using
C           XX standardized supplies, diversions, and retardations, 
C           XX assuming the outflow is constant over the month, and
C           XX rounded to 10m3s.  On rare occasions (such as when the 
C           XX lake level changes dramatically during a month), the
C           XX regulators use engineering judgement to come up with a
C           XX better estimate of the previous month's outflow.  When
C           XX this situation occurs and how the better flow is 
C           XX determined are not defined adequately to be modeled.
C           XX So we'll simply take last month's 'nominal regulated'
C           XX (i.e. 'plan') flow, and round it to 10m3s.
C    
C
            CALL Plan1977A( PeriodEnd(3,1), PeriodEnd(4,1), 
     &         Round(NominalStMarys, 1), PreviousGates, SPNWB, 
     &         SPEndofMonthL, MHEndofMonthL, SCEndofMonthL, 
     &         EREndofMonthL, NominalStMarys, SupGatesOpen)
C 
            CALL SupCriterions( PeriodEnd(4,1), PeriodEnd(3,1), 
     &         SPNWB, MHNWB, SCNWB, ERNWB, MHRet, SCRet, ERRet,
     &         SPEndofMonthL, MHEndofMonthL, SCEndofMonthL, 
     &         EREndofMonthL, LStoSWP, DBLE(SecsInMonth),NominalStMarys, 
     &         SWPier, ComputedNonGate, SupGatesOpen, SPRegCodes)
C              
C           At this point, the number of gates and non-gated flow
C           for the month are established.
            IF( SupVerbosity .GE. 8 )THEN
               WRITE( CTMP, '(A,A,I5,A,F7.2,A,D18.10,A1)') 
     &               ' Final Plan1977A results for month of ',
     &               MonthAbbrev(PeriodEnd(3,1)), PeriodEnd(4,1), 
     &               ': desired gates=', SupGatesOpen,
     &               ': desired regulated flow=', NominalStMarys, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF            
            SWPier = ZERO
            GateFlow = ZERO
            LastComputedNonGate = ComputedNonGate 
            IF( SupVerbosity .GE. 8 )THEN
               WRITE( CTMP, '(A,D18.10,A1)') 
     &               ' LastComputedNonGate= ', LastComputedNonGate, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF            
C
         ELSE
C           Running under Plan 1977A, but interval does not include a new month.  Keep previous settings
            IF( SupVerbosity .GE. 5 )THEN
               WRITE( CTMP, '(A,A,D18.10,A,F9.2,A1)') 
     &               ' Still in same month, so using previous ',
     &               'regulated flow=', LastComputedNonGate,
     &               ' and previous gates=', PreviousGates, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
            SupGatesOpen = PreviousGates 
            ComputedNonGate = LastComputedNonGate 
         ENDIF
      ENDIF
C     End of P77A block
C     IF(SupVerbosity .GE. 7 ) THEN
C        WRITE(CTMP,'(2X,3(A,F15.8),A,F15.3,A,F9.3,A,F9.3,A1)') 
C    &      ' B4FG SPBOP=', SPBOP, ' Superior NTS=', SPNTS,
C    &      ' Ice=', SMICE, ' dT=', dT(1), ' SupGatesOpen=', 
C    &      SupGatesOpen, ' ComputedNonGate=',ComputedNonGate, LLT
C        CALL LogIt(LOGD,0,0)
C        WRITE(CTMP,'(2X,5(A,F15.8),A1)')' B4FG SideChannel=',
C    &      SupSideChannel,' Forced Flow=',SupForceFlow,' ForcedGates=', 
C    &      NumSupForceGates, ' Deviation Flow=', SupDevFlow,
C    &      ' SupForceEOPLev=', SupForceEOPLev, LLT
C        CALL LogIt(LOGD,0,0)
C        IF(SupVerbosity .GE. 9 ) THEN
C           WRITE(CTMP,'(2X,3(A,F15.8),A1)') 
C    &         ' B4FG SPEOP=', SPEOP, ' SPMLV=', SPMLV, ' StMRivFlow=', 
C    &         StMRivFlow, LLT
C           CALL LogIt(LOGD,0,0)
C        ENDIF
C     ENDIF
C
C
C    -------------------------------------------------------------
C    Plan 2012
C    --------------------------------------------------------------  
      IF (StMarysSolution .EQ. 4) THEN
C        Plan 2012 
C        This plan was developed as Natural Plan 3 during the IUGLS.
C        
         CBUSSLevel = MissingVal
         IF(NewMonth) THEN
            CALL SetUpSupOperationalSupplies_P2012(P77ACompat,
     &         PeriodEnd(3,1),
     &         SPNWB, MHNWB, SCNWB, ERNWB, MHRet, SCRet, ERRet)
C
            CALL P2012_Regulation( PeriodEnd(3,1), SupVerbosity, 
     &         SPEndofMonthL, MHEndofMonthL, P77AWinterMax, 
     &         FlowRound, NominalStMarys, SupGatesOpen, ComputedNonGate,
     &         PreviousGates, IndexGatesSide)
 
            CALL P2012_Criterions( PeriodEnd(4,1), PeriodEnd(3,1), 
     &         SPNWB, MHNWB, SCNWB, ERNWB, MHRet, SCRet, ERRet,
     &         SPEndofMonthL, MHEndofMonthL, SCEndofMonthL, 
     &         EREndofMonthL, LStoSWP, DBLE(SecsInMonth),
     &         NominalStMarys, SWPier, ComputedNonGate, 
     &         SupGatesOpen, SPRegCodes, CBUSSLevel)          
     
C           At this point, the number of gates and non-gated flow
C           for the month are established.
            IF( SupVerbosity .GE. 8 )THEN
               WRITE( CTMP, '(A,A,I5,A,F7.2,A,D18.10,A1)') 
     &               ' Final Plan 2012 results for month of ',
     &               MonthAbbrev(PeriodEnd(3,1)), PeriodEnd(4,1), 
     &               ': desired gates=', SupGatesOpen,
     &               ': desired regulated flow=', NominalStMarys, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF            
            SWPier = ZERO
            GateFlow = ZERO
            LastComputedNonGate = ComputedNonGate 
            IF( SupVerbosity .GE. 8 )THEN
               WRITE( CTMP, '(A,D18.10,A1)') 
     &               ' LastComputedNonGate= ', LastComputedNonGate, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF            
C 
         ELSE
C           Running under Plan 2012, but interval does not include a new month.  Keep previous settings
            IF( SupVerbosity .GE. 5 )THEN
               WRITE( CTMP, '(A,A,D18.10,A,F9.2,A1)') 
     &               ' Still in same month, so using previous ',
     &               'regulated flow=', LastComputedNonGate,
     &               ' and previous gates=', PreviousGates, LLT
               CALL LogIt(LOGD,0,0)
            END IF
            
            SupGatesOpen = PreviousGates 
            ComputedNonGate = LastComputedNonGate 
         END IF
      
      END IF            
C
C
C------------------------------------------------------------------------
C     END OF REG PLANS  - should have plan outflows or gates, etc set now 
C
      IF( IDNINT(SupSideChannel) .NE. MissingVal) THEN
C     XX Side Channel forcing disabled because this requires more thought 
C     XX to implement properly nd we don't feel a compelling need for 
C     XX the capability right now
         STOP 'NoSideF'
      END IF
C
C     Check to see if outflows from a regulation plan jeopardize 
C     the structure. All plan flows are known and we now have to verify
C     these with physical capacities.  (See notes from 8 Oct 2008)
C------------------------------------------------------------------------
C Bookmark:PChecks
C
      IF(StMarysSolution .GE. 2 ) THEN
C
C        See if Lake Superior level causes overtopping Comp Works, or
C        violates sidechannel minimums or maximums.
C        Just use BOP level for simplicity.
         CWFreeBoard = CompWorksFreeboard(
     &      PeriodEnd(3,1), SPBOP, CWOvertopElev)
         IF( SupVerbosity .GE. 3 )THEN
           WRITE( CTMP, '(A,F10.6,A,F9.6,A,F9.4,A1)') 
     &         ' Superior level of ', SPBOP, 'm leaves ', CWFreeBoard,
     &         'm of freeboard below overtop at ', CWOvertopElev, LLT
           CALL LogIt(LOGD,0,0)
         ENDIF
C
C        Approximate sidechan max capacity for use in AdjustDueToCompOvertop 
C        based on SPBOP and desired flow
         CALL GetSideChanMaxCapacity( SPBOP, NominalStMarys, 
     &       PeriodEnd(3,1), SideChanMaxFlow)
        
C
         AdjOTFlow = MissingVal
         IF (CWFreeBoard .LT. CWMinimumFreeboard(PeriodEnd(3,1))) THEN
            WRITE( CTMP, '(A,I2,A3,I5,A,F10.6,A,F10.6,A,A1)') 
     &         ' Warning - on ', PeriodEnd(1,1),  
     &         MonthAbbrev(PeriodEnd(3,1)), PeriodEnd(4,1), 
     &         ' estimate only ',CWFreeBoard, 
     &         'm of freeboard but need ', 
     &         CWMinimumFreeboard(PeriodEnd(3,1)),
     &         'm, so overriding regulated flow', LLT
            CALL LogIt(LOGD,0,0)
C           should handle lots of these messages to console better by 
C           simply counting and issuing one at end.
            IF( SupVerbosity .GE. 7 ) CALL LogIt(LOGCE,0,0)
C
            CALL AdjustDueToCompOvertop(PeriodEnd(1,1), PeriodEnd(3,1), 
     &         PeriodEnd(4,1), dT(1), 3, SPBOP, DBLE(SPNWB), 
     &         SideChanMaxFlow, SupGatesOpen, 
     &         ComputedNonGate, AdjOTFlow, SPRegCodes)

C           Now update the final value of side channel flow
            SupSideChannel = ComputedNonGate
C  
         ELSE
            IF( SupVerbosity .GE. 7 )THEN
               WRITE( CTMP, '(4(A,F10.6),A1)') 
     &           ' Superior level of ', SPBOP, 'm leaves ', CWFreeBoard,
     &           'm of freeboard, which meets requirement of ', 
     &           CWMinimumFreeboard(PeriodEnd(3,1)), 
     &           'm before overtop at ', CWOvertopElev, LLT
               CALL LogIt(LOGD,0,1)
            ENDIF
            SupSideChannel = ComputedNonGate
         ENDIF
C
C
C        Verify that Sidechannel flows are within allowable range
         CALL CheckSideChannelLimits(SideChanMaxFlow, SupSideChannel, 
     &      PeriodEnd(3,1), PeriodEnd(4,1), PeriodEnd(1,1), SPRegCodes)
C 
         WRITE(SupSideChanFlow,'(I6)')IDNINT(SupSideChannel)
         IF( SupVerbosity .GE. 8 )THEN
            WRITE( CTMP, '(A,A,A,D18.10,A,A1)') 
     &         ' Setting SupSideChanFlow= ', SupSideChanFlow, ' (',
     &         SupSideChannel, ')',LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
      ENDIF
C
C     Check for forced gates condition specified in time-series data
C------------------------------------------------------------------------
C Bookmark:ForceG
C     Forced Gates 
      IF( IDNINT(NumSupForceGates) .NE. MissingVal) THEN
C
C     XX Disabled because this requires more thought to implement properly 
C     XX and we don't feel a compelling need for the capability right now
      STOP 'NoForcin'
C        Dump at run-time if computation not supported yet
C         IF( IDNINT(SupSideChannel) .NE. MissingVal .OR.
         IF( IDNINT(SupForceFlow) .NE. MissingVal ) THEN 
            WRITE(CTMP,'(A,I2,A3,I5,A,A1)') 
     &         ' Unsupported computation for ', PeriodEnd(1,1), 
     &         MonthAbbrev(PeriodEnd(3,1)), PeriodEnd(4,1), 
     &         ': Forcing gates and total flow.', LLT
C     &         ': Forcing gates and either total flow or sidechannel.',
C     &         LLT
            CALL LogIt(LOGCDE,1,0)
            WRITE(CTMP,'(A,F9.3,A3,I5,A,A1)') 
     &         ' Forced Sidechannel: ', SupSideChannel,
     &         MonthAbbrev(PeriodEnd(3,1)), PeriodEnd(4,1), 
     &         ': Forcing gates and either total flow or sidechannel.',
     &         LLT
            CALL LogIt(LOGCDE,1,0)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
C        Issue a warning or log entry if overriding a real value
         IF( NINT(SupGatesOpen) .NE. MissingVal ) THEN
            WRITE( CTMP, '(A,I2,A3,I5,A,F5.2,A,F5.2,A)') 
     &         ' Warning for ', PeriodEnd(1,1), 
     &         MonthAbbrev(PeriodEnd(3,1)), PeriodEnd(4,1), 
     &         ': forcing ', NumSupForceGates, 
     &         ' over reg plan setting of ', SupGatesOpen, LLT
            IF( SupVerbosity .GE. 3 )THEN
               CALL LogIt(LOGCD,1,1)
            ELSE
               CALL LogIt(LOGD,1,1)
            ENDIF
         ENDIF
         SPRegCodes = '   Forced Gates                         ' 
         SupGatesOpen = NumSupForceGates
      ENDIF
C
C
C
C     Check for forced flow condition specified in time-series data
C------------------------------------------------------------------------
C Bookmark:ForceQ
C
      IF(IDNINT(SupForceFlow) .NE. MissingVal)THEN
         STMRIVFLOW = SupForceFlow
C        XX Should we make a stab here regarding how many gates
C        XX are required to achieve the forced flow?  Not Yet.
C
C        Dump at run-time if computation not supported yet
C         IF( IDNINT(SupSideChannel) .NE. MissingVal .OR.
C     &      IDNINT(NumSupForceGates) .NE. MissingVal ) THEN 
          IF( IDNINT(NumSupForceGates) .NE. MissingVal ) THEN 
            WRITE(CTMP,'(A,I2,A3,I5,A,A1)') 
     &         ' Unsupported computation for ', PeriodEnd(1,1), 
     &         MonthAbbrev(PeriodEnd(3,1)), PeriodEnd(4,1), 
     &         ': Forcing total flow and either gates or sidechannel.',
     &         LLT
            CALL LogIt(LOGCDE,1,0)
            WRITE(CTMP,'(2(A,F9.2),A1)') ' NumSupForceGates=', 
     &         NumSupForceGates, '    SupSideChannel',SupSideChannel,LLT
            CALL LogIt(LOGCDE,1,0)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
         SupGatesOpen = -1.0
         ComputedNonGate = MissingVal
         SupSideChannel = MissingVal
         WRITE(SupSideChanFlow,'(I6)')IDNINT(SupSideChannel)
         SPRegCodes = '   Forced Flow                          ' 
         WRITE(SupGateFlow, '(I6)') MissingVal
         IF(SupVerbosity .GE. 3 ) THEN
            IF (StMarysSolution .GT. 1) THEN
               WRITE( CTMP, 2255) LLT
               CALL LogIt(LOGD,1,1)
            ENDIF
            WRITE( CTMP, 2256) SupForceFlow, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
      ENDIF
 2255 FORMAT( 'Now forget whatever Superior outflow may have been ',
     &   'computed using a regulation plan!', A1)
 2256 FORMAT( 'Forcing St Marys River flow to specified value:', 
     &   F12.3, ' m3s', A1) 
C
C     Beginning of ACTUAL routing through Superior
C---------------------------------------------------------------------
C Bookmark:RoutingSetup
C     See how flows are calculated.  
C
C     If outflow is allowed to vary within a routing period as a
C     function of lake level, then change-in-storage isn't linear,
C     and need to integrate over a number of increments while solving
C     simultaneous equations.
C
C     If assumed constant over the period (such as Plan77A/linear
C     certain gate setting will result in assumes a constant outflow 
C     over the month), then storage change is linear, and StmRivFlow 
C     equals NominalStMarys (or SupForceFlow).  
C
C     Forced flows are assumed to already be period-avergae, so must
C     only perform linear change in storage calcualtion.
C
C     Traditional Superior regulation practice was to make the
C     simplifying assumption that the regulated outflow remains
C     constant during the month.  With outflow not a function
C     of head, the change-in-storage over the month is a straight
C     line, so knowing the outflow we can plug in the storage 
C     equation and solve for EOP level, then mean level.  This 
C     case (i.e., 1 increment) follows through with the general
C     logic in this loop, but occasionally fails to reproduce the 
C     exact sequence of rounds necessary to reproduce the results of 
C     previous models.  So if the relevant parameters imply that the
C     modeler is trying to replicate the errors of the older models,
C     the computation is handled separately as follows:
C   
      IF ( BackCompat .OR. IDNINT(SupForceFlow) .NE. MissingVal) THEN
C        Assume that need full backward compatibility with earlier 
C        models.  Linear change-in-storage, just compute EOP and 
C        mean levels, no increment loop needed.   Also applies to
C        forced flows, since constant outflow.
C          
         IF( SupVerbosity .GE. 4) THEN
            WRITE( CTMP, '(A,A1)') ' Performing'//
     &         ' constant outflow calculation for Lake Superior', LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
C        Forget gate settings, monthly average outflow defined as
C        the plan flow, so StmRivFlow = NominalStMarys, unless
C        forced flow or gates.
         IF(IDNINT(SupForceFlow)     .EQ. MissingVal .AND. 
     &      IDNINT(NumSupForceGates) .EQ. MissingVal) 
     &      STMRIVFLOW = NominalStMarys
C
C        Sabotage significant digits in keeping with existing practice
         IF(SupVerbosity .GE. 7 ) THEN
            TmpA = GetArea(1,SPMLV)
            WRITE( CTMP, 6234) Round(SPBOP, BOMRound), SPNTS,
     &         STMRivFlow, RoutingIncrement(1), dT(1), TmpA, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
 6234    FORMAT( '   Linear change in storage: BOP=', F9.5,
     &      ' NTS=',F8.2, ' St Marys Flow=', F7.1, ' num inc=', I3, 
     &      ' dT=', F10.1, ' Area=', D27.13, A1)
C
         SPEOP=Round(SPBOP, BOMRound) + Round((SPNTS-STMRIVFLOW)* 
     &      RoutingIncrement(1) *dT(1)/GETAREA(1,SPBOP), BOMRound)
         SPMLV = Round((SPEOP + SPBOP) * HALF, BOMRound)

         IF ( IDNINT(SupForceFlow) .EQ. MissingVal) THEN
            CALL GetSWPier( SPMLV, PeriodEnd(3,1), SupGatesOpen, 
     &         IDNINT(SupSideChannel), BOMRound, LStoSWP, 
     &         8,PeriodEnd(4,1),'Superior Routing', SWPier, GateFlow)
C
            IF( SupVerbosity .GE. 5) THEN
               WRITE( CTMP, '(A,I5,A,A1)') ' Re-computed gate flow '//
     &            'with levels rounded to ', FlowRound, ' places:', LLT
               CALL LogIt(LOGD,0,0)
               WRITE( CTMP, 2305) GateFlow, SPMLV, SWPier, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
         ENDIF
 2305    FORMAT( '   Gated flow computed as ', F9.4, 
     &      ', based on Sup/SW Pier levels of ', 2F9.4, A1)
C
C        Update report fields
         WRITE(SWPierLevel, '(F8.3)') Round(SWPier, BOMRound)
         WRITE(SupGateFlow, '(I6)') IDNINT(Round(GateFlow, 0))
C
C        Do a continuity check to see if creating/destroying water
         SPSolutionCheck = ROUND(DABS( STMRIVFLOW - SPNTS +   
     &      ((SPEOP - SPBOP) * GETAREA(1,SPMLV))/dT(1)), FlowRound)
         IF(SPSolutionCheck .GT. SPErrorFlow)
     &      SPErrorFlow=SPSolutionCheck 
         IF(SupVerbosity .GT. 7 ) THEN
            WRITE(CTMP, 8033) PeriodEnd(5,1), 1,'SP', 
     &         SPSolutionCheck, SolutionTest, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C        Value of SolutionTest assigned in CONTROL.FOR 
         IF(SPSolutionCheck .GT. SolutionTest ) THEN
               SPSolveErrors = SPSolveErrors + 1
            IF(SupVerbosity .GT. 1) THEN
               WRITE(CTMP, 8072) PeriodEnd(5,1), 1, PeriodEnd(1,1), 
     &            PeriodEnd(3,1), PeriodEnd(4,1), LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
         ELSE
            IF(SupVerbosity .GT. 7 ) THEN
               WRITE(CTMP, 8071) 1, PeriodEnd(5,1), LLT
               CALL LogIt(LOGD,0,1)
            ENDIF
         ENDIF
C        End of constant Superior outflow condition
C        Short-circuit the looping through increments and jump in to
C        spot after SPEOP and SPMLV get computed for non-linear
C        change in storage.
         GOTO 1233
      ENDIF
C
C     At this point, the outflow varies with lake level.  So have 
C     non-linear change-in-storage, and need to integrate outflow 
C     over RoutingIncrement(1) sub-divisions of a routing interval.
C     Within each increment solve simultaneous equations for 
C     outflow, mean increment level, and end-of-increment level.
C
      IF( SupVerbosity .GE. 5) THEN
         WRITE(CTMP, '(A,I3,A,A1)') ' Lake Superior outflow varies '//
     &      'during an interval; integrating over ',RoutingIncrement(1), 
     &      ' increments.', LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
      SWPier = ZERO
      GateFlow = ZERO
C
C     Now loop through each increment of the interval to integrate
C     the change in storage
C------------------------------------------------------------------------
C Bookmark:TopofRoutingLoop
C
      DO 1222 I=1, RoutingIncrement(1)            
C
C        Will do iterative solution of 3 simultaneous eqns:
C           mean level for increment = (SPBOI + SPEOI)/2
C           outflow= f(mean level) or SPOutI = f(SPBOI + SPEOI)/2
C           (SPEOI-SPBOI)*Area/dT = Net Total Supply - outflow
C          
C        Initialize BOI levels as EOI levels from prev increment
         SPBOI = SPEOI
C        Initialize the outflows with bogus number
         SPOutI = MissingVal
C        Initialize the mean level value with the initial level
         SPIMLV = SPBOI
C
         EOILevels2Ago = -999
C
C        Dump some debugging info if desired:
         IF(SupVerbosity .GE. 8 ) THEN
            TmpA = GetArea(1,SPIMLV)
            WRITE(CTMP,'(A,I3,3(A,F15.8),A,D25.13,A1)') 
     &         ' Increment ', I, ' SPBOI=', SPBOI, ' SPEOI=', SPEOI,
     &         ' SPIMLV=', SPIMLV, ' Area=', TmpA, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
         IF( StMarysSolution .EQ. 1 )THEN 
C           Compute flow using outflow equation, such as 1887 relationship
C           Loop through so many iterations to solve simultaneous equations
            DO 1310, K = 1, MaxIterations
C              Save previous outflow for comparison to tolerance
               LastFlow = SPOutI
               LastSPEOI = SPEOI
C       
C              Approximate the outflow for the iteration
               SPOutI = Round( OUTFLOW( SPIMLV, ZERO, SPK, SPYM, SPA, 
     &                 SPB, SPWT, SPC)- SMICE, FlowRound)
C                 
C              Compute EOP levels for the iteration
               SPEOI = Round( SPBOI + (SPNTS-SPOutI) * dT(1)/
     &            GETAREA(1,SPIMLV), LevelRound )
C
C              Dump some debugging info if desired:
               IF(SupVerbosity .GE. 9 ) THEN
                  TmpA = GetArea(1,SPIMLV)
                  WRITE(CTMP,'(2(A,I3),3(A,F15.8),A1)') ' Increment ',
     &               I, ' Iteration ', K, ' Approx. mean=', SPIMLV, 
     &               ' Approx. outflow=', SPOutI, ' delta S=', 
     &               (SPNTS-SPOutI) * dT(1)/TmpA, LLT
                  CALL LogIt(LOGD,0,0)
               ENDIF
C
C              Compute mean level for the iteration
               SPIMLV = (SPEOI + SPBOI) * HALF
C
C              Dump some debugging info if desired:
               IF(SupVerbosity .GE. 9 ) THEN
C                 Prevent redundant logging within OUTFLOW()
                  SupVerbosity = 1
C                 Using FWRK(10,4) as scratch.
                  FWRK(10,4) = Round( OUTFLOW( SPIMLV, ZERO, SPK, SPYM, 
     &               SPA, SPB, SPWT, SPC) -SMICE, FlowRound)
                  WRITE(CTMP,'(2(A,I3),4(A,F15.8),A1)') ' Increment ',  
     &               I, ' Iteration ', K, ' SPBOI=', SPBOI, ' SPEOI=', 
     &               SPEOI,' SPIMLV=',SPIMLV,' Outflow=', FWRK(10,4),LLT
                  CALL LogIt(LOGD,0,0)
                  SupVerbosity = 9
               ENDIF
C
               IF( DABS(SPOutI - LastFlow) .LE. Tolerance) GOTO 1350
C
               IF (LevelRound .NE. -99) THEN
C                 Check to see if oscillating between two consecutive 
C                 values of rounded EOI levels.
C    
C 101             FORMAT(' Inc ', I3,' Iteration',I4,' SP 2 iter ago=', 
C    &               I10,' LastSPEOI= ', F11.5, ' This SPEOI=', I10,A1)
C                 WRITE(CTMP,101) I, K,  EOILevels2Ago, LastSPEOI, 
C    &               IDNINT(SPEOI*LevelToInteger), LLT
C                 CALL LogIt(LOGD,0,0)
                  IF(EOILevels2Ago .EQ.IDNINT(SPEOI*LevelToInteger))THEN
                     IF(SupVerbosity .GE. 7 ) THEN
                        WRITE(CTMP,'(4X,A,I4,A,A1)') 
     &                     '. . . solution converges at ', K,
     &                     ' iterations, but desired tolerance not met!'
     &                     ,LLT
                        CALL LogIt(LOGD,0,0)
                     ENDIF
                     GOTO 1350
                  ENDIF            
C            
C                 Convert to integers, with ones digit at rounding point, 
C                 and update EOILevels2Ago.
                  EOILevels2Ago = IDNINT(LastSPEOI*LevelToInteger)
               ENDIF            
C
 1310          CONTINUE
C
C           If fall to here, must have exceeded MaxIterations, so
C           failed to converge.
            WRITE(CTMP,614)MaxIterations, I, LLT
            CALL LogIt(LOGCD,0,0)
C
C           Re-compute mean outflow
 1350       SPOutI = Round(OUTFLOW( SPIMLV, ZERO, SPK, 
     &         SPYM, SPA, SPB, SPWT, SPC) - SMICE, FlowRound)
C
C           End of increment solution using an outflow equation
C
         ELSE IF(StMarysSolution .GE. 2)THEN 
C           This case routes outflows due to regulation plans (i.e.,
C           flow through gates and sidechannel)
C
C           The number of gates and the non-gated flow are established.
C           Not assuming that outflow during a period is constant.
C           Must compute actual gated flow within a period, and add to 
C           non-gated flow to get St Marys flow.
C
            IF(SupVerbosity .GE. 8 ) THEN
               WRITE(CTMP,'(4X,A,A,F5.2,A,I12,A,F12.4,A,I2,A1)') 
     &            '. . . Superior routing, calling SolveGatesSide(',
     &            'SupGatesOpen=', SupGatesOpen,
     &            ', IDNINT(SupSideChannel)=', IDNINT(SupSideChannel),
     &            ', SPNTS=', SPNTS,
     &            ', LStoSWP=', LStoSWP, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
            CALL SolveGatesSide(SupGatesOpen, IDNINT(SupSideChannel), 
     &         PeriodEnd(3,1), SPBOI, SPNTS, 
     &         dT(1), LStoSWP, SPEOI, SPIMLV, SWPI, GateFlowI, SPOutI)
C           Accumulate SWPier means for increment, for later averaging
            SWPier = SWPier + SWPI
            GateFlow = GateFlow + GateFlowI
            IF(SupVerbosity .GE. 8 ) THEN
               WRITE(CTMP,'(4X,A,A,F5.2,A,F12.4,A,F12.4,A,I2,A1)') 
     &            '. . . Superior routing, after SolveGatesSide(',
     &            'SupGatesOpen=', SupGatesOpen,
     &            ', spouti=', spouti,
     &            ', gateflow=', gateflow,
     &            ', LStoSWP=', LStoSWP, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
C
         ELSE
            STOP 'Chief Wiggins'
         ENDIF
C
C
         IF( SupVerbosity .GE. 8) THEN
            WRITE(CTMP,605)I, SPBOI, SPEOI, LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,606)I, SPOutI, SPIMLV, LLT
            CALL LogIt(LOGD,0,0)
            IF( StMarysSolution .EQ. 1) THEN
               WRITE(CTMP,612)I, K, Tolerance, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
         ENDIF
         IF( StMarysSolution .EQ. 1 .AND. K .GE. MaxIterations )THEN
            CALL ErrMess(1091,'SuperiorRoute'//LLT//'  ',LOGCDE)
            CALL BuildInfo
            CALL LogIt(LOGCDE,0,0)
            CALL SAVETOD
            CALL ADIOS()
         ENDIF
C
C        Accumulate outflow for increment, for later averaging
         STMRIVFLOW = STMRIVFLOW + SPOutI
C
C        Now check that the simultaneous equations were solved
C        correctly.  Verify that the inflow plus NBS, minus the
C        average rate of change in storage during the period,
C        nearly equals the outflow for the period.  Store the 
C        difference in the SolutionCheck variables for debugging
C        info.  Note that this implemenation uses surface area as
C        BOI value.  Can easily fix this by introducing a new 
C        increment mean level if want to check variable surface 
C        area runs.
C----&------------------------------------------------------------------
         SPSolutionCheck = ROUND(DABS( SPOutI -  SPNTS +   
     &      ((SPEOI - SPBOI) * GETAREA(1,SPEOI))/dT(1)), FlowRound)
         IF(SPSolutionCheck .GT. SPErrorFlow)SPErrorFlow=SPSolutionCheck
         IF(SupVerbosity .GT. 7 ) THEN
            WRITE(CTMP, 8033) PeriodEnd(5,1), I,'SP', 
     &         SPSolutionCheck, SolutionTest, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C        Value of SolutionTest assigned in CONTROL.FOR 
         IF(SPSolutionCheck .GT. SolutionTest ) THEN
            SPSolveErrors = SPSolveErrors + 1
            WRITE(CTMP, 8072) PeriodEnd(5,1), I, PeriodEnd(1,1), 
     &         PeriodEnd(3,1), PeriodEnd(4,1), LLT
            CALL LogIt(LOGD,0,0)
         ELSE
            IF(SupVerbosity .GT. 8 ) THEN
               WRITE(CTMP, 8071) I, PeriodEnd(5,1), LLT
               CALL LogIt(LOGD,0,1)
            ENDIF
         ENDIF
C
C        Accumulate the EOI levels into the mean level variables
 1222    SPMLV = SPMLV + SPEOI
C            
C     Ok, have looped through each increment of a period.  The EOP 
C     values are the final EOI levels.
      SPEOP = SPEOI
C           
C     Must adjust the mean level values to subtract half of the final
C     EOI values (was only supposed to add a half inside the
C     loop, but added a whole).  Then divide by number of increments.
      SPMLV = (SPMLV - (SPEOP * HALF)) / DBLE(RoutingIncrement(1))
C
C     Now calculate mean outflow for the period.
      STMRIVFLOW = Round(STMRIVFLOW / DBLE(RoutingIncrement(1)), 
     &   FlowRound)
C
      IF(StMarysSolution .GE. 2)THEN 
         SWPier = Round(SWPier / DBLE(RoutingIncrement(1)), LevelRound)
C        Since SupSideChannel is constant, these are mathematically equivalent
C        GateFlow = StMRivFlow - SupSideChannel
         GateFlow= Round(GateFlow / DBLE(RoutingIncrement(1)),FlowRound)
C        Update report fields
         WRITE(SWPierLevel, '(F8.3)') SWPier
         WRITE(SupGateFlow, '(I6)') IDNINT(GateFlow)
      ENDIF

C     Routing integration has finished.
C     Calculations for the entire routing period finished at this point.
C
C------------------------------------------------------------------------
C     Each St. Marys Solution method could have different intermediate 
C     results in outputs. 
C     (i.e, Plan 1977A vs outflow eqn vs Mod 1955 vs some fencepost have
C     such as pre-project flow)
C Bookmark:WriteOutput
C
C     The following line is a "stupid FORTRAN trick" to convert a 
C     character string to a double precision number.
      IF ( ( SupPlanFlow) .NE. '      ') then
          READ (SupPlanFlow, *) SupPlanFlow_Dble
      END IF
C     Send the results to the superior output file
 1233 IF (StMarysSolution .EQ. 1) THEN      
C        Equation results
C        Fit std format by compute SW Pier level, putting 0's in for 
C        gates and gateflow, and show flow as all sidechannel
         SupGatesOpen = 0.0
         CALL GetSWPier(SPMLV, PeriodEnd(3,1), SupGatesOpen,  
     &      IDNINT(Round(STMRIVFLOW,0)), -2, 3, 
     &      9,PeriodEnd(4,1),'Superior Output1',SWPier, GateFlow)

         WRITE(CTMP, 8080) PeriodEnd(1,1), PeriodEnd(3,1), 
     &      PeriodEnd(4,1), EOPType(1), SPBOP, SPMLV, SPEOP, SWPier, 
     &      SupGatesOpen, STMRIVFLOW, 0.0, STMRIVFLOW, LLT  
         CALL LogIt(LOGP,0,0)
C
      ELSE IF (StMarysSolution .EQ. 2) THEN 
C        Modified 1955 Rule of '49 plan results
       
         WRITE(CTMP, 8081) PeriodEnd(1,1), PeriodEnd(3,1),
     &      PeriodEnd(4,1), EOPType(1), SupBOMLevel, Round(SPMLV,-2), 
     &      Round(SPEOP,-2), SWPierLevel, SupGatesOpen,
     &      '        ', '        ', '        ', 
     &      IDNINT(Round(STMRIVFLOW,0)),  
     &      SupGateFlow, SupSideChanFlow, '        ', AdjOTFlow,
     &      SPRegCodes, LLT  
         CALL LogIt(LOGP,0,0)
C         
      ELSE IF (StMarysSolution .EQ. 3) THEN
C        Plan 1977A results
C
C        If not a monthly routing, some values have meaning only at BOM.
         IF(.NOT. NewMonth) THEN
            BalancingFlow = '        '
            SupPlanFlow = '      '
            SupPreProject = '      '
            SupBOMLevel = '        '
            IF(IDNINT(SupForceFlow)     .EQ. MissingVal .AND. 
     &         IDNINT(NumSupForceGates) .EQ. MissingVal) 
     &         SPRegCodes = '                                        '
         ENDIF

         IF(SupVerbosity .GE. 5 ) THEN
            WRITE(CTMP, '(A)') P77AHeader//LLT
            CALL LogIt(LOGD,1,0)
         ENDIF
         WRITE(CTMP, 8081) PeriodEnd(1,1), PeriodEnd(3,1),
     &      PeriodEnd(4,1), EOPType(1), SupBOMLevel, Round(SPMLV,-2),
     &      Round(SPEOP,-2), SWPierLevel, SupGatesOpen, BalancingFlow, 
     &      IDINT(Round(SupPlanFlow_Dble, 1)), SupPreProject, 
     &      IDNINT(Round(STMRIVFLOW,0)), SupGateFlow, SupSideChanFlow, 
     &      IDNINT(Round(NominalStMarys/10,0)*10), AdjOTFlow,
     &      SPRegCodes, LLT  

c         WRITE(CTMP, 8081) PeriodEnd(1,1), PeriodEnd(3,1),
c     &      PeriodEnd(4,1), EOPType(1), SupBOMLevel, Round(SPMLV,-2),
c     &      Round(SPEOP,-2), SWPierLevel, SupGatesOpen, BalancingFlow, 
c     &      SupPlanFlow, SupPreProject, 
c     &      IDNINT(Round(STMRIVFLOW,0)), SupGateFlow, SupSideChanFlow, 
c     &      IDNINT(Round(NominalStMarys/10,0)*10), AdjOTFlow,
c     &      SPRegCodes, LLT  

         CALL LogIt(LOGP,0,0)
         IF(SupVerbosity .GE. 3 ) CALL LogIt(LOGD,0,0)
C
      ELSE IF (StMarysSolution .EQ. 4) THEN 
C        Plan 2012 results
         IF(SupVerbosity .GE. 5 ) THEN
            WRITE(CTMP, '(A)') P2012Header//LLT
            CALL LogIt(LOGD,1,0)
         ENDIF
         WRITE(CTMP, 8082) PeriodEnd(1,1), PeriodEnd(3,1),
     &      PeriodEnd(4,1), EOPType(1), SupBOMLevel, Round(SPMLV,-2),
     &      Round(SPEOP,-2), SWPierLevel, SupGatesOpen, BalancingFlow, 
     &      IDINT(Round(SupPlanFlow_Dble, 1)), SupPreProject, 
     &      IDNINT(Round(STMRIVFLOW,0)), SupGateFlow, SupSideChanFlow,
     &      IDNINT(Round(NominalStMarys/10,0)*10), AdjOTFlow,
     &      CBUSSLevel, SPRegCodes, LLT  
c         WRITE(CTMP, 8082) PeriodEnd(1,1), PeriodEnd(3,1),
c     &      PeriodEnd(4,1), EOPType(1), SupBOMLevel, Round(SPMLV,-2),
c     &      Round(SPEOP,-2), SWPierLevel, SupGatesOpen, BalancingFlow, 
c     &      (SupPlanFlow), SupPreProject, 
c     &      IDNINT(Round(STMRIVFLOW,0)), SupGateFlow, SupSideChanFlow,
c     &      IDNINT(Round(NominalStMarys/10,0)*10), AdjOTFlow,
c     &      CBUSSLevel, SPRegCodes, LLT       
         CALL LogIt(LOGP,0,0)
         IF(SupVerbosity .GE. 3 ) CALL LogIt(LOGD,0,0)
C
      ELSE
         STOP 'Otto'
      ENDIF
      AdjOTFlow = MissingVal 
      IF(SupVerbosity .GE. 3 ) THEN
         WRITE(CTMP, 8089) PeriodEnd(5,1), Day(1), PeriodEnd(1,1), 
     &      PeriodEnd(3,1), PeriodEnd(4,1), SPBOP, SPMLV, 
     &      SPEOP, STMRIVFLOW, LLT  
         CALL LogIt(LOGD,0,0)
         CALL WriteHighPrecision(1, SPMLV, SPEOP, StMRivFlow)
      ENDIF
C
C     Force the EOP level, if necessary, for use as next period's BOP
C     level:
      IF( SupForceEOPLev .NE. MissingVal )THEN
         IF(SupVerbosity .GE. 3 ) THEN
            WRITE(CTMP,'(2X,A,F8.4,A,F12.7,A,A1)') 
     &         ' Replacing modeled Superior EOP level (', SPEOP, 
     &         ') with forced EOP level:', SupForceEOPLev, 'm.', LLT
            CALL LogIt(LOGD,1,1)
         ENDIF
         SPEOP = SupForceEOPLev
      ENDIF
C
 8080 FORMAT( 1X, 2(I3), I5, 1X, A1, 4F8.3, F6.2, 22X, 3F7.1, A1)

 8081 FORMAT( 1X, 2(I3), I5, 1X, A1, A8, F8.3, F8.3, A8, F6.2, A8, I7, 
     &    A7, I7, 2A7, I7, 1X, F8.0, 2X, A, A1)
 8082 FORMAT( 1X, 2(I3), I5, 1X, A1, A8, F8.3, F8.3, A8, F6.2, A8, I7, 
     &    A7, I7, 2A7, I7, 1X, F8.0, 2X,F7.2, A,  A1)
 8089 FORMAT(' Superior computation ending ', I10, 3(I3), I5,
     &   ' BOP/mean/EOP elevs: ',3F10.4,' mean outflow=',F9.1, A1)
C
  614 FORMAT(' Exceeded max of ', I5,' iterations for increment ',I2,A1)
  605 FORMAT(' Superior increment ',I5, 
     &            ' beginning- and end-of-increment levels: ',2F14.8,A1)
  606 FORMAT(' Superior increment ',I5, 
     &            '    mean outflow and mean level for inc: ',2F14.8,A1)
  612 FORMAT(I3, ' Required ', I5, 
     &               ' iterations to converge within ',F7.3,' m3s.',A1)
 8033 FORMAT(' End day: ',I9,' increment:', I3, A5,'Check=',F9.1,
     &   '  vs ',F10.2, A1)
 8071 FORMAT(' Successfully solved increment ', I5, 
     &   ' of routing interval ending on day number', I9, A1)
 8072 FORMAT(' Suspicious Superior routing for end day ',I9,
     &    ' increment', I5, ' (ending d/m/y= ',I2,'/',I2,'/',I5,')',A1)
C
C
      END
C----&------------------------------------------------------------------
C
C
      DOUBLE PRECISION FUNCTION NonCompWorksFlow( SupLevel, MHLevel, 
     &   MonthNum )
C
C     Returns the "Side Channel" flow (everthing that doesn't go
C     through the Compensating Works.
C
      IMPLICIT NONE 
C
C     Declarations regarding subroutine arguments:
C      
      DOUBLE PRECISION SupLevel, MHLevel
C                   (Input)  The levels of lakes Superior and Huron.
C
      INTEGER MonthNum
C             (Input) Month to check the side channel capacity for
C
C
C     Declarations regarding local variables or functions:
C
C
C
C     Declarations regarding COMMON block variables:
C
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      CHARACTER     CTMP*255,CWRK*80, LLT*1
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
      INTEGER SideChanMaxCapacity(12), SideChanMinCapacity(12)
      DOUBLE PRECISION CWOvertopElev, CWMinimumFreeboard(12)
C                   See descriptions in Block Data Genral
C
C  Common block statements.
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      
      COMMON /SPOTOP/ SideChanMaxCapacity, CWOvertopElev, 
     &                  CWMinimumFreeboard, SideChanMinCapacity        
C
C
C     Check to see if a maximum side channel has been specified and return it
C     XXX This does not use the Bruxer, at al. low lake level estimates
      IF (SideChanMaxCapacity(MonthNum) .LT. 9999) THEN
         NonCompWorksFlow = SideChanMaxCapacity(MonthNum)
      ELSE
C        No max side channel specified so use Plan 77A default
C        NOTE: This may be different from the current operational limit due to
C        changes at the hydropower plants or other structures, but is what was 
C        used for development of Plan 1977-A.  (Well, it was actually 
C        82,000cfs, which rounds to 2320cms.)
         NonCompWorksFlow = 2320.D0
      END IF
C
      IF( SupVerbosity .GE. 9) THEN
         WRITE( CTMP, '(2(A,F8.4), A, F8.3, A1)') ' With Superior at ',
     &      SupLevel, ' and Michigan-Huron at ', MHLevel, 
     &      ', non-gated flow= ', NonCompWorksFlow, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C  
      END
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
      SUBROUTINE SupClosestFlow( TargetFlow, SupBOPLevel, MHLevel, 
     &   NTS, dT, LevelRnd, FlowRnd, ClosestFlow, Gates, SupEOPLevel,
     &   MonthNum)
C
C     This module takes a specified flow (TargetFlow) and figures out
C     the gate setting that yields the closest period-average flow.
C
      IMPLICIT NONE 
C
C     Declarations regarding subroutine arguments:
C      
C
      DOUBLE PRECISION TargetFlow, SupBOPLevel, MHLevel
C                   (Input)  The desired flow and levels of lakes
C                   Superior and Michigan-Huron.  MHLevel isn't really
C                   used at this time, so put in anything you want
C                   when calling this function.  Just have it in place
C                   in case need to pass MH level to NonCompWorksFlow
C                   because made sidechannel flow a function of
C                   effective head.
C
      DOUBLE PRECISION NTS, dT
C                   (Input)  The Net Total Supply, and length of the
C                   period in seconds.
C
      INTEGER       LevelRnd, FlowRnd, MonthNum
C                   (Input)  What power of ten the levels and flows
C                   must be rounded to.  For compatibility with
C                   traditional practice, or looking at better ways.
C                   MonthNum is needed as an argument for P77ASolveSup.
C
      REAL          Gates
C                   (Output) The number of gates that results in flow
C                   closest to TargetFlow
C
      DOUBLE PRECISION ClosestFlow, SupEOPLevel
C                   (Output) The period-average flow and End-of-Period 
C                   level resulting from Gates.
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER       Setting
      DOUBLE PRECISION SettingFlow, MeanSupLevel
C                   Scratch variables
C
C      DOUBLE PRECISION NonCompWorksFlow
C                   Function that returns that flow for a given
C                   lake level, except for the amount going through
C                   the Compensating Works.
C
C
C     Declarations regarding COMMON block variables:
C
      REAL          SupGatesOpen
C                   The number of gates open in current routing period.
      DOUBLE PRECISION SPK,SPYM,SPA,SPB,SPWT,SPC
C                   See descriptions in CHECKPAR and OUTFLOW
C
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      CHARACTER     CTMP*255,CWRK*80, LLT*1
      INTEGER       GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &              SupVerbosity, MLVerbosity, OntVerbosity
C                   See descriptions in Block Data Genral
C
      DOUBLE PRECISION SGX(17), SGY(17), SGZ(17)
C                   SGX, SGY, SGZ are constants for gate equations.  
C
      INTEGER       MaxSGSetting
C                   The number of elements in SGX, SGY, and SGZ
C
C
C     Common block statements.
C
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /SPDAT/    SPK, SPYM, SPA, SPB, SPWT, SPC, 
     &                  SGX, SGY, SGZ, MaxSGSetting,
     &                  SupGatesOpen
C
C
      SettingFlow = 0.0D0
C
      DO  3923, Setting = 1, MaxSGSetting
C
C        Save flow from previous iteration in ClosestFlow
         ClosestFlow = SettingFlow
C
C        Figure out how many gates to pass to P77ASolveSup. Use variable
C        Gates for now, although its value will be reassigned at end.
         Gates = Setting - 1
         IF( Setting .EQ. 1) Gates = 0.5
C
C        Compute flow for this gate setting
         CALL P77ASolveSup(SupBOPLevel, NTS, dT, Gates, FlowRnd,   
     &      LevelRnd, SupEOPLevel, SettingFlow, MeanSupLevel, MonthNum)
C
C        Jump out if we get a higher flow 
 3923    IF( SettingFlow .GE. TargetFlow ) GOTO 5558
C
C
C     Setting is > MaxSGSetting (unable to pass TargetFlow)
      Setting = MaxSGSetting 
      GOTO 6232
C      
C     Pass through on Setting at minimum number (TargetFlow too small for gates)
C
 5558 IF( Setting .GT. 1 .AND. DABS(SettingFlow - TargetFlow)
     &   .GT. DABS(ClosestFlow - TargetFlow)) THEN
C        Since Setting > 1, Setting was incremented so that ClosestFlow
C        currently less than TargetFlow, and SettingFlow is greater than 
C        TargetFlow.  If ClosestFlow is closer, then need to decrement 
C        Setting and assign previous value (currently in ClosestFlow) to
C        SettingFlow.
C
         SettingFlow = ClosestFlow
         Setting = Setting - 1
C        Also need new EOP level 
         CALL P77ASolveSup(SupBOPLevel, NTS, dT, -1.0, FlowRnd,   
     &      LevelRnd, SupEOPLevel, ClosestFlow, MeanSupLevel, MonthNum)
      ENDIF
         WRITE( CTMP, 6557) TargetFlow, ClosestFlow, setting, LLT
         CALL LogIt(LOGD,0,0)
C
C     The flow and gate settings are figured out now.
 6232 ClosestFlow = SettingFlow
      Gates = Setting - 1
      IF(Setting .EQ. 1) Gates = 0.5
C
C
C
C      
      IF( SupVerbosity .GE. 5) THEN
         WRITE( CTMP, 6555) SupBOPLevel, MHLevel, LLT
         CALL LogIt(LOGD,0,0)
         WRITE( CTMP, 6556) TargetFlow, ClosestFlow, Gates, LLT
C         WRITE( CTMP, 6556) ClosestFlow, Gates, SGX(Setting) * 
C     &      (SGY(Setting) * SupLevel - SGZ(Setting)) ** 1.5D0,
C     &      NonCompWorksFlow(SupLevel, MHLevel),  LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
 6555 FORMAT(10X,'==> With Superior at ', F8.4, 
     &   ' and Michigan-Huron at ', F8.4, ',', A1)
 6556 FORMAT(10X,'==> the closest flow to ', F11.4,' m3s is ', F9.4, 
     &   ' m3s (', F5.2, ' gates)', A1)
 6557 FORMAT(10X,'==> the closest flow to ', F11.4,' m3s is ', F9.4, 
     &   ' m3s (', I5, ' gates setting)', A1)    
C 6556 FORMAT(10X,' is ', F9.4, ' (', F6.2, ' gates=', F11.6, 
C     &   ' sidechannel flow=', F8.3, ')', A1)
      END
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
C
      DOUBLE PRECISION FUNCTION USS (SUPFLOW, MHUMEAN)
C
C     Compute the level at the US Slip Gauge by method of bisection of GLERL
C
C     Lifted out of Env Canada Plan 77A model, with considerable mutation
C  
C
      INTEGER SUPFLOW
C     Must be to 10m3s if MHUMean is to cm, or won't converge.
C
      DOUBLE PRECISION MHUMEAN, UPLIMIT, LOLIMIT
      INTEGER VALUE, Tmp, I
C     DOUBLE PRECISION ROUND
C
C     initial limits are Uplimit = Superior max. Lolimit= current MHU
C     and USS = a typical US slip level
C
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      CHARACTER     CTMP*255,CWRK*80, LLT*1
      INTEGER       GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &              SupVerbosity, MLVerbosity, OntVerbosity
C                   See descriptions in Block Data Genral
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /WKCHAR/   CTMP,LLT,CWRK
C
      IF (MHUMEAN.LE.(172.30D0)) THEN
         USS = DBLE(-9999)
         RETURN
      ENDIF 
C
C     SUPFLOW = IDNINT(Round( DBLE(SupFlow), 1 ))
C     SupFlow already integer
      UPLIMIT = 183.8D0
      LOLIMIT = MHUMEAN
      USS = 176.5D0
C
C      Write(*, '(I8,4F9.3)') SUPFLOW, MHUMEAN, USS, UpLimit, LoLimit
C
C     Repeat to achieve tolerance of 1 cm. SUPFLOW in 10m3s
C
      DO 11, I = 1, 99
C
C         Tmp = IDNINT(Round(USS, -2) * 100)
         Tmp = IDNINT(USS * 100)
         USS = (UpLimit + LoLimit) / 2 
C
         IF( SupVerbosity .GE. 8) THEN
            WRITE(CTMP,'(A,I3,2(A,F12.6),A1)') 
     &         '  inside USS, iteration=', I,'  USS=', USS, 
     &         ' MHUMEAN=', MHUMEAN, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C        If equation goes invalid, just asusme USS = MH + 6cm
         IF( (USS-MHUMEAN-0.06D0) .LE. 0D0 .OR. 
     &      (USS - 172.30D0) .LE. 0D0 ) THEN
            USS = MHUMEAN + 0.06D0
            IF( SupVerbosity .GE. 3) THEN
               WRITE(CTMP,'(2(A,F12.6),A,I5,A1)') '  Assuming US Slip=',
     &            USS, ' because equation invalid with MH=', MHUMEAN, 
     &            ' and St. Marys=', SUPFLOW, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
            GOTO 12
         ENDIF
C
C         VALUE = IDNINT( Round(347.3D0 * 
C     &      ((USS-MHUMEAN-0.06D0)**0.4D0)*((USS - 172.30D0)**1.5D0), 1))
         VALUE = IDNINT( 347.3D0 * 
     &      ((USS-MHUMEAN-0.06D0)**0.4D0)*((USS - 172.30D0)**1.5D0))
C
C        Write(*, '(2I8,3F9.3)') Tmp, Value, USS, 34.73D0 * 
C    &      ((USS-MHUMEAN-0.06D0)**0.4D0)*((USS - 172.30D0)**1.5D0)
C
C         VALUE = IDNINT(Round( 34.73D0 * 
C     &      ((USS-MHUMEAN-0.06D0)**0.4D0)*((USS - 172.30D0)**1.5D0)
C     &      , 0)) - SUPFLOW
C
         IF (VALUE .GE. SUPFLOW) UPLIMIT=USS
         IF (VALUE .LT. SUPFLOW) LOLIMIT=USS
C
C   11    IF ( IDNINT(Round(USS,-2) * 100) .EQ. Tmp
C     &      .AND. VALUE .EQ. SUPFLOW) GOTO 12
   11    IF(IDNINT(USS * 100) .EQ. Tmp .AND. VALUE .EQ. SUPFLOW) GOTO 12
C
C  10 MDPNT = MIDPOINT
C     MIDPOINT=(UPLIMIT+LOLIMIT)*0.5
C       VALUE =  DBLE(34.73) * ((MIDPOINT-MHUMEAN-0.06)**0.4)*
C     &          ((MIDPOINT-DBLE(172.30))**1.5)-SUPFLOW
C       IF (VALUE.GE.DBLE(0.0)) UPLIMIT=MIDPOINT
C       IF (VALUE.LT.DBLE(0.0)) LOLIMIT=MIDPOINT
C       IF ((ROUND(MIDPOINT,DBLE(100)).NE.ROUND(MDPNT,DBLE(100))).AND.
C     &       (VALUE.NE.DBLE(0.))) GOTO 10
C      ENDIF 
C
C       USS = MIDPOINT
C
         STOP ' Gilligan'
C
   12 END
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
C
      SUBROUTINE GetSupPreProject( SupBOP, SupNTS, PerLen, Verb, Mon, 
     &   PreProjFlow, PreProjMean)
C
C
C
      DOUBLE PRECISION SupBOP, SupNTS, PerLen
C              (Input) The Superior beginning of period level, the net 
C              of the non-outflow water balance terms, and the number 
C              of seconds in the routing period.
C
      INTEGER  Mon
C              (Input) The month that the flow is computer for.
C
      INTEGER  Verb
C              (Input) The level of debugging output to write.
C
      DOUBLE PRECISION PreProjFlow, PreProjMean
C              (Output) The mean St. Mary's flow and Superior level 
C              for the routing period.
C
      INTEGER  retard(12)
C              monthly retardation values
C
      INTEGER  incs
C              Number of increments to integrate storage change over
C
      INTEGER  MaxIterations
C              Upper limit for solving simutaneous equations
C
      INTEGER  i, k
      DOUBLE PRECISION Qi, BOI, EOI, LastFlow, dT, DScratch
C
C
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      CHARACTER    CTMP*255, CWRK*80, LLT*1
C                  See descriptions in Block Data Genral
C
      DOUBLE PRECISION GetArea, tmpA
C                      Returns area of lake
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
C
C
C     From:           	David_Fay@pch.gc.ca
C     Date sent:      	Mon, 2 Oct 2000 14:03:12 -0400
C     Just talked to Chuck, who apparently did some of the work with Nanette. 
C     He pointed me to the June 1970 report of the Coordinating Committee on
C     Superior flows.   It gives an estimated average Jan-Apr ice retardation
C     for 1871 -1887 to be 4 tcfs or about 110 m3/s based on the differences
C     between the Marquette and SWP discharge equations.  They note that this
C     may include effects of change in slope between Marquette and SWP as well
C     as ice in the upper river.
C
C     From:	Woodruff, Carl L LRE 
C     Sent:	Tuesday, October 03, 2000 9:24 PM
C     Subject:	Ice retardation factor
C     As I now understand it after talking with Nanette the ice retardation / 
C     Criteria "c" concern arose from what was the ice retardation that existed
C     prior to 1887 (Pre-project).  With regulation ice retardation is not a 
C     concern as the outflow is deetermined for the coming month and it is 
C     flowed.  During preproject times there were no structures and natural 
C     conditions prevailed which could have included restrictions due to ice 
C     and weed retardation in the St. Marys River.
C
C     Quoting from  Paragraph 21, Page 9 of"Lake Superior Outflow 1860 - 1968" 
C     by the Coordinating Committee on Great lakes Basic Hydraulic and 
C     Hydrologic Data, June 1970:
C
C     "21.  Effects of ice and weed retardation on stage discharge 
C     relationships.   The stage-discharge relationships at Southwest Pier are 
C     considered to be unaffected by the presence of ice and weeds in the 
C     relatively short reach from the gauge to the head of the rapids.  For the 
C     1860-1870 period, an average winter flow retardation due to ice of 4,000 
C     cfs was applied to the January through April flows from the Marquette 
C     stage-discharge relationship.  The 4,000 cfs retardation is the average 
C     difference between winter discharges from recorded Marquette and Southwest 
C     Pier stages in the 1871-1887 period and the appropriate stage-discharge 
C     relationships.  Thus the retardation used includes the effects of such 
C     factors as slope between Marquette and the rapids as well as the effects 
C     of ice.  The May through December relationship between Marquette and the 
C     discharge, computed from the Southwest Pier equation, contains any weed 
C     effect that occurs."
C
C     But it seems to make a big difference using the retardation, so leave
C     set to zero for now.
      DATA retard / 12 * 0  /
C
      DATA MaxIterations / 999 /
C
C     Use DScratch just to hush up SparcCompiler caution about integer division
      DScratch =  PerLen / DBLE(43200) 
      K = IDINT( DScratch )
      IF(K .LE. 1) THEN
C        Trial-and-error experience shows that time slices finer than
C        about 12 hours don't improve the Superior change-in-storage 
C        computation.  But if interval is > than 12 hours, just use it.
         incs = 1          
         dT = PerLen
C        
      ELSE
C        PerLen given as argument is longer than 1/2 day.  Find out the 
C        least number of increments contained within PerLen so that they 
C        are no more than 12 hours each.          
         IF( DMOD( PerLen, 43200.0D0 ) .EQ. 0.0D0 ) THEN
C            PerLen is an whole multiple of 12 hours
             incs = k   
         ELSE
C            PerLen is not a whole multiple of 12 hours; round up
             incs = k + 1   
         ENDIF
         dT = PerLen / incs
      ENDIF
C
      IF( Verb .GE. 5) THEN
         WRITE(CTMP, '(A,I3,A,F13.1, A,A1)') ' Integrating Superior '//
     &      'pre-project outflow over ',incs, ' increments of ', dT, 
     &      ' seconds each', LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, '(A,F10.5,A,F10.5,A,F13.3,A1)')  ' BOP=', SupBOP, 
     &      '  net water balance=', SupNTS, '  period len=', PerLen, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
      PreProjFlow = 0.0D0
      PreProjMean = SupBOP / 2
C
      EOI = SupBOP
C
      DO 1350 I=1, incs            
C
C        Initialize BOI levels as EOI levels from prev increment
         BOI = EOI
C
C        Dump some debugging info if desired:
         IF(Verb .GE. 8 ) THEN
            WRITE(CTMP,'(A,I3,3(A,F15.8),A,F15.2,A1)') 
     &         ' Increment ', I, ' BOI=', BOI, ' EOI=', EOI,
     &         ' MLV=', (EOI + BOI)/2, ' Area=', 
     &         (EOI + BOI)/2 , LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
C        Will do iterative solution of 2 simultaneous eqns:
C           mean level for increment = (SPBOI + SPEOI)/2
C           SPOutI = f(SPBOI + SPEOI)/2
C           (SPEOI-SPBOI)*Area/dT = Net Total Supply - outflow
         DO 1310, K = 1, MaxIterations
C           Save previous outflow for comparison to tolerance
            LastFlow = Qi
C       
C           Approximate the outflow for the iteration
            Qi=824.7D0 * ((EOI + BOI)/2 - 181.43D0)**1.5D0 - retard(mon)
C                 
            tmpA = GetArea(1,(EOI + BOI)/2)
            EOI = BOI + (SupNTS - Qi) * dT / tmpA
C
C           Dump some debugging info if desired:
            IF(Verb .GE. 9 ) THEN
               WRITE(CTMP,'(2(A,I3),3(A,F15.8),A1)') ' Increment ',
     &            I, ' Iteration ', K, ' Approx. mean=', (EOI + BOI)/2, 
     &            ' Approx. outflow=', Qi, ' delta S=', 
     &            (SupNTS - Qi) * dT / tmpA, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
C
 1310       IF( DABS(Qi - LastFlow) .LE. 0.000001D0) GOTO 1330
C
C        If fall to here, must have exceeded MaxIterations, so
C        failed to converge.
         STOP 'Troy McClure'
C
C 
 1330    IF( Verb .GE. 8) THEN
            WRITE(CTMP,605)I, BOI, EOI, LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,606)I, Qi, (EOI + BOI)/2, LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,612)I, K, 0.000001D0, LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,607)I, PreProjFlow, PreProjMean, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C            
C        Accumulate incremental results
         PreProjMean = PreProjMean + EOI
 1350    PreProjFlow = PreProjFlow + Qi
C
      PreProjMean = (PreProjMean - EOI / 2) / incs
      PreProjFlow = PreProjFlow / incs
C           
      IF( Verb .GE. 5) THEN
         WRITE(CTMP,608)PreProjFlow, PreProjMean, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
  605 FORMAT(' Sup pre-project increment ',I3, 
     &            ' beginning- and end-of-increment levels: ',2F14.8,A1)
  606 FORMAT(' Sup pre-project increment ',I3, 
     &            '    mean outflow and mean level for inc: ',2F14.8,A1)
  607 FORMAT(' Sup pre-project increment ',I3, 
     &            '    accumulated outflow and mean level: ',2D21.14,A1)
  608 FORMAT(' Sup pre-project outflow and mean level: ', 2F14.8,A1)
  612 FORMAT(I3, ' Required ', I3, 
     &               ' iterations to converge within ',F7.3,' m3s.',A1)
C
      END
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
C
C
C
      SUBROUTINE GetSWPier(SupLev, Month, Gates, NonGatedFlow, SWPRnd, 
     &   LS2SWP, SquawkLev, Year, MsgPrefix, SWPierLev, GateFlow)
C      
C    This function takes a Lake Superior (ie, Pt Iroquois level) and
C    a month number, determines the corresponding level at S.W. Pier,
C    and applies it to the SW Pier gate equations.  The gate equations
C    are calibrated periodically, and fairly reliable.  The original
C    relationship between Pt. Iroquois and S.W. Pier is quite old, and
C    possibly based on very limited data.  Rob Caldwell at EC has come
C    up with an improved alternative.  It incorporates total flow, 
C    requiring that the non-gated flow be specified.
C
C    As far as Matt knows, the Pt Iroquois equations are derivatives of 
C    the old SW Pier relationship and gate equations, inheriting all 
C    their flaws, and then muddling things up a bit further.
C    
C
C
      IMPLICIT NONE 
C
C     Declarations regarding subroutine arguments:
C      
      DOUBLE PRECISION SupLev
C                      (input) The Lake Superior level.
C
      INTEGER          NonGatedFlow, Month, SWPRnd
C                      (input) The month of the lake level, and the 
C                      rounding position for SW Pier levels.
C
      INTEGER          LS2SWP
C                      (input) Indicates which LS:SWP relationship to use:     
C                         1 = old Pt Iroquois / SW Pier equation
C                         2 = Caldwell's new equation of 2000
C                         3 = Caldwell's newer equation of 2008
C
      INTEGER          SquawkLev 
C                      (input) Indicates when to write messages about SW 
C                      Pier computations.  Has no effect on calculations - 
C                      just issues warnings.  When calling this routine for 
C                      intermediate results, such as trying gate settings
C                      to achieve a target flow within a regulation plan's 
C                      logic, use a low value to avoid screenfuls of 
C                      warnings.  If computing final routing results, use
C                      a high value in order to flag any compromised 
C                      results. 
C
      INTEGER          Year 
C                      (input) Just used for logging and messages.  Has no 
C                      effect on calculations - 
C
      CHARACTER        MsgPrefix*16
C                      (input) String pre-pended to messages to provide 
C                      context, such as 'Superior routing'.  The string
C                      should be terminated with LLT.  
C                        
      REAL             Gates
C                      (input) The number of gates open in Month.
C
      DOUBLE PRECISION SWPierLev, GateFlow
C                      (output) Level at S.W. Pier, and flow through
C                      compensating works. 
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER          I, flowtest, Flow2ItsAgo, strlen
C                      Scratch variables
C     
      DOUBLE PRECISION ROUND
C                      Returns value rounded to a certain accuracy,
C                      according to the old Lake Survey method.
C
      DOUBLE PRECISION SWPierGateFlow
C                      Returnd flow through gates, given SWP level
C
C
C     Declarations regarding COMMON block variables:
C
      DOUBLE PRECISION OldSWPA(12), OldSWPB(12), OldSWPC
      DOUBLE PRECISION EC2000SWPA, EC2000SWPB, EC2000SWPC(12)
      DOUBLE PRECISION EC2008SWPA, EC2008SWPB, EC2008SWPC(12)
C                      Constants for SW Pier/Pt Iroquois relationship
C
      DOUBLE PRECISION SWPGateEqA(16), SWPGateEqB(16)  
C                      Constants for SW Pier gate equation, used for
C                      gated flow calculation.
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
      CHARACTER*3  MonthAbbrev(12)
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      CHARACTER    CTMP*255, CWRK*80, LLT*1
C                  See descriptions in Block Data Genral
C
      COMMON /SWPIER/   SWPGateEqA, SWPGateEqB, OldSWPA, OldSWPB, 
     &                  OldSWPC, EC2000SWPA, EC2000SWPB, EC2000SWPC,
     &                  EC2008SWPA, EC2008SWPB, EC2008SWPC
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /MONAME/   MonthAbbrev
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
C
      strlen = INDEX(MsgPrefix,LLT) - 1
      IF(strlen .LT. 1) strlen = 16
C
 6899 FORMAT('         ', A, ' Entering GetSWPier() ', A3,I5, 
     &    ' with Superior level of ',F14.9,', and ', F8.1,
     &    ' gates, and sidech of ', I6, A1)
      IF ( SquawkLev .GE. 6) THEN
         WRITE(CTMP,6899) MsgPrefix(1:strlen), MonthAbbrev(Month), Year, 
     &     SupLev, Gates, NonGatedFlow, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
C     XX TODO: make a nice error message here
      IF( SupLev .LT. 100) STOP ' Bogus Sup level!'
C
C
      IF( LS2SWP .EQ. 1) THEN 
C        Using original (old) LS:SWP relationship
         SWPierLev = Round(              
     &      (SupLev + OldSWPA(Month))/OldSWPB(Month) + OldSWPC, SWPRnd)
         IF( Gates .LT. 0.45 .OR. Gates .GT. 16.1) THEN
C           Bogus number.  No gates open, or don't care.
            GateFlow = 0.0D0
         ELSE
            GateFlow = SWPierGateFlow( SWPierLev, Gates )
         ENDIF
         IF ( SquawkLev .GE. 5) THEN
            WRITE(CTMP,'(A, F10.2, A1)') '    . . . gate '//
     &         ' flow of ',
     &         GateFlow, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C        
      ELSE IF( LS2SWP .EQ. 2) THEN 
C        Using EC 2000 LS:SWP relationship
C        Solve simultaneous eqn's iteratively because new SWP eqn is
C        f(gateflow) and gateflow is f(SWP)
         IF ( SquawkLev .GE. 6) THEN
            WRITE(CTMP,'(A, I5, A1)') '    . . . using 2000 EC '//
     &         'relationship for LS:SWP and non-gated flow of ',
     &         NonGatedFlow, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
         GateFlow = 100
         Flow2ItsAgo = 101
         IF( GateFlow + NonGatedFlow .LT. 0) STOP 'Kazoo1'
         flowtest = 102
         DO 233, I = 1, 9999
            SWPierLev = Round( EC2000SWPA * SupLev -              
     &         EC2000SWPB * LOG(GateFlow + NonGatedFlow) -
     &         EC2000SWPC(Month), SWPRnd)
            IF( Gates .LT. 0.45 .OR. Gates .GT. 16.1) THEN
C              Bogus number.  No gates open, or don't care.
               GateFlow = 0.0D0
            ELSE
               GateFlow = SWPierGateFlow( SWPierLev, Gates )
            ENDIF
            IF( flowtest .EQ. IDINT(GateFlow * 100000.D0) ) GOTO 5988
C           Check to see if oscillating between two values 
            IF(SWPRnd .NE. -99 )THEN
               IF(IDINT(GateFlow*100000.D0) .EQ. Flow2ItsAgo) GOTO 5988
               Flow2ItsAgo = flowtest
            ENDIF
C            
  233       flowtest = IDINT(GateFlow * 100000.D0)
C
 5988    CONTINUE
      ELSE IF( LS2SWP .EQ. 3) THEN 
C        Using EC 2008 LS:SWP relationship
C        Solve simultaneous eqn's iteratively because new SWP eqn is
C        f(gateflow) and gateflow is f(SWP)
         IF ( SquawkLev .GE. 6) THEN
            WRITE(CTMP,'(A, I5, A1)') '    . . . using 2008 EC '//
     &         'relationship for LS:SWP and non-gated flow of ',
     &         NonGatedFlow, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
         GateFlow = 100
         Flow2ItsAgo = 101
         IF( GateFlow + NonGatedFlow .LT. 0) STOP 'Kazoo2'
         flowtest = 102
         DO 234, I = 1, 9999
            SWPierLev = Round( EC2008SWPA * SupLev -              
     &         EC2008SWPB * LOG((GateFlow + NonGatedFlow)/10.0D0) -
     &         EC2008SWPC(Month), SWPRnd)
C           if SWP level is higher than Sup level, just set it equal to Sup level
            IF (SWPierLev .GT. SupLev)  THEN
               WRITE(CTMP,'(2(A,F8.3),A,A1)') 
     &            '   Warning: computed SWPier level ', SWPierLev,  
     &            'm being changed to Superior level ',SupLev, 
     &            'm.', LLT
               IF ( SquawkLev .GE. 6) THEN
                    CALL LogIt(LOGCDE,0,0)
               ELSE 
                  CALL LogIt(LOGD,1,0)
               ENDIF
               SWPierLev = SupLev
            ENDIF
            IF( Gates .LT. 0.45 .OR. Gates .GT. 16.1) THEN
C           Bogus number.  No gates open, or don't care.
               GateFlow = 0.0D0
            ELSE
               GateFlow = SWPierGateFlow( SWPierLev, Gates )
            ENDIF
            IF ( SquawkLev .GE. 6) THEN
               WRITE(CTMP,'(2(A,F8.3),A,I7,2(A,F16.9),A1)') 
     &            '            SupLev=',SupLev,' Gates=',Gates, 
     &            ' NonGatedFlow=',NonGatedFlow,' SWPierLev=', 
     &            SWPierLev, ' GateFlow=',GateFlow, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF

            IF( flowtest .EQ. IDINT(GateFlow * 100000.D0) ) GOTO 5989
C           Check to see if oscillating between two values 
            IF(SWPRnd .NE. -99 )THEN
               IF(IDINT(GateFlow*100000.D0) .EQ. Flow2ItsAgo) GOTO 5989
               Flow2ItsAgo = flowtest
            ENDIF
C            
  234       flowtest = IDINT(GateFlow * 100000.D0)
C
 5989    CONTINUE
C        re-purpose flowtest briefly to hold total outflow
         flowtest = GateFlow + NonGatedFlow
C From:  	 Caldwell,Rob [Ontario] [Rob.Caldwell@ec.gc.ca] 	 
C Sent:  	 Tue 10/19/2010 9:55 AM
C To:  	 Dahl, Travis A LRE; Fan,Yin [Ontario]; McPherson, Matthew M HEC; 
C        Koschik, John A LRE
C Cc:  	 Fay,David [Ontario]; Allis, John T LRE
C Subject:  	 RE: Latest CGLRRM
C Well, we may be "pushing it" a bit, but it should be fine to extrapolate
C and push the envelope a tad.  Of course, I can give no indication 
C whatsoever what sort of accuracy the relation will provide given these
C conditions.  Good to see that it significantly cuts down the flaggings,
C though.
C
C Rob Caldwell
C
C -----Original Message-----
C From: Dahl, Travis A LRE [mailto:Travis.A.Dahl@usace.army.mil]
C Sent: Monday, October 18, 2010 4:57 PM
C To: Fan,Yin [Ontario]; McPherson, Matthew M HEC; Koschik, John A LRE;
C     Caldwell,Rob [Ontario]
C Cc: Fay,David [Ontario]; Allis, John T LRE
C Subject: RE: Latest CGLRRM
C
C Yin (& Rob),
C    Attached is the latest code.  The main change we made last week was
C to the reporting for the range of validity of the SW Pier equation. 
C Rob had given us ranges of Lake Superior levels from 182.87 to 183.91
C and flows from 1250 to 3740 (based on what was used to create the
C equation).  We added a "tolerance" to the error message.  Now the
C output will only complain when the levels are outside of 182.50 and
C 184.00 for the level and 1000 and 4000 for the flow.  We chose these
C numbers because they looked nice, but would appreciate any commentary
C from Rob.  (It did have the desired effect of reducing the size of the
C Details file by almost an order of magnitude.)
C      -Travis         
c         IF (SupLev .LT. 182.87D0 .OR. SupLev .GT. 183.91D0 .OR. 
c     &        flowtest .LT. 1250 .OR. flowtest .GT. 3740) THEN
C           oops - out of valid equation range, build a warning message
c            WRITE(CWRK,'(A,A,1X,A3,I5,A,A1)') '    ', 
c     &         MsgPrefix(1:strlen), MonthAbbrev(Month), Year, 
c     &         ': EC2008 SW Pier eqn out of range: ', LLT
C           re-purpose strlen
c            strlen = INDEX(CWRK,LLT) - 1
c            IF (SupLev .LT. 182.87D0 .OR. SupLev .GT. 183.91D0 ) THEN
c                WRITE(CTMP,'(A,A,F8.3,A,A1)') CWRK(1:strlen),
c     &            ' Superior level of ', SupLev, 
c     &            ' m not between 182.87 and 183.91 m', LLT
c            ELSE 
c                WRITE(CTMP,'(A,A,I6,A,A1)') CWRK(1:strlen),
c     &            ' St. Marys flow of ', flowtest, 
c     &            ' m3/s not between 1250 and 3740 m3/s', LLT
c            ENDIF
c            IF ( SquawkLev .GE. 6) THEN
c              CALL LogIt(LOGCDE,0,0)
c            ELSE 
c               IF (SupLev .LT. 182.50D0 .OR. SupLev .GT. 184.00D0 .OR. 
c     &           flowtest .LT. 1000 .OR. flowtest .GT. 4000) THEN
C                 Numbers are outside of an (arbitrarily) "acceptable" tolerance
c                  CALL LogIt(LOGD,1,0)
c               END IF
c           ENDIF
C           end of case 3 (EC2008 eqn)
c         ENDIF

      ELSE
         
         WRITE(CTMP,'(A,I3,A1)') 
     &      '     Unsupported value of LS2SWP=', LS2SWP, LLT
         CALL LogIt(LOGCDE,1,0)
         STOP 'Alice!'
      ENDIF
C

      IF ( SquawkLev .GE. 6) THEN
         WRITE(CTMP,'(A,I3,A,2F16.10,A1)') 
     &      '            . . . with level rounding at ', 
     &      SWPRnd, ' SWP level and gated flow compute to ', 
     &      SWPierLev, GateFlow, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
C  From Rob in summer 2000:
C  
C The new relationship is a linear equation with monthly mean Lake Superior
C (LS) levels and the natural logarithm of Lake Superior monthly mean 
C outflow (ln(Q)) as the only independent variants.  There are constant 
C monthly fudge factors to the basic y-intercept (-8.125 m).  I'd enter 
C those using the SWPA array, though you might want to rename it SWPC 
C for clarity.  The basic form of equation is:
C 
C SWPierLev = 1.046385*LS - 0.10015*ln(Q) + SWPC, with everything in 
C metres and 10 cms, and SWPC as:
C 
C MONTH     SWPC (m)
C Jan       -8.152
C Feb       -8.161
C Mar       -8.163
C Apr       -8.158
C May       -8.129
C Jun       -8.129
C Jul       -8.125
C Aug       -8.124
C Sep       -8.125
C Oct       -8.126
C Nov       -8.118
C Dec       -8.127
C 
C The above new and improved relation was developed on data from 1970 to 
C 1998, since there was a funky discontinuity in the set prior to that 
C (your guess is as good as ours as to what caused it) for the open 
C water months from May to November.  Soooo, the above correction 
C factors mainly account for ice(, weeds, etc), and the equation 
C shouldn't be stretched beyond it's derivation limitations (183.08<LS<183.91 m & 125<Q<374 10 cms).

C   ---------------------------------------------------------------------------------   

C     From Rob's 2008 report: DEVELOPMENT OF A REVISED S.W. PIER WATER LEVEL RELATIONSHIP
     
C     The study described above resulted in the development of a new Lake Superior: S.W. Pier water level relationship.  
C     The authors developed an updated relationship [EC2008 henceforth] using May-November data from 1979-2007:  
     
C     SWP = aSUP+bln(Q)+c            (3)
     
C     Where:      SWP      = monthly mean water level at S.W. Pier gauge (m IGLD 1985),
C     a     = 1.033806,     
C     SUP      = monthly mean lake-wide average water level of Lake Superior (m IGLD 1985), 
C     b     = -0.08622,
C     ln(Q)     = natural logarithm of monthly mean Lake Superior outflow (10 m3/s), and
C          c     = base intercept (-5.892 m) and a monthly correction factor (see Table 3 below).
               
C     As this relationship was derived using open-water conditions, the following constants, "c", are applied to account for monthly variabilities due to ice, etc.:
     
C     Table 3: Monthly Constants for EC2008 Relationship
C     Mth.     Jan          Feb          Mar          Apr          May          Jun          Jul          Aug          Sep          Oct          Nov          Dec
C     c (m) -5.916  -5.923  -5.925  -5.923  -5.897  -5.897  -5.892  -5.891 -5.894   -5.890  -5.886  -5.884
C
C     So given a mean Lake Superior level of 184.133027m for October and mean outflow of 2443.0521 m3/s, SWP =183.7952256 m 
C
C From: Fan,Yin [Ontario] 
C Sent: October 1, 2009 11:20 AM
C To: 'McPherson, Matthew M HEC'; Koschik, John A LRE; Dahl, Travis A LRE
C Cc: Watkins, Dave W HEC Visiting Scholar; Bryan Tolson; Caldwell,Rob [Ontario]; Fay,David [Ontario]; Kropfreiter, Melissa A LRE; Bill Werick; Leger,Wendy [Burlington]
C Subject: RE: Updated CGLRRM
C 
C 
C Hi, Matt,
C  
C I talked to Rob this morning regarding the ranges of the SWP equation he developed. Here is the information he provided to me:
C  
C Sup level:  182.87m to 183.91 m  (IGLD 1985)
C Sup flow:  1250 to 3740 m3/s
C The dataset he used:  1979-2007
C  
C So if the flow outsides the ranges, we have to decide what we need to do.
C 
C David suggested that if SWP level is higher than Sup level, we may just set SWP level equals to Sup level. What do you think? If you are agree with that, would you please modify the code? Thanks.
C 
C Yin
C
C
      END
C


      DOUBLE PRECISION FUNCTION SWPierGateFlow( SWPLev, Gates )
C
C     Computes gate flow from SW Pier gate equations
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      DOUBLE PRECISION SWPLev
C                      (input) The SW Pier level.
C
      REAL             Gates
C                      (input) The number of gates open
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER          I, J
C                      Scratch variables
C     
C
C     Declarations regarding COMMON block variables:
C
      DOUBLE PRECISION OldSWPA(12), OldSWPB(12), OldSWPC
      DOUBLE PRECISION EC2000SWPA, EC2000SWPB, EC2000SWPC(12)
      DOUBLE PRECISION EC2008SWPA, EC2008SWPB, EC2008SWPC(12)
C                      Constants for SW Pier/Pt Iroquois relationship
C
      DOUBLE PRECISION SWPGateEqA(16), SWPGateEqB(16)  
C                      Constants for SW Pier gate equation, used for
C                      gated flow calculation.
C                      See descriptions in Block Data Genral
      LOGICAL          LOGC(8), LOGD(8), LOGS(8), 
     &                 LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &                 LOGM(8), LOGP(8), LOGO(8), 
     &                 LOGDM(8), LOGDP(8), LOGDO(8), 
     &                 LOGSM(8), LOGSP(8), LOGSO(8),
     &                 LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &                 LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      CHARACTER        CTMP*255, CWRK*80, LLT*1
C
C
      COMMON /SWPIER/  SWPGateEqA, SWPGateEqB, OldSWPA, OldSWPB, 
     &                 OldSWPC, EC2000SWPA, EC2000SWPB, EC2000SWPC,
     &                 EC2008SWPA, EC2008SWPB, EC2008SWPC
C
      COMMON /LOGCOD/  LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                 LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                 LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                 LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/  CTMP,LLT,CWRK
C
C     Get the index to the gates array (.5 gates when I=0)
      I = NINT(Gates - 0.1) 
C
C     Test for equation being arithmetically valid
      J=I
C     Use constant for first gate if setting is 0.5 gates
      IF(I .EQ. 0) J=1
      IF( SWPLev - SWPGateEqB(J) .LE. 0)THEN
C        SWPier level too low for gate equation
         WRITE(CTMP,'(A,F12.6,A,F12.6,A,A1)') 
     &      '     SWPier level of ',SWPLev,
     &      ' too low for gate equation (b=', SWPGateEqB(J), ')',LLT
         CALL LogIt(LOGCDE,1,0)
         STOP 'GATEEQ'
      ENDIF
C    
      IF(I .EQ. 0) THEN
C        Gates must be 0.5. Use half of diff between 1 and 2 gates
         SWPierGateFlow=
     &            (SWPGateEqA(2) * (SWPLev - SWPGateEqB(2))**1.5D0 -
     &            SWPGateEqA(1) * (SWPLev - SWPGateEqB(1))**1.5D0) /2
      ELSE
         SWPierGateFlow=SWPGateEqA(I) * (SWPLev-SWPGateEqB(I))**1.5D0 
      ENDIF
C
C     XXX TODO: do not stop simualtion - define a range for where gate 
C     equations are usuable.  When outside that range define a 
C     theoretical rating based on best engineering judgement (subject to 
C     available time and resources).  Write up a Word file documenting assumptions.
C
C
C
      END
C
C
      SUBROUTINE SolveGatesSide( Gates, SideChannel, Mon, BOP, NTS, 
     &   PeriodLength, LS2SWP, EOP, MeanLevel, SWPierLevel, GateFlow,
     &   MeanTotalFlow)
C      
C     This module solves the storage equation for Lake Superior 
C     for a specific sidechannel flow and gate setting.  Mean
C     outflow and lake level computed, along with mean monthly 
C     level.
C
      IMPLICIT NONE 
C
C     Declarations regarding subroutine arguments:
C      
      REAL             Gates
C                      (input) The number of gates open.
C
      INTEGER          SideChannel, Mon, LS2SWP
C                      (input) The constant non-gated outflow, the 
C                      month number, and which Superior:SWP eqn to use.
C
      DOUBLE PRECISION BOP,  NTS,  PeriodLength
C                      (input) The BOP Lake Superior level, sum
C                      of NBS and Long Lac/Ogoki, and number of seconds
C                      in the time slice.
C
      DOUBLE PRECISION EOP, MeanLevel,MeanTotalFlow,GateFlow,SWPierLevel
C                      (output) EOP and average Lake Superior levels,
C                      corresponding SWP level, gated flow, and total 
C                      St. Marys flow.
C
C     Declarations regarding local variables or functions:
C
      INTEGER          I
C                      Scratch variable
C
      INTEGER          MaxIterations, flowtest
C
      DOUBLE PRECISION GETAREA
C                      Function that returns surface area for a lake
C
      DATA MaxIterations / 100 /
      flowtest = 0
      MeanLevel = BOP
C
C     initiliaze to some reasonable value
      MeanTotalFlow = 2000
C
      DO 322 I=1, MaxIterations            
         EOP = BOP + (NTS - MeanTotalFlow) * 
     &      PeriodLength / GETAREA(1,MeanLevel)
         MeanLevel = (EOP + BOP) / 2.0D0
         CALL GetSWPier(MeanLevel, Mon, Gates, SideChannel, -99, 
     &      LS2SWP,3,-9999,'SolveGatesSide()', SWPierLevel,GateFlow)
         MeanTotalFlow = SideChannel + GateFlow
         IF( flowtest .EQ. IDINT(GateFlow * 100000.D0) ) GOTO 323
C
  322    flowtest = IDINT(GateFlow * 100000.D0)
C
C     Should never get here - there's no rounding, so should always
C     be able to converge to within 8-10 sigfigs 
      STOP ' Wilma!'
  323 END
C
C
C
C
      DOUBLE PRECISION FUNCTION CompWorksFreeboard( 
     &   Mon, SupLvl, CompWorksTopElev)
C
C     returns the amount of "freeboard" before the structure is threatened
C      
      IMPLICIT NONE 
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER          Mon
C                      (input) The month number
C
      DOUBLE PRECISION SupLvl, CompWorksTopElev
C                      (input) The Lake Superior level and the elevation
C                      at which the Compensating Works becomes threatened
C
C
C     Declarations regarding local variables or functions:
C
C      INTEGER          I
C                      Scratch variable
C
      INTEGER          SWPMethod, SWPRoundLvl
C
      DOUBLE PRECISION SWPierLevel
C                      Computed SW Pier level
C
      INTEGER          NonGatedFlow
      DOUBLE PRECISION GateFlow
      REAL             Gates
C                      Presently only required as dummy arguments
C                      to pass along to GetSWPier(), but may be
C                      used in future 
C
C     Presently hardcoded to use old relationship between SW Pier and 
C     Lake Superior instead of the one developed in 2000 by EC, since 
C     old one is simpler for this initial implementation.  Newer method 
C     is also function of St. Marys River flow, which requires 
C     decisions about the outflow, so leave that as a refinement.
      SWPMethod = 1
      Gates = -9999.0
      NonGatedFlow = -9999
C
C     Don't do any rounding of intermediate computations
      SWPRoundLvl = -99
C
      CALL GetSWPier(SupLvl, Mon, Gates, NonGatedFlow, SWPRoundLvl, 
     &   SWPMethod, 4,-9999,'CompWorksFreebrd', SWPierLevel, GateFlow)
C      
      CompWorksFreeboard = CompWorksTopElev - SWPierLevel
C
      END
C
C
      SUBROUTINE AdjustDueToCompOvertop(DayOfMon, Mon, Year, PeriodLen,
     &   Verb, SPBOP, SPNWBF, MaxSideChannel, Gates, SideChan,
     &   AdjOTFlow, CodeStr)
C
C     returns the adjusted gates and Side Channel flows in situation where 
C     Compensating Works may be overtopped
C      
      IMPLICIT NONE 
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER          DayOfMon, Mon, Year, Verb
C                      (input) The date, month number, year, and debug message level
C
      DOUBLE PRECISION  MaxSideChannel
C                      (input) Total capacity of all outflow not through gates
C
      DOUBLE PRECISION PeriodLen, SPBOP, SPNWBF
C                      (input) The number of seconds in routing period, Lake
C                      Superior beginning of period level and the sum of the
C                      water balance terms
C
      DOUBLE PRECISION SideChan
C                      (input/output) The total side channel flow when the 
C                      overtopping condition arose.  This value will be modified.
C
      REAL             Gates
C                      (input/output) The total gates open when the overtopping
C                      condition arose.  This value will be modified.
C
      CHARACTER*40     CodeStr
C                      (input/output) The codes used for this regulation decision, 
C                      which will be clobbered if reg plan is overridden
C
C     Declarations regarding local variables or functions:
C
      INTEGER          NumGates
C                      Scratch variable
C
      DOUBLE PRECISION PreProjFlow, tmplev, AdjOTFlow,
     &                 tmpEOP, tmpMLvl, tmpSWP, gatedFlow, totalFlow
C
C     COMMON Blocks
      CHARACTER*3  MonthAbbrev(12)
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      CHARACTER    CTMP*255, CWRK*80, LLT*1
      INTEGER       MaxSGSetting
      DOUBLE PRECISION SPK,SPYM,SPA,SPB,SPWT,SPC,SGX(17),SGY(17),SGZ(17)
      REAL      SupGatesOpen
C                  See descriptions in Block Data Genral
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /MONAME/   MonthAbbrev
      COMMON /SPDAT/    SPK, SPYM, SPA, SPB, SPWT, SPC, 
     &                  SGX, SGY, SGZ, MaxSGSetting,
     &                  SupGatesOpen
C
C     Dump input to log if desired  
      IF( Verb .GE. 3) THEN
         WRITE(CTMP,'(A,I2,A3,I5,A,I2, A1)') 
     &      '    AdjustDueToCompOvertop: ',DayOfMon,MonthAbbrev(Mon),
     &      Year, '  debug lvl=', Verb, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,'(A,F12.6,A,F11.4,A,F12.3,A,F12.3,A,A1)') 
     &      '         Superior BOP ',SPBOP, 'm,  net inflow: ',
     &      SPNWBF, ' m3/s,   routing period length: ', PeriodLen,
     &      ' seconds   Maximum Side Channel: ' , MaxSideChannel, 
     &      ' m3/s', LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,'(A,F6.2,A,F11.4,A,A1)') 
     &      '         beginning gates ',Gates, ',  sidechannel ',
     &      SideChan, ' m3/s', LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
C
C From: Fan,Yin [Ontario] [mailto:Yin.Fan@ec.gc.ca]
C Sent: Wednesday, December 05, 2007 1:25 PM
C  To: McPherson, Matthew M HEC
C  Cc: Fay,David [Ontario]; Kropfreiter, Melissa A LRE; Dahl, Travis A LRE
C  Subject: 
C
C Matt,
C  I tried some tests by changing P77AMHRangeFact, I got very high Lake Superior 
C  levels (over 185m) while Lake Superior outflow was less than 2000 cms. As we 
C  know the top of the lowest 4 gates of Comp Work is only 184.0 m, so the flow 
C  and Lake Superior level are not correct. I talked with David about this. We 
C  would like to know if it is possible to code something in CGLRRM like if Lake 
C  Superior level is greater than 184.0 m than the flow should be not less than 
C  preproject flow. To develop our fencepost plan we want to use Lake Superior as 
C  a reservoir for Lake MH, but the levels and flows still need to be physically 
C  real.
C
C  Experimental results indicate that using Criterion B  would allow a greater 
C  outflow than the pre-project flow used here. The implementation here based on 
C  pre-project flows is an arbitrary assumption. (from Yin Fan, 8 Oct 2008)
C
C     First determine pre-project flow.  
C  XX  todo: verify this is using same SPNWBF as criterion checks
C     Verbosity hardcoded because we really don't need all the messages in sub
      CALL GetSupPreProject( SPBOP, SPNWBF, PeriodLen, 1, Mon, 
     &   PreProjFlow, tmplev)
C
C     Establish value of mean outflow if already above 
C     XXHardcoded to use old SWP equations for now instead of yr 2008 eq
C     XXNeed to change this to use the function.
      CALL SolveGatesSide( Gates, IDNINT(SideChan), Mon, SPBOP, SPNWBF,
     &   PeriodLen, 3, tmpEOP, tmpMLvl, tmpSWP, gatedFlow, totalFlow)
C
      IF( Verb .GE. 3) THEN
         WRITE(CTMP,'(3(A,F11.6),A,F11.4,A,F11.4,A,A1)') 
     &      '        PreProjFlow: ',PreProjFlow, 
     &      ' m3/s, Sup Elev: ', tmpMLvl, ' m, S.W. Pier: ',tmpSWP, 
     &      ' m,  period-ave regulated flow ',totalFlow, 
     &      ' m3/s,  period-ave gated flow ',gatedFlow, ' m3/s',LLT
         CALL LogIt(LOGD,0,1)
      ENDIF
C
      AdjOTFlow = totalFlow
C     So if total expected flow is less than pre-project, raise sidechannel
C     to max and open gates until flow higher than pre-project.
C
C     See if already flowing greater than pre-project
      IF( PreProjFlow .LE. totalFlow) THEN
C
         WRITE(CTMP,'(2(A,F12.6),A,A1)') 
     &      '         Regulated flow ', totalFlow, 
     &      ' m3/s is already higher than Pre-project flow ', 
     &      PreProjFlow, 
     &      ' m3/s, so no adjustments made to avoid overtopping.', LLT
         IF( Verb .GE. 3) THEN
            CALL LogIt(LOGD,0,1)
            IF( Verb .GE. 7) CALL LogIt(LOGCE,0,0)
         ENDIF
C
      ELSE
C        Somehow need to raise outflow above Pre-Project
         CodeStr(32:34) = ' OT'
C
         IF( PreProjFlow .LE. (gatedFlow + MaxSideChannel) ) THEN
C           Simply increasing sidechannel is enough.  Leave
C           gates as before, but bump up SideChannel so that
C           total flow equals Pre-Project flow
            SideChan = PreProjFlow - gatedFlow    
            WRITE(CTMP,'(A,F6.2,A,F12.6,A,A1)') 
     &         '         kept ', gates, 
     &         ' gates but increased sidechannel flow to ', SideChan, 
     &         ' m3/s to avoid overtopping.', LLT
            IF( Verb .GE. 3) THEN
               CALL LogIt(LOGD,0,1)
               IF( Verb .GE. 7) CALL LogIt(LOGCE,0,0)
            ENDIF
C
            AdjOTFlow = SideChan + gatedFlow
         ELSE
C           Set to sidechannel to max and opening gates until get 
C           above preproj flow
            SideChan = MaxSideChannel
C           
C           0.5 gates becomes 0 for purposes of incrementing to next
            DO 352, NumGates = NINT(gates - .1) + 1, MaxSGSetting - 1
               CALL SolveGatesSide( REAL(NumGates), 
     &             IDNINT(MaxSideChannel), 
     &            Mon, SPBOP, SPNWBF, PeriodLen, 3, tmpEOP, tmpMLvl, 
     &            tmpSWP, gatedFlow, totalFlow )
C
               IF( Verb .GE. 3) THEN
                  WRITE(CTMP,'(A,I2,3(A,F11.4),A,A1)') 
     &                '         trying ',Numgates, 
     &                ' gates and sidechannel =', SideChan,
     &                ' m3/s:  gateflow=', gatedFlow, 
     &                ' m3/s,  totalflow=', totalFlow, ' m3/s',LLT
                  CALL LogIt(LOGD,0,0)
               ENDIF
C
               IF( totalFlow .GE. PreProjFlow ) THEN
C                 Opened enough, now set gates and jump out 
                  Gates = NumGates
                  GOTO 367      
               ENDIF
C
  352       CONTINUE
C              
C           At this point Comp works can no longer pass the pre-project flow.  The 
C           physical situation will be messy with power facilites shut down and large
C           areas inundated and drowned outlets.  The outflow under these conditions 
C           is beyond the ability of present information to effectively model. Per 
C           discussion among Matt, Rob, and Yin on 20Jan2011, CGLRRM will default the 
C           outflow to the greatest pre-project value still within the capacity of the 
C           gate ratings and full sidechannel.  Calculations show intersection of 
C           capacity curve for 16 gates plus sidechannel and 1887 equation to be
C           4630 m3/s at a Superior level of 184.59m (Rob by spreadsheet in 08Oct2008 
C           and by confirmed by hand in 2011 and also by Yin's model experimentation
C           in 2011). 
C           
C           Modelers should note this situation as a period of unreliable results. 
C
C           Note from Matt: oops - we can't actually set the outflow to an arbitrary 
C           value here - can only set gates and sidecahnnel.  So why don't we just 
C           set outflow to be the Comp Works capacity and stay as close to pre-project 
C           as possible?
C
            WRITE(CTMP,'(A,I2,A,I5,A,F8.3,2(A,F7.1),A,A1)') 
     &         '  Serious problem on ', 
     &         DayOfMon, MonthAbbrev(Mon),Year,
     &         ', Lake Superior overtopping at ', tmpMLvl,
     &         ' m, but only able to pass Comp Works capacity of ', 
     &         totalFlow, ' m3/s instead pre-project flow of ',
     &         PreProjFlow, ' m3/s',LLT
     &         
C
            CodeStr(32:34) = ' OC'

            IF( Verb .GE. 3) THEN
               CALL LogIt(LOGCDE,0,0)
            ELSE
               CALL LogIt(LOGD,0,0)
            ENDIF

C
  367       WRITE(CTMP,'(A,F6.2,A,F11.4,A,A1)') 
     &         '         adjusted to ',Gates, 
     &         ' gates and sidechannel=',SideChan, ' m3/s', LLT
            CALL LogIt(LOGD,0,1)
            AdjOTFlow = totalFlow
         ENDIF
C
      ENDIF
C     
      END
C
C
C
      SUBROUTINE CheckSideChannelLimits (MaxCap, 
     &   ProposedSideChannelFlow,  
     &   RegulatedMonth, year, day, CodeStr)
C
C Sidechannel capacity:  The total sidechannel flow used by the routing model
C included all outflow other than that through an open gate where flow 
C that is a function of head.  
C 
C Two considerations were discussed, but chose not to change value:
C 1) Leakage isn't known with any confidence, but suspected to be somewhat 
C     higher than the amount given due to locks, stop logs, valves, sandstone
C     formations. 
C 2) The total US capacity is  limited by Treaty to an equal
C     share with Canada.
C 
C Physical minimums reflect lowest ability to stay operating - not leakage in
C case of shutdown.
C 
C Minimum winter sidechannel (Dec-Apr):
C 565 BRP
C 305 ES
C 260 US Gov't plant
C 10 m3/s for leakage through gates (0.7m3/s*16)
C 10 m3/s for municipal, industrial
C =1150 m3/s minimum flow to keep hydropower open.  Track flows attempted below
C 1150 m3/s to determine how often the Treaty assumption comes into play.
C Winter = Dec-Apr assumed
C 
C Minimum summer sidechannel (May-Nov)
C 40 m3/s USGP for plant power
C 90 m3/s ES to justify operation
C 130 m3/s BRP by Treaty
C 10 m3/s for leakage through gates (0.7m3/s*16)
C 10 m3/s for municipal, industrial
C =280 m3/s summer min
C 
C Minimum reasonable flow - 1200 m3/s is historical low flow - flag any lower
C flows.
C 
C Mean Superior to SW Pier relationship 
C 
C In 1989, 1969, and 1936, 10 cm difference between Superior lake-wide average
C and SW Pier - should cause a warning that the equation is out of range. 
C 5 cm in 1929, but isostatic rebound affects.
C 
C EC 2008 SWP Relationship range of conditions was Q=1250-3740 m3/s ,
C SUP=182.87-183.91 m.


C     Variables
      DOUBLE PRECISION MaxCap
      DOUBLE PRECISION ProposedSideChannelFlow, RevisedSideChannelFlow
      INTEGER RegulatedMonth, L, year, day
      CHARACTER*40     CodeStr
C                      (input/output) The codes used for this regulation decision, 
C                      which will be clobbered if reg plan is overridden
C

      INTEGER      SideChanMaxCapacity(12), SideChanMinCapacity(12)
      DOUBLE PRECISION CWOvertopElev, CWMinimumFreeboard(12)
C                  See explanations in RouteSuperior() in Superior.for
      CHARACTER*3  MonthAbbrev(12)
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      CHARACTER    CTMP*255, CWRK*80, LLT*1
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
      INTEGER      numMaxSC, numMinSCWinter, numMinSCSummer, numJacSC
C     COMMON Blocks
      COMMON /SPOTOP/   SideChanMaxCapacity, CWOvertopElev, 
     &                  CWMinimumFreeboard, SideChanMinCapacity
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /MONAME/   MonthAbbrev
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /PHYCAP/   numMaxSC, numMinSCWinter, numMinSCSummer,
     &                  numJacSC
C
C     Dump entering information
C
      IF (SupVerbosity .GT. 6) THEN
       WRITE(CTMP,'(A,F8.1,A,I2,A1)')
     &  '  Entering CheckSideChannelLimits, ProposedSideChannelFlow= ',
     &  ProposedSideChannelFlow, ' RegulatedMonth=',RegulatedMonth,LLT
       CALL LogIt(LOGD,0,1)
      ENDIF

C
C     The Larson and Dumont models sometimes resulted in sidechannel flows 
C     up to 2330 m3/s, due to complexities involving rounding and routing 
C     assumptions while computing Plan flow and gated flow at different 
C     stages of the calculations.  Be aware that capping the sidechannel 
C     flows below 2330 m3/s will cause occasional minor departures from
C     the older results (and actual practice).
C
C
C       
C     initialize as if it is OK
      RevisedSideChannelFlow = ProposedSideChannelFlow
C
      CWRK(1:1) = LLT
C
C
      IF (ProposedSideChannelFlow .GT. MaxCap) THEN
         CWRK = 'physical maximum'//LLT
         IF (MaxCap .LT. SideChanMaxCapacity(RegulatedMonth)) THEN
           numJacSC = numJacSC +1
           CodeStr(35:40) = ' SCJAC'
C XX DEBUG
           IF ( SupVerbosity .GE. 2. ) THEN
            WRITE(CTMP,'(A,F9.3,A,I6,A,F9.3,A1)') 
     &        '  Proposed SideChannel: ',ProposedSideChannelFlow, 
     &        '  Sidechannel Max: ',SideChanMaxCapacity(RegulatedMonth)
     &        ,'  Jac: ',MaxCap,LLT
            IF ( SupVerbosity .GE. 3. ) CALL LogIt(LOGC,0,0) 
            IF ( SupVerbosity .GE. 2. ) CALL LogIt(LOGD,0,0) 
         ENDIF
C XX DEBUG
         ELSE
           numMaxSC = numMaxSC + 1
           CodeStr(35:40) = ' SCMAX'  
         END IF
         
         RevisedSideChannelFlow = MaxCap

      ELSE
C
         IF (ProposedSideChannelFlow .LT. 
     &      SideChanMinCapacity(RegulatedMonth)) THEN
C           Flow too low
            RevisedSideChannelFlow=SideChanMinCapacity(RegulatedMonth)
C
            IF(RegulatedMonth .GT. 4 .AND. RegulatedMonth .LT. 12)THEN
               numMinSCSummer = numMinSCSummer + 1
               CWRK = 'minimum summer'//LLT
               CodeStr(35:40) = ' SCSMN'
            ELSE
               numMinSCWinter = numMinSCWinter + 1
               CWRK = 'minimum winter'//LLT
               CodeStr(35:40) = ' SCWMN'
            ENDIF
C
C           If very low Superior levels, the MaxCap might be as low as 
C           800 m3/s.  So need to make sure that max governs over min,
C           and can only raise the proposed sidechannel up to the MaxCap.
C           Still violating minimum SC so leave numMinSCSummer/numMinSCWinter
C           counters incremmented, but also governed by MaxCap so increment
C           numMaxSC.
            IF( SideChanMinCapacity(RegulatedMonth) .GT. MaxCap) THEN
              RevisedSideChannelFlow=MaxCap
              numJacSC = numJacSC + 1
              RevisedSideChannelFlow = MaxCap
              CodeStr(35:40) = ' SCMXX'
            ENDIF
C
         ENDIF
      ENDIF     

      L = INDEX(CWRK,LLT) - 1
C     IF ( RevisedSideChannelFlow .NE. ProposedSideChannelFlow) THEN
C     Rather than compare floats, check to see if the CHAR got changed
      IF ( L .GE. 2 ) THEN
C        sidechannel flow was revised 
         IF ( SupVerbosity .GE. 2. ) THEN
            WRITE(CTMP,'(I3,A3,I5,A,F9.3,A,A,A,F9.3,A1)') 
     &        day, MonthAbbrev(RegulatedMonth), year, 
     &        '  Proposed SideChannel flow of ',ProposedSideChannelFlow, 
     &        ' changed to ',CWRK(1:L),' of ',RevisedSideChannelFlow,LLT
            IF ( SupVerbosity .GE. 3. ) CALL LogIt(LOGC,0,0) 
            IF ( SupVerbosity .GE. 2. ) CALL LogIt(LOGD,0,0) 
         ENDIF
         ProposedSideChannelFlow = RevisedSideChannelFlow
      ENDIF     

      END
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
      SUBROUTINE SetUpSupOperationalSupplies(CompatIndex, MonthIndex, 
     &    SPNTS, SPNWB, MHNWB, SCNWB, ERNWB, MHRet, SCRet, ERRet)
C      
C     This module assigns the desired "operational" supplies and
C     retardations to the variables used for Criterion B and C
C     checks. 
C
      IMPLICIT NONE 
C
C     input arguments
      INTEGER       CompatIndex
C                   aka P77ACompat
C                   1=outflow eqn, 2=1955 mod rule, 3=plan77a
C                   -1 = virtcarl, 0 = read-in, 1 = 6-month operational
      INTEGER       MonthIndex
C                   1=Jan
      DOUBLE PRECISION SPNTS
C                   sum of Superior water balance terms for current routing 
C                   period read from tsd file.
C
C     arguments used for output
      INTEGER       SPNWB, MHNWB, SCNWB, ERNWB, MHRet, SCRet, ERRet
C                   Values that hold routing data appropriate to setting
C                   of P77ACompat.  Not used after overtopping calls.
C
C
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
      INTEGER P77ASup05NBS(12), P77ASup50NBS(12), P77ASup95NBS(12), 
     &        P77ALLOgokiDiv(12),P77AStMarysRetrd(12), 
     &        P77AMH05NBS(12), P77AMH50NBS(12), P77AMH95NBS(12), 
     &        P77AChicagoDiv(12), P77AStClairRetrd(12), 
     &        P77AStC50NBS(12), P77ADetroitRetrd(12),
     &        P77AErie50NBS(12), P77AWellandFlow(12), 
     &        P77ANiagaraRetrd(12),
     &        P77AMidLakesInc, P77AMidLakesTol, 
     &        P77AMonthRound, P77AFlowRound, P77ALevelRound
      DOUBLE PRECISION  P77ACGIP, P77ALowTrigger, P77AHighTrigger, 
     &        P77ASPK, P77ASPYM, P77ASPA, P77ASPB, P77ASPWT, P77ASPC, 
     &        P77AMHK, P77AMHYM, P77AMHA, P77AMHB, P77AMHWT, P77AMHC, 
     &        P77ASCK, P77ASCYM, P77ASCA, P77ASCB, P77ASCWT, P77ASCC, 
     &        P77AERK, P77AERYM, P77AERA, P77AERB, P77AERWT, P77AERC
      INTEGER P77AWinterMax, P77AMaxChange, P77AMinimumFlow(12)
      DOUBLE PRECISION P77ASupTargetLev(12), P77ASupRangeFact(12),
     &        P77AMHTargetLev(12), P77AMHRangeFact(12), P77AAValue
      INTEGER P77ABaseFlow(12)
      INTEGER P77AMonthsToAve(12)
      LOGICAL P77AUniformMonth
C
      COMMON /SUPCOM/   SPOPNBS06, OPLLO, MHOPNBS06, OPChicago, SCRRet,
     &                  SCOPNBS06, DetRet, EROPNBS06, OPWelland, NiaRet,
     &                  SPOPNBS07, MHOPNBS07, SCOPNBS07, EROPNBS07,
     &                  SPOPNBS08, MHOPNBS08, SCOPNBS08, EROPNBS08,
     &                  OPWelland_P2012, SCRRet_P2012, DetRet_P2012,
     &                  NiaRet_P2012,     
     &                  ForecastMonths,
     &                  LStoSWP, SupCritRnd, BackCompat, MonthsChecked

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

C           Now sort out ways to remain compatible with existing
C           procedures.  
C
            IF ( CompatIndex .EQ. -1 ) THEN
C              Virtcarl-type supplies.  Use latest operational values for
C              month 1 of the operational values for final flow average, 
C              and criterion checks
               SPNWB=OPLLO(MonthIndex)+ SPOPNBS07(MonthIndex,1)
               MHNWB = OPChicago + MHOPNBS07( MonthIndex, 1)
               SCNWB = SCOPNBS07( MonthIndex, 1)
               ERNWB = OPWelland( MonthIndex ) + 
     &                   EROPNBS07( MonthIndex, 1)
               MHRet = SCRRet(MonthIndex)
               SCRet = DetRet(MonthIndex)
               ERRet = NiaRet(MonthIndex)
C
            ELSE IF( CompatIndex .EQ. 0 ) THEN
C              Replicating existing models.  Use observed supplies read
C              in for routing for final flow average, and criterion
C              checks.
               SPNWB = IDNINT(SPNTS)
C              We actually can't replicate Dumont/Larson model methods
C              without looking ahead to this month's supplies on the
C              middle lakes.  But this only comes into play regarding 
C              Criterion B, which is a mess anyway.  So we'll just
C              fake it here, and use the first month operationals.
               MHNWB = OPChicago + MHOPNBS07( MonthIndex, 1)
               SCNWB = SCOPNBS07( MonthIndex, 1)
               ERNWB = OPWelland( MonthIndex ) + 
     &                   EROPNBS07( MonthIndex, 1)
               MHRet = SCRRet(MonthIndex)
               SCRet = DetRet(MonthIndex)
               ERRet = NiaRet(MonthIndex)
C  
            ELSE IF( CompatIndex .EQ. 1 ) THEN
C              Using Plan 1977A Supplies instead of operational supplies
               SPNWB = P77ALLOgokiDiv(MonthIndex) + 
     &            P77ASup50NBS(MonthIndex)
               MHNWB = P77AChicagoDiv(MonthIndex) + 
     &            P77AMH50NBS(MonthIndex)
               SCNWB = P77AStC50NBS( MonthIndex)
               ERNWB = P77AWellandFlow( MonthIndex ) + 
     &                   P77AErie50NBS( MonthIndex)
               MHRet = P77AStClairRetrd(MonthIndex)
               SCRet = P77ADetroitRetrd(MonthIndex)
               ERRet = P77ANiagaraRetrd(MonthIndex)
C  
            ELSE IF( CompatIndex .EQ. 2 ) THEN
               STOP ' No P77A=2' 
cC              6-month regulation run.  For first six months, use 
cC              1989 operational values for final flow average, and 
cC              criterion checks.
c               MonthsChecked = MonthsChecked + 1
c               IF(MonthsChecked .GT. 6) THEN
cC                 Period ends more than 6 months into run, so unable
cC                 to use operational supply matrix.
c                  CALL ErrMess(1146,'SuperiorRoute'//LLT//'  ',LOGCDE)
c                  CALL SAVETOD
c                  CALL ADIOS()
c               ENDIF
cC  
c               SPNWB = OPLLO(StartMonth) + 
c     &            SPOPNBS07(StartMonth, MonthsChecked)
c               MHNWB = OPChicago + MHOPNBS07(StartMonth, MonthsChecked)
c               SCNWB = SCOPNBS07( StartMonth, MonthsChecked )
c               ERNWB = OPWelland( MonthIndex ) + 
c     &                 EROPNBS07( StartMonth, MonthsChecked )
c               MHRet = SCRRet(MonthIndex)
c               SCRet = DetRet(MonthIndex)
c               ERRet = NiaRet(MonthIndex)
c               IF(SupVerbosity .GE. 5 ) THEN
c                  WRITE( CTMP, '(3A,A1)') 
c     &               '   Using operational supply matrix for Criterion',
c     &               ' checks, like real practice but not old models', 
c     &               ' (P77ACompat=2)', LLT
c                  CALL LogIt(LOGD,0,0)
c                  WRITE( CTMP, '(3(A,I13),A1,A1)') 
c     &               '     SPNWB=', SPNWB,
c     &               ' (OPLLO     ', OPLLO(StartMonth),
c     &               ' + SPOPNBS07 ',
c     &               SPOPNBS07(StartMonth, MonthsChecked), ')', LLT
c                  CALL LogIt(LOGD,0,0)
c                  WRITE( CTMP, '(3(A,I13),A1,A1)') 
c     &               '     MHNWB=',MHNWB,
c     &               ' (OPChicago ',OPChicago,
c     &               ' + MHOPNBS07 ',
c     &               MHOPNBS07(StartMonth, MonthsChecked), ')', LLT
c                  CALL LogIt(LOGD,0,0)
c                  WRITE( CTMP, '(3(A,I13),A1,A1)') 
c     &               '     ERNWB=',ERNWB,
c     &               ' (OPWelland ',OPWelland( MonthIndex ),
c     &               ' + EROPNBS07 ',
c     &               EROPNBS07(StartMonth, MonthsChecked), ')', LLT
c                  CALL LogIt(LOGD,0,0)
c               ENDIF
C
            ELSE IF ( CompatIndex .EQ. 3 ) THEN
C              Virtcarl-type supplies.  Same as option -1, but use 
C              older operational values (only for month 1 of the  
C              6-month array for final flow average, and criterion checks).
               SPNWB=OPLLO(MonthIndex)+ SPOPNBS06(MonthIndex,1)
               MHNWB = OPChicago + MHOPNBS06( MonthIndex, 1)
               SCNWB = SCOPNBS06( MonthIndex, 1)
               ERNWB = OPWelland( MonthIndex ) + 
     &                   EROPNBS06( MonthIndex, 1)
               MHRet = SCRRet(MonthIndex)
               SCRet = DetRet(MonthIndex)
               ERRet = NiaRet(MonthIndex)
C
            ELSE
               STOP 'Datsyuk'
            ENDIF
      END
C
C----&------------------------------------------------------------------

      SUBROUTINE SetUpSupOperationalSupplies_P2012(CompatIndex, 
     &   MonthIndex, SPNWB, MHNWB, SCNWB, ERNWB, MHRet, SCRet, ERRet)
C      
C     This module assigns the desired "operational" supplies and
C     retardation values to the variables used for Criterion B 
C     check in Plan2012. 
C
      IMPLICIT NONE 
C
C     input arguments
      INTEGER       CompatIndex
C                   aka P77ACompat
C                   1=outflow eqn, 2=1955 mod rule, 3=plan77a
C                   -1 = virtcarl, 0 = read-in, 1 = 6-month operational
      INTEGER       MonthIndex
C                   1=Jan

C     arguments used for output
      INTEGER       SPNWB, MHNWB, SCNWB, ERNWB, MHRet, SCRet, ERRet
C                   Values that hold routing data appropriate to setting
C                   of P77ACompat.  Not used after overtopping calls.
C
C
      INTEGER   ForecastMonths, LStoSWP, SupCritRnd, MonthsChecked
         INTEGER SPOPNBS06(12,6), SPOPNBS07(12,6), OPLLO(12), 
     &           MHOPNBS06(12,6), MHOPNBS07(12,6), OPChicago, 
     &           SCRRet(12), SCOPNBS06(12,6), SCOPNBS07(12,6), 
     &           DetRet(12), EROPNBS06(12,6), EROPNBS07(12,6),
     &           OPWelland(12), NiaRet(12),
     &           SCRRet_P2012(12), SCOPNBS08(12,6), SPOPNBS08(12,6), 
     &           DetRet_P2012(12), EROPNBS08(12,6), MHOPNBS08(12,6),
     &           OPWelland_P2012(12), NiaRet_P2012(12)       
 
      LOGICAL   BackCompat

C
      COMMON /SUPCOM/   SPOPNBS06, OPLLO, MHOPNBS06, OPChicago, SCRRet,
     &                  SCOPNBS06, DetRet, EROPNBS06, OPWelland, NiaRet,
     &                  SPOPNBS07, MHOPNBS07, SCOPNBS07, EROPNBS07,
     &                  SPOPNBS08, MHOPNBS08, SCOPNBS08, EROPNBS08,
     &                  OPWelland_P2012, SCRRet_P2012, DetRet_P2012,
     &                  NiaRet_P2012,     
     &                  ForecastMonths,
     &                  LStoSWP, SupCritRnd, BackCompat, MonthsChecked


C
            IF ( CompatIndex .EQ. -1 ) THEN
C              Virtcarl-type supplies.  Use latest operational values for
C              month 1 of the operational values for final flow average, 
C              and criterion checks
               SPNWB = OPLLO(MonthIndex)+ SPOPNBS08(MonthIndex,1)
               MHNWB = OPChicago + MHOPNBS08( MonthIndex, 1)
               SCNWB = SCOPNBS08( MonthIndex, 1)
               ERNWB = OPWelland_P2012( MonthIndex ) + 
     &                   EROPNBS08( MonthIndex, 1)
               MHRet = SCRRet_P2012(MonthIndex)
               SCRet = DetRet_P2012(MonthIndex)
               ERRet = NiaRet_P2012(MonthIndex)
            ENDIF
      END
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
      SUBROUTINE GetInstantGate(SupLevel, TargetFlow, MonthNum, 
     &   RoundDownOrUp, PreviousGates, SupGates, GateFlow,
     &   SWPierLev, SideChan, NominalStMarys)
     
C        Added PreviousGates to implement the logic that      
C        the number of gate openning should not be changed from Dec to April.           
C
C        Purpose is to figure out the gates open and corresponding sidechannel 
C        given a instantaneous Superior level and target St. Marys River flow.
C        The total river flow is also returned in the last argument, since
C        sometimes the target flow is impossible to meet, and the routine 
C        returns the closest flow possible. 
C
C        This routine assumes steady-state conditions and relies on 
C        GetSWPier() to solve simultaneous equations for calculating SW Pier
C        level from Lake Superior elevation.  The gated flow represents an
C        instantaneous solution rather than a mean value over some interval.  
C        If mean values over a period are needed, use SolveGatesSide(), which
C        includes the change-in-storage while solving the simultaneous 
C        equations (look above Chief Wiggins for example).
C
C        Uses the SWPier-Lake Superior relationship indicated by LS2SWP held 
C        in a common block.
C
C        Hardcoded to do no rounding.
C       
C

         IMPLICIT NONE 
C
C        Input args
         DOUBLE PRECISION SupLevel, TargetFlow 
c         DOUBLE PRECISION leakage, fishery, totalleak          
         INTEGER MonthNum, RoundDownOrUp
C          RoundDownOrUp: -1 means always round gates down
C                         -9999 or 0 means choose whatever gate setting is  
C                             closest match to target flow (traditional)
C                         +1 means always round up

C        Output
         DOUBLE PRECISION GateFlow, SideChan, NominalStMarys, SWPierLev
         REAL SupGates, PreviousGates
C
C        local variables
         INTEGER I
         DOUBLE PRECISION  SideChanCapacity
         DOUBLE PRECISION TmpFlow, PrevGateFlowDiff
C         DOUBLE PRECISION ROUND 
C                      Returns value rounded to a certain accuracy,
C                      according to the old Lake Survey method.
C
C        commons
         INTEGER SideChanMaxCapacity(12), SideChanMinCapacity(12)
         DOUBLE PRECISION CWOvertopElev, CWMinimumFreeboard(12)

C               These variables included just beacuse they are
C               in same COMMON as LStoSWP and SupCritRnd
         INTEGER ForecastMonths, LStoSWP, SupCritRnd, MonthsChecked
         INTEGER SPOPNBS06(12,6), SPOPNBS07(12,6), OPLLO(12), 
     &           MHOPNBS06(12,6), MHOPNBS07(12,6), OPChicago, 
     &           SCRRet(12), SCOPNBS06(12,6), SCOPNBS07(12,6), 
     &           DetRet(12), EROPNBS06(12,6), EROPNBS07(12,6),
     &           OPWelland(12), NiaRet(12),
     &           SCRRet_P2012(12), SCOPNBS08(12,6), SPOPNBS08(12,6), 
     &           DetRet_P2012(12), EROPNBS08(12,6), MHOPNBS08(12,6),
     &           OPWelland_P2012(12), NiaRet_P2012(12)       
         LOGICAL BackCompat
         LOGICAL LOGC(8), LOGD(8), LOGS(8), 
     &           LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &           LOGM(8), LOGP(8), LOGO(8), 
     &           LOGDM(8), LOGDP(8), LOGDO(8), 
     &           LOGSM(8), LOGSP(8), LOGSO(8),
     &           LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &           LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
         INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &           SupVerbosity, MLVerbosity, OntVerbosity
         CHARACTER   CTMP*255, CWRK*80, LLT*1
         CHARACTER*3 MonthAbbrev(12)
      DOUBLE PRECISION OldSWPA(12), OldSWPB(12), OldSWPC
      DOUBLE PRECISION EC2000SWPA, EC2000SWPB, EC2000SWPC(12)
      DOUBLE PRECISION EC2008SWPA, EC2008SWPB, EC2008SWPC(12)
C                      Constants for SW Pier/Pt Iroquois relationship
C
      DOUBLE PRECISION SWPGateEqA(16), SWPGateEqB(16)
         COMMON /RCPARI/ GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
         COMMON /SUPCOM/ SPOPNBS06, OPLLO, MHOPNBS06, OPChicago, SCRRet,
     &                  SCOPNBS06, DetRet, EROPNBS06, OPWelland, NiaRet,
     &                  SPOPNBS07, MHOPNBS07, SCOPNBS07, EROPNBS07,
     &                  SPOPNBS08, MHOPNBS08, SCOPNBS08, EROPNBS08,
     &                  OPWelland_P2012, SCRRet_P2012, DetRet_P2012,
     &                  NiaRet_P2012,     
     &                  ForecastMonths,
     &                  LStoSWP, SupCritRnd, BackCompat, MonthsChecked
         COMMON /SPOTOP/ SideChanMaxCapacity, CWOvertopElev, 
     &                  CWMinimumFreeboard, SideChanMinCapacity        
         COMMON /MONAME/ MonthAbbrev
         COMMON /WKCHAR/ CTMP,LLT,CWRK
         COMMON /LOGCOD/ LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
         COMMON /SWPIER/SWPGateEqA, SWPGateEqB, OldSWPA, OldSWPB, 
     &                  OldSWPC, EC2000SWPA, EC2000SWPB, EC2000SWPC,
     &                  EC2008SWPA, EC2008SWPB, EC2008SWPC
      
C
         IF ( SupVerbosity .GE. 5. ) THEN
            WRITE(CTMP,'(A,A3,2(A,F9.3),A,I2, A,I6,A1)') 
     &        '  Entering GetInstantGate()', MonthAbbrev(MonthNum),
     &        ' with Sup= ', SupLevel,'m  TargetFlow= ', TargetFlow, 
     &        ' m3/s using SWP eq=', LStoSWP,'  RoundDownOrUp=', 
     &        RoundDownOrUp, LLT
            CALL LogIt(LOGD,0,0) 
         ENDIF
C
C         Make the conditional simpler later 
         IF (RoundDownOrUp .EQ. -9999) RoundDownOrUp = 0
C
         CALL GetSideChanMaxCapacity( SupLevel, TargetFlow, MonthNum, 
     &       SideChanCapacity)
     
         IF ( SupVerbosity .GE. 5. ) THEN
             WRITE(CTMP,'(A, F9.3, A, I5,A1)') 
     &        '  PreviousGates=', PreviousGates,
     &        ' month=', MonthNum, LLT
              CALL LogIt(LOGD,0,0)
         END IF  
C        See if Target flow at least fills up the sidechannel
C        If not, solve for SW Pier using Target flow
C        Calculate SWP level according to the Beginning Sup level and Target flow
c         TmpFlow = SideChanMaxCapacity(MonthNum)
         
         TmpFlow = SideChanCapacity
         IF( TargetFlow .LT. SideChanCapacity) 
     &      TmpFlow = TargetFlow               
         
         IF ((PreviousGates .GE. 0.5) .and. (MonthNum .LE. 4) .and. 
     &       (SupLevel .LE. CWOvertopElev))  THEN  

            CALL GetSWPier(SupLevel, MonthNum, PreviousGates, 
     &         IDNINT(TmpFlow),-99,LStoSWP,SupVerbosity,-9999,
     &         'GetInstantGate()',SWPierLev, GateFlow)         
               SupGates = PreviousGates
            IF (TargetFlow .LT. GateFlow) THEN
C              Exceed the target flow, but even .5 gate is too much.  
C              Set sidechannel to 0, and return riverflow = gate flow
               SideChan = 0.0D0
               NominalStMarys = GateFlow
            ELSE                 
               SideChan = TargetFlow - GateFlow
               IF (SideChan .GE. SideChanCapacity) THEN
                   SideChan = SideChanCapacity
               END IF    
               NominalStMarys = SideChan + GateFlow
            END IF  
            GOTO 3373  
         END IF  
                  
  
         PrevGateFlowDiff = 9999
C
C        Now loop through gates until pass target flow.  Save results of 
C        previous pass in case flow exceeding the target is farther than
C        flow short of target.
         DO 3353, I = 0, 16
            IF (I .EQ. 0 ) SupGates = 0.5
            CALL GetSWPier(SupLevel, MonthNum, SupGates, 
     &         IDNINT(TmpFlow),-99,LStoSWP,SupVerbosity,-9999,
     &         'GetInstantGate()',SWPierLev,GateFlow)

           IF ( SupVerbosity .GE. 5. ) THEN
             WRITE(CTMP,'(2(A,F9.3),A1)') 
     &        ' SideChanCapacity= ', SideChanCapacity, 
     &        '  GateFlow=', GateFlow, LLT
              CALL LogIt(LOGD,0,0)
            END IF   
C           See if total flow exceeds target yet
            IF((GateFlow + SideChanCapacity ) .LT. 
     &            TargetFlow) THEN 
C              Not there yet.  Save difference and continue  
               PrevGateFlowDiff = TargetFlow -
     &           (GateFlow+SideChanCapacity)
               IF (I .EQ. 0 ) THEN
                  SupGates = 1.0
               ELSE
                  SupGates = SupGates + 1
               ENDIF
            ELSE
C              Total outflow exceeded target flow.  Done looping, but
C              need to see which case.
C              Check for weird case:
               IF (TargetFlow .LT. GateFlow) THEN
C                 Exceed the target flow, but even .5 gate is too much.  
C                 Set sidechannel to 0, and return riverflow = gate flow
                  SideChan = 0.0D0
                  NominalStMarys = GateFlow
                  WRITE(CTMP,'(3(A,F9.3),A,A1)') 
     &               ' Unable to achieve target flow of ', TargetFlow, 
     &               ' m3/s because .5 gate outflow is ', 
     &               GateFlow, ' m3/s, so setting river flow to ',
     &               NominalStMarys, ' m3/s',LLT
                  IF ( SupVerbosity .GE. 5. )CALL LogIt(LOGC,0,0) 
                  IF ( SupVerbosity .GE. 3. )CALL LogIt(LOGD,0,0) 
C
C              If we always round up, or choosing closest match and higher 
C              gate flow from this pass best matches target: 
               ELSE IF( RoundDownOrUp .EQ. 1 .OR. 
     &           (RoundDownOrUp .EQ. 0 .AND. (Prevgateflowdiff .GT.
     &           (GateFlow+SideChanCapacity - TargetFlow)))
     &           ) THEN
C                This is closest possible setting.  Trim back sidechannel until
C                we hit the target flow.  (Or we simply always round up gates.)
                 SideChan = TargetFlow - GateFlow
                 NominalStMarys = TargetFlow
C                 WRITE(CTMP,'(4(A,F9.3),A,A1)') 
C     &             'Trim sidechannel flow, target flow is ',TargetFlow, 
C     &             ' m3/s  gate outflow is ', 
C     &               GateFlow, ' m3/s, so setting river flow to ',
C     &               NominalStMarys, ' m3/s, Prevgateflowdiff',
C     &               PrevGateFlowDiff,' m3/s', LLT
C                 CALL LogIt(LOGD,0,0) 
C
               ELSE
C              Flow was closer to target at previous gate setting, or
C              always round down anyway
                 SupGates = SupGates - 1
                 IF (NINT(SupGates - .1) .EQ. 0) SupGates = 0.5
                 GateFlow = TargetFlow - PrevGateFlowDiff - 
     &              SideChanCapacity
                 SideChan = SideChanCapacity
                 NominalStMarys = GateFlow + SideChan
C                                  WRITE(CTMP,'(4(A,F9.3),A,A1)') 
C     &               ' Round down, using target flow of ', TargetFlow, 
C     &            ' m3/s using previous gated flow is ', 
C     &               GateFlow ,' m3/s,',
C     &               NominalStMarys, ' m3/s, Prevgateflowdiff',
C     &               PrevGateFlowDiff,' m3/s', LLT
C                 CALL LogIt(LOGD,0,0) 
               ENDIF
C              OK have our answers
               GOTO 3373
            ENDIF
 3353    CONTINUE    
C              
C        If fall through here, then must be max gates, but 
C        couldn't reach targetflow
         SupGates = 16
         SideChan = SideChanCapacity
         NominalStMarys = SideChan + GateFlow
         WRITE(CTMP,'(3(A,F9.3),A,A1)') 
     &      ' Unable to achieve target flow of ', TargetFlow, 
     &      ' m3/s because 16 gate outflow is ', 
     &      GateFlow, ' m3/s, so setting river flow to ',
     &      NominalStMarys, ' m3/s',LLT
         IF ( SupVerbosity .GE. 5. )CALL LogIt(LOGC,0,0) 
         IF ( SupVerbosity .GE. 3. )CALL LogIt(LOGD,0,0) 
C
 3373    CONTINUE
 
 
         IF ( SupVerbosity .GE. 5. ) THEN
            WRITE(CTMP,'(A,F9.3,A,I2,A,F6.2,2(A,F9.3),A,A1)') 
     &        '  Leaving GetInstantGate() where TargetFlow= ', 
     &        TargetFlow, ' m3/s using SWP eq=', LStoSWP,
     &        ' gates=', SupGates, ' gateflow=', GateFlow, 
     &        ' river flow=', NominalStMarys, ' m3/s', LLT
            CALL LogIt(LOGD,0,0) 
         ENDIF

      END
C----&------------------------------------------------------------------      
      
      SUBROUTINE GetForcedGates(SupLevel, TargetFlow, MonthNum, 
     &   SupGates, GateFlow, SWPierLev, SideChan, NominalStMarys)
     
C
C        Purpose is to figure out the gated flow and corresponding sidechannel flow
C        given a forced gates condition (e.g., gate maintenance/limited capacity). 
C        The total river flow is also returned in the last argument, since
C        sometimes the target flow is impossible to meet. 
C
C        This routine assumes steady-state conditions and relies on 
C        GetSWPier() to solve simultaneous equations for calculating SW Pier
C        level from Lake Superior elevation.  The gated flow represents an
C        instantaneous solution rather than a mean value over some interval.  
C        If mean values over a period are needed, use SolveGatesSide(), which
C        includes the change-in-storage while solving the simultaneous 
C        equations (look above Chief Wiggins for example).
C
C        Uses the SWPier-Lake Superior relationship indicated by LS2SWP held 
C        in a common block.
C
C        Hardcoded to do no rounding.
C       
C

         IMPLICIT NONE 
C
C        Input args
         DOUBLE PRECISION SupLevel, TargetFlow 
c         DOUBLE PRECISION leakage, fishery, totalleak          
         INTEGER MonthNum

C        Output
         DOUBLE PRECISION GateFlow, SideChan, NominalStMarys, SWPierLev
         REAL SupGates
C
C        local variables
         INTEGER I
         DOUBLE PRECISION  SideChanCapacity
         DOUBLE PRECISION TmpFlow, PrevGateFlowDiff
C         DOUBLE PRECISION ROUND 
C                      Returns value rounded to a certain accuracy,
C                      according to the old Lake Survey method.
C
C        commons
         INTEGER SideChanMaxCapacity(12), SideChanMinCapacity(12)
         DOUBLE PRECISION CWOvertopElev, CWMinimumFreeboard(12)

C               These variables included just beacuse they are
C               in same COMMON as LStoSWP and SupCritRnd
         INTEGER ForecastMonths, LStoSWP, SupCritRnd, MonthsChecked
         INTEGER SPOPNBS06(12,6), SPOPNBS07(12,6), OPLLO(12), 
     &           MHOPNBS06(12,6), MHOPNBS07(12,6), OPChicago, 
     &           SCRRet(12), SCOPNBS06(12,6), SCOPNBS07(12,6), 
     &           DetRet(12), EROPNBS06(12,6), EROPNBS07(12,6),
     &           OPWelland(12), NiaRet(12),
     &           SCRRet_P2012(12), SCOPNBS08(12,6), SPOPNBS08(12,6), 
     &           DetRet_P2012(12), EROPNBS08(12,6), MHOPNBS08(12,6),
     &           OPWelland_P2012(12), NiaRet_P2012(12)       
         LOGICAL BackCompat
         LOGICAL LOGC(8), LOGD(8), LOGS(8), 
     &           LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &           LOGM(8), LOGP(8), LOGO(8), 
     &           LOGDM(8), LOGDP(8), LOGDO(8), 
     &           LOGSM(8), LOGSP(8), LOGSO(8),
     &           LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &           LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
         INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &           SupVerbosity, MLVerbosity, OntVerbosity
         CHARACTER   CTMP*255, CWRK*80, LLT*1
         CHARACTER*3 MonthAbbrev(12)
         DOUBLE PRECISION OldSWPA(12), OldSWPB(12), OldSWPC
         DOUBLE PRECISION EC2000SWPA, EC2000SWPB, EC2000SWPC(12)
         DOUBLE PRECISION EC2008SWPA, EC2008SWPB, EC2008SWPC(12)
C                      Constants for SW Pier/Pt Iroquois relationship
         DOUBLE PRECISION SWPGateEqA(16), SWPGateEqB(16)
         COMMON /RCPARI/ GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
         COMMON /SUPCOM/ SPOPNBS06, OPLLO, MHOPNBS06, OPChicago, SCRRet,
     &                  SCOPNBS06, DetRet, EROPNBS06, OPWelland, NiaRet,
     &                  SPOPNBS07, MHOPNBS07, SCOPNBS07, EROPNBS07,
     &                  SPOPNBS08, MHOPNBS08, SCOPNBS08, EROPNBS08,
     &                  OPWelland_P2012, SCRRet_P2012, DetRet_P2012,
     &                  NiaRet_P2012,     
     &                  ForecastMonths,
     &                  LStoSWP, SupCritRnd, BackCompat, MonthsChecked
         COMMON /SPOTOP/ SideChanMaxCapacity, CWOvertopElev, 
     &                  CWMinimumFreeboard, SideChanMinCapacity        
         COMMON /MONAME/ MonthAbbrev
         COMMON /WKCHAR/ CTMP,LLT,CWRK
         COMMON /LOGCOD/ LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
         COMMON /SWPIER/SWPGateEqA, SWPGateEqB, OldSWPA, OldSWPB, 
     &                  OldSWPC, EC2000SWPA, EC2000SWPB, EC2000SWPC,
     &                  EC2008SWPA, EC2008SWPB, EC2008SWPC

C
         IF ( SupVerbosity .GE. 5. ) THEN
            WRITE(CTMP,'(A,A3,2(A,F9.3),A,I2, A,F6.2,A1)') 
     &        '  Entering GetForcedGate()', MonthAbbrev(MonthNum),
     &        ' with Sup= ', SupLevel,'m  TargetFlow= ', TargetFlow, 
     &        ' m3/s using SWP eq=', LStoSWP,'  SupGates=', 
     &        SupGates, LLT
            CALL LogIt(LOGD,0,0) 
         ENDIF
C
C
         CALL GetSideChanMaxCapacity( SupLevel, TargetFlow, MonthNum, 
     &       SideChanCapacity)
         
         TmpFlow = TargetFlow
 
         CALL GetSWPier(SupLevel, MonthNum, SupGates, 
     &         IDNINT(TmpFlow),-99,LStoSWP,SupVerbosity,-9999,
     &         'GetForcedGate()',SWPierLev, GateFlow)         
               
         IF (TargetFlow .LT. GateFlow) THEN
C              Exceed the target flow, but even .5 gate is too much.  
C              Set sidechannel to 0, and return riverflow = gate flow
               SideChan = 0.0D0
               NominalStMarys = GateFlow
         ELSE                 
               SideChan = TargetFlow - GateFlow
               IF (SideChan .GE. SideChanCapacity) THEN
                   SideChan = SideChanCapacity
               END IF    
               NominalStMarys = SideChan + GateFlow
         END IF  
 
         IF ( SupVerbosity .GE. 5. ) THEN
            WRITE(CTMP,'(A,F9.3,A,I2,A,F6.2,4(A,F9.3),A,A1)') 
     &        '  Leaving GetForcedGates() where TargetFlow= ', 
     &        TargetFlow, ' m3/s using SWP eq=', LStoSWP,
     &        ' gates=', SupGates, ' gateflow=', GateFlow, 
     &        ' SideChan=', SideChan, ' m3/s  SWPierLev=',
     &         SWPierLev, ' m NominalStMarys=',
     &         NominalStMarys, ' m3/s', LLT
            CALL LogIt(LOGD,0,0) 
         ENDIF

      END      
C
C
C----&------------------------------------------------------------------


      SUBROUTINE GetForcedSideChan(SupLevel, TargetFlow, MonthNum, 
     &   RoundDownOrUp, PreviousGates, SupGates, GateFlow,
     &   SWPierLev, SideChanForced, NominalStMarys)   
C
C        Purpose is to figure out the gates open and corresponding gated flow 
C        given a forced side channel flow condition (e.g., during hydropower outages)
C        The total river flow is also returned in the last argument, since
C        sometimes the target flow is impossible to meet. 
C
C        This routine assumes steady-state conditions and relies on 
C        GetSWPier() to solve simultaneous equations for calculating SW Pier
C        level from Lake Superior elevation.  The gated flow represents an
C        instantaneous solution rather than a mean value over some interval.  
C        If mean values over a period are needed, use SolveGatesSide(), which
C        includes the change-in-storage while solving the simultaneous 
C        equations (look above Chief Wiggins for example).
C
C        Uses the SWPier-Lake Superior relationship indicated by LS2SWP held 
C        in a common block.
C
C        Hardcoded to do no rounding.
C       
C

         IMPLICIT NONE 
C
C        Input args
         DOUBLE PRECISION SupLevel, TargetFlow 
c         DOUBLE PRECISION leakage, fishery, totalleak          
         INTEGER MonthNum, RoundDownOrUp

C        Output
         DOUBLE PRECISION GateFlow, SideChan, NominalStMarys, SWPierLev
         REAL SupGates, PreviousGates
C
C        local variables
         INTEGER I
         DOUBLE PRECISION  SideChanCapacity, SideChanForced
         DOUBLE PRECISION TmpFlow, PrevGateFlowDiff
C         DOUBLE PRECISION ROUND 
C                      Returns value rounded to a certain accuracy,
C                      according to the old Lake Survey method.
C
C        commons
         INTEGER SideChanMaxCapacity(12), SideChanMinCapacity(12)
         DOUBLE PRECISION CWOvertopElev, CWMinimumFreeboard(12)

C               These variables included just beacuse they are
C               in same COMMON as LStoSWP and SupCritRnd
         INTEGER ForecastMonths, LStoSWP, SupCritRnd, MonthsChecked
         INTEGER SPOPNBS06(12,6), SPOPNBS07(12,6), OPLLO(12), 
     &           MHOPNBS06(12,6), MHOPNBS07(12,6), OPChicago, 
     &           SCRRet(12), SCOPNBS06(12,6), SCOPNBS07(12,6), 
     &           DetRet(12), EROPNBS06(12,6), EROPNBS07(12,6),
     &           OPWelland(12), NiaRet(12),
     &           SCRRet_P2012(12), SCOPNBS08(12,6), SPOPNBS08(12,6), 
     &           DetRet_P2012(12), EROPNBS08(12,6), MHOPNBS08(12,6),
     &           OPWelland_P2012(12), NiaRet_P2012(12)       
         LOGICAL BackCompat
         LOGICAL LOGC(8), LOGD(8), LOGS(8), 
     &           LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &           LOGM(8), LOGP(8), LOGO(8), 
     &           LOGDM(8), LOGDP(8), LOGDO(8), 
     &           LOGSM(8), LOGSP(8), LOGSO(8),
     &           LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &           LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
         INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &           SupVerbosity, MLVerbosity, OntVerbosity
         CHARACTER   CTMP*255, CWRK*80, LLT*1
         CHARACTER*3 MonthAbbrev(12)
         DOUBLE PRECISION OldSWPA(12), OldSWPB(12), OldSWPC
         DOUBLE PRECISION EC2000SWPA, EC2000SWPB, EC2000SWPC(12)
         DOUBLE PRECISION EC2008SWPA, EC2008SWPB, EC2008SWPC(12)
C                      Constants for SW Pier/Pt Iroquois relationship
         INTEGER ForcedSideChan (12) 
         DOUBLE PRECISION SWPGateEqA(16), SWPGateEqB(16)
         COMMON /RCPARI/ GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
         COMMON /SUPCOM/ SPOPNBS06, OPLLO, MHOPNBS06, OPChicago, SCRRet,
     &                  SCOPNBS06, DetRet, EROPNBS06, OPWelland, NiaRet,
     &                  SPOPNBS07, MHOPNBS07, SCOPNBS07, EROPNBS07,
     &                  SPOPNBS08, MHOPNBS08, SCOPNBS08, EROPNBS08,
     &                  OPWelland_P2012, SCRRet_P2012, DetRet_P2012,
     &                  NiaRet_P2012,     
     &                  ForecastMonths,
     &                  LStoSWP, SupCritRnd, BackCompat, MonthsChecked
         COMMON /SPOTOP/ SideChanMaxCapacity, CWOvertopElev, 
     &                  CWMinimumFreeboard, SideChanMinCapacity        
         COMMON /MONAME/ MonthAbbrev
         COMMON /WKCHAR/ CTMP,LLT,CWRK
         COMMON /LOGCOD/ LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
         COMMON /SWPIER/SWPGateEqA, SWPGateEqB, OldSWPA, OldSWPB, 
     &                  OldSWPC, EC2000SWPA, EC2000SWPB, EC2000SWPC,
     &                  EC2008SWPA, EC2008SWPB, EC2008SWPC
         COMMON /FORCESIDE/ ForcedSideChan  

C
         IF ( SupVerbosity .GE. 5. ) THEN
            WRITE(CTMP,'(A,A3,2(A,F9.3),A,I2, A,I6, A,F9.3, A1)') 
     &        '  Entering GetForcedSideChan()', MonthAbbrev(MonthNum),
     &        ' with Sup= ', SupLevel,'m  TargetFlow= ', TargetFlow, 
     &        ' m3/s using SWP eq=', LStoSWP,'  RoundDownOrUp=', 
     &        RoundDownOrUp,'SideChanForced= ', SideChanForced, LLT
            CALL LogIt(LOGD,0,0) 
         ENDIF
C
C
         IF (RoundDownOrUp .EQ. -9999) RoundDownOrUp = 0
C
         CALL GetSideChanMaxCapacity( SupLevel, TargetFlow, MonthNum, 
     &       SideChanCapacity)
     
c         IF ( SupVerbosity .GE. 5. ) THEN
c             WRITE(CTMP,'(A, F9.3, A, I5,A1)') 
c     &        '  PreviousGates=', PreviousGates,
c     &        ' month=', MonthNum, LLT
c              CALL LogIt(LOGD,0,0)
c         END IF  
C        See if Target flow at least fills up the sidechannel
C        If not, solve for SW Pier using Target flow
C        Calculate SWP level according to the Beginning Sup level and Target flow

         IF (SideChanForced .GE. SideChanCapacity) THEN
               SideChanForced = SideChanCapacity
         END IF 
            IF (SideChanForced .LE. SideChanMinCapacity(MonthNum)) THEN
            SideChanForced = SideChanCapacity
         END IF 
         
         TmpFlow = SideChanForced
       
         IF ((PreviousGates .GE. 0.5) .and. (MonthNum .LE. 4) .and. 
     &       (SupLevel .LE. CWOvertopElev))  THEN  

            CALL GetSWPier(SupLevel, MonthNum, PreviousGates, 
     &         IDNINT(TmpFlow),-99,LStoSWP,SupVerbosity,-9999,
     &         'GetInstantGate()',SWPierLev, GateFlow)         
               SupGates = PreviousGates
            IF (TargetFlow .LT. GateFlow) THEN
C              Exceed the target flow, but even .5 gate is too much.  
C              Set sidechannel to 0, and return riverflow = gate flow
               SideChan = 0.0D0
               NominalStMarys = GateFlow
            ELSE                 
c               SideChan = TargetFlow - GateFlow
c               IF (SideChan .GE. SideChanCapacity) THEN
c                   SideChan = SideChanCapacity
c               END IF    
               NominalStMarys = SideChanForced + GateFlow
            END IF  
            GOTO 3373  
         END IF  
                  
  
         PrevGateFlowDiff = 9999
C
C        Now loop through gates until pass target flow.  Save results of 
C        previous pass in case flow exceeding the target is farther than
C        flow short of target.
         DO 3353, I = 0, 16
            IF (I .EQ. 0 ) SupGates = 0.5
            CALL GetSWPier(SupLevel, MonthNum, SupGates, 
     &         IDNINT(TmpFlow),-99,LStoSWP,SupVerbosity,-9999,
     &         'GetInstantGate()',SWPierLev,GateFlow)

            IF ( SupVerbosity .GE. 5. ) THEN
             WRITE(CTMP,'(A, F9.3, (A,F9.3),A1)') 
     &        ' ForcedSideChan= ', SideChanForced, 
     &        '  GateFlow=', GateFlow, LLT
              CALL LogIt(LOGD,0,0)
            END IF   
C           See if total flow exceeds target yet
            IF((GateFlow + SideChanForced ) .LT. 
     &            TargetFlow) THEN 
C              Not there yet.  Save difference and continue  
               PrevGateFlowDiff = TargetFlow -
     &           (GateFlow + SideChanForced)
               IF (I .EQ. 0 ) THEN
                  SupGates = 1.0
               ELSE
                  SupGates = SupGates + 1
               ENDIF
            ELSE
C              Total outflow exceeded target flow.  Done looping, but
C              need to see which case.
C              Check for weird case:
               IF( RoundDownOrUp .EQ. 1 .OR. 
     &           (RoundDownOrUp .EQ. 0 .AND. (Prevgateflowdiff .GT.
     &           (GateFlow + SideChanForced - TargetFlow)))
     &           ) THEN
                 SideChan = SideChanForced
                 NominalStMarys = GateFlow + SideChan
C                 IF(NominalStMarys .GT. TargetFlow .AND.
C     &             NINT(SupGates*2) .EQ.1)THEN
C                 
C                    SideChan = TargetFlow - GateFlow
C                    NominalStMarys = TargetFlow
C                 END IF
                 WRITE(CTMP,'(4(A,F9.3),A,A1)') 
     &               ' round up target flow is ', TargetFlow, 
     &               ' m3/s  gate outflow is ', 
     &               GateFlow, ' m3/s, so setting river flow to ',
     &               NominalStMarys, ' m3/s, ForcedSideChan',
     &               SideChanForced,' m3/s', LLT
                    CALL LogIt(LOGD,0,0) 

C
               ELSE
C              Flow was closer to target at previous gate setting, or
C              always round down anyway
                 SupGates = SupGates - 1
                 IF (NINT(SupGates - .1) .EQ. 0) SupGates = 0.5
                 GateFlow = TargetFlow - PrevGateFlowDiff - 
     &              SideChanForced
                 SideChan = SideChanForced
                 NominalStMarys = GateFlow + SideChan
                 WRITE(CTMP,'(4(A,F9.3),A,A1)') 
     &               ' round down target flow is ', TargetFlow, 
     &            ' m3/s GateFlow  is ', 
     &               GateFlow ,' m3/s,',
     &               NominalStMarys, ' m3/s, ForcedSideChan',
     &               SideChanForced,' m3/s', LLT
                 CALL LogIt(LOGD,0,0) 
               ENDIF
C              OK have our answers
               GOTO 3373
            ENDIF
 3353    CONTINUE    
C              
C        If fall through here, then must be max gates, but 
C        couldn't reach targetflow
         SupGates = 16
         SideChan = SideChanForced
         NominalStMarys = SideChan + GateFlow
         WRITE(CTMP,'(3(A,F9.3),A,A1)') 
     &      ' Unable to achieve target flow of ', TargetFlow, 
     &      ' m3/s because 16 gate outflow is ', 
     &      GateFlow, ' m3/s, so setting river flow to ',
     &      NominalStMarys, ' m3/s',LLT
         IF ( SupVerbosity .GE. 5. )CALL LogIt(LOGC,0,0) 
         IF ( SupVerbosity .GE. 3. )CALL LogIt(LOGD,0,0) 
C
 3373    CONTINUE
 
 
         IF ( SupVerbosity .GE. 5. ) THEN
            WRITE(CTMP,'(A,F9.3,A,I2,A,F6.2,3(A,F9.3), A,A1)') 
     &        '  Leaving GetForcedSideChan() where TargetFlow= ', 
     &        TargetFlow, ' m3/s using SWP eq=', LStoSWP,
     &        ' gates=', SupGates, ' gateflow=', GateFlow, 
     &        ' river flow=', NominalStMarys, 'ForcedSideChan',
     &        SideChanForced,' m3/s', LLT
            CALL LogIt(LOGD,0,0) 
         ENDIF

      END
C
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------      
      
      SUBROUTINE GetForcedGatesSideChan(SupLevel, TargetFlow, MonthNum, 
     &   SupGates, GateFlow, SWPierLev, SideChan, NominalStMarys)
     
C
C       Forced gates and Side Channel flow (i.e., forced total flow, and therefore not Plan 2012 flow)
C        Purpose is to force total flow in the river rather than 
C        Determine the gate setting or side channel flow as is done in Plan 2012
C        This subroutine over-rides Plan 2012 flow and can be used to determine impacts of deviations  


         IMPLICIT NONE 
C
C        Input args
         DOUBLE PRECISION SupLevel, TargetFlow 
c         DOUBLE PRECISION leakage, fishery, totalleak          
         INTEGER MonthNum

C        Output
         DOUBLE PRECISION GateFlow, SideChan, NominalStMarys, SWPierLev
         REAL SupGates
C
C        local variables
         INTEGER I
         DOUBLE PRECISION  SideChanCapacity
         DOUBLE PRECISION TmpFlow, PrevGateFlowDiff
C         DOUBLE PRECISION ROUND 
C                      Returns value rounded to a certain accuracy,
C                      according to the old Lake Survey method.
C
C        commons
         INTEGER SideChanMaxCapacity(12), SideChanMinCapacity(12)
         DOUBLE PRECISION CWOvertopElev, CWMinimumFreeboard(12)

C               These variables included just beacuse they are
C               in same COMMON as LStoSWP and SupCritRnd
         INTEGER ForecastMonths, LStoSWP, SupCritRnd, MonthsChecked
         INTEGER SPOPNBS06(12,6), SPOPNBS07(12,6), OPLLO(12), 
     &           MHOPNBS06(12,6), MHOPNBS07(12,6), OPChicago, 
     &           SCRRet(12), SCOPNBS06(12,6), SCOPNBS07(12,6), 
     &           DetRet(12), EROPNBS06(12,6), EROPNBS07(12,6),
     &           OPWelland(12), NiaRet(12),
     &           SCRRet_P2012(12), SCOPNBS08(12,6), SPOPNBS08(12,6), 
     &           DetRet_P2012(12), EROPNBS08(12,6), MHOPNBS08(12,6),
     &           OPWelland_P2012(12), NiaRet_P2012(12)       
         LOGICAL BackCompat
         LOGICAL LOGC(8), LOGD(8), LOGS(8), 
     &           LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &           LOGM(8), LOGP(8), LOGO(8), 
     &           LOGDM(8), LOGDP(8), LOGDO(8), 
     &           LOGSM(8), LOGSP(8), LOGSO(8),
     &           LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &           LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
         INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &           SupVerbosity, MLVerbosity, OntVerbosity
         CHARACTER   CTMP*255, CWRK*80, LLT*1
         CHARACTER*3 MonthAbbrev(12)
         DOUBLE PRECISION OldSWPA(12), OldSWPB(12), OldSWPC
         DOUBLE PRECISION EC2000SWPA, EC2000SWPB, EC2000SWPC(12)
         DOUBLE PRECISION EC2008SWPA, EC2008SWPB, EC2008SWPC(12)
C                      Constants for SW Pier/Pt Iroquois relationship
         DOUBLE PRECISION ForcedSupGates(12)
         INTEGER          ForcedSideChan(12)
         DOUBLE PRECISION SWPGateEqA(16), SWPGateEqB(16)
         COMMON /RCPARI/ GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
         COMMON /SUPCOM/ SPOPNBS06, OPLLO, MHOPNBS06, OPChicago, SCRRet,
     &                  SCOPNBS06, DetRet, EROPNBS06, OPWelland, NiaRet,
     &                  SPOPNBS07, MHOPNBS07, SCOPNBS07, EROPNBS07,
     &                  SPOPNBS08, MHOPNBS08, SCOPNBS08, EROPNBS08,
     &                  OPWelland_P2012, SCRRet_P2012, DetRet_P2012,
     &                  NiaRet_P2012,     
     &                  ForecastMonths,
     &                  LStoSWP, SupCritRnd, BackCompat, MonthsChecked
         COMMON /SPOTOP/ SideChanMaxCapacity, CWOvertopElev, 
     &                  CWMinimumFreeboard, SideChanMinCapacity        
         COMMON /MONAME/ MonthAbbrev
         COMMON /WKCHAR/ CTMP,LLT,CWRK
         COMMON /LOGCOD/ LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
         COMMON /SWPIER/SWPGateEqA, SWPGateEqB, OldSWPA, OldSWPB, 
     &                  OldSWPC, EC2000SWPA, EC2000SWPB, EC2000SWPC,
     &                  EC2008SWPA, EC2008SWPB, EC2008SWPC
         COMMON /FORGATE/ ForcedSupGates 
         COMMON /FORSIDE/ ForcedSideChan  
C
         IF ( SupVerbosity .GE. 5. ) THEN
            WRITE(CTMP,'(A,A3,2(A,F9.3),A,I2, A,F6.2,A1)') 
     &        '  Entering GetForcedGateSideChan()', (MonthNum),
     &        ' with Sup= ', SupLevel,'m  TargetFlow= ', TargetFlow, 
     &        ' m3/s using SWP eq=', LStoSWP,'  SupGates=', 
     &        SupGates, LLT
            CALL LogIt(LOGD,0,0) 
         ENDIF
C
         CALL GetSideChanMaxCapacity( SupLevel, TargetFlow, MonthNum, 
     &       SideChanCapacity)
     
         IF (SideChan .GE. SideChanCapacity) THEN
               SideChan = SideChanCapacity
         END IF   
         IF (SideChan .LE. SideChanMinCapacity(MonthNum)) THEN
            SideChan = SideChanCapacity
         END IF 
         TmpFlow = SideChan
 
         CALL GetSWPier(SupLevel, MonthNum, SupGates, 
     &         IDNINT(TmpFlow),-99,LStoSWP,SupVerbosity,-9999,
     &         'GetForcedGate()',SWPierLev, GateFlow)         
               
         NominalStMarys = SideChan + GateFlow
 
 
         IF ( SupVerbosity .GE. 5. ) THEN
            WRITE(CTMP,'(A,F9.3,A,I2,A,F6.2,3(A,F9.3),A,A1)') 
     &        '  Leaving GetForcedGatesSideChan() where TargetFlow= ', 
     &        TargetFlow, ' m3/s using SWP eq=', LStoSWP,
     &        ' gates=', SupGates, ' gateflow=', GateFlow, 
     &        ' river flow=', NominalStMarys, 
     &        ' m3/s ForcedSideChan', SideChan,
     &        ' m3/s', LLT
            CALL LogIt(LOGD,0,0) 
         ENDIF

      END      
C
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
      SUBROUTINE GetSideChanMaxCapacity(SupLevel, TargetFlow, MonthNum, 
     &   SideChanMax)

         IMPLICIT NONE 
C
C        Input args
         DOUBLE PRECISION SupLevel, TargetFlow 
         INTEGER MonthNum
C
C        Local variables
         DOUBLE PRECISION  ApproxSWP, SideChanMax,
     &          JacobSideChanCapacity
C        Output
         DOUBLE PRECISION Round
         LOGICAL LOGC(8), LOGD(8), LOGS(8), 
     &           LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &           LOGM(8), LOGP(8), LOGO(8), 
     &           LOGDM(8), LOGDP(8), LOGDO(8), 
     &           LOGSM(8), LOGSP(8), LOGSO(8),
     &           LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &           LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
         INTEGER GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &           SupVerbosity, MLVerbosity, OntVerbosity
         CHARACTER   CTMP*255, CWRK*80, LLT*1
         CHARACTER*3 MonthAbbrev(12)
         DOUBLE PRECISION OldSWPA(12), OldSWPB(12), OldSWPC
         DOUBLE PRECISION EC2000SWPA, EC2000SWPB, EC2000SWPC(12)
         DOUBLE PRECISION EC2008SWPA, EC2008SWPB, EC2008SWPC(12)
         INTEGER SideChanMaxCapacity(12), SideChanMinCapacity(12)
         DOUBLE PRECISION CWOvertopElev, CWMinimumFreeboard(12)
         DOUBLE PRECISION SWPGateEqA(16), SWPGateEqB(16)
C
C        local variables
         COMMON /SPOTOP/ SideChanMaxCapacity, CWOvertopElev, 
     &                  CWMinimumFreeboard, SideChanMinCapacity        
         COMMON /MONAME/ MonthAbbrev
         COMMON /RCPARI/ GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
         COMMON /WKCHAR/ CTMP,LLT,CWRK
         COMMON /LOGCOD/ LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
         COMMON /SWPIER/SWPGateEqA, SWPGateEqB, OldSWPA, OldSWPB, 
     &                  OldSWPC, EC2000SWPA, EC2000SWPB, EC2000SWPC,
     &                  EC2008SWPA, EC2008SWPB, EC2008SWPC

C        Start by assuming that we are using our normal Side Channel Maximums
         SideChanMax = SideChanMaxCapacity(MonthNum)
C        Next we will check to see if we should be using the low water 
C        capacities based on the work of Bruxer, et al. 
         IF ( TargetFlow .LE. 0.D0) TargetFlow = 0.1D0     
         ApproxSWP = Round( EC2008SWPA * SupLevel -              
     &         EC2008SWPB * LOG(TargetFlow/10.0D0) -
     &         EC2008SWPC(MonthNum), -99)
C        The following code checks the maximum physical capacity based on twice
C        the Brookfield power capacity, as computed by Jacob's equations.  
C        Implicit in this check is the assumption that flows (NOT allocations) 
C        not to be equal on both sides of the border, per Treaty requirements.
C        Note: Much gnashing of teeth has accompanied the development of these 
C        equations.  For more details, see especially the following emails: 
C        from Jacob Bruxer dated 6/23/2011; from Rob Caldwell dated 7/6/2011; 
C        and from John Allis dated 7/7/2011.  It may also be helpful to read 
C        the previous report from Detroit District emailed out by Travis Dahl 
C        on 3/11/11.
         IF(MonthNum .GE. 5 .AND. MonthNum .LE. 11) THEN
C            "Ice-Free" equation
             JacobSideChanCapacity = 2 * (98.732D0 * ApproxSWP**2
     &                       - 35414 * ApproxSWP + 3175583) + 20
         ELSE 
C            "Ice-Affected" equation
             JacobSideChanCapacity = 2 * (95.3305D0 * ApproxSWP**2
     &                       - 34203.4D0 * ApproxSWP + 3067884) + 20
         END IF  
         IF(JacobSideChanCapacity .LE. SideChanMax) THEN
            SideChanMax = JacobSideChanCapacity
         END IF         

         IF ( SupVerbosity .GE. 5 ) THEN
            WRITE(CTMP,'(A,F9.3,2(A,F6.2),2(A,F9.3),A,A1)') 
     &        '  Leaving GetSideChanMaxCapacity() where TargetFlow= ',
     &        TargetFlow, ' m3/s using Lake Sup=', SupLevel,
     &        ' m SWP= ', ApproxSWP, ' m JacobSideChanCapacity=',
     &        JacobSideChanCapacity,' m3/s SideChanMax=', SideChanMax, 
     &        ' m3/s', LLT
            CALL LogIt(LOGD,0,0) 
         ENDIF
      END
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C