C
C       
C  XX     TODO:
C XXconvert monthlength to an array, so that can have superior in uniform, 
C XX  middle in actual, etc.
C XXwrite a 'nextjob' file, with EOP levels, and starting date of next run
C 
C 
      SUBROUTINE RUNIT()
C      
C
      IMPLICIT NONE 
C
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER I, J, K, L, M
      CHARACTER*512 ReadBuf
      CHARACTER     Unspecified*12
      CHARACTER*254 TmpStr
      DOUBLE PRECISION deltaStor, tmpD
C                   These are just buffers used throughout program.
C                   Scratch variables
C
      LOGICAL       FinishedLoop(2) 
C                   Get set to .TRUE. only after routing past the interval
C                   containing the run's ending date.  Element 1 refers
C                   to Lake Superior loop, second element is middle lakes
C                   loop.  If a downstream loop ends later than an upstream
C                   loop, the last flow from the upstream loop gets repeated
C                   until the the D/S loop gets past the run's end date.
C
      LOGICAL       NewWeek(3), NewQMonth(3), NewMonth(3)
C                   Indicates that the new routing interval begins
C                   a new week or month.  Note that weeks begin on the
C                   starting day for the simulation.  If routing weekly 
C                   for more than one lake then the periods for both 
C                   lakes must being on the same day of the week.  (You
C                   can't do weekly routing for Lake Superior starting
C                   on monday, and Lake Ontario starting on thursday).
C                   ReadLine indicates that a new summary input line 
C                   should be read in.
C
      INTEGER       OntarioPeriod, SuperiorPeriod, MiddlePeriod
C                   A counter of the number of routing periods elapsed
C                   since the start of the simulation.
C
      INTEGER       TotalDays
C                   A counter of the total number of days that will 
C                   during the routing, to be used eventually in kind 
C                   of progress estimate. Note that runs of
C                   greater than 2147483646  (or so) will blow this
C                   value up.  Fix it if becomes a problem; shouldn't
C                   cause model to puke, but the period number in  
C                   printouts will look goofy after a point.
C
      INTEGER       CompletionStage(2)
C                   Holds the day number appropriate to the next 20%
C                   of progress through the main loop, and the percent.
C
      INTEGER       SecsInRtPeriod(3)
C                   The SecsInRtPeriod value contains the number of
C                   seconds in the current routing interval for each
C                   lake, used for calculating dT.  It gets set once
C                   for weekly, daily, or uniform month/qmonth 
C                   routing, but gets updated each period for actual
C                   month/qmonth routing.
C
      CHARACTER*8   OutputFileNames(5)
C                   These are the root parts of the file names for 
C                   summary output.  They are passed to OutFileName()
C                   with length determined by ENDSTR().  If the max
C                   length changes from 8 chars, must update the call to 
C                   ENDSTR(), and revise padding in DATA statement for
C                   OutputFileNames.
C
      INTEGER       OutFileHandles(5)
C                   Unit specifiers for summary output files, 
C                   corresponding to OutputFileNames(5).  Must be 
C                   coordinated with values in LogIt().
C
      INTEGER       FileNamesIndex(5)
C                   Indices to FileNames corresponding to summary input
C                   filenames. 
C
      INTEGER       NumRtLoops
C                   Number of routing loops.  These are nested DO loops
C                   organizing the order of solution through the Great
C                   Lakes.   Currently only three loops.  Innermost is
C                   Superior, with Michigan-Erie wrapped around it.
C                   Outermost loop is Lake Ontario.  If get into three-
C                   lake regulation or something, will have to modify
C                   this value.  (and mess around quite a bit with the
C                   nesting logic, but the implementations can't be
C                   anticipated).
C
      INTEGER       RtLoop2LkSlvI(3)
C                   Lists which member of LakesToSolve to use for each 
C                   routing loop. 
C                   element       Routing Loop    Lake to Solve     
C                     1          Lake Superior                1
C                     2   Lakes Michigan-Huron                2
C                     3           Lake Ontario                5
C
      CHARACTER*19  ReadFormat(4)
C                   Format specifiers for reading the times-series 
C                   summary files.  First element is for no oddball
C                   columns (MH & SC); second element is for one 
C                   oddball column (ER), third is for 2 oddball
C                   columns (superior), third is for 10 oddball
C                   columns (ontario/st lawrence).
C
      DOUBLE PRECISION MHBOP, MHEOP, SCBOP, SCEOP, ERBOP, EREOP, 
     &              SPBOP, SPEOP, ONBOP, ONEOP 
C                   Beginning- and end-of-period levels for each lake.
C                   Note that these storages are actually aliased to
C                   memory in FWRK.
C
      DOUBLE PRECISION MHBOI, MHEOI, SCBOI, SCEOI, ERBOI, EREOI,
     &              SPBOI, SPEOI, ONBOI, ONEOI
C                   Beginning- and end-of-period levels for each lake.
C                   Note that these storages are actually aliased to
C                   memory in FWRK.
C
      DOUBLE PRECISION NIZ
C                   Variable used for two-gage equation on Niagara 
C                   River.  NIZ is mean level for previous increment.  
C                   Note that this storage is actually aliased to
C                   memory in FWRK.
C                   
      DOUBLE PRECISION NiagaraRead
C                   Variable used for reading Niagara flows from input
C                   data instead of being computed.  
C                   Note that this storage is actually aliased to
C                   memory in FWRK.
C                   
      DOUBLE PRECISION FlowIntoWelland
C                   Holds the Welland Canal flow for the current 
C                   Lake Erie routing period (converted to positve 
C                   number). Note that this storage is actually 
C                   aliased to memory in FWRK.
C
      DOUBLE PRECISION FlowFromWelland
C                   Holds the Welland Canal flow for the current 
C                   Lake Ontario routing period (converted to positve 
C                   number).  Data is from same source as for the 
C                   diversion out of Lake Erie, although different
C                   routing periods on Lake Erie and Lake Ontario may
C                   result in different values being used.  Note that 
C                   this storage is actually aliased to memory in FWRK.
C       
      DOUBLE PRECISION LakeOntarioFlow, StLouisFlow
C                   Hold the St. Lawrence River flow and outflow from
C                   Lake St. Louis.  I don't know what the distinction
C                   is at this time.  Note that 
C                   this storage is actually aliased to memory in FWRK.
C
      DOUBLE PRECISION StLouisRetard
C                   Hold the amount of retardation in the St. Lawrence
C                   River at Lake St. Louis?  I don't know what this
C                   is for, but it's m3s.  Note that this storage is 
C                   actually aliased to memory in FWRK.
C
      DOUBLE PRECISION OntDeviation
C                   Amount Lake Ontario deviation from regulation plan
C                   in m3s.  Note that this storage is 
C                   actually aliased to memory in FWRK.
C
      DOUBLE PRECISION LSRoughness
C                   Hold the roughness factor at Long Sault for the
C                   current Lake Ontario regulation period.  
C                   Note that this storage is 
C                   actually aliased to memory in FWRK.
C
      DOUBLE PRECISION DesPrFlow, MilleFlow
C                   The flow from Des Prairies-Mille Iles.  
C                   Note that this storage is 
C                   actually aliased to memory in FWRK.
C
      DOUBLE PRECISION 
c These are the variables representing the downstream roughnesses and flows 
c at Montreal and beyond (Sorel, Varennes, Batiscan, 3 Rivers, etc...)
     &              Qdpmi, Qstfran, Qstmau, QRich,
     &              RJetty1, RVarennes, RSorel, RPierre,
     &              R3Rivers, RBatiscan, Tidal, Onttnbs,
     &              Other06n, Other06Q
C
      DOUBLE PRECISION OntForceFlow
C                   Lake Ontario forced outflow
C                   Note that this storage is 
C                   actually aliased to memory in FWRK.
C
      DOUBLE PRECISION LouisOntFlowDiff
C                   Lake Ontario time-series data
C                   Note that this storage is 
C                   actually aliased to memory in FWRK.
C
      DOUBLE PRECISION NumSupForceGates, NumOntForceGates
C                   Variables used to store the number of gates open 
C                   at the Compensating Works or Long Sault Dam River,
C                   if gate settings are forced.  Note that this 
C                   storage is actually aliased to memory in FWRK.
C
      DOUBLE PRECISION SupSideChannel
C                   The Lake Superior sidechannel flow to use.  If you
C                   specify this, but you forced a gate setting for a
C                   period, the model assumes 2320 m3s.
C
      DOUBLE PRECISION SupForceFlow
C                   The St. Marys River flow to use if solving for
C                   Lake Superior levels, but forcing an outflow.
C
      DOUBLE PRECISION SupDevFlow
C                   Deviations from Lake Superior regulated discharge 
C                   .
C
      DOUBLE PRECISION SupplyMH, SupplyON
      DOUBLE PRECISION MHNTS, SCNBS, ERNBS, SPNTS, ONNTS
      DOUBLE PRECISION MHICE, SCICE, ERICE, SMICE, ONICE
C                   Values read from summary input files for net supply
C                   and ice/weed retardation.  Note that these storages
C                   are actually aliased to memory in FWRK.
C
      DOUBLE PRECISION SPMLV, MHMLV, SCMLV, ERMLV, CGMLV, ONMLV
C                   Period-average levels for each lake.  
C                   Note that these storages are actually aliased to
C                   memory in FWRK.
C
C      LOGICAL       MLSolveError       
C                   Indicates whether the simultaneous solution of the
C                   routing equations checks out regarding continuity.
C
      DOUBLE PRECISION MLErrorFlow, SPErrorFlow
C                   Maximum error (in terms of river flows) seen during
C                   continuity checks. 
C
      INTEGER       NumMLSolveError, NumSPSolveError
C                   Number of times that middle lakes routing solution
C                   fails continuity check.  And Superior.
C  
      DOUBLE PRECISION STMRIVFLOW, STCRIVFLOW, DETRIVFLOW, NIARIVFLOW,
     &              NominalStMarys
C                   Period-average flows for each connecting channel.
C                   NominalStMarys is what the regulation plan calls
C                   for (see description in subroutine RouteSuperior).
C                   Note that these storages are actually aliased
C                   to memory in FWRK.
C                   
      DOUBLE PRECISION MHOutI, SCOutI, EROutI, SPOutI, ONOutI
C                   Outflows computed from the mean level of the current
C                   increment for each lake.
C                   Note that these storages are actually aliased
C                   to memory in FWRK.
C                   
      DOUBLE PRECISION SolutionTest, MHSolutionCheck, SCSolutionCheck, 
     &              ERSolutionCheck, SPSolutionCheck
C                   Variables used to verify that the simultaneous
C                   solution of the  routing equations was accurate.
C                   Note that these storages are actually aliased to
C                   memory in FWRK.
C
      DOUBLE PRECISION StMarysFlowAccum, NiagaraFlowAccum
C                   Accumulation storage used for converting flows
C                   in the connecting channels from the upstream routing
C                   interval to the downstream routing interval.
C
      INTEGER       StMarysDaysAccum, NiagaraDaysAccum
C                   Number of days accumulated for converting flows
C                   in the connecting channels from the upstream routing
C                   interval to the downstream routing interval.
C
      DOUBLE PRECISION StMarysConverted, NiagaraConverted
C                   Connecting channel flows converted into the downtream
C                   routing interval.
C
      DOUBLE PRECISION MHEndofMonthL, SCEndofMonthL, EREndofMonthL,
     &              SPEndofMonthL, SPMonthlyMeanF, ERMonthlyMeanF
C                   The end-of-month levels and mean monthly flow, 
C                   computed by TallyMeans, that get used as beginning
C                   of month levels and previous month's outflow by 
C                   Superior regulation plans.              
C                   Note that these storages are actually aliased to
C                   memory in FWRK.
C----&------------------------------------------------------------------
C
C     DOUBLE PRECISION MHWeeklyMeanF, MHWeeklyAccumF, MHQmonthlyMeanF,
C    & MHQmonthlyAccumF, MHMonthlyMeanF, MHMonthlyAccumF,  
C    & SCWeeklyMeanF, SCWeeklyAccumF, SCQmonthlyMeanF, SCQmonthlyAccumF, 
C    & SCMonthlyMeanF, SCMonthlyAccumF, ERWeeklyMeanF, ERWeeklyAccumF, 
C    & ERQmonthlyMeanF,  ERQmonthlyAccumF, ERMonthlyMeanF, 
C    & ERMonthlyAccumF, SPWeeklyMeanF, SPWeeklyAccumF, SPQmonthlyMeanF,
C    & SPQmonthlyAccumF, SPMonthlyMeanF, SPMonthlyAccumF, ONWeeklyMeanF,
C    & ONWeeklyAccumF, ONQmonthlyMeanF, ONQmonthlyAccumF,ONMonthlyMeanF,
C    & ONMonthlyAccumF, MHWeeklyMeanL, MHWeeklyAccumL, MHQmonthlyMeanL,
C    & MHQmonthlyAccumL, MHMonthlyMeanL, MHMonthlyAccumL,  
C    & SCWeeklyMeanL, SCWeeklyAccumL, SCQmonthlyMeanL, SCQmonthlyAccumL, 
C    & SCMonthlyMeanL, SCMonthlyAccumL, ERWeeklyMeanL, ERWeeklyAccumL, 
C    & ERQmonthlyMeanL,  ERQmonthlyAccumL, ERMonthlyMeanL, 
C    & ERMonthlyAccumL, SPWeeklyMeanL, SPWeeklyAccumL, SPQmonthlyMeanL,
C    & SPQmonthlyAccumL, SPMonthlyMeanL, SPMonthlyAccumL, ONWeeklyMeanL,
C    & ONWeeklyAccumL, ONQmonthlyMeanL, ONQmonthlyAccumL,ONMonthlyMeanL,
C    & ONMonthlyAccumL
C
      DOUBLE PRECISION SupForceEOPLev, MHForceEOPLev, STCForceEOPLev,
     &              ErieForceEOPLev, OntForceEOPLev
C                   EOP levels that replace the modeled EOP levels, if 
C                   not missing.  Note that these values have no impact
C                   on the routing period for which are specified - 
C                   they just get uses as the BOP level for the next
C                   routing period.  Had to implement it this way
C                   (instead of just forcing BOP level for current 
C                   period) so that Plan 1977A would be aware of forced 
C                   BOP levels on the middle lakes.
C                   
C
      CHARACTER     CHNKST*16
C                   This is the calling subroutine name passed to
C                   ErrMess() in case of error messages.
C
      INTEGER       GETQMONNUM 
C                   Determines what quarter-month a date is in
C
      INTEGER       LASTKC
C                   Function returning the position of the last useful
C                   character in a string.
C
C     INTEGER       ENDSTR
C                   Function returning end position of current word.
C 
      LOGICAL       LEAPYEAR
C                   Function returns .TRUE. if argument is a leap year.
C 
      CHARACTER*(*) SumTitleStr
      CHARACTER*(*) UnsolvedMessage
C                   Just strings used repeatedly.
C
C      DOUBLE PRECISION OUTFLOW
C                   Returns connecting channel flow, given the
C                   parameters for the discharge equation.
C
      DOUBLE PRECISION GETAREA
C                   Function that returns surface area for a lake
C
      DOUBLE PRECISION ROUND
C                   Returns value rounded to a certain accuracy,
C                   according to the old Lake Survey method.
C
      CHARACTER*256 OutFileName
C                   Function that pre-pends the output directory and
C                   appends the file extension to a root-name given as
C                   an argument, creating an output filename.
C
      CHARACTER*16  MakeIntervalName
C                   Returns a string with the interval name.  For output.
C 
      LOGICAL       InitTSOutputFile
C                   Initializes time-series output files. 
C
C      INTEGER       TSOutputFileName
C                   Formulates filename and path for output files.
C
C
C
C     Declarations regarding COMMON block variables:
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
C                   month within the month (so Q(4,1) tells month length
C                   for Lake Superior routing).
C
      INTEGER       DaysSinceSaturday
C                   The start date's offset from the previous Saturday 
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
      INTEGER       LkSlv2RtLoopI(8)
C                   Lists which routing loop to use for each member
C                   of LakesToSolve. 
C                   element      Lake to Solve     Routing Loop 
C                     1          Lake Superior                1
C                     2   Lakes Michigan-Huron                2
C                     3         Lake St. Clair                2
C                     4              Lake Erie                2
C                     5           Lake Ontario                3
C                     6      Lake St. Lawrence                3
C                     7               (unused)                -
C                     8               (unused)                -
C
      CHARACTER     SDATE*9,STIME*8
C                   See descriptions in INITPARAMS()
C
      CHARACTER*80  RCTITLE(3), MLHEADER
C                   See descriptions in PARLINE()
C
      LOGICAL       LakesToSolve(8) 
C                   Describes which lakes to solve for:
C                     1          Lake Superior
C                     2   Lakes Michigan-Huron
C                     3         Lake St. Clair
C                     4              Lake Erie
C                     5           Lake Ontario
C                     6      Lake St. Lawrence
C                     7               (unused)
C                     8               (unused)
C
      INTEGER       StartYear, StartMonth, StartDay, 
     &              EndYear, EndMonth, EndDay
C                   See descriptions in Block Data Genral
C
      INTEGER       MidLakeSolution
C                   Mid-Lakes solution method (1=iterative, 2=matrix).
C
      INTEGER       StMarysSolution, P77ACompat
C                   How to compute Superior outflows (1=outflow eqn, 
C                   3=plan77A), and what "operational" suuplies to use
C
      INTEGER       UStLawrSolution
C                   How to compute Ontario outflows (1=outflow eqn, 
C                   2=plan 1568D, 3= IS Plan, 4=plan1998)
C
      INTEGER       RoutingInterval(3), RoutingIncrement(3)
C                   Routing interval to use on Superior, middle lakes,
C                   and Ontario (1=daily, 2=weekly, 3=qmonthly, 
C                   4=monthly), and the number of increments within each 
C                   interval.
C
      INTEGER       BOMRound, FlowRound, LevelRound
C                   Control how End-Of-Month levels, and intermediate
C                   flows and  levels are rounded during the 
C                   computations. See descriptions in Block Data Genral.
C
      DOUBLE PRECISION LevelToInteger 
C                   Aliased to FWRK, used in other functions for solving
C                   simultaneous equations iteratively.
C
      CHARACTER*9   MLSolutionName(2)
C                   Strings set to 'interative' and 'matrix', used for 
C                   output.
C
      CHARACTER*16  SMSolutionName(4), USLSolutionName(5)
C                   See description in Block Data Genral
C
      CHARACTER*8   IntervalName(4)
C                   an 8 element char array containign names of the
C                   routing  intervals supported by the model, for use
C                   in output.
C
      INTEGER       IntervalNameLen(4)
C                   Number of letters in each element of IntervalName
C
      INTEGER       NOPDIR, NOPEXT, TSOutputFileNum(4,2,8)
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
      CHARACTER*3   MonthAbbrev(12)
      CHARACTER*256 FileNames(100)
      CHARACTER*256 RepeatLines(80)
      CHARACTER*21  LAKENAME(16), RiverName(16)
      INTEGER       NLN(16), LkToSlv2LN(8), NRN(16)
      INTEGER       TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &              MissingVal, SecsInUniformMon
      DOUBLE PRECISION MPFT, HALF, ONE, ZERO
      DOUBLE PRECISION SGX(17), SGY(17), SGZ(17)
      INTEGER       MaxSGSetting
C                   See descriptions in Block Data Genral
C     
      REAL          SupGatesOpen
C                   The number of Superior gates open in current routing 
C                   period.
C
      DOUBLE PRECISION SPK,SPYM,SPA,SPB,SPWT,SPC
      DOUBLE PRECISION MHK,MHYM,MHA,MHB,MHWT,MHC
      DOUBLE PRECISION SCK,SCYM,SCA,SCB,SCWT,SCC
      DOUBLE PRECISION ERK,ERYM,ERA,ERB,ERWT,ERC
C                   See descriptions in CHECKPAR and OUTFLOW
C
      DOUBLE PRECISION SPSTLV, MHSTLV, SCSTLV, ERSTLV, CGSTLV, ONSTLV
C                   The initial lake levels.
C
      DOUBLE PRECISION SPBOMLev, MHBOMLev, SCBOMLev, ERBOMLev
C                   The beginning-of-month lake levels used by Plan 1977A
C
      DOUBLE PRECISION SupIndex1, SupIndex2, SupIndex3
C                   time-series of arbitrary indices used for Superior 
C                   regulation
C
      INTEGER       PrevStMarysFlow
C                   Initial Lake Superior outflow (average for month 0)
C
      REAL          PrevSupGates 
C                   Initial gates open at Comp Works
C 
      CHARACTER     CTMP*255, CWRK*80
C                   These are just buffers used throughout program.
C
      CHARACTER LLT*1
C                   This is the character that marks the end of a line
C                   sent to the LogIt subroutine.
C 
      DOUBLE PRECISION FWRK
C                   This is a common block used for 64-bit calculations.
      DIMENSION FWRK(225,12)
C                   Dimensioned it this way to provide a 225 element
C                   array for each pool.
C
      CHARACTER*101 MonthOutputBuf(10,2,8)
C                   Character buffer holding a decade of mean monthly
C                   levels and flows for each of the lakes.  Data gets 
C                   saved in char buffer to reduce memory requirements,
C                   and written in 10-year chunks to each file to help
C                   performance.
C
      CHARACTER*103 QMonthOutputBuf(4,10,2,8)
C                   Character buffer holding a decade of mean qmonthly
C                   levels and flows for each of the lakes.  Data gets 
C                   saved in char buffer to reduce memory requirements,
C                   and written in 10-year chunks to each file to help
C                   performance.
C
C     CHARACTER*90  WeekOutputBuf(11,2,8)
      CHARACTER*19  WeekOutputBuf(52,2,8)
C                   Character buffer holding two years of mean weekly
C                   levels and flows for each of the lakes.  Data gets 
C                   saved in char buffer to reduce memory requirements,
C                   and written in 1-year chunks to each file to help
C                   performance.
C
      CHARACTER*74  DailyOutputBuf(48,2,8)
C                   Character buffer holding one year of mean daily
C                   levels and flows for each of the lakes.  Data gets 
C                   saved in char buffer to reduce memory requirements,
C                   and written in 1-year chunks to each file to help
C                   performance.
C
      LOGICAL       WriteTSOutput(4,2,8)
C                   For each lake (designated by third index - 8), for
C                   each type of time-series output (flow or ave level),
C                   for each interval (d, w, q, or m), .TRUE. if should
C                   write out the data to a time-series output file. 
C       
      LOGICAL       WriteErieOutflow
C                   See description in Block Data Genral
C       
      LOGICAL       OnlyOneOpen
C                   .TRUE. if the time-series output files are to be 
C                   opened prior to EACH write, and closed afterward.  
C                   .FALSE. if files are to be left open throughout 
C                   routing calculations.
C
      CHARACTER     SupPreProject*6, BalancingFlow*8, SupGateFlow*6,
     &              SupPlanFlow*6, SupBOMLevel*8, SWPierLevel*8,
C    &              SupNavFlow*6, SupPowerFlow*6, 
     &              SupSideChanFlow*6, SPRegCodes*40, P77AHeader*120,
     &              P2012Header*120
C                   Lake Superior regulation output strings
C
      CHARACTER*(80) OutputHeader
C                   Used to identify columns in ontario output file.
C
      INTEGER      P58DStartDay, P58DStartMonth, P58DStartYear
C                  Indicate the date of the first Saturday of the simulation
C
      INTEGER      Extrapolations
      DOUBLE PRECISION LargestExtrapolation
C                  Variables used to describe inadequate area-elev data
C
      INTEGER      numMaxSC, numMinSCWinter, numMinSCSummer, numJacSC
C                   Counters for number of times the Superior sidechannel flow
C                   was overridden by physical capacity constraints
C
      CHARACTER     DirSep*5,OSName(5)*8,MessDB*256
      CHARACTER*64  ModelerName
C
      DOUBLE PRECISION WBCheckSP, WBCheckMH, WBCheckSC, WBCheckER
C                   WB variables keep a running summation of the water balance
C                   terms for each timestep.  The values are checked against  
C                   change in storage at each timestep.
C
      INTEGER       USSlipEquation
      DOUBLE PRECISION  USSIceWeedAdj
C
C
C  Common block statements.
C
      COMMON /USSLIP/   USSIceWeedAdj, USSlipEquation
      COMMON /GENRLC/   OSName,DirSep,MessDB, ModelerName
      COMMON /EXTRAP/   LargestExtrapolation, Extrapolations 
      COMMON /FSAT/     P58DStartDay, P58DStartMonth, P58DStartYear
      COMMON /ONTOUT/   OutputHeader
      COMMON /SUPOUT/   SupPreProject, BalancingFlow, SupGateFlow,
     &                  SupPlanFlow, SupBOMLevel, SWPierLevel, 
     &                  SupSideChanFlow, SPRegCodes, P77AHeader,
     &                  P2012Header
      COMMON /OUTNMI/   NOPDIR, NOPEXT, TSOutputFileNum, OnlyOneOpen
      COMMON /CLDATE/   dT, Day, Q, PeriodEnd, NumIncsVary,
     &                  DaysSinceSaturday 
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /STIME/    STIME,SDATE
      COMMON /TITLES/   RCTITLE, MLHEADER
      COMMON /RCRUNI/   LakesToSolve, 
     &                  IntervalNameLen, MonthLength,
     &                  StartYear, StartMonth, StartDay, 
     &                  EndYear, EndMonth, EndDay
      COMMON /RCRUNC/   RCRUNTYP, IntervalName, MonthLengthName,
     &                  MLSolutionName, SMSolutionName,
     &                  USLSolutionName
      COMMON /STAGE/    MHK,SCK,ERK,MHYM,SCYM,ERYM,MHA,SCA,ERA,
     &                  MHB,SCB,ERB,MHWT,SCWT,ERWT,MHC,SCC,ERC
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /MISINT/   MidLakeSolution, BOMRound, FlowRound, 
     &                  LevelRound, RoutingInterval, RoutingIncrement,
     &                  StMarysSolution, UStLawrSolution, P77ACompat
      COMMON /INITLV/   SPSTLV, MHSTLV, SCSTLV, ERSTLV, CGSTLV, ONSTLV,
     &                  PrevStMarysFlow, PrevSupGates, 
     &                  SPBOMLev, MHBOMLev, SCBOMLev, ERBOMLev 
      COMMON /SPDAT/    SPK, SPYM, SPA, SPB, SPWT, SPC,
     &                  SGX, SGY, SGZ, MaxSGSetting,
     &                  SupGatesOpen
      COMMON /FNAMEC/   FileNames, RepeatLines
      COMMON /FSPACE/   FWRK
      COMMON /LNAMEC/   LAKENAME, RiverName
      COMMON /LNAMEI/   NLN, LkToSlv2LN, NRN, LkSlv2RtLoopI
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
      COMMON /MONAME/   MonthAbbrev
      COMMON /OUTTXT/   MonthOutputBuf, WeekOutputBuf
      COMMON /OUTTX2/   DailyOutputBuf
      COMMON /OUTTX3/   QMonthOutputBuf
      COMMON /OUTSTF/   WriteTSOutput, WriteErieOutflow
      COMMON /PHYCAP/   numMaxSC, numMinSCWinter, numMinSCSummer,
     &                  numJacSC
      COMMON /SUPIDX/   SupIndex1, SupIndex2, SupIndex3
C
C
C     Assign convenient aliases for summary inflow values, which
C     will be copied into FRWK(1,?) in this subroutine.  Remember that
C     the M-H summary inflow may already include St. Marys River
C     contribution, if supplied as a given tiem-series.  Also recognize 
C     that consumptive use and diversions, etc, are already combined 
C     into each of these terms, which is a bit of a change from GLERL's 
C     original approach for MIDLAKES.
      EQUIVALENCE (MHNTS, FWRK(1,1))
      EQUIVALENCE (SCNBS, FWRK(1,2))
      EQUIVALENCE (ERNBS, FWRK(1,3))
C     No NBS applicable to CGIP, so borrow FWRK(1,4) for SupplyMH
      EQUIVALENCE (SupplyMH, FWRK(1,4))
      EQUIVALENCE (SPNTS, FWRK(1,5))
      EQUIVALENCE (ONNTS, FWRK(1,6))
C       
C     Assign convenient aliases for summary ice/weed retardation values, 
C     which will be copied into FRWK(2,?) in this subroutine.
      EQUIVALENCE (FWRK(2,1), MHICE)
      EQUIVALENCE (FWRK(2,2), SCICE)
      EQUIVALENCE (FWRK(2,3), ERICE)
C     No applicable term for CGIP, so borrow FWRK(2,4) for SupplyON
      EQUIVALENCE (SupplyON, FWRK(2,4))
      EQUIVALENCE (FWRK(2,5), SMICE)
      EQUIVALENCE (FWRK(2,6), ONICE)
C
C     Assign convenient aliases for beginning of increment levels,
C     which will be copied into FRWK(3,?) in this subroutine.
      EQUIVALENCE (FWRK(3,1), MHBOI)
      EQUIVALENCE (FWRK(3,2), SCBOI)
      EQUIVALENCE (FWRK(3,3), ERBOI)
C     No applicable term for CGIP; use the storage for NIZ
      EQUIVALENCE (FWRK(3,4), NIZ)
      EQUIVALENCE (FWRK(3,5), SPBOI)
      EQUIVALENCE (FWRK(3,6), ONBOI)
C        
C     Assign convenient aliases for end of increment levels, which will 
C     be written to FRWK(4,?) at the end of the routing subroutine.
      EQUIVALENCE (FWRK(4,1), MHEOI)
      EQUIVALENCE (FWRK(4,2), SCEOI)
      EQUIVALENCE (FWRK(4,3), EREOI)
C     No applicable term for CGIP; storage used for LSRoughness
      EQUIVALENCE (FWRK(4,4), LSRoughness)
      EQUIVALENCE (FWRK(4,5), SPEOI)
      EQUIVALENCE (FWRK(4,6), ONEOI)
C
C     Assign convenient aliases for the period average levels, using 
C     storage actually allocated to the common workspace.  The levels
C     will be initialized in FRWK(5,?) and computed by this subroutine.
      EQUIVALENCE (MHMLV, FWRK(5,1))
      EQUIVALENCE (SCMLV, FWRK(5,2))
      EQUIVALENCE (ERMLV, FWRK(5,3))
C     CGMLV is actually input variable; gets read instead of computed.
      EQUIVALENCE (CGMLV, FWRK(5,4))
      EQUIVALENCE (SPMLV, FWRK(5,5))
      EQUIVALENCE (ONMLV, FWRK(5,6))
C
C     Assign convenient aliases for the period final outflow terms, 
C     using storage actually allocated to the common workspace.  The 
C     flows will be copied into FRWK(6,?) in this subroutine.
      EQUIVALENCE (STCRIVFLOW, FWRK(6,1))
      EQUIVALENCE (DETRIVFLOW, FWRK(6,2))
      EQUIVALENCE (NIARIVFLOW, FWRK(6,3))
C     No applicable term for CGIP; use storage in FWRK(6,4) for MLErrorFlow.
      EQUIVALENCE (MLErrorFlow, FWRK(6,4))
      EQUIVALENCE (STMRIVFLOW, FWRK(6,5))
      EQUIVALENCE (LakeOntarioFlow, FWRK(6,6))
C        
C     Assign convenient aliases for beginning of period levels,
C     which will be copied into FRWK(7,?) in this subroutine.
      EQUIVALENCE (FWRK(7,1), MHBOP)
      EQUIVALENCE (FWRK(7,2), SCBOP)
      EQUIVALENCE (FWRK(7,3), ERBOP)
C     This implementation assumes that the CGIP levels are constant 
C     throughout the routing interval, so that the beginning-of-
C     period, the end-of-period, and period average values are all
C     equivalent.  So no need to specify a BOP variable.
C     No applicable term for CGIP; use storage in FWRK(7,4) for SPErrorFlow.
      EQUIVALENCE (SPErrorFlow, FWRK(7,4))
      EQUIVALENCE (FWRK(7,5), SPBOP)
      EQUIVALENCE (FWRK(7,6), ONBOP)
C        
C     Assign convenient aliases for end of period levels, which will 
C     be written to FRWK(8,?) after call to a routing subroutine.
      EQUIVALENCE (FWRK(8,1), MHEOP)
      EQUIVALENCE (FWRK(8,2), SCEOP)
      EQUIVALENCE (FWRK(8,3), EREOP)
C     This implementation assumes that the CGIP levels are constant 
C     throughout the routing interval, so that the beginning-of-
C     period, the end-of-period, and period average values are all
C     equivalent.  So no need to specify a EOP variable, and
C     leaving storage in FWRK(8,4) for FlowIntoWelland.
      EQUIVALENCE (FWRK(8,4), FlowIntoWelland)
      EQUIVALENCE (FWRK(8,5), SPEOP)
      EQUIVALENCE (FWRK(8,6), ONEOP)
C
C     Assign convenient aliases for solution check values, which will 
C     be written to FRWK(9,?) after call to a routing subroutine.
      EQUIVALENCE (FWRK(9,1), MHSolutionCheck)
      EQUIVALENCE (FWRK(9,2), SCSolutionCheck)
      EQUIVALENCE (FWRK(9,3), ERSolutionCheck)
C     No applicable term for CGIP, so use the storage for SolutionTest
      EQUIVALENCE (FWRK(9,4), SolutionTest)
      EQUIVALENCE (FWRK(9,5), SPSolutionCheck)
C     No applicable term for Ontario routing, so storage is available
C      EQUIVALENCE (FWRK(9,6), stillopen)
C
C     FWRK(10,?), FWRK(11,?), FWRK(12,?), FWRK(13,?), and FWRK(14,?)
C     get overwritten each increment in Middle Lakes routing.  
C     Use them for intermediate values if desired.
C
      EQUIVALENCE (FWRK(11,4), SupForceFlow)
      EQUIVALENCE (FWRK(11,5), NominalStMarys)
      EQUIVALENCE (FWRK(11,6), SupDevFlow)
C     EQUIVALENCE (FWRK(12,1), SWPI)      used in RouteSuperior() 
C     EQUIVALENCE (FWRK(12,2), GateFlowI) used in RouteSuperior()
C      EQUIVALENCE (FWRK(12,5), ComputedNonGate) used in RouteSuperior(), must persist
      EQUIVALENCE (FWRK(12,6), SupSideChannel)
C
C     Assign convenient aliases for the intermediate outflow terms, using 
C     storage actually allocated to the common workspace.  These 
C     represent outflows per increment.
      EQUIVALENCE (FWRK(15,1), MHOutI)
      EQUIVALENCE (FWRK(15,2), SCOutI)
      EQUIVALENCE (FWRK(15,3), EROutI)
C     No applicable term for CGIP; use FWRK(15,4) for LevelToInteger
      EQUIVALENCE (FWRK(15,4), LevelToInteger) 
      EQUIVALENCE (FWRK(15,5), SPOutI)
      EQUIVALENCE (FWRK(15,6), ONOutI)
C
C     Assign convenient aliases for the outflow accumulation terms, using 
C     storage actually allocated to the common workspace.  
      EQUIVALENCE (FWRK(16,5), StMarysFlowAccum)
      EQUIVALENCE (FWRK(16,3), NiagaraFlowAccum)
C     No applicable term for St. Clair River; storage used by DPMIFlow
      EQUIVALENCE (FWRK(16,1), MilleFlow)
C     No applicable term for Detroit River, storage for NumSupForceGates
      EQUIVALENCE (FWRK(16,2), NumSupForceGates)
C     No applicable term for CGIP, so use the storage for NumOntForceGates
      EQUIVALENCE (FWRK(16,4), NumOntForceGates)
C     No applicable term for St Lawrence River; use for Onteviation
      EQUIVALENCE (FWRK(16,6), OntDeviation)
C 
C     Assign convenient aliases for the converted outflow terms, using 
C     storage actually allocated to the common workspace.  
      EQUIVALENCE (FWRK(17,3), NiagaraConverted)
      EQUIVALENCE (FWRK(17,5), StMarysConverted)
C     No applicable term for St. Clair River; storage used by StLouisFlow
      EQUIVALENCE (FWRK(17,1), StLouisFlow)
C     No applicable term for Detroit River; storage used by FlowFromWelland
      EQUIVALENCE (FWRK(17,2), FlowFromWelland)
C     No applicable term for CGIP; storage used by StLouisRetard
      EQUIVALENCE (FWRK(17,4), StLouisRetard)
C     No applicable term for St Lawrence River; storage used by DPMIFlow
      EQUIVALENCE (FWRK(17,6), DesPRFlow)
C
C     Assign aliases for the terms required to compute mean outflows
C     and levels, using storage actually allocated to the common 
C     workspace.  Also save space for end-of-month levels.
C     Actually these aliases aren't even used; the data managed via
C     references directly to FWRK indices.  But they are listed here are
C     as comments to show what goes where, and as a reminder that the 
C     storage is being used.
C     EQUIVALENCE (FWRK(18,1), MHWeeklyMeanF)
C     EQUIVALENCE (FWRK(19,1), MHWeeklyAccumF)
C     EQUIVALENCE (FWRK(20,1), MHQmonthlyMeanF)
C     EQUIVALENCE (FWRK(21,1), MHQmonthlyAccumF)
C     EQUIVALENCE (FWRK(22,1), MHMonthlyMeanF)
C     EQUIVALENCE (FWRK(23,1), MHMonthlyAccumF)
C     EQUIVALENCE (FWRK(18,2), SCWeeklyMeanF)
C     EQUIVALENCE (FWRK(19,2), SCWeeklyAccumF)
C     EQUIVALENCE (FWRK(20,2), SCQmonthlyMeanF)
C     EQUIVALENCE (FWRK(21,2), SCQmonthlyAccumF)
C     EQUIVALENCE (FWRK(22,2), SCMonthlyMeanF)
C     EQUIVALENCE (FWRK(23,2), SCMonthlyAccumF)
C     EQUIVALENCE (FWRK(18,3), ERWeeklyMeanF)
C     EQUIVALENCE (FWRK(19,3), ERWeeklyAccumF)
C     EQUIVALENCE (FWRK(20,3), ERQmonthlyMeanF)
C     EQUIVALENCE (FWRK(21,3), ERQmonthlyAccumF)
      EQUIVALENCE (FWRK(22,3), ERMonthlyMeanF)
C     EQUIVALENCE (FWRK(23,3), ERMonthlyAccumF)
C     No applicable term for CGIP; storage in FWRK(18-23,4) unused.
C     EQUIVALENCE (FWRK(18,5), SPWeeklyMeanF)
C     EQUIVALENCE (FWRK(19,5), SPWeeklyAccumF)
C     EQUIVALENCE (FWRK(20,5), SPQmonthlyMeanF)
C     EQUIVALENCE (FWRK(21,5), SPQmonthlyAccumF)
      EQUIVALENCE (FWRK(22,5), SPMonthlyMeanF)
C     EQUIVALENCE (FWRK(23,5), SPMonthlyAccumF)
C     EQUIVALENCE (FWRK(18,6), ONWeeklyMeanF)
C     EQUIVALENCE (FWRK(19,6), ONWeeklyAccumF)
C     EQUIVALENCE (FWRK(20,6), ONQmonthlyMeanF)
C     EQUIVALENCE (FWRK(21,6), ONQmonthlyAccumF)
C     EQUIVALENCE (FWRK(22,6), ONMonthlyMeanF)
C     EQUIVALENCE (FWRK(23,6), ONMonthlyAccumF)
C     EQUIVALENCE (FWRK(24,1), MHWeeklyMeanL)
C     EQUIVALENCE (FWRK(25,1), MHWeeklyAccumL)
C     EQUIVALENCE (FWRK(26,1), MHQmonthlyMeanL)
C     EQUIVALENCE (FWRK(27,1), MHQmonthlyAccumL)
C     EQUIVALENCE (FWRK(28,1), MHMonthlyMeanL)
C     EQUIVALENCE (FWRK(29,1), MHMonthlyAccumL)
C     EQUIVALENCE (FWRK(24,2), SCWeeklyMeanL)
C     EQUIVALENCE (FWRK(25,2), SCWeeklyAccumL)
C     EQUIVALENCE (FWRK(26,2), SCQmonthlyMeanL)
C     EQUIVALENCE (FWRK(27,2), SCQmonthlyAccumL)
C     EQUIVALENCE (FWRK(28,2), SCMonthlyMeanL)
C     EQUIVALENCE (FWRK(29,2), SCMonthlyAccumL)
C     EQUIVALENCE (FWRK(24,3), ERWeeklyMeanL)
C     EQUIVALENCE (FWRK(25,3), ERWeeklyAccumL)
C     EQUIVALENCE (FWRK(26,3), ERQmonthlyMeanL) 
C     EQUIVALENCE (FWRK(27,3), ERQmonthlyAccumL)
C     EQUIVALENCE (FWRK(28,3), ERMonthlyMeanL)
C     EQUIVALENCE (FWRK(29,3), ERMonthlyAccumL)
C     No applicable term for CGIP; storage in FWRK(24-29,4) unused.
C     EQUIVALENCE (FWRK(24,5), SPWeeklyMeanL)
C     EQUIVALENCE (FWRK(25,5), SPWeeklyAccumL)
C     EQUIVALENCE (FWRK(26,5), SPQmonthlyMeanL)
C     EQUIVALENCE (FWRK(27,5), SPQmonthlyAccumL)
C     EQUIVALENCE (FWRK(28,5), SPMonthlyMeanL)
C     EQUIVALENCE (FWRK(29,5), SPMonthlyAccumL)
C     EQUIVALENCE (FWRK(24,6), ONWeeklyMeanL)
C     EQUIVALENCE (FWRK(25,6), ONWeeklyAccumL)
C     EQUIVALENCE (FWRK(26,6), ONQmonthlyMeanL)
C     EQUIVALENCE (FWRK(27,6), ONQmonthlyAccumL)
C     EQUIVALENCE (FWRK(28,6), ONMonthlyMeanL)
C     EQUIVALENCE (FWRK(29,6), ONMonthlyAccumL)
      EQUIVALENCE (FWRK(30,1), MHEndofMonthL)
      EQUIVALENCE (FWRK(30,2), SCEndofMonthL)
      EQUIVALENCE (FWRK(30,3), EREndofMonthL)
C     No applicable term for CGIP; storage in FWRK(30,4) used by OntForceFlow
      EQUIVALENCE (FWRK(30,4), OntForceFlow)
      EQUIVALENCE (FWRK(30,5), SPEndofMonthL)
C     EQUIVALENCE (FWRK(30,6), ONEndofMonthL)
      EQUIVALENCE (FWRK(31,1), MHForceEOPLev)
      EQUIVALENCE (FWRK(31,2), STCForceEOPLev)
      EQUIVALENCE (FWRK(31,3), ErieForceEOPLev)
C     No applicable term for CGIP; storage in FWRK(31,4) used by NiagaraRead
      EQUIVALENCE (FWRK(31,4), NiagaraRead)
      EQUIVALENCE (FWRK(31,5), SupForceEOPLev)
      EQUIVALENCE (FWRK(31,6), OntForceEOPLev)
C
C     These storages are used internally by the matrix solution to the 
C     middle lakes routing equations.  They can be re-used as long as 
C     you don't mind them being trashed each time the middle lakes
C     route.
C          
C      EQUIVALENCE (FWRK(32,1), OUTMH)
C      EQUIVALENCE (OUT1, FWRK(33,1))
C      EQUIVALENCE (FWRK(34,1), OUT2)
C      EQUIVALENCE (FWRK(32,2), OUTSC)
C      EQUIVALENCE (FWRK(33,2), OUT3)
C      EQUIVALENCE (FWRK(34,2), OUT4)
C      EQUIVALENCE (FWRK(32,3), OUTER)
C      EQUIVALENCE (FWRK(33,3), OUT5)
C      EQUIVALENCE (FWRK(34,3), OUT6)
C     No applicable term for CGIP; storage in FWRK(32-4,4) unused.
C      EQUIVALENCE (FWRK(32,4), ???)
C      EQUIVALENCE (FWRK(33,4), ???)
C      EQUIVALENCE (FWRK(34,4), ???)
C
C
      PARAMETER (SumTitleStr = 'Summary Output File')
      PARAMETER (UnsolvedMessage = 
     &  'NO MEANINGFUL OUTPUT FOLLOWS - Not specified in LakesToSolve!')
C
C     Pad file names up to eight characters with trailing blanks: 
      DATA OutputFileNames  / 'summary ', 'superior', 'micherie',
     &   'ontario ', 'lawrence'/
C
C     1=summary, 2=superior, 3=micherie, 4=ontario, 5=lawrence
      DATA OutFileHandles / 92, 94, 93, 95, 96 /
C
C     DATA SecsInUniformMon /2629746/
C
      DATA NumRtLoops /3/
      DATA RtLoop2LkSlvI / 1, 2, 5/
      DATA FileNamesIndex /41, 31, 32, 33, 57/
C
      DATA StMarysDaysAccum, NiagaraDaysAccum / 0, 0/
C
C      DATA MLSolveError / .FALSE. /
      DATA NumMLSolveError, NumSPSolveError / 2 * 0 /
C
      DATA OntarioPeriod, SuperiorPeriod, MiddlePeriod / 3 * 0 /
C
      DATA CompletionStage / 1, 2 / 
C
C     Initialize these to TRUE, so ready in first iteration of loop 
C     Don't worry about NewWeek, it's reset at start of each iteration.
      DATA NewQMonth, NewMonth / 3*.TRUE., 3*.TRUE. /
C
C     XX change form array into single string, replace 12th char with '2', etc
      DATA ReadFormat / '(22X,2F12.3,1F9.3)', '(22X,2F12.3,3F9.3)', 
     &     '(22X,2F12.3,8F9.3)', '(22X,2F12.3,34F9.3)'/ 
C
      DATA FinishedLoop / 2*.FALSE. /
C
      DATA WBCheckSP,WBCheckMH, WBCheckSC, WBCheckER / 4*0.D0 /
C
      CHNKST='RunIt'//LLT
      Unspecified='Unspecified'//LLT
C
C     Clear any trash in workspace
      DO 3, J=1,8
         DO 3, I=1,125
    3       FWRK(I, J) = ZERO

C
C     The following initializations are unnecessary because the storage
C     in FWRK has been zero'd.  But if the EQUIVALENCE approach goes the
C     way of the Buffalo, they'll need to be uncommented.
C
C     Initialize Lake Superior outflow to 0.  It gets added to SupplyMH to 
C     update MHNTS.  If Superior not being solved, doesn't change anything.
C     STMRIVFLOW = ZERO
C     SupplyMH = ZERO
C
C     StMarysFlowAccum = ZERO
C     StMarysConverted = ZERO
C
C     NiagaraFlowAccum = ZERO
C     NiagaraConverted = ZERO
C
C     Initialize MLErrorFlow
C     MLErrorFlow = ZERO
C     SPErrorFlow = ZERO
C
C    
C      CHARACTER*74  DailyOutputBuf(48,2,8)
C     Initialize daily output buffers
      DO 221, L=1,8
         DO 221, J=1,2
            DO 221, I=1,48
               DailyOutputBuf(I,J,L)(1:1) = '#' 
               DO 221, K=2,74
  221             DailyOutputBuf(I,J,L)(K:K) = ' ' 
C
C     Initialize quarter-monthly output buffers
      DO 223, L=1,8
         DO 223, J=1,2
            DO 223, M=1,4
               DO 223, I=1,10
                  QMonthOutputBuf(M,I,J,L)(1:1) = '#' 
                  DO 223, K=2,102
  223                QMonthOutputBuf(M,I,J,L)(K:K) = ' ' 
C
C    
C     Initialize monthly output buffers
      DO 224, L=1,8
         DO 224, J=1,2
            DO 224, I=1,10
               MonthOutputBuf(I,J,L)(1:1) = '#' 
               DO 224, K=2,100
  224             MonthOutputBuf(I,J,L)(K:K) = ' ' 
C
C----&------------------------------------------------------------------
C
C     Initialize summary output files in output directory:
C
      DO 2003, J=1,5
C
C
C        See if LOGCOD variables need to use the file.  The test is
C        very kludgy and relies on the the element number (out of 8) in 
C        the LOGCOD variable being different from the I/O unit numbers 
C        by 88.  Ugh.  Any re-arranging of the LOGCOD elements or 
C        changing of the unit numbers will break this scheme.  Rather
C        than cause cryptic errors when the indices bust, trap it here:
         IF( OutFileHandles(J) - 88 .LT. 4 .OR. 
     &      OutFileHandles(J) - 88 .GT. 8) STOP 'Milhouse'
C
C        Use LOGSMPO as test variable, since it requires all the files
         IF( LOGSMPO(OutFileHandles(J) - 88) ) THEN
C
C           Need to open the file.  Will be closed as part of normal
C           program termination.
C
C           Build file name, store it in TmpStr for now.  
            TmpStr = OutFileName(OutputFileNames(J), 
     &         LastKC(OutputFileNames(J),8), K ) 
            OPEN(UNIT=OutFileHandles(J), FILE=TmpStr(1:K), IOSTAT=I) 
            IF ( I .NE. 0 ) THEN
               CALL ErrMess(1097,CHNKST,LOGCDE)
               WRITE(CTMP,9015) I, TmpStr(1:K), LLT
               GOTO 39
            ENDIF      
C        ELSE
C          Nothing needs the file; don't open it
         ENDIF      
C     1=summary, 2=superior, 3=micherie, 4=ontario, 5=lawrence
C     DATA OutFileHandles / 92, 94, 93, 95, 96 /
C        Turn off all elements for Superior output 
C        LOGP(6) = .FALSE.
C        Turn off all elements for middle lakes output
C        LOGM(5) = .FALSE.
C        Turn off all elements for Ontario output 
C        LOGO(7) = .FALSE.
C                  col 1 - stdout (unit 6)
C                  col 2 - detailed log file (unit 90)
C                  col 3 - critical failure log (unit 91)
C                  col 4 - summary output file (unit 92)
C                  col 5 - Middle Lakes output file (unit 93)
C                  col 6 - Lake Superior output file (unit 94)
C                  col 7 - Lake Ontario output file (unit 95)
C                  col 8 - St. Lawrence output file (unit 96)
C
 2003    CONTINUE
C
C       
C       
C----&------------------------------------------------------------------
C
C
c     CALL InputSummary()
C
C
C     Now write titles to summary output files:
C
 4002 FORMAT (A,A1)
      WRITE(CTMP,4002) SumTitleStr, LLT
      CALL LogIt(LOGS,1,0)
      WRITE(CTMP,4002) 'Someday will write an executive summary of '//
     &   'the results for this run HERE.', LLT
      CALL LogIt(LOGS,0,1)
C     TODO: Probably make a subroutine to write attractive
C           header block for this file, re-cap the input
C           parameters, and generally describe processing
C           in friendly terms.
C
      WRITE(CTMP,4002)SumTitleStr//' for Lake Superior',LLT
      CALL LogIt(LOGP,1,1)
      IF( .NOT. LakesToSolve(1) )THEN
         CTMP = UnsolvedMessage//LLT      
         CALL LogIt(LOGP,1,1)
      ENDIF
C
      WRITE(CTMP,4002)
     &   SumTitleStr//' for Lakes Michigan-Erie',LLT
      CALL LogIt(LOGM,1,1)
      IF( .NOT. LakesToSolve(2) )THEN
         CTMP = UnsolvedMessage//LLT      
         CALL LogIt(LOGM,1,1)
      ENDIF
C
      WRITE(CTMP,4002)
     &   SumTitleStr//' for Lake Ontario',LLT
      CALL LogIt(LOGO,1,1)
      IF( .NOT. LakesToSolve(5) )THEN
         CTMP = UnsolvedMessage//LLT      
         CALL LogIt(LOGO,1,1)
      ENDIF
C
C
C     Write the title lines
C
      DO 733, K=1,3
         WRITE(CTMP,4002)RCTITLE(K)(1:LASTKC(RCTITLE(K),80)),LLT
  733    CALL LogIt(LOGDMPO,0,0)
C
 4003 FORMAT ('Increments per ', A,1X,A,' routing interval: ',A,A1)
 4007 FORMAT('Summary of ', A,' Routing',A1) 
      IF( LakesToSolve(1) ) THEN
         CWRK = MakeIntervalName( RoutingInterval(1), .TRUE., J)
         WRITE(CTMP, 4007)CWRK(1:J), LLT
         CALL LogIt(LOGDP,1,1)
         CWRK = MakeIntervalName( RoutingInterval(1), .FALSE., J)
         WRITE(TmpStr,'(I6)') RoutingIncrement(1)
         IF (RoutingIncrement(1) .EQ. 0) TmpStr = 'hourly'
         WRITE(CTMP,4003)LakeName(1)(1:NLN(1)),CWRK(1:J),TmpStr(1:6),LLT
         CALL LogIt(LOGDP,0,1)
      ENDIF
C
C
      IF( LakesToSolve(2) ) THEN
         WRITE(CTMP,4002) MLHEADER(1:LASTKC(MLHEADER,80)),LLT
         CALL LogIt(LOGDM,1,0)
         CWRK = MakeIntervalName( RoutingInterval(2), .TRUE., J)
         WRITE(CTMP, 4007)CWRK(1:J), LLT
         CALL LogIt(LOGDM,0,1)
         CWRK = MakeIntervalName( RoutingInterval(2), .FALSE., J)
         WRITE(TmpStr,'(I6)') RoutingIncrement(2)
         IF (RoutingIncrement(2) .EQ. 0) TmpStr = 'hourly'
         WRITE(CTMP,4003)LakeName(2)(1:NLN(2)),CWRK(1:J),TmpStr(1:6),LLT
         CALL LogIt(LOGDM,0,1)
         WRITE(CTMP,4002) 'Stage-discharge equation parameters:',LLT
         CALL LogIt(LOGDM,0,0)
         WRITE(CTMP,4500) LLT
         CALL LogIt(LOGDM,0,0)
         WRITE(CTMP,4004)MHK,MHYM,MHA,MHB,MHWT,MHC,LLT
         CALL LogIt(LOGDM,0,0)
         WRITE(CTMP,4005) SCK,SCYM,SCA,SCB,SCWT,SCC,LLT
         CALL LogIt(LOGDM,0,0)
         WRITE(CTMP,4006) ERK,ERYM,ERA,ERB,ERWT,ERC,LLT
         CALL LogIt(LOGDM,0,1)
 4004 FORMAT ('Michigan-Huron: ',
     &   F9.4,4X,F6.2,4X,F3.1,4X,F3.1,4X,F3.1,2X,F5.2,A1)
 4005 FORMAT ('St. Clair: ',
     &   5X,F9.4,4X,F6.2,4X,F3.1,4X,F3.1,4X,F3.1,2X,F5.2,A1)
 4006 FORMAT ('Erie: ',10X,F9.4,4X,F6.2,4X,F3.1,4X,F3.1,4X,F3.1,
     &   2X,F5.2,A1)
 4500 FORMAT(20X,'-K-',7X,'-YM-',5X,'-A-',4X,'-B-',3X,'-WT-',3X,
     &   '-ZC-',A1)
      ENDIF
C
      IF( LakesToSolve(5) ) THEN
         CWRK = MakeIntervalName( RoutingInterval(3), .TRUE., J)
         WRITE(CTMP, 4007)CWRK(1:J), LLT
         CALL LogIt(LOGDO,1,1)
         CWRK = MakeIntervalName( RoutingInterval(3), .FALSE., J)
         WRITE(CTMP,4003) LakeName(7)(1:NLN(7)), CWRK(1:J), 
     &      RoutingIncrement(3), LLT
         CALL LogIt(LOGDO,1,1)
      ENDIF
C
C 
C     Write line with modeler's name, date of run, and model version
C
      WRITE( CTMP, 4068)
     &   StartDay, MonthAbbrev(StartMonth), StartYear,
     &   EndDay, MonthAbbrev(EndMonth), EndYear, LLT
      CALL LogIt(LOGDMPO,0,1)
      WRITE(CTMP,4066)ModelerName(1:LASTKC(ModelerName,64)),
     &   STIME, SDATE, LLT
      CALL LogIt(LOGDMPO,0,0)
      CALL BuildInfo
      CALL LogIt(LOGDMPO,0,0)
 4066 FORMAT('Modeled by ', A, ' at ', A, ' on ', A, A1)
 4068 FORMAT('Data generated from simulation starting on ', 
     &   I2, A3, I5, ', and ending on ', I2, A3, I5, A1)
C 
C     No more titles, headers, etc.  
C
C       
C----&------------------------------------------------------------------
C
C       
C
C     Open inflow summary files:
      DO 33, K=1,5
C        Write(*,'(A,I3)')'idx=',INDEX(FileNames(FileNamesIndex(K)),LLT)
C        Write(*,'(A,A )')'   =',FileNames(FileNamesIndex(K))
C    &      (1:INDEX(FileNames(FileNamesIndex(K)),LLT) - 1)
   33    IF (LakesToSolve(K)) OPEN(UNIT=K + 80, 
     &      FILE=FileNames(FileNamesIndex(K))
     &      (1:INDEX(FileNames(FileNamesIndex(K)),LLT) - 1),
     &      IOSTAT=I,ERR=37) 
      GOTO 80
C        
   37 CALL ErrMess(1003,CHNKST,LOGCDE)
 2015 FORMAT(' Got IOSTAT=',I4,' when opening ', A,
     &   ' summary inflow/outflow file.',A1)
      WRITE(CTMP,2015) I, 
     &   LakeName(LkToSlv2LN(K))(1:NLN(LkToSlv2LN(K))),LLT
   38 CALL LogIt(LOGCDE,0,0)
      WRITE(CTMP,'(A,A,A1,A1)')' The filename is {',
     &   FileNames(FileNamesIndex(K))
     &   (1:INDEX(FileNames(FileNamesIndex(K)),LLT) - 1), '}',LLT
   39 CALL LogIt(LOGCDE,0,0)
      CALL BuildInfo
      CALL LogIt(LOGCDE,0,0)
      CALL SAVETOD
      CALL ADIOS()
C
   77 CALL ErrMess(1086,CHNKST,LOGCDE)
 2016 FORMAT(' Got IOSTAT=',I4,' when reading ', A,
     &   ' summary inflow/outflow file.',A1)
      WRITE(CTMP,2016) I, 
     &   LakeName(LkToSlv2LN(K))(1:NLN(LkToSlv2LN(K))),LLT
      GOTO 38
C----&------------------------------------------------------------------
C
   79 CALL ErrMess(1087,CHNKST,LOGCDE)
      WRITE(CTMP,'(A,I5,A,A,A1)') ' Got IOSTAT=', I,
     &   ' when reading following line as FORMAT ',
     &   ReadFormat(J), LLT
      CALL LogIt(LOGCDE,0,0)
      WRITE(CTMP,8653) ReadBuf(1:80)//LLT
      CALL LogIt(LOGCDE,0,0)
      WRITE(CTMP,8653) ReadBuf(81:160)//LLT
      CALL LogIt(LOGCDE,0,0)
      WRITE(CTMP,8653) ReadBuf(161:240)//LLT
      CALL LogIt(LOGCDE,0,0)
      WRITE(CTMP,8653) ReadBuf(241:320)//LLT
      CALL LogIt(LOGCDE,0,0)
      WRITE(CTMP,8653) ReadBuf(321:384)//LLT
      GOTO 38
C       
C 9455 CALL ErrMess(1126,CHNKST,LOGCDE)
C      WRITE(CTMP,2017) I, CWRK(1:K), LLT
C 2017 FORMAT(' Got IOSTAT=',I4,' when accessing file ',A,'.',A1)
C      GOTO 39
C
C       
C----&------------------------------------------------------------------
C
C     Do some data initializations
C
C     Initialize EOP levels - they will be rolled over to BOP levels
C     inside loop through each routing period
   80 SPEOP = SPSTLV
      MHEOP = MHSTLV
      SCEOP = SCSTLV
      EREOP = ERSTLV
      ONEOP = ONSTLV
C
C     Plan 1977A parameter initialization
C     SPMonthlyMeanF = DBLE(PrevStMarysFlow)
      NominalStMarys = DBLE(PrevStMarysFlow)
      SupGatesOpen = PrevSupGates
      SPEndofMonthL = SPBOMLev
      MHEndofMonthL = MHBOMLev
      SCEndofMonthL = SCBOMLev  
      EREndofMonthL = ERBOMLev
C
C     In case using two-gage eqn, initialize NIZ as level from previous 
C     routing interval, which is start level:
      NIZ = CGSTLV
      IF( NIZ .LT. ZERO) NIZ = ZERO
C
C     Set the tolerance for flagging failures to satisfy continuity
C     to the rounding point for intermediate flow calculations.
C     For example, if FlowRound==1, then message will pop up saying
C     simultaneous solution of equations was unsatisfactory for at
C     least one increment if the check on water balance terms comes
C     out more than 10 m3s different.  If FlowRound <= -2, use .1 m3s
C     
C     XX make SolutionTest user-defined (ie, avoid warning 
C     XX message if below 1tcfs, etc)
      IF( FlowRound .LE. -2) THEN
         SolutionTest = ONE / TEN
      ELSE
         SolutionTest = 10.0D0**FlowRound
      ENDIF
C    
      IF( LevelRound .NE. -99) THEN
         LevelToInteger = ONE/10.0D0**LevelRound
C         write(*,'(A,I3,A,F20.9)') ' LevelRound=',
C     &      LevelRound, ' LevelToInteger=', LevelToInteger
      ENDIF
C    
C----&------------------------------------------------------------------
C
C     Inner loops not incremented unless that lake is being solved.  
C     But the most downstream loop (Lake Ontario) is main control loop.
C     It always gets incremented, and computations finish when it hits
C     end date.  So if not solving for Lake Ontario, sync it in with
C     next closest loop that is being solved, so that end date is
C     reached for main control loop at same time as at the most
C     downstream lake that is being solved.
C
      IF( LakesToSolve(5) )THEN
         I = 5
      ELSE
         DO 81, I = 4, 1, -1
   81       IF( LakesToSolve(I) ) GOTO 181
         WRITE(CTMP, '(A,A1)') ' No lakes being solved!', LLT
         CALL LogIt(LOGCDE,1,0)
         STOP 'Lisa'         
  181    RoutingInterval(NumRtLoops) = RoutingInterval(LkSlv2RtLoopI(I))
      ENDIF
C      WRITE(*,'(3(A,I5))') ' I=', I, ' LkSlv2RtLoopI(I)=', 
C     &   LkSlv2RtLoopI(I), ' RoutingInterval(LkSlv2RtLoopI(I))=',
C     &   RoutingInterval(LkSlv2RtLoopI(I))
      CWRK = MakeIntervalName( RoutingInterval(NumRtLoops), .FALSE., J)
      WRITE(CTMP, '(A,A,A,A,A,A1)')
     &   ' Syncing main control loop date counters with ',
     &   LakeName(LkToSlv2LN(I))(1:NLN(LkToSlv2LN(I))), ' {', CWRK(1:J), 
     &   '}',LLT
C----&------------------------------------------------------------------
      CALL LogIt(LOGCD,1,1)
C
C
      DO 82, L = 1, NumRtLoops
C 
Cz         This line moved down a ways so that PeriodEnd registers always
Cz         get initialized.
Cz         IF( LakesToSolve(RtLoop2LkSlvI(L)) .OR. L .EQ. NumRtLoops )THEN
C           Routing through this lake, so must set PeriodEnd values
C           based on starting date.
C
C           Don't need to set Day(L)- gets reset based on PeriodEnd(1,L)
            PeriodEnd(4,L) = StartYear
            PeriodEnd(3,L) = StartMonth
C           Store the ending days of StartMonth's qmonths in Q
            CALL GetQtrMonBounds (StartYear, StartMonth, Q(1,L))
C           Initialize PeriodEnd(2,L) to the starting 1/4 month
            PeriodEnd(2,L)= GetQmonNum (StartYear, StartMonth, StartDay)
            PeriodEnd(1,L) = StartDay - 1
            PeriodEnd(5,L) = 0
C
         IF( LakesToSolve(RtLoop2LkSlvI(L)) .OR. L .EQ. NumRtLoops )THEN
C           Set the duration (in seconds) of the initial interval
C           for each lake in SecsInRtPeriod, and the duration of the
C           increments within .  
            IF( RoutingInterval(L) .EQ. 1 ) THEN
C              Daily routing
               SecsInRtPeriod(L) = SecsInDay
               dT(L) = DBLE(SecsInRtPeriod(L))/DBLE(RoutingIncrement(L))
C
            ELSE IF( RoutingInterval(L) .EQ. 2 ) THEN
C              Weekly routing
               SecsInRtPeriod(L) = 7 * SecsInDay
               dT(L) = DBLE(SecsInRtPeriod(L))/DBLE(RoutingIncrement(L))
C
            ELSE IF( RoutingInterval(L) .EQ. 3 ) THEN
C              Quarter-Monthly routing - see if actual length or uniform
               IF(MonthLength .EQ. 1) THEN
C                 Uniform quarter-month length
C                 Number of increments varies, but dT stays at 1 hour
                  IF( NumIncsVary(L) ) RoutingIncrement(L) = 183
C                 else number of increments stays same, but dT changes
                  SecsInRtPeriod(L) = SecsInUniformMon / 4
                  dT(L)= DBLE(SecsInUniformMon) /
     &               DBLE(4*RoutingIncrement(L))
               ELSE
C                 Actual quarter-month length
                  IF(PeriodEnd(2,L) .EQ. 1)THEN
                     SecsInRtPeriod(L) = Q(1,L) * SecsInDay
                  ELSE
                     SecsInRtPeriod(L) = (Q(PeriodEnd(2,L),L) -
     &                  Q(PeriodEnd(2,L) - 1,L)) * SecsInDay
                  ENDIF
C                 SecsInRtPeriod is multiple of 86400, so 3600 will go into it
                  IF( NumIncsVary(L) ) RoutingIncrement(L) = 
     &               SecsInRtPeriod(L) / 3600
                  dT(L) = DBLE(SecsInRtPeriod(L)) /
     &               DBLE(RoutingIncrement(L))
               ENDIF
               PeriodEnd(2,L) = PeriodEnd(2,L) - 1
C
            ELSE IF( RoutingInterval(L) .EQ. 4 ) THEN
C              Monthly routing - see if actual length or uniform
               IF(MonthLength .EQ. 1)  THEN
C                 Uniform month length
                  SecsInRtPeriod(L) = SecsInUniformMon
                  IF( NumIncsVary(L) ) RoutingIncrement(L) = 732
                  dT(L) = DBLE(SecsInRtPeriod(L)) /
     &               DBLE(RoutingIncrement(L))
               ELSE
C                 Actual Month length
                  SecsInRtPeriod(L) = Q(4,L) * SecsInDay
                  IF( NumIncsVary(L) ) RoutingIncrement(L) = 
     &               SecsInRtPeriod(L) / 3600
                  dT(L) = DBLE(SecsInRtPeriod(L)) /
     &               DBLE(RoutingIncrement(L))
               ENDIF
               PeriodEnd(2,L) = 0
               PeriodEnd(3,L) = PeriodEnd(3,L) - 1
C               PeriodEnd(6,L) = ?
            ELSE
               WRITE(CTMP, '(A,I3,A,I1,A1)') ' Bad routing interval {',
     &            RoutingInterval(L), '} for loop ',L,LLT
               CALL LogIt(LOGCDE,1,0)
               STOP 'Marge'            
            ENDIF
C
C           Set dT to 0 on master loop if not solving it
            IF( .NOT. LakesToSolve(RtLoop2LkSlvI(L)) .AND. L .EQ. 
     &         NumRtLoops) dT(L) = ZERO
C
         ELSE
C           Not doing this lake.
            dT(L) = ZERO
            RoutingInterval(L) = RoutingInterval(NumRtLoops) 
         ENDIF
C
C
C        IF(GenVerbosity .GE. 9) THEN
C           WRITE(CTMP,8976)  L, dT(L), Day(L), 
C    &         PeriodEnd(1,L), MonthAbbrev(PeriodEnd(3,L)), 
C    &         PeriodEnd(4,L), PeriodEnd(5,L), LLT
C           CALL LogIt(LOGD,1,1)
C8976       FORMAT(' Init loop #',I1, ' to dT=', F11.4,
C    &         ', start day=', I2, ', end date=',
C    &         I2, A3, I5, ', elapsed days=', I10, A1)
C        ENDIF
   82    CONTINUE
C
C
      IF( dT(2) .GT. 220000) THEN
C        Trial and error seems to indicate routing gets really 
C        risky with DT < 220000 seconds on middle lakes.  Haven't
C        looked into it on other lakes.  
         CALL ErrMess(1121,CHNKST,LOGCD)
      ENDIF
C
C

C     Lines commented with CX pertain to sliding ahead to end next Friday.
C     Paul believed strongly that is better to slide back to previous one.
C
CX    See if running Plan 1958D, which works off weeks defined 
CX    from Saturday at 0000 to Friday at 2400.  If so, then
CX    must slide PeriodEnd(?,3) ahead to the next Friday.
C
C     Must slide ending of initial routing period to previous Friday
      IF( LakesToSolve(5) .AND. UStLawrSolution  .EQ. 2 ) THEN
CX        I = 7 - DaysSinceSaturday
CX        IF( I .GT. 6) I = 0
C        WRITE(*,'(2X,I2,A,I2,A,I5,2I9)')  PeriodEnd(3,3), '/',
C    &   PeriodEnd(1,3),'/',PeriodEnd(4,3),PeriodEnd(5,3),PeriodEnd(6,3)
CX       PeriodEnd(5,3) = PeriodEnd(5,3) + I
CX       PeriodEnd(1,3) = PeriodEnd(1,3) + I
CX       IF( PeriodEnd(1,3) .GT. Q(4,3) ) THEN
CX          PeriodEnd(1,3) = PeriodEnd(1,3) - Q(4,3)
CX          PeriodEnd(2,3) = 1
CX          PeriodEnd(3,3) = PeriodEnd(3,3) + 1
CX          IF( PeriodEnd(3,3) .GT. 12) THEN
CX             PeriodEnd(3,3) = 1
CX             PeriodEnd(4,3) = PeriodEnd(4,3) + 1
CX          ENDIF
CX       ENDIF
         PeriodEnd(1,3) = P58DStartDay - 1
         IF( PeriodEnd(1,3) .LE. 0 ) THEN
            CALL GetQtrMonBounds(P58DStartYear, P58DStartMonth, Q(1,3))
            PeriodEnd(1,3) = PeriodEnd(1,3) + Q(4,3)
            PeriodEnd(2,3) = 4
            PeriodEnd(3,3) = P58DStartMonth - 1
            IF( PeriodEnd(3,3) .LE. 0 ) THEN 
               PeriodEnd(3,3) = 12
               PeriodEnd(4,3) = PeriodEnd(4,3) - 1
            ENDIF
         ENDIF
         WRITE(CTMP,3377)  PeriodEnd(1,3),
     &      MonthAbbrev(PeriodEnd(3,3)), PeriodEnd(4,3), LLT
         CALL LogIt(LOGCD,1,1)
      ENDIF
 3377 FORMAT(' First Ontario routing period starts after ',
     &   'midnight on Friday, ', I2, A3, I5, A1)
C
C
      DO 725, L=1,8
C
C        Initialize whatever output files are requested:
         IF( LakesToSolve(L) ) THEN
            DO 766, I = 1,2
               DO 766, J=RoutingInterval(LkSlv2RtLoopI(L)), 4
                  IF(WriteTSOutput(J,I,L)) THEN
                     IF(.NOT. InitTSOutputFile(J,I,L,     
     &                  RoutingInterval(LkSlv2RtLoopI(L)) )) GOTO 39
                  ENDIF
  766             CONTINUE
C    
C           Initialize daily, quarter-monthly and monthly output buffers.
C           Put year at start of first line (or block of lines) in the 
C           monthly and qmonthly output, etc.
C
C      CHARACTER*74  DailyOutputBuf(48,2,8)
C
C           Put year, month, qm at start of first line of daily output 
            DO 748, J = 1, 2
  748          WRITE(DailyOutputBuf(1,J,L)(1:9), '(I4,I3,I2)')  
     &            PeriodEnd(4,LkSlv2RtLoopI(L)), 
     &            PeriodEnd(3,LkSlv2RtLoopI(L)), 
     &            PeriodEnd(2,LkSlv2RtLoopI(L))
C
C           Write missing values for first 7 days in qmonth that run 
C           starts.  This will put a missing values in for the first
C           day(s) of run, but will be overwritten when output comes.
            DO 749, M = 1, 7
C              Figure out ending column for the day in I
               I = 9 + M * 8
               DO 749, J = 1, 2
  749             WRITE(DailyOutputBuf(1,J,L)(I-7:I),'(I8)') MissingVal 
C
            DO 753, J = 1, 2
               DO 753, M = 1, 4
  753             WRITE(QMonthOutputBuf(M,1,J,L)(1:7), '(I5,I2)')  
     &               PeriodEnd(4,LkSlv2RtLoopI(L)), M
C               WRITE(CTMP, '(3(6X,A,I2),A1)')'L=',L,'J=',J,'M=',M,LLT
C               CALL LogIt(LOGD,0,0)
C               WRITE(CTMP, '(A,A1)') QMonthOutputBuf(M,1,J,L)(1:6), LLT  
C               CALL LogIt(LOGD,0,0)
C  753          continue
C
            DO 752, M = 1, StartMonth
C              Figure out ending column for the month in I
               I = 7 + M * 8
               DO 752, K = 1, 4
                  DO 752, J = 1, 2
  752                WRITE( QMonthOutputBuf(K,1,J,L)(I-7:I), '(I8)') 
     &                  MissingVal  
C
C           Put year at start of first line in the monthly output 
            DO 754, J = 1, 2
  754          WRITE(MonthOutputBuf(1,J,L)(1:5), '(I5)')  
     &            PeriodEnd(4,LkSlv2RtLoopI(L))
C
C           Write missing value for leading months before run starts.
C           This will actually put a missing value in for the first
C           month, but it will be overwritten when the output comes.
            DO 755, M = 1, StartMonth
C              Figure out ending column for the month in I
               I = 5 + M * 8
               DO 755, J = 1, 2
  755             WRITE(MonthOutputBuf(1,J,L)(I-7:I),'(I8)') MissingVal 
C
         ENDIF
  725    CONTINUE
C
C     Figure out how many days for each year of run
      TotalDays = (EndYear - StartYear + 1) * 365
      DO 83, I=StartYear, EndYear
C        Add in leap years
   83    IF(LEAPYEAR(I)) TotalDays = TotalDays + 1
C     Subtract the months before StartMonth
      DO 84, I=1, StartMonth - 1
         CALL NUMDAYS (StartYear, I, K, J)
   84    TotalDays = TotalDays - J 
C     Subtract the days before StartDay
      TotalDays = TotalDays - StartDay + 1                 
C     Subtract the months after and including EndMonth
      DO 85, I=EndMonth, 12
         CALL NUMDAYS (EndYear, I, K, J)
   85    TotalDays = TotalDays - J
C     Add the days within EndMonth
      TotalDays = TotalDays + EndDay                 
      WRITE(CTMP, '(A,I11,A,A1)') 
     &   ' Routing through ',TotalDays,' days.',LLT
      CALL LogIt(LOGSDMPO,1,0)
C
C     TODO: experiment with VERY large runs.  loops should be
C     infinite, but test to see if some ints bust along the way.
      IF( TotalDays .GT. MaxInt .OR. TotalDays .LT. 1) THEN
         WRITE(CTMP,'(A, I9, A1)') 
     &      ' Too many days in simulation period: ', TotalDays, LLT
         CALL LogIt(LOGCDE,0,0)
         CALL SAVETOD
         CALL ADIOS
      ENDIF
C
C     
      IF( LakesToSolve(1) ) THEN
C        Write some column headers to Superior output file   
         WRITE(CTMP, 8082) LLT
         CALL LogIt(LOGP,1,0)
         IF (StMarysSolution .EQ. 3 ) THEN
            WRITE(CTMP, '(A)') P77AHeader//LLT
         ELSE IF (StMarysSolution .EQ. 4) THEN
            WRITE(CTMP, '(A)') P2012Header//LLT
         ELSE 
            WRITE(CTMP, '(A)') P77AHeader//LLT
C           Maybe add some alternate header if different plan
         ENDIF
         CALL LogIt(LOGP,0,0)
      ENDIF
 8082 FORMAT(15X, '____________Levels_____________', 7X, 
     &   '_____________________Flows________________', A1)
C             
C
C     
      IF( LakesToSolve(2) ) THEN
         WRITE(CTMP, 8069) LLT
         CALL LogIt(LOGM,0,0)
         WRITE(CTMP, 8070) LLT
         CALL LogIt(LOGM,0,0)
      ENDIF
C             
C     
      IF( LakesToSolve(5) ) THEN
         WRITE(CTMP, 8654) OutputHeader, LLT
         CALL LogIt(LOGO,0,0)
      ENDIF
C             
C       
      IF(GenVerbosity .GE. 2 ) THEN
         WRITE(CTMP,'(A,A1)')
     &      ' !!! Beginning routing computations . . .',LLT
         CALL LogIt(LOGCD,1,0)
      ENDIF
C
C     Initialize progress counter
      CompletionStage(1) = CompletionStage(2) * TotalDays / 10 
C
C
C     Begin the control (main) loops for routing calculations.  
C     Use PeriodEnd(5,?) as counter, if go all the way through MaxInt
C     then jumps back and starts control loop again.  Only
C     way out of loop is to reach end date, in which case flow
C     jumps out of loop.  So be darn sure that end date is legit!
C
C     The dates and IF's used here for loop control take a long
C     time to understand.  If you muck around with them, then you
C     will likely BREAK THE LOGIC.  Do a lot of testing with all
C     combinations of routing periods before locking in any changes.
C
C     This outer loop is the Lake Ontario Loop 
C
 8002 DO 1333, OntarioPeriod = OntarioPeriod + 1, MaxInt
C
         IF(GenVerbosity .GE. 4 ) THEN
            CTMP='               vvvvvvv   New master routing'//
     &         ' period starts here:  vvvvvvv'//LLT 
            CALL LogIt(LOGD,1,0)
         ENDIF
C
         CALL IncrementDates(3, SecsInRtPeriod(3),
     &      NewWeek(3), NewQMonth(3), NewMonth(3))
C
C        If we already have Niagara River flows for next Lake
C        Ontario loop, go there now:
         IF(PeriodEnd(5,3) .LE. PeriodEnd(5,2) )GOTO 1453
C
         DO 2333, MiddlePeriod = MiddlePeriod + 1, MaxInt
C
 8018       CALL IncrementDates(2, SecsInRtPeriod(2),
     &         NewWeek(2), NewQMonth(2), NewMonth(2))
C
            IF( LakesToSolve(1) .AND. 
     &         PeriodEnd(5,2) .GT. PeriodEnd(5,1) )THEN
C
C               WRITE(CTMP,'(A,I3,A,I3,A1)')
C     &            ' Inside 8021: PeriodEnd(5,2)=', PeriodEnd(5,2),
C     &             '  ?GT? PeriodEnd(5,1)=', PeriodEnd(5,1), LLT
C               CALL LogIt(LOGD,0,0)
C
               DO 3333, SuperiorPeriod = SuperiorPeriod + 1, MaxInt
C     
                  CALL IncrementDates(1, SecsInRtPeriod(1),
     &               NewWeek(1), NewQMonth(1), NewMonth(1))
C
C                 If passed ending date on Lake Superior loop, keep
C                 same St. Marys' flow 
                  IF(FinishedLoop(1)) GOTO 8013
C
C                 See if run ends during this Superior interval
                  IF( PeriodEnd(4,1) .GT. EndYear .OR.
     &               (PeriodEnd(4,1) .EQ. EndYear .AND.
     &                   PeriodEnd(3,1) .GT. EndMonth) .OR.
     &               (PeriodEnd(4,1) .EQ. EndYear .AND.
     &                   PeriodEnd(3,1) .EQ. EndMonth) .AND.
     &                   PeriodEnd(1,1) .GE. EndDay)
     &               FinishedLoop(1) = .TRUE.
C
C                 K indicates which summary file (lake) is being read,
C                 corresponds to index of LkToSlv2LN.  Set value prior 
C                 to call to READ 'cuz if problem will jump to 77 or 79,
C                 and want to know which lake went wrong.  Similar 
C                 mechanism for J, which tells which readformat was 
C                 attempted.
                  K = 1      
                  J = 3
8010              READ(81, 8653, IOSTAT=I, ERR=77) ReadBuf
                  IF(ReadBuf(1:1) .EQ. '#') GOTO 8010
                  READ(ReadBuf,ReadFormat(J),IOSTAT=I, ERR=79) SPNTS, 
     &               SMICE, NumSupForceGates, SupSideChannel, 
     &               SupForceFlow, SupDevFlow, SupForceEOPLev, 
     &               SupIndex1, SupIndex2, SupIndex3
                  CALL RouteSuperior(NewMonth(1), I)
                  IF (I .GT. 0) NumSPSolveError = NumSPSolveError + I
                  CALL TallyMeans(1, 1, STMRIVFLOW, SPMLV)
8013              IF( LakesToSolve(2) ) CALL ConChanConvert( 
     &               1, STMRIVFLOW, StMarysFlowAccum,
     &               StMarysDaysAccum, StMarysConverted)
C
                  WBCheckSP = WBCheckSP + 
     &               (SPNTS - STMRIVFLOW) * SecsInRtPeriod(1)
                  deltaStor = (SPEOP - SPSTLV) *
     &               ((GETAREA(1,SPEOP) + GETAREA(1,SPSTLV)))/2
                  IF(SupVerbosity .GE. 7 ) THEN
                     WRITE(CTMP, 8309)'Superior', PeriodEnd(1,1), 
     &                  MonthAbbrev(PeriodEnd(3,1)),
     &                  PeriodEnd(4,1), WBCheckSP, SPSTLV, SPEOP, 
     &                  deltaStor, deltaStor - WBCheckSP,
     &                  (deltaStor - WBCheckSP)/GETAREA(1,SPSTLV),
     &                  (deltaStor - WBCheckSP)/(86400*PeriodEnd(5,1)), 
     &                  LLT
                     CALL LogIt(LOGD,0,1)
                     IF(SupVerbosity .GE. 7 ) THEN
                        WRITE(CTMP,'(A,I10,A1)')
     &                     ' Superior outflow computed up to day ', 
     &                     PeriodEnd(5,1), LLT
                        CALL LogIt(LOGD,0,0)
                        WRITE(CTMP,'(A,I10,A,A1)')
     &                     ' Need Superior outflow computed up to day ', 
     &                     PeriodEnd(5,2), 
     &                     ' to route middle lakes.', LLT
                        CALL LogIt(LOGD,0,0)
                     ENDIF
                  ENDIF
C
C                 Jump out of loop #1 if have StMarys flows up to or 
C                 past the end of this Middle Lakes routing period.
 3333             IF(PeriodEnd(5,1) .GE. PeriodEnd(5,2) )GOTO 3353
C
            ENDIF
C----&------------------------------------------------------------------
C
 3353       IF( LakesToSolve(2) )THEN
C
               IF(MLVerbosity .GE. 7 ) THEN
                  WRITE(CTMP,'(A,I10,A1)')
     &              ' Beginning middle lake routing for period ending'//
     &              ' on day ', PeriodEnd(5,2), LLT
                  CALL LogIt(LOGD,0,0)
                  WRITE(CTMP,'(A,I10,A1)')
     &                  ' Superior outflows available up to day ', 
     &                  PeriodEnd(5,1), LLT
                  CALL LogIt(LOGD,0,0)
                  WRITE(CTMP,'(A,I10,A1)')
     &               ' Must loop through middle lakes until day ', 
     &               PeriodEnd(6,3), LLT
                  CALL LogIt(LOGD,0,0)
               ENDIF
C
C              If passed ending date on middle lakes loop, keep
C              same Niagara flow 
                  IF(FinishedLoop(2)) GOTO 8014
C
C              See if run ends during this middle lakes interval 
               IF( PeriodEnd(4,2) .GT. EndYear .OR.
     &            (PeriodEnd(4,2) .EQ. EndYear .AND.
     &                PeriodEnd(3,2) .GT. EndMonth) .OR.
     &            (PeriodEnd(4,2) .EQ. EndYear .AND.
     &                PeriodEnd(3,2) .EQ. EndMonth) .AND.
     &                PeriodEnd(1,2) .GE. EndDay ) 
     &               FinishedLoop(2) = .TRUE.
C
               K = 2      
               J = 2
 4011          READ(82, 8653, IOSTAT=I, ERR=77) ReadBuf
               IF(ReadBuf(1:1) .EQ. '#') GOTO 4011
               READ(ReadBuf,ReadFormat(J),IOSTAT=I,ERR=79) 
     &            SupplyMH, MHICE, MHForceEOPLev,ussiceweedadj,tmpD
               K = 3      
 4012          READ(83, 8653, IOSTAT=I, ERR=77) ReadBuf
               IF(ReadBuf(1:1) .EQ. '#') GOTO 4012
               J = 1
               READ(ReadBuf, ReadFormat(J), IOSTAT=I, ERR=79)
     &            SCNBS, SCICE, STCForceEOPLev
               K = 4      
 4013          READ(84, 8653, IOSTAT=I, ERR=77) ReadBuf
               IF(ReadBuf(1:1) .EQ. '#') GOTO 4013
               J = 2
               READ(ReadBuf,ReadFormat(J),IOSTAT=I,ERR=79) 
     &            ERNBS, ERICE, CGMLV, ErieForceEOPLev, FlowIntoWelland
C              Subtract Welland from ERNBS here for correct water balance
               ERNBS = ERNBS - FlowIntoWelland
C
C              Kludge that follows only for backward compatibility
C              with existing rounding practice.  For some reason, if doing
C              quarter-monthly or monthly routing, the end-of-month values
C              are rounded to the cm for the middle lakes.
               IF(BOMRound .GT. -5  .AND. NewMonth(2)  .AND. 
     &            RoutingInterval(2) .GE. 3)THEN
                  MHEOP = Round( MHEOP, BOMRound )
                  SCEOP = Round( SCEOP, BOMRound )
                  EREOP = Round( EREOP, BOMRound )
               ENDIF
C
C----&------------------------------------------------------------------
               IF( LakesToSolve(1) ) THEN
C                 Must add Superior outflows to MHNTS here.  If Superior
C                 is monthly routing and middle lake are weekly routing
C                 and the week is contained within the month, then use
C                 SPMonthlyMeanF as StMarysConverted.  Otherwise use the
C                 value computed by ConChanConvert for a week that 
C                 straddles two months. 
                  IF( RoutingInterval(2) .EQ. 2 .AND. 
     &               RoutingInterval(1) .EQ. 4 .AND. 
     &               PeriodEnd(1,2) .GT. 6) 
     &               StMarysConverted = SPMonthlyMeanF
                  IF(MLVerbosity .GE. 4 ) THEN
                     WRITE(CTMP,'(A,I9,A,I2,A3,I5,A,A1)') 
     &                  ' On day number', PeriodEnd(5,2),
     &                  ', set Michigan-Huron NTS for period ending ',
     &                  PeriodEnd(1,2), MonthAbbrev(PeriodEnd(3,2)), 
     &                  PeriodEnd(4,2), ' as', LLT
                     CALL LogIt(LOGD,1,0)
                     WRITE(CTMP,'(3(A,F18.9),A1)')' M-H NTS =',
     &                  StMarysConverted, ' + ', SupplyMH, ' = ', 
     &                  SupplyMH + StMarysConverted, LLT
                     CALL LogIt(LOGD,0,0)
                  ENDIF
                  MHNTS = SupplyMH + StMarysConverted
               ELSE
C                 Not computing St Marys flow; it has been read in and
C                 already lumped together with SupplyMH
                  MHNTS = SupplyMH
               ENDIF
C
               CALL RouteMiddleLakes(K)
               IF (K .GT. 0) NumMLSolveError = NumMLSolveError + K
               CALL TallyMeans(2, 2, STCRIVFLOW, MHMLV)
               CALL TallyMeans(3, 2, DETRIVFLOW, SCMLV)
               IF(WriteErieOutflow) THEN
                 CALL TallyMeans(4,2,NIARIVFLOW + FlowIntoWelland,ERMLV)
               ELSE
                 CALL TallyMeans(4, 2, NIARIVFLOW, ERMLV)
               ENDIF
 8014          IF( LakesToSolve(5) ) CALL ConChanConvert( 
     &            2, NIARIVFLOW, NiagaraFlowAccum,
     &            NiagaraDaysAccum, NiagaraConverted)
C
               WBCheckMH = WBCheckMH + 
     &            (MHNTS - STCRIVFLOW) * SecsInRtPeriod(2)
               WBCheckSC = WBCheckSC + (STCRIVFLOW +
     &            SCNBS - DETRIVFLOW) * SecsInRtPeriod(2)
C              Welland flows already deducted from ERNBS at this point
               WBCheckER = WBCheckER + (DETRIVFLOW + ERNBS - 
     &             NIARIVFLOW) * SecsInRtPeriod(2)
               IF(MLVerbosity .GE. 7 ) THEN
                  deltaStor = (MHEOP - MHSTLV) *
     &               ((GETAREA(2,MHEOP) + GETAREA(2,MHSTLV)))/2
                  WRITE(CTMP, 8309)'Mich-Huron',  
     &               PeriodEnd(1,2), MonthAbbrev(PeriodEnd(3,2)), 
     &               PeriodEnd(4,2), WBCheckMH, MHSTLV, MHEOP, 
     &               deltaStor, deltaStor - WBCheckMH,
     &               (deltaStor - WBCheckMH)/GETAREA(2,MHSTLV),
     &               (deltaStor - WBCheckMH)/(86400*PeriodEnd(5,2)), 
     &               LLT
                  CALL LogIt(LOGD,0,0)
                  deltaStor = (SCEOP - SCSTLV) *
     &               ((GETAREA(5,SCEOP) + GETAREA(5,SCSTLV)))/2
C                 WRITE(CTMP,'(4(A,D17.7),A1)') '  SCEOP=',SCEOP,
C    &               ' A2=',GETAREA(5,SCEOP),' SCSTLV=', SCSTLV,' A=1',
C    &               GETAREA(5,SCSTLV),LLT 
C                 CALL LogIt(LOGD,0,0)
                  WRITE(CTMP, 8309)'Lake St. Clair', 
     &               PeriodEnd(1,2), MonthAbbrev(PeriodEnd(3,2)), 
     &               PeriodEnd(4,2), WBCheckSC, SCSTLV, SCEOP, 
     &               deltaStor, deltaStor - WBCheckSC,
     &               (deltaStor - WBCheckSC)/GETAREA(5,SCSTLV),
     &               (deltaStor - WBCheckSC)/(86400*PeriodEnd(5,2)), 
     &               LLT
                  CALL LogIt(LOGD,0,0)
                  deltaStor = (EREOP - ERSTLV) *
     &               ((GETAREA(6,EREOP) + GETAREA(6,ERSTLV)))/2
                  WRITE(CTMP, 8309)'Erie',  
     &               PeriodEnd(1,2), MonthAbbrev(PeriodEnd(3,2)), 
     &               PeriodEnd(4,2), WBCheckER, ERSTLV, EREOP, 
     &               deltaStor, deltaStor - WBCheckER,
     &               (deltaStor - WBCheckER)/GETAREA(6,ERSTLV),
     &               (deltaStor - WBCheckER)/(86400*PeriodEnd(5,2)), 
     &               LLT
                  CALL LogIt(LOGD,0,0)
               ENDIF
            ENDIF
            IF(MLVerbosity .GE. 7 ) THEN
               WRITE(CTMP,'(A,I10,A1)')
     &            ' Middle Lakes outflow computed up to day ', 
     &            PeriodEnd(5,2), LLT
               CALL LogIt(LOGD,0,0)
               WRITE(CTMP,'(A,I10,A,A1)')
     &            ' Need middle lakes outflow computed up to day ', 
     &            PeriodEnd(5,3), ' to route next Ontario.', LLT
             CALL LogIt(LOGD,0,0)
            ENDIF
C           Jump out of loop #2 if have Niagara flows up to or past the
C           end of the Ontario routing period.
            IF(PeriodEnd(5,2) .GE. PeriodEnd(5,3) )GOTO 1453
C
            IF(MLVerbosity .GE. 7 ) THEN
               WRITE(CTMP,'(A,I10,A1)')
     &            ' StMarys flow available to day ', 
     &            PeriodEnd(5,1), LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
 2333       IF(PeriodEnd(6,2) .GT. PeriodEnd(5,1) )GOTO 8018
C
C1453    CONTINUE
 1453    IF( LakesToSolve(5) )THEN
C
            K = 5      
            J = 4
 8015       READ(85, 8653, IOSTAT=I, ERR=77) ReadBuf
            IF(ReadBuf(1:1) .EQ. '#') GOTO 8015
            READ(ReadBuf, ReadFormat(J), IOSTAT=I, ERR=79) SupplyON, 
     &          ONICE, NiagaraRead, FlowFromWelland, OntDeviation, 
     &          LakeOntarioFlow, StLouisRetard, StLouisFlow,LSRoughness,
     &          DesPrFlow, MilleFlow, NumOntForceGates, OntForceEOPLev, 
     &          OntForceFlow, LouisOntFlowDiff,
     &          Qdpmi, Qstfran, Qstmau, QRich,
     &          RJetty1, RVarennes, RSorel, RPierre,
     &          R3Rivers, RBatiscan, Tidal, Onttnbs,
     &          Other06n, Other06Q
            WRITE(CTMP, '(15F12.4,A1)', IOSTAT=I, ERR=79) SupplyON, 
     &          ONICE, NiagaraRead, FlowFromWelland, OntDeviation, 
     &          LakeOntarioFlow, StLouisRetard, StLouisFlow,LSRoughness, 
     &          DesPrFlow, MilleFlow, NumOntForceGates, OntForceEOPLev, 
     &          OntForceFlow, LouisOntFlowDiff, LLT
            CALL LogIt(LOGD,1,0)
            WRITE(CTMP, '(14F12.4,A1)', IOSTAT=I, ERR=79) 
     &          Qdpmi, Qstfran, Qstmau, QRich,
     &          RJetty1, RVarennes, RSorel, RPierre,
     &          R3Rivers, RBatiscan, Tidal, Onttnbs,
     &          Other06n, Other06Q, LLT
            CALL LogIt(LOGD,1,0)
C
            IF( LakesToSolve(4) ) THEN
C              Must add Niagara outflows to SupplyON here.  If Erie
C              is monthly routing and Ontario is weekly routing
C              and the week is contained within the month, then use
C              ERMonthlyMeanF as NiagaraConverted.  Otherwise 
C              NiagaraConverted remains set to the value computed by 
C              ConChanConvert (for a week that straddles two months).
C               
               IF( RoutingInterval(3) .EQ. 2 .AND. 
     &             RoutingInterval(2) .EQ. 4 .AND. 
     &             PeriodEnd(1,3) .GT. 6) THEN
                       NiagaraConverted = ERMonthlyMeanF
               ENDIF
            ELSE
C              Not computing Niagara flow; it has been read in
               NiagaraConverted = NiagaraRead
            ENDIF
C
C           Checking to for Lake Ontario Supplies data
            IF(filenames(89)(1:12).ne.Unspecified)THEN
C              Ont Total Supplies Data exists!
               ONNTS=Onttnbs
               SupplyON=Onttnbs - (NiagaraConverted + FlowFromWelland)
            ELSE
C              Ont LOCAL Supplies Data exists!
               ONNTS = SupplyON + NiagaraConverted + FlowFromWelland
            ENDIF
            IF(OntVerbosity .GE. 4 ) THEN
               WRITE(CTMP,'(A,I9,A,I2,A3,I5,A,A1)') 
     &             ' On day number', PeriodEnd(5,3),
     &             ', set Ontario NTS for period ending ',
     &             PeriodEnd(1,3), MonthAbbrev(PeriodEnd(3,3)), 
     &             PeriodEnd(4,3), ' as', LLT
               CALL LogIt(LOGD,1,0)
               IF(filenames(53).ne.Unspecified) THEN
                   WRITE(CTMP,'(4(A,F18.9),A1)')' ON NTS =',
     &                 NiagaraConverted, ' + ', SupplyON, ' + ', 
     &                 FlowFromWelland,' = ', ONNTS, LLT
               ELSE
                   WRITE(CTMP,'(A,F18.9,A1)')' ON NTS =',
     &                 Onttnbs, LLT
               ENDIF
               CALL LogIt(LOGD,0,0)
            ENDIF
C
            CALL OntarioPlans(Day(3), PeriodEnd(1,3), PeriodEnd(2,3), 
     &         PeriodEnd(3,3), PeriodEnd(4,3), PeriodEnd(5,3),
     &         NewWeek(3), NewQMonth(3), NewMonth(3), SecsInRtPeriod(3),
     &         NiagaraConverted, FlowFromWelland, SupplyON, ONNTS, 
     &         ONICE, OntDeviation, StLouisRetard, StLouisFlow, 
     &         LSRoughness, DesPrFlow, MilleFlow, NumOntForceGates, 
     &         OntForceEOPLev, OntForceFlow ,LouisOntFlowDiff, 
     &         Qdpmi, Qstfran, Qstmau, QRich, RJetty1, RVarennes, 
     &         RSorel, RPierre, R3Rivers, RBatiscan, Tidal, Onttnbs,
     &         Other06n, Other06Q,
     &         LakeOntarioFlow, ONBOP, ONEOP, ONMLV)
C
         ENDIF
C
C----&------------------------------------------------------------------
C
C        THIS IS THE ONLY WAY OUT!  DON'T BUNGLE ENDING DATE!
C        Jump out of loop if was last day of final period
         IF( PeriodEnd(4,NumRtLoops) .GT. EndYear .OR.
     &      (PeriodEnd(4,NumRtLoops) .EQ. EndYear .AND.
     &         PeriodEnd(3,NumRtLoops) .GT. EndMonth) .OR.
     &      (PeriodEnd(4,NumRtLoops) .EQ. EndYear .AND.
     &         PeriodEnd(3,NumRtLoops) .EQ. EndMonth) .AND.
     &         PeriodEnd(1,NumRtLoops) .GE. EndDay ) GOTO 8990
C
C
C            WRITE(CTMP,'(I12,A,A1)')PeriodEnd(5,NumRtLoops), 
C     &          ' days routed',        LLT
C            CALL LogIt(LOGCD,0,0)
         IF(GenVerbosity .GE. 2 .AND. PeriodEnd(5,NumRtLoops) .GE.
     &      CompletionStage(1)) THEN
            WRITE(CTMP,'(A,I2,A,A1)')'     . . . routing computations ',
     &         CompletionStage(2), '0% complete . . .', LLT
            CALL LogIt(LOGCD,0,0)
            CompletionStage(2) = CompletionStage(2) + 2
            CompletionStage(1) = CompletionStage(2) * TotalDays / 10 
         ENDIF
C
         IF(GenVerbosity .GE. 5 ) THEN
            WRITE(CTMP,8966)  PeriodEnd(1,NumRtLoops),
     &         MonthAbbrev(PeriodEnd(3,NumRtLoops)),
     &         PeriodEnd(4,NumRtLoops), LLT
 8966       FORMAT(15X, '^^^^^^^   Master period ending ', I2, A3, I5, 
     &         ' passes   ^^^^^^^', A1)
            CALL LogIt(LOGD,1,1)
         ENDIF
C
         IF( PeriodEnd(5,NumRtLoops) .GT. MAXINT ) THEN
C           Must recycle days elapsed counters.  If 
C           PeriodEnd(5,NumRtLoops) is past MAXINT, then all of the
C           elapsed days counters being used are over also.
            DO 72, I = 1, NumRtLoops
               IF( LakesToSolve(RtLoop2LkSlvI(I)) .OR. 
     &            I .EQ. NumRtLoops )THEN
                  PeriodEnd(5,I) = PeriodEnd(5,I) - MAXINT
                  PeriodEnd(6,I) = PeriodEnd(6,I) - MAXINT
               ENDIF
   72       CONTINUE
         ENDIF
C----&------------------------------------------------------------------
C
C
 1333    CONTINUE  
C            
C     Hmmm, if get here must be lotsa periods.  oh, well, maybe a 
C     humongous run.  jump back to beginning of loop.
      GOTO 8002 
C        
C            
 8990 CONTINUE    
C
C
C     Write out anything left in holding space.  
C
      DO 609, L = 1, 8
         DO 609, I = 1, 4
C           If routed through the lake, and want output, and output 
C           interval is greater than or equal to routing interval, then do it.
            IF( LakesToSolve(L)) J = LkSlv2RtLoopI(L)
  609       IF( LakesToSolve(L) .AND. WriteTSOutput(I,1,L) .AND.
     &         RoutingInterval(J) .LE. I) 
     &         CALL WriteToDisk( I, L )
C
C
C
C      IF(GenVerbosity .GE. 4) THEN
C         WRITE(CTMP,8654)'  . . . complete', LLT
C         CALL LogIt(LOGCD,0,0)
C      ENDIF
C
C     If got continuity-check errors during superior or middle-lakes routing,
C     issue warnings.
      IF(NumSPSolveError .GT. 0 .OR. NumMLSolveError .GT. 0) 
     &   CALL ErrMess(1090,CHNKST,LOGCD)
      IF(NumSPSolveError .GT. 0) THEN
         WRITE(CTMP,8352) 'Lake Superior', SolutionTest, 
     &      NumSPSolveError, LLT
         CALL LogIt(LOGCD,0,0)
         WRITE(CTMP, 8353) 'Lake Superior', SPErrorFlow, LLT
         CALL LogIt(LOGCD,0,1)
      ENDIF
      IF(NumMLSolveError .GT. 0) THEN
         WRITE(CTMP,8352) 'middle lakes', SolutionTest, 
     &      NumMLSolveError, LLT
         CALL LogIt(LOGCD,0,0)
         WRITE(CTMP, 8353) 'middle lakes', MLErrorFlow, LLT
         CALL LogIt(LOGCD,0,1)
      ENDIF
C
      IF(Extrapolations .GT. 0) THEN
         CALL ErrMess(1183,CHNKST,LOGCD)
         WRITE(CTMP,8354) Extrapolations, LargestExtrapolation, LLT
         CALL LogIt(LOGCD,0,0)
      ENDIF

      IF(numMaxSC .GT. 0 .OR. numMinSCWinter .GT. 0 .OR. 
     &   numMinSCSummer .GT. 0 .OR. numJacSC .GT. 0) THEN
         WRITE(CTMP,8367) 'physical maximum capacity', numMaxSC, LLT
         CALL LogIt(LOGCD,0,0)
         WRITE(CTMP,8367) 'winter minimum capacity', numMinSCWinter,LLT
         CALL LogIt(LOGCD,0,0)
         WRITE(CTMP,8367) 'summer minimum capacity', numMinSCSummer,LLT
         CALL LogIt(LOGCD,0,0)
         WRITE(CTMP,8367)'physical Jacob maximum capacity',numJacSC,LLT
         CALL LogIt(LOGCD,0,0)
      ENDIF
C
C     Put out cumulative continuity results regardless of verbosity
      deltaStor = (SPEOP - SPSTLV) *
     &    ((GETAREA(1,SPEOP) + GETAREA(1,SPSTLV)))/2
      WRITE(CTMP, 8309)'Superior', PeriodEnd(1,1),
     &   MonthAbbrev(PeriodEnd(3,1)), PeriodEnd(4,1), WBCheckSP, SPSTLV, 
     &   SPEOP, deltaStor, deltaStor - WBCheckSP,
     &   (deltaStor - WBCheckSP)/GETAREA(1,SPSTLV),
     &   (deltaStor - WBCheckSP)/(86400*PeriodEnd(5,1)), LLT
      CALL LogIt(LOGD,0,1)
C
      deltaStor = (MHEOP - MHSTLV) *
     &   ((GETAREA(2,MHEOP) + GETAREA(2,MHSTLV)))/2
      WRITE(CTMP, 8309)'MichHuron', PeriodEnd(1,2),
     &   MonthAbbrev(PeriodEnd(3,2)), PeriodEnd(4,2), WBCheckMH, MHSTLV, 
     &   MHEOP, deltaStor, deltaStor - WBCheckMH,
     &   (deltaStor - WBCheckMH)/GETAREA(2,MHSTLV),
     &   (deltaStor - WBCheckMH)/(86400*PeriodEnd(5,2)), LLT
      CALL LogIt(LOGD,0,0)
C
      deltaStor = (SCEOP - SCSTLV) *
     &   ((GETAREA(5,SCEOP) + GETAREA(5,SCSTLV)))/2
      WRITE(CTMP, 8309)'St. Clair', PeriodEnd(1,2), 
     &   MonthAbbrev(PeriodEnd(3,2)), PeriodEnd(4,2), WBCheckSC, SCSTLV,
     &   SCEOP, deltaStor, deltaStor - WBCheckSC,
     &   (deltaStor - WBCheckSC)/GETAREA(5,SCSTLV),
     &   (deltaStor - WBCheckSC)/(86400*PeriodEnd(5,2)), LLT
      CALL LogIt(LOGD,0,0)
C
      deltaStor = (EREOP - ERSTLV) *
     &   ((GETAREA(6,EREOP) + GETAREA(6,ERSTLV)))/2
      WRITE(CTMP, 8309)'Erie', PeriodEnd(1,2),
     &   MonthAbbrev(PeriodEnd(3,2)), PeriodEnd(4,2), WBCheckER, ERSTLV, 
     &   EREOP, deltaStor, deltaStor - WBCheckER,
     &   (deltaStor - WBCheckER)/GETAREA(6,ERSTLV),
     &   (deltaStor - WBCheckER)/(86400*PeriodEnd(5,2)), LLT
      CALL LogIt(LOGD,0,0)
C
C     Format statements.
 8069 FORMAT(32X,'Lakes Michigan-Huron', 17X, 
     &   'Lake St. Clair', 15X, 'Lake Erie', A1 )
 8070 FORMAT('   DayNum SD ED  M   Yr     US Slip      Mean      ', 
     &   'EOP StC Riv      Mean      EOP Det Riv      Mean      ',
     &   'EOP Nia Riv   CGIP Lvl', A1)
 8309 FORMAT(A, ' continuity check at end of ',I2,A3,I5,
     &   ' cum water bal=',D16.9, ' m3, init el=',F9.4, ' m, final el=',
     &   F9.4, ' m, stor change=',D15.9, ' m3, cum discontinuity=',
     &   D15.9, ' m3 =', F12.9,' m = ', F12.6, ' m3/s', A1)

 8352 FORMAT(' Error of ', A, ' solution exceeded ', F5.2, 
     &   ' m3s a total of ',I9,' times.', A1)
 8353 FORMAT(' Maximum error of solution in ', A, ':', F9.3,' m3s.', A1)
 8354 FORMAT(' Extrapolated a total of ',I9,
     &   ' times.  Greatest distance: ', F21.9,A1)
 8367 FORMAT(' Number of times Superior sidechannel overridden by ',
     &   A, ': ', I9, A1)
 8653 FORMAT(A)
 8654 FORMAT(A,A1)
 9015 FORMAT(' Got IOSTAT=',I5,' when accessing file ==>',A,'<==.',A1)
C
      END
C
C
C
C
      SUBROUTINE IncrementDates(Lake, SecsInPeriod, 
     &    NewWeek, NewQMonth, NewMonth)
C
C
      IMPLICIT NONE 
C
C     Declarations regarding function arguments:
C
      INTEGER       Lake
C                   1 for Lake Superior, 2 for middle lakes, 3 for
C                   Lake Ontario
C
      INTEGER       SecsInPeriod
C                   Number of seconds in the routng period
C
      LOGICAL       NewWeek, NewQMonth, NewMonth
C                   True if period starts a new one
C
C     Declarations regarding local variables or functions:
C
      INTEGER       I
C                   Scratch
C
      INTEGER       DaysToAdd(2)
C                   For daily and weekly routing. 
C
      CHARACTER*13  NameOfLake(3)
C
C
C
C     Declarations regarding COMMON block variables:
C
      INTEGER       PeriodEnd(6,3)
      INTEGER       Day(3), Q(4,3), DaysSinceSaturday
      LOGICAL       NumIncsVary(3) 
C                   See descriptions in subroutine RUNIT().
C
      DOUBLE PRECISION dT(3)
C                   "delta time" associated with routing interval (or
C                   increment).  Used in storage-change calculations.
C                   First element is for Lake Superior, second for
C                   middle lakes, third for Lake Ontario.
C
      INTEGER       RoutingInterval(3), RoutingIncrement(3)
C                   Routing interval to use on Superior, middle lakes,
C                   and Ontario (1=daily, 2=weekly, 3=qmonthly, 
C                   4=monthly), and the number of increments within each 
C                   interval.
C
      LOGICAL       LakesToSolve(8) 
      INTEGER       StartYear, StartMonth, StartDay, 
     &              EndYear, EndMonth, EndDay
      INTEGER       IntervalNameLen(4)
      INTEGER       MidLakeSolution, StMarysSolution, UStLawrSolution, 
     &              P77ACompat
      INTEGER       MonthLength
      INTEGER       BOMRound, FlowRound, LevelRound
C                   Not used here; just part of same COMMON 
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
      CHARACTER*3   MonthAbbrev(12)
      CHARACTER     CTMP*255, CWRK*80
      CHARACTER     LLT*1
C
      COMMON /CLDATE/   dT, Day, Q, PeriodEnd, NumIncsVary,
     &                  DaysSinceSaturday
      COMMON /RCRUNI/   LakesToSolve,  
     &                  IntervalNameLen, MonthLength,
     &                  StartYear, StartMonth, StartDay, 
     &                  EndYear, EndMonth, EndDay
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
      COMMON /MONAME/   MonthAbbrev
C        
      DATA DaysToAdd / 1, 7 /
C
      DATA NameOfLake /'Lake Superior','Middle Lakes ','Lake Ontario '/
C
C
C     Update day (beginning of new routing period) to be the day after 
C     current period ends
      Day(Lake) = PeriodEnd(1,Lake) + 1
C     If past end of a month, adjust date
      IF(Day(Lake) .GT. Q(4,Lake)) Day(Lake) = 1
C
      IF( RoutingInterval(Lake) .LE. 2 ) THEN
C
C        Daily or weekly routing - increment by DaysToAdd
C        Compute elapsed days at end of new routing interval
         PeriodEnd(5,Lake) = PeriodEnd(5,Lake) + 
     &      DaysToAdd(RoutingInterval(Lake))
C        Compute elapsed days at end of routing interval after new one
         PeriodEnd(6,Lake) = PeriodEnd(5,Lake) + 
     &      DaysToAdd(RoutingInterval(Lake))
C
C        Update day of month for end of new period
         PeriodEnd(1,Lake) = PeriodEnd(1,Lake) + 
     &      DaysToAdd(RoutingInterval(Lake))
C
C
         IF( RoutingInterval(Lake) .EQ. 1 ) THEN
C           If daily routing and after first period, turn off NewMonth, etc
            IF( PeriodEnd(5,Lake) .GT. 1 )THEN
               NewWeek = .FALSE.
               NewQMonth = .FALSE.
               NewMonth = .FALSE.
C              If were routing daily through Ontario under P58D, then
C              would need to add DaysSinceSaturday to PeriodEnd(5,?)
               IF( MOD(PeriodEnd(5,Lake), 7) .EQ. 0) NewWeek = .TRUE.
            ENDIF
         ELSE
            NewWeek = .TRUE.
         ENDIF
C
C        See if new ending date is greater than end of current
C        quartermonth number=PeriodEnd(2,Lake)
         IF(PeriodEnd(1,Lake) .GT. Q(PeriodEnd(2,Lake),Lake)) THEN
            NewQMonth = .TRUE.
            PeriodEnd(2,Lake) = PeriodEnd(2,Lake) + 1
C
C           Check to see if new month
            IF(PeriodEnd(2,Lake) .GT. 4) THEN
               NewMonth = .TRUE.
               PeriodEnd(2,Lake) = 1
C              Adjust date
               PeriodEnd(1,Lake) = PeriodEnd(1,Lake) - Q(4,Lake)
               PeriodEnd(3,Lake) = PeriodEnd(3,Lake) + 1
C
C              Check to see if new year
               IF(PeriodEnd(3,Lake) .GT. 12 ) THEN
                  PeriodEnd(3,Lake) = 1
                  PeriodEnd(4,Lake) = PeriodEnd(4,Lake) + 1
               ENDIF
C
C              Update Q array for the lake
               CALL GetQtrMonBounds( PeriodEnd(4,Lake), 
     &            PeriodEnd(3,Lake), Q(1,Lake))
            ELSE
               NewMonth = .FALSE.
            ENDIF
         ELSE
            NewQMonth = .FALSE.
            NewMonth = .FALSE.
         ENDIF
C        Oops!  Check to see if first time through:
         IF(PeriodEnd(5,Lake) .EQ. DaysToAdd(RoutingInterval(Lake)))THEN
            NewQMonth = .TRUE.
            NewMonth = .TRUE.
            NewWeek = .TRUE.
         ENDIF          
C
      ELSE IF( RoutingInterval(Lake) .EQ. 3 ) THEN
C        Quarter-monthly routing
         NewWeek = .TRUE.
         NewQMonth = .TRUE.
         PeriodEnd(2,Lake) = PeriodEnd(2,Lake) + 1
C
C        Check to see if new month (or 1st month of run)
         IF(PeriodEnd(2,Lake) .GT. 4 ) THEN
            NewMonth = .TRUE.
            PeriodEnd(2,Lake) = 1
            PeriodEnd(3,Lake) = PeriodEnd(3,Lake) + 1
C           Check to see if new year
            IF(PeriodEnd(3,Lake) .GT. 12 ) THEN
               PeriodEnd(3,Lake) = 1
               PeriodEnd(4,Lake) = PeriodEnd(4,Lake) + 1
            ENDIF
C           Update Q array for the lake
            CALL GetQtrMonBounds( PeriodEnd(4,Lake), 
     &            PeriodEnd(3,Lake), Q(1,Lake))
         ELSE 
            NewMonth = .FALSE.
         ENDIF
C
C        Set the date for end of quarter-month
         PeriodEnd(1,Lake) = Q( PeriodEnd(2,Lake) ,Lake)
C
C        Set NewMonth if first time 
         IF(PeriodEnd(5,Lake) .EQ. 0) NewMonth = .TRUE.
C
C        Adjust days elapsed to end of period
         PeriodEnd(5,Lake) = PeriodEnd(5,Lake) +
     &      (PeriodEnd(1,Lake) - Day(Lake) + 1)
C        Adjust days elapsed to end of next period.  Put length of next
C        period in I
         IF( PeriodEnd(2,Lake) .EQ. 4)THEN
C           Next quarter-month will be in new month.  First quarter of 
C           February is 7 days long; all other months first quarter is 
C           8 days long.  So see if next month will be Feb:
            IF( PeriodEnd(3,Lake) .EQ. 1) THEN
               I = 7
            ELSE
               I = 8
            ENDIF
         ELSE
C           Next quarter-month in same month, just subtract dates in Q
            I= Q(PeriodEnd(2,Lake) +1, Lake) - Q(PeriodEnd(2,Lake),Lake)
         ENDIF
         PeriodEnd(6,Lake) = PeriodEnd(5,Lake) + I
C
C
      ELSE IF( RoutingInterval(Lake) .EQ. 4 ) THEN
C        Monthly routing
         NewWeek = .TRUE.
         NewQMonth = .TRUE.
         NewMonth = .TRUE.
         PeriodEnd(2,Lake) = 4
         PeriodEnd(3,Lake) = PeriodEnd(3,Lake) + 1
C        Check to see if new year
         IF(PeriodEnd(3,Lake) .GT. 12 ) THEN
            PeriodEnd(3,Lake) = 1
            PeriodEnd(4,Lake) = PeriodEnd(4,Lake) + 1
         ENDIF
C        Update Q array for the lake
         CALL GetQtrMonBounds( PeriodEnd(4,Lake), 
     &            PeriodEnd(3,Lake), Q(1,Lake))
C        Set the date for end of month
         PeriodEnd(1,Lake) = Q( 4, Lake)
         PeriodEnd(5,Lake) = PeriodEnd(5,Lake) + 
     &      (PeriodEnd(1,Lake) - Day(Lake) + 1)
C
         IF(PeriodEnd(3,Lake) .EQ. 12 ) THEN
            PeriodEnd(6,Lake) = PeriodEnd(6,Lake) + 31
         ELSE         
C           Don't need number of days in year; use I as scratch
            CALL NUMDAYS(PeriodEnd(4,Lake), PeriodEnd(3,Lake) + 1,
     &         I, PeriodEnd(6,Lake))
            PeriodEnd(6,Lake) = PeriodEnd(6,Lake) + PeriodEnd(5,Lake)
         ENDIF
      ENDIF
C 
C
      IF(MonthLength .EQ. 2 .AND. RoutingInterval(Lake) .GE. 3) THEN
C        Actual quarter-month or month length
         SecsInPeriod = (PeriodEnd(1,Lake) - Day(Lake) + 1) * 86400
C        Adjust number of increments if keeping dT at 1 hour
         IF( NumIncsVary(Lake) ) RoutingIncrement(Lake) = 
     &      SecsInPeriod / 3600
         dT(Lake)= DBLE(SecsInPeriod)/DBLE(RoutingIncrement(Lake))
         I = 20
      ELSE
         CWRK = ' '
         I = 1
      ENDIF
C
      IF(GenVerbosity .GE. 5) THEN
         WRITE(CTMP,8966)  NameOfLake(Lake), CWRK(1:I), Day(Lake), 
     &      PeriodEnd(1,Lake), MonthAbbrev(PeriodEnd(3,Lake)), 
     &      PeriodEnd(4,Lake), PeriodEnd(5,Lake), NewWeek, NewQMonth, 
     &      NewMonth, LLT
         CALL LogIt(LOGD,1,1)
         WRITE(CTMP,9863) RoutingIncrement(Lake), dT(Lake), LLT
         CALL LogIt(LOGD,1,1)
      ENDIF
 9863 FORMAT(' Routing interval to have ', I5,
     &   ' increments, each of length dT= ', F12.4, ' secs', A1)
 8966 FORMAT(' Incrementing ', A, A, 'routing period:  start=',
     &    I2, ' end=', I2, A3, I5, ', elapsed days=', I10, 3(1X,L1), A1)
C
C           XX  Reset PeriodEnd(5,?) if GE MaxInt
      END
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------



