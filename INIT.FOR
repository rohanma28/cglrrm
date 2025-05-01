C-----------------------------------------------------------------------
C      
C     Steps for adding a new input parameter
C        1.  At first, hardcode the value for your new variable
C            with assignments or DATA statements.  Organize your 
C            variable and COMMON blocks according to what makes sense 
C            regardless of where the actual values come from.
C        2.  When your new variable is in place and performing its
C            function (in a very limited way), add it to the BLOCK DATA
C            GENRAL section at the top of INIT.FOR.  You can put it in
C            a new COMMON, or add it to an existing one.  If you add it
C            to an existing one, remember to update the declarations 
C            anywhere else that COMMON appears.  You should provide a 
C            short description of what the variable is for, and assign a
C            an obviously bogus value in a DATA statement.
C        3.  Make up a key phrase (16-char or less) and add it into an
C            empty space in the BLOCK DATA PLINIT section of INIT.FOR.
C            Empty slots contain a string of '?' characters.  Overtype
C            the '?' with the uppercase letters and numbers in your new
C            key phrase, and pad any extra characters with spaces.  Note
C            the number on the right, beyond column 72.  This is the 
C            order number of your parameter.
C        4.  Add the COMMON and declare the variable in subroutine 
C            ParLine() of INIT.FOR.  Then add 4000 to the previously 
C            mentioned order number and go to that label (i.e., if you
C            had parameter 14, look for label 4014).  This is where you
C            handle the input parsing, and assign the user-specified 
C            value to your variable.  There are functions established
C            to save you some work (for example, ReadIntegerValue), but
C            the details are too cumbersome to explain here.  The best
C            advice is to copy an existing section that parses the same
C            type of data (string, int array, float, filename, etc) and 
C            adapt it to your needs.  Some basics:
C               - the line being parsed is in PLINE
C               - the key phrase (left of the ':' char) is in PARCODE
C               - FKC and LKC refer to the position in PLINE of the 
C                 first and last 'keyboard' characters, excluding 
C                 whitespace, that come after the ':'. So, the value
C                 being parsed is in PLINE(FKC:LKC).
C               - If parsing is successful, jump to label 5000 (make
C                 sure EMessNo is less than 1000).
C               - If have an error during parsing, add a section to
C                 MESSBASE.TXT, assign the message number to
C                 EMessNo, and RETURN.  Or call ErrMess() and Adios()
C                 directly.
C               - Make your peace with LogIt() and use it.
C        5.  Add the COMMON and declare the variable in the CheckPar()
C            routine at the top of MISCUTL.FOR.  Find a good section
C            within the routine to do some checks on your variable.
C            At minimum, look to see if it still set at the bogus 
C            value and squawk about not specifying a value in the input.
C            Preferably, you would also make some value judgements about
C            whatever value was read in, within the context of the other
C            input.  Look at some of the neighboring checks for examples.
C        6.  Update InputSummary() to write the new parameter back out 
C            to the summary output file.  
C
C
      BLOCK DATA GENRAL
C
      IMPLICIT NONE
C
C     Just initializing some constants and strings used throughout the code.
C
C
      COMMON /PARMRD/   AlreadySpecified
      COMMON /ONTOUT/   OutputHeader
      COMMON /OUTSTF/   WriteTSOutput, WriteErieOutflow
      COMMON /FSPACE/   FWRK
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /OUTNAM/   RCOUTDIR, RCOUTEXT
      COMMON /OUTNMI/   NOPDIR,NOPEXT, TSOutputFileNum, OnlyOneOpen
      COMMON /MONAME/   MonthAbbrev
      COMMON /GENRLI/   OSTYPE, NUMPARAM
      COMMON /GENRLC/   OSName,DirSep,MessDB, ModelerName
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /RCRUNC/   RCRUNTYP, IntervalName, MonthLengthName,
     &                  MLSolutionName, SMSolutionName,
     &                  USLSolutionName
      COMMON /RCRUNI/   LakesToSolve, 
     &                  IntervalNameLen, MonthLength,
     &                  StartYear, StartMonth, StartDay, 
     &                  EndYear, EndMonth, EndDay
      COMMON /RCAVSZ/   NUMAVSZ, AVSZ
      COMMON /UNITSC/   UnitType
      COMMON /UNITSI/   UnitTypeLen, BaseUnitType, UnitConvFactor
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
      COMMON /PARSEC/   INDENT, PARCODE
      COMMON /LNAMEC/   LAKENAME, RiverName
      COMMON /LNAMEI/   NLN, LkToSlv2LN, NRN, LkSlv2RtLoopI
      COMMON /MISINT/   MidLakeSolution, BOMRound, FlowRound, 
     &                  LevelRound, RoutingInterval, RoutingIncrement,
     &                  StMarysSolution, UStLawrSolution, P77ACompat
      COMMON /EXTRAP/   LargestExtrapolation, Extrapolations 
      COMMON /STAGE/    MHK,SCK,ERK,MHYM,SCYM,ERYM,MHA,SCA,ERA,
     &                  MHB,SCB,ERB,MHWT,SCWT,ERWT,MHC,SCC,ERC
      COMMON /SPDAT/    SPK, SPYM, SPA, SPB, SPWT, SPC,
     &                  SGX, SGY, SGZ, MaxSGSetting,
     &                  SupGatesOpen
      COMMON /INITLV/   SPSTLV, MHSTLV, SCSTLV, ERSTLV, CGSTLV, ONSTLV,
     &                  PrevStMarysFlow, PrevSupGates, 
     &                  SPBOMLev, MHBOMLev, SCBOMLev, ERBOMLev  
      COMMON /FNAMEC/   FileNames, RepeatLines
      COMMON /FNAMEI/   RepeatIndex, NPFILZ
      COMMON /TITLES/   RCTITLE, MLHEADER
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
      COMMON /SUPOUT/   SupPreProject, BalancingFlow, SupGateFlow,
     &                  SupPlanFlow, SupBOMLevel, SWPierLevel, 
     &                  SupSideChanFlow, SPRegCodes, P77AHeader,
     &                  P2012Header
      COMMON /SUPCOM/   SPOPNBS06, OPLLO, MHOPNBS06, OPChicago, SCRRet,
     &                  SCOPNBS06, DetRet, EROPNBS06, OPWelland, NiaRet,
     &                  SPOPNBS07, MHOPNBS07, SCOPNBS07, EROPNBS07,
     &                  SPOPNBS08, MHOPNBS08, SCOPNBS08, EROPNBS08,
     &                  OPWelland_P2012, SCRRet_P2012, DetRet_P2012,
     &                  NiaRet_P2012,
     &                  ForecastMonths,
     &                  LStoSWP, SupCritRnd, BackCompat, MonthsChecked
      COMMON /SUPCRT/   CritBOMRound, CritFlowRound, CritLevelRound,
     &                  Q1887Method, CritBMethod, CritCThresholdEl
      COMMON /SWPIER/   SWPGateEqA, SWPGateEqB, OldSWPA, OldSWPB, 
     &                  OldSWPC, EC2000SWPA, EC2000SWPB, EC2000SWPC,
     &                  EC2008SWPA, EC2008SWPB, EC2008SWPC
      COMMON /USSLIP/   USSIceWeedAdj, USSlipEquation
      COMMON /SPOTOP/   SideChanMaxCapacity, CWOvertopElev, 
     &                  CWMinimumFreeboard, SideChanMinCapacity
      COMMON /PHYCAP/   numMaxSC, numMinSCWinter, numMinSCSummer,
     &                  numJacSC
      COMMON /INI58D/   CritDepths, P58DTotalSup, MaxILimit, MinPLimit,
     &                  P58DWeightedSupply, P58DRuleCurve, P58DAppLimit,
     &                  P58DRunType, P58DForecastYear, P58DForecastMon,
     &                  P58DForecastDay, SupplyIndicator, 
     &                  P58DPlanLevel, P58DActualLevel,
     &                  P58DSimulateFile, P58DSIDataFile, ElimDev, 
     &                  ForceFlow, SetFlow, PreProjConditions, 
     &                  P58DUseDatabase, P58DSiDatabase
      common /Ini58dc/  P58DAppLimitSym
      COMMON /INI58D1/  P58DAdjustment,P58DrcAdjustment,P58DActualFlow, 
     &                  P58DTotalDev,P58DLSLflow
      COMMON /INI58D2/  ComputeMontreal
      common /PreCond/  PreprojectLevel
      COMMON /CritK/    CritKfile

      COMMON /P2012/    P2012PARA, P2012SPTARGET, P2012MHTARGET, 
     &                  P2012SPHIGHLEVEL, P2012SPLOWLEVEL, 
     &                  P2012MHTARGETADJ 
      COMMON /STURGEON/ COUNTER_JUNEFLOW, Last_JuneFlow 
      COMMON /FORGATE/ ForcedSupGates
      COMMON /FORSIDE/ ForcedSideChan
      COMMON /INFORCE/ IndexGatesSide  
C
      CHARACTER*256 CritKfile
C
C     Seaway Alert Depths:
C   
      INTEGER      CritDepths
C                  0 = none, 1= alert, 2 = minimum
C
      INTEGER      P58DTotalSup
C
      INTEGER      MaxILimit
C                  623 (ie 220) or 680 (ie 240)
C
      INTEGER      MinPLimit
C                  0 = computed, 629 (ie 222) for all MINP
C
      INTEGER      P58DRunType
C                  0 = simulate only, 1= forecast only, 2 = both
C
      INTEGER      P58DForecastYear, P58DForecastMon, P58DForecastDay
C                  Required only if P58DRunType == 2
C
      CHARACTER*256 P58DSimulateFile
C                  Input file for SIMULATION (required only if P58DRunType==1) 
C
      CHARACTER*256 P58DSIDataFile
C                  database of supply indicators
C
      INTEGER      P58DWeightedSupply
C
      INTEGER      P58DRuleCurve
C
      INTEGER      P58DAppLimit
C
      CHARACTER*2  P58DAppLimitSym
C
      INTEGER      P58DAdjustment, P58DrcAdjustment, P58DActualFlow, 
     &             P58DTotalDev, P58DLSLflow
C
      INTEGER      SupplyIndicator(21)
C
      LOGICAL      ElimDev, ComputeMontreal
c
c     Should use the database for initial starting conditions
      logical      P58DUseDatabase, P58DSiDatabase
C
C                  set plan level = actual level
      LOGICAL      ForceFlow
C                  Ignore above levels for "forced" flows below
C
      LOGICAL      SetFlow
C                  Set flow = plan flow if total deviations < 0
C
      LOGICAL      PreProjConditions
C                  Apply Prep Conditions if 74.15 (242.77) > level > 75.37 (246.77)
C
      DOUBLE PRECISION P58DActualLevel, P58DPlanLevel
c
c     Preproject variables
c
      DOUBLE PRECISION PreprojectLevel
C                  
C                  
      DOUBLE PRECISION P2012SPTARGET(12), P2012MHTARGET(12),
     &                 P2012SPHIGHLEVEL, P2012SPLOWLEVEL,
     &                 P2012MHTARGETADJ   
      DOUBLE PRECISION P2012PARA(6) 
      DOUBLE PRECISION ForcedSupGates(12)
      INTEGER          ForcedSideChan(12)
      INTEGER          COUNTER_JUNEFLOW, IndexGatesSide 
C                      Counter since last June sturgeon spawning flow
      INTEGER          Last_JuneFlow
C                      Year of last June Sturgeon flow prior to start
C                      of simulation.  Currently (June 2013), a "sturgeon
C                      flow" is considered to be a June St. Marys flow of
C                      >~1700cms.  The year is compared to the starting
C                      year of the simulation (in XXX) to initialize 
C                      COUNTER_JUNEFLOW
      INTEGER          Extrapolations
      DOUBLE PRECISION LargestExtrapolation
C                  Variables used to describe inadequate area-elev data
C
      DOUBLE PRECISION FWRK(2700)
C                  This is a common block used for 64-bit calculations.
C                  This usage is subject to change pending better ideas.
C                  Idea is to re-use a large common block for doing
C                  double precision arithmetic, rather than allocating
C                  new variables all the time, since nothing ever seems
C                  to go out of scope, and can't allocate/deallocate on demand.
C
      CHARACTER*80 RCTITLE(3), MLHEADER
C                  These are descriptive lines for use in output.
C
      CHARACTER*3  MonthAbbrev(12)
C                  a 12 element char array of 3 chars per element.  Holds the
C                  abbreviated name of each month.  Used in various outputs. 
C
      CHARACTER*64 ModelerName
C                  Optional description for the person(?) running the model.
C                  Ends up in headers of the output files.
C
      CHARACTER*256 MessDB
C                  char array holding the the filespec for the message 
C                  database.  Initialized to default of "messbase.txt"
C                  in the pwd, but over-ridden by the parameter of the 
C                  same name.
C
      CHARACTER*128 RCOUTDIR
C                  char array holding the path to the directory in which
C                  to write the output files.  If not specified, the 
C                  present working directory is used.  Warning: be sure
C                  that the value for RCOUTDIR ends in the appropriate
C                  directory separator (ie, '\' for DOS, '/' for Unix,
C                  '*" for many mainframe OS's), because the value for 
C                  RCOUTDIR is simply pre-pended to the output filenames.
C                  Also be aware of the present working directory at run-
C                  time, if RCOUTDIR specifies a relative path.
C
      CHARACTER*16 RCOUTEXT
C                  char array holding the 1-16 Character
C                  filename extension to apply to output files.  If not
C                  specified, no extension is used.
C
C      INTEGER      OutFileHandles(5)
C                  Unit specifiers for summary output files, corresponding
C                  to OutputFileNames(5) in RunIt().  
C
      CHARACTER*8  OSName(5)
C                  a 5 element char array of 8 chars per element.  Holds 
C                  strings describing operating systems that might be
C                  encountered. This value is used as extra detail in output.
C
      CHARACTER*5  DirSep
C                  a 5 element char array listing the recognized directory
C                  separators corresponding to the operating systems listed 
C                  in OSName.  It is determined as the last character in 
C                  RCOUTDIR.  This value is used only to infer OSName.
C
      LOGICAL      LakesToSolve(8)
C                  an 8 element logical array listing which lakes to solve for.
C                     1          Lake Superior
C                     2   Lakes Michigan-Huron
C                     3         Lake St. Clair
C                     4              Lake Erie
C                     5           Lake Ontario
C                     6      Lake St. Lawrence
C                     7               (unused)
C                     8               (unused)
C
      CHARACTER*8  IntervalName(4)
C                  an 8 element char array containign names of the routing 
C                  intervals supported by the model, for use in output.
C
      INTEGER      IntervalNameLen(4)
C                  The number of letters in each element of IntervalName
C
      CHARACTER*8  RCRUNTYP
C                  an 8 element char array containing a code which activates
C                  a pre-determined set of parameter values.   (NOT USED)
C
      CHARACTER    PARCODE*16, INDENT*4
C                  See descriptions in PARLINE().
C
      LOGICAL      AlreadySpecified(420)
C                  Indicates whether a parameter has already been read in
C
      CHARACTER*21 LAKENAME(16)
C                  Array containing long names for each lake, used in output.
C                  The first eight names must correspond to the first eight
C                  entries in AVSZ and NUMAVSZ.
C
      CHARACTER*21 RiverName(16)
C                  Array containing long names for each river, used in output.
C                  The names correspond to the outlets of the lakes in LakeName.
C
      INTEGER      LkToSlv2LN(8)
C                  Gives element number for LakeName corresponding to each
C                  member of LakesToSolve 
C
      INTEGER      LkSlv2RtLoopI(8)
C                  Lists which routing loop to use for each member
C                  of LakesToSolve. 
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
      CHARACTER*256 FileNames(100)
C                  This is the array containing all the input filenames.
C                  The number of elements (100) and lengths of the array
C                  elements (80) are arbitrary, but FORTRAN complies poorly
C                  with 0-1-infinity rule, and this makes things simpler.
C                  It's a horribly expensive way to go about handling this,
C                  and will need to replace it with some clever scheme of 
C                  dodging FORTRAN 77's limitations when we eventually bump
C                  up against memory constraints for DOS implementations.  
C                  This array is also used for repeating time-series data.
C                  Instead of a file name, the element contains the string
C                  'Repeating Data'//LLT in the first 15 characters, 
C                  followed by a integer in the three characters, which
C                  tells how many lines of RepeatLines are used for this
C                  data.  The remainder of the FileNames entry for repeating
C                  data is a list of numbers in I3, each in format, which
C                  correspond to the elements in RepeatLines.
C
C                  The array entries are as follows:
C
C                     Lake Superior
C                         FileNames(42) = SUPNBSDATA
C                         FileNames(43) = SUPCONUSEDATA
C                         FileNames(44) = LONGLACDIVDATA
C                         FileNames(45) = OGOKIDIVDATA 
C                         FileNames(46) = COMBLLODIVDATA 
C                         FileNames(47) = SUPOTHERDIVDATA
C                         FileNames(48) = SUPOTHERSUPDATA
C                         FileNames(49) = SUPDEVIATIONDATA
C                         FileNames(40) = SUPGATESDATA
C                         FileNames(93) = SUPINDEX1DATA
C                         FileNames(94) = SUPINDEX2DATA
C                         FileNames(95) = SUPINDEX3DATA
C                         FileNames(68) = SUPOUTFLOWDATA
C                         FileNames(69) = SUPSIDECHANDATA 
C                         FileNames(70) = SUPEOPLEVELDATA 
C                         FileNames(41) = SUPSUMINPUTFILE 
C
C                     St. Marys River
C                         FileNames(19) = STMRIVFLOWDATA
C                         FileNames(50) = STMRIVRETRDATA
C                         FileNames(51) = STMRIVICEDATA
C                         FileNames(52) = STMRIVWEEDDATA
C
C                     Lakes Michigan-Huron
C                         FileNames(1)  = MHUNBSDATA
C                         FileNames(23) = MHUCONUSEDATA
C                         FileNames(24) = MHUOTHERDIVDATA
C                         FileNames(25) = MHUOTHERSUPDATA
C                         FileNames(4)  = MICCONUSEDATA
C                         FileNames(12) = MICOTHERDIVDATA
C                         FileNames(8)  = MICOTHERSUPDATA
C                         FileNames(5)  = HURCONUSEDATA
C                         FileNames(13) = HUROTHERDIVDATA
C                         FileNames(9)  = HUROTHERSUPDATA 
C                         FileNames(22) = CHICAGODIVDATA
C                         FileNames(71) = MHUEOPLEVELDATA 
C                         FileNames(31) = MHUSUMINPUTFILE
C                         FileNames(26) = USSIceWeedAdj 
C
C                     St. Clair River
C                         FileNames(28) = STCRIVFLOWDATA  
C                         FileNames(16) = STCRIVRETRDATA
C                         FileNames(34) = STCRIVICEDATA
C                         FileNames(35) = STCRIVWEEDDATA
C
C                     Lake St. Clair
C                         FileNames(2)  = STCNBSDATA      
C                         FileNames(6)  = STCCONUSEDATA   
C                         FileNames(14) = STCOTHERDIVDATA
C                         FileNames(10) = STCOTHERSUPDATA
C                         FileNames(72) = STCEOPLEVELDATA 
C                         FileNames(32) = STCSUMINPUTFILE 
C
C                     Detroit River
C                         FileNames(29) = DETRIVFLOWDATA  
C                         FileNames(17) = DETRIVRETRDATA
C                         FileNames(36) = DETRIVICEDATA  
C                         FileNames(37) = DETRIVWEEDDATA  
C
C                     Lake Erie
C                         FileNames(3)  = ERINBSDATA
C                         FileNames(7)  = ERICONUSEDATA
C                         FileNames(15) = ERIOTHERDIVDATA
C                         FileNames(11) = ERIOTHERSUPDATA
C                         FileNames(21) = WELLANDCANALDATA
C                         FileNames(73) = ERIEOPLEVELDATA 
C                         FileNames(33) = ERISUMINPUTFILE 
C
C                     Niagara River
C                         FileNames(30) = NIARIVFLOWDATA  
C                         FileNames(18) = NIARIVRETRDATA
C                         FileNames(38) = NIARIVICEDATA  
C                         FileNames(39) = NIARIVWEEDDATA  
C                         FileNames(20) = CGIPAVGLEVELDATA
C
C                     Lake Ontario
C                         FileNames(53) = ONTNBSDATA      
C                         FileNames(54) = ONTCONUSEDATA   
C                         FileNames(55) = ONTOTHERDIVDATA
C                         FileNames(56) = ONTOTHERSUPDATA
C                         FileNames(67) = ONTGATESDATA 
C                         FileNames(74) = ONTEOPLEVELDATA 
C                         FileNames(57) = ONTSUMINPUTFILE 
C
C                     St. Lawrence River
C                         FileNames(58) = STLRIVFLOWDATA  
C                         FileNames(62) = LSLRETRDATA  
C                         FileNames(66) = LSLFLOWDATA  
C                         FileNames(63) = LSROUGHNESSDATA  
C                         FileNames(64) = ONTDEVIATIONDATA
C                         FileNames(21) = WELLANDCANALDATA
C                         FileNames(65) = DPRRIVFLOWDATA  
C                         FileNames(77) = MIIRIVFLOWDATA  
C                         FileNames(59) = STLRIVPPRETRDATA
C                         FileNames(60) = STLRIVPPICEDATA  
C                         FileNames(61) = STLRIVPPWEEDDATA  
C                         FileNames(75) = ONTFORCEFLOWDATA  
C                         FileNames(76) = LSLONTDIFFDATA
C                         FileNames(78) = FLOWDPMIDATA
C                         FileNames(79) = FLOWSTFRANDATA
C                         FileNames(80) = FLOWSTMAURDATA
C                         FileNames(81) = FLOWRICHDATA
C                         FileNames(82) = ROUGHJETTY1DATA
C                         FileNames(83) = ROOUGHVARENDATA
C                         FileNames(84) = ROUGHSORELDATA
C                         FileNames(85) = ROUGHPIERREDATA
C                         FileNames(86) = ROUGH3RIVERSDATA
C                         FileNames(87) = ROUGHBATISDATA
C                         FileNames(88) = TIDALDATA
C                         FileNames(89) = ONTTOTALSUPDATA
C                         FileNames(90) = OTHER06ROUGHDATA
C                         FileNames(91) = OTHER06FLOWDATA
C 
C
C                     Unassigned
C                         FileNames(27) = 
C                         FileNames(92) =  
C                         FileNames(96) =  
C                         FileNames(97) =  
C                         FileNames(98) =  
C                         FileNames(99) =  
C                         FileNames(100) =  
C
      CHARACTER*256 RepeatLines(80)
C                  Stores lines of repeating time-series data (as opposed to
C                  time-series data being read from a file).  
C
      INTEGER      RepeatIndex
C                  Number of elements in RepeatLines
C 
      INTEGER      NLN(16), NRN(16)
C                  Array containing length of each name in LAKENAME. The
C                  first eight names must correspond to the first eight
C                  entries in AVSZ and NUMAVSZ.  Also number of characters 
C                  in RiverName.
C
      INTEGER      BOMRound 
C                  Specifies how to round the End-of-Month levels for 
C                  the middle lakes.  Typically set to -2 (for cm) in 
C                  interest of maintaining compatibility with traditional
C                  practice.  Matt feels strongly that there is no longer
C                  any reason to throw away the digits, and that it should
C                  be set to -99 (for no rounding).
C
      INTEGER      FlowRound, LevelRound
C                  These values determine how to round the flows and levels 
C                  during the internal routing and regulation computations.
C                  They correspond to powers of ten; for instance, a value
C                  of 3 means round to the thousand, while a value of -2 means
C                  round to the hundredth.  The allowable values are integers
C                  from -4 to 3.  Typically set to -4 (for 1/10 mm) for 
C                  levels, and 0 (for nearest m3s) for flows, in interest of 
C                  maintaining compatibility with existing models.  Matt
C                  feels strongly that there is no longer any reason to throw 
C                  away the digits, and that it should be set to -99 (for no 
C                  rounding).
C
      INTEGER      StMarysSolution
C                  How to compute Superior outflows. 1=outflow eqn, 2=1955 mod
C                  rule, 3=plan77a
C
      INTEGER      P77ACompat
C                  -1 = virtcarl supplies/diversions/retardations
C                   0 = existing supplies/diversions/retardations
C                   1 = operational supplies/diversions/retardations
C                   2 = 6-month reg run use of opp supply matrix
C
      INTEGER      MonthsChecked
C                  which month of OpSupply sequence we are in if p77Acompat=2.
C                    
      LOGICAL      BackCompat
C                  .TRUE. if program is trying to reproduce results from
C                  older upper lakes models.
C
      INTEGER      SupCritRnd
C                   1 = use P77A rounding for criterion calcs
C                   2 = use general settings for criterion calcs
C
      CHARACTER*16 SMSolutionName(4)
C                  Strings set to 'Outflow Equation', '1955 Mod. Rule',
C                  and 'Plan 1977A', etc., used for output.
C
      DOUBLE PRECISION SGX(17), SGY(17), SGZ(17)
C                  SGX, SGY, SGZ are constants for the Pt Iroquois gate equations.  
C
      INTEGER      MaxSGSetting
C                  The number of elements in SGX, SGY, and SGZ
C
      INTEGER      ForecastMonths
C                  The number of months in the future for which expected Lake 
C                  Superior supplies are provided or assumed by the regulation
C                  plan.  Should include diversions at Long Lac and Ogoki. 
C
      INTEGER      LStoSWP
C                  Which Lake Superior:SW Pier relationship to use for 
C                  calculating gated flow.  Used to determine sidechannel (i.e.
C                  non-gated flow in SupCriterions(), and to compute gated
C                  flow for routing in SuperiorRoute.
C
      DOUBLE PRECISION OldSWPA(12), OldSWPB(12), OldSWPC
      DOUBLE PRECISION EC2000SWPA, EC2000SWPB, EC2000SWPC(12)
      DOUBLE PRECISION EC2008SWPA, EC2008SWPB, EC2008SWPC(12)
C                  Constants for SW Pier/Pt Iroquois relationship
C
      DOUBLE PRECISION SWPGateEqA(16), SWPGateEqB(16)  
C                  Constants for SW Pier gate equations, used for gated
C                  flow calculation.
C
      INTEGER      USSlipEquation
C                  1=ancient eq (default), 2=2010 EC eq
      DOUBLE PRECISION USSIceWeedAdj             
C
      INTEGER      SideChanMaxCapacity(12), SideChanMinCapacity(12)
      DOUBLE PRECISION CWOvertopElev, CWMinimumFreeboard(12)
C                  See explanations in RouteSuperior() in Superior.for
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
C
      INTEGER      UStLawrSolution
C                  How to compute Ontario outflows. 1=outflow eqn, 2=1958D,
C                  3=IS Plan, 4=Plan 1998, 5=Preproject
C
      CHARACTER*16 USLSolutionName(5)
C                  Strings set to 'Outflow Equation', 'Plan 1958D', 'IS Plan',
C                  'Plan 1998' and 'Preproject'; used for output.
C
      INTEGER      MidLakeSolution
C                  The method used to solve the simulataneous equations for 
C                  routing through the middle lakes.  A value of one means
C                  interative solution, as described in the Plan 77A documents;
C                  a value of 2 means the matrix solution, as described in 
C                  GLERL's techmemo.
C
      CHARACTER*9  MLSolutionName(2)
C                  Strings set to 'interative' and 'matrix', used for output.
C
      INTEGER      RoutingInterval(3), RoutingIncrement(3)
C                  Routing interval to use on Superior, middle lakes,
C                  and Ontario (1=daily, 2=weekly, 3=qmonthly, 
C                  4=monthly), and the number of increments within each 
C                  interval.
C
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C                  8 element logcial arrays indicating different 
C                  combinations of log files to write messages. 
C
      INTEGER      StartYear, StartMonth, StartDay, 
     &             EndYear, EndMonth, EndDay
C                  Starting and ending year, month, day.  Read in as three
C                  sets of comma-delimited integers.
C                  EndYear must be >= StartYear, etc.
C
C----&------------------------------------------------------------------
      INTEGER      MonthLength
      CHARACTER*8  MonthLengthName(2)
C                  Used only for quarter-monthly and monthly routing.
C                  MonthLength is 1 for the traditional "uniform" month
C                  and quarter-month lengths, or 2 for the actual length
C                  of the period.  For actual quarter-month length, the 
C                  period is defined by the number of whole days assigned
C                  to each quarter-month, as defined in SUBROUTINE
C                  GetQtrMonBounds.  The uniform month length is taken
C                  as the average month length according to the Gregorian
C                  calendar, which calls for 97 leap days every 400 years.
C                  So the mean length of a month works out to:
C                    (400 * 365 + 97) days * (24 * 60 * 60) seconds/day
C                    --------------------------------------------------
C                                  (400 * 12) months    
C                    = 2629746 seconds/month 
C                  The length of a quarter-month is computed from that.
C                  MonthLengthName is two-element character array used
C                  for describing the value of MonthLength.  It is holds
C                  the strings "uniform" and "actual" (each padded with
C                  a trailing space for use with function EndStr).
C                  
      INTEGER      OSTYPE
C                  the index to the OSNAME array. 
C
      INTEGER      NPFILZ
C                  The number of parameter files read.  It is used for 
C                  detailed output.
C
      INTEGER      NUMPARAM
C                  the number of valid parameters read in. 
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C
C                  The degree of debugging information written out to the
C                  detailed log file.  A value of 1 is most terse, while
C                  9 means most verbose.  The default is 5.  The value
C                  used in any given module defaults to GenVerbosity,
C                  unless overridden by the appropriate one of TSDVerbosity,  
C                  ChkVerbosity, SupVerbosity, MLVerbosity, or OntVerbosity.
C
C                  These values control the output level during execution.
C                  GenVerbosity is a global value which sets the minimum 
C                  output for anywhere in the program. 
C                  The other values override the global value for output
C                  within specific program areas:
C                     TSDVerbosity: sets message level regarding 
C                                   processing of time series data
C                     ChkVerbosity: sets message level regarding 
C                                   the consistency check module
C                     SupVerbosity: sets message level regarding 
C                                   the Lake Superior modules
C                     MLVerbosity:  sets message level regarding 
C                                   the middle lakes modules
C                     OntVerbosity: sets message level regarding 
C                                   the Lake Ontario modules
C                  Note that only the GenVerbosity control is available 
C                  while input parameters are being parsed.  This value 
C                  may be set and re-set within the parameter files to 
C                  provide debugging info as desired, but the LAST value
C                  specified will remain in effect for the duration of the
C                  computations.
C
      INTEGER      NOPDIR, NOPEXT
C                  the length of variables RCOUTDIR and RCOUTEXT, (not 
C                  including trailing whitespace).  They are used in the 
C                  construction of output filenames.  
C
      LOGICAL      OnlyOneOpen
C                  .TRUE. if the time-series output files are to be 
C                  opened prior to EACH write, and closed afterward.  
C                  .FALSE. if files are to be left open throughout 
C                  routing calculations.
C
      INTEGER      TSOutputFileNum(4,2,8)
C                  Holds the filespecs to use for each of the potential
C                  time-series output files, such as Lake Erie monthly
C                  mean level, Lake Superior daily outflow, etc.  The
C                  first index (1-4) tells whether file is for daily,  
C                  weekly, quarter-monthly, or monthly data.  The second
C                  index in either level (1) or flow (2).  The third 
C                  index is for the lake (numbers correspond to names in
C                  LakesToSolve().
C
      LOGICAL      WriteTSOutput(4,2,8)
C                  For each lake (designated by third index - 8), for
C                  each type of time-series output (flow or ave level),
C                  for each interval (d, w, q, or m), .TRUE. if should
C                  write out the data to a time-series output file. 
C       
      LOGICAL      WriteErieOutflow
C                  .TRUE. if the time-series output files for Niagara River
C                  should be replaced by one for Lake Erie Outflow.  The 
C                  output file name will begin with 'er' instead of 'nr'.
C                  The values will be the Niagara flow plus the Welland flow.
C
      INTEGER      NUMAVSZ(8)
C                  The number of area/elevation points read in for each
C                  lake.  A zero means no such points were read.
C
      DOUBLE PRECISION AVSZ(8,5,2)
C                  This is the three-dimensional array containing the
C                  area-elevation curve for each lake.  The first dimension
C                  is the lake index, corresponding to the 8 names listed
C                  under NUMAVSZ.  The second dimension is the index to the
C                  area-elevation points held by the array; currently only
C                  5 elevation-area pairs are supported.  The third
C                  dimension refers to the individual value of the data pair
C                  ie, either the area or the elvation.  The elevations are
C                  ALWAYS expressed in meters IGLD-85, while the areas are
C                  ALWAYS expressed in  square meters.  The elevation
C                  may be '0' if keeping area constant.  For instance,
C                  values specified for various array members might be:
C                   AVSZ(2,3,1) = elev for 3rd point for M-H = 176.50D6
C                   AVSZ(2,3,2) = area for 3rd point for M-H = 117400D6
C                   AVSZ(6,1,2) = area for 1st point for Erie = 25700D6
C                   AVSZ(5,1,1) = elev for 1st point for St. Clair = 0
C                   AVSZ(5,1,2) = area for 1st point for St. Clair = 1114D6
C
C                  This order must be the same as the first eight lakes in 
C                  LAKENAME and NLN.
C
      CHARACTER*6  UnitType(20)
C                  The list of recognized unit types.  Most numeric 
C                  parameters must be accompanied by their units.  
C                   
      INTEGER      UnitTypeLen(20)
C                  The length of each member of UnitType.  
C                   
      INTEGER      BaseUnitType(20)
C                  The base member of UnitType that each member must be
C                  converted to on being read in (ie., the fourth element
C                  of UnitType is tcfs, which must be converted before
C                  being used into m3s, which is the second element, so
C                  BaseUnitType(4) = 2).  
C
      DOUBLE PRECISION UnitConvFactor(20)
C                  The conversion factor used to convert each unit type in
C                  UnitType to its base member.  Initialized in INITPARMS().
C
      INTEGER      TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &             MissingVal
      DOUBLE PRECISION MPFT, HALF, ONE, ZERO
C                  These are constants for the value 10, the number of
C                  feet per mile, the number of inches per foot, the number
C                  of seconds in a day, the number of meters per foot, 0.5D0,
C                  1.0D0., and 0.0D0.  And the minimum and maximum integer
C                  values permitted (used for loop control).
C
      INTEGER      SecsInUniformMon
C                  Used for calculating dT, and converting storage to
C                  flow. The number of seconds in a uniform month,  
C                  2629746, is taken from the average month length 
C                  according to the Gregorian calendar, which calls 
C                  for 97 leap days every 400 years. So the mean 
C                  length of a month works out to:
C                  (400 * 365 + 97) days * (24 * 60 * 60) seconds/day
C                  --------------------------------------------------
C                                  (400 * 12) months    
C                   = 2629746 seconds/month 
C                  Uniform quarter-month length taken as 1/4 of this.
C
      DOUBLE PRECISION SPSTLV, MHSTLV, SCSTLV, ERSTLV, CGSTLV, ONSTLV
C                  The initial lake levels.
C
      DOUBLE PRECISION SPBOMLev, MHBOMLev, SCBOMLev, ERBOMLev
C                  The beginning-of-month lake levels used by Plan 1977A
C
      INTEGER      PrevStMarysFlow
C                  Initial Lake Superior outflow (average for month 0)
C 
      REAL         PrevSupGates 
C                  Initial gates open at Comp Works
C 
      REAL         SupGatesOpen
C                  The number of Superior gates open in current routing 
C                  period.
C
      DOUBLE PRECISION SPK,SPYM,SPA,SPB,SPWT,SPC
      DOUBLE PRECISION MHK,MHYM,MHA,MHB,MHWT,MHC
      DOUBLE PRECISION SCK,SCYM,SCA,SCB,SCWT,SCC
      DOUBLE PRECISION ERK,ERYM,ERA,ERB,ERWT,ERC
C                   See descriptions in CHECKPAR and OUTFLOW
C
      CHARACTER*40     SPRegCodes
C                      Describes which Plan 1977A Requirements and 
C                      Criteria come into play, used for output.
C
      CHARACTER        BalancingFlow*8, SupBOMLevel*8
C                      The flow called for by the balancing equation,
C                      and beginning-of-month level.  Used for output.
C
      CHARACTER*120    P77AHeader
C                      Describes columns of Superior output
C
      CHARACTER*120    P2012Header
C                      Describes columns of Lake Superior Output under Plan 2012
C
      CHARACTER        SupPreProject*6, SupGateFlow*6,
     &                 SupSideChanFlow*6, SupPlanFlow*6, SWPierLevel*8
C                      More Lake Superior regulation output strings
C
      CHARACTER*(80)   OutputHeader
C                      Used to identify columns in ontario output file.
C
      INTEGER          MHOPNBS06(12,6), MHOPNBS07(12,6), OPChicago,  
     &                 SCOPNBS06(12,6), SCOPNBS07(12,6), 
     &                 EROPNBS06(12,6), EROPNBS07(12,6), OPWelland(12), 
     &                 SCRRet(12), DetRet(12), NiaRet(12),
     &                 OPLLO(12), SPOPNBS06(12,6), SPOPNBS07(12,6),
     &                 SPOPNBS08(12,6), MHOPNBS08(12,6),
     &                 SCOPNBS08(12,6), EROPNBS08(12,6),
     &                 OPWelland_P2012(12), SCRRet_P2012(12),
     &                 DetRet_P2012(12), NiaRet_P2012(12)      
C                      The 50% NBS and related values, potentially used for 
C                      Criterion B and C checking.
C
      INTEGER          CritBOMRound, CritFlowRound, CritLevelRound
C                      Rounding values for how the criterions get calculated.
C                      These are currently set to the same as Plan 1977A,
C                      to replicate existing operational practice.  There
C                      is a block of commented code for setting them to 
C                      the general parameters, if practice is changed.    
C
      INTEGER          Q1887Method, CritBMethod
C                      These specify which variation in computational 
C                      methods to use for Criterions B and C.
C                      If Q1887Method is 1, use the old way (get pre-project
C                      by popping the Plan 1977A forecasted level into 
C                      the 1887 equation.  If Q1887Method is 2, the sub
C                      GetSupPreProject solves for the Criterion C 
C                      iteratively, based on the BOP level and water 
C                      balance terms.  Either way, the trigger remains
C                      the Plan 1977A forecasted level.  If CritBMethod is
C                      1, then the Criterion B flow is allowed to fall below
C                      the pre-project flow (old way).  If CritBMethod is 2,
C                      then the Criterion B flow is not allowed to be less
C                      than the pre-project flow.
C
      DOUBLE PRECISION CritCThresholdEl
C                      Elevation at which minimum flow becomes pre-project - 
C                      specified as 183.40 in the Orders of Approval.
C
      INTEGER          numMaxSC, numMinSCWinter, numMinSCSummer,
     &                 numJacSC
C                      Counters for number of times the physical capacity 
C                      overides the desired sidechannel flow
C
      INTEGER          I, J, K
C                      Local scratch variables used for indices.
C
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
C
C     The following section of DATA statements fills in constants
C     that the program knows.
C                  
      DATA MonthAbbrev /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',
     &             'Sep','Oct','Nov','Dec'/
C
      DATA DirSep /'??/*:'/
C
      DATA MonthLengthName /'uniform ', 'actual  '/
C
      DATA OSName /'unknown ','DOS     ','Unix    ','Mac     ',
     &   'VMS     '/
C
      DATA MLSolutionName /'iterative', 'matrix   '/
C
      DATA SMSolutionName /   'outflow equation', '1955 Mod. Rule  ', 
     &    'Plan 1977A      ', 'Plan 2012       '/
C
      DATA USLSolutionName /'outflow equation', 'Plan 1958D      ', 
     &    'IS Plan         ','Plan 1998       ','Preproject      '/
C
      DATA OutputHeader(1:60) /
     &   '  Period End     BOP    Mean     EOP Outflow                '/
      DATA OutputHeader(61:80) / '                    ' /
C
C     Be careful how these are ordered, so that no spurious matches
C     when searching list (ie, "123.879mm" would match unit 'M' if
C     unit 'MM' didn't precede it in list.
      DATA UnitType /'10M3S','M3S','TCFS','CFS','MM','CM','M','KM2',
     &   'M2','M3','INCHES','FEET','FT2','FT3','MI2','LT10','MISS',
     &   'NONE', 'NA', '  --  ' /
C
      DATA UnitTypeLen /5,3,4,3,2,2,1,3,2,2,6,4,3,3,3,4,4,4,2,6/
C
      DATA BaseUnitType /2,2,2,2,7,7,7,9,9,10,7,7,9,10,9,16,17,18,19,20/
C
      DATA TEN, FTPMI, INPFT, SecsInDay /10, 5280, 12, 86400/
C
      DATA SecsInUniformMon / 2629746 /
C
      DATA MPFT, HALF, ONE, ZERO /0.3048D0, 0.5D0, 1.0D0, 0.0D0/
C 
C     The ANSI F77 standard calls for INTEGER data type to be a
C     signed, 32-bit value, so valid range is from -2147483648 to
C     -2147483648, plus or minus one, depending on the implementation.
C     And some compilers can act a little flaky with integers at the
C     max/min value, so say max/min is (-)2147483647.  But for the
C     mechanics of routing, we need to have at least 31 days of 
C     elbow room for knowing which routing step the current month
C     will end on.  So we're down to 2147483616.  Finally, we want
C     to cap the maximum at a number that is divisible by 7, for 
C     purposes of weekly routing.  So if we roll counters over at 
C     2147483611, then the change occurs on a week boundary.
      DATA MinInt, MaxInt / -2147483611, 2147483611/
C
      DATA MissingVal / -9999 /
C
      DATA Extrapolations / 0 /
      DATA LargestExtrapolation / 0.0D0 /
C
      DATA INDENT /'===>'/
C
      DATA IntervalName /'daily','weekly','qmonthly','monthly'/
C
      DATA IntervalNameLen /5,6,8,7/
C
      DATA LAKENAME /
     1          'Lake Superior',
     2   'Lakes Michigan-Huron',
     3          'Lake Michigan',
     4             'Lake Huron',
     5         'Lake St. Clair',
     6              'Lake Erie',
     7           'Lake Ontario',
     8         'Lake Champlain',
     9      'Grass Island Pool',
     A      'Lake St. Lawrence',
     B         'Lake St. Louis',
     C     'Lake 1 by Montreal',
     D     'Lake 2 by Montreal',
     E     'Lake 3 by Montreal',
     F     'Lake 4 by Montreal',
     G  'Lake Future Expansion'/
C
      DATA NLN / 13,20,13,10,14, 9,12,14,17,17,14,18,18,18,18,21 /
C
      DATA LkToSlv2LN /1,2,5,6,7,0,0,0/
C
      DATA LkSlv2RtLoopI / 1, 2, 2, 2, 3, 3, 0, 0/
C
      DATA RiverName /
     1        'St. Marys River',
     2        'St. Clair River',
     3                      ' ',
     4                      ' ',
     5          'Detroit River',
     6          'Niagara River',
     7           'Lake Ontario',
     8     'St. Lawrence River',
     9                      ' ',
     A                      ' ',
     B                      ' ',
     C                      ' ',
     D                      ' ',
     E                      ' ',
     F                      ' ',
     G                      ' '/
C
      DATA NRN / 15,15,0,0,13,13,12,18,0,0,0,0,0,0,0,0 /
C
      DATA P77AHeader( 1:40)/'    Period End  B-of-M    Mean     EOP S'/
      DATA P77AHeader(41:80)/'W Pier Gates  Bal Eq   Plan PrePrj Actua'/
      DATA P77AHeader(81:)  /'l  Gated SideCh   Crit AdjOTFlow Codes'/
C                  ____________Levels_____________       _____________________Flows_______________
C     Period End  B-of-M    Mean     EOP SW Pier Gates  Bal Eq   Plan PrePrj Actual  Gated SideCh   Crit AdjOTFlow Codes  
C    1  746814 D 183.440 183.440 183.440 183.257  4.00    5704   2902   2377   2899    573   2326   5704   -9999.  NS       RI          
C
C
      DATA P2012Header(1:40)/'    Period End  B-of-M    Mean     EOP S'/
      DATA P2012Header(41:79)/'W Pier Gates  Sturgn   Plan PrePrj Actu'/
      DATA P2012Header(80:115)  /'al  Gated SideCh   Crit AdjOTFlow  C'/
      DATA P2012Header(116:)  /'BUSS Codes'/
C                  ____________Levels_____________       _____________________Flows_______________
C     Period End  B-of-M    Mean     EOP SW Pier Gates  Sturgn   Plan PrePrj Actual  Gated SideCh   Crit AdjOTFlow Codes  
C    1  746814 D 183.440 183.440 183.440 183.257  4.00       1   2902   2377   2899    573   2326   5704   -9999.  NS       RI          
C
C
      DATA (AlreadySpecified(I), I = 1,420) / 420 * .FALSE. /
C
C     I/O unit numbers for output files
      DATA TSOutputFileNum /
     &   10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
     &   20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 
     &   30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 
     &   40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 
     &   50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 
     &   60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 
     &   70, 71, 72, 73 /
C
C     Definitions for the various log file possibilities
C                  col 1 - stdout (unit 6)
C                  col 2 - detailed log file (unit 90)
C                  col 3 - critical failure log (unit 91)
C                          (failure.log deactivated as of 12/31/01)
C                  col 4 - summary output file (unit 92)
C                  col 5 - Middle Lakes output file (unit 93)
C                  col 6 - Lake Superior output file (unit 94)
      DATA LOGC  / .TRUE., 7 * .FALSE. /
      DATA LOGD  / .FALSE., .TRUE., 6 * .FALSE. /
      DATA LOGS  / 3 * .FALSE., .TRUE., 4 * .FALSE. /
      DATA LOGCD  / 2 * .TRUE., 6 * .FALSE. /
      DATA LOGSD  / .FALSE., .TRUE., .FALSE., .TRUE., 4 * .FALSE. /
C      DATA LOGCE  / .TRUE., .FALSE., .TRUE., 5 * .FALSE.  /
      DATA LOGCE  / .TRUE., .FALSE., 6 * .FALSE.  /
C      DATA LOGCDE  / 3 * .TRUE., .FALSE., 4 * .FALSE.  /
      DATA LOGCDE  / 2 * .TRUE., 6 * .FALSE.  /
C     Messages for the middle lakes summary output:
      DATA LOGM  / 4 * .FALSE., .TRUE., 3 * .FALSE.  /
C     Messages for the detailed log and middle lakes summary output:
      DATA LOGDM / .FALSE., .TRUE., 2 * .FALSE., .TRUE., 3 * .FALSE. /
C     Messages for the Lake Superior summary output:
      DATA LOGP  / 5 * .FALSE., .TRUE., 2 * .FALSE.  /
C     Messages for the detailed log and Lake Superior summary output:
      DATA LOGDP / .FALSE., .TRUE., 3 * .FALSE., .TRUE., 2 * .FALSE. /
C     Messages for the Lake Ontario summary output:
      DATA LOGO  / 6 * .FALSE., .TRUE., .FALSE.  /
C     Messages for the detailed log and Lake Ontario summary output:
      DATA LOGDO / .FALSE., .TRUE., 4 * .FALSE., .TRUE., .FALSE. /
C     Messages for the detailed log, middle lakes summary output,
C        Lake Ontario, and Lake Superior summary output:
      DATA LOGDMPO / .FALSE., .TRUE., 2* .FALSE., 3 *.TRUE., .FALSE./
C     Messages for the summary, middle lakes summary output,
C        Lake Ontario, and Lake Superior summary output:
      DATA LOGSMPO / 3 * .FALSE., 4 *.TRUE., .FALSE./
C     Messages for the summary, detailed, middle lakes summary output,
C        and Lake Superior summary output:
      DATA LOGSDMPO / .FALSE., .TRUE., .FALSE., 4 *.TRUE., .FALSE. /
C     Messages for the summary, detailed, middle lakes summary output:
      DATA LOGSDM / .FALSE.,.TRUE.,.FALSE., 2 * .TRUE., 3 * .FALSE. /
C     Messages for the summary, detailed, lake Superior summary output:
      DATA LOGSDP / .FALSE., .TRUE., .FALSE., .TRUE., .FALSE., 
     &    .TRUE., 2 * .FALSE./
C     Messages for the summary, detailed, Lake Ontario summary output:
      DATA LOGSDO / .FALSE.,.TRUE.,.FALSE.,.TRUE., 2* .FALSE., .TRUE.,
     &    .FALSE./
C
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
C     These DATA statements initialize variables to bogus values.
C     Viable input for each of these must be received.  CheckPar()
C     looks to see if these are still set to MissingVal (or whatever), 
C     and aborts if so.
C 
      DATA LStoSWP, SupCritRnd /2 * -9999/
C
      DATA BOMRound , FlowRound, LevelRound /3 * -9999/
C
      DATA CritBOMRound, CritFlowRound, CritLevelRound /3 * -9999/
      DATA Q1887Method, CritBMethod /2 * -9999/
      DATA CritCThresholdEl / 183.40D0 /
C      
      DATA (NUMAVSZ(I), I = 1,8)     / 8 * -9999 /
C     Initialize all the AVSZ elevations and areas to zero:
      DATA (((AVSZ(I,J,K), J=1,5), I = 1,8), K = 1,2) / 80* -9999.D0 /
C
      DATA ModelerName /'(not given)                             '/ 
C
      DATA StartYear, StartMonth, StartDay  /3 * -9999/
C
      DATA EndYear, EndMonth, EndDay  /3 * -9999/ 
C
      DATA MonthLength / -9999/
      DATA StMarysSolution, P77ACompat, MidLakeSolution, 
     &     UStLawrSolution  / 4 * -9999/
C
      DATA RoutingIncrement /3 * -9999/
      DATA RoutingInterval /3 * -9999/
C
      DATA SPSTLV,MHSTLV,SCSTLV,ERSTLV,CGSTLV,ONSTLV / 6 * -9999.D0 /
C
      DATA SPK,SPYM,SPA,SPB,SPWT,SPC / 6 * -9999.D0 /      
      DATA MHK,MHYM,MHA,MHB,MHWT,MHC / 6 * -9999.D0 /      
      DATA SCK,SCYM,SCA,SCB,SCWT,SCC / 6 * -9999.D0 /      
      DATA ERK,ERYM,ERA,ERB,ERWT,ERC / 6 * -9999.D0 /      
C
      DATA PrevStMarysFlow / -9999 /
      DATA PrevSupGates  / -9999. /
      DATA SPBOMLev, MHBOMLev, SCBOMLev, ERBOMLev  / 4 * -9999.D0 /     
C 
C
C
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
C     These DATA statements initialize variables to 'empty' or default
C     values.  Input for these may be optional or determined 
C     automatically as a consequence of other inputs.
C     
C 
C 
C     DATA MessDB(  1:30 ) / '                             '/
C     DATA MessDB( 31:60 ) / '                             '/
C     DATA MessDB( 61:90 ) / '                             '/
C     DATA MessDB( 91:120) / '                             '/
C     DATA MessDB(121:150) / '                             '/
C     DATA MessDB(151:184) / '                                 '/
C
      DATA OSTYPE /1/
C
      DATA (((WriteTSOutput(I,J,K), I=1,4),J=1,2),K=1,8) /64 * .FALSE. /
      DATA WriteErieOutflow / .FALSE. /
C       
      DATA GenVerbosity, TSDVerbosity, ChkVerbosity, SupVerbosity, 
     &     MLVerbosity, OntVerbosity /5, 5 * -9999 /
C
      DATA LakesToSolve / 8 * .FALSE. /
C
      DATA RCRUNTYP /'        '/
C
      DATA NOPDIR, NOPEXT /0,0/
C
      DATA RepeatIndex /0/
C
      DATA NPFILZ /0/
      DATA NUMPARAM /0/
C
      DATA SPRegCodes  / '                                        ' /  
      DATA SupBOMLevel, SWPierLevel / 2*'        ' /
      DATA SupPreProject, SupGateFlow, SupPlanFlow, 
     &     SupSideChanFlow, BalancingFlow / 4*'      ','        ' /
C
C          SupGatesOpen is REAL, not DOUBLE PRECISION
      DATA SupGatesOpen    / -9999.0 /
C
C
      DATA numMaxSC, numMinSCWinter, numMinSCSummer, numJacSC /4 * 0/
C
      DATA BackCompat / .FALSE. /
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
C     These DATA statements initialize superior regulation variables.
C     Initializing them in this manner probably isn't consistent with
C     the preference for no defaults (to minimize confusion).  But it 
C     is expedient, and in the case of Plan 1977A data, justifiable on
C     the grounds that for something as fixed and cohesive as Plan 1977A
C     the model should 'know' on its own what the plan is.  For 
C     flexibility in modeling changes to the plan, many of them have 
C     already been implemented as inputs, which override these default
C     values.  The rest could be also be implemented as inputs by 
C     following steps 3-5 of the instructions for adding inputs at the
C     top of this file.
C     
      DATA P2012SPHIGHLEVEL /183.90D0/
      DATA P2012SPLOWLEVEL  /182.888D0/
      DATA P2012MHTARGETADJ /0.016982D0/
      DATA COUNTER_JUNEFLOW /0/   
      DATA Last_JuneFlow /-9999/
      DATA P2012SPTARGET   / 183.36D0, 183.29D0, 183.24D0, 183.25D0,
     &                       183.30D0, 183.39D0, 183.45D0, 183.50D0,
     &                       183.52D0, 183.51D0, 183.48D0, 183.43D0 / 
      DATA P2012MHTARGET   / 176.23D0, 176.20D0, 176.20D0, 176.25D0,
     &                       176.36D0, 176.44D0, 176.50D0, 176.51D0,
     &                       176.46D0, 176.41D0, 176.33D0, 176.31D0 / 
      DATA P2012PARA  /0.75D0, 0.007617D0,  0.075D0, 4.90D0, 183.90D0,
     &              3800.0D0/
      Data ForcedSupGates /0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
     &                     0.5, 0.5, 0.5, 0.5/ 
      Data ForcedSideChan /2320, 2320, 2320, 2320, 2320, 2320, 
     &                     2320, 2320, 2320, 2320, 2320, 2320/  
      Data IndexGatesSide /0/  
C     The 'operational' NBS and retardations are peculiar to one mode
C     of P77ACOMPAT, and will hopefully go away.  Until then, this 
C     seems to be the most convenient place for them.
C
      DATA P77AAValue      /  18580.0D0 /
      DATA P77AMonthsToAve / 12 * 5 /
      DATA P77AMinimumFlow / 12 * 1560 /
      DATA P77AMaxChange   / 850 /
      DATA P77AWinterMax    /2410/
      DATA P77AUniformMonth /.TRUE./
      DATA P77ALowTrigger, P77AHighTrigger / 183.19D0, 183.71D0 /
      DATA P77ASup05NBS    /   790,  1730,  3510,  6710,  8160,  7110,  
     &                        5750,  4900,  5040,  3340,  2550,   680 /
      DATA P77ASup50NBS    /  -400,   250,  1190,  4160,  5240,  4390,
     &                        3620,  2780,  2010,  1080,   420,  -650 /
      DATA P77ASup95NBS    / -1470,  -850,  -590,  2070,  2580,  2240,
     &                        2010,   960,  -200, -1020, -1420, -2070 /
      DATA P77AMH05NBS     /  4130,  4980,  9260, 12230, 11330,  9060,
     &                        6340,  4530,  4930,  3680,  4020,  4050 /
      DATA P77AMH50NBS     /  1420,  2460,  5100,  7990,  7020,  5690,
     &                        3540,  1500,   740,  -110,   930,   710 /
      DATA P77AMH95NBS     /  -650,   400,  1640,  4530,  3680,  3200,
     &                        1440, -1160, -2290, -2860, -2010, -2150 /
      DATA P77AStC50NBS    /   140,   170,   250,   230,   140,    80,   
     &                          60,    30,    30,    30,    60,   110 /
      DATA P77AErie50NBS   /   620,   960,  2010,  1900,  1250,   790,
     &                          80,  -370,  -570,  -680,  -200,   450 /
      DATA P77ALLOgokiDiv   / 12 *  159 /
      DATA P77AChicagoDiv         / 12 *  -91 /
      DATA P77AWellandFlow  / 12 * -261 /
      DATA P77AStMarysRetrd / 12 * 0 /
      DATA P77AStClairRetrd / 590, 480, 110, 110, 7 * 0, 110 /
      DATA P77ADetroitRetrd / 650, 510, 110, 8 * 0, 140 /
      DATA P77ANiagaraRetrd / 110, 140,  80, 140, 0, 60, 
     &                        140, 110,  80,  60, 0,  0 /
      DATA P77ABaseFlow     / 2010, 2010, 1980, 1980, 2180, 2350, 
     &                        2490, 2460, 2490, 2440, 2440, 2040  /
      DATA P77ASupTargetLev / 183.34D0, 183.27D0, 183.22D0, 183.20D0, 
     &                        183.28D0, 183.38D0, 183.46D0, 183.50D0, 
     &                        183.52D0, 183.51D0, 183.47D0, 183.42D0 /
      DATA P77ASupRangeFact /   0.13D0,   0.13D0,   0.13D0,   0.13D0, 
     &                          0.15D0,   0.14D0,   0.14D0,   0.13D0,
     &                          0.13D0,   0.12D0,   0.13D0,   0.12D0 /
      DATA P77AMHTargetLev /  176.38D0, 176.35D0, 176.35D0, 176.40D0, 
     &                        176.51D0, 176.59D0, 176.65D0, 176.66D0, 
     &                        176.62D0, 176.57D0, 176.49D0, 176.45D0 /
      DATA P77AMHRangeFact /    0.35D0,   0.35D0,   0.34D0,   0.36D0, 
     &                          0.36D0,   0.37D0,   0.37D0,   0.38D0,
     &                          0.38D0,   0.37D0,   0.37D0,   0.37D0 /
      DATA P77ASPK, P77ASPYM, P77ASPA, P77ASPB, P77ASPWT, P77ASPC /    
     &     824.7D0, 181.43D0,   1.5D0,    0.D0,     1.D0,    0.D0 /
      DATA P77AMHK, P77AMHYM, P77AMHA, P77AMHB, P77AMHWT, P77AMHC /
     &      46.4D0, 165.83D0,  2.00D0,  0.50D0,   0.50D0, -0.02D0 /
      DATA P77ASCK, P77ASCYM, P77ASCA, P77ASCB, P77ASCWT, P77ASCC /
     &      70.7D0, 165.94D0,  2.00D0,  0.50D0,   1.00D0, -0.01D0 /
      DATA P77AERK, P77AERYM, P77AERA, P77AERB, P77AERWT, P77AERC /
     &     616.7D0, 169.75D0,  1.50D0,  0.00D0,   1.00D0,    0.D0 /  
      DATA P77ACGIP / 171.16D0 /
      DATA P77AMidLakesInc, P77AMidLakesTol / 40, -2 /
      DATA P77AMonthRound, P77AFlowRound, P77ALevelRound / -2, 0, -4 /
C
C
C     SAVE X, Y, Z, MaxSetting
C      
C imperial IGLD 1955
C      DATA X /83.4, 196.8, 393.0, 588.6, 833, 1077, 1350, 1618,
C     1        1874, 2022, 2112, 2213, 2278, 2364, 2458, 2538, 2601/
C
C      DATA Y /1, 1, 1.0, 1, 1, 1, 1, 0.951, 0.954, 0.951, 0.946,
C     1        0.951, 0.949, 0.944, 0.946, 0.944, 0.944/
C
C      DATA Z /591.81, 591.81, 591.38, 591.62, 592.16, 592.65,
C     1        593.08, 564.048, 566.172, 564.478, 561.478, 564.498,
C     2        563.272, 560.262, 561.518, 560.322, 560.292/
C
C metric IGLD 1985 USING LAKE SUPERIOR IGLD55 TO IGLD85 CONVERSION
C                       OF .366 -.1 HYDRAULIC CORRECTOR
C
      DATA SGX /14.03D0, 33.12D0, 66.13D0, 99.05D0, 140.17D0, 181.23D0, 
     1        227.17D0, 272.27D0, 315.35D0, 340.25D0, 355.4D0, 372.39D0,
     2        383.33D0, 397.8D0, 413.62D0, 427.08D0, 437.69D0/
C
      DATA SGY /1.D0, 1.D0, 1.0D0, 1.D0, 1.D0, 1.D0, 1.D0, 0.951D0, 
     1        0.954D0, 0.951D0, 0.946D0, 0.951D0, 0.949D0, 0.944D0, 
     2        0.946D0, 0.944D0, 0.944D0/
C
      DATA SGZ /180.750D0, 180.750D0, 180.619D0, 180.692D0, 180.856D0,
     1        181.006D0, 181.137D0, 172.270D0, 172.918D0, 172.401D0,
     2        171.485D0, 172.407D0, 172.033D0, 171.113D0, 171.497D0,
     3        171.132D0, 171.122D0/
C
      DATA MaxSGSetting /17/
C
      DATA SWPGateEqA / 33.12D0, 66.13D0, 99.05D0, 140.17D0, 181.23D0, 
     &   227.17D0, 272.27D0, 315.35D0, 340.25D0, 355.40D0, 372.39D0, 
     &   383.33D0, 397.80D0, 413.62D0, 427.08D0, 437.69D0 /
C
      DATA SWPGateEqB / 180.59D0, 180.47D0, 180.54D0, 180.70D0, 
     &   180.84D0, 180.97D0, 181.09D0, 181.18D0, 181.21D0, 181.21D0, 
     &   181.22D0, 181.21D0, 181.21D0, 181.22D0, 181.22D0, 181.21D0 /
C
      DATA OldSWPA / 27.6027D0, 18.9144D0, 17.6510D0, 11.4080D0, 
     &   17.2385D0, 14.9971D0, 8.7544D0, 9.9005D0, 11.8851D0, 20.3023D0, 
     &   16.7230D0, 13.6718D0 /
      DATA OldSWPB / 1.1534D0, 1.1060D0, 1.0991D0, 1.0648D0, 1.0965D0, 
     &   1.0842D0, 1.0501D0, 1.0564D0, 1.0672D0, 1.1132D0, 1.0936D0, 
     &   1.0771D0 /
      DATA OldSWPC / 0.27D0 /
C
      DATA EC2000SWPA, EC2000SWPB / 1.046385D0, 0.10015D0 /
      DATA EC2000SWPC / 8.152D0, 8.161D0, 8.163D0, 8.158D0, 8.129D0, 
     &   8.129D0, 8.125D0, 8.124D0, 8.125D0, 8.126D0, 8.118D0, 8.127D0 /
C
C     The values below were used in the original implementation of the EC2000
C     relationship, but have been corrected to the values above.
C      DATA EC2000SWPC / 7.921D0, 7.930D0, 7.932D0, 7.927D0, 7.898D0, 
C     &   7.898D0, 7.894D0, 7.893D0, 7.894D0, 7.895D0, 7.887D0, 7.896D0 /
C
      DATA EC2008SWPA, EC2008SWPB / 1.033806D0, 0.08622D0 /
      DATA EC2008SWPC / 5.916D0, 5.923D0, 5.925D0, 5.923D0, 5.897D0, 
     &   5.897D0, 5.892D0, 5.891D0, 5.894D0, 5.890D0, 5.886D0, 5.884D0 /
C
      DATA USSlipEquation / 1 /     
C
C     CWOvertopElev is the top of the Compensating Works gates, when closed.
C     This also assumes that the Canadian berm (upstream of Whitefish Island)
C     stands firm.
      DATA CWOvertopElev /184.0D0/
C     Assumed from actual observed intermittent overtopping at a daily level 
C     of 183.91 (from Rob Caldwell/David Fay (2008))
      DATA CWMinimumFreeboard /12 * 0.10D0/

C     Use default sidechannel physical capacities such that it does not
C     change answers that the Larson and Dumont models would create (i.e., 
C     somtimes up to 2518 at levels around 183.7).  Override these defaults 
C     in input if wish to set a limit on Side Channel capacity. 
      DATA SideChanMaxCapacity /12 * 9999/
      DATA SideChanMinCapacity /12 * 0/
      
C     Data are updated on Nov 17 2014 for plan 2012 (only the first column), used in Criterion B check program      
C     The data are the average of coordinated NBS (1900-2008, actual month length, used in routing).  Only the first month 
C     is actually used in Plan 2012's Crit "b" (but six months specified as this was required in Plan 77A's operational forecast)
C     -330	130	1320	4280	5060	4400	3600	2650	1970	1170	480	-570
     
      DATA (SPOPNBS08( 1,I),I=1,6) / -330,  130, 1320, 4280, 5060, 4400/
      DATA (SPOPNBS08( 2,I),I=1,6) /  130, 1320, 4280, 5060, 4400, 3600/
      DATA (SPOPNBS08( 3,I),I=1,6) / 1320, 4280, 5060, 4400, 3600, 2650/
      DATA (SPOPNBS08( 4,I),I=1,6) / 4280, 5060, 4400, 3600, 2650, 1970/
      DATA (SPOPNBS08( 5,I),I=1,6) / 5060, 4400, 3600, 2650, 1970, 1170/
      DATA (SPOPNBS08( 6,I),I=1,6) / 4400, 3600, 2650, 1970, 1170,  480/
      DATA (SPOPNBS08( 7,I),I=1,6) / 3600, 2650, 1970, 1170,  480, -570/
      DATA (SPOPNBS08( 8,I),I=1,6) / 2650, 1970, 1170,  480, -570, -330/
      DATA (SPOPNBS08( 9,I),I=1,6) / 1970, 1170,  480, -570, -330,  130/
      DATA (SPOPNBS08(10,I),I=1,6) / 1170,  480, -570, -330,  130, 1320/
      DATA (SPOPNBS08(11,I),I=1,6) /  480, -570, -330,  130, 1320, 4280/
      DATA (SPOPNBS08(12,I),I=1,6) / -570, -330,  130, 1320, 4280, 5060/

C     1720	2470	5150	8100	7050	5800	3620	1530	730	90	1100	1000

      DATA (MHOPNBS08( 1,I),I=1,6) / 1720, 2470, 5150, 8100, 7050, 5800/
      DATA (MHOPNBS08( 2,I),I=1,6) / 2470, 5150, 8100, 7050, 5800, 3620/
      DATA (MHOPNBS08( 3,I),I=1,6) / 5150, 8100, 7050, 5800, 3620, 1530/
      DATA (MHOPNBS08( 4,I),I=1,6) / 8100, 7050, 5800, 3620, 1530,  730/
      DATA (MHOPNBS08( 5,I),I=1,6) / 7050, 5800, 3620, 1530,  730,   90/
      DATA (MHOPNBS08( 6,I),I=1,6) / 5800, 3620, 1530,  730,   90, 1100/
      DATA (MHOPNBS08( 7,I),I=1,6) / 3620, 1530,  730,   90, 1100, 1000/
      DATA (MHOPNBS08( 8,I),I=1,6) / 1530,  730,   90, 1100, 1000, 1720/
      DATA (MHOPNBS08( 9,I),I=1,6) /  730,   90, 1100, 1000, 1720, 2470/
      DATA (MHOPNBS08(10,I),I=1,6) /   90, 1100, 1000, 1720, 2470, 5150/
      DATA (MHOPNBS08(11,I),I=1,6) / 1100, 1000, 1720, 2470, 5150, 8100/
      DATA (MHOPNBS08(12,I),I=1,6) / 1000, 1720, 2470, 5150, 8100, 7050/

C    190	220	270	230	150	100	80	40	40	50	70	150
     
      DATA (SCOPNBS08( 1,I),I=1,6) /  190,  220,  270,  230,  150,  100/
      DATA (SCOPNBS08( 2,I),I=1,6) /  220,  270,  230,  150,  100,   80/
      DATA (SCOPNBS08( 3,I),I=1,6) /  270,  230,  150,  100,   80,   40/
      DATA (SCOPNBS08( 4,I),I=1,6) /  230,  150,  100,   80,   40,   40/
      DATA (SCOPNBS08( 5,I),I=1,6) /  150,  100,   80,   40,   40,   50/
      DATA (SCOPNBS08( 6,I),I=1,6) /  100,   80,   40,   40,   50,   70/
      DATA (SCOPNBS08( 7,I),I=1,6) /   80,   40,   40,   50,   70,  150/
      DATA (SCOPNBS08( 8,I),I=1,6) /   40,   40,   50,   70,  150,  190/
      DATA (SCOPNBS08( 9,I),I=1,6) /   40,   50,   70,  140,  190,  220/
      DATA (SCOPNBS08(10,I),I=1,6) /   50,   70,  140,  190,  220,  270/
      DATA (SCOPNBS08(11,I),I=1,6) /   70,  140,  190,  220,  270,  230/
      DATA (SCOPNBS08(12,I),I=1,6) /  153,  190,  220,  270,  230,  150/

C      820	1090  2010  1930  1340  900  190  -260  -500  -570  -80  610

      DATA (EROPNBS08( 1,I),I=1,6) /  820, 1090, 2010, 1930, 1340,  900/
      DATA (EROPNBS08( 2,I),I=1,6) / 1090, 2010, 1930, 1340,  900,  190/
      DATA (EROPNBS08( 3,I),I=1,6) / 2010, 1930, 1340,  900,  190, -260/
      DATA (EROPNBS08( 4,I),I=1,6) / 1930, 1340,  900,  190, -260, -500/
      DATA (EROPNBS08( 5,I),I=1,6) / 1340,  900,  190, -260, -500, -570/
      DATA (EROPNBS08( 6,I),I=1,6) /  900,  190, -260, -500, -570,  -80/
      DATA (EROPNBS08( 7,I),I=1,6) /  190, -260, -500, -570,  -80,  610/
      DATA (EROPNBS08( 8,I),I=1,6) / -260, -500, -570,  -80,  610,  820/
      DATA (EROPNBS08( 9,I),I=1,6) / -500, -570,  -80,  610,  820, 1090/
      DATA (EROPNBS08(10,I),I=1,6) / -570,  -80,  610,  820, 1090, 2010/
      DATA (EROPNBS08(11,I),I=1,6) /  -80,  610,  820, 1090, 2010, 1930/
      DATA (EROPNBS08(12,I),I=1,6) /  610,  820, 1090, 2010, 1930, 1340/
C
C     The operational NBS for Lake Superior, MH, SC and Erie
C
C     These values provided by Rob and added by Yin on Oct 1 2007. 
C     The value of P77ACompat determines whether to use these updated  
C     values or the older ones.
C
      DATA (SPOPNBS07( 1,I),I=1,6) / -360,  280, 1300, 4200, 5200, 4410/
      DATA (SPOPNBS07( 2,I),I=1,6) /  260, 1280, 4180, 5180, 4390, 3680/
      DATA (SPOPNBS07( 3,I),I=1,6) / 1300, 4200, 5200, 4410, 3700, 2780/
      DATA (SPOPNBS07( 4,I),I=1,6) / 4150, 5150, 4360, 3650, 2730, 1980/
      DATA (SPOPNBS07( 5,I),I=1,6) / 5240, 4450, 3740, 2820, 2070, 1140/
      DATA (SPOPNBS07( 6,I),I=1,6) / 4410, 3700, 2780, 2030, 1100,  510/
      DATA (SPOPNBS07( 7,I),I=1,6) / 3680, 2760, 2010, 1080,  490, -660/
      DATA (SPOPNBS07( 8,I),I=1,6) / 2750, 2000, 1070,  480, -670, -390/
      DATA (SPOPNBS07( 9,I),I=1,6) / 2010, 1080,  490, -660, -380,  260/
      DATA (SPOPNBS07(10,I),I=1,6) / 1080,  490, -660, -380,  260, 1280/
      DATA (SPOPNBS07(11,I),I=1,6) /  530, -620, -340,  300, 1320, 4220/
      DATA (SPOPNBS07(12,I),I=1,6) / -630, -350,  290, 1310, 4210, 5210/

      DATA (MHOPNBS07( 1,I),I=1,6) / 1630, 2470, 5220, 8040, 7070, 5760/
      DATA (MHOPNBS07( 2,I),I=1,6) / 2530, 5280, 8100, 7130, 5820, 3710/
      DATA (MHOPNBS07( 3,I),I=1,6) / 5220, 8040, 7070, 5760, 3650, 1560/
      DATA (MHOPNBS07( 4,I),I=1,6) / 8140, 7170, 5860, 3750, 1660,  910/
      DATA (MHOPNBS07( 5,I),I=1,6) / 7110, 5800, 3690, 1600,  850,   90/
      DATA (MHOPNBS07( 6,I),I=1,6) / 5760, 3650, 1560,  810,   50, 1110/
      DATA (MHOPNBS07( 7,I),I=1,6) / 3630, 1540,  790,   30, 1090,  830/
      DATA (MHOPNBS07( 8,I),I=1,6) / 1560,  810,   50, 1110,  850, 1630/
      DATA (MHOPNBS07( 9,I),I=1,6) /  810,   50, 1110,  850, 1630, 2470/
      DATA (MHOPNBS07(10,I),I=1,6) /   30, 1090,  830, 1610, 2450, 5200/
      DATA (MHOPNBS07(11,I),I=1,6) / 1110,  850, 1630, 2470, 5220, 8040/
      DATA (MHOPNBS07(12,I),I=1,6) /  860, 1640, 2480, 5230, 8050, 7080/

      DATA (SCOPNBS07( 1,I),I=1,6) /  190,  210,  270,  230,  150,  110/
      DATA (SCOPNBS07( 2,I),I=1,6) /  220,  280,  240,  160,  120,   90/
      DATA (SCOPNBS07( 3,I),I=1,6) /  270,  230,  150,  110,   80,   50/
      DATA (SCOPNBS07( 4,I),I=1,6) /  220,  140,  100,   70,   40,   30/
      DATA (SCOPNBS07( 5,I),I=1,6) /  150,  110,   80,   50,   40,   50/
      DATA (SCOPNBS07( 6,I),I=1,6) /  100,   70,   40,   30,   40,   60/
      DATA (SCOPNBS07( 7,I),I=1,6) /   80,   50,   40,   50,   70,  150/
      DATA (SCOPNBS07( 8,I),I=1,6) /   50,   40,   50,   70,  150,  190/
      DATA (SCOPNBS07( 9,I),I=1,6) /   30,   40,   60,  140,  180,  200/
      DATA (SCOPNBS07(10,I),I=1,6) /   40,   60,  140,  180,  200,  260/
      DATA (SCOPNBS07(11,I),I=1,6) /   60,  140,  180,  200,  260,  220/
      DATA (SCOPNBS07(12,I),I=1,6) /  150,  190,  210,  270,  230,  150/
C
      DATA (EROPNBS07( 1,I),I=1,6) /  770, 1040, 2030, 1880, 1300,  860/
      DATA (EROPNBS07( 2,I),I=1,6) / 1050, 2040, 1890, 1310,  870,  150/
      DATA (EROPNBS07( 3,I),I=1,6) / 2030, 1880, 1300,  860,  140, -330/
      DATA (EROPNBS07( 4,I),I=1,6) / 1860, 1280,  840,  120, -350, -530/
      DATA (EROPNBS07( 5,I),I=1,6) / 1270,  830,  110, -360, -540, -670/
      DATA (EROPNBS07( 6,I),I=1,6) /  840,  120, -350, -530, -660, -120/
      DATA (EROPNBS07( 7,I),I=1,6) /  110, -360, -540, -670, -130,  490/
      DATA (EROPNBS07( 8,I),I=1,6) / -350, -530, -660, -120,  500,  750/
      DATA (EROPNBS07( 9,I),I=1,6) / -520, -650, -110,  510,  760, 1030/
      DATA (EROPNBS07(10,I),I=1,6) / -660, -120,  500,  750, 1020, 2010/
      DATA (EROPNBS07(11,I),I=1,6) / -110,  510,  760, 1030, 2020, 1870/
      DATA (EROPNBS07(12,I),I=1,6) /  520,  770, 1040, 2030, 1880, 1300/
C
C
C     The following data corresponds to P77ACompat=3, in use from 199?-2007
      DATA (SPOPNBS06( 1,I),I=1,6) / -380,  300, 1280, 4210, 5260, 4470/
      DATA (SPOPNBS06( 2,I),I=1,6) /  280, 1260, 4190, 5240, 4450, 3650/
      DATA (SPOPNBS06( 3,I),I=1,6) / 1280, 4210, 5260, 4470, 3670, 2850/
      DATA (SPOPNBS06( 4,I),I=1,6) / 4210, 5260, 4470, 3670, 2850, 2080/
      DATA (SPOPNBS06( 5,I),I=1,6) / 5260, 4470, 3670, 2850, 2080, 1070/
      DATA (SPOPNBS06( 6,I),I=1,6) / 4460, 3660, 2840, 2070, 1060,  510/
      DATA (SPOPNBS06( 7,I),I=1,6) / 3640, 2820, 2050, 1040,  490, -700/
      DATA (SPOPNBS06( 8,I),I=1,6) / 2830, 2060, 1050,  500, -690, -400/
      DATA (SPOPNBS06( 9,I),I=1,6) / 2060, 1050,  500, -690, -400,  280/
      DATA (SPOPNBS06(10,I),I=1,6) / 1050,  500, -690 ,-400,  280, 1260/
      DATA (SPOPNBS06(11,I),I=1,6) /  520, -670, -380,  300, 1280, 4210/
      DATA (SPOPNBS06(12,I),I=1,6) / -670, -380,  300, 1280, 4210, 5260/
C
      DATA (MHOPNBS06( 1,I),I=1,6) / 1490, 2480, 5220, 8090, 7100, 5780/
      DATA (MHOPNBS06( 2,I),I=1,6) / 2560, 5300, 8170, 7180, 5860, 3680/
      DATA (MHOPNBS06( 3,I),I=1,6) / 5220, 8090, 7100, 5780, 3600, 1560/
      DATA (MHOPNBS06( 4,I),I=1,6) / 8090, 7100, 5780, 3600, 1560,  880/
      DATA (MHOPNBS06( 5,I),I=1,6) / 7150, 5830, 3650, 1610,  930,   80/
      DATA (MHOPNBS06( 6,I),I=1,6) / 5760, 3580, 1540,  860,   10, 1010/
      DATA (MHOPNBS06( 7,I),I=1,6) / 3560, 1520,  840,  -10,  990,  770/
      DATA (MHOPNBS06( 8,I),I=1,6) / 1550,  870,   20, 1020,  800, 1480/
      DATA (MHOPNBS06( 9,I),I=1,6) /  890,   40, 1040,  820, 1500, 2490/
      DATA (MHOPNBS06(10,I),I=1,6) /   10, 1010,  790, 1470, 2460, 5200/
      DATA (MHOPNBS06(11,I),I=1,6) / 1050,  830, 1510, 2500, 5240, 8110/
      DATA (MHOPNBS06(12,I),I=1,6) /  830, 1510, 2500, 5240, 8110, 7120/

C
      DATA (SCOPNBS06( 1,I),I=1,6) /  160,  190,  270,  220,  150,  100/
      DATA (SCOPNBS06( 2,I),I=1,6) /  190,  270,  220,  150,  100,   70/
      DATA (SCOPNBS06( 3,I),I=1,6) /  270,  220,  150,  100,   70,   40/
      DATA (SCOPNBS06( 4,I),I=1,6) /  220,  150,  100,   70,   40,   30/
      DATA (SCOPNBS06( 5,I),I=1,6) /  150,  100,   70,   40,   30,   40/
      DATA (SCOPNBS06( 6,I),I=1,6) /   90,   60,   30,   20,   30,   50/
      DATA (SCOPNBS06( 7,I),I=1,6) /   70,   40,   30,   40,   60,  140/
      DATA (SCOPNBS06( 8,I),I=1,6) /   30,   20,   30,   50,  130,  150/
      DATA (SCOPNBS06( 9,I),I=1,6) /   30,   40,   60,  140,  160,  190/
      DATA (SCOPNBS06(10,I),I=1,6) /   30,   50,  130,  150,  180,  260/
      DATA (SCOPNBS06(11,I),I=1,6) /   60,  140,  160,  190,  270,  220/
      DATA (SCOPNBS06(12,I),I=1,6) /  140,  160,  190,  270,  220,  150/
C
      DATA (EROPNBS06( 1,I),I=1,6) /  700, 1000, 2040, 1870, 1300,  840/
      DATA (EROPNBS06( 2,I),I=1,6) / 1010, 2050, 1880, 1310,  850,  130/
      DATA (EROPNBS06( 3,I),I=1,6) / 2040, 1870, 1300,  840,  120, -350/
      DATA (EROPNBS06( 4,I),I=1,6) / 1860, 1290,  830,  110, -360, -530/
      DATA (EROPNBS06( 5,I),I=1,6) / 1290,  830,  110, -360, -530, -660/
      DATA (EROPNBS06( 6,I),I=1,6) /  820,  100, -370, -540, -670, -160/
      DATA (EROPNBS06( 7,I),I=1,6) /  100, -370, -540, -670, -160,  460/
      DATA (EROPNBS06( 8,I),I=1,6) / -360, -530, -660, -150,  470,  690/
      DATA (EROPNBS06( 9,I),I=1,6) / -530, -660, -150,  470,  690,  990/
      DATA (EROPNBS06(10,I),I=1,6) / -670, -160,  460,  680,  980, 2020/
      DATA (EROPNBS06(11,I),I=1,6) / -140,  480,  700, 1000, 2040, 1870/
      DATA (EROPNBS06(12,I),I=1,6) /  490,  710, 1010, 2050, 1880, 1310/
C
C
C
C
      DATA (OPLLO(I),I= 1,12) / 113, 99, 85, 79, 181, 252, 207, 170, 
     &   133, 127, 122, 125 / 
C
      DATA (OPWelland(I),I= 1,12)  / -218, -215, -218, -252, -255, -255, 
     &   -249, -252, -258, -261, -258, -227 / 
C
      DATA (SCRRet(I),I= 1,12)  / 590,480,110,110,0,0,0,0,0,0,0,110 / 
C
      DATA (DETRet(I),I= 1,12)  / 650,510,110, 0, 0,0,0,0,0,0,0,140 / 
C
      DATA (NIARet(I),I= 1,12) / 110,140,80,140,0,60,230,140,80,60,0,0/


C     Added on Nov 19 2014 for Plan2012 (used in P2012_Criterion)
C     220,211,222,247,244,246,248,239,244,243,243,237
      DATA (OPWelland_P2012(I),I= 1,12)  / -220, -211, -222, -247, 
     &             -244, -246, -248, -239, -244, -243, -243, -237 / 
C     565	498	172	141	20	25	0	-19	-18	5	-8	85
      DATA (SCRRet_P2012(I),I= 1,12)  / 565, 498, 172, 141, 20, 25,
     &                                  0, -19,  -18, 5, -8 ,85 /   
C     503	371	66	-24	-9	139	255	278	246	104	-63	39
      DATA (DETRet_P2012(I),I= 1,12)  / 503, 371, 66, -24, -9, 139, 
     &                                  255, 278, 246, 104, -63, 39 /      
C     -127	-39	-18	38	-49	106	155	125	56	-79	-200	-239     
      DATA (NIARet_P2012(I),I= 1,12) / -127, -39, -18, 38, -49, 106, 
     &                                  155, 125, 56, 79, -200, -239/
              
      DATA OPChicago / -91 /
      DATA MonthsChecked / 0 /
      DATA ForecastMonths / 0 /
C
C
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
C     These DATA statements initialize ontario regulation variables.
C     Initializing them in this manner probably isn't consistent with
C     the preference for no defaults (to minimize confusion).  But that
C     is left for others to remedy   ;-)
C
      DATA CritDepths          / 1 /
      DATA P58DUseDatabase     / .FALSE. /
      DATA ElimDev             / .FALSE. /
      DATA SetFlow             / .FALSE. /
      DATA PreProjConditions   / .FALSE. /
      DATA P58DTotalSup        / 6250 /
      DATA MaxILimit           / 0/
      DATA MinPLimit           / 0 /
      DATA P58DWeightedSupply  / 6230 /
      DATA P58DRuleCurve       / 5750 /
      DATA P58DAppLimit        / 5950 /
      DATA P58DPlanLevel       / 74.55D0 /
      DATA P58DActualLevel     / 74.55D0 /
      DATA P58DRunType         / 2 / 
      DATA P58DTotaldev        / 0 /
      DATA P58DAdjustment      / 0 /
      DATA P58DrcAdjustment    / 0 /
      DATA P58DLSLflow         / 0 /
      DATA P58DActualflow      / 600 /
      DATA P58DForecastDay, P58DForecastMon, P58DForecastYear /3 * -1 /
      DATA P58DAppLimitSym     / '  '/
      DATA ForceFlow           / .FALSE. /
C     DATA P58DSimulateFile(1:34) / '   (none)                        '/
C     DATA P58DSimulateFile(35:65) /'                               '/
C     DATA P58DSimulateFile(66:80) /'               '/
C     DATA P58DSIDataFile(1:34) / '   (none)                        '/
C     DATA P58DSIDataFile(35:65) /'                               '/
C     DATA P58DSIDataFile(66:80) /'               '/
C     DATA SupplyIndicator     /  0, -68, -71, -68, -62, -59, -59, -59,
C    &                          -59, -62, -65, -68, -70, -59, -59, -59,
C    &                          -59, -59, -59, -59, -59/
C    Paul expects the first value to be 0
      DATA SupplyIndicator     /  0, 20 * -9999 /
      DATA PreprojectLevel     / 69.33D0 /
      DATA P58DSimulateFile(  1:30 ) / '                             '/
      DATA P58DSimulateFile( 31:60 ) / '                             '/
      DATA P58DSimulateFile( 61:90 ) / '                             '/
      DATA P58DSimulateFile( 91:120) / '                             '/
      DATA P58DSimulateFile(121:150) / '                             '/
      DATA P58DSimulateFile(151:180) / '                             '/
      DATA P58DSimulateFile(181:210) / '                             '/
      DATA P58DSimulateFile(211:240) / '                             '/
      DATA P58DSimulateFile(241:256) / '                '/
      DATA P58DSIDataFile(  1:30 ) / '                             '/
      DATA P58DSIDataFile( 31:60 ) / '                             '/
      DATA P58DSIDataFile( 61:90 ) / '                             '/
      DATA P58DSIDataFile( 91:120) / '                             '/
      DATA P58DSIDataFile(121:150) / '                             '/
      DATA P58DSIDataFile(151:180) / '                             '/
      DATA P58DSIDataFile(181:210) / '                             '/
      DATA P58DSIDataFile(211:240) / '                             '/
      DATA P58DSIDataFile(241:256) / '                '/
      DATA CritKfile(  1:30 ) / '                             '/
      DATA CritKfile( 31:60 ) / '                             '/
      DATA CritKfile( 61:90 ) / '                             '/
      DATA CritKfile( 91:120) / '                             '/
      DATA CritKfile(121:150) / '                             '/
      DATA CritKfile(151:180) / '                             '/
      DATA CritKfile(181:210) / '                             '/
      DATA CritKfile(211:240) / '                             '/
      DATA CritKfile(241:256) / '                '/
C
C
      END
C      
C
C-----------------------------------------------------------------------
C      
      BLOCK DATA PLINIT 
C
      IMPLICIT NONE
C
      COMMON /CPARID/   PARIDS
C
      CHARACTER*6400 PARIDS
C               This is a concatenated list of all the parameters that
C               the program knows to look for.  When a parameter ID is read,
C               the program searches for it in this list.  Since all the
C               parameters are sixteen chars, dividing the return value
C               from INDEX() by sixteen yields the order number of the 
C               parameter ID in the list.  This order number is used to
C               jump to the label appropriate for processing input lines
C               for that parameter.  
C
      DATA PARIDS(0001:0016) /'MODELERNAME     '/                          1
      DATA PARIDS(0017:0032) /'MESSAGEDATABASE '/                          2
      DATA PARIDS(0033:0048) /'OUTPUTDIRECTORY '/                          3
      DATA PARIDS(0049:0064) /'OUTPUTEXTENSION '/                          4
      DATA PARIDS(0065:0080) /'TITLE1          '/                          5
      DATA PARIDS(0081:0096) /'TITLE2          '/                          6
      DATA PARIDS(0097:0112) /'TITLE3          '/                          7
      DATA PARIDS(0113:0128) /'????????????????'/                          8
      DATA PARIDS(0129:0144) /'ONLYONEOPEN     '/                          9
      DATA PARIDS(0145:0160) /'GENERALVERBOSITY'/                         10
      DATA PARIDS(0161:0176) /'TSDATAVERBOSITY '/                         11
      DATA PARIDS(0177:0192) /'CHECKVERBOSITY  '/                         12
      DATA PARIDS(0193:0208) /'SUPVERBOSITY    '/                         13
      DATA PARIDS(0209:0224) /'MIDLAKEVERBOSITY'/                         14
      DATA PARIDS(0225:0240) /'ONTVERBOSITY    '/                         15
      DATA PARIDS(0241:0256) /'STARTDATE       '/                         16
      DATA PARIDS(0257:0272) /'ENDDATE         '/                         17
      DATA PARIDS(0273:0288) /'MONTHLENGTH     '/                         18
      DATA PARIDS(0289:0304) /'LAKESTOSOLVEFOR '/                         19
      DATA PARIDS(0305:0320) /'SUPSTARTLEVEL   '/                         20
      DATA PARIDS(0321:0336) /'SUPSIDECHANDATA '/                         21
      DATA PARIDS(0337:0352) /'SUPNBSDATA      '/                         22
      DATA PARIDS(0353:0368) /'SUPCONUSEDATA   '/                         23
      DATA PARIDS(0369:0384) /'LONGLACDIVDATA  '/                         24
      DATA PARIDS(0385:0400) /'OGOKIDIVDATA    '/                         25
      DATA PARIDS(0401:0416) /'COMBLLODIVDATA  '/                         26
      DATA PARIDS(0417:0432) /'SUPOTHERDIVDATA '/                         27
      DATA PARIDS(0433:0448) /'SUPOTHERSUPDATA '/                         28
      DATA PARIDS(0449:0464) /'SUPRUNOFFDATA   '/                         29
      DATA PARIDS(0465:0480) /'SUPPRECIPDATA   '/                         30
      DATA PARIDS(0481:0496) /'SUPEVAPDATA     '/                         31
      DATA PARIDS(0497:0512) /'SUPDEVIATIONDATA'/                         32
      DATA PARIDS(0513:0528) /'SUPEOPLEVELDATA '/                         33
      DATA PARIDS(0529:0544) /'SUPGATESDATA    '/                         34
      DATA PARIDS(0545:0560) /'SUPOUTFLOWDATA  '/                         35
      DATA PARIDS(0561:0576) /'MHUSTARTLEVEL   '/                         36
      DATA PARIDS(0577:0592) /'????????????????'/                         37
      DATA PARIDS(0593:0608) /'MHUNBSDATA      '/                         38
      DATA PARIDS(0609:0624) /'MHUCONUSEDATA   '/                         39
      DATA PARIDS(0625:0640) /'CHICAGODIVDATA  '/                         40
      DATA PARIDS(0641:0656) /'MHUOTHERDIVDATA '/                         41
      DATA PARIDS(0657:0672) /'MHUOTHERSUPDATA '/                         42
      DATA PARIDS(0673:0688) /'MHURUNOFFDATA   '/                         43
      DATA PARIDS(0689:0704) /'MHUPRECIPDATA   '/                         44
      DATA PARIDS(0705:0720) /'MHUEVAPDATA     '/                         45
      DATA PARIDS(0721:0736) /'MHUEOPLEVELDATA '/                         46
      DATA PARIDS(0737:0752) /'????????????????'/                         47
      DATA PARIDS(0753:0768) /'????????????????'/                         48
      DATA PARIDS(0769:0784) /'MICCONUSEDATA   '/                         49
      DATA PARIDS(0785:0800) /'MICOTHERDIVDATA '/                         50
      DATA PARIDS(0801:0816) /'MICOTHERSUPDATA '/                         51
      DATA PARIDS(0817:0832) /'MICRUNOFFDATA   '/                         52
      DATA PARIDS(0833:0848) /'MICPRECIPDATA   '/                         53
      DATA PARIDS(0849:0864) /'MICEVAPDATA     '/                         54
      DATA PARIDS(0865:0880) /'????????????????'/                         55
      DATA PARIDS(0881:0896) /'HURCONUSEDATA   '/                         56
      DATA PARIDS(0897:0912) /'HUROTHERDIVDATA '/                         57
      DATA PARIDS(0913:0928) /'HUROTHERSUPDATA '/                         58
      DATA PARIDS(0929:0944) /'HURRUNOFFDATA   '/                         59
      DATA PARIDS(0945:0960) /'HURPRECIPDATA   '/                         60
      DATA PARIDS(0961:0976) /'HUREVAPDATA     '/                         61
      DATA PARIDS(0977:0992) /'STCSTARTLEVEL   '/                         62
      DATA PARIDS(0993:1008) /'STCNBSDATA      '/                         63
      DATA PARIDS(1009:1024) /'STCCONUSEDATA   '/                         64
      DATA PARIDS(1025:1040) /'STCOTHERDIVDATA '/                         65
      DATA PARIDS(1041:1056) /'STCOTHERSUPDATA '/                         66
      DATA PARIDS(1057:1072) /'STCRUNOFFDATA   '/                         67
      DATA PARIDS(1073:1088) /'STCPRECIPDATA   '/                         68
      DATA PARIDS(1089:1104) /'STCEVAPDATA     '/                         69
      DATA PARIDS(1105:1120) /'STCEOPLEVELDATA '/                         70
      DATA PARIDS(1121:1136) /'ERISTARTLEVEL   '/                         71
      DATA PARIDS(1137:1152) /'ERINBSDATA      '/                         72
      DATA PARIDS(1153:1168) /'ERICONUSEDATA   '/                         73
      DATA PARIDS(1169:1184) /'WELLANDCANALDATA'/                         74
      DATA PARIDS(1185:1200) /'NYSBCANALDATA   '/                         75
      DATA PARIDS(1201:1216) /'ERIOTHERDIVDATA '/                         76
      DATA PARIDS(1217:1232) /'ERIOTHERSUPDATA '/                         77
      DATA PARIDS(1233:1248) /'ERIRUNOFFDATA   '/                         78
      DATA PARIDS(1249:1264) /'ERIPRECIPDATA   '/                         79
      DATA PARIDS(1265:1280) /'ERIEVAPDATA     '/                         80
      DATA PARIDS(1281:1296) /'ERIEOPLEVELDATA '/                         81
      DATA PARIDS(1297:1312) /'CGIPSTARTLEVEL  '/                         82
      DATA PARIDS(1313:1328) /'CGIPAVGLEVELDATA'/                         83
      DATA PARIDS(1329:1344) /'STMRIVFLOWDATA  '/                         84
      DATA PARIDS(1345:1360) /'STMRIVICEDATA   '/                         85
      DATA PARIDS(1361:1376) /'STMRIVWEEDDATA  '/                         86
      DATA PARIDS(1377:1392) /'STMRIVRETRDATA  '/                         87
      DATA PARIDS(1393:1408) /'????????????????'/                         88
      DATA PARIDS(1409:1424) /'STCRIVFLOWDATA  '/                         89
      DATA PARIDS(1425:1440) /'STCRIVICEDATA   '/                         90
      DATA PARIDS(1441:1456) /'STCRIVWEEDDATA  '/                         91
      DATA PARIDS(1457:1472) /'STCRIVRETRDATA  '/                         92
      DATA PARIDS(1473:1488) /'????????????????'/                         93
      DATA PARIDS(1489:1504) /'DETRIVFLOWDATA  '/                         94
      DATA PARIDS(1505:1520) /'DETRIVICEDATA   '/                         95
      DATA PARIDS(1521:1536) /'DETRIVWEEDDATA  '/                         96
      DATA PARIDS(1537:1552) /'DETRIVRETRDATA  '/                         97
      DATA PARIDS(1553:1568) /'????????????????'/                         98
      DATA PARIDS(1569:1584) /'NIARIVFLOWDATA  '/                         99
      DATA PARIDS(1585:1600) /'NIARIVICEDATA   '/                        100
      DATA PARIDS(1601:1616) /'NIARIVWEEDDATA  '/                        101
      DATA PARIDS(1617:1632) /'NIARIVRETRDATA  '/                        102
      DATA PARIDS(1633:1648) /'????????????????'/                        103
      DATA PARIDS(1649:1664) /'ONTSTARTLEVEL   '/                        104
      DATA PARIDS(1665:1680) /'ONTGATESDATA    '/                        105
      DATA PARIDS(1681:1696) /'ONTNBSDATA      '/                        106
      DATA PARIDS(1697:1712) /'ONTCONUSEDATA   '/                        107
      DATA PARIDS(1713:1728) /'ONTOTHERDIVDATA '/                        108
      DATA PARIDS(1729:1744) /'ONTOTHERSUPDATA '/                        109
      DATA PARIDS(1745:1760) /'ONTRUNOFFDATA   '/                        110
      DATA PARIDS(1761:1776) /'ONTPRECIPDATA   '/                        111
      DATA PARIDS(1777:1792) /'ONTEVAPDATA     '/                        112
      DATA PARIDS(1793:1808) /'ONTEOPLEVELDATA '/                        113
      DATA PARIDS(1809:1824) /'ONTDEVIATIONDATA'/                        114
      DATA PARIDS(1825:1840) /'ONTFORCEFLOWDATA'/                        115
      DATA PARIDS(1841:1856) /'SUPAREAELEV     '/                        116
      DATA PARIDS(1857:1872) /'MHUAREAELEV     '/                        117
      DATA PARIDS(1873:1888) /'MICAREAELEV     '/                        118
      DATA PARIDS(1889:1904) /'HURAREAELEV     '/                        119
      DATA PARIDS(1905:1920) /'STCAREAELEV     '/                        120
      DATA PARIDS(1921:1936) /'ERIAREAELEV     '/                        121
      DATA PARIDS(1937:1952) /'ONTAREAELEV     '/                        122
      DATA PARIDS(1953:1968) /'FORECASTMONTHS  '/                        123
      DATA PARIDS(1969:1984) /'SUPOUTPUTFILES  '/                        124
      DATA PARIDS(1985:2000) /'MHUOUTPUTFILES  '/                        125
      DATA PARIDS(2001:2016) /'STCOUTPUTFILES  '/                        126
      DATA PARIDS(2017:2032) /'ERIOUTPUTFILES  '/                        127
      DATA PARIDS(2033:2048) /'ONTOUTPUTFILES  '/                        128
C     DATA PARIDS(2049:2064) /'DPMIFLOWDATA    '/                        129
      DATA PARIDS(2049:2064) /'DPRRIVFLOWDATA  '/                        129
      DATA PARIDS(2065:2080) /'LSLFLOWDATA     '/                        130
      DATA PARIDS(2081:2096) /'STLRIVFLOWDATA  '/                        131
      DATA PARIDS(2097:2112) /'STLRIVPPRETRDATA'/                        132
      DATA PARIDS(2113:2128) /'STLRIVPPICEDATA '/                        133
      DATA PARIDS(2129:2144) /'STLRIVPPWEEDDATA'/                        134
      DATA PARIDS(2145:2160) /'LSLRETRDATA '/                            135
      DATA PARIDS(2161:2176) /'LSLROUGHNESSDATA'/                        136
      DATA PARIDS(2177:2192) /'LSLONTDIFFDATA  '/                        137
      DATA PARIDS(2193:2208) /'MIIRIVFLOWDATA  '/                        138
      DATA PARIDS(2209:2224) /'STMRIVSOLUTION  '/                        139
      DATA PARIDS(2225:2240) /'SUPINC          '/                        140
      DATA PARIDS(2241:2256) /'SUPTIMESTEP     '/                        141
      DATA PARIDS(2257:2272) /'STMRIVEQN       '/                        142
      DATA PARIDS(2273:2288) /'PREVSTMRIVFLOW  '/                        143
      DATA PARIDS(2289:2304) /'PREVSUPGATESOPEN'/                        144
      DATA PARIDS(2305:2320) /'P77ACOMPATMODE  '/                        145
      DATA PARIDS(2321:2336) /'????????????????'/                        146
      DATA PARIDS(2337:2352) /'BOMROUND        '/                        147
      DATA PARIDS(2353:2368) /'FLOWROUND       '/                        148
      DATA PARIDS(2369:2384) /'LEVELROUND      '/                        149
      DATA PARIDS(2385:2400) /'ROUTINGTITLE    '/                        150
      DATA PARIDS(2401:2416) /'STCRIVEQN       '/                        151
      DATA PARIDS(2417:2432) /'DETRIVEQN       '/                        152
      DATA PARIDS(2433:2448) /'NIARIVEQN       '/                        153
      DATA PARIDS(2449:2464) /'MIDLAKESOLUTION '/                        154
      DATA PARIDS(2465:2480) /'MIDLAKEINC      '/                        155
      DATA PARIDS(2481:2496) /'MIDLAKETIMESTEP '/                        156
      DATA PARIDS(2497:2512) /'SUPSUMINPUTFILE '/                        157
      DATA PARIDS(2513:2528) /'MHUSUMINPUTFILE '/                        158
      DATA PARIDS(2529:2544) /'STCSUMINPUTFILE '/                        159
      DATA PARIDS(2545:2560) /'ERISUMINPUTFILE '/                        160
      DATA PARIDS(2561:2576) /'ONTSUMINPUTFILE '/                        161
      DATA PARIDS(2577:2592) /'WRITEERIEOUTFLOW'/                        162
      DATA PARIDS(2593:2608) /'USTLRIVSOLUTION '/                        163
      DATA PARIDS(2609:2624) /'ONTINC          '/                        164
      DATA PARIDS(2625:2640) /'ONTTIMESTEP     '/                        165
      DATA PARIDS(2641:2656) /'USTLRIVEQN      '/                        166
      DATA PARIDS(2657:2672) /'SUPINDEX1DATA   '/                        167
      DATA PARIDS(2673:2688) /'SUPINDEX2DATA   '/                        168
      DATA PARIDS(2689:2704) /'SUPINDEX3DATA   '/                        169
      DATA PARIDS(2705:2720) /'P77ASUPTARGETLEV'/                        170
      DATA PARIDS(2721:2736) /'P77AMHUTARGETLEV'/                        171
      DATA PARIDS(2737:2752) /'P77ASUPRANGEFACT'/                        172
      DATA PARIDS(2753:2768) /'P77AMHURANGEFACT'/                        173
      DATA PARIDS(2769:2784) /'P77ABASEOUTFLOW '/                        174
      DATA PARIDS(2785:2800) /'P77AMINIMUMFLOW '/                        175
      DATA PARIDS(2801:2816) /'P77ALLOGOKIDIV  '/                        176
      DATA PARIDS(2817:2832) /'P77ACHICAGODIV  '/                        177
      DATA PARIDS(2833:2848) /'P77AWELLANDDIV  '/                        178
      DATA PARIDS(2849:2864) /'P77ASTCRIVRETR  '/                        179
      DATA PARIDS(2865:2880) /'P77ADETRIVRETR  '/                        180
      DATA PARIDS(2881:2896) /'P77ANIARIVRETR  '/                        181
      DATA PARIDS(2897:2912) /'P77ASUP5NBS     '/                        182
      DATA PARIDS(2913:2928) /'P77ASUP50NBS    '/                        183
      DATA PARIDS(2929:2944) /'P77ASUP95NBS    '/                        184
      DATA PARIDS(2945:2960) /'P77AMHU5NBS     '/                        185
      DATA PARIDS(2961:2976) /'P77AMHU50NBS    '/                        186
      DATA PARIDS(2977:2992) /'P77AMHU95NBS    '/                        187
      DATA PARIDS(2993:3008) /'P77ASTC50NBS    '/                        188
      DATA PARIDS(3009:3024) /'P77AERI50NBS    '/                        189
      DATA PARIDS(3025:3040) /'P77ANIARIVEQN   '/                        190
      DATA PARIDS(3041:3056) /'P77ADETRIVEQN   '/                        191
      DATA PARIDS(3057:3072) /'P77ASTCRIVEQN   '/                        192
      DATA PARIDS(3073:3088) /'P77ASTMRIVEQN   '/                        193
      DATA PARIDS(3089:3104) /'P77AMIDLAKEINC  '/                        194
      DATA PARIDS(3105:3120) /'P77AMIDLAKETOL  '/                        195
      DATA PARIDS(3121:3136) /'P77ACGIPLEVEL   '/                        196
      DATA PARIDS(3137:3152) /'P77AMONTHROUND  '/                        197
      DATA PARIDS(3153:3168) /'P77ALEVELROUND  '/                        198
      DATA PARIDS(3169:3184) /'P77AFLOWROUND   '/                        199
      DATA PARIDS(3185:3200) /'P77AMONTHSTOAVG '/                        200
      DATA PARIDS(3201:3216) /'P77AACONSTANT   '/                        201
      DATA PARIDS(3217:3232) /'P77AMAXCHANGE   '/                        202
      DATA PARIDS(3233:3248) /'P77ALOWTRIGGER  '/                        203
      DATA PARIDS(3249:3264) /'P77AHIGHTRIGGER '/                        204
      DATA PARIDS(3265:3280) /'SUPSWPIERFORMULA'/                        205
      DATA PARIDS(3281:3296) /'SUPCRITERIAROUND'/                        206
      DATA PARIDS(3297:3312) /'SUPBOMLEVEL     '/                        207
      DATA PARIDS(3313:3328) /'MHUBOMLEVEL     '/                        208
      DATA PARIDS(3329:3344) /'STCBOMLEVEL     '/                        209
      DATA PARIDS(3345:3360) /'ERIBOMLEVEL     '/                        210
      DATA PARIDS(3361:3376) /'SUPPREPROJMETHOD'/                        211
      DATA PARIDS(3377:3392) /'SUPCRITBMETHOD  '/                        212
      DATA PARIDS(3393:3408) /'SUPCRITCTHRESH  '/                        213
      DATA PARIDS(3409:3424) /'COMPUTEMONTREAL '/                        214
      DATA PARIDS(3425:3440) /'SEAWAYALERTDEPTH'/                        215
      DATA PARIDS(3441:3456) /'P58DELIMINATEDEV'/                        216
      DATA PARIDS(3457:3472) /'SEAWAYFORCEFLOW '/                        217
      DATA PARIDS(3473:3488) /'P58DTOTALSUPPLY '/                        218
      DATA PARIDS(3489:3504) /'P58DUSELACTLEV  '/                        219
      DATA PARIDS(3505:3520) /'P58DPREPROJCOND '/                        220
      DATA PARIDS(3521:3536) /'P58DMAXILIMIT   '/                        221
      DATA PARIDS(3537:3552) /'P58DMINPLIMIT   '/                        222
      DATA PARIDS(3553:3568) /'P58DACTUALLEVEL '/                        223
      DATA PARIDS(3569:3584) /'P58DKOFACTOR    '/                        224
      DATA PARIDS(3585:3600) /'CRITERIONKFILE  '/                        225
      DATA PARIDS(3601:3616) /'P58DAPPLIMIT    '/                        226
      DATA PARIDS(3617:3632) /'P58DSI          '/                        227
      DATA PARIDS(3633:3648) /'P58DRUNTYPE     '/                        228
      DATA PARIDS(3649:3664) /'P58DFORECASTDATE'/                        229
      DATA PARIDS(3665:3680) /'????????????????'/                        230
      DATA PARIDS(3681:3696) /'P58DAPPLIMITSYM '/                        231
      DATA PARIDS(3697:3712) /'P58DPLANLEVEL   '/                        232
      DATA PARIDS(3713:3728) /'P58DTOTALDEV    '/                        233
      DATA PARIDS(3729:3744) /'P58DADJUSTMENTSI'/                        234
      DATA PARIDS(3745:3760) /'P58DADJUSTEDSI  '/                        235 
      DATA PARIDS(3761:3776) /'P58DACTUALFLOW  '/                        236 
      DATA PARIDS(3777:3792) /'P58DUSEDATABASE '/                        237
      DATA PARIDS(3793:3808) /'P58DSIDATABASE  '/                        238
      DATA PARIDS(3809:3824) /'PREPROJECTLEVEL '/                        239
      DATA PARIDS(3825:3840) /'P58DSIDATAFILE  '/                        240
      DATA PARIDS(3841:3856) /'P58DLSLFLOW     '/                        241
      DATA PARIDS(3857:3872) /'FLOWDPMIDATA    '/                        242
      DATA PARIDS(3873:3888) /'FLOWSTFRANDATA  '/                        243
      DATA PARIDS(3889:3904) /'FLOWSTMAUDATA   '/                        244
      DATA PARIDS(3905:3920) /'FLOWRICHDATA    '/                        245
      DATA PARIDS(3921:3936) /'ROUGHJETTY1DATA '/                        246 
      DATA PARIDS(3937:3952) /'ROUGHVARENDATA  '/                        247
      DATA PARIDS(3953:3968) /'ROUGHSORELDATA  '/                        248
      DATA PARIDS(3969:3984) /'ROUGHPIERREDATA '/                        249
      DATA PARIDS(3985:4000) /'ROUGH3RIVERSDATA'/                        250
      DATA PARIDS(4001:4016) /'ROUGHBATISDATA  '/                        251
      DATA PARIDS(4017:4032) /'TIDALDATA       '/                        252
      DATA PARIDS(4033:4048) /'ONTTOTALSUPDATA '/                        253
      DATA PARIDS(4049:4064) /'OTHER06ROUGHDATA'/                        254
      DATA PARIDS(4065:4080) /'OTHER06FLOWDATA '/                        255
      DATA PARIDS(4081:4096) /'????????????????'/                        256
      DATA PARIDS(4097:4112) /'????????????????'/                        257
      DATA PARIDS(4113:4128) /'USSICEWEEDADJ   '/                        258
      DATA PARIDS(4129:4144) /'USSLIPEQUATION  '/                        259
      DATA PARIDS(4145:4160) /'SIDECHANNELMIN  '/                        260
      DATA PARIDS(4161:4176) /'SIDECHANNELMAX  '/                        261 
      DATA PARIDS(4177:4192) /'COMPWORKSOVERTOP'/                        262 
      DATA PARIDS(4193:4208) /'COMPWORKSFREEBRD'/                        263 
      DATA PARIDS(4209:4224) /'????????????????'/                        264
      DATA PARIDS(4225:4240) /'????????????????'/                        265
      DATA PARIDS(4241:4256) /'LASTSTURGEONYEAR'/                        266
      DATA PARIDS(4257:4272) /'P2012SPHIGHLEVEL'/                        267
      DATA PARIDS(4273:4288) /'FORCESUPSIDECHAN'/                        268
      DATA PARIDS(4289:4304) /'FORCESUPGATES   '/                        269
      DATA PARIDS(4305:4320) /'INDEXFORCEDGS   '/                        270
      DATA PARIDS(4321:4336) /'????????????????'/                        271
      DATA PARIDS(4337:4352) /'????????????????'/                        272
      DATA PARIDS(4353:4368) /'????????????????'/                        273
      DATA PARIDS(4369:4384) /'P2012SPTARGET   '/                        274
      DATA PARIDS(4385:4400) /'P2012MHTARGET   '/                        275
      DATA PARIDS(4401:4416) /'????????????????'/                        276
      DATA PARIDS(4417:4432) /'????????????????'/                        277
      DATA PARIDS(4433:4448) /'????????????????'/                        278
      DATA PARIDS(4449:4464) /'????????????????'/                        279
      DATA PARIDS(4465:4480) /'????????????????'/                        280
      DATA PARIDS(4481:4496) /'????????????????'/                        281
      DATA PARIDS(4497:4512) /'????????????????'/                        282
      DATA PARIDS(4513:4528) /'????????????????'/                        283
      DATA PARIDS(4529:4544) /'????????????????'/                        284
      DATA PARIDS(4545:4560) /'????????????????'/                        285
      DATA PARIDS(4561:4576) /'????????????????'/                        286
      DATA PARIDS(4577:4592) /'????????????????'/                        287
      DATA PARIDS(4593:4608) /'????????????????'/                        288
      DATA PARIDS(4609:4624) /'????????????????'/                        289
      DATA PARIDS(4625:4640) /'????????????????'/                        290
      DATA PARIDS(4641:4656) /'????????????????'/                        291
      DATA PARIDS(4657:4672) /'????????????????'/                        292
      DATA PARIDS(4673:4688) /'????????????????'/                        293
      DATA PARIDS(4689:4704) /'????????????????'/                        294
      DATA PARIDS(4705:4720) /'P2012PARA       '/                        295
      DATA PARIDS(4721:4736) /'????????????????'/                        296
      DATA PARIDS(4737:4752) /'????????????????'/                        297
      DATA PARIDS(4753:4768) /'????????????????'/                        298
      DATA PARIDS(4769:4784) /'????????????????'/                        299
      DATA PARIDS(4785:4800) /'P2012SPLOWLEVEL '/                        300
      DATA PARIDS(4801:4816) /'????????????????'/                        301
      DATA PARIDS(4817:4832) /'????????????????'/                        302
      DATA PARIDS(4833:4848) /'????????????????'/                        303
      DATA PARIDS(4849:4864) /'????????????????'/                        304
      DATA PARIDS(4865:4880) /'P2012MHTARGETADJ'/                        305
      DATA PARIDS(4881:4896) /'????????????????'/                        306
      DATA PARIDS(4897:4912) /'????????????????'/                        307
      DATA PARIDS(4913:4928) /'????????????????'/                        308
      DATA PARIDS(4929:4944) /'????????????????'/                        309
      DATA PARIDS(4945:4960) /'????????????????'/                        310
      DATA PARIDS(4961:4976) /'????????????????'/                        311
      DATA PARIDS(4977:4992) /'????????????????'/                        312
      DATA PARIDS(4993:5008) /'????????????????'/                        313
      DATA PARIDS(5009:5024) /'????????????????'/                        314
      DATA PARIDS(5025:5040) /'????????????????'/                        315
      DATA PARIDS(5041:5056) /'????????????????'/                        316
      DATA PARIDS(5057:5072) /'????????????????'/                        317
      DATA PARIDS(5073:5088) /'????????????????'/                        318
      DATA PARIDS(5089:5104) /'????????????????'/                        319
      DATA PARIDS(5105:5120) /'????????????????'/                        320
      DATA PARIDS(5121:5136) /'????????????????'/                        321
      DATA PARIDS(5137:5152) /'????????????????'/                        322
      DATA PARIDS(5153:5168) /'????????????????'/                        323
      DATA PARIDS(5169:5184) /'????????????????'/                        324
      DATA PARIDS(5185:5200) /'????????????????'/                        325
      DATA PARIDS(5201:5216) /'????????????????'/                        326
      DATA PARIDS(5217:5232) /'????????????????'/                        327
      DATA PARIDS(5233:5248) /'????????????????'/                        328
      DATA PARIDS(5249:5264) /'????????????????'/                        329
      DATA PARIDS(5265:5280) /'????????????????'/                        330
      DATA PARIDS(5281:5296) /'????????????????'/                        331
      DATA PARIDS(5297:5312) /'????????????????'/                        332
      DATA PARIDS(5313:5328) /'????????????????'/                        333
      DATA PARIDS(5329:5344) /'????????????????'/                        334
      DATA PARIDS(5345:5360) /'????????????????'/                        335
      DATA PARIDS(5361:5376) /'????????????????'/                        336
      DATA PARIDS(5377:5392) /'????????????????'/                        337
      DATA PARIDS(5393:5408) /'????????????????'/                        338
      DATA PARIDS(5409:5424) /'????????????????'/                        339
      DATA PARIDS(5425:5440) /'????????????????'/                        340
      DATA PARIDS(5441:5456) /'????????????????'/                        341
      DATA PARIDS(5457:5472) /'????????????????'/                        342
      DATA PARIDS(5473:5488) /'????????????????'/                        343
      DATA PARIDS(5489:5504) /'????????????????'/                        344
      DATA PARIDS(5505:5520) /'????????????????'/                        345
      DATA PARIDS(5521:5536) /'????????????????'/                        346
      DATA PARIDS(5537:5552) /'????????????????'/                        347
      DATA PARIDS(5553:5568) /'????????????????'/                        348
      DATA PARIDS(5569:5584) /'????????????????'/                        349
      DATA PARIDS(5585:5600) /'????????????????'/                        350
      DATA PARIDS(5601:5616) /'????????????????'/                        351
      DATA PARIDS(5617:5632) /'????????????????'/                        352
      DATA PARIDS(5633:5648) /'????????????????'/                        353
      DATA PARIDS(5649:5664) /'????????????????'/                        354
      DATA PARIDS(5665:5680) /'????????????????'/                        355
      DATA PARIDS(5681:5696) /'????????????????'/                        356
      DATA PARIDS(5697:5712) /'????????????????'/                        357
      DATA PARIDS(5713:5728) /'????????????????'/                        358
      DATA PARIDS(5729:5744) /'????????????????'/                        359
      DATA PARIDS(5745:5760) /'????????????????'/                        360
      DATA PARIDS(5761:5776) /'????????????????'/                        361
      DATA PARIDS(5777:5792) /'????????????????'/                        362
      DATA PARIDS(5793:5808) /'????????????????'/                        363
      DATA PARIDS(5809:5824) /'????????????????'/                        364
      DATA PARIDS(5825:5840) /'????????????????'/                        365
      DATA PARIDS(5841:5856) /'????????????????'/                        366
      DATA PARIDS(5857:5872) /'????????????????'/                        367
      DATA PARIDS(5873:5888) /'????????????????'/                        368
      DATA PARIDS(5889:5904) /'????????????????'/                        369
      DATA PARIDS(5905:5920) /'????????????????'/                        370
      DATA PARIDS(5921:5936) /'????????????????'/                        371
      DATA PARIDS(5937:5952) /'????????????????'/                        372
      DATA PARIDS(5953:5968) /'????????????????'/                        373
      DATA PARIDS(5969:5984) /'????????????????'/                        374
      DATA PARIDS(5985:6000) /'????????????????'/                        375
      DATA PARIDS(6001:6016) /'????????????????'/                        376
      DATA PARIDS(6017:6032) /'????????????????'/                        377
      DATA PARIDS(6033:6048) /'????????????????'/                        378
      DATA PARIDS(6049:6064) /'????????????????'/                        379
      DATA PARIDS(6065:6080) /'????????????????'/                        380
      DATA PARIDS(6081:6096) /'????????????????'/                        381
      DATA PARIDS(6097:6112) /'????????????????'/                        382
      DATA PARIDS(6113:6128) /'????????????????'/                        383
      DATA PARIDS(6129:6144) /'????????????????'/                        384
      DATA PARIDS(6145:6160) /'????????????????'/                        385
      DATA PARIDS(6161:6176) /'????????????????'/                        386
      DATA PARIDS(6177:6192) /'????????????????'/                        387
      DATA PARIDS(6193:6208) /'????????????????'/                        388
      DATA PARIDS(6209:6224) /'????????????????'/                        389
      DATA PARIDS(6225:6240) /'????????????????'/                        390
      DATA PARIDS(6241:6256) /'????????????????'/                        391
      DATA PARIDS(6257:6272) /'????????????????'/                        392
      DATA PARIDS(6273:6288) /'????????????????'/                        393
      DATA PARIDS(6289:6304) /'????????????????'/                        394
      DATA PARIDS(6305:6320) /'????????????????'/                        395
      DATA PARIDS(6321:6336) /'????????????????'/                        396
      DATA PARIDS(6337:6352) /'????????????????'/                        397
      DATA PARIDS(6353:6368) /'????????????????'/                        398
      DATA PARIDS(6369:6384) /'P77AWINTERMAX   '/                        399
      DATA PARIDS(6385:6400) /'????????????????'/                        400
C
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE INITPARAMS()
C     -------------------      
C
C     Initialize stuff, read in model parameters, open data files
C
      IMPLICIT NONE
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER I, J
C               These are scratch variables for use in
C               loop control and for intermediate results.
C
      LOGICAL LOGVAR
C               This is a scratch variable used here and there.
C
      CHARACTER CHNKST*16
C               This is the calling subroutine name passed to ErrMess() in
C               case of error messages.
C
      CHARACTER*(*) InitializationLog, MasterInputFile
C               These are the hardcoded names of the log for input parsing 
C               and the master input file. 
C
      CHARACTER INLINE*256
C               This is a scratch variable used for hold lines from the
C               master input file and pass them to MIFLINE().
C
      CHARACTER ComputeLogName*256
      INTEGER   ComputeLogNameLen
C               This is the name of the file for the detailed processing log.
C
      CHARACTER*256 OutFileName
C                   Function tha builds complete filespec.
C     
      INTEGER      LastKC
C                  Function that returns the position of last 'keyboard char'
C 
C     Declarations regarding COMMON block variables:
C
      CHARACTER*256  FileNames(100)
C                   See descriptions in Block Data FileNames
C     
      CHARACTER*256 RepeatLines(80)
C                  See description in Block Data Genral
C
      CHARACTER CTMP*255, CWRK*80
C               These are just buffers used throughout program.
C
      CHARACTER LLT*1
C               This is the character that marks the end of a line sent
C               to the LogIt subroutine.
C
      CHARACTER STIME*8,SDATE*9
C               These hold the starting date and time, used for detail 
C               in logs and outputs
C
      CHARACTER MonthAbbrev(12)*3
C               See description in Block Data Genral
C
      CHARACTER     DirSep*5,OSName(5)*8,MessDB*256
      CHARACTER*64  ModelerName
C               See descriptions in Block Data Genral
C
      LOGICAL   LOGC(8), LOGD(8), LOGS(8), 
     &          LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &          LOGM(8), LOGP(8), LOGO(8), 
     &          LOGDM(8), LOGDP(8), LOGDO(8), 
     &          LOGSM(8), LOGSP(8), LOGSO(8),
     &          LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &          LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C               See descriptions in Block Data Genral
C
      CHARACTER RCOUTDIR*128, RCOUTEXT*16
C               See descriptions in Block Data Genral
C
      LOGICAL   OnlyOneOpen
      INTEGER   NOPDIR,NOPEXT, TSOutputFileNum(4,2,8)
C               See descriptions in Block Data Genral
C
      INTEGER   OSTYPE, NUMPARAM
C               See descriptions in Block Data Genral
C
      INTEGER   UnitTypeLen(20), BaseUnitType(20)
C               Not used in this function; only declared because in 
C               same COMMON block as UnitConvFactor.
C                   
      INTEGER   TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &          MissingVal, SecsInUniformMon
      DOUBLE PRECISION MPFT, UnitConvFactor(20), HALF, ONE, ZERO
C               See descriptions in Block Data Genral
C                   
      CHARACTER*80 RCTITLE(3), MLHEADER
      CHARACTER*256 CritKfile
C
C
C
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /MONAME/   MonthAbbrev
      COMMON /GENRLI/   OSTYPE, NUMPARAM
      COMMON /GENRLC/   OSName,DirSep,MessDB, ModelerName
      COMMON /STIME/    STIME,SDATE
      COMMON /OUTNAM/   RCOUTDIR, RCOUTEXT
      COMMON /OUTNMI/   NOPDIR,NOPEXT, TSOutputFileNum, OnlyOneOpen
C
      COMMON /FNAMEC/   FileNames, RepeatLines
      COMMON /UNITSI/   UnitTypeLen, BaseUnitType, UnitConvFactor
      COMMON /TITLES/   RCTITLE, MLHEADER
      COMMON /CritK/    CritKfile
C
C
      PARAMETER (InitializationLog='init.log')
      PARAMETER (MasterInputFile='cglrrm.ini')
C
C  Initialize some values.
C
C
C  Store starting date and time for later use.
C
      CALL MyTime(STIME)
      CALL MyDate(SDATE)
C
C     DOS '\' character in DirSep(2) initialized here because some compilers
C     think you're escaping the next character in a DATA statement.
      DirSep(2:2)=CHAR(92)
C
C     Line terminator initialized here because I like using the NULL, which
C     also can cause compiler fits if in a DATA or PARAMETER statement.
      LLT=CHAR(0)
C
C     Initialize the units conversion constants:
C     For converting 10m3s to m3s:
      UnitConvFactor(1) = TEN
C     For converting m3s to m3s:
      UnitConvFactor(2) = ONE
C     For converting tcfs to m3s:
      UnitConvFactor(3) = ( MPFT  * TEN )**3
C     For converting cfs to m3s:
      UnitConvFactor(4) = ( MPFT )**3
C     For converting mm to m:
      UnitConvFactor(5) = ( TEN )**3
C     For converting cm to m:
      UnitConvFactor(6) = ( TEN )**2
C     For converting m to m:
      UnitConvFactor(7) = ONE
C     For converting km2 to m2:
      UnitConvFactor(8) = ( TEN )**6
C     For converting m2 to m2:
      UnitConvFactor(9) = ONE
C     For converting m3 to m3:
      UnitConvFactor(10)= ONE
C     For converting inches to m:
      UnitConvFactor(11)= MPFT / INPFT
C     For converting feet to m:
      UnitConvFactor(12)= MPFT
C     For converting ft2 to m2:
      UnitConvFactor(13)= ( MPFT )**2
C     For converting ft3 to m3:
      UnitConvFactor(14)= ( MPFT )**3
C     For converting mi2 to m2:
      UnitConvFactor(15)= ( FTPMI * MPFT )**2
C     For converting Less Than 10 to LT10
      UnitConvFactor(16)= ONE
C     For converting Missing to MISS
      UnitConvFactor(17)= ONE
C     For converting Not Applicable to NA
      UnitConvFactor(18)= ONE
C     For converting None to NONE
      UnitConvFactor(19)= ONE
C     For converting Not Applicable to '  --  '
      UnitConvFactor(20)= ONE
C
C
C     This identifies which subroutine calls for a message.
      CHNKST='INITPARAMS'//LLT
C
C     Initialize time-series input filespecs
      DO 161 J=12,256
C        The next three are initialized to spaces in Block Data Genral
C        P58DSIDataFile(J:J)= LLT
C        P58DSimulateFile(J:J)= LLT
C        CritKfile(J:J)= LLT
  161    MessDB(J:J)= LLT
      DO 163, I=1,100
         DO 162 J=12,256
  162       FileNames(I)(J:J)= LLT
  163    FileNames(I)(1:11) = 'Unspecified'

C
      DO 164, I=1,80
  164    MLHEADER(I:I) = ' '
      DO 165, I=1,3
  165    RCTITLE(I) = MLHEADER
C
C
C     
C  Open the log file for critical errors in current directory.
C
C     OPEN(UNIT=91,FILE=LF91,IOSTAT=I,ACCESS='APPEND') 
C     IF ( I .NE. 0 ) THEN
C        CALL ErrMess(1000,CHNKST,LOGC)
C        WRITE(CTMP,2015)I,LF91,LLT
C        CALL LogIt(LOGC,1,1)
C        CALL BuildInfo
C        CALL LogIt(LOGC,0,0)
C        CALL SAVETOD
C        CALL ADIOS
C     ENDIF
C     
C  Open the Initialization log file in current directory.
C
      OPEN(UNIT=90,FILE=InitializationLog,IOSTAT=I) 
      IF ( I .NE. 0 ) THEN
         CALL ErrMess(1000,CHNKST,LOGCE)
         WRITE(CTMP,2015)I,InitializationLog,LLT
         CALL LogIt(LOGCE,1,1)
         CALL BuildInfo
         CALL LogIt(LOGCE,0,0)
         CALL SAVETOD
         CALL ADIOS
      ENDIF
      WRITE(CTMP,'(A,A1)')' Detailed log of CGLRRM processing.', LLT
      CALL LogIt(LOGD,1,1)
      WRITE(CTMP,'(A,A,A,A,A1,A1)')' Run started at ',STIME,' on ',
     &   SDATE,LLT,'.'
      CALL LogIt(LOGD,1,1)
      CALL BuildInfo
      CALL LogIt(LOGD,0,1)
C     
C  See if master input file exists.  If not, squawk and abort.
C
      INQUIRE(FILE=MasterInputFile,EXIST=LOGVAR)
      IF(.NOT. LOGVAR) THEN
         CALL ErrMess(1001,CHNKST,LOGCDE)
         WRITE(CTMP,'(A,A,A1,A1)')' Expected master input filename: ',
     &      MasterInputFile,'.',LLT
         CALL LogIt(LOGCDE,1,0)
         CALL SAVETOD
         CALL ADIOS
      ENDIF
C     
C  Open the master input file.
C
      OPEN(UNIT=96,FILE=MasterInputFile,STATUS='OLD',IOSTAT=I) 
      IF ( I .NE. 0 ) THEN
         CALL ErrMess(1002,CHNKST,LOGCDE)
         WRITE(CTMP,2015) I,MasterInputFile,LLT
         CALL LogIt(LOGCDE,1,1)
         CALL BuildInfo
         CALL LogIt(LOGCE,0,0)
         CALL SAVETOD
         CALL ADIOS()
      ENDIF
C
C
      DO 1231 J= -32766, 32766, 1
C
C        read in a line at a time, and call processing routine until EOF
C
         READ(96, 2014, END=1233, ERR=1232, IOSTAT=I) INLINE
 1231    CALL MIFLINE(INLINE)
      GOTO 1233
C
 1232 CALL ErrMess(1003,CHNKST,LOGCDE)
      WRITE(CTMP,2015) J,MasterInputFile,LLT
      CALL LogIt(LOGCDE,1,0)
      CALL SAVETOD
      CALL ADIOS
C
 1233 CLOSE(96)
C
C     Build the compute log filespec, now that the output dir and extension 
C     should be known and validated.
      ComputeLogName = OutFileName( 'details', 7, ComputeLogNameLen)
C
C     Now close Unit 90 (the initialization log).  Open Unit 90 again as the 
C     detailed compute log. 
      CLOSE(90)
      OPEN(UNIT=90,FILE=ComputeLogName(1:ComputeLogNameLen),IOSTAT=I) 
      IF ( I .NE. 0 ) THEN
         CALL ErrMess(1000,CHNKST,LOGCE)
         WRITE(CTMP,2015)I,ComputeLogName(1:ComputeLogNameLen),LLT
         CALL LogIt(LOGCE,1,1)
         CALL BuildInfo
         CALL LogIt(LOGCE,0,0)
         CALL SAVETOD
         CALL ADIOS
      ENDIF
      WRITE(CTMP, 2014) 'This file is named '//
     &   ComputeLogName(1:ComputeLogNameLen)//LLT
      CALL LogIt(LOGD,1,1)
C
C     Open Unit 97 as the initialization log, and copy the contents to the new log.
      OPEN(UNIT=97,FILE=InitializationLog,IOSTAT=I) 
      IF ( I .NE. 0 ) THEN
         CALL ErrMess(1000,CHNKST,LOGCE)
         WRITE(CTMP,2015)I,InitializationLog,LLT
         CALL LogIt(LOGCE,1,1)
         CALL BuildInfo
         CALL LogIt(LOGCE,0,0)
         CALL SAVETOD
         CALL ADIOS
      ENDIF
      DO 1235 J= -32766, 32766, 1
         READ(97, 2014, END=1238, ERR=1236, IOSTAT=I) INLINE
 1235    WRITE(90,2014,ERR=1237,IOSTAT=I) InLine(1:LastKC(InLine,256))
      GOTO 1238
C
 1236 CALL ErrMess(1003,CHNKST,LOGCDE)
      WRITE(CTMP,2015) J,InitializationLog,LLT
      CALL LogIt(LOGCDE,1,0)
      CALL SAVETOD
      CALL ADIOS

 1237 CALL ErrMess(1127,CHNKST,LOGCDE)
      WRITE(CTMP,2015) J,ComputeLogName(1:ComputeLogNameLen),LLT
      CALL LogIt(LOGCDE,1,0)
      CALL SAVETOD
      CALL ADIOS
C
 1238 CLOSE(97, STATUS = 'DELETE')
C
 2014 FORMAT(A)
 2015 FORMAT(' Got IOSTAT=',I4,' when accessing file ',A,'.',A1)
C2024 FORMAT(A255)
C
      END
C
C
C
C-----------------------------------------------------------------------
C
      LOGICAL FUNCTION ReadAvsZ( LAKENUM, STRING, SLEN )
C
C     This function processes a surface area input line. 
C
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER STRING*(*)
C                  The character-array to parse.
C
      INTEGER      SLEN
C                  The length of the string.
C
      INTEGER      LAKENUM
C                  The index to the first dimension of array AVSZ, 
C                  corresponding to which lake the data refers.
C
C     Declarations regarding local variables or functions:
C
      INTEGER      I, J, FS, FE
C                  These are scratch variables.
C
      DOUBLE PRECISION AREAVAL, ELEVVAL
C                  The area and elevation values on the line.
C
      CHARACTER    CHNKST*16
C                  This is the calling subroutine name passed to ErrMess() in
C                  case of error messages.
C
      LOGICAL      FIELDBND, GETFLOAT
C                  Functions that helps parse comma-delimited input lines
C                  and floating-point input values.
C
C
C     Declarations regarding COMMON block variables:
C
C
      INTEGER      NUMAVSZ(8)
C                  See descriptions in Block Data Genral
C
      DOUBLE PRECISION AVSZ(8,5,2)
C                  See descriptions in Block Data Genral
C
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C                  See description in Block Data Genral
C
      CHARACTER*21 LAKENAME(16), RiverName(16)
C                  See description in Block Data Genral
C
      INTEGER      NLN(16), LkToSlv2LN(8), NRN(16), LkSlv2RtLoopI(8)
C                  See description in Block Data Genral
C
      CHARACTER    PARCODE*16, INDENT*4
C                  See descriptions in PARLINE().
C
      CHARACTER    CTMP*255, CWRK*80
C                  These are just buffers used throughout program.
C
      CHARACTER    LLT*1
C                  This is the character that marks the end of a line sent
C                  to the LogIt subroutine.
C
C
      COMMON /RCAVSZ/   NUMAVSZ, AVSZ
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /PARSEC/   INDENT, PARCODE
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /LNAMEC/   LAKENAME, RiverName
      COMMON /LNAMEI/   NLN, LkToSlv2LN, NRN, LkSlv2RtLoopI
C
C
      CHNKST='ReadAvsZ'//LLT
C
      ReadAvsZ = .FALSE.
C
C     Get field bounds of expected area value
      IF( .NOT. FIELDBND(STRING, SLEN, 1, FS, FE, I ) ) THEN 
         CALL ErrMess(1039,CHNKST,LOGCDE)
         GOTO 9999
      ENDIF 
C     Save the beginning of the next field, if any.
      I = I + 1
C
C     Read/convert the value and store in AREAVAL
      IF(.NOT. GETFLOAT( STRING(FS:FE),FE-FS+1,9,AREAVAL)) GOTO 9999
C
C     Get field bounds of expected elevation value.
      IF( FIELDBND(STRING, SLEN, I, FS, FE, J ) ) THEN 
C    
C        Looks like we're specifying a point.
C        Read/convert the value and store in ELEVVAL
         IF(.NOT. GETFLOAT(STRING(FS:FE),FE-FS+1,7,ELEVVAL))GOTO 9999
C     
C        See if contradicts earlier parameter calling for constant area.
C        The first elevation will be -1.0 is so.
         IF(NUMAVSZ(LAKENUM).EQ.1 .AND. AVSZ(LAKENUM,1,1).LT. 0.D0)THEN
            CALL ErrMess(1047,CHNKST,LOGCDE)
            GOTO 9999
         ENDIF
C     
C        Initialize count if this is first area record (will be MissingVal)
         IF(NUMAVSZ(LAKENUM) .LT. 0) NUMAVSZ(LAKENUM) = 0
         NUMAVSZ(LAKENUM) = NUMAVSZ(LAKENUM) + 1
C     
C        Make sure not too many points.
         IF(NUMAVSZ(LAKENUM) .GT. 5) THEN
            CALL ErrMess(1045,CHNKST,LOGCDE)
            GOTO 9999
         ENDIF
C
C        assign the values.  they will get sorted after init() closes.
         AVSZ(LAKENUM, NUMAVSZ(LAKENUM), 1) = ELEVVAL
         AVSZ(LAKENUM, NUMAVSZ(LAKENUM), 2) = AREAVAL
C
         WRITE(CTMP,9995)INDENT,PARCODE,
     &      LAKENAME(LAKENUM)(1:NLN(LAKENUM)), NUMAVSZ(LAKENUM), 
     &      AVSZ(LAKENUM, NUMAVSZ(LAKENUM), 2),AVSZ(LAKENUM, 
     &      NUMAVSZ(LAKENUM), 1), LLT     
      ELSE
C        found empty elevation field, or no data.  assume constant area.
C        See if contradicts earlier parameter calling for area vs Z curve
         IF(NUMAVSZ(LAKENUM) .GE.1 .AND. AVSZ(LAKENUM,1,1).GT. 0.D0)THEN
            CALL ErrMess(1048,CHNKST,LOGCDE)
            GOTO 9999
         ENDIF
         NUMAVSZ(LAKENUM) = 1
         AVSZ(LAKENUM,1,1) = -1.D0
         AVSZ(LAKENUM,1,2) = AREAVAL
C
         WRITE(CTMP,9994)INDENT,PARCODE,
     &      LAKENAME(LAKENUM)(1:NLN(LAKENUM)),AVSZ(LAKENUM,1,2),LLT     
      ENDIF 
C
      ReadAvsZ = .TRUE.
C
 9994 FORMAT(A4,A16,': read ',A,' constant area of ',F18.3,' m2.',A1)
     &   
 9995 FORMAT(A4,A16,': read ',A,' area-elevation pair number ',I1,
     &   ' as  ',F16.4,' m2 at ',F7.3,' m.',A1)
 9999 RETURN
      END
C
C-----------------------------------------------------------------------
C
      LOGICAL FUNCTION GETFLOAT( STRING, STRLEN, UNITNUM, OUTFLOAT)
C
C     This function attempts to:
C         1) determine which units are included at the end of STRING
C         2) read the preceding floating-point number into OUTFLOAT
C         3) convert (if necessary) the value in OUTFLOAT to the units
C               called for by UNITNUM.
C
C     It returns .TRUE. if all the above goes successfully.
C
C                   
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER    STRING*(*)
C                  The character-array to parse.
C
      INTEGER      STRLEN
C                  The length of the array.  
C
      INTEGER      UNITNUM
C                  The index to array UnitType that corresponds to the
C                  desired units to return the float in.  
C
      DOUBLE PRECISION OUTFLOAT
C                  The value of the floating point number read in, after
C                  conversion to the proper units.
C
C      
C     Declarations regarding local variables or functions:
C
      INTEGER      I, J, K, L
C                  These are scratch variables.
C
      CHARACTER    CHNKST*16
C                  This is the calling subroutine name passed to ErrMess() in
C                  case of error messages.
C
      INTEGER      LASTKC
C                  Function returning last non-whitespace "keyboard" 
C                  character in a string.
C
C
C     Declarations regarding COMMON block variables:
C
      CHARACTER*6  UnitType(20)
C                  See description in Block Data Genral
C                   
      INTEGER      UnitTypeLen(20)
C                  See description in Block Data Genral
C                   
      INTEGER      BaseUnitType(20)
C                  See description in Block Data Genral
C                   
      INTEGER      TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &             MissingVal, SecsInUniformMon
      DOUBLE PRECISION MPFT, UnitConvFactor(20), HALF, ONE, ZERO
C                  See descriptions in Block Data Genral
C
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C                  See description in Block Data Genral
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                  See description in Block Data Genral
C
      CHARACTER    CTMP*255, CWRK*80
C                  These are just buffers used throughout program.
C
      CHARACTER    LLT*1
C                  This is the character that marks the end of a line sent
C                  to the LogIt subroutine.
C
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /UNITSC/   UnitType
      COMMON /UNITSI/   UnitTypeLen, BaseUnitType, UnitConvFactor
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
C
C
      CHNKST='GETFLOAT'//LLT
C
      GETFLOAT = .FALSE.
C
C
      IF( STRLEN .GE. 32) THEN
         CALL ErrMess(1040,CHNKST,LOGCDE)
         WRITE(CTMP,940)' Unable to read float from:',STRING,LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
C     Now copy to CWRK, making sure everything is uppercase:
C
      DO 40, J = 1, STRLEN
         I = ICHAR( STRING(J:J) )
         IF( I .GE. 97 .AND. I .LE. 122) THEN
            CWRK(J:J) = CHAR(I-32)
         ELSE
            CWRK(J:J) = STRING(J:J)
         ENDIF
 40      CONTINUE
C
C     First determine what units accompany the data (if any).  Search
C     through the list of recognized units
C     
      IF(UNITNUM .EQ. 0)THEN
C        No units to read.  Position J to the space after the number ends:
         J=STRLEN + 1
         I=0
      ELSE
C        Note that only expecting units for data types up to 16.  
C        Shouldn't be specifying a conversion for screwy ones.
         DO 311, I=1,16
            IF( STRLEN .GE. UnitTypeLen(I) )THEN
               J = STRLEN - UnitTypeLen(I) + 1
               IF(CWRK(J:STRLEN) .EQ. UnitType(I)(1:UnitTypeLen(I)))THEN
C                  WRITE(*,'(6X,A,A,A)')CWRK(J:),' matches ',
C     &               UnitType(I)(1:UnitTypeLen(I))
                  GOTO 312
C              ELSE
C                  WRITE(*,'(6X,A,A,A)')CWRK(J:),' does not match ',
C     &               UnitType(I)(1:UnitTypeLen(I))
               ENDIF
            ENDIF
  311       J = 0
  312    IF( J .EQ. 0 )THEN
C           Didn't match anything in UnitType!
            CALL ErrMess(1041,CHNKST,LOGCDE)
            WRITE(CTMP,940)
     &         ' Unable to recognize units in:',CWRK(1:STRLEN),LLT
            CALL LogIt(LOGCDE,0,0)
            RETURN
         ENDIF
         IF( J .LE. 1 )THEN
C           The only thing in the field was the units!
            CALL ErrMess(1042,CHNKST,LOGCDE)
            WRITE(CTMP,940)
     &         ' Unable to find a number in:',CWRK(1:STRLEN),LLT
            CALL LogIt(LOGCDE,0,0)
            RETURN
         ENDIF
      ENDIF
C
C
C     So now I contains the index of the given units corresponding to
C     the list of units in UnitType, and J contains the position of the
C     first character of the units label within STRING.  Next try to read
C     the number.  It seems that the two compilers I'm testing with have 
C     slightly varying nuances regarding interpretation of floating point
C     numbers using the same format specification.  But they both
C     do things the same way if a decimal point is found.  So, let's
C     make the string as readable into REAL as possible.
C     Get position just before start of units label
      K = J - 1
C     Trim any trailing whitespace that may confuse the READ function
      J = LASTKC(CWRK(1:K), K)
C
C     If units are 10m3s, check to see if some separation between number
C     and units.  For instance, 510m3s may mean 5 10m3s (50 m3s), or 
C     possibly 510 m3s.
      IF(I .EQ. 1 .AND. J .EQ. K)THEN
         CALL ErrMess(1046,CHNKST,LOGCDE)
         RETURN
      ENDIF
C
C     See if there is a '.' specified.  If not, supply one.
      IF(INDEX(CWRK(1:J), '.') .EQ. 0 ) THEN
C        No '.' found.  See if exponential form.
         K = INDEX(CWRK(1:J), 'E')
         IF( K .EQ. 0 ) THEN
            K = INDEX(CWRK(1:J), 'D')
         ENDIF
         IF( K .EQ. 0 ) THEN
C           Not exponent format, so add one after digits
            CWRK(J+1:J+1) = '.'
         ELSE
C           Exponent format, so shift everything from 'E' or 'D' to right
C           one position and add '.' where 'E' or 'D' used to be
            DO 333, L=J, K, -1
  333          CWRK(L+1:L+1) = CWRK(L:L)
            CWRK(K:K) = '.'
         ENDIF
         J = J + 1
      ENDIF
C
      READ( CWRK(1:J), '(D62.31)', ERR=2015, IOSTAT=K ) OUTFLOAT
      IF(GenVerbosity .GE. 8) THEN
         WRITE(CTMP, '(A,D21.9,A,A,A,A1)') ' Read a value of ',
     &      OUTFLOAT,' from string ::',CWRK(1:J),'::',LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
C      
      IF( UNITNUM .EQ. 0) THEN     
C        All done.  leave now.
         GETFLOAT = .TRUE. 
         RETURN
      ENDIF
C
C     Next see if any conversion from the given units is required:
C             Index  Units   base unit
C               0     none    -   -
C               1    10m3s    2 (m3s)    
C               2      m3s    2   - 
C               3     tcfs    2 (m3s)
C               4      cfs    2 (m3s)
C               5       mm    7  (m) 
C               6       cm    7  (m)
C               7        m    7   -
C               8      km2    9 (m2) 
C               9       m2    9   -
C              10       m3   10   -
C              11   inches    7  (m)
C              12     feet    7  (m)
C              13      ft2    9 (m2)
C              14      ft3   10 (m3)
C              15      mi2    9 (m2)
C              16     LT10   16 (LT10)
C              17     MISS   17 (MISS)
C              18     NONE   18 (NONE)
C              19       NA   19 (NA)
C              20     --     20 (  --  )
C               
C     See if the base unit for the units given match up with what we expect
C
      IF( BaseUnitType(I) .NE. UNITNUM ) THEN
         CALL ErrMess(1044,CHNKST,LOGCDE)
         WRITE(CTMP,942)' Unable to convert ',
     &      UnitType(I)(1:UnitTypeLen(I)),' into ',
     &      UnitType(UNITNUM)(1:UnitTypeLen(UNITNUM)),
     &      '.  Input field: ', STRING(1:STRLEN),LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
      IF( I .NE. UNITNUM ) THEN
C        Conversion is necessary.  UnitConvFactor(I) contains conversion factor 
C        from UnitType(I) to the base units for UnitType(I).
         OUTFLOAT = OUTFLOAT * UnitConvFactor(I)
C
C        Yet another place where a case statement would be nice.  ugh.
C         IF( I .EQ. 1 ) THEN
C           Convert 10m3s to m3s
C            OUTFLOAT = OUTFLOAT * TEN        
C         ELSE IF( I .EQ. 3 ) THEN
C           Convert tcfs to m3s
C            OUTFLOAT = OUTFLOAT * ( MPFT * TEN )**3
C         ELSE IF( I .EQ. 4 ) THEN
C           Convert cfs to m3s
C            OUTFLOAT = OUTFLOAT * ( MPFT )**3
C         ELSE IF( I .EQ. 5 ) THEN
C           Convert mm to m
C            OUTFLOAT = OUTFLOAT * ( TEN )**3
C         ELSE IF( I .EQ. 6 ) THEN
C           Convert cm to m
C            OUTFLOAT = OUTFLOAT * ( TEN )**2
C         ELSE IF( I .EQ. 8 ) THEN
C           Convert km2 to m2
C            OUTFLOAT = OUTFLOAT * ( TEN )**6
C         ELSE IF( I .EQ. 11 ) THEN
C           Convert inches to m
C            OUTFLOAT = OUTFLOAT * ( MPFT ) / INPFT
C         ELSE IF( I .EQ. 12 ) THEN
C           Convert feet to m
C            OUTFLOAT = OUTFLOAT * ( MPFT )
C         ELSE IF( I .EQ. 13 ) THEN
C           Convert ft2 to m2
C            OUTFLOAT = OUTFLOAT * ( MPFT )**2
C         ELSE IF( I .EQ. 14 ) THEN
C           Convert ft3 to m3
C            OUTFLOAT = OUTFLOAT * ( MPFT )**3
C         ELSE IF( I .EQ. 15 ) THEN
C           Convert mi2 to m2
C            OUTFLOAT = OUTFLOAT * ( FTPMI * MPFT )**2
C         ELSE
C            STOP 'You cant get here from there!'
C         ENDIF
      ENDIF
C
C  
      IF(GenVerbosity .GE. 2) THEN
         WRITE(CTMP, '(A,A,1X,A,A,G25.16,1X,A,A1)') ' Converted ',
     &      CWRK(1:J), UnitType(I)(1:UnitTypeLen(I)), ' into ',OUTFLOAT,
     &     UnitType(BaseUnitType(I))(1:UnitTypeLen(BaseUnitType(I))),LLT
         IF( I .NE. UNITNUM .AND. GenVerbosity .GE. 8) THEN
            CALL LogIt(LOGC,0,0)
         ENDIF
         CALL LogIt(LOGD,0,0)
      ENDIF
C
      GETFLOAT = .TRUE. 
      RETURN
C
 2015 CALL ErrMess(1043,CHNKST,LOGCDE)
      WRITE(CTMP,941)' Unable to read "', CWRK(1:J),
     &   '" according to format (F62.31).  Got IOSTAT=',K,'.',LLT
      CALL LogIt(LOGCDE,0,0)
      CALL BuildInfo
      CALL LogIt(LOGCDE,0,0)
      RETURN
C
  940 FORMAT(A,A,A1) 
  941 FORMAT(A,A,A,I6,A,A1) 
  942 FORMAT(A,A,A,A,A,A,A1) 
      END
C
C
C
C-----------------------------------------------------------------------
C
      LOGICAL FUNCTION FIELDBND( STRING, STRLEN, FLDST,
     &                           FNWKC, LNWKC, FLDEN )
C
C     This function takes a string, its length, and the starting position 
C     of a filed within the string, and returns the positions of the next 
C     non-whitespace "keyboard" character in the field, the last 
C     non-whitespace "keyboard" character in the field, and the postion
C     of the next comma (field delimiter).  For example, given
C             123456789012345678901234567890     (this for easy reference)
C     STRING="RCSTDATE:  1918, 12  ,18" and STRLEN=24
C     With STPOS<=0, function returns 'F', with FNWKC=-1, LNWKC=-1, FLDEN=-1
C     With STPOS= 1, function returns 'T', with FNWKC= 1, LNWKC=15, FLDEN=16
C     With STPOS= 2, function returns 'T', with FNWKC= 2, LNWKC=15, FLDEN=16
C     With STPOS= 3, function returns 'T', with FNWKC= 3, LNWKC=15, FLDEN=16
C     With STPOS= 4, function returns 'T', with FNWKC= 4, LNWKC=15, FLDEN=16
C     With STPOS= 5, function returns 'T', with FNWKC= 5, LNWKC=15, FLDEN=16
C     With STPOS= 6, function returns 'T', with FNWKC= 6, LNWKC=15, FLDEN=16
C     With STPOS= 7, function returns 'T', with FNWKC= 7, LNWKC=15, FLDEN=16
C     With STPOS= 8, function returns 'T', with FNWKC= 8, LNWKC=15, FLDEN=16
C     With STPOS= 9, function returns 'T', with FNWKC= 9, LNWKC=15, FLDEN=16
C     With STPOS=10, function returns 'T', with FNWKC=12, LNWKC=15, FLDEN=16
C     With STPOS=11, function returns 'T', with FNWKC=12, LNWKC=15, FLDEN=16
C     With STPOS=12, function returns 'T', with FNWKC=12, LNWKC=15, FLDEN=16
C     With STPOS=13, function returns 'T', with FNWKC=13, LNWKC=15, FLDEN=16
C     With STPOS=14, function returns 'T', with FNWKC=14, LNWKC=15, FLDEN=16
C     With STPOS=15, function returns 'T', with FNWKC=15, LNWKC=15, FLDEN=16
C     With STPOS=16, function returns 'F', with FNWKC= 0, LNWKC= 0, FLDEN=16
C     With STPOS=17, function returns 'T', with FNWKC=18, LNWKC=19, FLDEN=22
C     With STPOS=18, function returns 'T', with FNWKC=18, LNWKC=19, FLDEN=22
C     With STPOS=19, function returns 'T', with FNWKC=19, LNWKC=19, FLDEN=22
C     With STPOS=20, function returns 'F', with FNWKC= 0, LNWKC= 0, FLDEN=22
C     With STPOS=21, function returns 'F', with FNWKC= 0, LNWKC= 0, FLDEN=22
C     With STPOS=22, function returns 'F', with FNWKC= 0, LNWKC= 0, FLDEN=22
C     With STPOS=23, function returns 'T', with FNWKC=23, LNWKC=24, FLDEN= 0
C     With STPOS=24, function returns 'T', with FNWKC=24, LNWKC=24, FLDEN= 0
C     With STPOS=25, function returns 'F', with FNWKC=-1, LNWKC=-1, FLDEN=-1
C     With STPOS>25, function returns 'F', with FNWKC=-1, LNWKC=-1, FLDEN=-1
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER STRING*(*)
C               Input - The character-array to parse.
C
      INTEGER   STRLEN
C               Input - The length of the array.  
C
      INTEGER   FLDST
C               Input - The first position of the field (must be at least 1).  
C
      INTEGER   FLDEN
C               Output - The last position of the field (location of comma).  
C               If there is no field delimiter (comma) between 1 and
C               STRLEN, then FLDEN is assigned 0.
C
      INTEGER   FNWKC, LNWKC
C               Output - The first and last non-whitespace "keyboard"
C                charactersof the field.  If there are no non-whitespace 
C               "keyboard"characters between FLDST and FLDEN, the values are 0.
C
C     Declarations regarding local variables or functions:
C
C
      INTEGER   NEXTDLM, FIRSTKC, LASTKC, I
C
C     Initialize return value to .FALSE., others to -1:
      FIELDBND = .FALSE.
      FNWKC = -1
      LNWKC = -1
      FLDEN = -1
C
C
C     See if the arguments make sense:
      IF( STRLEN .LE. 0 .OR. FLDST .GT. STRLEN)THEN
C         GOTO 111
         RETURN
      ENDIF
C
C     write(*,'(A,A,/,A,A,A1)')
C    & '  123456789012345678901234567890123456789012345678901234567890',
C    &   '12345678901234567890',' {',
C    &   STRING(1:STRLEN),'}'
C 
C     write(*,'(A,I2,A,I2,A,I2,A,I2,A,I2,A,I2,A,L1)')
C    &   ' 1STRLEN:',STRLEN,'    FLDST:',FLDST,
C    &   '    FNWKC:',FNWKC,'    LNWKC:',LNWKC,'    FLDEN:',FLDEN,
C    &   '    I:',I,'    FIELDBND:', FIELDBND
C
C     Look for next delimiter:
      I = NEXTDLM(STRING(FLDST:STRLEN),STRLEN-FLDST+1)
C
C     write(*,'(A,I2,A,I2,A,I2,A,I2,A,I2,A,I2,A,L1)')
C    &   ' 2STRLEN:',STRLEN,'    FLDST:',FLDST,
C    &   '    FNWKC:',FNWKC,'    LNWKC:',LNWKC,'    FLDEN:',FLDEN,
C    &   '    I:',I,'    FIELDBND:', FIELDBND
C     Save the value for search limits
      IF( I .EQ. 0) THEN
C        No delimiter on line, must be last field
         FLDEN = STRLEN
         I = STRLEN
      ELSE
         FLDEN = I+FLDST-1
         I = FLDEN - 1
      ENDIF
C
      FNWKC = FIRSTKC(STRING(FLDST:FLDEN), FLDEN-FLDST+1) +FLDST-1
      IF( FNWKC .EQ. 0 .OR. STRING(FNWKC:FNWKC) .EQ. ',' ) THEN
         FNWKC = 0
         LNWKC = 0
      ELSE
         LNWKC = LASTKC(STRING(1:I), I) 
         FIELDBND = .TRUE.
      ENDIF
C
C111  write(*,'(A,A,/,A,A,A1,/,A,I2,A,I2,A,I2,A,I2,A,I2,A,I2,A,L1)')
C    & '  123456789012345678901234567890123456789012345678901234567890',
C    & '12345678901234567890',' {',
C    & STRING(1:STRLEN),'}','    STRLEN:',STRLEN,'    FLDST:',FLDST,
C    & '    FNWKC:',FNWKC,'    LNWKC:',LNWKC,'    FLDEN:',FLDEN,
C    & '    I:',I,'    FIELDBND:', FIELDBND
C
      RETURN 
C
      END
C
C
C
C
C
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION NEXTDLM( STRING, STRLEN )
C
C     Figures out postion of next comma, returns 0 if no comma.
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER STRLEN, I
      CHARACTER STRING*(*)
C
C
      DO 12, I = 1, STRLEN
   12    IF( STRING(I:I) .EQ. ',') GOTO 13

   13 IF(I .GT. STRLEN .OR. STRING(I:I) .NE. ',') THEN
         I = 0
      ENDIF
C
C      write(*,'(A,I2,A,A,A)')' NEXTDLM=',NEXTDLM,'{',STRING(1,NEXTDLM),
C     &   '}'
C
      NextDLM = I
      RETURN 
C
      END
C
C-----------------------------------------------------------------------
C      
      SUBROUTINE MIFLINE(ParameterFile)
C     -------------------      
C
C     Process line from master input file (ie, read in and act on 
C     arameter file names).
C
      IMPLICIT NONE
C                                                                             
C                                                                             
C     Declarations regarding subroutine arguments:                            
C                                                                             
      CHARACTER ParameterFile*256
C               This is line read from the master input file which may
C               contain a parameter filename.
C
C
C     Declarations regarding local variables or functions:
C
C
      LOGICAL   LOGVAR
      INTEGER   I, J, BFN, EFN
C               These are some scratch variables.
C
      INTEGER   FIRSTKC, LASTKC
C               Functions that return positions of the first and last useful
C               chars in a string.
C
      CHARACTER CHNKST*16
C               This is the calling subroutine name passed to ErrMess() in
C               case of error messages.
C
      CHARACTER CSTR*256
C               This is a buffer used to hold lines read from parameter files.
C     
C
C
C     Declarations regarding COMMON block variables:
C
      CHARACTER    CTMP*255, CWRK*80
C                  These are just buffers used throughout program.
C
      CHARACTER    LLT*1
C                  This is the character that marks the end of a line sent
C                  to the LogIt subroutine.
C
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C                  See description in Block Data Genral
C     
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                  See description in Block Data Genral
C  
C
C     INTEGER    NUMPARM
C
C     
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
C     
C     IF( NPFILZ .GT. MPFILZ) THEN
C        CALL ADIOS()
C     ENDIF
C
C
C     WRITE(*,'(1X,A)')ParameterFile
C
      CHNKST='MIFLINE'//LLT
C
C     Find the first useful character on line:
      BFN = FIRSTKC(ParameterFile, 256)
C
      IF(BFN .EQ. 0)THEN
C        Nothing to do here; move on to next line
         RETURN
      ENDIF
C
      IF(ParameterFile(BFN:BFN) .EQ. '#')THEN
C        Comment line; skip it
         RETURN
      ENDIF
C
C     Ok, must process the line.  Supposed to be a filename.
C     Determine end of filename:
      EFN=LastKC( ParameterFile(BFN:), 256 - BFN + 1) + BFN - 1
C
C     Now see if the file exists:
C
C     WRITE(*,'(A,A,A1)')' Looking for file ',ParameterFile(BFN:EFN),'!'
      INQUIRE(FILE=ParameterFile(BFN:EFN),EXIST=LOGVAR)
      IF(.NOT. LOGVAR) THEN
         CALL ErrMess(1004,CHNKST,LOGCDE)
         WRITE(CTMP,2014)ParameterFile(BFN:EFN),LLT
         CALL LogIt(LOGCDE,0,1)
         CALL SAVETOD
         CALL ADIOS
      ENDIF
C     WRITE(*,'(A,A)')' Found file ',ParameterFile(BFN:EFN)
C     
C     Open the parameter file.
C
      OPEN(UNIT=98,FILE=ParameterFile(BFN:EFN),STATUS='OLD',IOSTAT=I) 
      IF ( I .NE. 0 ) THEN                            
         CALL ErrMess(1005,CHNKST,LOGCDE)
         WRITE(CTMP,2015) I,ParameterFile(BFN:EFN),LLT
         CALL LogIt(LOGCDE,1,1)
         CALL BuildInfo
         CALL LogIt(LOGCDE,0,0)
         CALL SAVETOD
         CALL ADIOS
      ENDIF
C
      IF(GenVerbosity .GE. 2 ) THEN
         WRITE(CTMP,'(A,A,A,A)') ' Reading parameters from file ',
     &      ParameterFile(BFN:EFN),' . . .',LLT
         CALL LogIt(LOGCD,1,0)
      ENDIF
C
      DO 1231 J= 1, 32766, 1
C
C        Read in a line at a time, and call processing routine until EOF
C
         READ(98, '(A127)', END=1233, ERR=1232, IOSTAT=I) CSTR
C
C        Write the input line to the detailed log file:
C
         IF( GenVerbosity .GE. 5 ) THEN
            WRITE(CTMP,2014) ParameterFile(BFN:EFN)//
     &         ': '//CSTR(1:LASTKC(CSTR,125)),LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
         CALL PARLINE(CSTR, I)
C
         IF(I .GE. 1000)THEN
C           Function returned an error number               
            WRITE(CTMP,2014) CSTR(1:LASTKC(CSTR,125)), LLT
            CALL LogIt(LOGCDE,0,1)
            WRITE(CTMP,2016) J,ParameterFile(BFN:EFN),LLT
            CALL LogIt(LOGCDE,0,1)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
 1231    CONTINUE
C
      IF (J .GE. 32766)THEN
C        Something is wrong with the file.  Shouldn't have read so may lines.
         CALL ErrMess(1006,CHNKST,LOGCDE)
         WRITE(CTMP,2014) ParameterFile(BFN:EFN),LLT
         CALL LogIt(LOGCDE,1,1)
         CALL SAVETOD
         CALL ADIOS
      ENDIF
C
 1232 CALL ErrMess(1003,CHNKST,LOGCDE)
      WRITE(CTMP,2015) I,ParameterFile(BFN:EFN),LLT
      CALL LogIt(LOGCDE,1,1)
      CALL BuildInfo
      CALL LogIt(LOGCDE,0,0)
      CALL SAVETOD
      CALL ADIOS
C
 1233 CLOSE(98)
C
C     
 2014 FORMAT(1X,A,A1)
 2015 FORMAT(' Got IOSTAT=',I4,' when accessing file ',A,'.',A1)
 2016 FORMAT(' Problem occurred on line ',I4,' of file ',A,'.',A1)
C
      RETURN
      END
C      
C
C
C-----------------------------------------------------------------------
C      
      SUBROUTINE PARLINE(PLINE, EMESSNO)
C     -------------------      
C      
C
C     Process line from a parameter file (ie, read in and act on 
C     parameter id codes).
C
      IMPLICIT NONE
C                                                                             
C                                                                             
C     Declarations regarding subroutine arguments:                            
C                                                                             
      CHARACTER PLINE*127
C               This is line read from the parameter file which supposedly
C               contains parameter data.
C
      INTEGER   EMESSNO
C               This return value tells the calling routine which error
C               message to display.  Returning 0 means no error.
C
C
C     Declarations regarding local variables or functions:
C
      DOUBLE PRECISION TMPDBL
      INTEGER   I, J, K, FKC, LKC, SFNW, EFNW, BPID, PNUM
      CHARACTER PROBENAME*14
      LOGICAL   LOGVAR
C               These are some scratch variables.
C
      LOGICAL   AlreadySpecified(420)
C               Indicates whether a parameter has already been read in
C
      INTEGER   FIRSTKC, LASTKC, ENDSTR
C               Functions that return positions of the first and last useful
C               chars in a string, the end position of the current word, and
C               the position of the next delimiter (usually a comma).
C
      LOGICAL   FIELDBND, READQEQN, ReadValue, UseTSData, 
     &          ReadIntegerList, ReadFloatList, ReadIntegerValue 
C               Functions that help parse comma-delimited input lines,
C               discharge equation lines, elevation values, data file names.
C               
C      DOUBLE PRECISION Round
C
      LOGICAL    ReadAvsZ
C               Function that processes area-elevation parameter records.
C     
      LOGICAL   WhichOutputFiles
C               Function that sets the flags for which output files to make. 
C
      CHARACTER CHNKST*16
C               This is the calling subroutine name passed to ErrMess() in
C               case of error messages.
C
C     Declarations regarding COMMON block variables:
C
      CHARACTER*6400 PARIDS
C               This is a concatenated list of all the parameters that
C               the program knows to look for.  When a parameter ID is read,
C               the program searches for it in this list.  Since all the
C               parameters are eight chars, dividing the return value
C               from INDEX() by eight yields the order number of the 
C               parameter ID in the list.  This order number is used to
C               jump to the label appropriate for processing input
C               lines for that parameter.  
C
      CHARACTER PARCODE*16, INDENT*4
C               These are output-formatting strings used by subroutines
C               specialized for parsing various types of parameters.
C
      CHARACTER CTMP*255, CWRK*80
C               These are just buffers used throughout program.
C
      CHARACTER LLT*1
C               This is the character that marks the end of a line sent
C               to the LogIt subroutine.
C
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C               See description in Block Data Genral
C
      CHARACTER    DirSep*5,OSName(5)*8,MessDB*256
      CHARACTER*64 ModelerName
C               See descriptions in Block Data Genral
C
      CHARACTER*9  MLSolutionName(2)
C               Strings set to 'interative' and 'matrix', used for output.
C
      CHARACTER*16 SMSolutionName(4), USLSolutionName(5)
C               See description in Block Data Genral
C
      CHARACTER*80 RCTITLE(3), MLHEADER
C               These are descriptive lines for use in output.
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
C               See descriptions in SUBROUTINE Plan1977A()
C
      CHARACTER RCOUTDIR*128, RCOUTEXT*16
C               See descriptions in Block Data Genral
C
      CHARACTER*8 RCRUNTYP, IntervalName(4)
C               See description in Block Data Genral
C
      INTEGER   IntervalNameLen(4)
C               See description in Block Data Genral
C
      LOGICAL   LakesToSolve(8)
C               See description in Block Data Genral
C
      INTEGER   StartYear, StartMonth, StartDay, 
     &          EndYear, EndMonth, EndDay
C               See descriptions in Block Data Genral
C
      INTEGER   OSTYPE, NUMPARAM
C               See descriptions in Block Data Genral
C
      DOUBLE PRECISION SGX(17), SGY(17), SGZ(17)
      INTEGER       MaxSGSetting
      INTEGER   BOMRound , FlowRound, LevelRound,
     &          MidLakeSolution, StMarysSolution, UStLawrSolution, 
     &          P77ACompat
C               See descriptions in Block Data Genral
C
      INTEGER   RoutingInterval(3), RoutingIncrement(3)
C               Routing interval to use on Superior, middle lakes,
C               and Ontario (1=daily, 2=weekly, 3=qmonthly, 
C               4=monthly), and the number of increments within each 
C               interval.
C
      INTEGER   MonthLength
      CHARACTER*8  MonthLengthName(2)
C                  See description in Block Data Genral
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
      INTEGER   NOPDIR, NOPEXT, TSOutputFileNum(4,2,8)
      LOGICAL   OnlyOneOpen
C               See descriptions in Block Data Genral
C
      DOUBLE PRECISION SPSTLV, MHSTLV, SCSTLV, ERSTLV, CGSTLV, ONSTLV
C               See descriptions in Block Data Genral
C     
      DOUBLE PRECISION SPBOMLev, MHBOMLev, SCBOMLev, ERBOMLev
C               The beginning-of-month lake levels used by Plan 1977A
C
      INTEGER   PrevStMarysFlow
C               Initial Lake Superior outflow (average for month 0)
C 
      REAL      PrevSupGates 
C               Initial gates open at Comp Works
C 
      REAL      SupGatesOpen
C               The number of Superior gates open in current routing 
C               period.
C
      DOUBLE PRECISION SPK,SPYM,SPA,SPB,SPWT,SPC
      DOUBLE PRECISION MHK,MHYM,MHA,MHB,MHWT,MHC
      DOUBLE PRECISION SCK,SCYM,SCA,SCB,SCWT,SCC
      DOUBLE PRECISION ERK,ERYM,ERA,ERB,ERWT,ERC
      DOUBLE PRECISION ONK,ONYM,ONA,ONB,ONWT,ONC
C               See descriptions in CHECKPAR and OUTFLOW
C
      INTEGER       SideChanMaxCapacity(12), SideChanMinCapacity(12)
      DOUBLE PRECISION CWOvertopElev, CWMinimumFreeboard(12)
C               See descriptions in RouteSuperior
C
      LOGICAL   DoDSMontrealComp
      INTEGER   SupIndReadArray(20)
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
      CHARACTER MonthAbbrev(12)*3
      INTEGER   SPOPNBS06(12,6), SPOPNBS07(12,6), OPLLO(12), 
     &          MHOPNBS06(12,6), MHOPNBS07(12,6), OPChicago, 
     &          SCRRet(12), SCOPNBS06(12,6), SCOPNBS07(12,6), 
     &          DetRet(12), EROPNBS06(12,6), EROPNBS07(12,6),
     &          OPWelland(12), NiaRet(12),
     &          SPOPNBS08(12,6), MHOPNBS08(12,6),  
     &          SCRRet_P2012(12), SCOPNBS08(12,6), 
     &          DetRet_P2012(12), EROPNBS08(12,6),
     &          OPWelland_P2012(12), NiaRet_P2012(12)
      INTEGER   ForecastMonths, LStoSWP, SupCritRnd, MonthsChecked
C                   See descriptions in Block Data Genral
C
      INTEGER   CritBOMRound, CritFlowRound, CritLevelRound, 
     &          Q1887Method, CritBMethod
      INTEGER      USSlipEquation
C                  1=ancient eq (default), 2=2010 EC eq
      DOUBLE PRECISION USSIceWeedAdj
C

      DOUBLE PRECISION CritCThresholdEl
      LOGICAL   BackCompat
      LOGICAL   WriteTSOutput(4,2,8), WriteErieOutflow
C               See descriptions in Block Data General.  
C
      CHARACTER*256 CritKfile
      DOUBLE PRECISION P2012SPTARGET(12), P2012MHTARGET(12),
     &                 P2012SPHIGHLEVEL, P2012SPLOWLEVEL,
     &                 P2012MHTARGETADJ   
      DOUBLE PRECISION P2012PARA(6) 
      INTEGER          COUNTER_JUNEFLOW, Last_JuneFlow
      INTEGER          IndexGatesSide
      DOUBLE PRECISION ForcedSupGates(12)
      INTEGER          ForcedSideChan(12)       
C       
C  Common block statements.
C 
      COMMON /P2012/    P2012PARA, P2012SPTARGET, P2012MHTARGET, 
     &                  P2012SPHIGHLEVEL, P2012SPLOWLEVEL,
     &                  P2012MHTARGETADJ      
      COMMON /STURGEON/ COUNTER_JUNEFLOW, Last_JuneFlow 
      COMMON /FORGATE/ ForcedSupGates  
      COMMON /FORSIDE/ ForcedSideChan
      COMMON /INFORCE/ IndexGatesSide
         
      COMMON /CRITK/    CritKfile
C
      COMMON /OUTSTF/   WriteTSOutput, WriteErieOutflow
      COMMON /PARMRD/   AlreadySpecified
      COMMON /MONAME/   MonthAbbrev
      COMMON /INI58D/   CritDepths, P58DTotalSup, MaxILimit, MinPLimit,
     &                  P58DWeightedSupply, P58DRuleCurve, P58DAppLimit,
     &                  P58DRunType, P58DForecastYear, P58DForecastMon,
     &                  P58DForecastDay, SupplyIndicator, 
     &                  P58DPlanLevel, P58DActualLevel,
     &                  P58DSimulateFile, P58DSIDataFile, ElimDev, 
     &                  ForceFlow, SetFlow, PreProjConditions,
     &                  P58DUseDatabase, P58DSiDatabase
      common /Ini58dc/  P58DAppLimitSym
      COMMON /INI58D1/  P58DAdjustment,P58DrcAdjustment,P58DActualFlow, 
     &                  P58DTotalDev,P58DLSLflow
      COMMON /INI58D2/  ComputeMontreal
      common /PreCond/  PreprojectLevel
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /OUTNAM/   RCOUTDIR, RCOUTEXT
      COMMON /OUTNMI/   NOPDIR, NOPEXT, TSOutputFileNum, OnlyOneOpen
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /RCRUNC/   RCRUNTYP, IntervalName, MonthLengthName,
     &                  MLSolutionName, SMSolutionName,
     &                  USLSolutionName
      COMMON /RCRUNI/   LakesToSolve, 
     &                  IntervalNameLen,  MonthLength,
     &                  StartYear, StartMonth, StartDay, 
     &                  EndYear, EndMonth, EndDay
      COMMON /CPARID/   PARIDS
      COMMON /TITLES/   RCTITLE, MLHEADER
      COMMON /GENRLC/   OSName,DirSep,MessDB, ModelerName
      COMMON /GENRLI/   OSTYPE, NUMPARAM
      COMMON /PARSEC/   INDENT, PARCODE
      COMMON /STAGE/    MHK,SCK,ERK,MHYM,SCYM,ERYM,MHA,SCA,ERA,
     &                  MHB,SCB,ERB,MHWT,SCWT,ERWT,MHC,SCC,ERC
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
      COMMON /SPDAT/    SPK,SPYM,SPA,SPB,SPWT,SPC,
     &                  SGX, SGY, SGZ, MaxSGSetting,
     &                  SupGatesOpen
      COMMON /SPOTOP/   SideChanMaxCapacity, CWOvertopElev, 
     &                  CWMinimumFreeboard, SideChanMinCapacity
      COMMON /ONDAT/    ONK,ONYM,ONA,ONB,ONWT,ONC
      COMMON /MISINT/   MidLakeSolution, BOMRound , FlowRound, 
     &                  LevelRound, RoutingInterval, RoutingIncrement,
     &                  StMarysSolution, UStLawrSolution, P77ACompat
      COMMON /INITLV/   SPSTLV, MHSTLV, SCSTLV, ERSTLV, CGSTLV, ONSTLV,
     &                  PrevStMarysFlow, PrevSupGates, 
     &                  SPBOMLev, MHBOMLev, SCBOMLev, ERBOMLev 
      COMMON /SUPCOM/   SPOPNBS06, OPLLO, MHOPNBS06, OPChicago, SCRRet,
     &                  SCOPNBS06, DetRet, EROPNBS06, OPWelland, NiaRet,
     &                  SPOPNBS07, MHOPNBS07, SCOPNBS07, EROPNBS07,
     &                  SPOPNBS08, MHOPNBS08, SCOPNBS08, EROPNBS08,
     &                  OPWelland_P2012, SCRRet_P2012, DetRet_P2012,
     &                  NiaRet_P2012,
     &                  ForecastMonths,
     &                  LStoSWP, SupCritRnd, BackCompat, MonthsChecked
      COMMON /SUPCRT/   CritBOMRound, CritFlowRound, CritLevelRound,
     &                  Q1887Method, CritBMethod, CritCThresholdEl
      COMMON /SLDS/     DoDSMontrealComp
      COMMON /USSLIP/   USSIceWeedAdj, USSlipEquation

C
      EQUIVALENCE (SupplyIndicator(2), SupIndReadArray)
C
      CHNKST='PARLINE'//LLT
C
      EMESSNO = 0
C
C     Locate first useful char
      BPID = FIRSTKC(PLINE, 126)
C
      IF(BPID .LE. 0 .OR. PLINE(BPID:BPID) .EQ. '#')THEN
C        It's a comment or empty line, return without doing anything
         RETURN
      ENDIF
C
C     Ok, must process the line.  Supposed to be parameter id,
C     followed by a semi-colon.  Find first ':' on line:
C
      I = INDEX(PLINE(BPID:), ':')
      IF (I .LE. 6)THEN
         IF( I .LT. 1)THEN
            EMESSNO = 1007
         ELSE IF( I .EQ. 1)THEN
            EMESSNO = 1066
         ELSE
            EMESSNO = 1008
         ENDIF
      CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
      RETURN
      ENDIF
C
C
C     First find the end of the string:
      LKC = LASTKC(PLINE, 127)
C
C     Get start of data pieces after ':' 
      J = BPID + I
      FKC = FIRSTKC(PLINE(J:127), 127 - J) 
      IF(FKC .LE. 0 )THEN
C        No data!
         EMESSNO = 1009
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF
      FKC = FKC + J - 1
C
C
C     Ok, established that have an ID code, with data following a ':'.
C     Copy it to PARCODE, making sure everything uppercase:
C
C      
C     Using SFNW as scratch here, gets re-defined with a more specific
C     meaning later on.  
      SFNW = 0
      DO 2566, J = 1, I - 1
         LOGVAR= .FALSE.
         I = ICHAR( PLINE(BPID+J-1:BPID+J-1) )
         IF( I .GE. 97 .AND. I .LE. 122) THEN
            I = I-32
            LOGVAR = .TRUE.
         ELSE IF( (I .GE. 48 .AND. I .LE. 57) .OR. 
     &             I .GE. 65 .AND. I .LE. 90) THEN
C            PARCODE(SFNW:SFNW) = PLINE(BPID+J-1:BPID+J-1)
            LOGVAR = .TRUE.
         ENDIF
         IF(LOGVAR)THEN
            SFNW = SFNW + 1 
            IF (SFNW .GT. 16) THEN
              EMESSNO = 1066
              CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
              RETURN
            ENDIF
            PARCODE(SFNW:SFNW) = CHAR(I)
         ENDIF
 2566 CONTINUE
C
C     Pad any remaining spots in PARCODE with spaces      
      DO 2567, J = SFNW+1, 16
 2567    PARCODE(J:J) = ' '
C 
C     Time to get messy and tedious.  Must go through and figure out which
C     parameter it is and then how to read the data list.
C
      J = INDEX(PARIDS, PARCODE)
      IF(J .LE. 0 )THEN
C        Not in list!
         EMESSNO = 1010
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF
C
C     Adjust to last character of the ID code in PARIDS
      J = J + 15
C
C     Make sure not a strange coincidence match.
C     The parids are padded to fillup all the 16 characters, so if ends on
C     complete id code, will be evenly divisble.
      I=MOD(J,16)
      IF(I .NE. 0 )THEN
         EMESSNO = 1011
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF
C
C     Get the index number of the array
      PNUM=J/16
C
C     Now see if we already read it
      IF( AlreadySpecified(PNUM) ) THEN
         EMESSNO = 1173
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
C        WRITE(CTMP,'(A,A,A,A1)') 'Warning!  This line overrides'//
C    &      ' previously specified value for ', PARCODE, ':', LLT
C        CALL LogIt(LOGCD,0,0)
C        WRITE(CTMP,'(A,A1)') PLINE(1:LKC), LLT
C        CALL LogIt(LOGCD,0,0)
      ELSE
         AlreadySpecified(PNUM) = .TRUE.
      ENDIF
C
C
C     WRITE(*,'(A,A,A,A,A,I3,A,A,A)')'     ',PLINE(BPID:BPID+15),
C    &   ' (',PARCODE,', index #',PNUM,') ===>',PLINE(FKC:LKC),'<==='
C
      IF ( PNUM .LE. 100 ) THEN
         GOTO (4001,4002,4003,4004,4005,4005,4005,4008,4009,4010,
     1         4011,4012,4013,4014,4015,4016,4017,4018,4019,4020,
     2         4021,4022,4023,4024,4025,4026,4027,4028,4029,4030,
     3         4031,4032,4033,4034,4035,4036,4037,4038,4039,4040,
     4         4041,4042,4043,4044,4045,4046,4047,4048,4049,4050,
     5         4051,4052,4053,4054,4055,4056,4057,4058,4059,4060,
     6         4061,4062,4063,4064,4065,4066,4067,4068,4069,4070,
     7         4071,4072,4073,4074,4075,4076,4077,4078,4079,4080,
     8         4081,4082,4083,4084,4085,4086,4087,4088,4089,4090,
     9         4091,4092,4093,4094,4095,4096,4097,4098,4099,4100) PNUM
      ELSEIF ( PNUM .LE. 200 ) THEN
         I = PNUM - 100
         GOTO (4101,4102,4103,4104,4105,4106,4107,4108,4109,4110,
     1         4111,4112,4113,4114,4115,4116,4117,4118,4119,4120,
     2         4121,4122,4123,4124,4125,4126,4127,4128,4129,4130,
     3         4131,4132,4133,4134,4135,4136,4137,4138,4139,4140,
     4         4141,4142,4143,4144,4145,4146,4147,4148,4149,4150,
     5         4151,4152,4153,4154,4155,4156,4157,4158,4159,4160,
     6         4161,4162,4163,4164,4165,4166,4167,4168,4169,4170,
     7         4171,4172,4173,4174,4175,4176,4177,4178,4179,4180,
     8         4181,4182,4183,4184,4185,4186,4187,4188,4189,4190,
     9         4191,4192,4193,4194,4195,4196,4197,4198,4199,4200) I
      ELSEIF ( PNUM .LE. 300 ) THEN
         I = PNUM - 200
         GOTO (4201,4202,4203,4204,4205,4206,4207,4208,4209,4210,
     1         4211,4212,4213,4214,4215,4216,4217,4218,4219,4220,
     2         4221,4222,4223,4224,4225,4226,4227,4228,4229,4230,
     3         4231,4232,4233,4234,4235,4236,4237,4238,4239,4240,
     4         4241,4242,4243,4244,4245,4246,4247,4248,4249,4250,
     5         4251,4252,4253,4254,4255,4256,4257,4258,4259,4260,
     6         4261,4262,4263,4264,4265,4266,4267,4268,4269,4270,
     7         4271,4272,4273,4274,4275,4276,4277,4278,4279,4280,
     8         4281,4282,4283,4284,4285,4286,4287,4288,4289,4290,
     9         4291,4292,4293,4294,4295,4296,4297,4298,4299,4300) I
      ELSE
         I = PNUM - 300
         GOTO (4301,4302,4303,4304,4305,4306,4307,4308,4309,4310,
     1         4311,4312,4313,4314,4315,4316,4317,4318,4319,4320,
     2         4321,4322,4323,4324,4325,4326,4327,4328,4329,4330,
     3         4331,4332,4333,4334,4335,4336,4337,4338,4339,4340,
     4         4341,4342,4343,4344,4345,4346,4347,4348,4349,4350,
     5         4351,4352,4353,4354,4355,4356,4357,4358,4359,4360,
     6         4361,4362,4363,4364,4365,4366,4367,4368,4369,4370,
     7         4371,4372,4373,4374,4375,4376,4377,4378,4379,4380,
     8         4381,4382,4383,4384,4385,4386,4387,4388,4389,4390,
     9         4391,4392,4393,4394,4395,4396,4397,4398,4399,4400) I
      ENDIF
C
      EMESSNO = 1012
      CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
      WRITE(CTMP,'(A,A,A,I3,A,A1)')' Parameter ID "',PARCODE,
     &   '" matches index number ',PNUM,'.  Read from line:',LLT
      CALL LogIt(LOGCDE,0,0)
      RETURN
C
C
C     OK, here's the truly ugly part.  
C   
C     -------------------      
C
C     Parse MODELERNAME line  (parameter #1)
C     
C     Assign contents of line to ModelerName
 4001 AlreadySpecified(PNUM) = .FALSE.
      ModelerName = PLINE(FKC:LKC)
C     write(*, '(1X, I2,A,A,A,I2)') 
C
C     Pad the rest of ModelerName with LLT
      DO 3001, I = LKC - FKC + 2, 64
 3001    ModelerName(I:I) = LLT
C
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5017) INDENT, PARCODE, ': Modeler Name set to {',
     &      ModelerName(1:LASTKC(ModelerName,64)), '}', LLT
      ENDIF
      GOTO 5000
C
C
C     -------------------      
C
C
C     Parse MESSAGEDATABASE line  (parameter #2)
C     
C     Determine end of filename
 4002 I=LastKC( PLINE(FKC:LKC), LKC - FKC + 1) - 1
C
C     WRITE(*,'(A,A,A1)')' Looking for file {',PLINE(FKC:FKC+I),'}'
C     Now see if the file exists:
      INQUIRE(FILE=PLINE(FKC:FKC+I),EXIST=LOGVAR,ERR=0803)
 0803 IF(.NOT. LOGVAR) THEN
         EMESSNO = 1013
         WRITE(CTMP,5018)' Could not find the message data',
     &      'base file specified by parameter ', LLT
         CALL LogIt(LOGCDE,0,0)
         WRITE(CTMP,5019)' MESSAGEDATABASE.  The filename was ',
     &      'given as {',PLINE(FKC:FKC+I), '}', LLT 
         CALL LogIt(LOGCDE,0,0)
         WRITE(CTMP,5018)' Input ', 'line:',LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C     
C     Ok, file is there (we'll assume that it's readable and intact)
C     Erase current contents of MessDB and put in new name:
      DO 3002, J=1,256
 3002    MessDB(J:J) = LLT
      MessDB(1:I+1)=PLINE(FKC:FKC+I)
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5017) INDENT,PARCODE,
     &      ': set message database filename to {',
     &      MessDB(1:I+1),'}',LLT
      ENDIF
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse OUTPUTDIRECTORY line  (parameter #3)
C     
C     Determine end of directory name
 4003 NOPDIR=LastKC( PLINE(FKC:LKC), LKC - FKC + 1)
C
C     Erase current contents of RCOUTDIR and put in new name:
      DO 3003, J=1,128
 3003    RCOUTDIR(J:J) = LLT
      RCOUTDIR=PLINE(FKC:FKC+NOPDIR-1)
C
C     Now see if directory exists and writable by opening/writing to a file.
      PROBENAME='rootname.12345'
C
C     If extension already defined, add it to existing 8-char root name:
      IF(NOPEXT .GT. 0)THEN
         PROBENAME(9:8+NOPEXT)=RCOUTEXT
         I = 8 + NOPEXT
      ELSE
         I = 12
      ENDIF
C
      OPEN(UNIT=97,FILE=RCOUTDIR(1:NOPDIR)//PROBENAME(1:I),IOSTAT=J) 
      IF ( J .NE. 0 ) THEN
         EMESSNO = 1015
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5015) J,RCOUTDIR(1:NOPDIR)//PROBENAME,LLT
         CALL LogIt(LOGCDE,1,0)
         CALL BuildInfo
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
      WRITE(97,'(A)',ERR=2003, IOSTAT=J)' Hi!'
      CLOSE(97,ERR=1003,IOSTAT=J,STATUS='DELETE')
C     
C     Ok, established that we can open/write files in RCOUTDIR.
C     now figure out OS and DirSep
C   
      OSTYPE=INDEX(DirSep,RCOUTDIR(NOPDIR:NOPDIR))
      IF ( OSTYPE .LE. 0) THEN
C        uh-oh, not in list, but not fatal error.  squawk and continue.
         CALL ErrMess(1019,CHNKST,LOGCDE)
         WRITE(CTMP,5019) ' Last character in ',RCOUTDIR(1:NOPDIR),
     &       ' not found in known list: ',DirSep,LLT
         CALL LogIt(LOGCDE,1,0)
         OSTYPE=1
C     ELSE
      ENDIF
      WRITE(CTMP,5018)
     &   '===>Probable operating system: ',OSName(OSTYPE),LLT
      CALL LogIt(LOGD,0,0)
C
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5017) INDENT,PARCODE,': set output directory to {',
     &      RCOUTDIR(1:NOPDIR),'}',LLT
      ENDIF
C
      GOTO 5000 
C
 2003 EMESSNO = 1016
      GOTO 0903
 1003 EMESSNO = 1017
 0903 CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
      WRITE(CTMP,5015) J,RCOUTDIR(1:NOPDIR)//PROBENAME,LLT
      CALL LogIt(LOGCDE,1,0)
      CALL BuildInfo
      CALL LogIt(LOGCDE,0,0)
      RETURN
C
C     -------------------      
C
C
C     Parse OUTPUTEXTENSION line  (parameter #4)
C     
 4004 NOPEXT=LastKC( PLINE(FKC:LKC), LKC - FKC + 1) 
      IF( NOPEXT .GT. 9) THEN
         EMESSNO = 1018
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5021) ' Attempted length of',NOPEXT, 
     &      ' characters in line:',LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
      RCOUTEXT=PLINE(FKC:FKC+NOPEXT-1)//LLT
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5017)INDENT,PARCODE,
     &      ': set output file extension to {',
     &      RCOUTEXT(1:NOPEXT),'}',LLT
      ENDIF
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse TITLE1, TITLE2, TITLE3 lines  (parameters #5,6,7)
C     
C     Convert parameter number to an index from 1-3
 4005 I=PNUM - 4
C
      IF( LKC .GT. 79 + FKC - 1 ) THEN
         WRITE(CTMP,5005) I, LLT
         CALL LogIt(LOGCD,1,0)
         LKC = 79 + FKC - 1
      ENDIF
 5005 FORMAT( ' Only first 79 characters in Title ', I1, 
     &   ' will be used.', A1 )
C      
C     Assign contents of line to RCTITLE(I)
      RCTITLE(I)=PLINE(FKC:LKC)
C
C     Pad extra length in RCTITLE(I) with line terminating characters
      DO 3005 J=80,LKC+1,-1
 3005    RCTITLE(I)(J:J) = LLT
C
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5016) INDENT, PARCODE, ': Job title #',I,
     &      ' set to {', RCTITLE(I)(1:LKC-FKC+1),'}', LLT
      ENDIF
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ??????? line  (parameter #8)
C     
C     Not Implemented!
 4008 GOTO 4990
C
C     -------------------      
C
C     Parse ONLYONEOPEN line  (parameter #9)
C     
 4009 IF( PLINE(FKC:FKC) .EQ. 'F' .OR. PLINE(FKC:FKC) .EQ. 'f') THEN 
C        If explicitly given as FALSE, set it FALSE.  This is 
C        default anyway.
         OnlyOneOpen = .FALSE.
      ELSE
C        If find this parameter and doesn't start with 'F', assume 
C        they want it set to TRUE.
         OnlyOneOpen = .TRUE.
      ENDIF
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5023) INDENT, PARCODE, 
     &      ': Only one output file open at a time set to {',
     &      OnlyOneOpen, '}', LLT
      ENDIF
      GOTO 5000
C
C
C     -------------------      
C     Parse GENERALVERBOSITY line  (parameter #10)
C     
C     Unset the 'overriding' flag for this parameter, so the
C     doesn't squawk about it.
 4010 AlreadySpecified(PNUM) = .FALSE.
C
C     Parse the line:
      IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': set general message level to '//LLT,
     &   ' '//LLT, 1, 9, 1014, GenVerbosity)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse TSDATAVERBOSITY line  (parameter #11)
C     
 4011 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': set time-series data message level to '//LLT,
     &   ' '//LLT, 1, 9, 1014, TSDVerbosity)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse CHECKVERBOSITY line  (parameter #12)
C     
 4012 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': set consistency-check message level to '//LLT,
     &   ' '//LLT, 1, 9, 1014, ChkVerbosity)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse SUPVERBOSITY line  (parameter #13)
C     
 4013 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': set Superior message level to '//LLT,
     &   ' '//LLT, 1, 9, 1014, SupVerbosity)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse MIDLAKEVERBOSITY line  (parameter #14)
C     
 4014 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': set middle lakes message level to '//LLT,
     &   ' '//LLT, 1, 9, 1014, MLVerbosity)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse ONTVERBOSITY line  (parameter #15)
C     
 4015 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': set Ontario/St. Lawrence message level to '//LLT,
     &   ' '//LLT, 1, 9, 1014, OntVerbosity)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse STARTDATE line  (parameter #16)
C     
C     Read starting year into variable StartYear
 4016 IF( .NOT. FIELDBND(PLINE, LKC, FKC, SFNW, EFNW, I ) ) THEN 
C        Since there should be a useful character on line to make it this 
C        far, the only way to achieve this condition should be to specify 
C        a ',' as the next charcater after ':'.  dump if it happens anyway.  
         EMESSNO = 1027
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF 
C     WRITE(*,'(A,A,A)')' trimmed year field ::', PLINE(SFNW:EFNW),'::'
C     Save the begnning of the next field, if any, in J.
      J = I + 1
C
      READ( PLINE(SFNW:EFNW), 5012, ERR=3016, IOSTAT=I ) StartYear
      IF( StartYear .LE. 0 .OR. StartYear .GE. 99999) THEN
         EMESSNO = 1021
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5024) StartYear, LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
      IF( .NOT. FIELDBND(PLINE, LKC, J, SFNW, EFNW, I ) ) THEN 
C        No delimiter, or no more data.  Presume only year specified on line.
         StartMonth = 1         
         StartDay = 1
         IF( GenVerbosity .GE. 4 ) THEN
C        Spit out a warning of our assumption
            IF( GenVerbosity .GE. 6 ) THEN
               CALL ErrMess(1032,CHNKST,LOGC)
            ENDIF
            CALL ErrMess(1022,CHNKST,LOGD)
         ENDIF
         GOTO 5516
      ENDIF 
C     WRITE(*,'(A,A,A)')' trimmed month field ::', PLINE(SFNW:EFNW),'::'
C     Save the begnning of the next field, if any, in J.
      J = I + 1
C
      READ( PLINE(SFNW:EFNW), 5012, ERR=2016, IOSTAT=I ) StartMonth
      IF( StartMonth .LE. 0 .OR. StartMonth .GT. 12) THEN
         EMESSNO = 1024
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5024) StartMonth, LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
      IF( .NOT. FIELDBND(PLINE, LKC, J, SFNW, EFNW, I ) ) THEN 
C        No delimiter, or no more data.  Presume no day specified on line.
         StartDay = 1
         IF( GenVerbosity .GE. 4 ) THEN
C        Spit out a warning of our assumption
            IF( GenVerbosity .GE. 6 ) THEN
               CALL ErrMess(1033,CHNKST,LOGC)
            ENDIF
            CALL ErrMess(1031,CHNKST,LOGD)
         ENDIF
         GOTO 5516
      ENDIF
C     
C     Find # of days in month, store in J.  Use I for scratch.
      CALL NUMDAYS (StartYear, StartMonth, I, J)
C
      READ( PLINE(SFNW:EFNW), 5012, ERR=1016, IOSTAT=I ) StartDay
      IF( StartDay .LE. 0 .OR. StartDay .GT. J) THEN
         EMESSNO = 1026
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5025) StartDay, StartMonth, StartYear, J, LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
 5516 IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5030) INDENT, PARCODE, 'Start', StartDay, 
     &      MonthAbbrev(StartMonth), StartYear, LLT
      ELSE
         CTMP(1:1)=LLT
      ENDIF
      GOTO 5000 

 3016 EMESSNO = 1020
      GOTO 0016
 2016 EMESSNO = 1023
      GOTO 0016
 1016 EMESSNO = 1037
 0016 CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
      WRITE(CTMP,5022)I,LLT
      CALL LogIt(LOGCDE,0,0)
      CALL BuildInfo
      CALL LogIt(LOGCDE,0,1)
      RETURN
C
C     -------------------      
C
C
C     Parse ENDDATE line  (parameter #17)
C     
C     Read ending year into variable EndYear
 4017 IF( .NOT. FIELDBND(PLINE, LKC, FKC, SFNW, EFNW, I ) ) THEN 
C        Since there should be a useful character on line to make it this 
C        far, the only way to achieve this condition should be to specify 
C        a ',' as the next charcater after ':'.  dump if it happens anyway.  
         EMESSNO = 1027
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF 
C     WRITE(*,'(A,A,A)')' trimmed year field ::', PLINE(SFNW:EFNW),'::'
C     Save the begnning of the next field, if any, in J.
      J = I + 1
C
      READ( PLINE(SFNW:EFNW), 5012, ERR=3017, IOSTAT=I ) EndYear
      IF( EndYear .LE. 0 .OR. EndYear .GE. 99999) THEN
         EMESSNO = 1021
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5024) EndYear, LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
      IF( .NOT. FIELDBND(PLINE, LKC, J, SFNW, EFNW, I ) ) THEN 
C        No delimiter, or no more data.  Presume only year specified on line.
         EndMonth = 12         
         EndDay = 31
         IF( GenVerbosity .GE. 4 ) THEN
C        Spit out a warning of our assumption
            IF( GenVerbosity .GE. 6 ) THEN
               CALL ErrMess(1034,CHNKST,LOGC)
            ENDIF
            CALL ErrMess(1030,CHNKST,LOGD)
         ENDIF
         GOTO 5517
      ENDIF 
C     WRITE(*,'(A,A,A)')' trimmed month field ::', PLINE(SFNW:EFNW),'::'
C     Save the begnning of the next field, if any, in J.
      J = I + 1
C
      READ( PLINE(SFNW:EFNW), 5012, ERR=2017, IOSTAT=I ) EndMonth
      IF( EndMonth .LE. 0 .OR. EndMonth .GT. 12) THEN
         EMESSNO = 1024
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5024) EndMonth, LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
C     Find # of days in month, store in FKC.  Overwrite values in FKC and I
C     because I'm too lazy to declare some more integer variables, and I'm
C     done with them anyway.
      CALL NUMDAYS (EndYear, EndMonth, I, FKC)
C
      IF( .NOT. FIELDBND(PLINE, LKC, J, SFNW, EFNW, I ) ) THEN 
C        No delimiter, or no more data.  Presume no day specified on line.
         EndDay = FKC
         IF( GenVerbosity .GE. 4 ) THEN
C        Spit out a warning of our assumption
            IF( GenVerbosity .GE. 6 ) THEN
               CALL ErrMess(1035,CHNKST,LOGC)
            ENDIF
            CALL ErrMess(1031,CHNKST,LOGD)
         ENDIF
         GOTO 5517
      ENDIF
C     
      READ( PLINE(SFNW:EFNW), 5012, ERR=1017, IOSTAT=I ) EndDay
      IF( EndDay .LE. 0 .OR. EndDay .GT. FKC) THEN
         EMESSNO = 1026
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5025) EndDay, EndMonth, EndYear, FKC, LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
 5517 IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5030) INDENT, PARCODE, 'End', EndDay, 
     &      MonthAbbrev(EndMonth), EndYear, LLT
      ELSE
         CTMP(1:1)=LLT
      ENDIF
      GOTO 5000 
C
 3017 EMESSNO = 1029
      GOTO 0017
 2017 EMESSNO = 1036
      GOTO 0017
 1017 EMESSNO = 1038
 0017 CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
      WRITE(CTMP,5022)I,LLT
      CALL LogIt(LOGCDE,0,0)
      CALL BuildInfo
      CALL LogIt(LOGCDE,0,1)
      RETURN
C
C     
C     -------------------      
C
C
C     Parse MonthLength line  (parameter #18)
C     
C     Only look at first character
 4018 IF( PLINE(FKC:FKC) .EQ. 'u' .OR. PLINE(FKC:FKC).EQ. 'U') THEN
         MonthLength = 1
      ELSE IF( PLINE(FKC:FKC) .EQ. 'a' .OR. PLINE(FKC:FKC).EQ. 'A') THEN
         MonthLength = 2
      ELSE
         EMESSNO = 1098
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5027) PLINE(FKC:FKC), LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5028) INDENT, PARCODE, 
     &      ': Month or quarter-month length set to {',
     &      MonthLengthName(MonthLength)
     &      (1:EndStr(MonthLengthName(MonthLength),8)), '}',LLT
      ENDIF
      GOTO 5000 
C
C     -------------------      
C
C     Parse LAKESTOSOLVEFOR line  (parameter #19)
C     The valid combinations of s (Superior), m (Michigan-Huron),
C     c (St. Clair), e (Erie), o (Ontario)
C
C     Determine end of string
 4019 I=ENDSTR( PLINE(FKC:LKC), LKC - FKC + 1)
      IF( I .GT. 8) THEN
         STOP 'bad string'
      ENDIF
C
C     Copy to CWRK, and make sure it's uppercase
      CWRK = PLINE(FKC: FKC + I - 1)
      DO 3019, J = 1, I
 3019    IF(CWRK(J:J) .GE. 'a' .AND. CWRK(J:J) .LE. 'z' )
     &      CWRK(J:J)=CHAR( ICHAR(PLINE(FKC+J-1:FKC+J-1)) - 32 )
C
C     See if there is an 's'; if so, need to solve for Superior level
      IF( INDEX( CWRK(1:I), 'S') .GT. 0 )LakesToSolve(1) = .TRUE.
C
C     See if there is an 'm'; if so, need to solve for MicHuron level
      IF( INDEX( CWRK(1:I), 'M') .GT. 0 )LakesToSolve(2) = .TRUE.
C
C     See if there is an 'c'; if so, need to solve for St. Clair level
      IF( INDEX( CWRK(1:I), 'C') .GT. 0 )LakesToSolve(3) = .TRUE.
C
C     See if there is an 'e'; if so, need to solve for Erie level
      IF( INDEX( CWRK(1:I), 'E') .GT. 0 )LakesToSolve(4) = .TRUE.
C
C     See if there is an 'o'; if so, need to solve for Ontario level
      IF( INDEX( CWRK(1:I), 'O') .GT. 0 )LakesToSolve(5) = .TRUE.
C
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5017) INDENT,PARCODE,
     &      ': will solve for lakes indicated by {',
     &       CWRK(1:I),'}',LLT
C        CALL LogIt(LOGD,1,0)
      ELSE
C        Avoid repeating contents of CTMP
         CTMP(1:1)=LLT
      ENDIF
C
      GOTO 5000
C
C     -------------------      
C
C
C     Parse SUPSTARTLEVEL line  (parameter #20)
C     
C     Parse the line:
 4020 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,SPSTLV))GOTO 5000
      EMESSNO = 9999
      RETURN
C
C
C     -------------------      
C
C
C     Parse SUPSIDECHANDATA line  (parameter #21)
C     
C     Parse the line:
 4021 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 69)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse SUPNBSDATA line  (parameter #22)
C     
C     Parse the line:
 4022 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 42)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse SUPCONUSEDATA line  (parameter #23)
C     
C     Parse the line:
 4023 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 43)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse LONGLACDIVDATA line  (parameter #24)
C     
C     Parse the line:
 4024 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 44)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse OGOKIDIVDATA line  (parameter #25)
C     
C     Parse the line:
 4025 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 45)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse COMBLLODIVDATA line  (parameter #26)
C     
C     Parse the line:
 4026 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 46)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse SUPOTHERDIVDATA line  (parameter #27)
C     
C     Parse the line:
 4027 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 47)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse SUPOTHERSUPDATA line  (parameter #28)
C     
C     Parse the line:
 4028 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 48)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse SUPRUNOFF line  (parameter #29)
C     
C     Not Implemented!
 4029 GOTO 4990 
C
C     -------------------      
C
C
C     Parse SUPPRECIPDATA line  (parameter #30)
C     
C     Not Implemented!
 4030 GOTO 4990 
C
C     -------------------      
C
C
C     Parse SUPEVAPDATA line  (parameter #31)
C     
C     Not Implemented!
 4031 GOTO 4990 
C
C     -------------------      
C
C
C     Parse SUPDEVIATIONDATA line  (parameter #32)
C     
C     Parse the line:
 4032 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 49)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse SUPEOPLEVELDATA line  (parameter #33)
C     
C     Parse the line:
 4033 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 70)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse SUPGATESDATA line  (parameter #34)
C     
C4034 EMESSNO = 1174
C     WRITE(CTMP,'(A)') ' ==> Forced gate settings not implemented '//
C    &      'yet for Lake Superior. <=='//LLT
C     CALL LogIt(LOGCDE,1,0)
C     CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
C     RETURN
C     Parse the line:
 4034 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 40)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse SUPOUTFLOWDATA line  (parameter #35)
C     
C     Parse the line:
 4035 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 68)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse MHUSTARTLEVEL line  (parameter #36)
C
C     Parse the line:
 4036 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,MHSTLV))GOTO 5000
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse ???????? line  (parameter #37)
C     
C     Determine end of filename
 4037 I=ENDSTR( PLINE(FKC:LKC), LKC - FKC + 1) - 1
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse MHUNBSDATA line  (parameter #38)
C     
C     Parse the line:
 4038 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 1)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse MHUCONUSEDATA line  (parameter #39)
C
C     Determine end of filename
 4039 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 23)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse CHICAGODIVDATA line  (parameter #40)
C     
C     Parse the line:
 4040 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 22)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse MHUOTHERDIVDATA line  (parameter #41)
C     
C     Parse the line:
 4041 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 24)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse MHUOTHERSUPDATA line  (parameter #42)
C     
C     Parse the line:
 4042 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 25)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse MHURUNOFFDATA line  (parameter #43)
C     Not Implemented!
 4043 GOTO 4990 
C
C     -------------------      
C
C
C     Parse MHUPRECIPDATA line  (parameter #44)
C     Not Implemented!
 4044 GOTO 4990 
C
C     -------------------      
C
C
C     Parse MHUEVAPDATA line  (parameter #45)
C     Not Implemented!
 4045 GOTO 4990 
C
C     -------------------      
C
C
C     Parse MHUEOPLEVELDATA line  (parameter #46)
C 
C     Parse the line:
 4046 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 71)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ???????? line  (parameter #47)
C     
C     Determine end of filename
 4047 I=ENDSTR( PLINE(FKC:LKC), LKC - FKC + 1) - 1
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ???????? line  (parameter #48)
C     
 4048 I=ENDSTR( PLINE(FKC:LKC), LKC - FKC + 1) - 1
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse MICCONUSEDATA line  (parameter #49)
C     
C     Parse the line:
 4049 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 4)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse MICOTHERDIVDATA line  (parameter #50)
C     
C     Parse the line:
 4050 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 12)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse MICOTHERSUPDATA line  (parameter #51)
C     
C     Parse the line:
 4051 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 8)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse MICRUNOFFDATA line  (parameter #52)
C     Not Implemented!
 4052 GOTO 4990 
C
C     -------------------      
C
C
C     Parse MICPRECIPDATA line  (parameter #53)
C     Not Implemented!
 4053 GOTO 4990 
C
C     -------------------      
C
C
C     Parse MICEVAPDATA line  (parameter #54)
C     Not Implemented!
 4054 GOTO 4990 
C
C     -------------------      
C
C
C     Parse ???????? line  (parameter #55)
C     
 4055 I=ENDSTR( PLINE(FKC:LKC), LKC - FKC + 1) - 1
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse HURCONUSEDATA line  (parameter #56)
C     
C     Parse the line:
 4056 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 5)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse HUROTHERDIVDATA line  (parameter #57)
C     
C     Parse the line:
 4057 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 13)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse HUROTHERSUPDATA line  (parameter #58)
C     
C     Parse the line:
 4058 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 9)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse HURRUNOFFDATA line  (parameter #59)
C     Not Implemented!
 4059 GOTO 4990 
C
C     -------------------      
C
C
C     Parse HURPRECIPDATA line  (parameter #60)
C     Not Implemented!
 4060 GOTO 4990 
C
C     -------------------      
C
C
C     Parse HUREVAPDATA line  (parameter #61)
C     Not Implemented!
 4061 GOTO 4990 
C
C     -------------------      
C
C
C     Parse STCSTARTLEVEL line  (parameter #62)
C     
C     Parse the line:
 4062 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,SCSTLV))GOTO 5000
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse STCNBSDATA line  (parameter #63)
C     
C     Parse the line:
 4063 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 2)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STCCONUSEDATA line  (parameter #64)
C     
C     Parse the line:
 4064 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 6)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STCOTHERDIVDATA line  (parameter #65)
C     
C     Parse the line:
 4065 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 14)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STCOTHERSUPDATA line  (parameter #66)
C     
C     Parse the line:
 4066 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 10)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STCRUNOFFDATA line  (parameter #67)
C     Not Implemented!
 4067 GOTO 4990 
C
C     -------------------      
C
C
C     Parse STCPRECIPDATA line  (parameter #68)
C     Not Implemented!
 4068 GOTO 4990 
C
C     -------------------      
C
C
C     Parse STCEVAPDATA line  (parameter #69)
C     Not Implemented!
 4069 GOTO 4990 
C
C     -------------------      
C
C
C     Parse STCEOPLEVELDATA line  (parameter #70)
C
C     Parse the line:
 4070 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 72)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ERISTARTLEVEL line  (parameter #71)
C     
C     Parse the line:
 4071 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,ERSTLV))GOTO 5000
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse ERINBSDATA line  (parameter #72)
C     
C     Parse the line:
 4072 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 3)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ERICONUSEDATA line  (parameter #73)
C     
C     Parse the line:
 4073 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 7)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C     Parse WELLANDCANALDATA line  (parameter #74)
C
C     Parse the line:
 4074 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 21)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse NYSBCANALDATA line  (parameter #75)
C     Not Implemented!
 4075 GOTO 4990 
C
C     -------------------      
C
C
C     Parse ERIOTHERDIVDATA line  (parameter #76)
C     
C     Parse the line:
 4076 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 15)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ERIOTHERSUPDATA line  (parameter #77)
C     
C     Parse the line:
 4077 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 11)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ERIRUNOFFDATA line  (parameter #78)
C     Not Implemented!
 4078 GOTO 4990 
C
C     -------------------      
C
C
C     Parse ERIPRECIPDATA line  (parameter #79)
C     Not Implemented!
 4079 GOTO 4990 
C
C     -------------------      
C
C
C     Parse ERIEVAPDATA line  (parameter #80)
C     Not Implemented!
 4080 GOTO 4990 
C
C     -------------------      
C 
C
C     Parse ERIEOPLEVELDATA line  (parameter #81)
C     
C     Parse the line:
 4081 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 73)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse CGIPSTARTLEVEL line  (parameter #82)
C     
C     Parse the line:
 4082 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,CGSTLV))GOTO 5000
      EMESSNO = 9999
      RETURN
C     
C     -------------------      
C
C
C     Parse CGIPAVGLEVELDATA line  (parameter #83)
C     
C     Parse the line:
 4083 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 20)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STMRIVFLOWDATA line  (parameter #84)
C     
C     Parse the line:
 4084 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 19)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STMRIVICEDATA line  (parameter #85)
C     
C     Parse the line:
 4085 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 51)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STMRIVWEEDDATA line  (parameter #86)
C     
C     Parse the line:
 4086 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 52)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STMRIVRETRDATA line  (parameter #87)
C     
C     Parse the line:
 4087 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 50)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ???????? line  (parameter #88)
C     
C     Determine end of filename
 4088 I=ENDSTR( PLINE(FKC:LKC), LKC - FKC + 1) - 1
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STCRIVFLOWDATA line  (parameter #89)
C     
C     Parse the line:
 4089 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 28)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STCRIVICEDATA line  (parameter #90)
C     
C     Parse the line:
 4090 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 34)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STCRIVWEEDDATA line  (parameter #91)
C     
C     Parse the line:
 4091 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 35)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STCRIVRETRDATA line  (parameter #92)
C     
C     Parse the line:
 4092 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 16)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ???????? line  (parameter #93)
C     
C     Determine end of filename
 4093 I=ENDSTR( PLINE(FKC:LKC), LKC - FKC + 1) - 1
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse DETRIVFLOWDATA line  (parameter #94)
C     
C     Parse the line:
 4094 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 29)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse DETRIVICEDATA line  (parameter #95)
C     
C     Parse the line:
 4095 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 36)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse DETRIVWEEDDATA line  (parameter #96)
C     
C     Parse the line:
 4096 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 37)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse DETRIVRETRDATA line  (parameter #97)
C     
C     Parse the line:
 4097 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 17)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ???????? line  (parameter #98)
C     
C     Determine end of filename
 4098 I=ENDSTR( PLINE(FKC:LKC), LKC - FKC + 1) - 1
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse NIARIVFLOWDATA line  (parameter #99)
C     
C     Parse the line:
 4099 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 30)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse NIARIVICEDATA line  (parameter #100)
C     
C     Parse the line:
 4100 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 38)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse NIARIVWEEDDATA line  (parameter #101)
C     
C     Parse the line:
 4101 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 39)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse NIARIVRETRDATA line  (parameter #102)
C     
C     Parse the line:
 4102 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 18)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ???????? line  (parameter #103)
C     
C     Determine end of filename
 4103 I=ENDSTR( PLINE(FKC:LKC), LKC - FKC + 1) - 1
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ONTSTARTLEVEL line  (parameter #104)
C     
C     Parse the line:
 4104 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,ONSTLV))GOTO 5000
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse ONTGATESDATA line  (parameter #67)
C     
C     Parse the line:
 4105 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 67)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ONTNBSDATA  line  (parameter #106)
C     
C     Parse the line:
 4106 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 53)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ONTCONUSEDATA line  (parameter #107)
C     
C     Parse the line:
 4107 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 54)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ONTOTHERDIVDATA line  (parameter #108)
C     
C     Parse the line:
 4108 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 55)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ONTOTHERSUPDATA line  (parameter #109)
C     
C     Parse the line:
 4109 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 56)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
 
C     -------------------      
C
C
C     Parse ONTRUNOFFDATA line  (parameter #110)
C     Not Implemented!
 4110 GOTO 4990 
C
C     -------------------      
C
C
C     Parse ONTPRECIPDATA line  (parameter #111)
C     Not Implemented!
 4111 GOTO 4990 
C
C     -------------------      
C
C
C     Parse ONTEVAPDATA line  (parameter #112)
C     Not Implemented!
 4112 GOTO 4990 
C
C     -------------------      
C
C
C     Parse ONTEOPLEVELDATA line  (parameter #113)
C     
C     Parse the line:
 4113 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 74)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ONTDEVIATIONDATA line  (parameter #114)
C     
C     Parse the line:
 4114 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 64)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ONTFORCEFLOWDATA line  (parameter #115)
C     
C     Determine end of filename
 4115 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 75)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse SUPAREAELEV line  (parameter #116)
C     
C     Parse the line:
 4116 AlreadySpecified(PNUM) = .FALSE. 
      IF(ReadAvsZ( 1, PLINE(FKC:LKC), LKC - FKC + 1)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse MHUAREAELEV line  (parameter #117)
C     
C     Parse the line:
 4117 AlreadySpecified(PNUM) = .FALSE. 
      IF(ReadAvsZ( 2, PLINE(FKC:LKC), LKC - FKC + 1)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse MICAREAELEV line  (parameter #118)
C     
C     Not Implemented!
 4118 GOTO 4990 
C     Parse the line:
C4118 IF(ReadAvsZ( 3, PLINE(FKC:LKC), LKC - FKC + 1)) GOTO 5000 
C
C     EMESSNO = 9999
C     RETURN
C
C     -------------------      
C
C
C     Parse HURAREAELEV line  (parameter #119)
C     
C     Not Implemented!
 4119 GOTO 4990 
C     Parse the line:
C4119 IF(ReadAvsZ( 4, PLINE(FKC:LKC), LKC - FKC + 1)) GOTO 5000 
C
C     EMESSNO = 9999
C     RETURN
C
C     -------------------      
C
C
C     Parse STCAREAELEV line  (parameter #120)
C     
C     Parse the line:
 4120 AlreadySpecified(PNUM) = .FALSE. 
      IF(ReadAvsZ( 5, PLINE(FKC:LKC), LKC - FKC + 1)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse ERIAREAELEV line  (parameter #121)
C     
C     Parse the line:
 4121 AlreadySpecified(PNUM) = .FALSE. 
      IF(ReadAvsZ( 6, PLINE(FKC:LKC), LKC - FKC + 1)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse ONTAREAELEV line  (parameter #122)
C     
C     Parse the line:
 4122 AlreadySpecified(PNUM) = .FALSE. 
      IF(ReadAvsZ( 7, PLINE(FKC:LKC), LKC - FKC + 1)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse FORECASTMONTHS line  (parameter #123)
C     
C     Parse the line:
 4123 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': months to forecast Lake Superior supplies '//LLT, 
     &   ' (0=no forecast)'//LLT,
     &    0, 12, 1183, ForecastMonths)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C           
C     Parse SUPOUTPUTFILES line  (parameter #124)
C     
 4124 IF(WhichOutputFiles( 1, PLINE(FKC:LKC), LKC - FKC + 1)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse MHUOUTPUTFILES line  (parameter #125)
C     
 4125 IF(WhichOutputFiles( 2, PLINE(FKC:LKC), LKC - FKC + 1)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse STCOUTPUTFILES line  (parameter #126)
C     
 4126 IF(WhichOutputFiles( 3, PLINE(FKC:LKC), LKC - FKC + 1)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse ERIOUTPUTFILES line  (parameter #127)
C     
 4127 IF(WhichOutputFiles( 4, PLINE(FKC:LKC), LKC - FKC + 1)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse ONTOUTPUTFILES line  (parameter #128)
C     
 4128 IF(WhichOutputFiles( 5, PLINE(FKC:LKC), LKC - FKC + 1)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse DPRRIVFLOWDATA line  (parameter #129)
C     
C     Parse the line:
 4129 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 65)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse LSLFLOWDATA line  (parameter #130)
C     
C     Parse the line:
 4130 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 66)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STLRIVFLOWDATA line  (parameter #131)
C     
C     Parse the line:
 4131 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 58)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STLRIVRETRDATA line  (parameter #132)
C     
C     Parse the line:
 4132 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 59)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STLRIVICEDATA line  (parameter #133)
C     
C     Parse the line:
 4133 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 60)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STLRIVWEEDDATA line  (parameter #134)
C     
C     Parse the line:
 4134 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 61)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse LSLRETRDATA line  (parameter #135)
C     
C     Parse the line:
 4135 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 62)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse LSLROUGHNESSDATA line  (parameter #136)
C     
C     Parse the line:
 4136 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 63)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse LSLONTDIFFDATA line  (parameter #137)
C     
C     Parse the line:
 4137 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 76)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse MIIRIVFLOWDATA line  (parameter #138)
C     
C     Parse the line:
 4138 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 77)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STMRIVSOLUTION line  (parameter #139)
C     
C     Determine end of string
 4139 IF( .NOT. FIELDBND(PLINE, LKC, FKC, SFNW, EFNW, I ) ) THEN 
C        Since there should be a useful character on line to make it this 
C        far, the only way to achieve this condition should be to specify 
C        a ',' as the next charcater after ':'.  dump if it happens anyway.  
         EMESSNO = 1059
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF 
C        
C     Strip out everything but letters and numbers and make uppercase to
C     simplify string comparisons (use CWRK to hold this version of 
C     parameter value).
      K=0
      DO 4566, J = 0, LKC - FKC
         I = ICHAR( PLINE(SFNW+J:SFNW+J) )
         IF( I .GE. 97 .AND. I .LE. 122) THEN
            K=K+1
            CWRK(K:K) = CHAR(I-32)
         ELSE IF( (I .GE. 48 .AND. I .LE. 57) .OR. 
     &             I .GE. 65 .AND. I .LE. 90) THEN
            K=K+1
            CWRK(K:K) = PLINE(SFNW+J:SFNW+J)
         ENDIF
 4566 CONTINUE
C     WRITE(*,'(A,A,A,I3)') "===>",CWRK(1:K),"<===   K=",K  
C
      IF(CWRK(1:4) .EQ. '1955' )THEN
         StMarysSolution = 2
      ELSE IF(CWRK(1:K) .EQ. '1977A' )THEN
         StMarysSolution = 3
      ELSE IF(CWRK(1:K) .EQ. 'PLAN2012' )THEN
         StMarysSolution = 4
      ELSE IF(CWRK(1:1) .EQ. 'O' )THEN
C        Using some type of outflow equation
         StMarysSolution = 1
      ELSE
         EMESSNO = 1117
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5029) PLINE(SFNW:EFNW), LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
      IF( GenVerbosity .GE. 3 ) THEN
         I = LASTKC(SMSolutionName(StMarysSolution), 16)
         WRITE(CTMP,5028) INDENT,PARCODE,
     &      ': Lake Superior outflow computed according to {',
     &      SMSolutionName(StMarysSolution)(1:I),'}',LLT
C     DATA SMSolutionName /   'outflow equation', '1955 Mod. Rule  ', 
C    &    'Plan 1977A      ', 'Plan 2012       '/
      ENDIF
C
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse SUPINC line  (parameter #140)
C     
C     Parse the line:
 4140 IF( PLINE(FKC:FKC) .EQ. 'H' .OR. PLINE(FKC:FKC) .EQ. 'h') THEN
C        Want it to be hourly; set to 0 for now and handle in CHeckPar()
         RoutingIncrement(1) = 0
         WRITE(CTMP,5028) INDENT,PARCODE,
     &      ': routing increments through Lake Superior set to {',
     &      'hourly','}',LLT
         GOTO 5000 
      ELSE 
         IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &      ': set number of increments within a period to '//LLT,
     &      ' for Lake Superior'//LLT, 1, 999,1060,RoutingIncrement(1))) 
     &      GOTO 5000 
      ENDIF
C
      EMESSNO = 9999
      RETURN
C     
C     -------------------      
C
C
C     Parse SUPTIMESTEP line  (parameter #141)
C     
C     Determine end of string
 4141 IF( .NOT. FIELDBND(PLINE, LKC, FKC, SFNW, EFNW, I ) ) THEN 
C        Since there should be a useful character on line to make it this 
C        far, the only way to achieve this condition should be to specify 
C        a ',' as the next charcater after ':'.  dump if it happens anyway.  
         EMESSNO = 1059
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF 
C         
      IF(PLINE(SFNW:SFNW) .EQ.'d' .OR. PLINE(SFNW:SFNW) .EQ.'D')THEN
         RoutingInterval(1) = 1
      ELSE IF(PLINE(SFNW:SFNW).EQ.'W' .OR. PLINE(SFNW:SFNW).EQ. 'w')THEN
         RoutingInterval(1) = 2
      ELSE IF(PLINE(SFNW:SFNW).EQ.'Q' .OR. PLINE(SFNW:SFNW).EQ. 'q')THEN
         RoutingInterval(1) = 3
      ELSE IF(PLINE(SFNW:SFNW).EQ.'m' .OR. PLINE(SFNW:SFNW).EQ. 'M')THEN
         RoutingInterval(1) = 4
      ELSE
         EMESSNO = 1088
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5027) PLINE(SFNW:SFNW), LLT
C         WRITE(CTMP,5027) PLINE(SFNW:SFNW), ''''//LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5028) INDENT,PARCODE,
     &      ': routing interval through Lake Superior set to {',
     &      IntervalName(RoutingInterval(1))
     &      (1:IntervalNameLen(RoutingInterval(1))), '}', LLT
      ENDIF
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STMRIVEQN line  (parameter #142)
C     
C     Parse the line:
 4142 IF(READQEQN( PLINE(FKC:LKC), LKC - FKC + 1, SPK, SPYM, SPA, 
     &   SPB, SPWT, SPC)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse PREVSTMRIVFLOW line  (parameter #143)
C     
C     Parse the line:
 4143 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1094, I, 2, TMPDBL)) THEN
         PrevStMarysFlow = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse PREVSUPGATESOPEN line  (parameter #144)
C     
C     Parse the line:
 4144 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1094, I, 0, TMPDBL)) THEN
         PrevSupGates = SNGL(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77ACOMPATMODE line  (parameter #145)
C     
 4145 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': compatibility with existing Plan 1977A practices set to '//
     &   'mode '//LLT, ' '//LLT, -1, 3, 1148, P77ACompat)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse ?????????????? line  (parameter #146)
C     
C     Determine end of filename
 4146 I=ENDSTR( PLINE(FKC:LKC), LKC - FKC + 1) - 1
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse BOMROUND line  (parameter #147)
C
 4147 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': during computations will round beginning-of-month levels'//
     &   ' to ten to the '//LLT, ' power'//LLT, -6, 3, 1116, 
     &   BOMRound)) 
     &   GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse FLOWROUND line  (parameter #148)
C     
 4148 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': during computations will round intermediate flows to ten'//
     &   ' to the '//LLT, ' power'//LLT, -6, 3, 1116, FlowRound)) 
     &   GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C     -------------------      
C
C
C     Parse LEVELROUND line  (parameter #149)
C     
 4149 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': during computations will round intermediate levels to ten'//
     &   ' to the '//LLT, ' power'//LLT, -6, 3, 1116, LevelRound)) 
     &   GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse ROUTINGTITLE line  (parameter #150)
C     
      
 4150 MLHEADER=PLINE(FKC:LKC)//LLT
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5017) INDENT, PARCODE, 
     &      ': Middle Lakes Routing Title {',
     &      MLHEADER(1:LKC-FKC+1),'}', LLT
      ENDIF
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STCRIVEQN line  (parameter #151)
C     
C     Parse the line:
 4151 IF(READQEQN( PLINE(FKC:LKC), LKC - FKC + 1, MHK, MHYM, MHA, 
     &   MHB, MHWT, MHC)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse DETRIVEQN line  (parameter #152)
C     
C     Parse the line:
 4152 IF(READQEQN( PLINE(FKC:LKC), LKC - FKC + 1, SCK, SCYM, SCA, 
     &   SCB, SCWT, SCC)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse NIARIVEQN line  (parameter #153)
C     
C     Parse the line:
 4153 IF(READQEQN( PLINE(FKC:LKC), LKC - FKC + 1, ERK, ERYM, ERA, 
     &   ERB, ERWT, ERC)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse MIDLAKESOLUTION line  (parameter #154)
C      
C     Determine end of string
 4154 IF( .NOT. FIELDBND(PLINE, LKC, FKC, SFNW, EFNW, I ) ) THEN 
C        Since there should be a useful character on line to make it this 
C        far, the only way to achieve this condition should be to specify 
C        a ',' as the next charcater after ':'.  dump if it happens anyway.  
         EMESSNO = 1059
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF 
C         
      IF(PLINE(SFNW:SFNW) .EQ.'m' .OR. PLINE(SFNW:SFNW) .EQ.'M')THEN
         MidLakeSolution = 2
         I = 6
      ELSE IF(PLINE(SFNW:SFNW).EQ.'I' .OR. PLINE(SFNW:SFNW).EQ. 'i')THEN
         MidLakeSolution = 1
         I = 9
      ELSE
         EMESSNO = 1092
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5027) PLINE(SFNW:SFNW), LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5028) INDENT,PARCODE,
     &      ': routing method through middle lakes set to {',
     &      MLSolutionName(MidLakeSolution)(1:I),'}',LLT
      ENDIF
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse MIDLAKEINC line  (parameter #155)
C     
C     Parse the line:
 4155 IF( PLINE(FKC:FKC) .EQ. 'H' .OR. PLINE(FKC:FKC) .EQ. 'h') THEN
C        Want it to be hourly; set to 0 for now and handle in CHeckPar()
         RoutingIncrement(2) = 0
         WRITE(CTMP,5028) INDENT,PARCODE,
     &      ': routing increments through middle lakes set to {',
     &      'hourly','}',LLT
         GOTO 5000 
      ELSE 
         IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &      ': set number of increments within a period to '//LLT,
     &      ' for Lakes Michigan-Erie'//LLT, 1, 999, 1060, 
     &      RoutingIncrement(2))) 
     &      GOTO 5000 
      ENDIF
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse MIDLAKETIMESTEP line  (parameter #156)
C     
C     Determine end of string
 4156 IF( .NOT. FIELDBND(PLINE, LKC, FKC, SFNW, EFNW, I ) ) THEN 
C        Since there should be a useful character on line to make it this 
C        far, the only way to achieve this condition should be to specify 
C        a ',' as the next charcater after ':'.  dump if it happens anyway.  
         EMESSNO = 1059
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF 
C         
      IF(PLINE(SFNW:SFNW) .EQ.'d' .OR. PLINE(SFNW:SFNW) .EQ.'D')THEN
         RoutingInterval(2) = 1
      ELSE IF(PLINE(SFNW:SFNW).EQ.'W' .OR. PLINE(SFNW:SFNW).EQ. 'w')THEN
         RoutingInterval(2) = 2
      ELSE IF(PLINE(SFNW:SFNW).EQ.'Q' .OR. PLINE(SFNW:SFNW).EQ. 'q')THEN
         RoutingInterval(2) = 3
      ELSE IF(PLINE(SFNW:SFNW).EQ.'m' .OR. PLINE(SFNW:SFNW).EQ. 'M')THEN
         RoutingInterval(2) = 4
      ELSE
         EMESSNO = 1088
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5027) PLINE(SFNW:SFNW), LLT
C         WRITE(CTMP,5027) PLINE(SFNW:SFNW), ''''//LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5028) INDENT,PARCODE,
     &      ': routing interval through middle lakes set to {',
     &      IntervalName(RoutingInterval(2))
     &      (1:IntervalNameLen(RoutingInterval(2))), '}', LLT
      ENDIF
      GOTO 5000 
C
C
C     -------------------      
C
C
C     Parse SUPSUMINPUTFILE line  (parameter #157)
C     
C     Parse the line:
 4157 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 41)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse MHUSUMINPUTFILE line  (parameter #158)
C     
C     Parse the line:
 4158 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 31)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse STCSUMINPUTFILE line  (parameter #159)
C     
C     Determine end of filename
 4159 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 32)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ERISUMINPUTFILE line  (parameter #160)
C     
C     Determine end of filename
 4160 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 33)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ONTSUMINPUTFILE line  (parameter #161)
C     
C     Determine end of filename
 4161 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 57)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse WRITEERIEOUTFLOW line  (parameter #162)
C     
C     Determine end of filename
 4162 IF( PLINE(FKC:FKC) .EQ. 'F' .OR. PLINE(FKC:FKC) .EQ. 'f') THEN 
C        If explicitly given as FALSE, set it FALSE.  This is 
C        default anyway.
         WriteErieOutflow = .FALSE.
      ELSE
C        If find this parameter and doesn't start with 'F', assume 
C        they want it set to TRUE.
         WriteErieOutflow = .TRUE.
      ENDIF
C
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5023) INDENT, PARCODE, 
     &      ': Generate Erie Outflow output file instead of '// 
     &      'Niagara River Flow file set to {',WriteErieOutflow,'}',LLT
      ENDIF
      GOTO 5000
C     -------------------      
C
C
C     Parse USTLRIVSOLUTION line  (parameter #163)
C      DATA PARIDS(2593:2608) /'  '/                        163
C     
C     Determine end of string
 4163 IF( .NOT. FIELDBND(PLINE, LKC, FKC, SFNW, EFNW, I ) ) THEN 
C        Since there should be a useful character on line to make it this 
C        far, the only way to achieve this condition should be to specify 
C        a ',' as the next charcater after ':'.  dump if it happens anyway.  
         EMESSNO = 1059
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF 
C         
      IF(PLINE(SFNW:SFNW) .EQ.'o' .OR. PLINE(SFNW:SFNW) .EQ.'O')THEN
C        Using some type of outflow equation for Lake Ontario outflow
         UStLawrSolution = 1
         I = 16
      ELSE IF(PLINE(SFNW:SFNW+4).EQ.'1958D' .OR. 
     &   PLINE(SFNW:SFNW+4).EQ. '1958d')THEN
C        Using 1958D for Lake Ontario outflow
         UStLawrSolution = 2
         I = 10
      ELSE IF(PLINE(SFNW:SFNW) .EQ.'i'.OR. PLINE(SFNW:SFNW) .EQ.'I')THEN
C        Using IS Plan for Lake Ontario outflow
         UStLawrSolution = 3
         I = 7
      ELSE IF(PLINE(SFNW:SFNW+3).EQ.'1998' )THEN 
C        Using 1998 for Lake Ontario outflow
         UStLawrSolution = 4
         I = 9
      ELSE IF(PLINE(SFNW:SFNW+2).EQ.'Pre' .OR. PLINE(SFNW:SFNW+2)
     &   .EQ.'pre' .OR. PLINE(SFNW:SFNW+2).EQ. 'PRE')THEN
C        Using "Preproject" conditions for Lake Ontario outflow
         UStLawrSolution = 5
         I = 10
      ELSE
         EMESSNO = 1134
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5029) PLINE(SFNW:EFNW), LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5028) INDENT,PARCODE,
     &      ': Lake Ontario outflow computed according to {',
     &      USLSolutionName(UStLawrSolution)(1:I),'}',LLT
      ENDIF
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse ONTINC line  (parameter #164)
C     
C     Parse the line:
 4164 IF( PLINE(FKC:FKC) .EQ. 'H' .OR. PLINE(FKC:FKC) .EQ. 'h') THEN
C        Want it to be hourly; set to 0 for now and handle in CHeckPar()
         RoutingIncrement(3) = 0
         GOTO 5000 
      ELSE 
         IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &      ': set number of increments within a period to '//LLT,
     &      ' for Lake Ontario'//LLT, 1, 999,1060, RoutingIncrement(3))) 
     &      GOTO 5000 
      ENDIF
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse ONTTIMESTEP line  (parameter #165)
C     
C     Determine end of string
 4165 IF( .NOT. FIELDBND(PLINE, LKC, FKC, SFNW, EFNW, I ) ) THEN 
C        Since there should be a useful character on line to make it this 
C        far, the only way to achieve this condition should be to specify 
C        a ',' as the next charcater after ':'.  dump if it happens anyway.  
         EMESSNO = 1059
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF 
C         
      IF(PLINE(SFNW:SFNW) .EQ.'d' .OR. PLINE(SFNW:SFNW) .EQ.'D')THEN
         RoutingInterval(3) = 1
      ELSE IF(PLINE(SFNW:SFNW).EQ.'W' .OR. PLINE(SFNW:SFNW).EQ. 'w')THEN
         RoutingInterval(3) = 2
      ELSE IF(PLINE(SFNW:SFNW).EQ.'Q' .OR. PLINE(SFNW:SFNW).EQ. 'q')THEN
         RoutingInterval(3) = 3
      ELSE IF(PLINE(SFNW:SFNW).EQ.'m' .OR. PLINE(SFNW:SFNW).EQ. 'M')THEN
         RoutingInterval(3) = 4
      ELSE
         EMESSNO = 1088
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5027) PLINE(SFNW:SFNW), LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5028) INDENT,PARCODE,
     &      ': routing interval through Lake Ontario set to {',
     &      IntervalName(RoutingInterval(3))
     &      (1:IntervalNameLen(RoutingInterval(3))), '}', LLT
      ENDIF
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse USTLRIVEQN line  (parameter #166)
C     
C     Parse the line:
 4166 IF(READQEQN( PLINE(FKC:LKC), LKC - FKC + 1, ONK, ONYM, ONA, 
     &   ONB, ONWT, ONC)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse SUPINDEX1DATA line  (parameter #167)
C     
C     Determine end of filename
 4167 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 93)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse SUPINDEX2DATA line  (parameter #168)
C     
C     Determine end of filename
 4168 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 94)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse SUPINDEX3DATA line  (parameter #169)
C     
C     Determine end of filename
 4169 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 95)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse P77ASUPTARGETLEV line  (parameter #170)
C     
 4170 IF(ReadFloatList(PLINE(FKC:LKC), LKC - FKC + 1, 12, 
     &   P77ASupTargetLev)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77AMHUTARGETLEV line  (parameter #171)
C     
 4171 IF(ReadFloatList(PLINE(FKC:LKC), LKC - FKC + 1, 12, 
     &   P77AMHTargetLev)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77ASUPRANGEFACT line  (parameter #172)
C     
 4172 IF(ReadFloatList(PLINE(FKC:LKC), LKC - FKC + 1, 12, 
     &   P77ASupRangeFact)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77AMHURANGEFACT line  (parameter #173)
C     
 4173 IF(ReadFloatList(PLINE(FKC:LKC), LKC - FKC + 1, 12, 
     &   P77AMHRangeFact)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77ABASEOUTFLOW line  (parameter #174)
C     
 4174 IF(ReadIntegerList(PLINE(FKC:LKC), LKC - FKC + 1, 12, 
     &   P77ABaseFlow)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     SUPERIOR LOWER LIMIT FLOWS        ISL   P77AMINIMUMFLOW   175
C     Parse P77AMINIMUMFLOW line  (parameter #175)
C     
C      LOGICAL FUNCTION ReadIntegerList(InLine, LineLength, 
C     &   ValuesExpected, IntArray)
C     Parse the line:
 4175 IF(ReadIntegerList(PLINE(FKC:LKC), LKC - FKC + 1, 12, 
     &   P77AMinimumFlow)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77ALLOGOKIDIV line  (parameter #176)
C     
C     Not Implemented!
 4176 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77ACHICAGODIV line  (parameter #177)
C     
C     Not Implemented!
 4177 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77AWELLANDDIV line  (parameter #178)
C     
C     Not Implemented!
 4178 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77ASTCRIVRETR line  (parameter #179)
C     
C     Not Implemented!
 4179 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77ADETRIVRETR line  (parameter #180)
C     
C     Not Implemented!
 4180 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77ANIARIVRETR line  (parameter #181)
C     
C     Not Implemented!
 4181 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77ASUP5NBS line  (parameter #182)
C     
C     Not Implemented!
 4182 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77ASUP50NBS line  (parameter #183)
C     
C     Not Implemented!
 4183 GOTO 4990 
C
C     -------------------      
C 681 FORMAT ( 12F6.0 )                                                     4020
C     SUPERIOR TARGET STAGES            PSM   P77ASUPTARGETLEV  170
C     SUPERIOR STANDARD DEVATIONS       PSSD  P77AMHTARGETLEV   171
C     MICHIGAN-HURON TARGET STAGES      PHM   P77ASUPRANGEFACT  172
C     MICHIGAN-HURON STANDARD DEVATIONS PHSD  P77AMHRANGEFACT   173
C     SUPERIOR MEAN FLOWS               QB    P77ABASEOUTFLOW   174
C6001 FORMAT ( 9X,I4,7X,12I5 )                                              4260
C     SUPERIOR LOWER LIMIT FLOWS        ISL   P77AMINIMUMFLOW   175
C     LLO DIV                           JNSUP P77ALLOGOKIDIV    176
C     CHICAGO DIV                       JCD   P77ACHICAGODIV    177
C     WELLAND DIV                       JWEL  P77AWELLANDDIV    178
C     STMARYS RETARDATION               JICEM P77ASTMARYSRETRD  179
C     STCLAIR RETARDATION               JICEC P77ASTCLAIRRETRD  180
C     NIAGARA RETARDATION               JICEE P77ANIAGARARETRD  181
C     SUPERIOR 5% NBS                   NS5   P77ASUP5NBS       182
C     SUPERIOR 50% NBS                  NS50  P77ASUP50NBS      183
C     SUPERIOR 95% NBS                  NS95  P77ASUP95NBS      184
C     MICHIGAN-HURON 5% NBS             NH5   P77AMH5NBS        185
C     MICHIGAN-HURON 50% NBS            NH50  P77AMH50NBS       186
C     MICHIGAN-HURON 95% NBS            NH95  P77AMH95NBS       187
C     ST CLAIR 50% NBS                  NC50  P77ASTC50NBS      188
C     ERIE 50% NBS                      NE50  P77AERIE50NBS     189
C
C
C     Parse P77ASUP95NBS line  (parameter #184)
C     
C     Not Implemented!
 4184 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77AMHU5NBS line  (parameter #185)
C     
C     Not Implemented!
 4185 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77AMHU50NBS line  (parameter #186)
C     
C     Not Implemented!
 4186 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77AMHU95NBS line  (parameter #187)
C     
C     Not Implemented!
 4187 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77ASTC50NBS line  (parameter #188)
C     
C     Not Implemented!
 4188 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77AERI50NBS line  (parameter #189)
C     
C     Not Implemented!
 4189 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77ANIARIVEQN line  (parameter #190)
C     
C     Parse the line:
 4190 IF(READQEQN( PLINE(FKC:LKC), LKC - FKC + 1, P77AERK, P77AERYM, 
     &   P77AERA, P77AERB, P77AERWT, P77AERC)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77ADETRIVEQN line  (parameter #191)
C     
C     Parse the line:
 4191 IF(READQEQN( PLINE(FKC:LKC), LKC - FKC + 1, P77ASCK, P77ASCYM,
     &   P77ASCA, P77ASCB, P77ASCWT, P77ASCC)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77ASTCRIVEQN line  (parameter #192)
C     
C     Parse the line:
 4192 IF(READQEQN( PLINE(FKC:LKC), LKC - FKC + 1, P77AMHK, P77AMHYM, 
     &   P77AMHA, P77AMHB, P77AMHWT, P77AMHC)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77ASTMRIVEQN line  (parameter #193)
C     
C     Parse the line:
 4193 IF(READQEQN( PLINE(FKC:LKC), LKC - FKC + 1, P77ASPK, P77ASPYM, 
     &   P77ASPA, P77ASPB, P77ASPWT, P77ASPC)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77AMIDLAKEINC line  (parameter #194)
C     
 4194 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': set number of increments within a period to '//LLT,
     &   ' for middle lakes under Plan 1977A'//LLT, 
     &   1, 999, 1060, P77AMidLakesInc)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C     -------------------      
C
C
C     Parse P77AMIDLAKETOL line  (parameter #195)
C     
C     Not Implemented!
 4195 GOTO 4990 
C
C      DATA P77AMidLakesInc, P77AMidLakesTol / 40, -2 /
C     -------------------      
C
C
C     Parse P77ACGIPLEVEL line  (parameter #196)
C      DATA P77ACGIP / 171.16D0 /
C     
C     Not Implemented!
 4196 GOTO 4990 
C
C     -------------------      
C
C
C     Parse P77AMONTHROUND line  (parameter #197)
C     
C     Parse the line:
 4197 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': during Plan 1977A computations will round '//
     &   'end-of-month lake levels to ten to the '//LLT, 
     &   ' power'//LLT, -6, 3, 1116, P77AMonthRound)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77ALEVELROUND line  (parameter #198)
C
C     Parse the line:
 4198 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': during Plan 1977A computations will round '//
     &   'intermediate lake levels to ten to the '//LLT, 
     &   ' power'//LLT, -6, 3, 1116, P77ALevelRound)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77AFLOWROUND line  (parameter #199)
C     
C     Parse the line:
 4199 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': during Plan 1977A computations will round '//
     &   'intermediate outflows to ten to the '//LLT, 
     &   ' power'//LLT, -6, 3, 1116, P77AFlowRound)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77AMONTHSTOAVG line  (parameter #200)
C     
C     Read values into MonthsToAverage
C
 4200 FKC = FKC - 1
C
C     Read each field 
      DO 3200, I = 1, 12
C
C        Save the begnning of the next field in J.
         J = FKC + 1
C
         IF( .NOT. FIELDBND(PLINE, LKC, J, SFNW, EFNW, FKC ) ) GOTO 2200
C
         READ( PLINE(SFNW:EFNW), 5012, ERR=2200, IOSTAT=K ) 
     &      P77AMonthsToAve(I)
         IF( P77AMonthsToAve(I) .LE. 0 .OR. 
     &       P77AMonthsToAve(I) .GT. 13) THEN
            IF(K .NE. 0) GOTO 2200
            EMESSNO = 1136
            GOTO 1200
         ENDIF
 3200    CONTINUE
C
      IF( GenVerbosity .GE. 5 ) THEN
         DO 5200, K=1,12
 5200       WRITE(CWRK((K-1)*3+1:K*3),'(I2,A1)') P77AMonthsToAve(K),','
         WRITE(CTMP,5028) INDENT,PARCODE,
     &      ': Plan1977A months to average set to {',   
     &      CWRK(1:35), '}', LLT
      ENDIF
C
      GOTO 5000 
C
 2200 EMESSNO = 1135
 1200 CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
      IF( K. NE. 0) THEN 
         WRITE(CTMP,'(A,I8,A,A1)')' Call to READ returned ', 
     &       K, '.', LLT
         CALL LogIt(LOGCDE,0,0)
      ENDIF
      WRITE(CTMP,'(A,I2,A,A1)')' Problem occurred reading field ',
     &    I, '.', LLT
      CALL LogIt(LOGCDE,0,0)
      RETURN
C
C
C     -------------------      
C
C
C     Parse P77AACONSTANT line  (parameter #201)
C
 4201 IF(ReadValue(
     &   PLINE(FKC:LKC),LKC-FKC+1,1059,I,0,P77AAValue))GOTO 5000
      EMESSNO = 9999
      RETURN
C     
C     -------------------      
C
C
C     Parse P77AMAXCHANGE line  (parameter #202)
C     
 4202 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 2, TMPDBL)) THEN
         P77AMaxChange = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77ALOWTRIGGER line (parameter #203)
C     
 4203 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,
     &   P77ALowTrigger))GOTO 5000
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P77AHIGHTRIGGER line (parameter #204)
C     
 4204 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,
     &   P77AHighTrigger))GOTO 5000
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse SUPSWPIERFORMULA line  (parameter #205)
C     
C     Parse the line:
 4205 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': using Lake Superior:SW Pier relationship '//LLT, 
     &   ' (1=old method,2=2000 EC equation,3=2008 EC equation)'//LLT,
     &    1, 3, 1170, LStoSWP)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse SUPCRITERIAROUND line  (parameter #206)
C     
C     Parse the line:
 4206 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': rounding for computations regarding Superior criterions '//
     &   'set to '//LLT, ' (1=like Plan 1977A, 2=general setting)'//LLT,
     &    1, 2, 1172, SupCritRnd)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse SUPBOMLEVEL line  (parameter #207)
C     
 4207 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,SPBOMLev))GOTO 5000
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse MHUBOMLEVEL line  (parameter #208)
C     
 4208 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,MHBOMLev))GOTO 5000
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse STCBOMLEVEL line  (parameter #209)
C     
 4209 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,SCBOMLev))GOTO 5000
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse ERIBOMLEVEL line  (parameter #210)
C     
 4210 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,ERBOMLev))GOTO 5000
      EMESSNO = 9999
      RETURN
C
C
C     -------------------      
C
C
C     Parse SUPPREPROJMETHOD line  (parameter #211)
C     
 4211 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': Superior pre-project flow determined using method'//LLT,
     &   ' (1=old way, 2=new way)'//LLT,1,2,1177,Q1887Method)) GOTO 5000
     &    
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse SUPCRITBMETHOD line  (parameter #212)
C     
 4212 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': Superior Criterion B flow determined using method'//LLT,
     &   ' (1=old way, 2=new way)'//LLT,1,2,1178,CritBMethod)) GOTO 5000
     &    
C
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse SUPCRITCTHRESH line  (parameter #213)
C     
C     Determine end of filename
 4213 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,CritCThresholdEl))
     &    GOTO 5000
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse Compute Montreal line  (parameter #214)
C     
 4214 IF( PLINE(FKC:FKC) .EQ. 'F' .OR. PLINE(FKC:FKC) .EQ. 'f') THEN 
C        If explicitly given as FALSE, set it FALSE.  This is 
C        default anyway, but what the heck.
         ComputeMontreal = .FALSE.
      ELSE
C        If find this parameter and doesn't start with 'F', assume 
C        they want it set to TRUE.
         ComputeMontreal = .TRUE.
      ENDIF
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5023) INDENT, PARCODE, 
     &      ': Perform St. Lawrence computations downstream of '//
     &      'Montreal set to {',DoDSMontrealComp, '}', LLT
      ENDIF
      GOTO 5000
C
C     -------------------      
C
C
C     Parse SEAWAYALERTDEPTH line  (parameter #215)
C     
 4215 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 0, TMPDBL)) THEN
         CritDepths = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P58DELIMINATEDEV line  (parameter #216)
C     
 4216 IF( PLINE(FKC:FKC) .EQ. 'F' .OR. PLINE(FKC:FKC) .EQ. 'f') THEN 
C        If explicitly given as FALSE, set it FALSE.  This is 
C        default anyway, but what the heck.
         ElimDev = .FALSE.
      ELSE
C        If find this parameter and doesn't start with 'F', assume 
C        they want it set to TRUE.
         ElimDev = .TRUE.
      ENDIF
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5023) INDENT, PARCODE, 
     &      ': Plan 58D eliminate deviations set to {',
     &      ElimDev, '}', LLT
      ENDIF
      GOTO 5000
C
C     -------------------      
C
C
C     Parse P58DFORCEFLOW line  (parameter #217)
C     
 4217 IF( PLINE(FKC:FKC) .EQ. 'F' .OR. PLINE(FKC:FKC) .EQ. 'f') THEN 
C        If explicitly given as FALSE, set it FALSE.  This is 
C        default anyway, but what the heck.
         ForceFlow = .FALSE.
      ELSE
C        If find this parameter and doesn't start with 'F', assume 
C        they want it set to TRUE.
         ForceFlow = .TRUE.
      ENDIF
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5023) INDENT, PARCODE, 
     &      ': ignore above levels for forced flows below set to {',
     &      ForceFlow, '}', LLT
      ENDIF
      GOTO 5000
C
C     -------------------      
C
C
C     Parse P58DTOTALSUPPLY line  (parameter #218)
C     
 4218 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 2, TMPDBL)) THEN
         P58DTotalSup = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse USELACTLEV line  (parameter #219)
C     
 4219 IF( PLINE(FKC:FKC) .EQ. 'F' .OR. PLINE(FKC:FKC) .EQ. 'f') THEN 
C        If explicitly given as FALSE, set it FALSE.  This is 
C        default anyway, but what the heck.
         SetFlow = .FALSE.
      ELSE
C        If find this parameter and doesn't start with 'F', assume 
C        they want it set to TRUE.
         SetFlow = .TRUE.
      ENDIF
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5023) INDENT, PARCODE, 
     &      ': Plan 58D set flow=plan flow if total devs < 0 set to {',
     &      SetFlow, '}', LLT
      ENDIF
      GOTO 5000
C
C     -------------------      
C
C
C     Parse P58DPREPROJCOND line  (parameter #220)
C
 4220 IF( PLINE(FKC:FKC) .EQ. 'F' .OR. PLINE(FKC:FKC) .EQ. 'f') THEN 
C        If explicitly given as FALSE, set it FALSE.  This is 
C        default anyway, but what the heck.
         PreProjConditions = .FALSE.
      ELSE
C        If find this parameter and doesn't start with 'F', assume 
C        they want it set to TRUE.
         PreProjConditions = .TRUE.
      ENDIF
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5023) INDENT, PARCODE, 
     &      ': Plan 58D apply pre-project conditions (if 74.15 < '//
     &      'level < 75.37) set to {', PreProjConditions, '}', LLT
      ENDIF
      GOTO 5000
C
C
C     -------------------      
C
C
C     Parse P58DMAXILIMIT line  (parameter #221)
C     
 4221 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 2, TMPDBL)) THEN
         MaxILimit = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P58DMINPLIMIT line  (parameter #222)
C     
 4222 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 2, TMPDBL)) THEN
         MinPLimit = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P58DACTUALLEVEL line  (parameter #223)
C     
 4223 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,
     &   P58DActualLevel))GOTO 5000
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P58DKOFACTOR line  (parameter #224)
C     
 4224 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 2, TMPDBL)) THEN
         P58DWeightedSupply = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse CRITERIONKFILE line  (parameter #225)
C     
C     Determine end of filename
 4225 I=LastKC(PLINE(FKC:LKC), LKC - FKC + 1) - 1
C
C     Now see if the file exists:
      INQUIRE(FILE=PLINE(FKC:FKC+I),EXIST=LOGVAR)
      IF(.NOT. LOGVAR) THEN
         EMESSNO = 1182
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5020)
     &      ' Filename "',PLINE(FKC:FKC+I),'" read from line:',LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C     
C     Ok, file is there (we'll assume that it's readable and intact)
C     Erase current contents of CritKFile and put in new name:
      DO 3225, J=1,256
 3225    CritKFile(J:J) = ' ' 
      CritKFile(1:I+1)=PLINE(FKC:FKC+I)
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5017) INDENT,PARCODE,
     &      ': set Criterion K File name to {',
     &      CritKFile(1:I+1),'}',LLT
      ENDIF
      GOTO 5000 
C
C
C     -------------------      
C
C
C     Parse P58DAPPLIMIT line  (parameter #226)
C     
 4226 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 2, TMPDBL)) THEN
         P58DAppLimit = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P58DSUPINDICATOR line  (parameter #227)
C     
c4227 IF(ReadIntegerList(PLINE(FKC:LKC), LKC - FKC + 1, 20,
c    &   SupIndReadArray)) GOTO 5000 
 4227 DO 3227, K = 1, 20 
 3227    IF(SupIndReadArray(K) .LE. -9999) GOTO 2227
      EMESSNO = 1179
      CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
      RETURN     
 2227 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 0, TMPDBL)) THEN
         J = IDNINT(TMPDBL)
C        First value (J) must be number of values on the line
         IF( J .LT. 1 .OR. J .GT. 20) THEN
            EMESSNO = 1180
            CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
            WRITE(CTMP,'(A,I7,A1)') ' Bogus number of values on line: ',
     &          J, LLT
            CALL LogIt(LOGCDE,0,0)
            RETURN     
         ENDIF
         IF(ReadIntegerList(PLINE(FKC + I:LKC), LKC - FKC - I + 1, J,
     &      SupIndReadArray(K))) THEN
            IF( GenVerbosity .GE. 5 ) THEN
               WRITE(CTMP,5026) INDENT, PARCODE, 
     &            ': Plan 58D supply indicators set up to element {',
     &            J + K - 1, '}', LLT
            ENDIF
            AlreadySpecified(PNUM) = .FALSE.
            GOTO 5000 
         ENDIF
      ENDIF
C
      EMESSNO = 9999
      RETURN
C
C
C     -------------------      
C
C
C     Parse P58DRUNTYPE line  (parameter #228)
C     
C     Determine end of filename
 4228 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 0, TMPDBL)) THEN
         P58DRunType = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C     Parse P58DFORECASTDATE line  (parameter #229)
C     
C     Read ending year into variable EndYear
 4229 IF( .NOT. FIELDBND(PLINE, LKC, FKC, SFNW, EFNW, I ) ) THEN 
C        Since there should be a useful character on line to make it this 
C        far, the only way to achieve this condition should be to specify 
C        a ',' as the next charcater after ':'.  dump if it happens anyway.  
         EMESSNO = 1027
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF 
C     WRITE(*,'(A,A,A)')' trimmed year field ::', PLINE(SFNW:EFNW),'::'
C     Save the begnning of the next field, if any, in J.
      J = I + 1
C
      READ( PLINE(SFNW:EFNW), 5012, ERR=3016, IOSTAT=I) P58DForecastYear
      IF( P58DForecastYear .LE. 0 .OR. P58DForecastYear .GE. 99999) THEN
         EMESSNO = 1021
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5024) P58DForecastYear, LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
C
      IF( .NOT. FIELDBND(PLINE, LKC, J, SFNW, EFNW, I ) ) THEN 
C        No delimiter, or no more data.
         EMESSNO = 1149
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF 
C     WRITE(*,'(A,A,A)')' trimmed month field ::', PLINE(SFNW:EFNW),'::'
C     Save the begnning of the next field, if any, in J.
      J = I + 1
C
      READ( PLINE(SFNW:EFNW), 5012, ERR=2016, IOSTAT=I ) P58DForecastMon
      IF( P58DForecastMon .LE. 0 .OR. P58DForecastMon .GT. 12) THEN
         EMESSNO = 1024
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5024) P58DForecastMon, LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
C     Find # of days in month, store in FKC.  Use FKC and I as scratch
      CALL NUMDAYS (P58DForecastYear, P58DForecastMon, I, FKC)
C
      IF( .NOT. FIELDBND(PLINE, LKC, J, SFNW, EFNW, I ) ) THEN 
C        No delimiter, or no more data.  Presume no day specified on line.
         EMESSNO = 1149
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         RETURN
      ENDIF
C
C     
      READ( PLINE(SFNW:EFNW), 5012, ERR=1016, IOSTAT=I) P58DForecastDay 
      IF( P58DForecastDay  .LE. 0 .OR. P58DForecastDay .GT. FKC) THEN
         EMESSNO = 1026
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5025) P58DForecastDay , P58DForecastMon, 
     &       P58DForecastYear, FKC, LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5030) INDENT, PARCODE, 'P58DForecastDay', 
     &      P58DForecastDay, MonthAbbrev(P58DForecastMon), 
     &      P58DForecastYear, LLT
      ELSE
         CTMP(1:1)=LLT
      ENDIF
      GOTO 5000 
C
      EMESSNO = 1029
      GOTO 0029
C2229 EMESSNO = 1036
      GOTO 0029
C1229 EMESSNO = 1038
 0029 CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
      WRITE(CTMP,5022)I,LLT
      CALL LogIt(LOGCDE,0,0)
      CALL BuildInfo
      CALL LogIt(LOGCDE,0,1)
      RETURN





C     -------------------      
C
C
C     Parse P58DSIMULATEFILE line  (parameter #230)
C     
C     Determine end of filename
 4230 I=LastKC( PLINE(FKC:LKC), LKC - FKC + 1) - 1
C
C     Now see if the file exists:
      INQUIRE(FILE=PLINE(FKC:FKC+I),EXIST=LOGVAR)
      IF(.NOT. LOGVAR) THEN
         EMESSNO = 1150
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5020)
     &      ' Filename "',PLINE(FKC:FKC+I),'" read from line:',LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C     
C     Ok, file is there (we'll assume that it's readable and intact)
C     Erase current contents of P58DSIMULATEFILE and put in new name:
      DO 3230, J=1,256
 3230    P58DSIMULATEFILE(J:J) = ' '
      P58DSIMULATEFILE(1:I+1)=PLINE(FKC:FKC+I)
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5017) INDENT,PARCODE,
     &      ': set Plan 1958D Simulate filename to {',
     &      P58DSIMULATEFILE(1:I+1),'}',LLT
      ENDIF
      GOTO 5000 
C
C
C     -------------------      
C
C
C     Parse P58DAPPLIMITSYM line  (parameter #231)
C     
 4231 P58DAppLimitSym = PLINE(FKC:FKC + 1)
C
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5017) INDENT, PARCODE, 
     &      ': Plan 1958D App Limit Sym set to {',
     &      P58DAppLimitSym(1:LASTKC(P58DAppLimitSym,2)), '}', LLT
      ENDIF
      GOTO 5000

C     -------------------      
C
C
C     Parse P58DPLANLEVEL line  (parameter #232)
C     
 4232 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,
     &   P58DPlanLevel))GOTO 5000
      EMESSNO = 9999
      RETURN
C4232 IF( PLINE(FKC:FKC) .EQ. 'F' .OR. PLINE(FKC:FKC) .EQ. 'f') THEN 
C        If explicitly given as FALSE, set it FALSE.  This is 
C        default anyway, but what the heck.
C        IgnoreDeviations = .FALSE.
C     ELSE
C        If find this parameter and doesn't start with 'F', assume 
C        they want it set to TRUE.
C        IgnoreDeviations = .TRUE.
C     ENDIF
C     IF( GenVerbosity .GE. 5 ) THEN
C        WRITE(CTMP,5023) INDENT, PARCODE, 
C    &      ': Plan 58D Ignore deviations set to {',
C    &      IgnoreDeviations, '}', LLT
C     ENDIF
C     GOTO 5000
C
C     -------------------      
C
C
C     Parse P58DTotalDev line  (parameter #233)
C     
 4233 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 2, TMPDBL)) THEN
         P58DTotalDev = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P58DAdjustment line  (parameter #234)
C     
 4234 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 2, TMPDBL)) THEN
         P58DAdjustment = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C     Parse P58DrcAdjustment line  (parameter #235)
C     
 4235 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 2, TMPDBL)) THEN
         P58DrcAdjustment = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C     Parse P58DActualflow line  (parameter #236)
C     
 4236 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 2, TMPDBL)) THEN
         P58DActualflow = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
c
C     Parse P58DUseDatabase line  (parameter #237)
C     
 4237 IF( PLINE(FKC:FKC) .EQ. 'F' .OR. PLINE(FKC:FKC) .EQ. 'f') THEN 
C        If explicitly given as FALSE, set it FALSE.  This is 
C        default anyway, but what the heck.
         P58DUseDatabase = .FALSE.
      ELSE
C        If find this parameter and doesn't start with 'F', assume 
C        they want it set to TRUE.
         P58DUseDatabase = .TRUE.
      ENDIF
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5023) INDENT, PARCODE, 
     &      ': Use the database for starting INITIAL CONDITIONS {',
     &      P58DUseDatabase, '}', LLT
      ENDIF
      GOTO 5000
C
C     -------------------      
C
C
c
C     Parse P58DUseDatabase line  (parameter #238)
C     
 4238 IF( PLINE(FKC:FKC) .EQ. 'F' .OR. PLINE(FKC:FKC) .EQ. 'f') THEN 
C        If explicitly given as FALSE, set it FALSE.  This is 
C        default anyway, but what the heck.
         P58DSiDatabase = .FALSE.
      ELSE
C        If find this parameter and doesn't start with 'F', assume 
C        they want it set to TRUE.
         P58DSiDatabase = .TRUE.
      ENDIF
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5023) INDENT, PARCODE, 
     &      ': Use the database for starting INITIAL CONDITIONS {',
     &      P58DSiDatabase, '}', LLT
      ENDIF
      GOTO 5000 
C
C     -------------------      
C
C     Parse PREPROJECTLEVEL line  (parameter #239)
C     
 4239 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,
     &   PreprojectLevel))GOTO 5000
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse P58DSIDATAFILE line  (parameter #240)
C     
C     Determine end of filename
 4240 I=LastKC( PLINE(FKC:LKC), LKC - FKC + 1) - 1
C
C     Now see if the file exists:
      INQUIRE(FILE=PLINE(FKC:FKC+I),EXIST=LOGVAR)
      IF(.NOT. LOGVAR) THEN
         EMESSNO = 1150
         CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
         WRITE(CTMP,5020)
     &      ' Filename "',PLINE(FKC:FKC+I),'" read from line:',LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C     
C     Ok, file is there (we'll assume that it's readable and intact)
C     Erase current contents of P58DSIDataFile and put in new name:
      DO 3240, J=1,256
 3240    P58DSIDataFile(J:J) = ' '
      P58DSIDataFile(1:I+1)=PLINE(FKC:FKC+I)
      IF( GenVerbosity .GE. 5 ) THEN
         WRITE(CTMP,5017) INDENT,PARCODE,
     &      ': set Plan 1958D Supply Indicator DataBase filename to {',
     &      P58DSIDataFile(1:I+1),'}',LLT
      ENDIF
      GOTO 5000 
C
C     -------------------      
C
C     Parse P58DLSLflow line  (parameter #241)
C     
 4241 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 2, TMPDBL)) THEN
         P58DLSLflow = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C
C     -------------------      
C
C
C     Parse Des Prairies & Mille Iles Flow line  (parameter #242)
C     
C     Parse the line:
 4242 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 78)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse St. Francois Flow line  (parameter #243)
C
C     Parse the line:
 4243 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 79)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse St. Maurice Flow line  (parameter #244)
C     
C     Parse the line:
 4244 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 80)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse Richelieu Flow line  (parameter #245)
C     
C     Parse the line:
 4245 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 81)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse Roughness at Jetty1 line  (parameter #246)
C     
C     Parse the line:
 4246 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 82)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse Roughness at Varennes line  (parameter #247)
C     
C     Parse the line:
 4247 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 83)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse Roughness at Sorel line  (parameter #248)
C     
C     Parse the line:
 4248 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 84)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse Roughness at Lac St. Pierre line  (parameter #249)
C     
C     Parse the line:
 4249 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 85)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C
C     Parse Roughness at 3 Rivers line  (parameter #250)
C     
C     Parse the line:
 4250 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 86)) GOTO 5000
      EMESSNO = 9999
      GOTO 5000 
C
C     -------------------      
C
C 
C 
C     Parse Roughness at Batiscan line  (parameter #251) 
C     
C     Parse the line:
 4251 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 87)) GOTO 5000
      EMESSNO = 9999
C 
C     -------------------      
C 
C 
C     Parse Tidal Wave line  (parameter #252) 
C     
C     Parse the line:
 4252 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 88)) GOTO 5000
      EMESSNO = 9999
C 
C     -------------------      
C 
C 
C     Parse ONTTOTALSUPDATA line  (parameter #253) 
C     
C     Parse the line:
 4253 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 89)) GOTO 5000
      EMESSNO = 9999
C 
C     -------------------      
C 
C 
C     Parse OTHER06ROUGHDATA line  (parameter #254) 
C     
C     Parse the line:
 4254 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 90)) GOTO 5000
      EMESSNO = 9999
C 
C     -------------------      
C 
C 
C     Parse OTHER06FLOWDATA line  (parameter #255) 
C     
C     Parse the line:
 4255 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 91)) GOTO 5000
      EMESSNO = 9999
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #256) 
C     
C     Not Implemented! 
 4256 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #257) 
C     
C     Not Implemented! 
 4257 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse USSICEWEEDADJ line  (parameter #258) 
 4258 IF(UseTSData(PLINE(FKC:LKC), LKC-FKC+1, 1064,
     &   AlreadySpecified(PNUM), 26)) GOTO 5000
      EMESSNO = 9999
      RETURN 
C 
C     -------------------      
C 
C 
C     Parse USSLIPEQUATION line  (parameter #259) 
C     
 4259 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1, 
     &   ': US Slip Equation to be determined using equation '//LLT,
     &   ' (1=old eqn, 2=2010)'//LLT,
     &   1, 2, 1184, USSlipEquation)) GOTO 5000
C
      EMESSNO = 9999
      RETURN
C     -------------------      
C 
C 
C     Parse SIDECHANNELMIN line  (parameter #260) 
C     
 4260 IF(ReadIntegerList(PLINE(FKC:LKC), LKC - FKC + 1, 12, 
     &   SideChanMinCapacity)) GOTO 5000  
C
      EMESSNO = 9999
      RETURN
C     -------------------      
C 
C 
C     Parse SIDECHANNELMAX line  (parameter #261) 
C     
 4261 IF(ReadIntegerList(PLINE(FKC:LKC), LKC - FKC + 1, 12, 
     &   SideChanMaxCapacity)) GOTO 5000  
C
      EMESSNO = 9999
      RETURN
C     -------------------      
C 
C 
C     Parse COMPWORKSOVERTOP line  (parameter #262) 
C     
 4262 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,7,CWOvertopElev))
     &   GOTO 5000
C
      EMESSNO = 9999
      RETURN
C     -------------------      
C
C     Parse COMPWORKSFREEBRD line  (parameter #263) 
C     
 4263 IF(ReadFloatList(PLINE(FKC:LKC), LKC - FKC + 1, 12, 
     &   CWMinimumFreeboard)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #264) 
C     
C     Not Implemented! 
 4264 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #265) 
C     
C     Not Implemented! 
 4265 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse LASTSTURGEONYEAR line  (parameter #266) 
C     
 4266 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 4,
     &   ': set year of last Sturgeon flow to '//LLT,
     &   ' '//LLT, 0, 99999, 1021, Last_JuneFlow)) 
     &      GOTO 5000 
      
      EMESSNO = 9999
      RETURN
C 
C     -------------------  
C 
C 
C     PDATA PARIDS(4337:4352) /'P2012SPHIGHLEVEL '/                           267
C     
 4267 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,0,P2012SPHIGHLEVEL))
     &  GOTO 5000  
 
      EMESSNO = 9999
      RETURN
C 
C     -------------------      
C 
C 
C     Parse ForcedSideChan line  (parameter #268) 
C     
 4268 IF(ReadIntegerList(PLINE(FKC:LKC), LKC - FKC + 1, 12, 
     &   ForcedSideChan)) GOTO 5000  
C
      EMESSNO = 9999
      RETURN
C 
C     -------------------      
C 
C 
C     Parse ForcedSupGates line  (parameter #269) 
C     
C     
 4269 IF(ReadFloatList(PLINE(FKC:LKC), LKC - FKC + 1, 12, 
     &   ForcedSupGates)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN 
C 
C     -------------------      
C 
C 
C     Parse IndexGatesSide line  (parameter #270) 
C     
C     Not Implemented! 
 4270 IF(ReadIntegerValue(PLINE(FKC:LKC), LKC - FKC + 1,
     &   ' '//LLT,
     &   ' '//LLT, -1, 5, 1021, IndexGatesSide)) 
     &      GOTO 5000 
      
      EMESSNO = 9999
      RETURN
C 
C     -------------------      
C 
C 
C     
C     Parse ???????? line  (parameter #271) 
C     
C     Not Implemented!      
 4271 GOTO 5000 
C 
C     -------------------  
C 
C 
C     Parse ???????? line  (parameter #272) 
C     
C     Not Implemented! 
 4272 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #273) 
C     
C     Not Implemented! 
 4273 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse P2012SPTARGET line  (parameter #274) 
C     DATA PARIDS(4353:4368) /'P2012SPTARGET '/                        274
 4274 IF(ReadFloatList(PLINE(FKC:LKC), LKC - FKC + 1, 12, 
     &   P2012SPTARGET)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C 
C     -------------------      
C 
C 
C     Parse P2012MHTARGET line  (parameter #275) 
C     DATA PARIDS(4353:4368) /'P2012MHTARGET '/                        275
 4275 IF(ReadFloatList(PLINE(FKC:LKC), LKC - FKC + 1, 12, 
     &   P2012MHTARGET)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #276) 
C     
C     Not Implemented! 
 4276 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #277) 
C     
C     Not Implemented! 
 4277 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #278) 
C     
C     Not Implemented! 
 4278 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #279) 
C     
C     Not Implemented! 
 4279 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #280) 
C     
C     Not Implemented! 
 4280 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #281) 
C     
C     Not Implemented! 
 4281 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #282) 
C     
C     Not Implemented! 
 4282 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #283) 
C     
C     Not Implemented! 
 4283 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #284) 
C     
C     Not Implemented! 
 4284 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #285) 
C     
C     Not Implemented! 
 4285 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #286) 
C     
C     Not Implemented! 
 4286 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #287) 
C     
C     Not Implemented! 
 4287 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #288) 
C     
C     Not Implemented! 
 4288 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #289) 
C     
C     Not Implemented! 
 4289 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #290) 
C     
C     Not Implemented! 
 4290 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #291) 
C     
C     Not Implemented! 
 4291 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #292) 
C     
C     Not Implemented! 
 4292 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #293) 
C     
C     Not Implemented! 
 4293 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #294) 
C     
C     Not Implemented! 
 4294 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse P2012PARA line  (parameter #295) 
C     
C     DATA PARIDS(4721:4736) /'P2012PARA '/
 4295 IF(ReadFloatList(PLINE(FKC:LKC), LKC - FKC + 1, 6, 
     &   P2012PARA)) GOTO 5000 
C
      EMESSNO = 9999
      RETURN 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #296) 
C     
C     Not Implemented! 
 4296 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #297) 
C     
C     Not Implemented! 
 4297 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #298) 
C     
C     Not Implemented! 
 4298 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse ???????? line  (parameter #299) 
C     
C     Not Implemented! 
 4299 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     
C     DATA PARIDS(4769:4784) /'P2012SPLOWLEVEL    '/                         300 
 4300 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,0,P2012SPLOWLEVEL))
     &  GOTO 5000  
 
      EMESSNO = 9999
      RETURN 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #301) 
C     
C     Not Implemented! 
 4301 GOTO 5000 
C 
C     ------------------- 

C 
C 
C     Parse ???????? line  (parameter #302) 
C     
C     Not Implemented! 
 4302 GOTO 5000 
C 
C     ------------------- 

C 
C 
C     Parse ???????? line  (parameter #303) 
C     
C     Not Implemented! 
 4303 GOTO 5000 
C 
C     ------------------- 

C 
C 
C     Parse ???????? line  (parameter #304) 
C     
C     Not Implemented! 
 4304 GOTO 5000 
C 
C     ------------------- 
C 
C 
C     Parse PARIDS(4769:4784) /'P2012MHTARGETADJ'/  (parameter #305) 
C     
 4305 IF(ReadValue(PLINE(FKC:LKC),LKC-FKC+1,1059,I,0,P2012MHTARGETADJ))
     &  GOTO 5000  
 
      EMESSNO = 9999
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #306) 
C     
C     Not Implemented! 
 4306 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #307) 
C     
C     Not Implemented! 
 4307 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #308) 
C     
C     Not Implemented! 
 4308 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #309) 
C     
C     Not Implemented! 
 4309 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #310) 
C     
C     Not Implemented! 
 4310 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #311) 
C     
C     Not Implemented! 
 4311 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #312) 
C     
C     Not Implemented! 
 4312 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #313) 
C     
C     Not Implemented! 
 4313 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #314) 
C     
C     Not Implemented! 
 4314 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #315) 
C     
C     Not Implemented! 
 4315 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #316) 
C     
C     Not Implemented! 
 4316 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #317) 
C     
C     Not Implemented! 
 4317 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #318) 
C     
C     Not Implemented! 
 4318 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #319) 
C     
C     Not Implemented! 
 4319 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #310) 
C     
C     Not Implemented! 
 4320 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #311) 
C     
C     Not Implemented! 
 4321 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #312) 
C     
C     Not Implemented! 
 4322 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #323) 
C     
C     Not Implemented! 
 4323 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #324) 
C     
C     Not Implemented! 
 4324 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #325) 
C     
C     Not Implemented! 
 4325 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #326) 
C     
C     Not Implemented! 
 4326 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #327) 
C     
C     Not Implemented! 
 4327 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #328) 
C     
C     Not Implemented! 
 4328 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #329) 
C     
C     Not Implemented! 
 4329 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #330) 
C     
C     Not Implemented! 
 4330 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #331) 
C     
C     Not Implemented! 
 4331 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #332) 
C     
C     Not Implemented! 
 4332 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #333) 
C     
C     Not Implemented! 
 4333 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #334) 
C     
C     Not Implemented! 
 4334 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #335) 
C     
C     Not Implemented! 
 4335 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #336) 
C     
C     Not Implemented! 
 4336 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #337) 
C     
C     Not Implemented! 
 4337 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #338) 
C     
C     Not Implemented! 
 4338 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #339) 
C     
C     Not Implemented! 
 4339 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #340) 
C  x  
C     Not Implemented! 
 4340 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #341) 
C     
C     Not Implemented! 
 4341 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #342) 
C     
C     Not Implemented! 
 4342 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #343) 
C     
C     Not Implemented! 
 4343 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #344) 
C     
C     Not Implemented! 
 4344 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #345) 
C     
C     Not Implemented! 
 4345 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #346) 
C     
C     Not Implemented! 
 4346 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #347) 
C     
C     Not Implemented! 
 4347 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #348) 
C     
C     Not Implemented! 
 4348 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #349) 
C     
C     Not Implemented! 
 4349 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #350) 
C     
C     Not Implemented! 
 4350 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #351) 
C     
C     Not Implemented! 
 4351 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #352) 
C     
C     Not Implemented! 
 4352 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #353) 
C     
C     Not Implemented! 
 4353 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #354) 
C     
C     Not Implemented! 
 4354 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #355) 
C     
C     Not Implemented! 
 4355 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #356) 
C     
C     Not Implemented! 
 4356 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #357) 
C     
C     Not Implemented! 
 4357 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #358) 
C     
C     Not Implemented! 
 4358 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #359) 
C     
C     Not Implemented! 
 4359 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #360) 
C     
C     Not Implemented! 
 4360 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #361) 
C     
C     Not Implemented! 
 4361 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #362) 
C     
C     Not Implemented! 
 4362 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #363) 
C     
C     Not Implemented! 
 4363 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #364) 
C     
C     Not Implemented! 
 4364 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #365) 
C     
C     Not Implemented! 
 4365 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #366) 
C     
C     Not Implemented! 
 4366 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #367) 
C     
C     Not Implemented! 
 4367 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #368) 
C     
C     Not Implemented! 
 4368 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #369) 
C     
C     Not Implemented! 
 4369 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #370) 
C     
C     Not Implemented! 
 4370 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #371) 
C     
C     Not Implemented! 
 4371 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #372) 
C     
C     Not Implemented! 
 4372 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #373) 
C     
C     Not Implemented! 
 4373 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #374) 
C     
C     Not Implemented! 
 4374 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #375) 
C     
C     Not Implemented! 
 4375 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #376) 
C     
C     Not Implemented! 
 4376 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #377) 
C     
C     Not Implemented! 
 4377 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #378) 
C     
C     Not Implemented! 
 4378 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #379) 
C     
C     Not Implemented! 
 4379 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #380) 
C     
C     Not Implemented! 
 4380 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #381) 
C     
C     Not Implemented! 
 4381 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #382) 
C     
C     Not Implemented! 
 4382 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #383) 
C     
C     Not Implemented! 
 4383 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #384) 
C     
C     Not Implemented! 
 4384 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #385) 
C     
C     Not Implemented! 
 4385 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #386) 
C     
C     Not Implemented! 
 4386 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #387) 
C     
C     Not Implemented! 
 4387 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #388) 
C     
C     Not Implemented! 
 4388 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #389) 
C     
C     Not Implemented! 
 4389 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #390) 
C     
C     Not Implemented! 
 4390 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #391) 
C     
C     Not Implemented! 
 4391 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #392) 
C     
C     Not Implemented! 
 4392 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #393) 
C     
C     Not Implemented! 
 4393 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #394) 
C     
C     Not Implemented! 
 4394 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #395) 
C     
C     Not Implemented! 
 4395 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #396) 
C     
C     Not Implemented! 
 4396 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #397) 
C     
C     Not Implemented! 
 4397 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #398) 
C     
C     Not Implemented! 
 4398 GOTO 5000 
C 
C     -------------------      
C 
C 
C     Parse P77AWinterMax line  (parameter #399) 
C     DATA PARIDS(4705:4720) /'P77AWINTERMAX   '/                        399
 4399 IF(ReadValue(PLINE(FKC:LKC), LKC-FKC+1, 1059, I, 2, TMPDBL)) THEN
         P77AWinterMax = IDNINT(TMPDBL)
         GOTO 5000
      ENDIF
      EMESSNO = 9999
      RETURN
C 
C     -------------------      
C 
C 
C     Parse ???????? line  (parameter #400) 
C     
C     Not Implemented! 
 4400 GOTO 5000 
C 
C     -------------------      
C
C
C     Jumped here because specified a keyphrase that had been
C     recognized in PARIDS, but not implemented here in the 
C     parsing module.  This is likely to occur due to some 
C     confusion of 'present' and 'future' tense in the 
C     documentation.  In many cases regarding Plan 1977A or
C     Ontario regulation, future input variables have been
C     identified and incorporated into PARIDS, but lacking
C     an immediate need, haven't been further processed.  
C     The various 'precip', 'evap', and 'runoff' data 
C     inputs likewise have spots reserved, but not used.
C
 4990 EMESSNO = 1165
      CALL ErrMess(EMESSNO,CHNKST,LOGCDE)
      WRITE(CTMP,4991) PARCODE, PNUM, LLT
      CALL LogIt(LOGCDE,0,1)
      CTMP(1:1) = LLT
 4991 FORMAT(' Input recognized but not implemented: ', A, 
     &   ' (parameter #',I3, ')', A1)
C      
C
C
C
C
 5000 NUMPARAM = NUMPARAM + 1
      IF( GenVerbosity .GE. 5 .AND. CTMP(1:1) .NE. LLT) THEN
         CALL LogIt(LOGD,0,1)
      ENDIF
C     CTMP(1:1) = LLT
C     CALL LogIt(LOGD,0,0)
C
 5012 FORMAT(I6)
 5015 FORMAT(' Got IOSTAT=',I4,' when accessing file ',A,'.',A1)
 5016 FORMAT(A4,A16,A,I1,A,A,A1,A1)
 5017 FORMAT(A4,A16,A,A,A1,A1)
 5018 FORMAT(A,A,A1)
 5019 FORMAT(A,A,A,A,A1)
 5020 FORMAT(A,A,A,A1)
 5021 FORMAT(A,I5,A,A1)
 5022 FORMAT(' IOSTAT=',I5,
     &   ', with build info and offending parameter line as:',A1)
 5023 FORMAT(A4,A16,A,L1,A1,A1)
 5024 FORMAT(' The invalid value is ',I5,A1)
 5025 FORMAT(' Value: ',I2,' (month ',I2,' in year ',I6,
     &   ' has ', I2,' days).',A1)
 5026 FORMAT(A4,A16,A,I3,A,A1)
 5027 FORMAT(' The invalid character is ===>',A1,'<===',A1)
 5028 FORMAT(A4,A16,A,A,A1,A1)
 5029 FORMAT(' The invalid string is ===>',A,'<===',A1)
 5030 FORMAT(A4,A16,': ', A, ' date set to {', I2, A3, I5, '}', A1)
C
      END
C
C      
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION FIRSTKC (STRING, STRLEN)
C
      IMPLICIT NONE
C
      INTEGER STRLEN, ASEQNO, I
      CHARACTER STRING*(*)
C      CHARACTER STRING
C      DIMENSION STRING (1,*)
C
C     Search from beginning of string looking for a regular "keyboard"  
C     character (ASCII char between decimal 33 and 126, not including 
C     the space character).
C
C     If the string contains no "keyboard" characters, 0 is returned.
C     DON'T USE THE RESULT OF THIS FUNCTION IN A STRING LIMIT UNLESS
C     YOU POSITIVELY KNOW THAT THE STRING CONTAINS A "KEYBOARD" CHARACTER!
C     For instance, SAMPLE(FIRSTKC(SAMPLE, 127):) will cause unexpected
C     results if SAMPLE is all blanks or something.  So safe use of this
C     for parsing lines read in from data files would be two-step, as in:
C          I= FIRSTKC(SAMPLE, 127)
C          IF(I .NE. 0)THEN
C              WRITE(90,2015)SAMPLE(I:LASTKC(SAMPLE,127-I))
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
   10    IF(ASEQNO .GE. 33 .AND. ASEQNO .LE. 126) GOTO 20

   20 IF (ASEQNO .LT. 33 .OR. ASEQNO .GT. 126) THEN
         I=0
      ENDIF
      FIRSTKC = I
      END
C
C-----------------------------------------------------------------------
C
      LOGICAL FUNCTION READQEQN(STRING, SLEN, K, YM, A, B, WT, ZC)
C
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER    STRING*(*)
C                  The character-array to parse.
C
      INTEGER      SLEN
C                  The length of the string.
C
      DOUBLE PRECISION K, YM, A, B, WT, ZC
C                  The variables to store the discharge equation
C                  coefficients in.  See descriptions in function
C                  OUTFLOW and subroutine CHECKPAR.
C                   
C      
C     Declarations regarding local variables or functions:
C
      DOUBLE PRECISION X
      INTEGER      I, J, M, FS, FE
C                  These are scratch variables.
C
      LOGICAL      FIELDBND, GETFLOAT
C                  Functions that helps parse comma-delimited input lines
C                  and floating-point input values.
C
C
C
C     Declarations regarding COMMON block variables:
C
C
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C                  See description in Block Data Genral
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                  See description in Block Data Genral
C
      CHARACTER    PARCODE*16, INDENT*4
C                  See descriptions in PARLINE().
C
      CHARACTER    CTMP*255, CWRK*80
C                  These are just buffers used throughout program.
C
      CHARACTER    LLT*1
C                  This is the character that marks the end of a line sent
C                  to the LogIt subroutine.
C
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /PARSEC/   INDENT, PARCODE
      COMMON /WKCHAR/   CTMP,LLT,CWRK
C
C
      READQEQN = .FALSE.
C
      J = 1
C
      DO 233, M=1,6
C        Get field bounds of expected value
         IF( .NOT. FIELDBND(STRING, SLEN, J, FS, FE, I ) ) THEN 
            CALL ErrMess(1052 + M, 'READQEQN'//LLT//'       ', LOGCDE)
            RETURN
         ENDIF 
C        Save the beginning of the next field, if any.
         J = I + 1
C
C        Read/convert the unitless value and store it in appropriate variable
         IF(GETFLOAT( STRING(FS:FE),FE-FS+1,0,X)) THEN
            GOTO (201, 202, 203, 204, 205, 206) M
  201          K = X
               GOTO 233
  202          YM = X
               GOTO 233
  203          A = X
               GOTO 233
  204          B = X
               GOTO 233
  205          WT = X
               GOTO 233
  206          ZC = X
         ELSE
            RETURN
         ENDIF
  233 CONTINUE
C
      IF(GenVerbosity .GE. 5) THEN
         WRITE(CTMP,'(A,A,6(A,F12.6),A1)')INDENT,PARCODE,
     &      ' K=',K,' YM=',YM,' A',A,' B=',B,' WT=',WT,' ZC=',ZC,LLT  
C        CALL LogIt(LOGD,0,0)
      ENDIF
C
      READQEQN = .TRUE.
C
      END
C
C-----------------------------------------------------------------------
C
      LOGICAL FUNCTION ReadValue(STRING, SLEN, ERRNUM, NEXTFLD, Units, 
     &   Value)
C
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER    STRING*(*)
C                  The character-array to parse.
C
      INTEGER      SLEN
C                  The length of the string.
C
      INTEGER      NEXTFLD
C                  The position of the next field delimiter plus one.
C
      INTEGER      ERRNUM
C                  The error number to pass to ErrMess() in case no value to read.
C
      INTEGER      Units
C                  The units that they value should be expressed returned in.
C
      DOUBLE PRECISION Value
C                  The variable to store the distance value in.
C                   
C      
C     Declarations regarding local variables or functions:
C
      INTEGER      FS, FE
C                  These are scratch variables.
C
      LOGICAL      FIELDBND, GETFLOAT
C                  Functions that helps parse comma-delimited input lines
C                  and floating-point input values.
C
C
C
C     Declarations regarding COMMON block variables:
C
C
      CHARACTER    PARCODE*16, INDENT*4
C                  See descriptions in PARLINE().
C
      CHARACTER    CTMP*255, CWRK*80
      CHARACTER    LLT*1
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      CHARACTER*6  UnitType(20)
      INTEGER      UnitTypeLen(20)
      INTEGER      BaseUnitType(20)
      DOUBLE PRECISION UnitConvFactor(20)
C                  See descriptions in Block Data Genral
C
C
C
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /PARSEC/   INDENT, PARCODE
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /UNITSC/   UnitType
      COMMON /UNITSI/   UnitTypeLen, BaseUnitType, UnitConvFactor
C
C
C      WRITE(CTMP, '(A,A,A,A1)')'  in=',STRING(1:SLEN),'=',LLT
C      CALL LogIt(LOGD,1,0)
      ReadValue = .FALSE.
C
C     Get field bounds of expected value
      IF( .NOT. FIELDBND(STRING, SLEN, 1, FS, FE, NEXTFLD) ) THEN 
         CALL ErrMess(ERRNUM, 'ReadValue'//LLT//'      ', LOGCDE)
         RETURN
      ENDIF 
C        Save the beginning of the next field, if any.
      NEXTFLD = NEXTFLD + 1
C
C     Read/convert the value and store it in Value
C     Units = 7 for m, 2 for m3s, 0 for no conversion
      IF(GETFLOAT( STRING(FS:FE),FE-FS+1,Units,Value)) THEN
         ReadValue = .TRUE.
         IF(GenVerbosity .GE. 5) THEN
            IF( Units .EQ. 0) THEN 
               CWRK(1:1) = LLT
            ELSE
               WRITE(CWRK, '(A,A1)') 
     &            UnitType(Units)(1:UnitTypeLen(Units)), LLT  
            ENDIF
            WRITE(CTMP,'(A,A,A,F12.6,1X,A,A1)') INDENT, PARCODE,
     &         ': read value of ', Value, CWRK(1:INDEX(CWRK, LLT)), LLT
C           CALL LogIt(LOGD,0,0)
         ENDIF
      ENDIF

C
      END
C
C----&------------------------------------------------------------------
C
      LOGICAL FUNCTION ReadFloatList(InLine, LineLength, 
     &   ValuesExpected, DoubleArray)
C
C
      IMPLICIT NONE
C
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER*(*) InLine
C                   (Input) The line to be parsed.
C
      INTEGER       LineLength
C                   (Input) The number of characters in InLine.
C
      INTEGER       ValuesExpected
C                   (Input) The number of floating point values expected to be read.
C
      DOUBLE PRECISION DoubleArray
      DIMENSION     DoubleArray(*)
C                   (Output) Floating point values read from InLine.
C
C
C
C     Declarations regarding local variables or functions:
C
      CHARACTER     CTMP2*256
      INTEGER       I, J, K, M, FS, FE
C                   These are scratch variables.
C
      LOGICAL       FIELDBND
C      LOGICAL       FIELDBND, GETFLOAT
C                   Functions that helps parse comma-delimited input 
C                   lines and floating-point input values.
C
C
C
C     Declarations regarding COMMON block variables:
C
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
      CHARACTER    PARCODE*16, INDENT*4
C                  See descriptions in PARLINE().
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                  See description in Block Data Genral
C
      CHARACTER    CTMP*255, CWRK*80
C                  These are just buffers used throughout program.
C
      CHARACTER    LLT*1
C                  This is the character that marks the end of a line 
C                  sent to the LogIt subroutine.
C
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /PARSEC/   INDENT, PARCODE
C
C
C
      ReadFloatList = .FALSE. 
      J = 1
C
C     Read each field 
      DO 233, M = 1, ValuesExpected
C        Get field bounds of expected value
         IF( .NOT. FIELDBND(InLine, LineLength, J, FS, FE, I ) ) THEN 
            CALL ErrMess(1137, 'ReadFloatList'//LLT, LOGCDE)
            WRITE(CTMP,'(A,I2,A,I2,A,A1)')' Expecting ', ValuesExpected,
     &         ' values on line, but only found ', M - 1, '!', LLT
            CALL LogIt(LOGCDE,0,0)
            RETURN
         ENDIF 
C        Save the beginning of the next field, if any.
         J = I + 1
C
         READ(InLine(FS:FE),'(F21.10)',ERR=521,IOSTAT=K) DoubleArray(M)
  233    CONTINUE
C
      IF(GenVerbosity .GE. 4) THEN
         DO 5143, K = 1, ValuesExpected
 5143       WRITE(CTMP2((K-1)*11+1:K*11),'(F10.4,A1)')DoubleArray(K),','
         WRITE(CTMP,'(A4,A16,A,A,A1,A1)') INDENT,PARCODE,
     &     ': Float list set to {',CTMP2(1:ValuesExpected*11-1),'}',LLT
C        CALL LogIt(LOGD,0,0)
      ENDIF
C
C
      ReadFloatList = .TRUE.
      RETURN
C
  521 CALL ErrMess(1138, 'ReadFloatList'//LLT, LOGCDE)
      WRITE(CTMP,'(A,I6,A,A,A1)') ' IOSTAT=',K,': Non-float field:', 
     &   InLine(FS:FE), LLT
      CALL LogIt(LOGCDE,1,0)
      END
C
C----&------------------------------------------------------------------
C
      LOGICAL FUNCTION UseTSData(STRING, SLEN, ERRNUM, SquawkIfAgain, 
     &   TSIndex)
C
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER    STRING*(*)
C                  (input) The character-array to parse.
C
      INTEGER      SLEN
C                  (input) The length of the string.
C
      INTEGER      ERRNUM
C                  (input) The error number to pass to ErrMess() in case file not found.
C
      LOGICAL      SquawkIfAgain
C                  (input) .TRUE. if don't expect more lines of this ID.
C                  .FALSE. if additonal values expected.
C
      INTEGER      TSIndex
C                  (input) The index to FileNames where the filename should be stored.
C
C
C     Declarations regarding local variables or functions:
C
C
      INTEGER      ENDSTR
C                  Function that returns end position of the current word.
C
      INTEGER      LastKC
C                  Function that returns the position of last 'keyboard char'
C 
C      LOGICAL      FIELDBND 
      INTEGER      I,J,K
      LOGICAL      LOGVAR
C                  Scratch variables
C
C
C
C     Declarations regarding COMMON block variables:
C
C
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C                  See description in Block Data Genral
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                  See description in Block Data Genral
C
      CHARACTER    PARCODE*16, INDENT*4
C                  See descriptions in PARLINE().
C
      CHARACTER    CTMP*255, CWRK*80
C                  These are just buffers used throughout program.
C
      CHARACTER    LLT*1
C                  This is the character that marks the end of a line sent
C                  to the LogIt subroutine.
C
      CHARACTER*256 FileNames(100)
C                  See description in Block Data Genral
C
      CHARACTER*256 RepeatLines(80)
C                  See description in Block Data Genral
C
      INTEGER      RepeatIndex
C                  Number of elements in RepeatLines
C 
      INTEGER      NPFILZ
C                  See description in Block Data Genral
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /PARSEC/   INDENT, PARCODE
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /FNAMEC/   FileNames, RepeatLines
      COMMON /FNAMEI/   RepeatIndex, NPFILZ
C
C
      UseTSData = .FALSE.
C
C
      IF( GenVerbosity .GE. 7 ) THEN
         WRITE(CTMP,'(2A,3(A,I4),A1)') 
     &      ' Running UseTSData with args String=', String, ', Slen=',
     &      Slen, ', ErrNum=', ErrNum, ' TSIndex=', TSIndex, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
c     IF( GenVerbosity .GE. 7 ) THEN
c        WRITE(CTMP,'(2A,I3,3(A,I4),A,L1,A1)') 
c    &      ' Running UseTSData with args String=', String,', Slen=',
c    &      Slen, ', ErrNum=', ErrNum, ' TSIndex=', TSIndex, 
c    &      ' SquawkIfAgain=', SquawkIfAgain, LLT
c        CALL LogIt(LOGD,0,0)
c     ENDIF
C
C     See if line refers to a file, or repeated data.  If repeating
C     data, should have at least one field delimiters (commas).
      IF( INDEX( String(1:Slen), ',' )  .GE. 1 ) THEN
C        There are commas.  So assume it's not a filename.  Process
C        line as repeating data.
C
C
C        Increment RepeatIndex, and store the data in RepeatLines:
         RepeatIndex = RepeatIndex + 1
C        XX TODO: make a nice error message for this:
         IF( RepeatIndex .GT. 80 ) STOP ' RepeatIndex >80'
         RepeatLines(RepeatIndex) = String(1:Slen)//LLT
         IF( GenVerbosity .GE. 3 ) THEN
            WRITE(CTMP,'(A,I2,A,A1)') 
     &         ' Assigned the following string as RepeatLines(',
     &         RepeatIndex, '):', LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,'(4X,A1,A,A1,A1)') '"', 
     &         RepeatLines(RepeatIndex)(1:Slen), '"', LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
C        Initialize K to 0.  This will be used later as number of  
C        related elements in RepeatLines.  
         K = 0
C
C        Write 3-digit (inc. leading spaces) number of element in array 
C        to unused portion of CTMP, for later use.
         WRITE(CTMP(250:252), '(I3)') RepeatIndex
C
C        See if already read a line for this input data.  
         IF(FileNames(TSIndex)(1:12) .NE. 'Unspecified'//LLT)THEN
C           Yep.  This is not first time parameter was supplied.
C           See if it is already defined as a filename, or 'Repeated'
            IF(FileNames(TSIndex)(1:15) .EQ. 
     &         'Repeating Data'//LLT)THEN
C              So this line gives additional repeating data.
C              Read number of previous lines read in, and their
C              entries in array RepeatData:
               READ(FileNames(TSIndex)(16:18), '(I3)') K
               CWRK = FileNames(TSIndex)(19:18+3*K)//CTMP(250:252)
C      WRITE(*, '(A,A,A)') '   recycled cwrk=', CWRK(1:55), '==55'
            ELSE
C              XXX TODO:message about re-defining from filename to repeating
               WRITE(*,'(A)') ' Hey! You tried to use both filename'//
     &            ' and repeating data for the same parameter!'
               STOP ' Stop that!'
            ENDIF
         ELSE
            CWRK = CTMP(250:252)
         ENDIF
C        Increment K to include the array element being added now
         K = K + 1        
C
C      WRITE(*, '(A,A,A)') '   new cwrk=', CWRK(1:55), '== up to 55'
C
C        Now write info to FileNames(TSIndex).  First 15 chars are
C        keyphrase 'Repeating Data'//LLT, which identifies data as not
C        coming from file, but repeating.  Next is the 3-digit number of 
C        elements in RepeatLines used for this data, followed by 3-digit
C        list of the indexes to the appropriate elements in RepeatLines.
C
         WRITE(FileNames(TSIndex), '(A,A1,I3,A,A1)') 
     &      'Repeating Data', LLT, K, CWRK(1:K*3), LLT
         IF( GenVerbosity .GE. 3 ) THEN
            WRITE(CTMP,5027) INDENT,PARCODE,
     &          ': set data filename number ',TSIndex,' to {',
     &          FileNames(TSIndex)(1:14)//
     &          FileNames(TSIndex)(16:(K+1)*3+15),'}',LLT
         ENDIF
         UseTSData = .TRUE.
C        Set the AlreadyRead flag to .FALSE. because additional lines 
C        might come.  This is sorta loose, because if all the values 
C        were specified on this record, a subsequent record with same 
C        parameter and using repeating data wouldn't be caught until
C        repeating data gets processed.  But this should handle most
C        user errors.
         SquawkIfAgain = .FALSE.
         RETURN
      ENDIF
C   
C     If made it here, then line is presumed to contain a filename.
C
C      WRITE(*,'(A)') '  UseTSData 04'
C
C     Determine end of filename
      I=LastKC( STRING, SLEN )
C
C      WRITE(*,'(A)') '  UseTSData 05'
C
C     Calling routine guarantees that I > 0, so no need to check that 
C     condition.
      IF( I .GT. 255 )THEN
C        Not enough room in FileNames element for STRING + LLT
C        TODO:  issue message via LogIt and ErrMess
         STOP 'fname too long'
      ELSE IF (I .EQ. 4) THEN
C        If four characters, could be 'null', so convert it to 
C        uppercase now before testing it later.
         CTMP = '     '
         DO 40, J = 1, 4
            I = ICHAR( STRING(J:J) )
            IF( I .GE. 97 .AND. I .LE. 122) THEN
               CTMP(J:J) = CHAR(I-32)
            ELSE
               CTMP(J:J) = STRING(J:J)
            ENDIF
 40      CONTINUE
         I=4
      ENDIF
C           
      IF ( I .EQ. 4 .AND. CTMP(1:4) .EQ. 'NULL') THEN
         STRING(1:4) = 'NULL'
      ELSE
C        See if the file exists:
C        WRITE(*,'(A,A,A1)')' Looking for file ',ParameterFile(BFN:EFN),'!'
         INQUIRE(FILE=STRING(1:I),EXIST=LOGVAR)
         IF(.NOT. LOGVAR) THEN
            CALL ErrMess(ERRNUM,'UseTSData'//LLT//'      ',LOGCDE)
            WRITE(CTMP,5020)
     &         ' Filename "',STRING(1:I),'" read from line:',LLT
            CALL LogIt(LOGCDE,0,0)
            RETURN
         ENDIF
      ENDIF
C     
      IF(FileNames(TSIndex)(1:12) .NE. 'Unspecified'//LLT)THEN 
         CALL ErrMess(1173,'UseTSData'//LLT//'      ',LOGCDE)
         WRITE(CTMP,'(6A)')' Re-defining value of ',
     &      PARCODE(1:EndStr(PARCODE,16)), ' from "',
     &      FileNames(TSIndex)(1:INDEX(FileNames(TSIndex),LLT) - 1), 
     &      '"', LLT
         CALL LogIt(LOGCDE,0,0)
         WRITE(CTMP,'(4A)') ' to file "', String(1:I), '"', LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
C     Ok, file is there (we'll assume that it's readable and intact) or
C     else that it's set to "NULL".
C     Overwrite current contents of FileNames(TSIndex) with put in new name:
      WRITE(FileNames(TSIndex), '(A,A1)') STRING(1:I), LLT
      IF( GenVerbosity .GE. 3 ) THEN
         WRITE(CTMP,5027) INDENT,PARCODE,
     &      ': set data filename number ',TSIndex,' to {',
     &      FileNames(TSIndex)(1:I),'}',LLT
      ENDIF
C
 5027 FORMAT(A4,A16,A,I2,A,A,A1,A1)
 5020 FORMAT(A,A,A,A1)
C
      SquawkIfAgain = .TRUE.
      UseTSData = .TRUE.
      END
C
C
C----&------------------------------------------------------------------
C 
      LOGICAL FUNCTION WhichOutputFiles( LakeNum, String, Slen )
C
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER    String*(*)
C                  The character-array to parse.
C
      INTEGER      Slen
C                  The length of the string.
C
      INTEGER      LakeNum
C                  Refers to the lake in question, corresponding to
C                  the elements in LakesToSolve.
C
C
C     Declarations regarding local variables or functions:
C
C
C
      INTEGER      I,J
C                  Scratch variables
C
      CHARACTER    Intervals*4
C                  Holds letters 'DWQM' for testing against String
C
      CHARACTER*16  MakeIntervalName
C                   Returns a string with the interval name.  For output.
C
C 
C     Declarations regarding COMMON block variables:
C
      LOGICAL      WriteTSOutput(4,2,8)
C                  For each lake (designated by third index - 8), for
C                  each type of time-series output (flow or ave level),
C                  for each interval (d, w, q, or m), .TRUE. if should
C                  write out the data to a time-series output file. 
C       
      LOGICAL      WriteErieOutflow
C                  See description in Block Data Genral
C
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
      CHARACTER*21 LAKENAME(16), RiverName(16)
      INTEGER      NLN(16), LkToSlv2LN(8), NRN(16), LkSlv2RtLoopI(8)
      CHARACTER    CTMP*255, CWRK*80, LLT*1
C                  See descriptions in Block Data Genral
C
      CHARACTER    PARCODE*16, INDENT*4
C                  See descriptions in PARLINE().
C
C
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /PARSEC/   INDENT, PARCODE
      COMMON /LNAMEC/   LAKENAME, RiverName
      COMMON /LNAMEI/   NLN, LkToSlv2LN, NRN, LkSlv2RtLoopI
      COMMON /OUTSTF/   WriteTSOutput, WriteErieOutflow
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
C
      DATA Intervals /'DWQM'/
C
C     Make sure string is uppercase
      DO 40, J = 1, Slen
         I = ICHAR( String(J:J) )
         IF( I .GE. 97 .AND. I .LE. 122) THEN
            CWRK(J:J) = CHAR(I-32)
         ELSE
            CWRK(J:J) = String(J:J)
         ENDIF
 40      CONTINUE
C
C
C     Turn on flag for each output interval, as specified
      DO 50, I = 1, 4
C
C        Get intervalname for use in logging.  Should be lotsa
C        room at end of CWRK, since String should only needs to
C        be four characters long.  
         CWRK(60:) = MakeIntervalName( I, .FALSE., J)
C
         IF( INDEX( CWRK(1:Slen), Intervals(I:I) )  .GE. 1 ) THEN
C           Found the letter we're looking for.  Turn on flags for
C           both levels and flows.
C
            IF( GenVerbosity .GE. 3 ) THEN
               WRITE(CTMP,5027) INDENT,PARCODE,  CWRK(60:60+J-1), 
     &            LakeName(LkToSlv2LN(LakeNum))
     &            (1:NLN(LkToSlv2LN(LakeNum))), LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
C          
            DO 60, J = 1, 2
 60            WriteTSOutput(I, J, LakeNum) = .TRUE.
C
C
C         ELSE
C           Didn't find the letter, so leave flag at it's default
C           setting of off.
C
         ENDIF
 50      CONTINUE
C
C     Kinda superfluous.  Never gets set to .false.
      WhichOutputFiles = .TRUE.
C
C     Avoid repeating message already sent
      CTMP(1:1)=LLT
C
 5027 FORMAT( A4, A16, ': requesting output files with mean ', A, 
     &   ' levels and flows for ', A, A1)
C
      END
C
C----&------------------------------------------------------------------
C
      LOGICAL FUNCTION ReadIntegerList(InLine, LineLength, 
     &   ValuesExpected, IntArray)
C
C
      IMPLICIT NONE
C
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER*(*) InLine
C                   (Input) The line to be parsed.
C
      INTEGER       LineLength
C                   (Input) The number of characters in InLine.
C
      INTEGER       ValuesExpected
C                   (Input) The number of integers expected to be read.
C
      INTEGER       IntArray
      DIMENSION     IntArray(*)
C                   (Output) Integers read from InLine.
C
C
C
C     Declarations regarding local variables or functions:
C
      CHARACTER     CTMP2*256
      INTEGER       I, J, K, M, FS, FE
C                   These are scratch variables.
C
      LOGICAL       FIELDBND
C      LOGICAL       FIELDBND, GETFLOAT
C                   Functions that helps parse comma-delimited input 
C                   lines and floating-point input values.
C
C
C
C     Declarations regarding COMMON block variables:
C
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
      CHARACTER    PARCODE*16, INDENT*4
C                  See descriptions in PARLINE().
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                  See description in Block Data Genral
C
      CHARACTER    CTMP*255, CWRK*80
C                  These are just buffers used throughout program.
C
      CHARACTER    LLT*1
C                  This is the character that marks the end of a line 
C                  sent to the LogIt subroutine.
C
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /PARSEC/   INDENT, PARCODE
C
C
C
      ReadIntegerList = .FALSE. 
      J = 1
C
C     Read each field 
      DO 233, M = 1, ValuesExpected
C        Get field bounds of expected value
         IF( .NOT. FIELDBND(InLine, LineLength, J, FS, FE, I ) ) THEN 
            CALL ErrMess(1137, 'ReadIntegerList'//LLT, LOGCDE)
            WRITE(CTMP,'(A,I2,A,I2,A,A1)')' Expecting ', ValuesExpected,
     &         ' values on line, but only found ', M - 1, '!', LLT
            CALL LogIt(LOGCDE,0,0)
            RETURN
         ENDIF 
C        Save the beginning of the next field, if any.
         J = I + 1
C
         READ( InLine(FS:FE), '(I9)', ERR=521, IOSTAT=K) IntArray(M)
  233    CONTINUE
C
      IF(GenVerbosity .GE. 5) THEN
         DO 5143, K = 1, ValuesExpected
 5143       WRITE(CTMP2((K-1)*7+1:K*7), '(I6,A1)') IntArray(K), ',' 
         WRITE(CTMP,'(A4,A16,A,A,A1,A1)') INDENT,PARCODE,
     &     ': Integer list set to {',CTMP2(1:ValuesExpected*7-1),'}',LLT
C        CALL LogIt(LOGD,0,0)
      ENDIF
C
C
      ReadIntegerList = .TRUE.
      RETURN
C
  521 CALL ErrMess(1138, 'ReadIntegerList'//LLT, LOGCDE)
      WRITE(CTMP,'(A,I6,A,A,A1)') ' IOSTAT=',K,': Non-integer field:', 
     &   InLine(FS:FE), LLT
      CALL LogIt(LOGCDE,1,0)
      END
C
C----&------------------------------------------------------------------
C
      LOGICAL FUNCTION ReadIntegerValue(InLine, LineLength, TextFrag1,
     &   TextFrag2, Min, Max, ErrNum, Value)
C
C
      IMPLICIT NONE
C
C
C     Declarations regarding subroutine arguments:
C      
      CHARACTER*(*) InLine
C                   (Input) The line to be parsed.
C
      INTEGER       LineLength
C                   (Input) The number of characters in InLine.
C
      CHARACTER*(*) TextFrag1, TextFrag2
C                   (Input) Pieces of a message to write to detailed
C                   log if number successfully read.  Must be 
C                   terminated by LLT.
C
      INTEGER       Min, Max, ErrNum
C                   (Input) The minimum and maximum values permitted,
C                   and the error message number for values out of
C                   this range.
C
      INTEGER       Value
C                   (Output) Integer read from InLine.
C
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER       J, FS, FE
C                   These are scratch variables.
C
      LOGICAL       FIELDBND
C                   Function that helps parse comma-delimited input 
C
      INTEGER       ENDSTR
C                   Function that returns end position of the current word.
C
C
C     Declarations regarding COMMON block variables:
C
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
      CHARACTER    PARCODE*16, INDENT*4
C                  See descriptions in PARLINE().
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                  See description in Block Data Genral
C
      CHARACTER    CTMP*255, CWRK*80
C                  These are just buffers used throughout program.
C
      CHARACTER    LLT*1
C                  This is the character that marks the end of a line 
C                  sent to the LogIt subroutine.
C
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /PARSEC/   INDENT, PARCODE
C
C
C
      ReadIntegerValue = .FALSE. 
C
C     WRITE(*,'(3A)') ' In ReadIntegerValue1:', InLine(1:LineLength), '#'
C
C     Get field bounds of expected value
      IF( .NOT. FIELDBND(InLine, LineLength, 1, FS, FE, J ) ) THEN 
         CALL ErrMess(1059, 'ReadIntegerValu'//LLT, LOGCDE)
         RETURN
      ENDIF 
C
C      WRITE(*,'(3A)') ' In ReadIntegerValue2:', InLine(1:LineLength), '#'
C
      READ( InLine(FS:ENDSTR( InLine(FS:FE), FE - FS + 1)), '(I9)', 
     &   ERR=521, IOSTAT=J) Value
C
  444 FORMAT( I24)
C
      WRITE(CWRK, 444) Value
      DO 455, J = 1, 24
  455    IF( CWRK(J:J) .NE. ' ') GOTO 454
C 454 J = J - 1
C
C      WRITE(*,'(3A)') ' In ReadIntegerValue3:', InLine(1:LineLength), '#'
C
  454 IF( Value .NE. -99 .AND. 
     &   (Value .LT. Min .OR. Value .GT. Max)) THEN
         CALL ErrMess(ErrNum,'ReadIntegerValu'//LLT,LOGCDE)
         WRITE(CWRK(55:), 444) Min
         DO 456, FS = 55, 78
  456       IF( CWRK(FS:FS) .NE. ' ') GOTO 457
  457    FS = 78 - FS + 1
         CWRK(25:25 + FS + 4) = CWRK(78 - FS:78)//' to '
C
         WRITE(CWRK(55:), 444) Max
         DO 458, FE = 55, 78
  458       IF( CWRK(FE:FE) .NE. ' ') GOTO 459
  459    FE = 78 - FE + 1
         CWRK(25 + FS + 4:25 + FS + 4 + FE) = CWRK(78 - FE:78)
C        WRITE(*, '(2(6X,A,I3), /, A, A, A)') 'FS=', FS, 'FE=', FE,
C    &        '  CWRK=', CWRK(1:70), '###'
 5024    FORMAT(' Value read: ',A, '; allowable range: ',A,A1)
         WRITE(CTMP,5024) CWRK(J:24), 
     &      CWRK(25:25 + FS + 4 + FE), LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
C      WRITE(*,'(3A)') ' In ReadIntegerValue4:',InLine(1:LineLength),'#'
C
      IF(GenVerbosity .GE. 5) THEN
         WRITE(CTMP,'(A4,A16,A,A1,A,A1,A,A1)') INDENT,PARCODE,
     &      TextFrag1(1:INDEX(TextFrag1,LLT) - 1), '{', CWRK(J:24),'}', 
     &      TextFrag2(1:INDEX(TextFrag2,LLT) - 1), LLT   
C        CALL LogIt(LOGD,0,0)
      ENDIF
C
C
      ReadIntegerValue = .TRUE.
      RETURN
C
  521 CALL ErrMess(1140, 'ReadIntegerValu'//LLT, LOGCDE)
      WRITE(CTMP,'(A,I6,A,A,A1)') ' IOSTAT=', J,
     &   ': Non-integer characters:', InLine(FS:FE), LLT
      CALL LogIt(LOGCDE,1,0)
      END
