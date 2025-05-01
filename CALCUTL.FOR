C
C-----------------------------------------------------------------------
C
C
      DOUBLE PRECISION FUNCTION Round( Value, PowerOfTen )
C 
C     Returns Value rounded to the PowerOfTen position.
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      DOUBLE PRECISION Value
C                  The 64-bit number to be rounded.
C
      INTEGER      PowerOfTen
C                  The digit which Value is to be rounded to.  3 means
C                  round the thousand, 2 rounds to hundredth, 1 to
C                  the tens, 0 to the ones position, -1 is to the tenth,
C                  -2 to the hundredth, -3 to the thousandth, -4 to 
C                  the ten-thousandth, -5 to the hundred-thousandth, and
C                  -6 to the millionth.  No other numbers are allowed, 
C                  since it refers to an element of array ExpFactor.
C
C     Declarations regarding local variables or functions:
C
      DOUBLE PRECISION tmp
C                  scratch variable
C
      DOUBLE PRECISION Constant(0:1)
C                  The number to add for rounding purposes.  Even
C                  numbers add first element, odd numbers add 2nd.
C
      DOUBLE PRECISION ExpFactor(-6:3)
C                  This is the array containing the position which is
C                  to be rounded to.  Array index set up to correspond
C                  to values of PowerOfTen (ie, 10E-3 to 10E+3).
C
C     Declarations regarding COMMON block variables:
C
C
C     
      DATA Constant / 0.4999999D0, 0.5000001D0/
C
C     These values seem in reverse order, because we multiply by them
      DATA ExpFactor /1000000.D0, 100000.D0, 10000.D0, 1000.D0, 100.D0,
     &    10.D0, 1.D0, 0.1D0, 0.01D0, 0.001D0/
C
C
      IF( PowerOfTen .GE. -6 .AND. PowerOfTen .LE. 3) THEN
C
C        Use Lake Survey rounding technique 
         tmp = DABS(Value * ExpFactor(PowerOfTen))
C
         Round = DINT( Constant(  MOD( IDINT(tmp),2 ) ) + tmp) /
     &      DSIGN(ExpFactor(PowerOfTen), Value)
C
      ELSE IF ( PowerOfTen .EQ. -99 ) THEN
C        Do no rounding 
         Round = Value
      ELSE
         STOP 'Round out of range'
      ENDIF
      END
C      WRITE(*, 161) 567.344999D0, Round(567.344999D0, -2)
C      WRITE(*, 161) 567.3449999D0, Round(567.3449999D0, -2)
C      WRITE(*, 161) 567.34499999D0, Round(567.34499999D0, -2)
C      WRITE(*, 161) 567.344999999D0, Round(567.344999999D0, -2)
C      WRITE(*, 161) 567.3449999999D0, Round(567.3449999999D0, -2)
C      WRITE(*, 161) 567.34499999999D0, Round(567.34499999999D0, -2)
C      WRITE(*, 161) 567.34500000000D0, Round(567.34500000000D0, -2)
C      WRITE(*, 161) 567.34500000001D0, Round(567.34500000001D0, -2)
C      WRITE(*, 161) 567.3450000001D0, Round(567.3450000001D0, -2)
C      WRITE(*, 161) 567.345000001D0, Round(567.345000001D0, -2)
C      WRITE(*, 161) 567.34500001D0, Round(567.34500001D0, -2)
C      WRITE(*, 161) 567.3450001D0, Round(567.3450001D0, -2)
C      WRITE(*, 161) 567.345001D0, Round(567.345001D0, -2)
C      WRITE(*, 161) 567.354999D0, Round(567.354999D0, -2)
C      WRITE(*, 161) 567.3549999D0, Round(567.3549999D0, -2)
C      WRITE(*, 161) 567.35499999D0, Round(567.35499999D0, -2)
C      WRITE(*, 161) 567.354999999D0, Round(567.354999999D0, -2)
C      WRITE(*, 161) 567.3549999999D0, Round(567.3549999999D0, -2)
C      WRITE(*, 161) 567.35499999999D0, Round(567.35499999999D0, -2)
C      WRITE(*, 161) 567.35500000000D0, Round(567.35500000000D0, -2)
C      WRITE(*, 161) 567.35500000001D0, Round(567.35500000001D0, -2)
C      WRITE(*, 161) 567.3550000001D0, Round(567.3550000001D0, -2)
C      WRITE(*, 161) 567.355000001D0, Round(567.355000001D0, -2)
C      WRITE(*, 161) 567.35500001D0, Round(567.35500001D0, -2)
C      WRITE(*, 161) 567.3550001D0, Round(567.3550001D0, -2)
C      WRITE(*, 161) 567.355001D0, Round(567.355001D0, -2)
C  161 FORMAT(' Double ', F18.12, ' rounds to ', F18.3)
C      STOP
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION GETAREA( LAKENUM, Z )
C
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      DOUBLE PRECISION Z
C                  The water surface elevation of the lake.
C
      INTEGER      LAKENUM
C                  The lake id number to be used as the first element of
C                  array AVSZ.
C
C     Declarations regarding local variables or functions:
C
      INTEGER      I, J
      DOUBLE PRECISION TMPX1, TMPX2
      DOUBLE PRECISION POLATE
C                  scratch variable
C
C     Declarations regarding COMMON block variables:
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
      INTEGER      NUMAVSZ(8)
C                  See descriptions in Block Data Genral
C
      DOUBLE PRECISION AVSZ(8,5,2)
C                  See descriptions in Block Data Genral
C
C
      COMMON /RCAVSZ/   NUMAVSZ, AVSZ
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
C
      IF ( NUMAVSZ(LAKENUM) .LE. 1 ) THEN
C        constant area:
         GETAREA=AVSZ(LAKENUM,1,2) 
      ELSE
C
C        inter/extra polate from curve.  find closest X values to use:
         J = 2
         TMPX1 = DABS((Z - AVSZ(LAKENUM,1,1))+(Z - AVSZ(LAKENUM,2,1))) 
         DO 311, I = 3, NUMAVSZ(LAKENUM)
            TMPX2=DABS((Z-AVSZ(LAKENUM,I-1,1))+(Z-AVSZ(LAKENUM,I,1)))
            IF( TMPX2 .LT. TMPX1 ) THEN
               J = I
               TMPX1 = TMPX2
            ENDIF
  311       CONTINUE
C
C        Figure out which message level to use
         IF( LakeNum .EQ. 1 ) THEN
            I = SupVerbosity
         ELSE IF ( LakeNum .EQ. 5) THEN
            I = OntVerbosity
         ELSE 
            I = MLVerbosity
         ENDIF
         GETAREA=POLATE(AVSZ(LAKENUM,J-1,1),AVSZ(LAKENUM,J-1,2),
     &      AVSZ(LAKENUM,J,1), AVSZ(LAKENUM,J,2), Z, I)
      ENDIF
      END        
C
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION POLATE( X1, Y1, X2, Y2, X, Verbosity )
C
C     This funtion returns the value corresponding to X that is inter- 
C     or extrapolated from the line segment defined by (X1,Y1) and 
C     (X2, Y2).  It tries to catch undefined slopes and avoid division
C     by zero, but relies on the calling procedure to provide the 
C     messaging information that describes how the situation came about.
C
C
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      DOUBLE PRECISION X1, Y1, X2, Y2
C                  Floating point numbers defining the line to 
C                  interpolate within or extrapolate from.
C
      DOUBLE PRECISION X
C                  The floating point ordinate that we are
C                  interpolating an abscissa for.
C
      INTEGER      Verbosity
C
C     Declarations regarding local variables or functions:
C
      INTEGER      I, J
      DOUBLE PRECISION ExtrapDist
      LOGICAL      LOGUSE(8)
C                  scratch variables
C
      CHARACTER    CHNKST*16
C                  This is the calling subroutine name passed to
C                  ErrMess() in case of error messages.
C
C     Declarations regarding COMMON block variables:
C
      INTEGER      Extrapolations
      DOUBLE PRECISION LargestExtrapolation
C                  Variables used to describe inadequate area-elev data
C
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C                  See description in Block Data Genral
C
      CHARACTER    CTMP*255,CWRK*80
C                  These are just buffers used throughout program.
C
      CHARACTER    LLT*1
C                  This is the character that marks the end of a line 
C                  sent to the LogIt subroutine.
C
C
      COMMON /EXTRAP/   LargestExtrapolation, Extrapolations 
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
C
C
      CHNKST = 'POLATE'//LLT
C
      I = 0
C     First make sure we don't divide by zero.  Should get 15 
C     significant digits in a DOUBLE, so test at that level:
C
      IF( (X2 - X1)*1.D14 .LT. 1.D0) THEN
         CALL ErrMess(1051,CHNKST,LOGCDE)
         I = 1051
      ELSE
C        Should be able to compute a result, so do it:
         POLATE = Y1 + (X-X1) * (Y2-Y1)/(X2-X1)
C
C        If extrapolating, squawk about it, but do it anyway.
         IF((X .LT. X1 .AND. X .LT. X2) .OR. 
     &      (X .GT. X2 .AND. X .GT. X1)) THEN
            Extrapolations = Extrapolations + 1
            IF(X .LT. X1 .AND. X .LT. X2) THEN
               IF ( X2 .LT. X1 ) THEN 
                  ExtrapDist = X2 - X
               ELSE
                  ExtrapDist = X1 - X
               ENDIF
            ELSE 
               IF ( X2 .GT. X1 ) THEN 
                  ExtrapDist = X - X2
               ELSE
                  ExtrapDist = X - X1
               ENDIF
            ENDIF
            IF ( ExtrapDist > LargestExtrapolation ) THEN
               LargestExtrapolation = ExtrapDist
            ENDIF
            CALL ErrMess(1052,CHNKST,LOGD)
            I = 1052
         ENDIF
      ENDIF 


C
      IF(I .GE. 1051 .OR. Verbosity .GE. 9) THEN
C        Send messages to console if fatal error, else just send to
C        detailed log.
         IF(I .EQ. 1051) THEN
            DO 6255, J=1,8
 6255          LOGUSE(J) = LOGCDE(J)
         ELSE
            DO 6256, J=1,8
 6256          LOGUSE(J) = LOGD(J)
         ENDIF
C
         WRITE(CTMP,5) ' Given X1,Y1=', X1, ',',Y1, LLT
         CALL LogIt(LOGUSE,1,0)
         WRITE(CTMP,5) '   and X2,Y2=', X2, ',',Y2, LLT
         CALL LogIt(LOGUSE,0,0)
         WRITE(CTMP,5) '  then  X,Y =', X, ',', POLATE, LLT
         CALL LogIt(LOGUSE,0,0)
         WRITE(CTMP,'(A,I19,A1)') '  I =', I, LLT
         CALL LogIt(LOGUSE,0,0)
         IF(I .EQ. 1051)THEN
            CALL SAVETOD
            CALL ADIOS
         ENDIF
      ENDIF
C
    5 FORMAT(A,D27.17,A,D27.17,A1)
      END
C
C

C----------------------------------------------------------------------
C
C
      DOUBLE PRECISION FUNCTION OUTFLOW(ZU,ZD,K,YM,A,B,WT,ZC)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION ZU,ZD,ZC,K,YM,WT,A,B
C      DOUBLE PRECISION WTDEPTH, FALL
C
C  Given constants k, ym, a, and b, and variables zu and zd, 
C  calculate channel outflow.  Terms include:
C         zu = upstream lake elevation
C         zd = downstream lake elevation
C         ZC = fall constant
C         k = outflow equation coefficient
C         ym = outlet bottom elevation
C         a = depth exponent
C         b = fall exponent
C         wt = relative weight between upstream and downstream lake levels
C  The constants k, ym, a, and b are dictated by the stage-discharge 
C  relationship for each connecting channel.  This function was written by
C  Anne H. Clites, NOAA/GLERL, September, 1995.  The ZC variable was added
C  by Matt in order to accommodate the metric forms of the equations.
C
C
      INTEGER      Verbosity
C     Declarations regarding COMMON block variables:
C
C
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C                   See descriptions in Block Data Genral
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
      INTEGER      TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &             MissingVal, SecsInUniformMon

      DOUBLE PRECISION MPFT, HALF, ONE, ZERO
C                   See descriptions in Block Data General.  
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
C
C
C
C     Figure out which message level to use (that of US body)
      IF( ZU .GT. 180 ) THEN
         Verbosity = SupVerbosity
      ELSE IF( ZU .GT. 100 ) THEN
         Verbosity = MLVerbosity
      ELSE 
         Verbosity = OntVerbosity
      ENDIF
C
C      INTEGER COUNT
C
C      DATA COUNT /0/
C      COUNT = COUNT + 1
C
C     Error checking.
C
C     TODO:  read in arange for which the equation is valid, test 
C            and issue warning if outside calibrated range.
      IF (ZD .GE. (ZU+ZC) .OR. YM .GT. ZU) THEN
         CALL ErrMess(1065,'OUTFLOW'//LLT//'        ',LOGCDE)
         WRITE(CTMP,5387) ZU, YM, ZD, K, LLT
 5387    FORMAT(' Upstream level:',F8.4,' Channel invert:',F8.4,
     &      '  Downstream level:',F8.4,'  K:', F12.7, A1)
         CALL LogIt(LOGCDE,1,0)
         CALL SAVETOD
         CALL ADIOS
      ENDIF
C
C
      OUTFLOW= K* (WT*ZU+(ONE-WT)*ZD-YM)**A * (ZU-ZD+ZC)**B
C
C
      IF( Verbosity .GE. 8) THEN                          
         WRITE(CTMP,5444) K, WT, YM, ZC, A, B, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,5445) ZU, ZD, OUTFLOW, (WT*ZU+(ONE-WT)*ZD-YM)**A,LLT
         CALL LogIt(LOGD,0,0)
C        WRITE(CTMP,5445) ZC, ONE, WT*ZU + (ONE-WT)*ZD - YM, ZU-ZD+ZC, 
C     &     OUTFLOW, LLT
 5444   FORMAT(5X, 'Outflow Eqn: ', ' K=', F12.7, ' WT=', F5.3, ' YM=', 
     &     F12.7, ' ZC=', F5.3,' A=', F15.12, ' B=', F15.12, A1)
 5445   FORMAT(5X, 'Input/Output: ', ' ZU=',F15.11,' ZD=', F21.11, 
     &     ' Outflow=', F15.8, ' Depth term=', F15.11,A1)
      ENDIF
C
C      WTDEPTH = WT*ZU + (ONE-WT)*ZD - YM
C      FALL = ZU-ZD+ZC
C      OUTFLOW = K*(WTDEPTH)**A*(FALL)**B
C
      END
C
C
C----&------------------------------------------------------------------
C
C
      SUBROUTINE ConChanConvert( USLoop, USFlow, AccumulatedFlow,
     &    DaysInAccum, WeightedFlow)
C
C     Call this routine after computing a connecting channel flow.  Does
C     the arithmetic for converting the flow to the routing interval
C     used in the next lake (routing loop).  This routine handles all
C     possible translations completely, except the case of monthly flows
C     US heading into weekly routing downstream.  In this case, the
C     routine calculates the weighted weekly flow for the pending DS
C     routing interval in WeightedFlow.  Subsequent weeks wihin the 
C     month must draw upon the monthly average US flow (FWRK(22,?)).
C
C
      IMPLICIT NONE 
C
C     Declarations regarding function arguments:
C
      INTEGER       USLoop
C                   (Input) Loop number for 1 for Lake Superior, 2 for 
C                   middle lakes, 3 for Lake Ontario.  
C
      DOUBLE PRECISION USFlow
C                   (Input) The outflow coming from the US loop.
C
      DOUBLE PRECISION WeightedFlow
C                   (Output) The outflow converted into the routing
C                   period used in the DS loop.
C
      DOUBLE PRECISION AccumulatedFlow
C                   (Input/Output) Accumulation of outflow since last
C                   time WeightedFlow was completed, weighted according
C                   to the number of days spent at each flowrate.
C
      INTEGER       DaysInAccum
C                   (Input/Output) Number of days accumulated since
C                   last time WeightedFlow was completed.
C
C     Declarations regarding local variables or functions:
C
      CHARACTER*18  NameOfRiver(3)      
      INTEGER       NameOfRiverNum(3)
C
      INTEGER       DaysTooFar
C                   Number of days that the current outflow from
C                   USLoop extends beyond the current DS routing 
C                   interval.  (frinstance, if weekly outflow from
C                   28nov-4dec goes to monthly period downstream,
C                   finishing off november, then DaysTooFar = 4.)
C
      INTEGER       Weight
C                   The number of days in US routing period.
C
      INTEGER       Verbosity
C
C
C     Declarations regarding COMMON block variables:
C
      INTEGER       PeriodEnd(6,3)
      INTEGER       Day(3), Q(4,3), DaysSinceSaturday
      DOUBLE PRECISION dT(3)
      LOGICAL       NumIncsVary(3) 
C                   See descriptions in subroutine RUNIT().
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
      INTEGER       MidLakeSolution, StMarysSolution, UStLawrSolution
      INTEGER       MonthLength
      INTEGER       BOMRound, FlowRound, LevelRound, P77ACompat
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
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
C        
      DATA NameOfRiver /'St. Marys','Niagara','Upper St. Lawrence'/
      DATA NameOfRiverNum / 9, 8, 12/
C
C     Make absolutely sure that USLoop + 1 is valid
      IF(USLoop .LT. 1 .AND. USLoop .GT. 2) 
     &   STOP 'what loop downstream?'
C
C     Figure out which message level to use (that of receiving body)
      IF( USLoop + 1 .EQ. 2 ) THEN
         Verbosity = MLVerbosity
      ELSE 
         Verbosity = OntVerbosity
      ENDIF
C
C
      IF(Verbosity .GE. 8) THEN
         WRITE(CTMP,8956) NameOfRiver(USLoop)(1:NameOfRiverNum(USLoop)),
     &      RoutingInterval(USLoop), RoutingInterval(USLoop + 1), LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,8959)'entering', AccumulatedFlow, DaysInAccum, 
     &      WeightedFlow, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,8957) USFlow, PeriodEnd(5,USLoop), Weight, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
C     If US period ending before DS period ends, DaysTooFar < 0.
C     If US period ending on same day as DS period, DaysTooFar == 0.
C     If US period ending after DS period ends, DaysTooFar > 0.
C     
      Weight = -999
      DaysTooFar = PeriodEnd(5,USLoop) - PeriodEnd(5,USLoop + 1)
C     USFlow, AccumulatedFlow, DaysInAccum, WeightedFlow
C
      IF( RoutingInterval(USLoop) .EQ. 1 ) THEN
C        US outflow is daily.  Weight is one.
C         WRITE(*,'(A,I4)')  '  Weight=', Weight 
         Weight = 1
         WRITE(CTMP,'(A,I4)')  '  Weight=', Weight 
C
      ELSE IF( RoutingInterval(USLoop) .EQ. 2 ) THEN
C        US outflow is weekly.  Weight is 7.
         Weight = 7
         IF( RoutingInterval(USLoop+1) .EQ. 1) DaysTooFar = 0
C           Weekly US and daily DS.  Use the weekly flow for
C           all downstream periods, until get next weekly flow.
C
      ELSE IF( RoutingInterval(USLoop) .EQ. 3 ) THEN
C        US outflow is quarter-monthly.  Weight is 
C        day-of-month at end of quarter-month less day-of-
C        month at start of quarter-month, plus one.
         Weight = PeriodEnd(1,USLoop) - Day(USLoop) + 1
         IF( RoutingInterval(USLoop+1) .EQ. 1) DaysTooFar = 0
C           Quarter-Monthly US and daily DS.  Use the qmonthly 
C           flow for all downstream periods, until get next 
C           qmonthly flow.
C
      ELSE IF( RoutingInterval(USLoop) .EQ. 4 ) THEN
C        Special handling for monthly routing, just use monthly
C        mean flow, unless downstream interval is weekly.
         IF( RoutingInterval(USLoop+1) .EQ. 2) THEN
C           Monthly US and weekly DS.  Calculate the weighted weekly 
C           flow for the pending DS routing interval in WeightedFlow.
C           Subsequent weeks within the month must draw upon the 
C           monthly average US flow (FWRK(22,?)).
            WeightedFlow= (AccumulatedFlow + USFlow*(7-DaysInAccum)) / 7
C           If the DS week ends before the US month ends, which 
C           should be all the time except if the first week straddles
C           a month, then need to look ahead beyond current DS week
C           to see how many days the final week of the month will have. 
            DaysTooFar = PeriodEnd(5,USLoop) - PeriodEnd(5,USLoop+1)
            IF( DaysTooFar .GT. 0) THEN
               Weight = MOD(DaysTooFar, 7) 
            ELSE
C              I think? this only possible for first month
               Weight = 7 - PeriodEnd(1,USLoop+1)
            ENDIF 
C
C           Set up for next time
            AccumulatedFlow = USFlow * Weight
            DaysInAccum = Weight
         ELSE
            WeightedFlow = USFlow
         ENDIF
         GOTO 7654
      ELSE
         STOP 'Ethel!'
      ENDIF
C
      IF(Verbosity .GE. 8) THEN
         WRITE(CTMP,8958) PeriodEnd(5,USLoop+1), DaysTooFar, Weight, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
C     Update the accumulated value
      AccumulatedFlow = AccumulatedFlow + USFlow * Weight
      DaysInAccum = DaysInAccum + Weight
C
      IF( DaysTooFar .GE. 0)THEN
         IF(DaysInAccum .EQ. DaysTooFar) STOP 'Fred!'
C        Finished off a DS routing period, so need to put a value
C        in WeightedFlow, and reset accumulators.
         WeightedFlow = (AccumulatedFlow - USFlow * DaysTooFar) /
     &      (DaysInAccum - DaysTooFar)
C
C        If DaysTooFar == 0, these get set to 0.  Otherwise hold
C        the weighted overlap.  
         AccumulatedFlow = USFlow * DaysTooFar
         DaysInAccum = DaysTooFar
      ELSE
C        Shouldn't need to do anything else if DaysTooFar < 0, because 
C        control loops should keep cycling through USLoop until DaysTooFar
C        >= 0.  But put a bogus flow in WeightedFlow just in case the
C        control logic fails.  This should result in catastrophic error,
C        and will hopefully crash the model.  It sounds horrible, but in
C        truth spectacular failure is more user-friendly than an insidious,
C        hard to detect bug.
         WeightedFlow = -99999999999.D0
      ENDIF
C
 7654 IF(Verbosity .GE. 8) THEN
         WRITE(CTMP,8959)'leaving', AccumulatedFlow, DaysInAccum, 
     &      WeightedFlow, LLT
         CALL LogIt(LOGD,0,1)
      ENDIF
C
 8956 FORMAT(' Converting ', A, ' River flow from routing interval ',
     &   I1, ' to ', I1, A1)
 8957 FORMAT(10X, ' flow = ', F9.2, '  upstream period ends on day ', 
     &   I10, ' Weight=', I5, A1)
 8958 FORMAT(10X, ' downstream period ends on day ', I10, 
     &   ' DaysTooFar=', I3, ' weight=', I4, A1)
 8959 FORMAT(10X, A, ' accum flow:',  F15.2, ' days:', I2, 
     &   ' weighted:', F15.2, A1)
      END
C
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
C
      SUBROUTINE TallyMeans(Lake, Loop, Flow, Level)
C
C     Call this routine after computing a connecting channel flow.  Does
C     the arithmetic for calculating the mean flows and levels, and 
C     writes them to the output files.
C
      IMPLICIT NONE 
C
C     Declarations regarding function arguments:
C
      INTEGER       Lake, Loop
C                   (Input) Lake number corresponding to element in
C                   LakesToSolve, and routing loop number according to
C                   PeriodEnd.
C                                 Name   Lake   Loop   FWRK Column
C                        Lake Superior      1      1          5
C                 Lakes Michigan-Huron      2      2          1
C                       Lake St. Clair      3      2          2
C                            Lake Erie      4      2          3
C                         Lake Ontario      5      3          6
C
      DOUBLE PRECISION Flow
C                   (Input) The outflow coming from the lake.
C
      DOUBLE PRECISION Level
C                   (Input) The lake level.
C
C     Declarations regarding local variables or functions:
C
      INTEGER       I, J, K
C                   Scratch
C
      INTEGER       Weight 
C                   The number of days in routing period of input data.
C
      INTEGER       DaysTooFar
C                   Number of days that the current level/outflow from
C                   Lake extends beyond the averaging period.  Frinstance, 
C                   if weekly flow (from 28nov-4dec goes to monthly
C                   averaging period, finishing off november, then 
C                   DaysTooFar = 4.  If qmonthly routing loop, then not
C                   days, but qmonths.
C
      INTEGER       LakeToFWRK(8)
C                   Gives the column number of the FWRK array that 
C                   corresponds to Lake.  Frinstance, for Ontario, Lake 
C                   = 5.  LakeToFWRK(5) = 6, which can be used to refer
C                   to the Ontario column in FWRK.  So 
C                   FWRK(26, LakeToFWRK(5)) references the Lake Ontario
C                   quarter-monthly mean level value.  
C
      INTEGER       DaysInAccum
      DIMENSION     DaysInAccum(2:4,8)
C                   Number of days accumulated for each mean for each
C                   lake since last time the mean was computed.  The
C                   x-subscript is for the mean (2=weekly, 3=qmonthly,
C                   4=monthly) while the y-subscript is for the lake.
C
      INTEGER       MBufRow(8)
C                   Row of MonthlyOutputBuf currently being filled out
C                   for each Lake.  When row 10 is complete, write data 
C                   to disk and flush MonthlyOutputBuf.
C                   
      INTEGER       QBufBlok(8)
C                   Which block of 4 rows of QMonthlyOutputBuf is 
C                   currently being filled out for each Lake.  When row 
C                   10 is complete, write data to disk and flush 
C                   QMonthlyOutputBuf.
C                   
      INTEGER       WBufRow(8)
C                   Row of WeeklyOutputBuf currently being filled out
C                   for each Lake.  When row 10 is complete, write data 
C                   to disk and flush WeeklyOutputBuf.
C                   
      INTEGER       DBufRow(8)
C                   Row of DailyOutputBuf currently being filled out
C                   for each Lake.  When row 48 is complete, write data 
C                   to disk and flush DailyOutputBuf.
C                   
      INTEGER       QMonth, Month, Year
C                   Quarter-Month, Month and Year of the monthly 
C                   average computed.
C
      INTEGER       DailyField
C                   Keeps track of how many values are on a line of 
C                   daily output for each lake, so know what 
C                   columns to write to.  Unneeded for other output 
C                   intervals, because can be figured out from the 
C                   qmonth and month.
C
      CHARACTER*8   LevelForm, FlowForm
C                   FORTRAN format specifiers to use for levels and flows.  
C
      DOUBLE PRECISION ROUND
C
      INTEGER       Verbosity
C
C
C     Declarations regarding COMMON block variables:
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
C                   Character buffer holding two years of mean weekly
C                   levels and flows for each of the lakes.  Data gets 
C                   saved in char buffer to reduce memory requirements,
C                   and written in 2-year chunks to each file to help
C                   performance.
C
      CHARACTER*19  WeekOutputBuf(52,2,8)
C                   Character buffer holding two years of mean weekly
C                   levels and flows for each of the lakes.  Data gets 
C                   saved in char buffer to reduce memory requirements,
C                   and written in (about) 1-year chunks to each file 
C                   to help performance.
C
      CHARACTER*74  DailyOutputBuf(48,2,8)
C                   Character buffer holding one year of mean daily
C                   levels and flows for each of the lakes.  Data gets 
C                   saved in char buffer to reduce memory requirements,
C                   and written in 1-year chunks to each file to help
C                   performance.
C
C     Don't really need these aliases.  Program deals directly with the
C     values in FWRK.  But it's nice to refer to them to see what is 
C     where.  So just comment them out for now.
C      DOUBLE PRECISION MHWeeklyMeanF, MHWeeklyAccumF, MHQmonthlyMeanF,
C     & MHQmonthlyAccumF, MHMonthlyMeanF, MHMonthlyAccumF,  
C     & SCWeeklyMeanF, SCWeeklyAccumF, SCQmonthlyMeanF, SCQmonthlyAccumF, 
C     & SCMonthlyMeanF, SCMonthlyAccumF, ERWeeklyMeanF, ERWeeklyAccumF, 
C     & ERQmonthlyMeanF,  ERQmonthlyAccumF, ERMonthlyMeanF, 
C     & ERMonthlyAccumF, SPWeeklyMeanF, SPWeeklyAccumF, SPQmonthlyMeanF,
C     & SPQmonthlyAccumF, SPMonthlyMeanF, SPMonthlyAccumF, ONWeeklyMeanF,
C     & ONWeeklyAccumF, ONQmonthlyMeanF, ONQmonthlyAccumF,ONMonthlyMeanF,
C     & ONMonthlyAccumF, MHWeeklyMeanL, MHWeeklyAccumL, MHQmonthlyMeanL,
C     & MHQmonthlyAccumL, MHMonthlyMeanL, MHMonthlyAccumL,  
C     & SCWeeklyMeanL, SCWeeklyAccumL, SCQmonthlyMeanL, SCQmonthlyAccumL, 
C     & SCMonthlyMeanL, SCMonthlyAccumL, ERWeeklyMeanL, ERWeeklyAccumL, 
C     & ERQmonthlyMeanL,  ERQmonthlyAccumL, ERMonthlyMeanL, 
C     & ERMonthlyAccumL, SPWeeklyMeanL, SPWeeklyAccumL, SPQmonthlyMeanL,
C     & SPQmonthlyAccumL, SPMonthlyMeanL, SPMonthlyAccumL, ONWeeklyMeanL,
C     & ONWeeklyAccumL, ONQmonthlyMeanL, ONQmonthlyAccumL,ONMonthlyMeanL,
C     & ONMonthlyAccumL
C 
C
      DOUBLE PRECISION FWRK
C                   This is a common block used for 64-bit calculations.
      DIMENSION FWRK(225,12)
C                   Dimensioned it this way to provide a 225 element
C                   array for each pool.
C       
      INTEGER       PeriodEnd(6,3)
      INTEGER       Day(3), Q(4,3), DaysSinceSaturday
      DOUBLE PRECISION dT(3)
      LOGICAL       NumIncsVary(3) 
C                   See descriptions in subroutine RUNIT().
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
      INTEGER       BOMRound , FlowRound, LevelRound
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
      CHARACTER*21  LAKENAME(16), RiverName(16)
      INTEGER       NLN(16), LkToSlv2LN(8), NRN(16), LkSlv2RtLoopI(8)
C
      COMMON /CLDATE/   dT, Day, Q, PeriodEnd, NumIncsVary,
     &                  DaysSinceSaturday
      COMMON /RCRUNI/   LakesToSolve, 
     &                  IntervalNameLen, MonthLength,
     &                  StartYear, StartMonth, StartDay, 
     &                  EndYear, EndMonth, EndDay
      COMMON /MISINT/   MidLakeSolution, BOMRound , FlowRound, 
     &                  LevelRound, RoutingInterval, RoutingIncrement,
     &                  StMarysSolution, UStLawrSolution, P77ACompat
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /FSPACE/   FWRK
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /MONAME/   MonthAbbrev
      COMMON /LNAMEC/   LAKENAME, RiverName
      COMMON /LNAMEI/   NLN, LkToSlv2LN, NRN, LkSlv2RtLoopI
      COMMON /OUTTXT/   MonthOutputBuf, WeekOutputBuf
      COMMON /OUTTX2/   DailyOutputBuf
      COMMON /OUTTX3/   QMonthOutputBuf
C        
C----&------------------------------------------------------------------
C
C     Assign convenient aliases for the period average levels, using 
C     storage actually allocated to the common workspace.  The levels
C     will be initialized in FRWK(5,?) and computed by this subroutine.
C     EQUIVALENCE (MHMLV, FWRK(5,1))
C     EQUIVALENCE (SCMLV, FWRK(5,2))
C     EQUIVALENCE (ERMLV, FWRK(5,3))
C     CGMLV is actually input variable; gets read instead of computed.
C     EQUIVALENCE (CGMLV, FWRK(5,4))
C     EQUIVALENCE (SPMLV, FWRK(5,5))
C
C     Assign aliases for the terms required to compute mean outflows
C     and levels, using storage actually allocated to the common workspace.  
C     Don't really need these aliases.  Program deals directly with the
C     values in FWRK.  But it's nice to refer to them to see what is 
C     where.  So just comment them out for now.
C      EQUIVALENCE (FWRK(18,1), MHWeeklyMeanF)
C      EQUIVALENCE (FWRK(19,1), MHWeeklyAccumF)
C      EQUIVALENCE (FWRK(20,1), MHQmonthlyMeanF)
C      EQUIVALENCE (FWRK(21,1), MHQmonthlyAccumF)
C      EQUIVALENCE (FWRK(22,1), MHMonthlyMeanF)
C      EQUIVALENCE (FWRK(23,1), MHMonthlyAccumF)
C      EQUIVALENCE (FWRK(18,2), SCWeeklyMeanF)
C      EQUIVALENCE (FWRK(19,2), SCWeeklyAccumF)
C      EQUIVALENCE (FWRK(20,2), SCQmonthlyMeanF)
C      EQUIVALENCE (FWRK(21,2), SCQmonthlyAccumF)
C      EQUIVALENCE (FWRK(22,2), SCMonthlyMeanF)
C      EQUIVALENCE (FWRK(23,2), SCMonthlyAccumF)
C      EQUIVALENCE (FWRK(18,3), ERWeeklyMeanF)
C      EQUIVALENCE (FWRK(19,3), ERWeeklyAccumF)
C      EQUIVALENCE (FWRK(20,3), ERQmonthlyMeanF)
C      EQUIVALENCE (FWRK(21,3), ERQmonthlyAccumF)
C      EQUIVALENCE (FWRK(22,3), ERMonthlyMeanF)
C      EQUIVALENCE (FWRK(23,3), ERMonthlyAccumF)
C     No applicable term for CGIP; storage in FWRK(18-23,4) unused.
C      EQUIVALENCE (FWRK(18,5), SPWeeklyMeanF)
C      EQUIVALENCE (FWRK(19,5), SPWeeklyAccumF)
C      EQUIVALENCE (FWRK(20,5), SPQmonthlyMeanF)
C      EQUIVALENCE (FWRK(21,5), SPQmonthlyAccumF)
C      EQUIVALENCE (FWRK(22,5), SPMonthlyMeanF)
C      EQUIVALENCE (FWRK(23,5), SPMonthlyAccumF)
C      EQUIVALENCE (FWRK(18,6), ONWeeklyMeanF)
C      EQUIVALENCE (FWRK(19,6), ONWeeklyAccumF)
C      EQUIVALENCE (FWRK(20,6), ONQmonthlyMeanF)
C      EQUIVALENCE (FWRK(21,6), ONQmonthlyAccumF)
C      EQUIVALENCE (FWRK(22,6), ONMonthlyMeanF)
C      EQUIVALENCE (FWRK(23,6), ONMonthlyAccumF)
C      EQUIVALENCE (FWRK(24,1), MHWeeklyMeanL)
C      EQUIVALENCE (FWRK(25,1), MHWeeklyAccumL)
C      EQUIVALENCE (FWRK(26,1), MHQmonthlyMeanL)
C      EQUIVALENCE (FWRK(27,1), MHQmonthlyAccumL)
C      EQUIVALENCE (FWRK(28,1), MHMonthlyMeanL)
C      EQUIVALENCE (FWRK(29,1), MHMonthlyAccumL)
C      EQUIVALENCE (FWRK(24,2), SCWeeklyMeanL)
C      EQUIVALENCE (FWRK(25,2), SCWeeklyAccumL)
C      EQUIVALENCE (FWRK(26,2), SCQmonthlyMeanL)
C      EQUIVALENCE (FWRK(27,2), SCQmonthlyAccumL)
C      EQUIVALENCE (FWRK(28,2), SCMonthlyMeanL)
C      EQUIVALENCE (FWRK(29,2), SCMonthlyAccumL)
C      EQUIVALENCE (FWRK(24,3), ERWeeklyMeanL)
C      EQUIVALENCE (FWRK(25,3), ERWeeklyAccumL)
C      EQUIVALENCE (FWRK(26,3), ERQmonthlyMeanL)
C      EQUIVALENCE (FWRK(27,3), ERQmonthlyAccumL)
C      EQUIVALENCE (FWRK(28,3), ERMonthlyMeanL)
C      EQUIVALENCE (FWRK(29,3), ERMonthlyAccumL)
C     No applicable term for CGIP; storage in FWRK(24-29,4) unused.
C      EQUIVALENCE (FWRK(24,5), SPWeeklyMeanL)
C      EQUIVALENCE (FWRK(25,5), SPWeeklyAccumL)
C      EQUIVALENCE (FWRK(26,5), SPQmonthlyMeanL)
C      EQUIVALENCE (FWRK(27,5), SPQmonthlyAccumL)
C      EQUIVALENCE (FWRK(28,5), SPMonthlyMeanL)
C      EQUIVALENCE (FWRK(29,5), SPMonthlyAccumL)
C      EQUIVALENCE (FWRK(24,6), ONWeeklyMeanL)
C      EQUIVALENCE (FWRK(25,6), ONWeeklyAccumL)
C      EQUIVALENCE (FWRK(26,6), ONQmonthlyMeanL)
C      EQUIVALENCE (FWRK(27,6), ONQmonthlyAccumL)
C      EQUIVALENCE (FWRK(28,6), ONMonthlyMeanL)
C      EQUIVALENCE (FWRK(29,6), ONMonthlyAccumL)
C      EQUIVALENCE (FWRK(30,1), MHEndofMonthL)
C      EQUIVALENCE (FWRK(30,2), SCEndofMonthL)
C      EQUIVALENCE (FWRK(30,3), EREndofMonthL)
C      No applicable term for CGIP; storage in FWRK(30,4) unused.
C      EQUIVALENCE (FWRK(30,4), ??EndofMonthL)
C      EQUIVALENCE (FWRK(30,5), SPEndofMonthL)
C      EQUIVALENCE (FWRK(30,6), ONEndofMonthL)
C
C     This isn't done at all.  Just kept on in case want to see how
C     to alias arrays to FWRK.
C     DIMENSION MHMonthlyLevels(12), SCMonthlyLevels(12), 
C    &   ERMonthlyLevels(12), SPMonthlyLevels(12), ONMonthlyLevels(12),
C    &   MHMonthlyFlows(12), SCMonthlyFlows(12), ERMonthlyFlows(12),
C    &   SPMonthlyFlows(12), ONMonthlyFlows(12)
C
C     Storage in FWRK(30-53,?) used as holding space to get 12 months
C     worth of data before writing levels or flows out.
C     EQUIVALENCE (FWRK(30,1), MHMonthlyLevels)
C     EQUIVALENCE (FWRK(42,1), MHMonthlyFlows)
C     EQUIVALENCE (FWRK(30,2), SCMonthlyLevels)
C     EQUIVALENCE (FWRK(42,2), SCMonthlyFlows)
C     EQUIVALENCE (FWRK(30,3), ERMonthlyLevels)
C     EQUIVALENCE (FWRK(42,3), ERMonthlyFlows)
C     Storage in FWRK(30,4) to FWRK(53,4) unused
C     EQUIVALENCE (FWRK(30,5), SPMonthlyLevels)
C     EQUIVALENCE (FWRK(42,5), SPMonthlyFlows)
C     EQUIVALENCE (FWRK(30,6), ONMonthlyLevels)
C     EQUIVALENCE (FWRK(42,6), ONMonthlyFlows)
C
C
      DATA LakeToFWRK / 5, 1, 2, 3, 6, 0, 0, 0 /
C
C
      DATA MBufRow / 8 * 1 / 
      DATA QBufBlok  / 8 * 1 / 
      DATA WBufRow  / 8 * 1 / 
      DATA DBufRow  / 8 * 1 / 
C
      DATA DaysInAccum / 24 * 0 / 
C
      DATA LevelForm /'(F8.2)'/
      DATA FlowForm /'(F8.0)'/
C
C     Put bogus values into monthly ave flow and level variables
      FWRK(22, LakeToFWRK(Lake)) = -999999999.D0
      FWRK(28, LakeToFWRK(Lake)) = -999999999.D0
C
C     Put bogus values into quarter-monthly ave  flow and level variables
      FWRK(20, LakeToFWRK(Lake)) = -999999999.D0
      FWRK(26, LakeToFWRK(Lake)) = -999999999.D0
C
C     Put bogus values into weekly ave  flow and level variables
      FWRK(18, LakeToFWRK(Lake)) = -999999999.D0
      FWRK(24, LakeToFWRK(Lake)) = -999999999.D0
C
      DaysTooFar = -999
C
C     Figure out which message level to use
      IF( Loop .EQ. 1 ) THEN
         Verbosity = SupVerbosity
      ELSE IF ( Loop .EQ. 2 ) THEN
         Verbosity = MLVerbosity
      ELSE 
         Verbosity = OntVerbosity
      ENDIF
C
      IF( RoutingInterval(Loop) .EQ. 1 ) THEN
C        Routed level/flow are daily.  Weight is one.
         Weight = 1
C        IF(Verbosity .GE. 8) THEN
C           WRITE(CTMP, '(3(A, I3), A1)') ' PeriodEnd(1,Loop)=', 
C    &         PeriodEnd(1,Loop), ' Q(4,Loop)=', Q(4,Loop), 
C    &         ' DaysTooFar=', DaysTooFar, LLT
C           CALL LogIt(LOGD,0,0)
C        ENDIF
C
      ELSE IF( RoutingInterval(Loop) .EQ. 2 ) THEN
C        Routed level/flow are weekly.  Weight is 7.
         Weight = 7
C
      ELSE IF( RoutingInterval(Loop) .EQ. 3 ) THEN
C        Routed level/flow are quarter-monthly.  Weight is day-of-month 
C        at end of quarter-month less day-of-month at start of 
C        quarter-month, plus one.
         Weight = PeriodEnd(1,Loop) - Day(Loop) + 1
C        Period ave levels are quarter-monthly mean levels:
C         FWRK(26, LakeToFWRK(Lake)) = FWRK(5, LakeToFWRK(Lake))
C        Period ave flows are quarter-monthly mean flows:
C         FWRK(20, LakeToFWRK(Lake)) = FWRK(6, LakeToFWRK(Lake))
C
      ELSE
C        Routed level/flow are monthly.  Just assign period average 
C        values, and move on to write-out portion.  Weight shouldn't 
C        values, used (and since assigning 0, will puke if used in   
C        denom).  
         Weight = 0
C        Period ave levels are monthly mean levels:
C         FWRK(28, LakeToFWRK(Lake)) = FWRK(5, LakeToFWRK(Lake))
C        Period ave flows are monthly mean flows:
C         FWRK(22, LakeToFWRK(Lake)) = FWRK(6, LakeToFWRK(Lake))
      ENDIF
C
C     Store provisional values for qmonth, month, and year
C     These will later be incremented backward if DaysTooFar > 0
      QMonth = PeriodEnd(2,Loop)
      Month = PeriodEnd(3,Loop)
      Year = PeriodEnd(4,Loop)
C
      IF(Verbosity .GE. 8) THEN
         WRITE(CTMP,'(10X, A,A,A,I2,A3,I5,A,I1,A,I2,A1)') 
     &      'Entering TallyMeans for ', 
     &      LakeName(LkToSlv2LN(Lake))(1:NLN(LkToSlv2LN(Lake))), 
     &      ' on ', PeriodEnd(1,Loop), 
     &      MonthAbbrev(PeriodEnd(3,Loop)), PeriodEnd(4,Loop), 
     &      ' with routing interval ', RoutingInterval(Loop), 
     &      ' and weight ', weight, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
      IF( RoutingInterval(Loop) .EQ. 4)THEN
C        Routed level/flow are monthly.  Just assign period average 
C        values, and move on to write-out portion. 
C        Period ave levels are monthly mean levels:
         FWRK(28, LakeToFWRK(Lake)) = Level
C        Period ave flows are monthly mean flows:
         FWRK(22, LakeToFWRK(Lake)) = Flow
C      EQUIVALENCE (FWRK(30,5), SPEndofMonthL)
C      EQUIVALENCE (FWRK(8,1), MHEOP)
         FWRK(30, LakeToFWRK(Lake)) = FWRK(8, LakeToFWRK(Lake))
      ELSE
C        So incoming level/flow is daily, weekly, or qmonthly.
C        Must compute monthly average.
C
         IF(RoutingInterval(Loop) .EQ. 3 .AND. MonthLength .EQ. 1) THEN
C           Concept of DaysTooFar applies, but qmons instead of days.
C           So if the quarter-month == 4, then DaysTooFar == 0
            DaysTooFar = PeriodEnd(2,Loop) - 4
            WRITE(CTMP,'(10X,A,3(A,I3),A1)')'@@',' PeriodEnd(2,Loop)=', 
     &         PeriodEnd(2,Loop), ' Q(4,Loop)=', Q(4,Loop), 
     &         ' DaysTooFar=', DaysTooFar, LLT
C           Uniform quarter-month length, so weight doesn't matter, as
C           long as it's equal each time.
            Weight = 1 
         ELSE IF(RoutingInterval(Loop) .EQ. 2 .AND. 
     &      PeriodEnd(1,Loop) .LE. 6 .AND. PeriodEnd(5,Loop) .GT. 7)THEN
C           If weekly routing with day-of-month of 6 or less, then 
C           starting a new month, unless initial period.  
            DaysTooFar = PeriodEnd(1,Loop)
            WRITE(CTMP,'(A,I3)')'@@',DaysTooFar
C 
C        XXX not sure this logic handles odd start dates all that well
         ELSE
C           Daily or Qmonthly, or weekly after first week in month.  If 
C           day-of-month is last day in month, then DaysTooFar == 0,
C           else < 0.
C
C           Something very odd going on in this routine.  Assigment of
C           DaysTooFar fails to happen (can be observed in debugger) unless
C           do something with the variable.  
            DaysTooFar = PeriodEnd(1,Loop) - Q(4,Loop)
            WRITE(CTMP,'(A,I3)')'@@',DaysTooFar
C
C           Now see if weekly routing on Lake Superior using Plan 1977A,
C           which needs to know the previous month's outflow.  This 
C           wouldn't be available except when week ends on last day of
C           month.  Need to estimate the monthly outflow when get to
C           second-to-last week of month, so it's available when Plan 
C           1977A runs in the next week.  Just presume that remainder
C           of month is same flow as this week. 
            IF(DaysTooFar .GT. -7 .AND. Loop .EQ. 1 .AND.
     &         RoutingInterval(Loop) .EQ. 2 )THEN
C              EQUIVALENCE (FWRK(22,5), SPMonthlyMeanF)
               FWRK(22,5)= (FWRK(23,5) + Flow * (Weight - DaysTooFar)) /
     &            (DaysInAccum(4,Lake) + Weight - DaysTooFar)
C             
            ENDIF
C
C
         ENDIF
C
         IF( DaysTooFar .GT. 0)THEN
C           This should only happen if data from weekly routing
            IF(RoutingInterval(Loop) .NE. 2) STOP ' HTH?'
            Month = Month - 1
            IF(Month .EQ. 0 ) THEN 
               Month = 12
               Year = Year - 1
            ENDIF
         ENDIF
C
         IF(Verbosity .GE. 8) THEN
            WRITE(CTMP,8959)'entering', ' monthly', 'level', 
     &         FWRK(29,LakeToFWRK(Lake)), FWRK(28,LakeToFWRK(Lake)),LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,8959)'entering', ' monthly', ' flow',
     &         FWRK(23,LakeToFWRK(Lake)), FWRK(22,LakeToFWRK(Lake)),LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,8958) DaysInAccum(4,Lake), 'month',DaysTooFar,LLT
            CALL LogIt(LOGD,0,0)
C           WRITE(CTMP, '(10X, 2(A, I3), A1)') ' PeriodEnd(1,Loop)=', 
C    &         PeriodEnd(1,Loop), ' Q(4,Loop)=', Q(4,Loop), LLT
C           CALL LogIt(LOGD,0,0)
         ENDIF
C
C        Update the monthly accumulated values for the lake.
         FWRK(23, LakeToFWRK(Lake)) = FWRK(23, LakeToFWRK(Lake)) + 
     &      Flow * Weight
         FWRK(29, LakeToFWRK(Lake)) = FWRK(29, LakeToFWRK(Lake)) + 
     &      Level * Weight
         DaysInAccum(4,Lake) = DaysInAccum(4,Lake) + Weight
C
         IF( DaysTooFar .GE. 0)THEN
C           Finished off a DS routing period, so need to put a value
C           in monthly mean flow and level, and reset accumulators.
            IF(DaysTooFar .EQ. DaysInAccum(4,Lake)) STOP 'Maggie'
            FWRK(22, LakeToFWRK(Lake)) = (FWRK(23, LakeToFWRK(Lake)) 
     &         - Flow * DaysTooFar) / (DaysInAccum(4,Lake) - DaysTooFar)
            FWRK(28, LakeToFWRK(Lake)) = (FWRK(29, LakeToFWRK(Lake)) 
     &         - Level* DaysTooFar) / (DaysInAccum(4,Lake) - DaysTooFar)
C
C
C      EQUIVALENCE (FWRK(30,5), SPEndofMonthL)
C      EQUIVALENCE (FWRK(7,1), MHBOP)
C      EQUIVALENCE (FWRK(8,1), MHEOP)
C           Any forced EOP levels should already be in FWRK(8,?).
C      XX   Must also determine End-of-month level for use with
C      XX   regulation plans.  Should work out correctly for qmonthly
C      XX   daily routing data, but linear interpolation of weekly
C      XX   might introduce some error.  Maybe reduce error by 
C      XX   interpolating from period average (weekly mean) level
C      XX   to the BOP or EOP levels.   
C           
C           If weekly routing, back off from EOP level according to 
C           ratio of distance from BOP level and DaysTooFar / 7.  If 
C           daily or qmonthly routing, then DaysTooFar == 0, and
C           End-of-Month level works out to EOP level, so the equation 
C           still works.  
C           XX OTOH, why bother at all?  The only reason
C           XX we care about BOM levels is for Plan 1977A, right?  
C           XX And Plan 1977A will have already run by the time the
C           XX overlapping week routes and execution arrives here.  
C           XX We can leave it in this way for now, but unless need 
C           XX BOM levels for Ontario, should probably just assign
C           XX FWRK(30,?) to FWRK(8,?) sometime in the future.
C           XX Better place to set up weekly exception is in 
C           XX RouteSuperior, before call to Plan 1977A.         
C           XX PS  It's there, but didn't disable here.  Not hurting
C           XX anything - I think.  Need some testing with this!
            FWRK(30, LakeToFWRK(Lake)) = FWRK(8, LakeToFWRK(Lake)) -
     &         (FWRK(8, LakeToFWRK(Lake)) - FWRK(7, LakeToFWRK(Lake)))
     &         * DaysTooFar / 7
C
C           If DaysTooFar == 0, these get set to 0.  Otherwise they hold
C           the weighted overlap.  
            FWRK(23, LakeToFWRK(Lake)) = Flow * DaysTooFar
            FWRK(29, LakeToFWRK(Lake)) = Level * DaysTooFar
            DaysInAccum(4,Lake) = DaysTooFar
         ENDIF
C
C----&------------------------------------------------------------------
C        That does it for accumulating monthly values.
C        Now do the same for qmonthly.    Flow and level already set if 
C        incoming (routed)
C        data is quarter-monthly, so only do accumulating if routing 
C        is weekly or daily.
         IF( RoutingInterval(Loop) .EQ. 3)THEN
C           Routed level/flow are quarter-monthly.  
C           Period ave levels are quarter-monthly mean levels:
            FWRK(26, LakeToFWRK(Lake)) = Level
C           Period ave flows are quarter-monthly mean flows:
            FWRK(20, LakeToFWRK(Lake)) = Flow
C            WRITE(CTMP,'(A, F9.4, A1)') ' here i am assigning ', 
C     &         FWRK(26, LakeToFWRK(Lake)), LLT    
C            CALL LogIt(LOGD,0,0)
C            DO 833, I=1,8 
C               WRITE(CTMP,'(A, I1, A, F9.4, A1)') '    FWRK(5,', I, 
C     &         ')= ', FWRK(5,I), LLT
C  833          CALL LogIt(LOGD,0,0)
C
         ELSE 
C           Must compute quarter-monthly average from weekly or daily.
C
            IF(RoutingInterval(Loop) .EQ. 2 )THEN
C
               IF( PeriodEnd(2,Loop) .GT. 1) THEN
C                 Routed data ends in 2-4 qtr of same month, so Q() 
C                 dates still applicable to averaging period.  End of
C                 previous qtrmon given by Q(PeriodEnd(2,Loop)-1,Loop).
                  DaysTooFar = PeriodEnd(1,Loop) - 
     &               Q(PeriodEnd(2,Loop) - 1,Loop)
                  WRITE(CTMP,'(A,I3)')'@@',DaysTooFar
               ELSE
C                 Routed data ends in first qmonth of month
                  IF( PeriodEnd(1,Loop) .LE. 6 ) THEN
C                    Bridged over to new month, DaysTooFar simply = date
                     DaysTooFar = PeriodEnd(1,Loop)
                     WRITE(CTMP,'(A,I3)')'@@',DaysTooFar
                  ELSE
C                    Still in same qmonth; if week ends on last day of
C                    qtrmonth then DaysTooFar == 0, else < 0.
                     DaysTooFar= PeriodEnd(1,Loop) - Q(1,Loop)
                     WRITE(CTMP,'(A,I3)')'@@',DaysTooFar
                  ENDIF
               ENDIF
C
C              Now add a final check in case hit a week ending on same
C              day as quarter-month:
               IF( DaysTooFar .GE. 7 ) DaysTooFar = PeriodEnd(1,Loop) - 
     &            Q(PeriodEnd(2,Loop),Loop)
               WRITE(CTMP,'(A,I3)')'@@',DaysTooFar
C
               IF( DaysTooFar .GT. 0)THEN
                  QMonth = QMonth - 1
                  IF(QMonth .EQ. 0 ) THEN 
                     QMonth = 4
C                    Don't adjust month; should have been done during 
C                    monthly accumulating arithmetic.
C                     Month = Month - 1
C                     IF(Month .EQ. 0 ) THEN 
C                        Month = 12
C                        Year = Year - 1
C                     ENDIF
                  ENDIF
               ENDIF
C
            ELSE 
C              Daily - if day-of-month is last day in qmonth, then 
C              DaysTooFar == 0, else < 0.
            DaysTooFar= PeriodEnd(1,Loop) - Q(PeriodEnd(2,Loop),Loop)
            WRITE(CTMP,'(A,I3)')'@@',DaysTooFar
            ENDIF
C
            IF(Verbosity .GE. 8) THEN
               WRITE(CTMP,8959)'entering', ' qmonthly', 'level', 
     &            FWRK(27,LakeToFWRK(Lake)), 
     &            FWRK(26,LakeToFWRK(Lake)), LLT
                CALL LogIt(LOGD,0,0)
               WRITE(CTMP,8959)'entering', ' qmonthly', ' flow',
     &            FWRK(21,LakeToFWRK(Lake)), 
     &            FWRK(20,LakeToFWRK(Lake)), LLT
               CALL LogIt(LOGD,0,0)
               WRITE(CTMP,8958) DaysInAccum(3,Lake), 'qmonth',
     &            DaysTooFar, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
C
C           Update the quarter-monthly accumulated values for the lake.
            FWRK(21, LakeToFWRK(Lake)) = FWRK(21, LakeToFWRK(Lake)) + 
     &         Flow * Weight
            FWRK(27, LakeToFWRK(Lake)) = FWRK(27, LakeToFWRK(Lake)) + 
     &         Level * Weight
            DaysInAccum(3,Lake) = DaysInAccum(3,Lake) + Weight
C
            IF( DaysTooFar .GE. 0)THEN
C              Finished off a DS routing period, so need to put a value
C              in qmonthly mean flow and level, and reset accumulators.
               FWRK(20, LakeToFWRK(Lake)) = (FWRK(21, LakeToFWRK(Lake)) 
     &            - Flow * DaysTooFar) / 
     &            (DaysInAccum(3,Lake) - DaysTooFar)
               FWRK(26, LakeToFWRK(Lake)) = (FWRK(27, LakeToFWRK(Lake)) 
     &            - Level* DaysTooFar) / 
     &            (DaysInAccum(3,Lake) - DaysTooFar)
C
C              If DaysTooFar == 0, these get set to 0.  Otherwise they 
C              hold the weighted overlap.  
               FWRK(21, LakeToFWRK(Lake)) = Flow * DaysTooFar
               FWRK(27, LakeToFWRK(Lake)) = Level * DaysTooFar
               DaysInAccum(3,Lake) = DaysTooFar
            ENDIF
C
C           That does it for accumulating quarter-monthly values.
C----&------------------------------------------------------------------
C
C           Now do the same for weekly.  Flow and level already set if 
C           incoming (routed)  data is quarter-monthly, so only do 
C           accumulating if routing is daily.
            IF( RoutingInterval(Loop) .EQ. 2)THEN
C              Routed level/flow are weekly.  
C              Period ave levels are weekly mean levels:
               FWRK(24, LakeToFWRK(Lake)) = Level
C              Period ave flows are weekly mean flows:
               FWRK(18, LakeToFWRK(Lake)) = Flow
C               WRITE(CTMP,'(A, F9.4, A1)') ' here i am assigning ', 
C     &            FWRK(24, LakeToFWRK(Lake)), LLT    
C               CALL LogIt(LOGD,0,0)
C
            ELSE 
C              Must compute weekly average from daily.  Really don't 
C              even have to bother with DaysTooFar and weight, since 
C              weight is always 1 and can only increment exactly to a
C              week boundary, not beyond it.
C              If already have 6 days accumulated, this will be 7th, so
C              DaysTooFar == 0, else < 0.
               DaysTooFar = 6 - DaysInAccum(2,Lake)
C
               IF(Verbosity .GE. 8) THEN
                  WRITE(CTMP,8959)'entering', ' weekly', 'level', 
     &               FWRK(25,LakeToFWRK(Lake)), 
     &               FWRK(24,LakeToFWRK(Lake)), LLT
                  CALL LogIt(LOGD,0,0)
                  WRITE(CTMP,8959)'entering', ' weekly', ' flow',
     &               FWRK(19,LakeToFWRK(Lake)), 
     &               FWRK(18,LakeToFWRK(Lake)), LLT
                  CALL LogIt(LOGD,0,0)
                  WRITE(CTMP,8958) DaysInAccum(2,Lake), 'week',
     &               DaysTooFar, LLT
                  CALL LogIt(LOGD,0,0)
               ENDIF
C
C              Update the weekly accumulated values for the lake.
               FWRK(19,LakeToFWRK(Lake))=FWRK(19,LakeToFWRK(Lake))+Flow 
               FWRK(25,LakeToFWRK(Lake))=FWRK(25,LakeToFWRK(Lake))+Level
               DaysInAccum(2,Lake) = DaysInAccum(2,Lake) + Weight
C
               IF( DaysInAccum(2,Lake) .EQ. 7)THEN
C                 Finished off a DS routing period, so must put a value
C                 in weekly mean flow and level, and reset accumulators.
                  FWRK(18,LakeToFWRK(Lake))= FWRK(19,LakeToFWRK(Lake))/7
                  FWRK(24,LakeToFWRK(Lake))= FWRK(25,LakeToFWRK(Lake))/7
C
                  FWRK(19, LakeToFWRK(Lake)) = 0.0D0
                  FWRK(25, LakeToFWRK(Lake)) = 0.0D0
                  DaysInAccum(2,Lake) = 0
               ENDIF
C
C              That does it for accumulating weekly values.
            ENDIF
         ENDIF
      ENDIF
C----&------------------------------------------------------------------
C
C
C     
C
      IF( FWRK(28, LakeToFWRK(Lake)) .GT. 0 )THEN
C        Now we have a monthly average
C
C        XX Keep track of the extrema, and mean (std dev?)
C
C        Write to a string; hold a bunch of strings in memory, then dump 
C        to file in a chunk for faster disk access.  Use a RAM disk?
C        Calculate the end column of field to hold this months' level:
         I = 5 + Month * 8
         WRITE(MonthOutputBuf(MBufRow(Lake),1,Lake)(I-7:I), LevelForm)  
     &      Round(FWRK(28, LakeToFWRK(Lake)), -2)
         WRITE(MonthOutputBuf(MBufRow(Lake), 2,Lake)(I-7:I), FlowForm)  
     &      Round(FWRK(22, LakeToFWRK(Lake)), 0)
C
C
         IF(Verbosity .GE. 8) THEN
            WRITE(CTMP,'(A, A, I2, A, I2, A, I3, A, A1)')    
     &         LakeName(LkToSlv2LN(Lake))(1:NLN(LkToSlv2LN(Lake))),
     &         ' monthly output row number',MBufRow(Lake), 
     &          ' columns ',I-7,' to ', I,' built to: ', LLT
            CALL LogIt(LOGD,0,0)
            DO 588, J=1,2
               WRITE(CTMP,8654) 
     &            MonthOutputBuf(MBufRow(Lake), J, Lake), LLT
  588          CALL LogIt(LOGD,0,0)
         ENDIF
C
         IF(Month .EQ. 1)THEN
C           Put the year at start of line when hit January.
            DO 589, J=1,2
  589          WRITE( MonthOutputBuf(MBufRow(Lake), J, Lake)(1:5), 
     &            '(I5)')Year  
C
         ELSE IF(Month .EQ. 12)THEN
            IF( MBufRow(Lake) .EQ. 10 )THEN
C              Monhtly output buffers full; save the data.
               CALL WriteToDisk(4, Lake) 
C              Flush out the buffers
               DO 490, J=1,2
                  DO 490, I=1,10
                     MonthOutputBuf(I,J,Lake)(1:1) = '#' 
                     DO 490, K=2,100
  490                   MonthOutputBuf(I,J,Lake)(K:K) = ' ' 
C              Reset the line counter
               MBufRow(Lake) = 1
            ENDIF
            MBufRow(Lake) = MBufRow(Lake) + 1
         ENDIF
      ENDIF
C     End of processing monthly averages
C----&------------------------------------------------------------------
C
      IF( FWRK(26, LakeToFWRK(Lake)) .GT. 0 )THEN
C        Now we have a qmonthly average
C
C        XX Keep track of the extrema, and mean (std dev?)
C
C        Write to a string; hold a bunch of strings in memory, then dump 
C        to file in a chunk for faster disk access.  Use a RAM disk?
C        Calculate the end column of field to hold this qmonths' level:
         I = 7 + Month * 8
C
         WRITE(QMonthOutputBuf(Qmonth, QBufBlok(Lake),1,Lake)(I-7:I), 
     &      LevelForm) Round(FWRK(26, LakeToFWRK(Lake)), -2)
         WRITE(QMonthOutputBuf(Qmonth, QBufBlok(Lake), 2,Lake)(I-7:I), 
     &      FlowForm) Round(FWRK(20, LakeToFWRK(Lake)), 0)
C         WRITE(CTMP,'(A)') ' here i am storing'//LLT    
C         CALL LogIt(LOGD,0,0)
C
         IF(Verbosity .GE. 8) THEN
            WRITE(CTMP,'(A, A, I1, A, I2, A, I2, A, I3, A, A1)')    
     &         LakeName(LkToSlv2LN(Lake))(1:NLN(LkToSlv2LN(Lake))),
     &         ' qmonth #', QMonth, ' block number',QBufBlok(Lake), 
     &         ' columns ',I-7,' to ', I,' built to: ', LLT
            CALL LogIt(LOGD,0,0)
            DO 688, J=1,2
               WRITE(CTMP,8654) 
     &            QMonthOutputBuf(QMonth, QBufBlok(Lake), J, Lake), LLT
  688          CALL LogIt(LOGD,0,0)
         ENDIF
C
         IF(Month .EQ. 1 .AND. QMonth .EQ. 1)THEN
C           Put the year at start of line when hit January.
            DO 689, J=1,2
               DO 689, I = 1,4
  689             WRITE(QMonthOutputBuf(I,QBufBlok(Lake), J, Lake)(1:7), 
     &            '(I5, I2)') Year, I
C
         ELSE IF(Month .EQ. 12 .AND. QMonth .EQ. 4)THEN
            IF( QBufBlok(Lake) .EQ. 10 )THEN
C              QMonthly output buffers full; save the data.
               CALL WriteToDisk(3, Lake) 
C              Flush out the buffer
               DO 690, J=1,2
                  DO 690, QMonth=1,4
                     DO 690, I=1,10
                        QMonthOutputBuf(QMonth,I,J,Lake)(1:1) = '#' 
                        DO 690, K=2,100
  690                      QMonthOutputBuf(QMonth,I,J,Lake)(K:K) = ' ' 
C              Reset the line counter
               QBufBlok(Lake) = 1
               QMonth = 4
            ENDIF
            QBufBlok(Lake) = QBufBlok(Lake) + 1
         ENDIF
      ENDIF
C----&------------------------------------------------------------------
C
C
      IF( FWRK(24, LakeToFWRK(Lake)) .GT. 0 )THEN
C        We have a weekly average
C
C        XX Keep track of the extrema, and mean (std dev?)
C
C        Put the year at start when beginning new line
         DO 789, J=1,2
  789       WRITE( WeekOutputBuf(WBufRow(Lake), J, Lake)(1:11), 
     &         '(I5,I3,I3)')Year, PeriodEnd(3,Loop), PeriodEnd(1,Loop)
C
C        CHARACTER*18  WeekOutputBuf(52,2,8)
C        Write to a string; hold a bunch of strings in memory, then dump 
C        to file in a chunk for faster disk access. 
         WRITE(WeekOutputBuf(WBufRow(Lake),1,Lake)(12:19), LevelForm)  
     &      Round(FWRK(24, LakeToFWRK(Lake)), -2)
         WRITE(WeekOutputBuf(WBufRow(Lake), 2,Lake)(12:19), FlowForm)  
     &      Round(FWRK(18, LakeToFWRK(Lake)), 0)
C
         IF(Verbosity .GE. 8) THEN
            WRITE(CTMP,'(A, A, I2, A, A1)')    
     &         LakeName(LkToSlv2LN(Lake))(1:NLN(LkToSlv2LN(Lake))),
     &         ' weekly outputs for row number',WBufRow(Lake), ': ', LLT
            CALL LogIt(LOGD,0,0)
            DO 788, J=1,2
               WRITE(CTMP,8654) 
     &            WeekOutputBuf(WBufRow(Lake), J, Lake), LLT
  788          CALL LogIt(LOGD,0,0)
         ENDIF
C
         IF( WBufRow(Lake) .EQ. 52 )THEN
C           Weekly output buffers full; save the data.
            CALL WriteToDisk(2, Lake) 
C           Flush out the buffers
            DO 790, J=1,2
               DO 790, I=1,52
                  WeekOutputBuf(I,J,Lake)(1:1) = '#' 
                  DO 790, K=2,18
  790                WeekOutputBuf(I,J,Lake)(K:K) = ' ' 
C           Reset the line counter
            WBufRow(Lake) = 1
         ELSE
            WBufRow(Lake) = WBufRow(Lake) + 1
         ENDIF
      ENDIF
C     End of processing weekly averages
C----&------------------------------------------------------------------
C
      IF( RoutingInterval(Loop) .EQ. 1 )THEN
C        We have a daily average input
C
C        Determine which field it goes in
         IF( PeriodEnd(2,Loop) .EQ. 1 ) THEN
            I = 0
         ELSE
            I = Q(PeriodEnd(2,Loop)-1,Loop) 
         ENDIF
         DailyField = PeriodEnd(1,Loop) - I
         IF(Verbosity .GE. 8) THEN
            WRITE(CTMP,'(A, 4(A, I2), A1)')    
     &         LakeName(LkToSlv2LN(Lake))(1:NLN(LkToSlv2LN(Lake))),
     &         ' DailyField:', DailyField,
     &         ' PeriodEnd(1,Loop):', PeriodEnd(1,Loop),
     &         ' PeriodEnd(2,Loop):', PeriodEnd(2,Loop),  
     &         ' Q(PeriodEnd(2,Loop)-1,Loop)', I, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
C        more than 8 fields blows 74-column line length
         IF(DailyField .GT. 8) STOP ' Ralph'
C
C        XX Keep track of the extrema, and mean (std dev?)
C
C        Write to a string; hold a bunch of strings in memory, then dump 
C        to file in a chunk for faster disk access.  Use a RAM disk?
C        Calculate the end column of field to hold this months' level:
         I = 9 + DailyField * 8
         WRITE(DailyOutputBuf(DBufRow(Lake),1,Lake)(I-7:I), LevelForm)  
     &      Round(Level, -2)
         WRITE(DailyOutputBuf(DBufRow(Lake), 2,Lake)(I-7:I), FlowForm)  
     &      Round(Flow, 0)
C
         IF(Verbosity .GE. 8) THEN
            WRITE(CTMP,'(A, A, I2, A, I2, A, I2, A, I3, A, A1)')    
     &         LakeName(LkToSlv2LN(Lake))(1:NLN(LkToSlv2LN(Lake))),
     &         ' daily output row number',DBufRow(Lake), ' field ',
     &         DailyField, ' columns ',I-7,' to ', I,' built to: ', LLT
            CALL LogIt(LOGD,0,0)
            DO 988, J=1,2
               WRITE(CTMP,8654) 
     &            DailyOutputBuf(DBufRow(Lake), J, Lake), LLT
  988          CALL LogIt(LOGD,0,0)
         ENDIF
C
         IF(DailyField .EQ. 1)THEN
C           Put the year, month, and quarter-month at start when 
C           beginning new line:
            DO 989, J=1,2
  989          WRITE( DailyOutputBuf(DBufRow(Lake), J, Lake)(1:9), 
     &            '(I4, I3, I2)')Year, PeriodEnd(3,Loop), 
     &            PeriodEnd(2,Loop)  
C
C         ELSE IF(DailyField(Lake) .EQ. 10)THEN
         ELSE IF(PeriodEnd(1,Loop) .EQ. Q(PeriodEnd(2,Loop),Loop))THEN
C           Need to start a new line next time
            IF( DBufRow(Lake) .EQ. 48 )THEN
C              Daily output buffers full; save the data.
               CALL WriteToDisk(1, Lake) 
C              Flush out the buffers
               DO 990, J=1,2
                  DO 990, I=1,48
                     DailyOutputBuf(I,J,Lake)(1:1) = '#' 
                     DO 990, K=2,74
  990                   DailyOutputBuf(I,J,Lake)(K:K) = ' ' 
C              Reset the line counter
               DBufRow(Lake) = 1
            ENDIF
            DBufRow(Lake) = DBufRow(Lake) + 1
         ENDIF
C     End of processing daily averages
      ENDIF
C
C
 8654 FORMAT(A,A1)
 8958 FORMAT(10X, I3, ' days accumulated for ', A, ' (', I3, 
     &   ' days too far)', A1)
 8959 FORMAT(10X,A,A,' accum ', A, ':',  F15.8, ' weighted:', F21.8,A1)
      END
