      SUBROUTINE P77ASolveSup(BOPLevel, NTS, dT, Gates, FlowRound, 
     &   LevelRound, EOPLevel, MeanOutflow, MeanLevel, MonthNum)
C     
C     Solves simultaneous equations for Lake Superior outflow,
C     mean level, and change-in-storage, where outflow consists of
C     flow through the compensating works and sidechannels.  Beefed-up 
C     version of EC function "FLOW" or Corps subroutine "MNLVL". 
C
C     Carries out steps 6 - 10 on pages 27 - 28 of Plan 77.  A comment
C     in the old Environment Canada model claims "A fixed number of 
C     iterations is used for coordination with other programs."
C     As far as Matt can tell, this means that four iterations are used
C     because the blue book shows an example doing it that way.  Matt 
C     thinks we should eliminate the rounding in this routine, and 
C     iterate to a tolerance.  
C
C     Any changes in this routine will affect backward-compatibility
C     with the older models.
C     
C     SUBROUTINE MNLVL(S,SM,SN,O,A1,A2,A3)
C     Revised to metric and IGLD 85, Apr 1993
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     O=(A1*(A2*S-A3)**1.5d0)*0.1d0+232.0d0
C     DO 13 L=1,4
C     CALL ROUND(((SN-O)/31.38d0),ISM)
C     CALL ROUND(S*100.0d0+ISM*0.5d0,ISM2)
C     SM=ISM2*0.01d0
C     O=(A1*(A2*SM-A3)**1.5d0)*0.1d0+232.0d0
C     CALL ROUND(O*10.0d0,ISO)
C     O=ISO*0.1d0
C  13 CONTINUE
C     RETURN
C     NOTE: MonthNum added 8/16/2013 to allow NonCompWorksFlow to access 
C     the SideChanMaxCapacity(Month) array.  It is an input variable to 
C     make sure that the month is always provided.
C
C
      IMPLICIT NONE
C
C
C     Declarations regarding subroutine arguments:
C      
C
      DOUBLE PRECISION BOPLevel, NTS, dT 
C                   (Input)  The Beginning-of-Period Superior level, 
C                   Net Total Supply (including diversions), and number
C                   of seconds in the month.  
C
      REAL          Gates
C                   (Input)  The number of gates open.  If gates given 
C                   as a negative value, then MeanOutflow is an input.
C                    
      INTEGER       FlowRound, LevelRound, MonthNum
C                   (Input)  The decimal place to round flows and levels
C                   to ( -2 means hundredths, 1 means 10's).
C                   MonthNum added as pass-through for NonCompWorksFlow
C
      DOUBLE PRECISION MeanOutflow
C                   (Input/Output) The mean outflow for period. If gates
C                   is less than zero, then determine EOP and mean level
C                   for Lake Superior based on MeanOutflow.  If gates is
C                   positive, determine EOP and mean level based on gate
C                   equations, and return average outflow in MeanOutflow
C
      DOUBLE PRECISION EOPLevel, MeanLevel
C                   (Output) The End-Of-Period and period-average level.
C            
C
C
C     Declarations regarding local variables or functions:
C
C
      INTEGER       I, G, MaxIterations
      DOUBLE PRECISION TmpF, TmpA
C                   Scratch
C
      DOUBLE PRECISION NonCompWorksFlow, Round, GetArea
C                   Function returns sidechannel flow given lake levels,
C                   function to round, and function to return area of
C                   lake given its elevation.
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
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
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
C     XX This hardcoded limit should be changed if have the chance to 
C     XX tweak Plan 1977A.  Better to iterate to a tolerance.  That 
C     XX itself is easier to do if no rounding.
      DATA MaxIterations /4/
C
C     Init MeanLevel
      MeanLevel = BOPLevel
C
C     This initialization doesn't affect the calculation, but it avoids
C     confusion in the first print statement.
      EOPLevel = -9999.0D0
C
C     Compute first approximation
C
      IF(Gates .GT. 0.0) THEN 
C        Save index to SGX, SGY, SGZ in G 
         G = NINT(Gates + 0.51)
         IF(SGY(G) * MeanLevel - SGZ(G) .LT. 0.0D0) THEN
            CALL ErrMess(1062,'P77ASolveSup'//LLT//'   ',LOGCDE)
            WRITE( CTMP, '(A,F7.3,A,F6.2,A1)') ' Superior Level=', 
     &         MeanLevel, '   Gates=', Gates, LLT
            CALL LogIt(LOGCDE,1,0)
            CALL SAVETOD
            CALL ADIOS
         ENDIF
C        original
         MeanOutflow = Round(NonCompWorksFlow(MeanLevel, 0.D0, MonthNum)
     &      +SGX(G) * (SGY(G) * MeanLevel - SGZ(G)) ** 1.5D0, FlowRound)
      ENDIF
C
C
      IF( SupVerbosity .GE. 8) THEN
         TmpA = GetArea(1,MeanLevel)
         Write(CTMP, '(A,F10.5,A,D21.14,A,F13.2,A,F9.2,A,F6.2,A1)')
     &      ' Entering P77ASolveSup, BOP Level= ', BOPLevel, '  Area= ',
     &      TmpA, ' dT= ',dT, '   NTS=', NTS, '   Gates=',  Gates, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C      
C     Repeat steps 6 - 10
      DO 10, I = 1, MaxIterations
         IF( SupVerbosity .GE. 8) THEN
            Write(CTMP, '(2(A,I2),3(A,F11.5),A1)')' Iteration ', 
     &          I, ' of ',MaxIterations,
     &          '  EOPLevel=', EOPLevel, '  MeanLevel=', MeanLevel, 
     &          ' MeanOutflow=', MeanOutflow, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
         EOPLevel= Round(dT * (NTS - MeanOutflow)/(GetArea(1,MeanLevel))
     &      + BOPLevel, LevelRound)
         MeanLevel = Round((BOPLevel + EOPLevel) * 0.5D0, LevelRound)
         IF( SupVerbosity .GE. 8) THEN
C           Write(CTMP, '(2(A,F15.8),A,F16.10,A1)')
C    &          '      NTS - outflow: ', NTS - MeanOutflow, 
C    &          ' ( ', Round((NTS - MeanOutflow), FlowRound),
C    &          ')  dT * (NTS - MeanOutflow)/Area',
C    &          Round(dT *(NTS - MeanOutflow)/(GetArea(1,MeanLevel)), 
C    &          LevelRound),LLT
C           CALL LogIt(LOGD,0,0)
            TmpA = GetArea(1,MeanLevel)
            Write(CTMP,'(A,F16.10,A)')'      Unrounded stage change= ',
     &         dT * (NTS - MeanOutflow)/(TmpA),LLT
            CALL LogIt(LOGD,0,0)
            Write(CTMP,'(A,F16.10,A)')'      New EOPLevel=',EOPLevel,LLT
            CALL LogIt(LOGD,0,0)
            Write(CTMP, '(A,F16.10,A1)')
     &         '      New MeanLevel=', MeanLevel, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
C        Use Pt Iroquois Eq's because that's what Plan 1977A does.
C        XX Maybe we can add an additional argument to this function
C        XX for optionally using SW Pier, is calling this function 
C        XX for calculations outside Plan1977A.
C 
C
         IF(Gates .GT. 0.0)MeanOutflow= Round(
     &      NonCompWorksFlow(MeanLevel, 0.D0, MonthNum) +
     &      SGX(G) * (SGY(G) * MeanLevel - SGZ(G)) ** 1.5D0, FlowRound)
C
         IF( SupVerbosity .GE. 8) THEN
            IF(Gates .GT. 0.6) THEN
               WRITE(CWRK, '(A,3F9.4,A,F12.7,A1)') 
     &           ' X,Y,Z= ',SGX(G),SGY(G),SGZ(G),  
     &          ' SGX*(SGY*MeanLevel - SGZ)= ',
     &          SGX(G) * (SGY(G) * MeanLevel-SGZ(G)), LLT
            ELSE IF(Gates .LT. 0.1) THEN
               WRITE(CWRK, '(A,A1)') ' ', LLT
            ELSE 
               WRITE(CWRK, '(A,A1)') ', 1/2 gate open, ', LLT
            ENDIF
C           Mysterious behavior here.  Embedding the call to NonCompWorksFlow
C           directly into the following WRITE corrupted the string.  
C           Assigning to a temporary value first works OK.  Is there some
C           missing typecasting?  A compiler bug (MS Powerstation)?  Or is 
C           this just the manifestation of a subtle bug somewhere else?
C           XX Try getting rid of TmpF and embedding call in WRITE again.
            TmpF = NonCompWorksFlow(MeanLevel, 0.D0, MonthNum)
            WRITE(CTMP, '(A,F10.5,A,A,F10.5,A1)') 
     &         '      NonCompWorksFlow= ', TmpF, 
C     &      '      NonCompWorksFlow= ',NonCompWorksFlow(MeanLevel,0.D0), 
     &         CWRK(1:INDEX(CWRK,LLT)-1),        
     &         ' New Mean Outflow=', MeanOutflow, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
   10 CONTINUE
C      
      IF( SupVerbosity .GE. 6) THEN
         WRITE( CTMP, 6556) Gates, MeanOutflow, MeanLevel, LLT
         CALL LogIt(LOGD,0,0)
         IF( SupVerbosity .GE. 8) THEN
            WRITE( CTMP, 6555) BOPLevel, EOPLevel, NTS, dT, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
      ENDIF
 6555 FORMAT(10X,'Superior BOP level: ', F8.4, ' Superior EOP level: ',
     &    F8.4, ' net total supply: ', F9.3, ' seconds in period: ', 
     &    F9.1, A1)
 6556 FORMAT(10X, 'Gates open: ', F6.2,
     &   ' mean outflow:', F8.3, ' mean level: ', F8.4, A1)
C  
      END
C
C----&------------------------------------------------------------------
C
C
      SUBROUTINE Plan1977A(Month, Year, PreviousFlow, PreviousGates,
     &   SupNTSfromTS, SupBOM, MHBOM, StClairBOM, ErieBOM, 
     &   Plan1977AFlow, GatesOpen)
C      
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER          Month, Year
C                      (input) The month number and year to compute 
C                      outflow for.
C
      DOUBLE PRECISION PreviousFlow
      REAL             PreviousGates
C                      (input) The number of gates open last month, and
C                      the outflow.
C
      INTEGER          SupNTSfromTS
C                      (input) Net total supply for the month, from the
C                      time-series input.
C    
      DOUBLE PRECISION SupBOM, MHBOM, StClairBOM, ErieBOM
C                      (input) Beginning-of-Month lake levels 
C    
      DOUBLE PRECISION Plan1977AFlow 
C                      (output) The flow for Month
C
      REAL             GatesOpen
C                      (output) The number of gates open in Month.
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER          I, L
C                      Scratch variables
C
      CHARACTER*(*)    IDString
C                      Used as prefix on lines of detailed.log 
C
      INTEGER          ForecastPeriod, ForecastMonth(2)
C                      Identify the month for which forecasted 
C                      outflow is being computed.  ForecastPeriod
C                      refers to sequence in P77AMonthsToAve, while
C                      ForecastMonth(1) contains the month number from
C                      1-12, ForecastMonth(2) contains the year.  So if
C                      computing Plan 1977A for October 1998, with 
C                      5-month averaging, values get updated 5 times: 
C                       
C                          ForecastPeriod   ForecastMonth 
C                               1             10   1998
C                               2             11   1998
C                               3             12   1998
C                               4              1   1999
C                               5              2   1999
C
C      DOUBLE PRECISION MHNTS 
C                      The Net Total Supply to use for Lakes 
C                      Michigan-Huron (considering trigger levels, 
C                      and diverions)
C
      DOUBLE PRECISION ForecastSPEOP, ForecastSPBOP
      DOUBLE PRECISION ForecastMHEOP, ForecastMHBOP
      DOUBLE PRECISION ForecastSCEOP, ForecastSCBOP
      DOUBLE PRECISION ForecastEREOP, ForecastERBOP
      DOUBLE PRECISION ForecastSTMFlow, ForecastPrevFlow,ForecastBalance
C  ForecastSPMLV, ForecastMHMLV, ForecastERMLV, ForecastSCMLV,
C   ForecastSCRFlow, ForecastDETFlow, ForecastNIAFlow,
      DOUBLE PRECISION SecsInPeriod
C                      The number of seconds in routing period (month).
C
C      DOUBLE PRECISION GETAREA
C                      Function that returns surface area for a lake
C
C      DOUBLE PRECISION OUTFLOW
C                      Function that returns connecting channel flow, 
C                      given the parameters for the discharge equation.
C
      DOUBLE PRECISION Round
C                      Function for rounding.
C
C
C
C     Declarations regarding COMMON block variables:
C
C
      CHARACTER     SupPreProject*6, BalancingFlow*8, SupGateFlow*6,
     &              SupPlanFlow*6, SupBOMLevel*8, SWPierLevel*8,
C    &              SupNavFlow*6, SupPowerFlow*6, 
     &              SupSideChanFlow*6, SPRegCodes*40, P77AHeader*120,
     &              P2012Header
C                   Lake Superior regulation output strings
C
      INTEGER P77AWinterMax, P77AMaxChange, P77AMinimumFlow(12)
C     Plan 1977A parameters for bounds-checking
C
      DOUBLE PRECISION P77ASupTargetLev(12), P77ASupRangeFact(12),
     &        P77AMHTargetLev(12), P77AMHRangeFact(12), P77AAValue
      INTEGER P77ABaseFlow(12)
C     Plan 1977A parameters for Balancing Equation
C
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
C     Plan 1977A routing parameters
C
      INTEGER P77AMonthsToAve(12)
      LOGICAL P77AUniformMonth
C
C
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
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
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
      COMMON /MONAME/   MonthAbbrev
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /SUPOUT/   SupPreProject, BalancingFlow, SupGateFlow,
     &                  SupPlanFlow, SupBOMLevel, SWPierLevel, 
     &                  SupSideChanFlow, SPRegCodes, P77AHeader,
     &                  P2012Header
C
      PARAMETER (IDString='    Plan 1977A==> ')
C
      DATA SecsInPeriod     /2629746.D0/
C
C
      SPRegCodes= '                                        '
C
C     Log the data entering Plan 1977A    
C
      IF( SupVerbosity .GE. 5) THEN
         WRITE(CTMP, 6523) IDString, MonthAbbrev(Month), Year, LLT
         CALL LogIt(LOGD,1,0)
         WRITE(CTMP, 6524) IDString,SupBOM,MHBOM,StClairBOM,ErieBOM, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6525) IDString, PreviousFlow, PreviousGates, 
     &      SupNTSfromTS, LLT
         CALL LogIt(LOGD,0,0)
         IF( SupVerbosity .GE. 8) THEN
            WRITE(CTMP, 6526) IDString, P77AMonthsToAve(Month), LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
      ENDIF
C
      Plan1977AFlow = 0.0D0
C
C     Initialize BOP/EOP levels, and month/year counters
C
      ForecastSPEOP = Round(SupBOM,P77AMonthRound)
      ForecastMHEOP = Round(MHBOM,P77AMonthRound)
      ForecastSCEOP = Round(StClairBOM,P77AMonthRound)
      ForecastEREOP = Round(ErieBOM,P77AMonthRound)
      ForecastSTMFlow = PreviousFlow
      ForecastMonth(1) = Month - 1
      ForecastMonth(2) = Year
C
C
C     Write interesting data to character buffers for later output
C
C
      WRITE(BalancingFlow, '(I8)') IDNINT(Round( P77ABaseFlow(Month) + 
     &      P77AAValue * 
     &         (ForecastSPEOP - 
     &             (P77ASupTargetLev(Month) + 
     &                (ForecastMHEOP - P77AMHTargetLev(Month))
     &                * P77ASupRangeFact(Month) / 
     &                P77AMHRangeFact(Month))), 0))
C
C
C     Loop through each month of the Plan 1977A forecast    
C
      DO 5555, ForecastPeriod = 1, P77AMonthsToAve(Month)
C
C
C        Increment month/year counters
C
         ForecastMonth(1) = ForecastMonth(1) + 1
         IF( ForecastMonth(1) .GT. 12) THEN
            ForecastMonth(1) = 1
            ForecastMonth(2) = ForecastMonth(2) + 1
         ENDIF
C
C
C        Determine length of timestep (unless uniform month)
C
         IF( .NOT. P77AUniformMonth )THEN
C           Get number of days in month in L (I is unused)
            CALL NumDays(ForecastMonth(2), ForecastMonth(1), I, L)
            SecsInPeriod = L * 86400
         ENDIF
C
C
C        Roll over the EOP/BOP values for this forecast period
C
         ForecastSPBOP = Round(ForecastSPEOP, P77AMonthRound)
         ForecastMHBOP = Round(ForecastMHEOP, P77AMonthRound)
         ForecastSCBOP = Round(ForecastSCEOP, P77AMonthRound)
         ForecastERBOP = Round(ForecastEREOP, P77AMonthRound) 
         ForecastPrevFlow = ForecastSTMFlow
C
         IF( SupVerbosity .GE. 5) THEN
            WRITE(CTMP, 6527) IDString, ForecastPeriod, 
     &         P77AMonthsToAve(Month), MonthAbbrev(ForecastMonth(1)),
     &         ForecastMonth(2), LLT
            CALL LogIt(LOGD,1,0)
         ENDIF
C
C
C        Compute balancing equation flow:
C
C
         ForecastBalance = Round( P77ABaseFlow(ForecastMonth(1)) + 
     &      P77AAValue * 
     &         (ForecastSPBOP - 
     &             (P77ASupTargetLev(ForecastMonth(1)) + 
     &                (ForecastMHBOP -P77AMHTargetLev(ForecastMonth(1)))
     &                * P77ASupRangeFact(ForecastMonth(1)) / 
     &                P77AMHRangeFact(ForecastMonth(1)))),
     &      P77AFlowRound)
C
         IF( SupVerbosity .EQ. 9) THEN
            WRITE(CTMP, 2456) IDString, ForecastMonth(1), P77AAValue, 
     &         ForecastMHBOP, ForecastSPBOP, 
     &         P77ABaseFlow(ForecastMonth(1)), 
     &         P77AMHTargetLev(ForecastMonth(1)), 
     &         P77AMHRangeFact(ForecastMonth(1)),
     &         P77ASupTargetLev(ForecastMonth(1)), 
     &         P77ASupRangeFact(ForecastMonth(1)), LLT
            CALL LogIt(LOGD,0,0)
 2456       FORMAT(A,' Month=', I2, ' A=',F9.2,' MHBOP=', F9.4, 
     &         ' SupBOP=', F9.4, ' Qavg=', I7, ' MH stats=', F9.4,
     &         F7.4, ' Sup stats=', F9.4, F7.4, A1)
         ENDIF
C 
         ForecastSTMFlow = ForecastBalance
C
         SPRegCodes(1:3) = 'BE '
C
         CALL P77ALimitations( SecsInPeriod, ForecastPrevFlow, -9999,
     &      ForecastSPBOP, ForecastMHBOP, ForecastSCBOP, ForecastERBOP,
     &      ForecastSPEOP, ForecastMHEOP, ForecastSCEOP, ForecastEREOP,
     &      ForecastSTMFlow, GatesOpen, ForecastMonth(1), 
     &      SPRegCodes, IDString//' (fcast flow adjust) '//LLT)
C
C
         IF( SupVerbosity .GE. 5) THEN
            WRITE(CTMP, 6528) IDString, ForecastBalance, LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP, 6522) IDString, 'forecast', ForecastSTMFlow, 
     &         GatesOpen, LLT 
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP, 6521) IDString, ForecastSPBOP, ForecastMHBOP, 
     &         ForecastSCBOP, ForecastERBOP, LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP, 6520) IDString, ForecastSPEOP, ForecastMHEOP, 
     &         ForecastSCEOP, ForecastEREOP, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF 
C
C
C
C
C
C
C
C        Accumulate forecast flows in Plan1977AFlow for later averaging
 5555    Plan1977AFlow = Plan1977AFlow + ForecastSTMFlow
C
C
C     Compute mean of forecast outflows
      Plan1977AFlow = Round(Plan1977AFlow / (ForecastPeriod - 1), 
     &   P77AFlowRound)
C
C
C
      IF( SupVerbosity .GE. 5) THEN
         WRITE(CTMP, 6529) IDString, Plan1977AFlow, LLT
         CALL LogIt(LOGD,1,0)
      ENDIF
C  
C
      IF( .NOT. P77AUniformMonth )THEN
C        Get number of days in month in L (I is unused)
         CALL NumDays(Year, Month, I, L)
         SecsInPeriod = L * 86400
      ENDIF
C
C
C     Do final adjustments of averaged forecast flow.  
C     ForecastSPEOP and ForecastSPMLV are not used (Superior actually takes
C     place in SuperiorRoute).  
C     
C
      SPRegCodes(1:3) = 'AF '
C
      CALL P77ALimitations( SecsInPeriod, PreviousFlow, SupNTSfromTS,
C     XXXX Another dumb hack for backward compatibility.  Round levels
C     XXXX here to 2 places instead of usual four.  Would like to change
C     XXXX this!  Consider it a bug in Carl's model
     &   Round(SupBOM, P77AMonthRound), Round(MHBOM, P77AMonthRound), 
     &   Round(StClairBOM,P77AMonthRound),Round(ErieBOM,P77AMonthRound),
     &   ForecastSPEOP, ForecastMHEOP, ForecastSCEOP, ForecastEREOP,
     &   Plan1977AFlow, GatesOpen, Month, 
     &   SPRegCodes, IDString//' (final flow adjust) '//LLT)
C
C
      WRITE(SupPlanFlow,'(I6)') IDNINT(Plan1977AFlow) 
C
      IF( SupVerbosity .GE. 5) THEN
         WRITE(CTMP, 6522) IDString, 'Plan 1977A', Plan1977AFlow, 
     &      GatesOpen, LLT 
         CALL LogIt(LOGD,0,0)
C        WRITE(CTMP, 6521) IDString, ForecastSPBOP, ForecastMHBOP, 
C    &      ForecastSCBOP, ForecastERBOP, LLT
C        CALL LogIt(LOGD,0,0)
C        WRITE(CTMP, 6520) IDString, ForecastSPEOP, ForecastMHEOP, 
C    &      ForecastSCEOP, ForecastEREOP, LLT
C        CALL LogIt(LOGD,0,1)
      ENDIF 
C
C
C
 6522 FORMAT( A,' Computed St. Marys ', A,' flow of ', F6.1, 
     &   ' m3s (',  F5.2, ' gates)', A1)
 6521 FORMAT( A,' Beginning-of-Period levels:', 4F11.5, A1)
 6520 FORMAT( A,'       End-of-Period levels:', 4F11.5, A1)
 6523 FORMAT( A,'Entering Plan1977A to compute nominal outflow for ', 
     &   A3, I5, '.', A1) 
 6524 FORMAT( A,'Beginning of month levels:', 4F16.11, A1)
 6525 FORMAT( A,'Previous month outflow=', F20.3, ', gates open=', 
     &   F8.2, '   Sup NTS=',I13, A1) 
 6526 FORMAT( A,'Averaging over ', I2, ' forecasted months.', A1)
 6527 FORMAT( A,' Calculating forecasted outflow for month ', I2,
     &   ' of ', I2, ' (', A3, I5, ').', A1)
 6528 FORMAT( A,' Balancing equation flow: ', F12.3, A1 )
 6529 FORMAT( A,' Average of forecasted outflows: ', F12.3, A1 )
C 
      END
C
C
C    
C
C
C----&------------------------------------------------------------------
C
C
C
      SUBROUTINE P77AMidLakeRoute (SecondsInPeriod, MHNTS, SCNTS, ERNTS, 
     &   MHBOP, SCBOP, ERBOP, SCRRetard, DETRetard, NIARetard, 
     &   Frnd, Lrnd, Mrnd, NumInc, MHEOP, SCEOP, EREOP, MHMLV)
C      
C     This function does middle lakes rounding as described by Plan
C     1977A.  It is a simpler, less robust method than IterativeRoute,
C     but doesn't rely as much on COMMON block data.  Is used for both
C     Plan 1977A routing, and Criterion checks.  Only the Plan 1977A
C     discharge equation constants are used from the Plan 1977A COMMON
C     block; all else must be specified as parameters.  This is a 
C     confusing approach, but that's because the existing practice for
C     Criterion checks has very messy logic.  Believe it or not, this
C     is actually a much clearer computer implementation than existing 
C     models.
C     
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      DOUBLE PRECISION SecondsInPeriod     
C                      (input) The length of the routing interval.
C
      DOUBLE PRECISION MHBOP, SCBOP, ERBOP
C                      (input) Beginning-of-Period lake levels 
C    
      INTEGER          MHNTS, SCNTS, ERNTS
C                      (Input) Net Total Supply to use for the lakes 
C                      (including diverions).
C
      INTEGER          SCRRetard, DETRetard, NIARetard
C                      (input) Retardation amounts to use
C    
      INTEGER          Frnd, Lrnd, Mrnd, NumInc
C                      (Input) Rounding levels for flows and levels,
C                      expessed as powers of ten (so that -2 means cm),
C                      and end-of-month level rounding factor,
C                      and the number of increments to use in the month.
C
      DOUBLE PRECISION MHEOP, SCEOP, EREOP, MHMLV
C                      (output) End-of-Period lake levels, and average
C                      MH level for the period.
C
C                       
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER          I, K
C                      Scratch variables
C
      DOUBLE PRECISION SCRIncFlow, DETIncFlow, NIAIncFlow, LMHIncLevel, 
     &                 LSCIncLevel, LERIncLevel, MHBOI, SCBOI, ERBOI
C                      Intermediate routing terms (mean levels and
C                      flows for an increment, beginning-of-increment
C                      levels)
C
      CHARACTER*(*)    Prefix
C                      Used as prefix on lines of detailed.log 
C
      DOUBLE PRECISION Round
C                      Function for rounding.
C
      DOUBLE PRECISION GETAREA
C                      Function that returns surface area for a lake
C
      DOUBLE PRECISION OUTFLOW
C                      Function that returns connecting channel flow, 
C                      given the parameters for the discharge equation.
C
C
C     Declarations regarding COMMON block variables:
C
      INTEGER P77AWinterMax, P77AMaxChange, P77AMinimumFlow(12)
C     Plan 1977A parameters for bounds-checking
C
      DOUBLE PRECISION P77ASupTargetLev(12), P77ASupRangeFact(12),
     &        P77AMHTargetLev(12), P77AMHRangeFact(12), P77AAValue
      INTEGER P77ABaseFlow(12)
C     Plan 1977A parameters for Balancing Equation
C
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
C     Plan 1977A routing parameters
C
      INTEGER P77AMonthsToAve(12)
      LOGICAL P77AUniformMonth
C
C
C     
C
C     Declarations regarding COMMON block variables:
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
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
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
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
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
C
      PARAMETER (Prefix = 'P77A MidLake Routing')
C
      IF( SupVerbosity .GE. 7) THEN
         WRITE(CTMP,11116) Prefix, ' MHBOP=',MHBOP, 
     &      ' SCBOP=', SCBOP, ' ERBOP=', ERBOP,LLT
         CALL LogIt(LOGD,1,0)
         WRITE(CTMP,11117) Prefix, ' MH NTS=', MHNTS, 
     &      ' StC NTS=', SCNTS, ' Erie NTS=', ERNTS, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,11117) Prefix, ' IceMH=', SCRRetard,
     &      ' IceSC=', DETRetard, ' IceER=', NIARetard, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
C
      SCRIncFlow = 0.D0
      DETIncFlow = SCRIncFlow 
      NIAIncFlow = SCRIncFlow 
C
      MHEOP = MHBOP
      SCEOP = SCBOP
      EREOP = ERBOP
C
      MHMLV = 0.D0
C
      DO 5210,  I = 1, NumInc
C
C
         MHBOI = MHEOP
         SCBOI = SCEOP
         ERBOI = EREOP
C
C
C        Approximate the mean level for the increment as BOI level
         LMHIncLevel = MHBOI
         LSCIncLevel = SCBOI 
         LERIncLevel = ERBOI 
C
C        Iterate five times (maybe in future remove rounding, and
C        iterate to a tolerance?)
         DO 5230, K = 1, 5            
C
            SCRIncFlow = Round( OUTFLOW( LMHIncLevel, LSCIncLevel, 
     &         P77AMHK, P77AMHYM, P77AMHA,P77AMHB, P77AMHWT, P77AMHC) 
     &         - SCRRetard, Frnd ) 
            DETIncFlow = Round( OUTFLOW( LSCIncLevel, LERIncLevel, 
     &         P77ASCK, P77ASCYM, P77ASCA, P77ASCB, P77ASCWT, P77ASCC) 
     &         - DETRetard, Frnd ) 
            NIAIncFlow = Round( OUTFLOW( LERIncLevel, P77ACGIP, P77AERK, 
     &         P77AERYM, P77AERA, P77AERB, P77AERWT, P77AERC)
     &         - NIARetard, Frnd ) 
C
C           Compute EOP levels for the iteration
C            MHEOP = MHBOI + (MHNTS - SCRIncFlow) * 
C     &         SecondsInPeriod/(P77AMidLakesInc*GetArea(2,LMHIncLevel)) 
C            WRITE(CTMP,'(22X,2(3X,A,F9.2),A,F14.12,A1)')'MH NTS=',MHNTS, 
C     &         'SCRIncFlow=',SCRIncFlow, '  dT/A=', SecondsInPeriod
C     &         /(P77AMidLakesInc*GetArea(2,LMHIncLevel)), LLT
C            CALL LogIt(LOGD,0,0)
C            SCEOP = SCBOI + 
C     &         (P77AStC50NBS(Mon) + SCRIncFlow - DETIncFlow) * 
C     &         SecondsInPeriod/(P77AMidLakesInc*GetArea(5,LSCIncLevel))
C            EREOP = ERBOI + (P77AErie50NBS(Mon) 
C     &         + DETIncFlow - NIAIncFlow + P77AWellandFlow(Mon)) * 
C     &         SecondsInPeriod/(P77AMidLakesInc*GetArea(6,LERIncLevel))
C
            MHEOP = Round( MHBOI + (MHNTS - SCRIncFlow) * 
     &         SecondsInPeriod/(NumInc*GetArea(2,LMHIncLevel)), Lrnd )
            SCEOP = Round( SCBOI + (SCNTS + SCRIncFlow - DETIncFlow) *  
     &         SecondsInPeriod/(NumInc*GetArea(5,LSCIncLevel)), Lrnd)
            EREOP = Round( ERBOI + (ERNTS + DETIncFlow - NIAIncFlow) * 
     &         SecondsInPeriod/(NumInc*GetArea(6,LERIncLevel)), Lrnd )
C
C           Compute mean levels for the iteration
            LMHIncLevel = ( MHBOI + MHEOP ) / 2
            LSCIncLevel = ( SCBOI + SCEOP ) / 2
            LERIncLevel = ( ERBOI + EREOP ) / 2 
 5230       CONTINUE
C
C 
C 
C
         IF( SupVerbosity .GE. 9) THEN
            WRITE(CTMP,11111) Prefix, I, ' OutMH=',SCRIncFlow, 
     &         ' OutSC=', DETIncFlow,' OutER=', NIAIncFlow, LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,11111) Prefix, I, ' MHBOI=',MHBOI, 
     &         ' SCBOI=', SCBOI, ' ERBOI=', ERBOI,LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,11111) Prefix, I, ' MHEOI=',MHEOP, 
     &         ' SCEOI=', SCEOP, ' EREOI=', EREOP,LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,11111) Prefix, I, ' MHMLV=',LMHIncLevel, 
     &         ' SCMLV=',LSCIncLevel,' ERMLV=', LERIncLevel, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
11111    FORMAT(A,' Inc ', I2, 3(A,F15.8),A1)
11117    FORMAT(A,3(A,I8),A1)
11116    FORMAT(A,3(A,F15.8),A1)
C
C
 5210    MHMLV = MHMLV + LMHIncLevel
C
      MHEOP = Round( MHEOP, Mrnd )
      SCEOP = Round( SCEOP, Mrnd )
      EREOP = Round( EREOP, Mrnd )
C
      MHMLV = Round( MHMLV / NumInc, Mrnd )
C
      IF( SupVerbosity .GE. 7) THEN
         WRITE(CTMP,11116) Prefix, ' MHEOP=',MHEOP, 
     &      ' SCEOP=', SCEOP, ' EREOP=', EREOP,LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,'(A,A,F15.8,A1)') Prefix, ' MHMLV=', MHMLV, LLT 
         CALL LogIt(LOGD,0,1)
      ENDIF
C
C
      END
C
C---------------------------------------------------------------------------------------

      SUBROUTINE P2012MidLakeRoute (SecondsInPeriod,MHNTS,SCNTS,ERNTS, 
     &   MHBOP, SCBOP, ERBOP, SCRRetard, DETRetard, NIARetard, 
     &   Frnd, Lrnd, Mrnd, NumInc, MHEOP, SCEOP, EREOP, MHMLV)
C      
C     This function does middle lakes borrowed from P77AMidLakeRoute but using the updated lake to lake equations
C     
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      DOUBLE PRECISION SecondsInPeriod     
C                      (input) The length of the routing interval.
C
      DOUBLE PRECISION MHBOP, SCBOP, ERBOP
C                      (input) Beginning-of-Period lake levels 
C    
      INTEGER          MHNTS, SCNTS, ERNTS
C                      (Input) Net Total Supply to use for the lakes 
C                      (including diverions).
C
      INTEGER          SCRRetard, DETRetard, NIARetard
C                      (input) Retardation amounts to use
C    
      INTEGER          Frnd, Lrnd, Mrnd, NumInc
C                      (Input) Rounding levels for flows and levels,
C                      expessed as powers of ten (so that -2 means cm),
C                      and end-of-month level rounding factor,
C                      and the number of increments to use in the month.
C
      DOUBLE PRECISION MHEOP, SCEOP, EREOP, MHMLV
C                      (output) End-of-Period lake levels, and average
C                      MH level for the period.
C
      DOUBLE PRECISION MHK,MHYM,MHA,MHB,MHWT,MHC
      DOUBLE PRECISION SCK,SCYM,SCA,SCB,SCWT,SCC
      DOUBLE PRECISION ERK,ERYM,ERA,ERB,ERWT,ERC
C     Declarations regarding local variables or functions:
C
      INTEGER          I, K
C                      Scratch variables
C
      DOUBLE PRECISION SCRIncFlow, DETIncFlow, NIAIncFlow, LMHIncLevel, 
     &                 LSCIncLevel, LERIncLevel, MHBOI, SCBOI, ERBOI
C                      Intermediate routing terms (mean levels and
C                      flows for an increment, beginning-of-increment
C                      levels)
C
      CHARACTER*(*)    Prefix
C                      Used as prefix on lines of detailed.log 
C
      DOUBLE PRECISION Round
C                      Function for rounding.
C
      DOUBLE PRECISION GETAREA
C                      Function that returns surface area for a lake
C
      DOUBLE PRECISION OUTFLOW
C                      Function that returns connecting channel flow, 
C                      given the parameters for the discharge equation.
C
C     Declarations regarding COMMON block variables:
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
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
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK

      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
C
      COMMON /STAGE/    MHK,SCK,ERK,MHYM,SCYM,ERYM,MHA,SCA,ERA,
     &                  MHB,SCB,ERB,MHWT,SCWT,ERWT,MHC,SCC,ERC
      PARAMETER (Prefix = 'P2012 MidLake Routing')
C
      IF( SupVerbosity .GE. 7) THEN
         WRITE(CTMP,11116) Prefix, ' MHBOP=',MHBOP, 
     &      ' SCBOP=', SCBOP, ' ERBOP=', ERBOP,LLT
         CALL LogIt(LOGD,1,0)
         WRITE(CTMP,11117) Prefix, ' MH NTS=', MHNTS, 
     &      ' StC NTS=', SCNTS, ' Erie NTS=', ERNTS, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,11117) Prefix, ' IceMH=', SCRRetard,
     &      ' IceSC=', DETRetard, ' IceER=', NIARetard, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
C
      SCRIncFlow = 0.D0
      DETIncFlow = SCRIncFlow 
      NIAIncFlow = SCRIncFlow 
C
      MHEOP = MHBOP
      SCEOP = SCBOP
      EREOP = ERBOP
C
      MHMLV = 0.D0
C
      DO 5210,  I = 1, NumInc
C
C
         MHBOI = MHEOP
         SCBOI = SCEOP
         ERBOI = EREOP
C
C
C        Approximate the mean level for the increment as BOI level
         LMHIncLevel = MHBOI
         LSCIncLevel = SCBOI 
         LERIncLevel = ERBOI 
C
C        Iterate five times (maybe in future remove rounding, and
C        iterate to a tolerance?)
         DO 5230, K = 1, 5            
C
            SCRIncFlow = Round( OUTFLOW( LMHIncLevel, LSCIncLevel, 
     &         MHK, MHYM, MHA,MHB, MHWT, MHC) 
     &         - SCRRetard, Frnd ) 
            DETIncFlow = Round( OUTFLOW( LSCIncLevel, LERIncLevel, 
     &         SCK, SCYM, SCA, SCB, SCWT, SCC) 
     &         - DETRetard, Frnd ) 
            NIAIncFlow = Round( OUTFLOW( LERIncLevel, 0.D0, ERK, 
     &         ERYM, ERA, ERB, ERWT, ERC)
     &         - NIARetard, Frnd ) 
     
C
C           Compute EOP levels for the iteration
C            MHEOP = MHBOI + (MHNTS - SCRIncFlow) * 
C     &         SecondsInPeriod/(P77AMidLakesInc*GetArea(2,LMHIncLevel)) 
C            WRITE(CTMP,'(22X,2(3X,A,F9.2),A,F14.12,A1)')'MH NTS=',MHNTS, 
C     &         'SCRIncFlow=',SCRIncFlow, '  dT/A=', SecondsInPeriod
C     &         /(P77AMidLakesInc*GetArea(2,LMHIncLevel)), LLT
C            CALL LogIt(LOGD,0,0)
C            SCEOP = SCBOI + 
C     &         (P77AStC50NBS(Mon) + SCRIncFlow - DETIncFlow) * 
C     &         SecondsInPeriod/(P77AMidLakesInc*GetArea(5,LSCIncLevel))
C            EREOP = ERBOI + (P77AErie50NBS(Mon) 
C     &         + DETIncFlow - NIAIncFlow + P77AWellandFlow(Mon)) * 
C     &         SecondsInPeriod/(P77AMidLakesInc*GetArea(6,LERIncLevel))
C
            MHEOP = Round( MHBOI + (MHNTS - SCRIncFlow) * 
     &         SecondsInPeriod/(NumInc*GetArea(2,LMHIncLevel)), Lrnd )
            SCEOP = Round( SCBOI + (SCNTS + SCRIncFlow - DETIncFlow) *  
     &         SecondsInPeriod/(NumInc*GetArea(5,LSCIncLevel)), Lrnd)
            EREOP = Round( ERBOI + (ERNTS + DETIncFlow - NIAIncFlow) * 
     &         SecondsInPeriod/(NumInc*GetArea(6,LERIncLevel)), Lrnd )
C
C           Compute mean levels for the iteration
            LMHIncLevel = ( MHBOI + MHEOP ) / 2
            LSCIncLevel = ( SCBOI + SCEOP ) / 2
            LERIncLevel = ( ERBOI + EREOP ) / 2 
 5230       CONTINUE
C
C 
C 
C
         IF( SupVerbosity .GE. 7) THEN
            WRITE(CTMP,11111) Prefix, I, ' OutMH=',SCRIncFlow, 
     &         ' OutSC=', DETIncFlow,' OutER=', NIAIncFlow, LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,11111) Prefix, I, ' MHBOI=',MHBOI, 
     &         ' SCBOI=', SCBOI, ' ERBOI=', ERBOI,LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,11111) Prefix, I, ' MHEOI=',MHEOP, 
     &         ' SCEOI=', SCEOP, ' EREOI=', EREOP,LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,11111) Prefix, I, ' MHMLV=',LMHIncLevel, 
     &         ' SCMLV=',LSCIncLevel,' ERMLV=', LERIncLevel, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
11111    FORMAT(A,' Inc ', I2, 3(A,F15.8),A1)
11117    FORMAT(A,3(A,I8),A1)
11116    FORMAT(A,3(A,F15.8),A1)
C
C
 5210    MHMLV = MHMLV + LMHIncLevel
C
      MHEOP = Round( MHEOP, Mrnd )
      SCEOP = Round( SCEOP, Mrnd )
      EREOP = Round( EREOP, Mrnd )
C
      MHMLV = Round( MHMLV / NumInc, Mrnd )
C
      IF( SupVerbosity .GE. 7) THEN
         WRITE(CTMP,11116) Prefix, ' MHEOP=',MHEOP, 
     &      ' SCEOP=', SCEOP, ' EREOP=', EREOP,LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,'(A,A,F15.8,A1)') Prefix, ' MHMLV=', MHMLV, LLT 
         CALL LogIt(LOGD,0,1)
      ENDIF
C
C
      END
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
      SUBROUTINE P77ALimitations( SecsInPeriod, PreviousFlow, InSupNTS,
     &   SuperiorBOPLevel,MicHuronBOPLevel,StClairBOPLevel,ErieBOPLevel,  
     &   SuperiorEOPLevel,MicHuronEOPLevel,StClairEOPLevel,ErieEOPLevel,
     &   StMarysFlow, Gates, Month, AdjustmentCodes, OutputPrefix)
C
C     This function takes a preliminary flow, from either the balancing
C     equation or the average of forecasted plan flows, and applies the 
C     various outflow limitations described in Plan 1977A.  These are:
C
C        Requirement A:  Maximum outflow May-Nov = 2320 m3s + 16 gates
C        Requirement B:  Maximum outflow Dec-Apr = 2410 m3s
C        Requirement C:  Minimum outflow = 1560 m3s
C        Requirement D:  (not applicable in this context)
C        Requirement E:  (not applicable in this context)
C        Requirement F:  (not applicable in this context)
C        Requirement G:  (not applicable in this context)
C        Requirement H:  i don't think there is one
C        Requirement I:  Maximum change between months = 850 m3s
C        Requirement J:  Minimum gates open = 0.5 m3s
C        Requirement K:  Maximum outflow Oct = 4110 m3s, Nov = 3260 m3s
C
C
C
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER          Month
C                      (input) The month number to compute outflow for.
C
      INTEGER          InSupNTS
C                      (input) The Superior net total supply to use 
C                      for computing mean level/outflow for a gate
C                      setting.  If -9999, indicates use the Plan 1977A
C                      supplies based on trigger level.  Otherwise, 
C                      call to this function was made after averaging 
C                      outflows from Plan77A forecast months.
C
      DOUBLE PRECISION SecsInPeriod     
C                      (input) The length of the routing interval.
C
      DOUBLE PRECISION PreviousFlow
C                      (input) The outflow last routing interval.
C
      DOUBLE PRECISION SuperiorBOPLevel, MicHuronBOPLevel, 
     &                 StClairBOPLevel, ErieBOPLevel
C                      (input) Beginning-of-Period lake levels 
C    
      CHARACTER*(*)    OutputPrefix
C                      (input) String to use at beginning of lines 
C                      written to detailed.log, including indent.
C    
      DOUBLE PRECISION StMarysFlow 
C                      (input/output) The Superior outflow for the 
C                      routing interval.  Initially given as some kind
C                      of target, either from the balancing equation or
C                      forecast averaging.  This routine adjusts and 
C                      modifies this flow according to the various Plan 
C                      77A limitations, and returns the adjusted value.
C
      REAL             Gates
C                      (output) The number of gates open in Month.
C
      CHARACTER*(40)   AdjustmentCodes
C                      (input/output) String containing code letters 
C                      describing what adjustments were made to come
C                      with the Plan flow.
C    
      DOUBLE PRECISION SuperiorEOPLevel, MicHuronEOPLevel, 
     &                 StClairEOPLevel, ErieEOPLevel
C                      (output) End-of-Period lake levels 
C
C                       
C
C
C     Declarations regarding local variables or functions:
C
      CHARACTER*(3)    CodeChar
      DOUBLE PRECISION tmpflow, HalfGateFlow, FullGateFlow
      INTEGER          I, K, L 
C                      Scratch variables
C
      INTEGER          PL
C                      Length of OutputPrefix
C
      DOUBLE PRECISION SupNTS, MHNTS
C                      The Net Total Supply to use for Lakes Superior
C                      and Michigan-Huron (considering trigger levels, 
C                      and diverions)
C
      DOUBLE PRECISION SuperiorAveLevel, MicHuronAveLevel
C                      Mean levels for the routing period, as computed
C                      within Plan 1977A.  NOT necessarily same results
C                      as routing through these lakes after Plan 1977A 
C                      outflow is determined.
C
      DOUBLE PRECISION Round, GetArea
C                      Function for rounding
C
C
C     Declarations regarding COMMON block variables:
C
      INTEGER P77AWinterMax, P77AMaxChange, P77AMinimumFlow(12)
C                      Plan 1977A parameters for bounds-checking
C
      DOUBLE PRECISION P77ASupTargetLev(12), P77ASupRangeFact(12),
     &        P77AMHTargetLev(12), P77AMHRangeFact(12), P77AAValue
      INTEGER P77ABaseFlow(12)
C                      Plan 1977A parameters for Balancing Equation
C
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
C                      Plan 1977A routing parameters
C
      INTEGER P77AMonthsToAve(12)
      LOGICAL P77AUniformMonth
C
C
C
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
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
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
      COMMON /MONAME/   MonthAbbrev
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
C
C
C 
      PL = INDEX(OutputPrefix, LLT) - 1
      IF( PL .LT. 1) STOP ' Norton!'
C
C
C     Log the entering data
C
      IF( SupVerbosity .GE. 3) THEN
         WRITE(CTMP, 6543) OutputPrefix(1:PL), MonthAbbrev(Month), LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6544) OutputPrefix(1:PL), StMarysFlow, 
     &      AdjustmentCodes, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6545) OutputPrefix(1:PL), PreviousFlow, 
     &      SuperiorBOPLevel, MicHuronBOPLevel, StClairBOPLevel,
     &      ErieBOPLevel, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
 6543 FORMAT(A,'Checking and adjusting Plan 1977A outflow for ',A3,A1) 
 6544 FORMAT(A,'Unadjusted flow=', F9.1, '  AdjustmentCodes=', A, A1) 
 6545 FORMAT(A,'Previous flow=', F9.1,' BOP Levels=', 4F8.3, A1) 
C
C
      IF(InSupNTS .EQ. -9999) THEN
C        Check trigger levels and determine the supplies to use
C        for routing this period.  Just put the NBS values in the
C        NTS variables for now; will add in the diversions after
C        write to logs.  Also record an indicator letter into 2nd
C        digit of AdjustmentCodes for which supplies being used.
C
         IF( SuperiorBOPLevel .LE. P77ALowTrigger ) THEN 
            SupNTS = P77ASup95NBS(Month) 
            MHNTS  = P77AMH95NBS(Month) 
            WRITE(CTMP, 6529) OutputPrefix(1:PL), SuperiorBOPLevel,
     &         'less', P77ALowTrigger, '95', LLT 
            AdjustmentCodes(4:6) = 'LS '
         ELSE IF( SuperiorBOPLevel .GE. P77AHighTrigger ) THEN 
            SupNTS = P77ASup05NBS(Month) 
            MHNTS  = P77AMH05NBS(Month) 
            WRITE(CTMP, 6529) OutputPrefix(1:PL), SuperiorBOPLevel,
     &         'greater', P77AHighTrigger, '5', LLT
            AdjustmentCodes(4:6) = 'HS '
         ELSE
            SupNTS = P77ASup50NBS(Month) 
            MHNTS  = P77AMH50NBS(Month) 
            WRITE(CTMP, 6530) OutputPrefix(1:PL), SuperiorBOPLevel, 
     &         P77ALowTrigger, P77AHighTrigger, LLT
            AdjustmentCodes(4:6) = 'NS '
         ENDIF
C
C        Write info to log
C                               `
         IF( SupVerbosity .GE. 5) THEN
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP, 6531) OutputPrefix(1:PL), 'Superior', 
     &         SupNTS + P77ALLOgokiDiv(Month), 
     &         IDNINT(SupNTS), P77ALLOgokiDiv(Month), LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP, 6531) OutputPrefix(1:PL), 'Michigan-Huron', 
     &         MHNTS + P77AChicagoDiv(Month), 
     &         IDNINT(MHNTS), P77AChicagoDiv(Month), LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
 6529    FORMAT(A,'BOP level ',F7.3,1X,A,' than trigger level (', F7.3, 
     &      ') so using ', A, '% supplies.', A1 )
 6530    FORMAT(A,'BOP level ',F7.3,' greater than low trigger level (',
     &      F7.3, '), and less than high trigger level (', F7.3,
     &      '), so 50% supplies.', A1 )
 6531    FORMAT(A, A,' NTS=', F8.1, ' (', I5, ' + ', I4, ')', A1)
C
C
C        Add in diversions.
C
         SupNTS = SupNTS + P77ALLOgokiDiv(Month) 
         MHNTS  = MHNTS + P77AChicagoDiv(Month)
C
C
C
      ELSE
         SupNTS = InSupNTS
         IF( SupVerbosity .GE. 5) THEN
            WRITE(CTMP, 6528) OutputPrefix(1:PL), SupNTS, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
 6528    FORMAT(A,'Using operational Superior net total supply: ',F17.3, 
     &      A1 )
      ENDIF
C
C
C     ===> Requirements C, J, and A
C
C     Check to see if incoming flow (from balancing equation or 
C     average of forecast flows) is less than minimum.  
C
      IF( SupVerbosity .GE. 5) THEN
         WRITE(CTMP, 6541) OutputPrefix(1:PL), LLT
 6541    FORMAT(A,'Checking preliminary flow against minimum:', A1)
         CALL LogIt(LOGD,0,0)
      ENDIF 
C
      IF( IDNINT(StMarysFlow) .LE. P77AMinimumFlow(Month)) THEN
C
C        Flow is too low.
         IF( SupVerbosity .GE. 3) THEN
            WRITE(CTMP, 6533) OutputPrefix(1:PL), StMarysFlow, 
     &         P77AMinimumFlow(Month), LLT
 6533       FORMAT(A,'Preliminary flow of ', F9.1,   
     &         ' m3s below minimum; setting it to ', I6,' m3s', A1)
            CALL LogIt(LOGD,0,0)
         ENDIF 
C    
         AdjustmentCodes(7:12) = 'RC RJ '
         StMarysFlow = P77AMinimumFlow(Month)
         Gates = 0.5
         CALL P77ASolveSup(SuperiorBOPLevel, SupNTS, SecsInPeriod,
     &      -1.0, P77AFlowRound, P77AMonthRound, SuperiorEOPLevel, 
     &      StMarysFlow, SuperiorAveLevel, Month)
      ELSE
C
C        Flow is above minimum.  If desired, log a message.
         IF( SupVerbosity .GE. 7) THEN
            WRITE(CTMP, 6537) OutputPrefix(1:PL), LLT
            CALL LogIt(LOGD,0,0)
         ENDIF 
C     
C        Find out gate position and flow (assuming maximum sidechannel) 
C        closest to entering flow, and return it in tmpflow.
         CALL SupClosestFlow( StMarysFlow, SuperiorBOPLevel, 
     &      MicHuronBOPLevel, SupNTS, SecsInPeriod, P77AMonthRound, 
     &      P77AFlowRound, tmpflow, Gates, SuperiorEOPLevel, Month)
C
         IF( SupVerbosity .GE. 5) THEN
            WRITE(CTMP, 6532) OutputPrefix(1:PL), tmpflow, Gates, LLT
 6532       FORMAT(A,'Closest flow (assuming maximum sidechannel)=', 
     &         F8.1, ' (', F5.2,' gates)', A1)
            CALL LogIt(LOGD,0,0)
         ENDIF 
C       
C
C        XXX Carl-compatibility fix.  Somehow Carl's model sets the plan
C        XXX flow to permit full sidechannel flow at 1/2 gate, if the 
C        XXX entering flow is within 1 m3s of the 1/2 gate full 
C        XXX sidechannel flow.  Nuke this line and uncomment the next
C        XXX one when backward compatibility not an issue regrading P77A
         IF(tmpflow .GT. (StMarysFlow +1) .AND. NINT(Gates*2) .EQ.1)THEN
C        XXX Matt's preference:
C        IF( tmpflow .GT. StMarysFlow .AND. NINT(Gates * 2) .EQ. 1) THEN
C           
C           If closest flow is greater than entering flow, and only half
C           gate open, then can't have maximum sidechannel flow.  Since
C           entering flow is above the minimum, and less than 1/2 gate
C           plus full sidechannel, it is a valid flow itself, so just
C           leave it unadjusted in StMarysFlow.  
C
            IF( SupVerbosity .GE. 5) THEN
              WRITE(CTMP, 6542) OutputPrefix(1:PL), LLT
 6542         FORMAT(A,'Unnecessary to adjust flow according to even ',
     &           'gate setting (sidechannel less than maximum).',A1)
              CALL LogIt(LOGD,0,0)
            ENDIF 
C
            CALL P77ASolveSup(SuperiorBOPLevel, SupNTS, SecsInPeriod,
     &         -1.0, P77AFlowRound, P77AMonthRound, SuperiorEOPLevel, 
     &         StMarysFlow, SuperiorAveLevel, Month)
            AdjustmentCodes(7:12) = '   RJ '
C
C        Note that older model simply presumes that total flow for 
C        1/2 gate is always 2390 m3s.  So to be compatible with older
C        answers, must make same assumption.  If 2390 is greater than 
C        entering flow, and only half gate open, then can't have
C        maximum sidechannel flow, so just leave StMarysFlow unadjusted.  
C        XXX  comment this clause out when don't need backward
C        XXX  compatibility with old bugs!
         ELSE IF( StMarysFlow .LT. 2390 ) THEN
            IF( SupVerbosity .GE. 3) THEN
               WRITE(CTMP, 6601) OutputPrefix(1:PL), LLT
               CALL LogIt(LOGD,0,0)
               WRITE(CTMP, 6602) OutputPrefix(1:PL), LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
 6601       FORMAT(A,'Not adjusting flow according to nearest gate ',
     &         'setting in order to replicate bug in the old model.',A1)
 6602       FORMAT(A,'For more information see the implementation of ',
     &         'Requirements C, J, and A in sub P77ALimitations of ',
     &         'SUPERIOR.FOR.',A1)
            CALL P77ASolveSup(SuperiorBOPLevel, SupNTS, SecsInPeriod,
     &         -1.0, P77AFlowRound, P77AMonthRound, SuperiorEOPLevel, 
     &         StMarysFlow, SuperiorAveLevel, Month)
            AdjustmentCodes(7:12) = '   RJ '
C        XXX End of clause to remove someday
C 
         ELSE
C            
            IF( NINT(Gates) .EQ. 16 .AND. StMarysFlow .GT. tmpflow)THEN
C              Entering flow is greater than capacity! (Requirement A) 
               AdjustmentCodes(16:18) = 'RA '
            ELSE
               AdjustmentCodes(16:18) = '   '
            ENDIF
C              
C           Entering flow high enough for maximum sidechannel flow with
C           half a gate.  Use the closest flow.
C           
            StMarysFlow = tmpflow
            AdjustmentCodes(7:12) = '      '
C
C           XXX Another older model compatibility fix.  In order to replicate
C           XXX rounding in older model, must run subroutine again.
C           XXX Nuke this block when backward compatibility not an issue.
            CALL P77ASolveSup(SuperiorBOPLevel, SupNTS, SecsInPeriod,
     &         -1.0, P77AFlowRound, P77AMonthRound, SuperiorEOPLevel, 
     &         StMarysFlow, SuperiorAveLevel, Month)
C           XXX end of fix
C
         ENDIF
      ENDIF
C
C        
C
C
C
C     Check rate-of-change limitation (Requirement I)
C
      IF( SupVerbosity .GE. 5) THEN
         WRITE(CTMP, 6539) OutputPrefix(1:PL), LLT
         CALL LogIt(LOGD,0,0)
      ENDIF 
C
C     XXX Hardcoded rounding in Plan1977A:
      IF( IABS(IDNINT(StMarysFlow) - IDNINT(PreviousFlow)) 
     &   .GE. P77AMaxChange) THEN
C
         IF( SupVerbosity .GE. 3) THEN
            WRITE(CTMP, 6536) OutputPrefix(1:PL), StMarysFlow, Gates, 
     &         P77AMaxChange, PreviousFlow, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF 
C        Flow change is too great.  Loop around until jump out.
         DO 5734, I=1,20
C
            IF( SupVerbosity .GE. 3) THEN
               WRITE(CTMP, '(A, A, F5.2, A1)') OutputPrefix(1:PL),
     &             ' Checking with gates= ', Gates, LLT
               CALL LogIt(LOGD,0,0)
            ENDIF 
C           Increment/decrement gate setting
            IF(StMarysFlow .GT. PreviousFlow) THEN
C              Flow rising too much this month, back off a gate
               Gates= Gates - 1  
               IF(NINT(Gates) .EQ. 0) Gates = 0.5
C 
C              Really bizarre numbers mean a bust somewhere
               IF(StMarysFlow.LT.50 .OR. PreviousFlow.LT.50)STOP 'Lucy!'
C
               IF(Gates .LT. 0.0) THEN
C               
C                 Hmmm.  Max sidechannel plus 1/2 gate is more than 850 
C                 higher than previous flow.  So can't have max 
C                 sidechannel.  Set flow to 850 plus previous flow, and
C                 recompute EOP level.
                  Gates = 0.5
                  StMarysFlow = IDNINT(PreviousFlow) + P77AMaxChange
C                 Check for inability to reconcile Requirements C, I, 
C                 and J.  Per Nanette's advice, let C control over I 
C                 (stick with 1560 m3s, don't worry about 850 change).   
                  IF( IDNINT(StMarysFlow) .LT. P77AMinimumFlow(Month))
     &               StMarysFlow = P77AMinimumFlow(Month)
                  CALL P77ASolveSup(SuperiorBOPLevel, SupNTS, 
     &               SecsInPeriod, -1.0, P77AFlowRound, P77AMonthRound,
     &               SuperiorEOPLevel, StMarysFlow, SuperiorAveLevel, 
     &               Month)
                  IF( SupVerbosity .GE. 3) THEN
                     WRITE(CTMP,6102) OutputPrefix(1:PL), 
     &                  IDNINT(StMarysFlow), LLT
                     CALL LogIt(LOGD,0,0)
                  ENDIF
C                 Now leave loop                    
                  GOTO 5735 
               ENDIF
            ELSE
C              Flow dropping too much this month
C              
               IF(NINT(Gates - 0.1) .EQ. 0) THEN
C
C                 See if previousflow - 850 is between 1/2 and 1 gate
C                 If so, set it at previousflow - 850 
C
                  IF( SupVerbosity .GE. 6) THEN
                     WRITE(CTMP, 6298) OutputPrefix(1:PL), '0.5', LLT
                     CALL LogIt(LOGD,0,0)
                  ENDIF 
                  CALL P77ASolveSup(SuperiorBOPLevel, SupNTS, 
     &               SecsInPeriod, 0.5, P77AFlowRound, P77AMonthRound, 
     &               SuperiorEOPLevel, HalfGateFlow, SuperiorAveLevel, 
     &               Month)
C
                  IF( SupVerbosity .GE. 6) THEN
                     WRITE(CTMP, 6298) OutputPrefix(1:PL), '1.0', LLT
                     CALL LogIt(LOGD,0,0)
                  ENDIF 
                  CALL P77ASolveSup(SuperiorBOPLevel, SupNTS, 
     &               SecsInPeriod, 1.0, P77AFlowRound, P77AMonthRound, 
     &               SuperiorEOPLevel, FullGateFlow, SuperiorAveLevel, 
     &               Month)
                  IF( SupVerbosity .GT. 2) THEN
                     WRITE(CTMP,6104) OutputPrefix(1:PL), HalfGateFlow, 
     &                  FullGateFlow, LLT
                     CALL LogIt(LOGD,0,0)
                  ENDIF
C
C                 XXX Fix to accomodate bug in older model, which
C                 XXX allows the flow to change more than 850 m3s
C                 XXX if 1/2 gate with full sidechannel is closer to
C                 XXX target (PreviousFlow -P77AMaxChange) than 1 gate
C                 XXX with full sidechannel, when the target flow is 
C                 XXX the max for the month, rounded to the nearest 10m3s
C                 XXX Simply remove this block to remove the "fix":
C
                  IF( Month .LE. 4 .OR. Month .EQ. 12) THEN
                     L = P77AWinterMax
                  ELSE
                     L = (12 - Month) * P77AMaxChange + P77AWinterMax
                  ENDIF
                  IF(ABS( (IDNINT(PreviousFlow) - P77AMaxChange) -
     &               IDNINT( HalfGateFlow ) )  .LT. 
     &               ABS( (IDNINT(PreviousFlow) - P77AMaxChange) -
     &               IDNINT( FullGateFlow ) ) .AND. 
     &               IDNINT( (PreviousFlow - P77AMaxChange) / 10.0d0 ) 
     &               * 10 .EQ. L ) THEN
C                    1/2 gate flow with full sidechannel is closer to 
C                    last month's flow minus 850, and full gate flow
C                    exceeds max for month.                       
                     StMarysFlow = HalfGateFlow
                     IF( SupVerbosity .GT. 5) THEN
                        WRITE(CTMP, 6108) OutputPrefix(1:PL),
     &                     ' Running P77ASolveSup in Req I compat fix'//
     &                     ' to determine new mean level', LLT
                        CALL LogIt(LOGD,0,0)
                     ENDIF
                     CALL P77ASolveSup(SuperiorBOPLevel, SupNTS, 
     &                  SecsInPeriod, -1., P77AFlowRound,P77AMonthRound, 
     &                  SuperiorEOPLevel, StMarysFlow, 
     &                  SuperiorAveLevel, Month)
C                    Now leave loop                    
C                    Leave gates = 0.5
                     GOTO 5735 
                  ELSE
                     IF( SupVerbosity .GT. 5) THEN
                        WRITE(CTMP, 6108 ) OutputPrefix(1:PL), 
     &                     ' Req I compat fix not required', LLT
                        CALL LogIt(LOGD,0,0)
C                     IF( NINT(Gates) .LE. 1 .AND. HalfGateFlow < 2380) 
C     &                  GOTO 5727
                     ENDIF
                  ENDIF
C                 XXX End of "fix" block
C
C                 Say that StMarysFlow 
C
                  StMarysFlow = FullGateFlow
                  Gates = 1
C
C
C
                  IF((PreviousFlow - P77AMaxChange).LT.StMarysFlow)THEN
C
                     IF( SupVerbosity .GT. 5) THEN
                        WRITE(CTMP, 6834) OutputPrefix(1:PL),
     &                     StMarysFlow, Gates, 'now', P77AMaxChange,
     &                     IDNINT(PreviousFlow), LLT
                        CALL LogIt(LOGD,0,0)
                        WRITE(CTMP, '(A,A,I6,A1)' ) OutputPrefix(1:PL), 
     &                     ' Setting StMarysFlow to PreviousFlow - ',
     &                     P77AMaxChange, LLT
                        CALL LogIt(LOGD,0,0)
                     ENDIF
                     StMarysFlow = IDNINT(PreviousFlow) - P77AMaxChange
C
C                    XXX Yet another odd bit logic intended to replicate 
C                    XXX answers from the older models.
                     IF( StMarysFlow .LT. 2390 .AND. 
     &                  StMarysFlow .LE. HalfGateFlow) THEN 
                        Gates = 0.5
                     ELSE
                        Gates = 1.0
                        IF( SupVerbosity .GT. 5) THEN
                           WRITE(CTMP, 6108 ) OutputPrefix(1:PL), 
     &                        ' Running P77ASolveSup after setting '//
     &                        '1.0 gate for 2390 m3/s compat fix ', LLT
                           CALL LogIt(LOGD,0,0)
                        ENDIF
                        CALL P77ASolveSup(SuperiorBOPLevel, SupNTS, 
     &                     SecsInPeriod, Gates, P77AFlowRound,
     &                     P77AMonthRound, SuperiorEOPLevel, 
     &                     StMarysFlow, SuperiorAveLevel, Month)
                     ENDIF
C                    XXX End of "fix" block
                     IF( SupVerbosity .GT. 5) THEN
                        WRITE(CTMP, '(A,A,F5.2,A)') OutputPrefix(1:PL), 
     &                     ' Gates set to ', gates, LLT
                        CALL LogIt(LOGD,0,0)
                     ENDIF
                     IF( SupVerbosity .GT. 5) THEN
                        WRITE(CTMP, 6108 ) OutputPrefix(1:PL), 
     &                     ' Running P77ASolveSup after Req I compat '//
     &                     'fix to determine new mean level', LLT
                        CALL LogIt(LOGD,0,0)
                     ENDIF
                     CALL P77ASolveSup(SuperiorBOPLevel, SupNTS, 
     &                  SecsInPeriod, -1., P77AFlowRound,P77AMonthRound,
     &                  SuperiorEOPLevel, StMarysFlow, SuperiorAveLevel,
     &                  Month)
C                    Now leave loop                    
                     GOTO 5735 
                  ELSE
                     IF( SupVerbosity .GT. 5) THEN
                        WRITE(CTMP, 6834) OutputPrefix(1:PL),
     &                     StMarysFlow, Gates, 'not', P77AMaxChange,
     &                     IDNINT(PreviousFlow), LLT
                        CALL LogIt(LOGD,0,0)
                     ENDIF
                  ENDIF
C
               ENDIF
C               
C              Increase a gate and try again
               Gates= NINT(Gates - 0.1) + 1  
               IF(Gates .GT. 16.0) THEN
                  CALL ErrMess(1171,'P77ALimitations'//LLT,LOGCDE)
                  WRITE(CTMP, 6519)OutputPrefix(1:PL),
     &               MonthAbbrev(Month), SuperiorAveLevel, Gates,
     &               IDNINT(StMarysFlow), IDNINT(PreviousFlow), LLT
                  CALL LogIt(LOGCDE,0,0)
                  CALL SAVETOD
                  CALL ADIOS
               ENDIF
            ENDIF
C
C           Get new mean outflow for adjusted gate setting
            IF( SupVerbosity .GT. 5) THEN
               WRITE(CTMP, 6108 ) OutputPrefix(1:PL), 
     &            ' Running P77ASolveSup to compute new outflow for '//
     &            ' adjusted gate setting . . .', LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
            CALL P77ASolveSup(SuperiorBOPLevel, SupNTS, SecsInPeriod, 
     &         Gates, P77AFlowRound, P77AMonthRound, 
     &         SuperiorEOPLevel, StMarysFlow, SuperiorAveLevel, Month)
C
C           XXX Another older model compatibility fix.  In order to replicate
C           XXX rounding in older model, must run subroutine again. 
C           XXX Nuke this block when backward compatibility with earlier
C           XXX models regarding P77A is not an issue.
            IF( SupVerbosity .GT. 8) THEN
               WRITE(CTMP, 6108 ) OutputPrefix(1:PL), 
     &            ' Running P77ASolveSup again to re-compute mean '//
     &            ' level, as done in old code . . .', LLT
               CALL LogIt(LOGD,0,0) 
            ENDIF
            CALL P77ASolveSup(SuperiorBOPLevel, SupNTS, 
     &         SecsInPeriod, -1., P77AFlowRound,P77AMonthRound, 
     &         SuperiorEOPLevel, StMarysFlow, SuperiorAveLevel, Month)
C           XXX End of "fix" block
C
C
C        XXX hardcoded rounding to nearest integer flow replicates 
C        XXX previous model, but is not flexible.
 5734    IF( IABS(IDNINT(StMarysFlow) - IDNINT(PreviousFlow)) 
     &      .LE. P77AMaxChange) GOTO 5735
C
         STOP 'Little Ricky!'
 5735    IF( SupVerbosity .GE. 3) THEN
            WRITE(CTMP,6535) OutputPrefix(1:PL), StMarysFlow, Gates,LLT
            CALL LogIt(LOGD,0,0)
         ENDIF 
         AdjustmentCodes(13:15) = 'RI '
      ELSE
         AdjustmentCodes(13:15) = '   '
         IF( SupVerbosity .GE. 7) THEN
            WRITE(CTMP, 6537) OutputPrefix(1:PL), LLT
            CALL LogIt(LOGD,0,0)
         ENDIF 
      ENDIF
 6102 FORMAT( A, 'Unable to allow maximum sidechannel flow;',
     &   ' Req. I limits flow to ', I7, ' m3s', A1)      
 6104 FORMAT( A, 'Flow for .5 and 1 gates with max sidechannel flow: ',
     &   2F8.1, A1)      
 6108 FORMAT( A, A, A1)      
 6298 FORMAT (A, ' Running P77ASolveSup() to determine ', A, 
     &   ' gate flow (full sidechan) . . .', A1)

 6536 FORMAT(A,'Preliminary flow of ', F6.1,' m3s (', F5.2,
     &   ' gates) changed more than ', I4,
     &   ' from previous period (', F6.1, ')', A1)
 6539 FORMAT(A,
     &   'Checking preliminary flow against rate-of-change limitation:',
     &    A1)
 6834 FORMAT(A,' StMarysFlow of ',F10.3,' with ',F5.2,' gates open is ',
     &       A,' within ',I4,' m3/s of previous flow (',I4,')',A1) 
C
C
C     Check maxmimum flow limitations (Requirements B, K) 
C     (Requirement A checked earlier)
C
C        Requirement B:  Maximum outflow Dec-Apr = 2410 m3s
C        Requirement K:  Maximum outflow Oct = 4110 m3s, Nov = 3260 m3s
C
      IF( SupVerbosity .GE. 5) THEN
         WRITE(CTMP, 6540) OutputPrefix(1:PL), LLT
         CALL LogIt(LOGD,0,0)
      ENDIF 
C
      IF( Month .LE. 4 .OR. Month .EQ. 12) THEN
C        It's winter.  Set maximum to P77AWinterMax
         L = P77AWinterMax
         CodeChar = 'RB '
      ELSE
C        It's not winter.  Set maximum to whatever we have to stay
C        under in order to get to P77AWinterMax by December, without
C        violating P77AMaxChange.
         L = (12 - Month) * P77AMaxChange + P77AWinterMax
         CodeChar = 'RK '
C         WRITE(*,'(A,I4)') 'Max = ', L
      ENDIF
      IF( SupVerbosity .GE. 7) THEN
         WRITE(CTMP,6538)OutputPrefix(1:PL), MonthAbbrev(Month), L, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF 
 6538 FORMAT(A,'Maximum flow for ', A3, ': ', I8, ' m3s', A1)
C
      IF( IDNINT(StMarysFlow) .GT. L) THEN
C        Flow is too high.  Loop around until jump out.
         IF( SupVerbosity .GE. 3) THEN
            WRITE(CTMP, 6534)OutputPrefix(1:PL), StMarysFlow,Gates,L,LLT
            CALL LogIt(LOGD,0,0)
         ENDIF 
         K = 0
         DO 5434, I=1,99
C           Decrement gate setting
            Gates= Gates - 1  
            IF(NINT(Gates) .LE. 0) THEN
                Gates = 0.5
C               Set K so that we will know if loop back here again
                K = K + 1
                IF ( K .EQ. 2 ) THEN
C                  Last time through was half gate also, and flow was too
C                  high for Requirements B or K.  So have conflict between
C                  Req J and B/K.  Per David Fay on 7Aug03, let the 1/2 gate
C                  Requirement control.
C                  
                   WRITE(CTMP, 6546)OutputPrefix(1:PL),
     &                MonthAbbrev(Month), SuperiorAveLevel, 
     &                IDNINT(StMarysFlow), L, LLT
                   CALL LogIt(LOGD,0,0)
                   WRITE(CTMP, '(A,A,A,A1)')OutputPrefix(1:PL),
     &                ' Requirement J controls; ', 
     &                'maximum winter flow disregarded.', LLT
                   CALL LogIt(LOGD,0,0)
C                  Now jump back in at end of Requirement checking
                   GOTO 5435
                ENDIF
            ENDIF
C
C           Get new mean outflow for adjusted gate setting
            CALL P77ASolveSup(SuperiorBOPLevel, SupNTS, SecsInPeriod, 
     &         Gates, P77AFlowRound, P77AMonthRound, 
     &         SuperiorEOPLevel, StMarysFlow, SuperiorAveLevel, Month)
C
C           If make it below outflow ceiling jump out of gate-closing loop 
 5434       IF( IDNINT(StMarysFlow) .LE. L) GOTO 5435
C
C        Should never get here - either get below outflow ceiling or hit 
C        1/2 gate twice in a row.
         STOP 'Trixie'
C
 5435    IF( SupVerbosity .GE. 3) THEN
            WRITE(CTMP,6535) OutputPrefix(1:PL), StMarysFlow, Gates,LLT
            CALL LogIt(LOGD,0,0)
         ENDIF 
         AdjustmentCodes(16:18) = CodeChar
      ELSE
         IF( SupVerbosity .GE. 7) THEN
            WRITE(CTMP, 6537) OutputPrefix(1:PL), LLT
            CALL LogIt(LOGD,0,0)
         ENDIF 
      ENDIF
 6534 FORMAT(A,'Preliminary flow of ', F6.1, ' m3s (', F5.2,
     &   ' gates) above maximum of ', I8,' m3s', A1)
 6540 FORMAT(A,'Checking preliminary flow against maximum:', A1)
 6546 FORMAT(A,' Month: ', A3, ' working Sup level: ', F6.2, 
     &   ' Req J (1/2 gate) flow:', I5, 
     &   ' in conflict with Req B/K max flow:', I5, A1)
 6519 FORMAT(A,' Month: ', A3, ' working Sup level: ', F8.2, 
     &   ' Gates: ', F5.2, ' St. Marys flow:', I9,
     &   ' Req I (prev) flow:', I9, A1)
C
C
C     XXX Another ugly re-compute necessary to exactly match rounding
C     XXX practices in older model.  Someday just remove it.
      SuperiorEOPLevel = Round( SuperiorBOPLevel, P77AMonthRound ) + 
     &   Round( (SupNTS - StMarysFlow) * SecsInPeriod /
     &   GetArea(1,SuperiorBOPLevel), P77AMonthRound )
C     XXX End of fix clause 
C
C
C     Now route the flow downstream, so can have a MH level for the
C     next period of forecast, using the method defined in Plan1977A.
C     If called to do check on average of forecasted flows, don't
C     need to do this.
      IF(InSupNTS .EQ. -9999) CALL P77AMidLakeRoute ( SecsInPeriod,
     &   IDNINT(Round(MHNTS + StMarysFlow,P77AFlowRound)), 
     &   P77AStC50NBS(Month),
     &   P77AErie50NBS(Month) + P77AWellandFlow(Month), 
     &   MicHuronBOPLevel,StClairBOPLevel,ErieBOPLevel,
     &   P77AStClairRetrd(Month), 
     &   P77ADetroitRetrd(Month), P77ANiagaraRetrd(Month), 
     &   P77AFlowRound, P77ALevelRound, P77AMonthRound, P77AMidLakesInc,
     &   MicHuronEOPLevel,StClairEOPLevel,ErieEOPLevel,MicHuronAveLevel)
C 
C
C 6525 FORMAT( A,'Previous month outflow=', F20.3, ', gates open=', 
C     &   F5.2, '.', A1) 
 6535 FORMAT(A,'Adjusted preliminary flow to ', F6.1, ' m3s (', F5.2,
     &   ' gates)', A1)
 6537 FORMAT(A, ' . . . passed check OK.', A1)
      END
C
C
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
C
C
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
C
      SUBROUTINE SupCriterions( Year, Month, SuperiorNWB, MichuronNWB,
     &   StClairNWB, ErieNWB, SCRivRet, DetRivRet, NiaRivRet,
     &   SuperiorBOPLevel, MicHuronBOPLevel, StClairBOPLevel,
     &   ErieBOPLevel, LStoSWP, SecsInPeriod, StMarysFlow, SWPier, 
     &   NonGatedFlow, Gates, AdjustmentCodes )
C
C     This function takes a preliminary flow, from either Plan 1977A (or
C     potentially in the future from the 1955 Modified Rule of '49), and 
C     checks/modifies the flow according to Criterions B and C of the 
C     Orders of Approval (at least according to the operational 
C     practices in effect in October 2000).
C
C          Criterion A:  (not applicable in this context)
C          Criterion B:  Prevent levels at US Slip greater than 177.94m
C          Criterion C:  Outflows capped at pre-project below 183.40m
C
C
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER          Year, Month
C                      (input) The date of the outflow.
C
      INTEGER          SuperiorNWB
C                      (input) The net water balance (NBS + LLO) to 
C                      use for Lake Superior change-in-storage 
C                      calculations.  Could be determined in three 
C                      different ways, depending on value of P77ACompat.
      INTEGER          MichuronNWB
C                      (input) The net water balance (NBS + ChiDiv), 
C                      calculated according to value of P77ACompat.
      INTEGER          StClairNWB
C                      (input) The net water balance (NBS), 
C                      calculated according to value of P77ACompat.
      INTEGER          ErieNWB
C                      (input) The net water balance (NBS - Welland), 
C                      calculated according to value of P77ACompat.
      INTEGER          SCRivRet, DetRivRet, NiaRivRet
C                      (input) Flow retardation values to use, 
C                      calculated according to value of P77ACompat.
C
      DOUBLE PRECISION SuperiorBOPLevel, MicHuronBOPLevel, 
     &                 StClairBOPLevel, ErieBOPLevel
C                      (input) Beginning-of-Period lake levels 
C    
      INTEGER          LStoSWP
C                      (input) Which Lake Superior to SW Pier 
C                      relationship to use.
C
      DOUBLE PRECISION SecsInPeriod     
C                      (input) The length of the routing interval.
C
      DOUBLE PRECISION StMarysFlow
C                      (input/output) The Superior outflow for the 
C                      routing interval.  Initially given as some kind
C                      of target, either from the balancing equation or
C                      forecast averaging.  This routine adjusts and 
C                      modifies this flow according to Criterions B and
C                      C, and returns the adjusted value in StMarysFlow.
C
      DOUBLE PRECISION SWPier, NonGatedFlow
C                      (output) The SWPier level as the old models would
C                      compute it, which yields a "gated flow" which
C                      yields the 'sidechannel' flow.
C
      REAL             Gates
C                      (input/output) The number of gates open in Month.
C
      CHARACTER*(40)   AdjustmentCodes
C                      (input/output) String containing code letters 
C                      describing what adjustments were made.
C    
C
C                       
C
C
C     Declarations regarding local variables or functions:
C
      DOUBLE PRECISION tmplev, gatedflow
      INTEGER          I, L, Flow2ItsAgo
C                      Scratch variables
C
      DOUBLE PRECISION SuperiorAveLevel, MicHuronAveLevel
C                      Mean levels for the routing period, as computed
C                      within Plan 1977A.  NOT necessarily same results
C                      as routing through these lakes after Plan 1977A 
C                      outflow is determined.
C
      DOUBLE PRECISION SuperiorEOPLevel, MicHuronEOPLevel, 
     &                 StClairEOPLevel, ErieEOPLevel
C                      End-of-Period lake levels 
C    
C
      CHARACTER*(*)    OutputPrefix
C                      Used as prefix on lines of detailed.log 
C
      DOUBLE PRECISION Round
C                      Function for rounding.
C
      DOUBLE PRECISION GetArea
C                      Returns area of lake
C
      DOUBLE PRECISION USS
C                      Function that estimates level at US Slip based
C                      on St Marys flow and period average level of
C                      Michigan-Huron
C
      DOUBLE PRECISION GetLevelAtUSSlip 
C                   Returns US Slip level corresponding to given 
C                   St Marys flow and MH elevation
C
      DOUBLE PRECISION Q1887
C                      Pre-project flow
C
C
C
C     
C
C     Declarations regarding COMMON block variables:
C
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
      INTEGER      MonthLength
C                  1 if uniform, 2 if actual
C
      LOGICAL      LakesToSolve(8)
      INTEGER      StartYear, StartMonth, StartDay,
     &             EndYear, EndMonth, EndDay, IntervalNameLen(4)
C                  Just declared because in same COMMON as MonthLength
C
      CHARACTER    SupPreProject*6, BalancingFlow*8, SupGateFlow*6,
     &             SupPlanFlow*6, SupBOMLevel*8, SWPierLevel*8,
     &             SupSideChanFlow*6, SPRegCodes*40, P77AHeader*120,
     &             P2012Header
C                  Lake Superior regulation output strings
      INTEGER      CritBOMRound, CritFlowRound, CritLevelRound, 
     &                  Q1887Method, CritBMethod
      INTEGER      USSlipEquation
C                  1=ancient eq (default), 2=2010 EC eq
C
      DOUBLE PRECISION CritCThresholdEl, USSIceWeedAdj
C                  Rounding values for how the criterions get calculated,
C                  and some choices regarding methodology
C
C
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /MONAME/   MonthAbbrev
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /RCRUNI/   LakesToSolve, 
     &                  IntervalNameLen, MonthLength,
     &                  StartYear, StartMonth, StartDay, 
     &                  EndYear, EndMonth, EndDay
      COMMON /SUPOUT/   SupPreProject, BalancingFlow, SupGateFlow,
     &                  SupPlanFlow, SupBOMLevel, SWPierLevel, 
     &                  SupSideChanFlow, SPRegCodes, P77AHeader,
     &                  P2012Header
      COMMON /SUPCRT/   CritBOMRound, CritFlowRound, CritLevelRound, 
     &                  Q1887Method, CritBMethod, CritCThresholdEl
      COMMON /USSLIP/   USSIceWeedAdj, USSlipEquation
C
C
      PARAMETER (OutputPrefix = 'Superior Criterion check: ')
C 
C
C     Log the entering data
C
      IF( SupVerbosity .GE. 5) THEN
         WRITE(CTMP, 6546) OutputPrefix, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6542) OutputPrefix,MonthAbbrev(Month),Year, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6544) OutputPrefix, StMarysFlow, Gates,
     &      AdjustmentCodes, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6545) OutputPrefix, SuperiorBOPLevel, 
     &      MicHuronBOPLevel, StClairBOPLevel, ErieBOPLevel, 
     &      SecsInPeriod, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6541) OutputPrefix, SuperiorNWB, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6543) OutputPrefix, MichuronNWB,
     &      StClairNWB, ErieNWB, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6540) OutputPrefix, SCRivRet, DetRivRet, 
     &      NiaRivRet, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
 6541 FORMAT(A,'Superior NTS=', I6, A1) 
 6542 FORMAT(A,'Evaluating plan flow for ',A3,I5,A1) 
 6543 FORMAT(A,'MH/SC/ER net water balances: ',3I7,A1) 
 6540 FORMAT(A,'MH/SC/ER outflow retardations: ',3I7,A1) 
 6544 FORMAT(A,'Plan flow=', F9.1, ' (', F5.2, ' gates)', 
     &   '  AdjustmentCodes=', A, A1) 
 6545 FORMAT(A,'BOP Levels=', 4F8.3, ' dT=', F14.3, A1) 
 6546 FORMAT(A,'Checking preliminary flow against level at US Slip:',A1)
C
C
      SuperiorBOPLevel = Round( SuperiorBOPLevel, CritBOMRound )
C
C     Check Criterion C
C
C     Pre-existing practice is to estimate the pre-project flow using the
C     mean monthly Superior level associated with the outflow from Plan
C     1977A/Criterion B.  
C
C     More proper computation would be to solve 3 simultaneous equations: 
C     pre-project eqn for outflow, mean level,
C     and change-in-storage.  
C
C
C     Env Canada code:
C
C     Apply criterion C of IJC Orders here.  If SUPBOP < 183.40m
C     then compute 1887 relation and test against SUPFLOW.
C     Use the Point Iroquois mtric IGLD 1985 relation
C
C     IF(SUPMEAN.LT.DBLE(183.40))THEN
C
C        if in range of equation, else set flow=10 and set flag=F
C
C        IF(SUPMEAN.GT.DBLE(181.43))THEN
C           Q1887=ROUND((DBLE(82.47)*(SUPMEAN-DBLE(181.43))**DBLE(1.5)),DBLE(10))
C        ELSE
C           Q1887=10.
C           CRTRNC=3
C        ENDIF
C
C        IF(SUPFLOW.GT.Q1887)THEN
C           criterion C applies
C           CRTRNC=1
C           CALL QSUP(Q1887, SUPBOP, SUPNTS, X, Y, Z, QFLOW,
C    1         SIDECHAN, GATESET, MONTH, 1)
C
C           compute final criterion C values for Lake Superior
C
C           SUPFLOW = QFLOW
C           SUPEOP = ROUND(SUPBOP + (SUPNTS - SUPFLOW)/3138,DBLE(100))
C           SUPMEAN = ROUND((SUPBOP + SUPEOP)/2, DBLE(100))
C           DATUM  = SUPMEAN - SUPDATUM
C        ENDIF
C     ENDIF
C
C
C     First we need to figure out Superior mean level (and possibly 
C     outflow) based on Plan 1977A.  
C
C     XX TODO: see if this flag for determining whether full sidechannel capacity,
C     XX causes logic to be dependent on using Plan 1977A 
C     15Oct2010 Yes - it modifies plan flow if fails to find 'RJ', which means 
C     need to fix this some way for other plans to use routine to do crit c
C     checks without inadvertently changing plan flow
      IF( INDEX(AdjustmentCodes, 'RJ') .EQ. 0) THEN
C        Full sidechannel capacity.  Calculate P77A mean level
C        and outflow based on gates.
         CALL P77ASolveSup(SuperiorBOPLevel, 
     &      DBLE(SuperiorNWB), SecsInPeriod, Gates, CritFlowRound, 
     &      CritBOMRound, SuperiorEOPLevel, StMarysFlow, 
     &      SuperiorAveLevel, Month)
      ELSE
C        Partial sidechannel capacity.  Calculate P77A mean level
C        based on StMarysFlow
         CALL P77ASolveSup(SuperiorBOPLevel,
     &      DBLE(SuperiorNWB), SecsInPeriod, -1.0, CritFlowRound, 
     &      CritBOMRound, SuperiorEOPLevel, StMarysFlow, 
     &      SuperiorAveLevel, Month)
      ENDIF  
C
      IF( SupVerbosity .GE. 5) THEN
         WRITE(CTMP,6512) OutputPrefix, SuperiorAveLevel, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF 
 6512 FORMAT( A,'trying Criterion C, with mean Sup level= ', F12.6, A1)
C
C     XX TODO: create logic to punt and return some flow if elev too 
C     XX low for 1887 eqn
      IF(SuperiorAveLevel .LT. 181.43D0) STOP ' Mr. Howell!'
C
      IF( Q1887Method .EQ. 1) THEN
         Q1887= Round(824.7D0*(SuperiorAveLevel-181.43D0)**1.5D0,
     &      CritFlowRound)
      ELSE IF( Q1887Method .EQ. 2) THEN
         CALL GetSupPreProject( SuperiorBOPLevel, DBLE(SuperiorNWB), 
     &      SecsInPeriod, SupVerbosity, Month, Q1887, tmplev)
         Q1887= Round(Q1887, CritFlowRound)
      ENDIF 
C
      IF( SupVerbosity .GE. 4) THEN
         WRITE(CTMP,'(A, A3,I5, A, F15.10,A1)')' Pre-project flow for ', 
     &      MonthAbbrev(Month), Year, ': ',  Q1887, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF 
C
      WRITE(SupPreProject, '(I6)') IDNINT(Q1887) 
C
      IF( SuperiorAveLevel .LT. CritCThresholdEl) THEN
C
         IF( SupVerbosity .GE. 3) THEN
            WRITE(CTMP,'(A, A3,I5, A, F15.10,A,F10.5,A,A1)')
     &         '   Must consider Crit C for ', MonthAbbrev(Month), Year,
     &         ' because Superior level for Crit C test (', 
     &         SuperiorAveLevel, 'm) is below threshold of ',
     &         CritCThresholdEl, 'm.', LLT
            CALL LogIt(LOGD,0,0)
         ENDIF 
C
         DO 845, I = 1, 999
C
C           Jump out if not a problem 
            IF(StMarysFlow .LE. Q1887) GOTO 9763
C
C           Violated Criterion C!
            AdjustmentCodes(22:24) = 'CC '
            IF( SupVerbosity .GT. 2) THEN
               WRITE(CTMP,6550)OutputPrefix,MonthAbbrev(Month),Year,LLT
               CALL LogIt(LOGD,0,0)
            ENDIF
C
            IF( NINT(Gates + 0.51) .EQ. 1 ) THEN
C              Only half a gate open!  Can't close any more! 
               StMarysFlow = Q1887
               IF( SupVerbosity .GT. 1) THEN
                  WRITE(CTMP,'(A, A3,I5,A)') OutputPrefix, 
     &                MonthAbbrev(Month), Year,
     &                '  Re-computing mean Superior '//
     &                'level changing regulated flow to '//
     &                'pre-project flow.'//LLT
                  CALL LogIt(LOGD,0,0)
               ENDIF
               CALL P77ASolveSup( SuperiorBOPLevel, DBLE(SuperiorNWB), 
     &            SecsInPeriod, -1.0, CritFlowRound, CritBOMRound, 
     &            SuperiorEOPLevel, StMarysFlow, SuperiorAveLevel, 
     &            Month)
               IF( SupVerbosity .GT. 1) THEN
                  WRITE(CTMP,6551) OutputPrefix, MonthAbbrev(Month),
     &               Year,StMarysFlow, LLT
                  CALL LogIt(LOGD,0,0)
               ENDIF
               GOTO 9764
            ENDIF
C
C           Decrement gate setting
            Gates = Gates - 1  
            IF(NINT(Gates) .EQ. 0) Gates = 0.5
            IF(Gates .LT. 0.0) STOP ' Professor'
C
C           Get new mean outflow for adjusted gate setting, loop again
  845       CALL P77ASolveSup(SuperiorBOPLevel,
     &         DBLE(SuperiorNWB), SecsInPeriod, Gates, CritFlowRound, 
     &         CritBOMRound, SuperiorEOPLevel, StMarysFlow, 
     &         SuperiorAveLevel, Month)
C
         STOP ' Skipper'   
C
         ELSE 
C           Superior Level > CritCThresholdEl
            IF( SupVerbosity .GE. 3) THEN
            WRITE(CTMP,'(A, A3,I5, A, F15.10,A,F10.5,A,A1)')
     &         '   Not considering Crit C for ', MonthAbbrev(Month), 
     &         Year, ' because Superior level for Crit C test (', 
     &         SuperiorAveLevel, 'm) is above threshold of ',
     &         CritCThresholdEl, 'm.', LLT
            CALL LogIt(LOGD,0,0)
         ENDIF 
C
      ENDIF
 9763 IF(SupVerbosity .GE. 4 .AND. AdjustmentCodes(22:24) .NE.'CC ')THEN
         WRITE(CTMP,6537) OutputPrefix, 'C', LLT
         CALL LogIt(LOGD,0,0)
      ENDIF 

C     XX Re-compute mean so that can reproduce rounding effects of 
C     XX earlier models.
C     XX Don't bother if no rounding, or if using new Crit C method,
C     XX or if didn't close a gate because of Criterion C.   
C     XX Remove if compatibility with older models not required.
      IF(AdjustmentCodes(22:24) .EQ. 'CC ' .AND. Q1887Method .EQ. 1
     &   .AND. CritBOMRound .NE. -99) THEN
         SuperiorAveLevel = Round(SuperiorBOPLevel, CritBOMRound) + 
     &      Round((DBLE(SuperiorNWB) - StMarysFlow) * SecsInPeriod * 
     &      0.5D0 / GetArea(1, SuperiorAveLevel), CritBOMRound)
      ENDIF
C     XX End of fix
C
C
 9764 CONTINUE
C     Check Criterion B level
C
C        There is a problem replicating existing practice for 
C        implementing criterions B and C.  The existing implementations
C        use monthly mean lake levels for Superior and Michigan-Huron,
C        based on supplies for the next month.  This is handled by 
C        the parameter P77A Compat.
C
C         =========================================================
C                      EC code
C
C            USSLIP(YEAR,MONTH) = ROUND(USS(SUPFLOW, MHUMEAN),
C    1                                        DBLE(100))
C
CC                If the US Slip level is too high, decrease the Lake
CC                Superior gate setting and compute a new plan flow
C
C             IF (USSLIP(YEAR,MONTH) .GT. DBLE(177.94)) THEN
C                GATESET = GATE(GATESET, -1)
C                CALL QSUP(DBLE(0), DBLE(0), DBLE(0), X, Y, Z,
C    1              DBLE(0), SIDECHAN, GATESET, MONTH, 2)
C                PLANFLOW = ROUND(FLOW(SUPBOP, SUPNTS, SIDECHAN,
C    1              X, Y, Z), DBLE(10))
C                IF (ROUND(PLANFLOW, DBLE(1000)) .EQ.
C    1              ROUND(SIDECHAN, DBLE(1000))) THEN
C                   PLANFLOW = ROUND(SUPFLOW - 1, DBLE(1))
C                END IF
C                QFLOW = PLANFLOW
C             END IF
C
C             IF (ROUND(USSLIP(YEAR,MONTH), DBLE(1000)) .GT.
CC   1           ROUND(DBLE(582.9), DBLE(1000))) GOTO 10
C    1        ROUND(DBLE(177.94), DBLE(1000))) GOTO 10
CC            UNTIL USSLIP <= 177.94
C
CC              Adjust the gate setting
C
C             IF (ROUND(GATESET, DBLE(100)) .EQ. DBLE(0)) THEN
C                GATESET = DBLE(0.5)
C             END IF
C
C         =========================================================
C
C
C
C
C        Get the estimated MH mean for planned St Marys flow, using
C        the operational NBS's etc.
C
      DO 855, I = 1, 999
C
C        Route St Marys flow thru middle lakes, to get ave MH level
C
C        XX 40 increments hardcoded here
         CALL P77AMidLakeRoute ( SecsInPeriod,
     &      IDNINT(Round(StMarysFlow, CritFlowRound)) + 
     &      MichuronNWB, StClairNWB, ErieNWB, 
     &      Round(MicHuronBOPLevel,CritBOMRound),
     &      Round(StClairBOPLevel,CritBOMRound),
     &      Round(ErieBOPLevel,CritBOMRound),
     &      SCRivRet, DetRivRet, NiaRivRet,
     &      CritFlowRound, CritLevelRound, CritBOMRound, 40, 
     &      MicHuronEOPLevel, StClairEOPLevel, ErieEOPLevel, 
     &      MicHuronAveLevel)
C 
C
         IF( CritBMethod .EQ. 1) THEN
C           XX Old models round to 10 m3s here.  Go figure.
            tmplev = Round(USS ( IDNINT(Round(StMarysFlow, 1)), 
     &         MicHuronAveLevel), CritBOMRound)
         ELSE
C           XX If using new Crit B method, don't perpetuate the 
C           XX inconsistency.  Find US Slip level using the
C           XX GetLevelAtUSSlip function which, in turn, uses
C           XX the equation specified by the state of USSlipEquation
            tmplev = Round( 
     &         GetLevelAtUSSlip(StMarysFlow, MicHuronAveLevel,
     &            USSlipEquation, month, .TRUE.)
     &         , CritBOMRound)
         ENDIF 
C
         IF( SupVerbosity .GE. 5) THEN
            WRITE(CTMP,6547) OutputPrefix, tmplev, 
     &         MicHuronAveLevel, StMarysFlow, Gates, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF 
 6547    FORMAT(A,'US Slip= ',F7.3, ' with MH level= ', F7.3, 
     &      ' and St Marys flow= ', F6.1, ' (', F4.1,')', A1)
C
C
C        If not a problem with US Slip, jump out
         IF( tmplev .LE. 177.94D0 ) THEN
            IF( SupVerbosity .GE. 5) THEN
               WRITE(CTMP,6537) OutputPrefix, 'B', LLT
               CALL LogIt(LOGD,0,0)
            ENDIF 
            GOTO 9753
         ENDIF
C
C        Violated Criterion B!
         AdjustmentCodes(19:21) = 'CB '
         IF( SupVerbosity .GT. 1) THEN
            WRITE(CTMP,6548) OutputPrefix, tmplev, LLT
            CALL LogIt(LOGD,0,0)
 6548       FORMAT(A,'Criterion B violation!  US Slip level (',F7.3, 
     &         ') is greater than 177.94. Must close a gate . . .',A1)
         ENDIF
C             
C
         IF( CritBMethod .EQ. 2) THEN
C           Don't allow discharge less than pre-project flow
            IF( StMarysFlow .LE. Q1887 ) THEN
               IF( SupVerbosity .GT. 1) THEN
                  IF( I .EQ. 1 ) THEN 
                     CWRK(1:1) = ' '
                     L = 1
                  ELSE
                     CWRK = ' further '
                     L = 9
                  ENDIF
                  WRITE(*,6529)CWRK(1:L),MonthAbbrev(Month),Year,LLT
                  WRITE(CTMP,6529)CWRK(1:L),MonthAbbrev(Month),Year,LLT
                  CALL LogIt(LOGCDE,0,0)
                  WRITE(CTMP,6530) Q1887,StMarysFlow, Gates, LLT 
                  CALL LogIt(LOGCDE,0,0) 
               ENDIF
               GOTO 9753
            ENDIF
         ENDIF
 6529    FORMAT(' Unable to', A, 'reduce ',A,I6,
     &      ' Crit B outflow and remain above pre-project flow!', A1)
 6530    FORMAT(' Pre-project flow: ', F12.1, '   violating flow: ', 
     &      F12.1, '  violating gates: ', F5.2, A1)
C
         IF( NINT(Gates + 0.51) .EQ. 1 ) THEN
C           Only half a gate open!  Can't close any more!  Logic
C           undefined for existing practice, so just squawk and 
C           move on.  XX Maybe should cut sidechannel flows to 0?
            WRITE(CTMP,6549) MonthAbbrev(Month),Year, LLT
            CALL LogIt(LOGCDE,0,0)
 6549       FORMAT(' Unable to comply with Criterion B during ',A,
     &         I5,'!',A1)
            GOTO 9753
         ENDIF
C
C        Decrement gate setting
         Gates = Gates - 1  
         IF(NINT(Gates) .EQ. 0) Gates = 0.5
C
C        XX This may be possible in weird scenario with MH very high,
C        XX so should come up with logic someday to handle it.
         IF(Gates .LT. 0.0) STOP ' Ginger'
C
C        Get new mean outflow for adjusted gate setting, loop again.
C        Put it straight into StMarysFlow, if were at 1/2 gate 
C        shouldn't have gotten here anyway.  Be aware that this 
C        leaves a hole in the logic regarding minimum flow situations.
         CALL P77ASolveSup(SuperiorBOPLevel, DBLE(SuperiorNWB), 
     &      SecsInPeriod, Gates, CritFlowRound, CritBOMRound, 
     &      SuperiorEOPLevel, StMarysFlow, SuperiorAveLevel, Month)
C
C        Knock accuracy down to cm or whatever
  855    SuperiorEOPLevel = Round( SuperiorEOPLevel, CritBOMRound )
C
C     Get here only by falling out of loop.  How?
      STOP ' Mary Ann' 
C
C
 9753 CONTINUE 
C
C
C     Sidechannel (i.e., non-gated) flow determination
C 
C     Compute non-gated flow for the month.  Need to determine GateFlow,
C     and subtract it from the nominal total outflow.  This requires a
C     SW Pier level, computed from the current estimate of the Superior 
C     monthly mean level.  Existing practice was to perform these  
C     calculations by hand, so all levels are rounded to cm.  
C
C     This is set up for iterative solution of simultaneous equations in
C     if we want to use the new EC equation for LS:SWP.  The equations are 
C     simultaneous because SWPier level is a function of both Superior
C     level and outflow, while outflow is a function of Superior level.
C     So solve for SWP, get NonGatedFlow, get SWP, get NonGatedFlow, etc
C     until converge.  
C
C     Traditional practice uses the old SWP:LS relationship, which is not
C     a function of outflow, so there is no need to solve simultaneous 
C     equations.  When the old SWP:LS relationship is selected, the 
C     loop cycles twice (need at least two iterations to test for 
C     convergence), computing the same answer both times.
C              
C     Seed NonGatedFlow with some number
      NonGatedFlow = 2320
C
C
C
      IF(SupVerbosity .GE. 6 ) THEN
         WRITE(CTMP,'(A, A1)') ' Solving for sidechannel flow . . .',LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
      L = -9999
      Flow2ItsAgo = L
C
      DO 4778, I=1,999
         CALL GetSWPier( SuperiorAveLevel, Month, Gates, 
     &      IDNINT(NonGatedFlow), CritBOMRound, LStoSWP, 5, Year,
     &      'SupCriterions'//LLT, SWPier, gatedflow)
C
C        This is the method that Matt proposes, not the way it
C        is calculated in "existing practices". 
         NonGatedFlow = IDNINT(StMarysFlow - gatedflow)
C        XX Backward compat fix:
C        XX This following line replicates the goofy mixed rounding of
C        XX of existing practices.  If we're rounding flows, we'll do
C        XX it the old fashioned way. Remove it when allowed.
         IF(CritFlowRound .NE. -99 ) NonGatedFlow = 
     &      Round(StMarysFlow,1) - Round(gatedflow,CritFlowRound)
C        XX End of fix block
C
         IF( L .EQ. IDINT(gatedflow * 100000.D0) ) GOTO 4779
C        Check to see if oscillating between two values 
         IF(CritBOMRound .NE. -99 )THEN
            IF(IDINT(gatedflow*100000.D0) .EQ. Flow2ItsAgo) GOTO 4779
            Flow2ItsAgo = L
         ENDIF
 4778    L = IDINT(gatedflow * 100000.D0)
C
 4779 IF(SupVerbosity .GE. 3 ) THEN
         WRITE( CTMP, 2310) NonGatedFlow,LLT
         CALL LogIt(LOGD,0,0)
         WRITE( CTMP, 2309) StMarysFlow, SuperiorAveLevel, 
     &      SWPier, gatedflow, LLT
         CALL LogIt(LOGD,0,1)
      ENDIF
C
      WRITE(SupSideChanFlow,'(I6)') 
     &   IDNINT(Round(NonGatedFlow,CritFlowRound))
      IF(SupVerbosity .GE. 7 ) THEN
         WRITE( CTMP, '(A,A1)') ' finished SupCriterions',LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
 2309 FORMAT( '   (based on nominal total flow of ', F9.4, 
     &   ' Sup/SW Pier levels of ', 2F9.4,', and gate flow of', F9.4, 
     &   ')', A1)
 2310 FORMAT( '   Computed "Official" non-gated flow as ', F9.4, A1)
C
 6550 FORMAT(A, A3,I5,
     &   ' Criterion C violation!  Will try to close a gate . . .',A1)
 6551 FORMAT(A, A3,I5,
     &   ' Cannot close any more gates, setting regulated flow to',
     &   F8.2, ' m3s', A1)
 6537 FORMAT(A, ' . . . passed Criterion ', A1, ' check OK.', A1)
      END
C----&------------------------------------------------------------------
C  -----------------------------------------------------------------------
C
C
      SUBROUTINE P55MR49(Month, Year, PrevSupMean, SPBOP, P55MR49Q, 
     &                   P55Gates, ComputedNonGate)
C     
C     This data is interpreted from "A Report on Lake Superior 
C     Regulation Operational Guides to the 1955 Modified Rule of 1949" 
C     dated January 1973. 
C
C     Underlying this code, as it stands now, is the assumption that 
C     this regulation plan will only be run to comply with the Orders
C     of Approval that the plan "will not exceed xxx.xx more than the
C     55 Modified Rule of 49".
C 
C     NOTE: If the first month of the simulation is January-April, this 
C     code will assume the previous winter months had the same gate 
C     setting.
C
C     (13 September 2010) Although Criterion C (limiting outflow to pre-
C     project when Superior is below 183.40m) will only rarely be an issue 
C     for the '55 Modified Rule of '49, there is now an explicit check for 
C     Criterion C hardcoded into the logic.
C     NOTE: This is coded differently than the current practice for Plan 77A.
C     Plan77A currently selects a proposed Lake Superior outflow and then 
C     estimates the coming month's mean level (based on an average supply) 
C     and uses this level for checking Criterion C.  The method used here 
C     looks at the previous month's mean level which is the same method that 
C     was used historically for making the regulation decision.  It is 
C     believed that this method more accurately represents the method 
C     actually used for performing Criterion C checks with the '55 Modified 
C     Rule of '49.  
C     Pre-Project outflow here is determined using the current (end of month) 
C     level.  It is unclear how this lines up with historical practices for 
C     the '55 Modified Rule of '49, although elsewhere in this code there is 
C     a statement with regard to Plan 77A that "Pre-existing practice is to 
C     estimate the pre-project flow using the mean monthly Superior level 
C     associated with the outflow from Plan 1977A/Criterion B."
C
C     (26 December 2009) Changed to using previous month's Superior mean 
C     level to determine regulation for coming month.  This included 
C     initializing the mean level at the start of the simulation 
C     (control.for: ~1294) and making sure it gets stored every regulation 
C     interval for use by this routine (superior.for: ~500).
C
C     (30 November 2009) Modified code to correctly assign 16 gates plus 
C     65,000 cfs for maximum flow months.
C
C     (27 October 2009) This subroutine now includes logic for limiting gate
C     changes in winter, based on the following information:
C     A program that was used to run this regulation plan for the 1977A 
C     development was obtained from Rob Caldwell & David Fay of Environment 
C     Canada.  This code may not represent the actual, historical procedure 
C     properly, given that it was developed 15+ years after the plan was last 
C     used.  However, a careful comparison of actual winter gate changes 
C     during the time period this regulation plan was in effect show little 
C     internal consistency (and none with any interpretation of the 
C     regulation plan).  Therefore, the logic used in the code from Plan 
C     1977A is implemented here to at least maintain consistency between 
C     regulation plan studies.
C
C -----Original Message-----
C From: Dahl, Travis A LRE 
C Sent: Wednesday, October 28, 2009 12:29 PM (Davis)
C To: Fay,David [Ontario]; btolson@civmail.uwaterloo.ca; Watkins, Dave W HEC 
C Visiting Scholar; Bill Werick; dhburn@civmail.uwaterloo.ca; Fan,Yin 
C [Ontario]; Leger,Wendy [Burlington]; McPherson, Matthew M HEC; Kropfreiter,
C Melissa A LRE; Thomas, Richard J LRB; Rob.Caldwell@ec.gc.ca; Thieme, Scott 
C J LRE; moins@ottawa.ijc.org; Eberhardt, Anthony J IWR
C Cc: Allis, John T LRE; Koschik, John A LRE
C Subject: RE: Plan formulation status
C
C All,
C   I have completed debugging the code for the '55 Modified Rule of '49.  In
C coding the winter gate change logic, I started by looking at the actual 
C winter gate changes during the time period this regulation plan was in 
C effect and found no consistency with any interpretation of the rules based 
C on the existing documentation.  Therefore, the implemented procedure for 
C winter gate changes is based on the code provided by David & Rob that had,
C presumably, been used for Plan 77A development.
C   As a check, I confirmed that Plan 77A met Criteria "a" (“no greater 
C probability of of exceeding elevation (183.86 m IGLD 1985) than would have 
C occurred using the 1955 Modified Rule of 1949”) for the 1900-2006 
C historical dataset.  (Plan 77A had 0 violations while the '55 Modified 
C Rule of '49 had 7.)
C   Documentation for the plan still needs to be created so that future 
C studies will not be as hampered as we were.  This is still on my plate.
C     -Travis
C   

      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER          Month, Year
C                      (input) The month number and year to compute 
C                      outflow for.
C
      DOUBLE PRECISION PrevSupMean
C                      (input) The Superior mean level of previous month
      DOUBLE PRECISION SPBOP
C                      (input) The Superior beginning of month level 
C                      (used by Gate setting routine)
C
      DOUBLE PRECISION P55MR49Q
C                      (output) The plan flow 
C
      REAL             P55Gates
C                      (output) The plan gates open 
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER          I
C                      Scratch variables
C
      CHARACTER*(*)    IDString
C                      Used as prefix on lines of detailed.log 
C
      REAL             BoMayS(6), BoJunS(6), BoJulS(6), BoAugS(6), 
     &                 BoSepS(6), BoOctS(6), BoNovS(6), BoDecS(4), 
     &                 BoJanS(4), BoFebS(4), BoMarS(4), BoAprS(4)   
      INTEGER          BoMayQ(6), BoJunQ(6), BoJulQ(6), BoAugQ(6), 
     &                 BoSepQ(6), BoOctQ(6), BoNovQ(6), BoDecQ(4), 
     &                 BoJanQ(4), BoFebQ(4), BoMarQ(4), BoAprQ(4)   
      DOUBLE PRECISION PrevSup55ft
C
      DOUBLE PRECISION NominalStMarys
      DOUBLE PRECISION ComputedNonGate, SWPierLev, GateFlow
C
      DOUBLE PRECISION SupGates_LastMonth
C                   This variable is used to track the previous 
C                   month's winter gate settings
      DOUBLE PRECISION AvgSupNBS(12)
C        This variable contains the long-term average (1900-2006) monthly
C        supplies for Lake Superior and functions as a crude 50% forecast to
C        be used when determining maximum Criterion C outflow.
      CHARACTER*(*)    OutputPrefix
C                      Used as prefix on lines of detailed.log 
      DOUBLE PRECISION Q1887
C                      Pre-project flow
      CHARACTER    SupPreProject*6, BalancingFlow*8, SupGateFlow*6,
     &             SupPlanFlow*6, SupBOMLevel*8, SWPierLevel*8,
     &             SupSideChanFlow*6, SPRegCodes*40, P77AHeader*120,
     &             P2012Header
C                  Lake Superior regulation output strings
      DOUBLE PRECISION Round
C                      Function for rounding.

C
C     Declarations regarding COMMON block variables:
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
      INTEGER      CritBOMRound, CritFlowRound, CritLevelRound, 
     &             Q1887Method, CritBMethod
      DOUBLE PRECISION CritCThresholdEl
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
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /MONAME/   MonthAbbrev
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /P55M49/   SupGates_LastMonth
      COMMON /SUPCRT/   CritBOMRound, CritFlowRound, CritLevelRound,
     &                  Q1887Method, CritBMethod, CritCThresholdEl
      COMMON /SUPOUT/   SupPreProject, BalancingFlow, SupGateFlow,
     &                  SupPlanFlow, SupBOMLevel, SWPierLevel, 
     &                  SupSideChanFlow, SPRegCodes, P77AHeader,
     &                  P2012Header
C
      PARAMETER (IDString='    Plan 55MR49==> ')
      PARAMETER (OutputPrefix = 'Plan 55MR49: ')
      SPRegCodes = '                                        '
C                      (input/output) String containing code letters 
C                      describing what adjustments were made.

C
      DATA BoMayS / 600.44, 600.31, 600.16, 600.04, 599.89, 599.75 /
      DATA BoMayQ / 98000, 89000, 80000, 70000, 67000, 58000 /
      DATA BoJunS / 600.62, 600.51, 600.40, 600.27, 600.15, 600.04 /
      DATA BoJunQ / 98000, 89000, 80000, 70000, 67000, 58000 /
      DATA BoJulS / 600.78, 600.69, 600.62, 600.53, 600.41, 600.30 /
      DATA BoJulQ / 100000, 90000, 81000, 70000, 68000, 58000 /
      DATA BoAugS / 601.00, 600.91, 600.81, 600.73, 600.61, 600.49 /
      DATA BoAugQ / 103000, 90000, 81000, 70000, 68000, 58000 /
      DATA BoSepS / 601.09, 600.99, 600.89, 600.80, 600.70, 600.60 /
      DATA BoSepQ / 103000, 91000, 81000, 70000, 68000, 58000 /
      DATA BoOctS / 601.16, 601.05, 600.94, 600.84, 600.72, 600.60 /
      DATA BoOctQ / 104000, 91000, 81000, 70000, 68000, 58000 /
      DATA BoNovS / 601.09, 600.95, 600.83, 600.70, 600.58, 600.46 /
      DATA BoNovQ / 103000, 90000, 81000, 70000, 68000, 58000 /
C
      DATA BoDecS / 600.90, 600.70, 600.51, 600.31 /
      DATA BoDecQ / 76000, 70000, 68000, 55000 /
      DATA BoJanS / 600.71, 600.53, 600.33, 600.13 /
      DATA BoJanQ / 76000, 70000, 67000, 55000 /
      DATA BoFebS / 600.50, 600.31, 600.11, 599.92 /
      DATA BoFebQ / 75000, 70000, 67000, 55000 /
      DATA BoMarS / 600.41, 600.22, 599.96, 599.70 /
      DATA BoMarQ / 75000, 70000, 67000, 55000 /
      DATA BoAprS / 600.34, 600.11, 599.86, 599.61 /
      DATA BoAprQ / 75000, 70000, 67000, 55000 /      
C      
      DATA AvgSupNBS / -330.0, 130.0, 1320, 4280, 5060, 4400, 3600, 
     &                2650, 1970, 1170, 480, -570 /
C 
C
C     Log the data entering Plan  
C
      IF( SupVerbosity .GE. 5) THEN
         WRITE(CTMP, '(A,A3,I6,A,F8.3,A1)') IDString,MonthAbbrev(Month),
     &    Year, '  previous Superior mean level=', PrevSupMean, LLT
         CALL LogIt(LOGD,1,0)
      END IF
C
C     0.366 is lakewide adjustment from IGLD 85 IGLD 55
      PrevSup55ft = (PrevSupMean - 0.366) / 0.3048

CXX   Check to see how this affects the even-zero hundredths in table
C      PrevSup55ft = Round (PrevSup55ft, 2)

      P55MR49Q = 0.0D0
      P55Gates = 0.0
C
C
      IF( Month .EQ. 1) THEN 
         IF( PrevSup55ft .GT. BoJanS(1)) THEN
               P55MR49Q = 85000D0
         ELSE
            DO 5501, I = 2, 4
C
               IF( PrevSup55ft .GT. BoJanS(I)) THEN
                  P55MR49Q = BoJanQ(I-1)
                  GOTO 7777
               END IF
            P55MR49Q = BoJanQ(4)
 5501       CONTINUE
         END IF
C
      ELSE IF( Month .EQ. 2) THEN 
         IF( PrevSup55ft .GT. BoFebS(1)) THEN
               P55MR49Q = 85000D0
         ELSE
            DO 5502, I = 2, 4
C
               IF( PrevSup55ft .GT. BoFebS(I)) THEN
                  P55MR49Q = BoFebQ(I-1)
                  GOTO 7777
               END IF
            P55MR49Q = BoFebQ(4)
 5502       CONTINUE
         END IF
C
      ELSE IF( Month .EQ. 3) THEN 
         IF( PrevSup55ft .GT. BoMarS(1)) THEN
               P55MR49Q = 85000D0
         ELSE
            DO 5503, I = 2, 4
C
               IF( PrevSup55ft .GT. BoMarS(I)) THEN
                  P55MR49Q = BoMarQ(I-1)
                  GOTO 7777
               END IF
            P55MR49Q = BoMarQ(4)
 5503       CONTINUE
         END IF
C
      ELSE IF( Month .EQ. 4) THEN 
         IF( PrevSup55ft .GT. BoAprS(1)) THEN
               P55MR49Q = 85000D0
         ELSE
            DO 5504, I = 2, 4
C
               IF( PrevSup55ft .GT. BoAprS(I)) THEN
                  P55MR49Q = BoAprQ(I-1)
                  GOTO 7777
               END IF
            P55MR49Q = BoAprQ(4)
 5504       CONTINUE
         END IF
C
      ELSE IF( Month .EQ. 5) THEN 
         IF( PrevSup55ft .GT. BoMayS(1)) THEN
C              XX Will need to make sure that calling routine know to use 16
C                 gates + 65k
C              999997D0 used as a flag to indicate 16 Gates + 65k cfs 
               P55MR49Q = 999997D0
         ELSE
            DO 5505, I = 2, 6
C
               IF( PrevSup55ft .GT. BoMayS(I)) THEN
                  P55MR49Q = BoMayQ(I-1)
                  GOTO 7777
               END IF
            P55MR49Q = BoMayQ(6)
 5505       CONTINUE
         END IF
C
      ELSE IF( Month .EQ. 6) THEN 
         IF( PrevSup55ft .GT. BoJunS(1)) THEN
C              999997D0 used as a flag to indicate 16 Gates + 65k cfs 
               P55MR49Q = 999997D0
         ELSE
            DO 5506, I = 2, 6
C
               IF( PrevSup55ft .GT. BoJunS(I)) THEN
                  P55MR49Q = BoJunQ(I-1)
                  GOTO 7777
               END IF
            P55MR49Q = BoJunQ(6)
 5506       CONTINUE
         END IF
C
      ELSE IF( Month .EQ. 7) THEN 
         IF( PrevSup55ft .GT. BoJulS(1)) THEN
C              999997D0 used as a flag to indicate 16 Gates + 65k cfs 
               P55MR49Q = 999997D0
         ELSE
            DO 5507, I = 2, 6
C
               IF( PrevSup55ft .GT. BoJulS(I)) THEN
                  P55MR49Q = BoJulQ(I-1)
                  GOTO 7777
               END IF
            P55MR49Q = BoJulQ(6)
 5507       CONTINUE
         END IF
C
      ELSE IF( Month .EQ. 8) THEN 
         IF( PrevSup55ft .GT. BoAugS(1)) THEN
C              999997D0 used as a flag to indicate 16 Gates + 65k cfs 
               P55MR49Q = 999997D0
         ELSE
            DO 5508, I = 2, 6
C
               IF( PrevSup55ft .GT. BoAugS(I)) THEN
                  P55MR49Q = BoAugQ(I-1)
                  GOTO 7777
               END IF
            P55MR49Q = BoAugQ(6)
 5508       CONTINUE
         END IF
C
      ELSE IF( Month .EQ. 9) THEN 
         IF( PrevSup55ft .GT. BoSepS(1)) THEN
C              999997D0 used as a flag to indicate 16 Gates + 65k cfs 
               P55MR49Q = 999997D0
         ELSE
            DO 5509, I = 2, 6
C
               IF( PrevSup55ft .GT. BoSepS(I)) THEN
                  P55MR49Q = BoSepQ(I-1)
                  GOTO 7777
               END IF
            P55MR49Q = BoSepQ(6)
 5509       CONTINUE
         END IF
C
      ELSE IF( Month .EQ. 10) THEN 
         IF( PrevSup55ft .GT. BoOctS(1)) THEN
C              999997D0 used as a flag to indicate 16 Gates + 65k cfs 
               P55MR49Q = 999997D0
         ELSE
            DO 5510, I = 2, 6
C
               IF( PrevSup55ft .GT. BoOctS(I)) THEN
                  P55MR49Q = BoOctQ(I-1)
                  GOTO 7777
               END IF
            P55MR49Q = BoOctQ(6)
 5510       CONTINUE
         END IF
C
      ELSE IF( Month .EQ. 11) THEN 
         IF( PrevSup55ft .GT. BoNovS(1)) THEN
C              999997D0 used as a flag to indicate 16 Gates + 65k cfs 
               P55MR49Q = 999997D0
         ELSE
            DO 5511, I = 2, 6
C
               IF( PrevSup55ft .GT. BoNovS(I)) THEN
                  P55MR49Q = BoNovQ(I-1)
                  GOTO 7777
               END IF
            P55MR49Q = BoNovQ(6)
 5511       CONTINUE
         END IF
C
      ELSE IF( Month .EQ. 12) THEN 
         IF( PrevSup55ft .GT. BoMayS(1)) THEN
               P55MR49Q = 85000D0
         ELSE
            DO 5512, I = 2, 4
C
               IF( PrevSup55ft .GT. BoDecS(I)) THEN
                  P55MR49Q = BoDecQ(I-1)
                  GOTO 7777
               END IF
            P55MR49Q = BoDecQ(4)
 5512       CONTINUE
         END IF
      END IF
      
 7777 CONTINUE

      IF (P55MR49Q .EQ. 999997D0) THEN
C        Set flow to 16 gates + 65tcfs      
         P55Gates = 16
C        65,000 cfs = 1840 cms
C        XX Verify the 65tcfs number, since this does not equal the 
C        XX current maximum side channel flow
         ComputedNonGate = 1840
C        XX NominalStMarys?  Use GetSWPier?
C            CALL GetSWPier(SupLevel, MonthNum, SupGates, 
C     &         IDNINT(TmpFlow),-99,LStoSWP,7,-9999,'GetInstantGate()',
C     &         SWPierLev,GateFlow)

      ELSE
C        Convert from CFS to CMS     
         P55MR49Q = (0.3048D0**3D0) * P55MR49Q
C        Calculate gate setting/Nominal St Marys flow
         CALL GetInstantGate(SPBOP, P55MR49Q, Month, 0, 
     &      SupGates_LastMonth, P55Gates, GateFlow, SWPierLev, 
     &      ComputedNonGate, NominalStMarys)
      END IF      
C
C     Criterion C Check
C     (Outflows capped at Pre-Project if Lake Superior < 183.40)
C     NOTE: This was not explicitly part of the original plan and was likely 
C     checked by hand outside of the regulation plan.  
C     This implementation is different from the current practice used for 
C     Plan 77A.  Rather than estimating what next month's mean level might be
C     with the proposed outflow, it checks last month's mean level aginst the
C     threshold.  Pre-Project outflow is determined using the BOP level.
      IF (PrevSupMean .LT. 183.40) THEN
         IF( SupVerbosity .GE. 3) THEN
            WRITE(CTMP,'(A, A3,I5, A, F15.10, A, A1)')
     &         '   Must consider Crit C for ', MonthAbbrev(Month), Year,
     &         ' because Superior level for Crit C test (', 
     &         PrevSupMean, 'm) is below threshold of 183.40m.', LLT
            CALL LogIt(LOGD,0,0)
         END IF 
C        Estimate Superior's mean level next month if pre-project flow is 
C        used. This is based on the long-term (1900-2008) average supply.
         
C        NOTE: This assumes that Q1887Method is "1" (Use expected mean level
C           of Lake Superior rather than route)
         Q1887= Round(824.7D0*(SPBOP-181.43D0)**1.5D0,
     &      CritFlowRound)
C
         IF( SupVerbosity .GE. 4) THEN
            IF (Month .GT. 1) THEN
               WRITE(CTMP,
     &            '(A, A3,I5, A, F15.10,A1)')' Pre-project flow for ', 
     &            MonthAbbrev(Month-1), Year, ': ',  Q1887, LLT
            ELSE
               WRITE(CTMP,
     &            '(A, A3,I5, A, F15.10,A1)')' Pre-project flow for ', 
     &            MonthAbbrev(12), Year, ': ',  Q1887, LLT
            END IF
            CALL LogIt(LOGD,0,0)
         END IF 
         WRITE(SupPreProject, '(I6)') IDNINT(Q1887) 
C        Now Check to see if the proposed flow is less than Pre-Project.  If 
C        it isn't, reduce the flow...
         IF(NominalStMarys .GT. Q1887) THEN
C           Violated Criterion C!
C           Write Violation Code to log
            SPRegCodes(22:24) = 'CC '
            IF( SupVerbosity .GT. 2) THEN
               WRITE(CTMP,6553)OutputPrefix,MonthAbbrev(Month),Year,LLT
               CALL LogIt(LOGD,0,0)
            END IF
C           Set flow to Pre-Project (1887 Equation)
            P55MR49Q = Q1887
C           Calculate new gate setting
            CALL GetInstantGate(SPBOP,P55MR49Q,Month,0, 
     &         SupGates_LastMonth, P55Gates, GateFlow, SWPierLev, 
     &         ComputedNonGate, NominalStMarys)
         END IF
      END IF
C
C XX NOTE: This section should never be triggered because changes in the code
C XX(in GetInstantGate) prevent gate changes during the "winter".  This means 
C XX that there could be small differences between the implementation in this 
C XX version of the code and the documented version of Plan 55MR49.  Testing 
C XX with the 13 NBS sequences used for the IUGLS, however, showed no 
C XX differences.
C     Check against winter flow gate change constraints
      IF (Month .LE. 4) THEN
C        Once gates are set during December regulation, they can only
C        be changed during winter (December-April) under certain conditions
C        The logic used for winter gate changes in Plan77A is as follows:
C        - Based on the proposed gate setting, call Qsup.
C        - If a gate change is proposed
C           - IF 0.5<= proposed gate setting <=2
C              - IF 0.5<= previous gate setting <=2
C                 - Both old and new gate settings are between 0.5 and 2
C                 - Call Qsup again
C                 - Don't change gates!
C              - ELSE
C                 - Change gates
C           - ELSE IF proposed gate setting > 2
C              - IF previous gate setting  > 2
C                 - Both old and new gate settings are greater than 2
C                 - Call Qsup again
C                 - Don't change gates!
C              - ELSE
C                 - Change gates
         IF (SupGates_LastMonth .GE. 0) THEN
C           Gate setting present for last month, so proceed with checking to 
C           see if a change is allowed.  This code assumes that whatever gate 
C           setting is proposed for the first month of a simulation is the 
C           same as month before it and is therefore allowed.
C        - If a gate change is proposed
           IF (P55Gates .NE. SupGates_LastMonth) THEN
C             - IF 0.5<= proposed gate setting <=2
               IF ((P55Gates .GE. 0.5) .AND. (P55Gates .LE. 2)) THEN
C                 - IF 0.5<= previous gate setting <=2
                  IF ((SupGates_LastMonth .GE. 0.5) .AND. 
     &               (SupGates_LastMonth .LE. 2)) THEN
C                    - Both old and new gate settings are between 0.5 and 2
C                    - Don't change gates!
                     P55Gates = SupGates_LastMonth
C                    NOTE: This will clobber the gate setting, potentially 
C                    leading to significant discrepancies between the 
C                    proposed and actual flows.
                  ELSE
C Winter Gate change called for by plan, but verboten under current (2011) practice.  
C Die here, in order to ensure that this is noticed.
C XX Ideally, this would just get flagged as a violation code ('WG'?)
                     STOP 'Winter gate change called for by 55MR49!'
                  END IF
C              - ELSE IF proposed gate setting > 2
               ELSE IF (P55Gates .GT. 2) THEN
C                 - IF previous gate setting  > 2
                  IF (SupGates_LastMonth .GT. 2) THEN
C                    - Both old and new gate settings are greater than 2
C                    - Don't change gates!
                     P55Gates = SupGates_LastMonth
C                    NOTE: This will clobber the gate setting, potentially 
C                    leading to significant discrepancies between the 
C                    proposed and actual flows.
                  ELSE
C Winter Gate change called for by plan, but verboten under current (2011) practice.  
C Die here, in order to ensure that this is noticed.
C XX Ideally, this would just get flagged as a violation code ('WG'?)
                     STOP 'Winter gate change called for by 55MR49!'
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C     This bit of code will hopefully keep track of the previous 
C     month's gate setting, for use with winter month gate changes
      SupGates_LastMonth = P55Gates
C
      WRITE( CTMP, '(A,F5.1,A,F9.1,A,F9.1,A1)') 
     &   ' Following 55 Modified Rule of 49, Gates: ', P55Gates,
     &   ' Side Channel: ', ComputedNonGate, ' & St Marys Flow: ',
     &   P55MR49Q, LLT
      CALL LogIt(LOGD,1,1)
C
      IF( SupVerbosity .GE. 5) THEN
         WRITE(CTMP, '(A,A3,I6,A,F8.3,A,F8.1, A1)') IDString,
     &    MonthAbbrev(Month), Year, 
     &    '  previous Superior mean level=', PrevSup55ft, 
     &    " plan flow=", P55MR49Q, LLT
         CALL LogIt(LOGD,1,0)
      END IF
C
 6553 FORMAT(A, A3,I5,
     &   ' Criterion C violation!  Will try to close a gate . . .',A1)
C
      END
C
C     ---------------------------------------------------------------------------
      SUBROUTINE P2012_Criterions( Year, Month, SuperiorNWB,
     &   MichuronNWB, StClairNWB, ErieNWB, SCRivRet, DetRivRet,
     &   NiaRivRet, SuperiorBOPLevel, MicHuronBOPLevel,StClairBOPLevel,
     &   ErieBOPLevel, LStoSWP, SecsInPeriod, StMarysFlow, SWPier,
     &   NonGatedFlow, Gates, AdjustmentCodes, tmplev)
C
C     This function takes a preliminary flow from Plan 2012 and 
C     checks/modifies the flow according to Criterion B of the 
C     Orders of Approval (at least according to the operational 
C     practices expected to be implemented in ~September 2013).
C
C          Criterion A:  (not applicable in this context)
C          Criterion B:  Prevent flows greater than pre-project when levels at US Slip greater than 177.94 m
C          Criterion C:  (Checked during plan development for the 
C                         International Upper Great Lakes Study)
C
C
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER          Year, Month
C                      (input) The date of the outflow.
C
      INTEGER          SuperiorNWB
C                      (input) The net water balance (NBS + LLO) to 
C                      use for Lake Superior change-in-storage 
C                      calculations.  Could be determined in three 
C                      different ways, depending on value of P77ACompat.
      INTEGER          MichuronNWB
C                      (input) The net water balance (NBS + ChiDiv), 
C                      calculated according to value of P77ACompat.
      INTEGER          StClairNWB
C                      (input) The net water balance (NBS), 
C                      calculated according to value of P77ACompat.
      INTEGER          ErieNWB
C                      (input) The net water balance (NBS - Welland), 
C                      calculated according to value of P77ACompat.
      INTEGER          SCRivRet, DetRivRet, NiaRivRet
C                      (input) Flow retardation values to use, 
C                      calculated according to value of P77ACompat.
C
      DOUBLE PRECISION SuperiorBOPLevel, MicHuronBOPLevel, 
     &                 StClairBOPLevel, ErieBOPLevel
C                      (input) Beginning-of-Period lake levels 
C    
      INTEGER          LStoSWP
C                      (input) Which Lake Superior to SW Pier 
C                      relationship to use.
C
      DOUBLE PRECISION SecsInPeriod     
C                      (input) The length of the routing interval.
C
      DOUBLE PRECISION StMarysFlow
C                      (input/output) The Superior outflow for the 
C                      routing interval.  Initially given as some kind
C                      of target, either from the balancing equation or
C                      forecast averaging.  This routine adjusts and 
C                      modifies this flow according to Criterions B and
C                      C, and returns the adjusted value in StMarysFlow.
C
      DOUBLE PRECISION SWPier
C                      (output) The SWPier level as the old models would
C                      compute it, which yields a "gated flow" which
C                      yields the 'sidechannel' flow.
      DOUBLE PRECISION NonGatedFlow
C                      (input/output) the sidechannel flow
C
      REAL             Gates
C                      (input/output) The number of gates open in Month.
C
      CHARACTER*(40)   AdjustmentCodes
C                      (input/output) String containing code letters 
C                      describing what adjustments were made.
C    
C
C                       
C
C
C     Declarations regarding local variables or functions:
C
      DOUBLE PRECISION tmplev, gatedflow
      INTEGER          I, L, Flow2ItsAgo
C                      Scratch variables
C
      DOUBLE PRECISION SuperiorAveLevel, MicHuronAveLevel
C                      Mean levels for the routing period, as computed
C                      within Plan 1977A.  NOT necessarily same results
C                      as routing through these lakes after Plan 1977A 
C                      outflow is determined.
C
      DOUBLE PRECISION SuperiorEOPLevel, MicHuronEOPLevel, 
     &                 StClairEOPLevel, ErieEOPLevel
C                      End-of-Period lake levels 
C    
C
      CHARACTER*(*)    OutputPrefix
C                      Used as prefix on lines of detailed.log 
C
      DOUBLE PRECISION Round
C                      Function for rounding.
C
      DOUBLE PRECISION USS
C                      Function that estimates level at US Slip based
C                      on St Marys flow and period average level of
C                      Michigan-Huron
C
      DOUBLE PRECISION GetLevelAtUSSlip 
C                   Returns US Slip level corresponding to given 
C                   St Marys flow and MH elevation
C
      DOUBLE PRECISION Q1887
C                      Pre-project flow
      DOUBLE PRECISION MeanTotalFlow
C                      Total flow in St. Marys River returned from SolveGatesSide
C
C
C
C     
C
C     Declarations regarding COMMON block variables:
C
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
      INTEGER      MonthLength
C                  1 if uniform, 2 if actual
C
      LOGICAL      LakesToSolve(8)
      INTEGER      StartYear, StartMonth, StartDay,
     &             EndYear, EndMonth, EndDay, IntervalNameLen(4)
C                  Just declared because in same COMMON as MonthLength
C
      CHARACTER    SupPreProject*6, BalancingFlow*8, SupGateFlow*6,
     &             SupPlanFlow*6, SupBOMLevel*8, SWPierLevel*8,
     &             SupSideChanFlow*6, SPRegCodes*40, P77AHeader*120,
     &             P2012Header
C                  Lake Superior regulation output strings
      INTEGER      CritBOMRound, CritFlowRound, CritLevelRound, 
     &                  Q1887Method, CritBMethod
      INTEGER      USSlipEquation
C                  1=ancient eq (default), 2=2010 EC eq
C
      DOUBLE PRECISION CritCThresholdEl, USSIceWeedAdj
C
C
C
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /MONAME/   MonthAbbrev
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /RCRUNI/   LakesToSolve, 
     &                  IntervalNameLen, MonthLength,
     &                  StartYear, StartMonth, StartDay, 
     &                  EndYear, EndMonth, EndDay
      COMMON /SUPOUT/   SupPreProject, BalancingFlow, SupGateFlow,
     &                  SupPlanFlow, SupBOMLevel, SWPierLevel, 
     &                  SupSideChanFlow, SPRegCodes, P77AHeader,
     &                  P2012Header
      COMMON /SUPCRT/   CritBOMRound, CritFlowRound, CritLevelRound,
     &                  Q1887Method, CritBMethod, CritCThresholdEl
      COMMON /USSLIP/   USSIceWeedAdj, USSlipEquation
C
C
      PARAMETER (OutputPrefix = 'Superior Criterion B check:')
C 
C
C     Log the entering data
C
      IF( SupVerbosity .GE. 2) THEN
         WRITE(CTMP, 6546) OutputPrefix, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6542) OutputPrefix,MonthAbbrev(Month),Year, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6544) OutputPrefix, StMarysFlow, Gates,
     &      NonGatedFlow, AdjustmentCodes, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6545) OutputPrefix, SuperiorBOPLevel, 
     &      MicHuronBOPLevel, StClairBOPLevel, ErieBOPLevel, 
     &      SecsInPeriod, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6541) OutputPrefix, SuperiorNWB, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6543) OutputPrefix, MichuronNWB,
     &      StClairNWB, ErieNWB, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP, 6540) OutputPrefix, SCRivRet, DetRivRet, 
     &      NiaRivRet, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
 6541 FORMAT(A,'Superior NTS=', I6, A1) 
 6542 FORMAT(A,'Evaluating plan flow for ',A3,I5,A1) 
 6543 FORMAT(A,'MH/SC/ER net water balances: ',3I7,A1) 
 6540 FORMAT(A,'MH/SC/ER outflow retardations: ',3I7,A1) 
 6544 FORMAT(A,'Plan flow=', F9.1, ' (', F5.2, ' gates)  sidech=', 
     &   F7.2,'  AdjustmentCodes=', A, A1) 
 6545 FORMAT(A,'BOP Levels=', 4F8.3, ' dT=', F14.3, A1) 
 6546 FORMAT(A,'Checking preliminary flow against level at US Slip:',A1)
C
C
      SuperiorBOPLevel = Round( SuperiorBOPLevel, CritBOMRound )
C     The next two bits of code (SolveGatesSide and assignment of Q1887 value)
C     are necessary in preperation for the Criterion B calculations.
C
C     XX TODO: see if this flag for determining whether full sidechannel capacity,
C     XX causes logic to be dependent on using Plan 1977A 
C     15Oct2010 Yes - it modifies plan flow if fails to find 'RJ', which means 
C     need to fix this some way for other plans to use routine to do crit c
C     checks without inadvertently changing plan flow
      CALL SolveGatesSide( Gates, IDNINT(NonGatedFlow), Month, 
     &   SuperiorBOPLevel, 
     &   DBLE(SuperiorNWB), SecsInPeriod, LStoSWP, SuperiorEOPLevel, 
     &   SuperiorAveLevel, SWPier, gatedflow, MeanTotalFlow)      
C     XX TODO: create logic to punt and return some flow if elev too 
C     XX low for 1887 eqn
      IF(SuperiorAveLevel .LT. 181.43D0) STOP ' Mr. Howell II!'
C
      IF( Q1887Method .EQ. 1) THEN
         Q1887= Round(824.7D0*(SuperiorAveLevel-181.43D0)**1.5D0,
     &      CritFlowRound)
      ELSE IF( Q1887Method .EQ. 2) THEN
         CALL GetSupPreProject( SuperiorBOPLevel, DBLE(SuperiorNWB), 
     &      SecsInPeriod, SupVerbosity, Month, Q1887, tmplev)
         Q1887= Round(Q1887, CritFlowRound)
      ENDIF 
C
      IF( SupVerbosity .GE. 2) THEN
         WRITE(CTMP,'(A, A3,I5, A, F15.10,A1)')' Pre-project flow for ', 
     &      MonthAbbrev(Month), Year, ': ',  Q1887, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF 
C
      WRITE(SupPreProject, '(I6)') IDNINT(Q1887) 
C     Check Criterion B level
C
C        There is a problem replicating existing practice for 
C        implementing criterions B and C.  The existing implementations
C        use monthly mean lake levels for Superior and Michigan-Huron,
C        based on supplies for the next month.  This is handled by 
C        the parameter P77A Compat.
C
C         =========================================================
C                      EC code
C
C            USSLIP(YEAR,MONTH) = ROUND(USS(SUPFLOW, MHUMEAN),
C    1                                        DBLE(100))
C
CC                If the US Slip level is too high, decrease the Lake
CC                Superior gate setting and compute a new plan flow
C
C             IF (USSLIP(YEAR,MONTH) .GT. DBLE(177.94)) THEN
C                GATESET = GATE(GATESET, -1)
C                CALL QSUP(DBLE(0), DBLE(0), DBLE(0), X, Y, Z,
C    1              DBLE(0), SIDECHAN, GATESET, MONTH, 2)
C                PLANFLOW = ROUND(FLOW(SUPBOP, SUPNTS, SIDECHAN,
C    1              X, Y, Z), DBLE(10))
C                IF (ROUND(PLANFLOW, DBLE(1000)) .EQ.
C    1              ROUND(SIDECHAN, DBLE(1000))) THEN
C                   PLANFLOW = ROUND(SUPFLOW - 1, DBLE(1))
C                END IF
C                QFLOW = PLANFLOW
C             END IF
C
C             IF (ROUND(USSLIP(YEAR,MONTH), DBLE(1000)) .GT.
CC   1           ROUND(DBLE(582.9), DBLE(1000))) GOTO 10
C    1        ROUND(DBLE(177.94), DBLE(1000))) GOTO 10
CC            UNTIL USSLIP <= 177.94
C
CC              Adjust the gate setting
C
C             IF (ROUND(GATESET, DBLE(100)) .EQ. DBLE(0)) THEN
C                GATESET = DBLE(0.5)
C             END IF
C
C         =========================================================
C
C
C
C
C        Get the estimated MH mean for planned St Marys flow, using
C        the operational NBS's etc.
C
      DO 855, I = 1, 999
C
C
C        Route St Marys flow thru middle lakes, to get ave MH level
C
C        XX 40 increments hardcoded here
         CALL P2012MidLakeRoute ( SecsInPeriod,
     &      IDNINT(Round(StMarysFlow, CritFlowRound)) + 
     &      MichuronNWB, StClairNWB, ErieNWB, 
     &      Round(MicHuronBOPLevel,CritBOMRound),
     &      Round(StClairBOPLevel,CritBOMRound),
     &      Round(ErieBOPLevel,CritBOMRound),
     &      SCRivRet, DetRivRet, NiaRivRet,
     &      CritFlowRound, CritLevelRound, CritBOMRound, 40, 
     &      MicHuronEOPLevel, StClairEOPLevel, ErieEOPLevel, 
     &      MicHuronAveLevel)
C 
C           WRITE(CTMP, '(A,A,I2,A4,I5,2(A,F12.5),A,I3,A1)') 
C    &         OutputPrefix,' after P77AMidLakeRoute CritBMethod=', 
C    &         CritBMethod, MonthAbbrev(Month),Year, 
C    &         '  StMarysFlow=', StMarysFlow, 
C    &         '  MicHuronAveLevel=', MicHuronAveLevel, 
C    &         '  CritBOMRound=', CritBOMRound, LLT
C           CALL LogIt(LOGD,0,0)
C
         IF( CritBMethod .EQ. 1) THEN
C           XX Old models round to 10 m3s here.  Go figure.
            tmplev = Round(USS ( IDNINT(Round(StMarysFlow, 1)), 
     &         MicHuronAveLevel), CritBOMRound)
         ELSE
C           XX If using new Crit B method, don't perpetuate the 
C           XX inconsistency.  Find US Slip level using the
C           XX GetLevelAtUSSlip function which, in turn, uses
C           XX the equation specified by the state of USSlipEquation
            tmplev = Round( 
     &         GetLevelAtUSSlip(StMarysFlow, MicHuronAveLevel,
     &            USSlipEquation, month, .TRUE.)
     &         , CritBOMRound)
         ENDIF 
C
C
         IF( SupVerbosity .GE. 5) THEN
            WRITE(CTMP,6547) OutputPrefix, tmplev, 
     &         MicHuronAveLevel, StMarysFlow, Gates, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF 
 6547    FORMAT(A,'US Slip= ',F7.3, ' with MH level= ', F7.3, 
     &      ' and St Marys flow= ', F6.1, ' (', F4.1,')', A1)
C
C
C        If not a problem with US Slip, jump out
         IF( tmplev .LE. 177.94D0 ) THEN
            IF ( tmplev .GE. 177.77D0 ) THEN
                AdjustmentCodes(19:21) = 'BW '
                IF( SupVerbosity .GT. 1) THEN
                  WRITE(CTMP,6550) OutputPrefix, tmplev, LLT
                  CALL LogIt(LOGD,0,0)
 6550             FORMAT(A,'Criterion B Warning! US Slip level (',F7.3, 
     &                  ') is greater than 177.77 m',A1)
                ENDIF
            END IF
            IF( SupVerbosity .GE. 5) THEN
               WRITE(CTMP,6537) OutputPrefix, 'B', LLT
               CALL LogIt(LOGD,0,0)
            ENDIF 
            GOTO 4779
         ENDIF
C
C        Violated Criterion B!
         AdjustmentCodes(19:21) = 'CB '
         IF( SupVerbosity .GT. 1) THEN
            WRITE(CTMP,6548) OutputPrefix, tmplev, LLT
            CALL LogIt(LOGD,0,0)
 6548       FORMAT(A,'Criterion B violation!  US Slip level (',F7.3, 
     &         ') is greater than 177.94. Must close a gate . . .',A1)
         ENDIF
C             
C
         IF( CritBMethod .EQ. 2) THEN
C           Don't allow discharge more than pre-project flow
            IF( StMarysFlow .LE. Q1887 ) THEN
               IF( SupVerbosity .GT. 1) THEN
                  IF( I .EQ. 1 ) THEN 
                     CWRK(1:1) = ' '
                     L = 1
                  ELSE
                     CWRK =  ' '
                     L = 1
                  ENDIF
                  WRITE(*,6529)CWRK(1:L),MonthAbbrev(Month),Year,LLT
                  WRITE(CTMP,6529)CWRK(1:L),MonthAbbrev(Month),Year,LLT
                  CALL LogIt(LOGCDE,0,0)
                  WRITE(CTMP,6530) Q1887,StMarysFlow, Gates, LLT 
                  CALL LogIt(LOGCDE,0,0) 
               ENDIF
               GOTO 9753
            ENDIF
         ENDIF
 6529    FORMAT(' Crit B outflow is', A, 'below ',A,I6,
     &      ' pre-project flow', A1)
 6530    FORMAT(' Pre-project flow: ', F12.1, '  Crit B outflow flow: ', 
     &      F12.1, '  Crit B outflow gates: ', F5.2, A1)
C
         IF( NINT(Gates + 0.51) .EQ. 1 ) THEN
C           Only half a gate open!  Can't close any more!  Logic
C           undefined for existing practice, so just squawk and 
C           move on.  XX Maybe should cut sidechannel flows to 0?
            WRITE(CTMP,6549) MonthAbbrev(Month),Year, LLT
            CALL LogIt(LOGCDE,0,0)
 6549       FORMAT(' Unable to comply with Criterion B during ',A,
     &         I5,'!',A1)
            GOTO 9753
         ENDIF
C
C        Decrement gate setting
         Gates = Gates - 1  
         IF(NINT(Gates) .EQ. 0) Gates = 0.5
C
C        XX This may be possible in weird scenario with MH very high,
C        XX so should come up with logic someday to handle it.
         IF(Gates .LT. 0.0) STOP ' Ginger'
C
C        Get new mean outflow for adjusted gate setting, loop again.
C        Put it straight into StMarysFlow, if were at 1/2 gate 
C        shouldn't have gotten here anyway.  Be aware that this 
C        leaves a hole in the logic regarding minimum flow situations.
         CALL P77ASolveSup(SuperiorBOPLevel, DBLE(SuperiorNWB), 
     &      SecsInPeriod, Gates, CritFlowRound, CritBOMRound, 
     &      SuperiorEOPLevel, StMarysFlow, SuperiorAveLevel, Month)
C
C        Knock accuracy down to cm or whatever
  855    SuperiorEOPLevel = Round( SuperiorEOPLevel, CritBOMRound )
C
C     Get here only by falling out of loop.  How?
      STOP ' Mary Ann' 
C
C
 9753 CONTINUE 
C
C
C     Sidechannel (i.e., non-gated) flow determination
C 
C     Compute non-gated flow for the month.  Need to determine GateFlow,
C     and subtract it from the nominal total outflow.  This requires a
C     SW Pier level, computed from the current estimate of the Superior 
C     monthly mean level.  Existing practice was to perform these  
C     calculations by hand, so all levels are rounded to cm.  
C
C     This is set up for iterative solution of simultaneous equations in
C     if we want to use the new EC equation for LS:SWP.  The equations are 
C     simultaneous because SWPier level is a function of both Superior
C     level and outflow, while outflow is a function of Superior level.
C     So solve for SWP, get NonGatedFlow, get SWP, get NonGatedFlow, etc
C     until converge.  
C
C     Traditional practice uses the old SWP:LS relationship, which is not
C     a function of outflow, so there is no need to solve simultaneous 
C     equations.  When the old SWP:LS relationship is selected, the 
C     loop cycles twice (need at least two iterations to test for 
C     convergence), computing the same answer both times.
C              
C     Seed NonGatedFlow with some number
      NonGatedFlow = 2320
C
C
C
      IF(SupVerbosity .GE. 6 ) THEN
         WRITE(CTMP,'(A, A1)') ' Solving for sidechannel flow . . .',LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
      L = -9999
      Flow2ItsAgo = L
C
      DO 4778, I=1,999
         CALL GetSWPier( SuperiorAveLevel, Month, Gates, 
     &      IDNINT(NonGatedFlow), CritBOMRound, LStoSWP, 6, Year,
     &      'P2012_Criterions'//LLT, SWPier, gatedflow)
C
C        This is the method that Matt proposes, not the way it
C        is calculated in "existing practices". 
         NonGatedFlow = IDNINT(StMarysFlow - gatedflow)
C        XX Backward compat fix:
C        XX This following line replicates the goofy mixed rounding of
C        XX of existing practices.  If we're rounding flows, we'll do
C        XX it the old fashioned way. Remove it when allowed.
         IF(CritFlowRound .NE. -99 ) NonGatedFlow = 
     &      Round(StMarysFlow,1) - Round(gatedflow,CritFlowRound)
C        XX End of fix block

C         WRITE(CTMP,'(A,I4,A,I15,A,F12.5,A1)') 
C     &      '         . . . CritFlowRound=', CritFlowRound,
C     &      '  L=', L, '  new NonGatedFlow=', NonGatedFlow, LLT
C         CALL LogIt(LOGD,0,0)

C
         IF( L .EQ. IDINT(gatedflow * 100000.D0) ) GOTO 4779
C        Check to see if oscillating between two values 
         IF(CritBOMRound .NE. -99 )THEN
            IF(IDINT(gatedflow*100000.D0) .EQ. Flow2ItsAgo) GOTO 4779
            Flow2ItsAgo = L
         ENDIF

C         WRITE( CTMP, '(A,I15,A1)')  
C     &       '         . . . IDINT(gatedflow * 100000.D0)=', 
C     &       IDINT(gatedflow * 100000.D0), LLT
C         CALL LogIt(LOGD,0,1)

 4778    L = IDINT(gatedflow * 100000.D0)
C
 4779 IF(SupVerbosity .GE. 3 ) THEN
         WRITE( CTMP, 2310) NonGatedFlow,LLT
         CALL LogIt(LOGD,0,0)
         WRITE( CTMP, 2309) StMarysFlow, SuperiorAveLevel, 
     &      SWPier, gatedflow, LLT
         CALL LogIt(LOGD,0,1)
      ENDIF
C
      WRITE(SupSideChanFlow,'(I6)') 
     &   IDNINT(Round(NonGatedFlow,CritFlowRound))
      IF(SupVerbosity .GE. 7 ) THEN
         WRITE( CTMP, '(A,A1)') ' finished SupCriterions',LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
 2309 FORMAT( '   (based on nominal total flow of ', F9.4, 
     &   ' Sup/SW Pier levels of ', 2F9.4,', and gate flow of', F9.4, 
     &   ')', A1)
 2310 FORMAT( '   Computed "Official" non-gated flow as ', F9.4, A1)
C
 6537 FORMAT(A, ' . . . passed Criterion ', A1, ' check OK.', A1)

      END      
C
C
C
C    ************************************************************************
C    Plan 2012
C    ************************************************************************
      SUBROUTINE P2012_Regulation( MonthIndex, SupVerbosity, SPBOP,
     &   MHBOP, P77AWinterMax, FlowRound, NominalStMarys, SupGatesopen,
     &   ComputedNonGate, PreviousGates, IndexGatesSide)
      IMPLICIT NONE

      DOUBLE PRECISION SPBOP, MHBOP, ComputedNonGate, NominalStMarys

      INTEGER  MonthIndex, SupVerbosity, FlowRound, P77AWinterMax
C
      DOUBLE PRECISION Tolerance 
      DOUBLE PRECISION SWPIERLVL
      DOUBLE PRECISION Spoutk, SPout1, MHSLOPE, NovMax, SPSLOPE
      DOUBLE PRECISION SPBOk, 
     &                 SPIMLVK, EOILevels2Ago, gatedflow
      DOUBLE PRECISION P2012SPTARGET(12), P2012MHTARGET(12),
     &                 P2012SPHIGHLEVEL, P2012SPLOWLEVEL,
     &                 P2012MHTARGETADJ
      DOUBLE PRECISION P2012PARA(6),AlterWinterMax   
      DOUBLE PRECISION ForcedSupGates(12)
      INTEGER          ForcedSideChan(12)
      INTEGER       MaxIterations, IndexGatesSide
      CHARACTER     CTMP*255, CWRK*80
      CHARACTER     LLT*1
      LOGICAL       LOGC(8), LOGD(8), LOGS(8), 
     &              LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &              LOGM(8), LOGP(8), LOGO(8), 
     &              LOGDM(8), LOGDP(8), LOGDO(8), 
     &              LOGSM(8), LOGSP(8), LOGSO(8),
     &              LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &              LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      DOUBLE PRECISION JUNEFLOW
C
      REAL MinGates(12)
      INTEGER       SideChanMaxCapacity(12), SideChanMinCapacity(12)
      DOUBLE PRECISION CWOvertopElev, CWMinimumFreeboard(12)     
      INTEGER       MaxSGSetting, COUNTER_JUNEFLOW, Last_JuneFlow
C                   The number of elements in SGX, SGY, and SGZ
C
      REAL          Supgatesopen, PreviousGates
C                   The number of gates open in current routing period.
C
      CHARACTER     SupPreProject*6, BalancingFlow*8, SupGateFlow*6,
     &              SupPlanFlow*6, SupBOMLevel*8, SWPierLevel*8,
     &              SupSideChanFlow*6, SPRegCodes*40, P77AHeader*120,
     &              P2012Header
C
      DOUBLE PRECISION Round
C
C                   Lake Superior regulation output strings
      COMMON /SUPOUT/   SupPreProject, BalancingFlow, SupGateFlow,
     &                  SupPlanFlow, SupBOMLevel, SWPierLevel, 
     &                  SupSideChanFlow, SPRegCodes, P77AHeader,
     &                  P2012Header
      COMMON /SPOTOP/   SideChanMaxCapacity, CWOvertopElev, 
     &                  CWMinimumFreeboard, SideChanMinCapacity      
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /P2012/    P2012PARA, P2012SPTARGET, P2012MHTARGET, 
     &                  P2012SPHIGHLEVEL, P2012SPLOWLEVEL,
     &                  P2012MHTARGETADJ 
      COMMON /FORGATE/ ForcedSupGates     
      COMMON /FORSIDE/ ForcedSideChan 

      COMMON /STURGEON/ COUNTER_JUNEFLOW, Last_JuneFlow 
      DATA MaxIterations /9999/
      DATA MinGates / 12*0.5 /
      DATA NovMax / 3260 / 
      DATA AlterWinterMax / 2690 / 
      DATA JUNEFLOW /1700.0D0/               
     
      DATA MaxSGSetting /17/
 
C     Log the entering data for debugging porpoises
      IF( SupVerbosity .GE. 6) THEN
         WRITE(CTMP,'(A,I7, 6(A,F12.6),A, I6, A1)'),
     &      'p2012: IndexGatesSide= ', IndexGatesSide,' SPBOP=',SPBOP,
     &      ' m   MHBOP=', MHBOP, 'SupGatesopen=', SupGatesopen,
     &      'ComputedNonGate=', ComputedNonGate, 'PreviousGates= ',
     &      PreviousGates, '   ForcedSupGates=', 
     &      ForcedSupGates(MonthIndex),'   ForcedSupSide=',
     &      ForcedSideChan(MonthIndex), LLT
         CALL LogIt(LOGD,0,0)
      ENDIF

      Tolerance = 10.0D0**FlowRound
      IF( Tolerance .LT. 0.1D0 ) Tolerance = 0.1D0
      SPBOk = SPBOP
C     Initialize the outflows with bogus number
      SPOutk = -9999.D0
C     Initialize the mean level value with the initial level
      SPIMLVk = SPBOk
C
      EOILevels2Ago = -999            
      
C     using beginning month level to calculate outflow

      SPout1 = 824.7*(SPbop-181.43 + P2012PARA(2))**1.5D0 
      SPSLOPE = (1+(SPbop - P2012SPTARGET(MonthIndex))*P2012PARA(1))
      IF (SPbop .LE. P2012SPLOWLEVEL) THEN
          MHSLOPE =(1-P2012PARA(3)*(MHbop-P2012MHTARGET(MonthIndex)
     &            + P2012MHTARGETADJ+0.20*(P2012SPLOWLEVEL-SPbop)))
     &            **P2012PARA(4)
      ELSE
          MHSLOPE =(1-P2012PARA(3)*(MHbop-P2012MHTARGET(MonthIndex)
     &            ))**P2012PARA(4)
      END IF           
      SPoutk = SPout1 * MHSLOPE * SPSLOPE
      IF(MonthIndex .EQ. 11 ) THEN
         IF (Spoutk .GE. NovMax) THEN
            Spoutk = NovMax
         END IF           
      END IF   
      IF(MonthIndex .GE. 5 .AND. MonthIndex .LE. 11) THEN         
         IF(SPbop .GE. P2012PARA(5) ) THEN
            SPoutk = P2012PARA(6)
         END IF   
      END IF 
      IF(MonthIndex .GE. 12 .or. MonthIndex .LE. 4) THEN
         IF (SPbop .GE. P2012SPHIGHLEVEL) THEN
            IF (Spoutk .GE. AlterWinterMax) THEN
                Spoutk = AlterWinterMax
            END IF
         ELSE
            IF (Spoutk .GE. P77AWinterMax) THEN
                Spoutk = P77AWinterMax
            END IF
         END IF      
      END IF 
 
c      SPRegCodes(1:3) = 'PP'
      WRITE(CTMP,'(6(A,F12.6),A1)'), 
     &      ' SPBOP=',SPBOP,' m   Suptarg=',P2012SPTARGET(MonthIndex), 
     &      ' m   Spout1', spout1,
     &      '  MHbop=', MHbop,'  MHslope', mhslope, 
     &      ' Spoutk=', SPoutk, LLT
      CALL LogIt(LOGD,0,0)
  
C     Modified on Sept 24, 2014
      
     
      IF (MonthIndex .EQ. 6 .AND. Round(SPoutK/10,0)*10 .LT. JUNEFLOW)
     &   THEN   
            COUNTER_JUNEFLOW = COUNTER_JUNEFLOW + 1
      END IF
      IF (MonthIndex .EQ. 6 .AND. Round(SPoutK/10,0)*10 .GE. JUNEFLOW)
     &  THEN   
            COUNTER_JUNEFLOW = 0
      END IF 
  
        
      IF (MonthIndex .EQ. 6 .AND. COUNTER_JUNEFLOW .GE. 5) THEN
          Call GetInstantGate(SPBOP, JUNEFLOW, MonthIndex,1, 
     &         PreviousGates, SupGatesopen, GatedFlow, SWPIERLVL,
     &         ComputedNonGate, NominalStMarys)        
          IF (NominalStMarys .GE. JUNEFLOW) THEN
               COUNTER_JUNEFLOW = 0
          END IF    
      ELSE 
          IF(IndexGatesSide .EQ. 0) THEN
              Call GetInstantGate(SPBOP, Round(SPoutK/10,0)*10, 
     &        MonthIndex,0, PreviousGates, SupGatesopen, GatedFlow, 
     &        SWPIERLVL, ComputedNonGate, NominalStMarys) 
          ELSE IF (IndexGatesSide .EQ. 1) THEN
C         Forced Gates          
              IF(SNGL(ForcedSupGates(MonthIndex)) .EQ. -9999) THEN
                 Call GetInstantGate(SPBOP, Round(SPoutK/10,0)*10, 
     &           MonthIndex,0, PreviousGates, SupGatesopen, GatedFlow,  
     &           SWPIERLVL, ComputedNonGate, NominalStMarys) 
              ELSE             
                 Call GetForcedGates(SPBOP, Round(SPoutK/10,0)*10, 
     &         MonthIndex, SNGL(ForcedSupGates(MonthIndex)), GatedFlow, 
     &         SWPIERLVL, ComputedNonGate, NominalStMarys)
               SupGatesOpen = SNGL(ForcedSupGates(MonthIndex))
              END IF  
          ELSE IF (IndexGatesSide .EQ. 2) THEN
C         Forced SideChan  
              IF ( DBLE(ForcedSideChan(MonthIndex))  .EQ. -9999) THEN 
                 Call GetInstantGate(SPBOP, Round(SPoutK/10,0)*10, 
     &           MonthIndex,0, PreviousGates, SupGatesopen, GatedFlow,  
     &           SWPIERLVL, ComputedNonGate, NominalStMarys) 
              ELSE                    
                 Call GetForcedSideChan(SPBOP, Round(SPoutK/10,0)*10,
     &                MonthIndex, 0, PreviousGates, SupGatesopen, 
     &                GatedFlow, SWPIERLVL,
     &                DBLE(ForcedSideChan(MonthIndex)), NominalStMarys)
                  ComputedNonGate = DBLE(ForcedSideChan(MonthIndex))
              END IF    
          ELSE IF (IndexGatesSide .EQ. 3) THEN
C         Forced Gates and SideChan       
             IF(SNGL(ForcedSupGates(MonthIndex)) .EQ. -9999 .AND.
     &          DBLE(ForcedSideChan(MonthIndex)) .EQ. -9999  ) THEN   
                Call GetInstantGate(SPBOP, Round(SPoutK/10,0)*10, 
     &            MonthIndex,0, PreviousGates, SupGatesopen, GatedFlow,  
     &            SWPIERLVL, ComputedNonGate, NominalStMarys) 
             ELSE IF (SNGL(ForcedSupGates(MonthIndex)) .EQ. -9999 )THEN
                Call GetForcedSideChan(SPBOP, Round(SPoutK/10,0)*10,
     &                MonthIndex, 0, PreviousGates, SupGatesopen, 
     &                GatedFlow, SWPIERLVL,
     &                DBLE(ForcedSideChan(MonthIndex)), NominalStMarys)
                ComputedNonGate = DBLE(ForcedSideChan(MonthIndex))
             ELSE IF ( DBLE(ForcedSideChan(MonthIndex)) .EQ. -9999)THEN
                  Call GetForcedGates(SPBOP, Round(SPoutK/10,0)*10, 
     &         MonthIndex, SNGL(ForcedSupGates(MonthIndex)), GatedFlow, 
     &            SWPIERLVL, ComputedNonGate, NominalStMarys)
               SupGatesOpen = SNGL(ForcedSupGates(MonthIndex)) 
             ELSE                                 
              Call GetForcedGatesSideChan(SPBOP, Round(SPoutK/10,0)*10, 
     &        MonthIndex, SNGL(ForcedSupGates(MonthIndex)), GatedFlow, 
     &        SWPIERLVL, DBLE(ForcedSideChan(MonthIndex)),
     &        NominalStMarys)
              SupGatesOpen = SNGL(ForcedSupGates(MonthIndex))
              ComputedNonGate = DBLE(ForcedSideChan(MonthIndex))
             END IF   
          END IF
  
      END IF      

      
C     Store some intermediate results in character variables
      WRITE(SupPlanFlow,'(I6)') IDNINT(NominalStMarys) 
      WRITE(SupPreProject,'(F6.3)') (MHSLOPE)           
      WRITE(BalancingFlow,'(I8)') (COUNTER_JUNEFLOW)      

      IF( SupVerbosity .GE. 2) THEN
         WRITE(CTMP,'(A,I5, A, I5, 4(A,F12.6),A,F5.2,A1)') 
     &      '  leaving P2012, month=', (MonthIndex),
     &      '  COUNTER_JUNEFLOW=',COUNTER_JUNEFLOW, 
     &         'SPBOP=', SPBOP, ' SIMLVk=', SPIMLVk,  
     &      ' NominalStMarys=', NominalStMarys, '  ComputedNonGate=',
     &      ComputedNonGate, '  SupGates=', SupGatesOpen, LLT
         CALL LogIt(LOGD,0,0)
      END IF

      END