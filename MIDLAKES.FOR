C----------------------------------------------------------------------
C       
C
C      This is the header block from the original MIDLAKES.FOR from GLERL,
C      retained here for reference.  Only the last two sentences remain
C      applicable; although this model's code was the basis for the CGLRRM,
C      after instituting user-selected routing intervals, infinite looping
C      through non-fixed number of periods, purely double-precision
C      arithmetic, extended error messaging, different output formats,
C      pre-processing of time-series input data, variable-surface area,
C      and other design objectives, almost every line was either re-written
C      or discarded.  The functions DATESEQUENCE and OUTFLOW remain very
C      close to their original implementation.  Although now combined into
C      a single function and apparently quite extensively modified, the 
C      changes to the COEFF and MATRIX subroutines are largely superficial,
C      preserving the heart of the computational mechanism as close as 
C      possible to the version in GLERL's tech memo.
C      ----------------------------------------------------------------      
C                NOAA/Great Lakes Environmental Research Laboratory
C                            Ann Arbor, Michigan
C      -----------------------------------------------------------------
C      -----------------------------------------------------------------
C
C     Program Name      Midlakes
C     Language          FORTRAN
C     Compiler          Lahey Revision B
C
C     Source File       midlakes.for
C     Object File       midlakes.o(bj)
C
C     Subroutines       COEFF, MATRIX, STATS, STATSUB, OUTDATA,
C                       MONOUTLINE, STATOUTLINE, WRITEMONLINE,
C                       MONSTATSLINE
C                       
C      
C       This program is a hydrologic response model of the unregulated
C       portion of the Great Lakes (Michigan-Huron, St. Clair, and Erie,
C       and their connecting channels).  Quarter-monthly supplies are 
C       input either in the form of Net Basin Supplies or as hydrologic
C       components.  Input includes overlake precipitation, 
C       basin runoff, lake evaporation, diversions, groundwater, 
C       consumptive use, initial outflows and lake levels
C       for each lake and ice retardation rates and discharge equation
C       parameters for each of the connecting channels.  Outputs include
C       monthly mean lake levels, beginning-of-month levels, and
C       monthly flow rates in the connecting channels.  Quarter-monthly 
C       Lake Ontario inflows are also calculated.  The Midlakes routing
C       model is comprised of continuity equations for each lake, solved
C       using a second-order finite-difference technique (Quinn, F.H.,
C       1978. Hydrologic Response Model of the North American Great 
C       Lakes, Journal of Hydrology 37(3) pp. 295-307.).
C  
C       For further information contact:
C             Anne H. Clites
C             Great Lakes Environmental Research Laboratory
C             2205 Commonwealth Blvd.
C             Ann Arbor, Michigan  48105-1593
C             clites@glerl.noaa.gov
C    
C       Last revision: September 1997
C
      SUBROUTINE COEFFMATRIX(DT, IncNum)
C
C     This subroutine calculates the 9 matrix coefficients for a specified
C     combination of gauges. When b=0, there is a single gauge; when b is 
C     greater than 0, there is more than one gauge involved.  This routine
C     written by Anne H. Clites, NOAA/GLERL, October, 1995.
C
C     Ruthlessly butchered by Matt during the spring of 1998.
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      DOUBLE PRECISION DT
C                   (Input) The number of seconds per increment.  Intrinsic
C                   part of the calculations, DT stands for "delta t". 
C
      INTEGER       IncNum
C                   (Input) The number of increments within each interval.
C                   Needed here only for logging purposes.
C
C
C     Declarations regarding local variables or functions:
C
      DOUBLE PRECISION CDET,Z1DET,Z2DET,Z3DET
      DOUBLE PRECISION C11,C12,C13,C21,C22,C23,C31,C32,C33
      DOUBLE PRECISION MHWB, SCWB, ERWB
      DOUBLE PRECISION OUT1,OUT2,OUT3,OUT4,OUT5,OUT6
      DOUBLE PRECISION OUTMH, OUTSC, OUTER
C                   Coefficients, determinants, etc, used in the
C                   matrix solution of the simultaneous equations.
C                   Aliased to values in FWRK for economy of storage
C                   and passing values to/from calling routine.
C
      DOUBLE PRECISION MHBOI, SCBOI, ERBOI
      DOUBLE PRECISION MHNTS, SCNBS, ERNBS
      DOUBLE PRECISION MHICE, SCICE, ERICE
C                   Input values assigned by calling routine.
C                   Aliased to values in FWRK for economy of storage
C                   and passing values to/from calling routine.
C                   MHBOI, etc are beginning of increment levels.
C                   MHICE, etc are retardation values for the period.
C                   MHNTS, SCNBS, ERNBS are summarized inflow/outflow
C                   for each lake.  MH includes St. Marys River flow, 
C                   hence named NTS instead of NBS.
C
      DOUBLE PRECISION MHEOI, SCEOI, EREOI
C                   Output end-of-increment values calculated by this 
C                   routine.  Aliased to values in FWRK for economy of
C                   storage and passing values to/from calling routine.
C
      DOUBLE PRECISION NIZ, DeltaZ4, CGMLV
C                   Variables used for two-gage equation on Niagara 
C                   River.  NIZ is mean level for the previous 
C                   routing interval.  DeltaZ4 is difference between
C                   NIZ and mean level for current routing interval.
C                   CGMLV is mean level for current routing interval.
C                   Note that these storages are actually aliased to
C                   memory in FWRK.
C                   
      DOUBLE PRECISION GETAREA
C                   Function that returns surface area for a lake
C
      DOUBLE PRECISION OUTFLOW
C                   Function that returns connecting channel flow, given 
C                   the parameters for the discharge equation.
C
C 
C
C     Declarations regarding COMMON block variables:
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                   See description in Block Data Genral
C
      CHARACTER     LLT*1
C                   This is the character that marks the end of a line 
C                   sent to the LogIt subroutine.
C 
      CHARACTER     CTMP*255, CWRK*80
C                   These are just buffers used throughout program.
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
      INTEGER       TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &              MissingVal, SecsInUniformMon
      DOUBLE PRECISION MPFT, HALF, ONE, ZERO
C                   See descriptions in Block Data General.  
C
      DOUBLE PRECISION MHK,MHYM,MHA,MHB,MHWT,MHC
      DOUBLE PRECISION SCK,SCYM,SCA,SCB,SCWT,SCC
      DOUBLE PRECISION ERK,ERYM,ERA,ERB,ERWT,ERC
C                   See descriptions in CHECKPAR and OUTFLOW
C
      DOUBLE PRECISION FWRK
C                   This is a common block used for 64-bit calculations.
      DIMENSION FWRK(225,12)
C                   Dimensioned it this way to provide a 225 element
C                   array for each pool lake, and some scratch.
C       
C
C     Common block statements.

      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /STAGE/    MHK,SCK,ERK,MHYM,SCYM,ERYM,MHA,SCA,ERA,
     &                  MHB,SCB,ERB,MHWT,SCWT,ERWT,MHC,SCC,ERC
      COMMON /FSPACE/   FWRK
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
C        
C         
C     Assign convenient aliases for summary inflow values, which
C     where copied into FRWK(1,?) in calling subroutine.  Remember that
C     the M-H summary inflow already includes St. Marys River 
C     contribution.  Also recognize that consumptive use and diversions,
C     etc, are already combined into each of these terms, which is a bit
C     of a change from GLERL's original approach.
      EQUIVALENCE (MHNTS, FWRK(1,1))
      EQUIVALENCE (SCNBS, FWRK(1,2))
      EQUIVALENCE (ERNBS, FWRK(1,3))
C       
C     Assign convenient aliases for summary ice/weed retardation values, 
C     which where copied into FRWK(2,?) in calling subroutine.
      EQUIVALENCE (FWRK(2,1), MHICE)
      EQUIVALENCE (FWRK(2,2), SCICE)
      EQUIVALENCE (FWRK(2,3), ERICE)
C        
C     Assign convenient aliases for beginning of increment levels,
C     which where copied into FRWK(3,?) in calling subroutine.
      EQUIVALENCE (FWRK(3,1), MHBOI)
      EQUIVALENCE (FWRK(3,2), SCBOI)
      EQUIVALENCE (FWRK(3,3), ERBOI)
      EQUIVALENCE (FWRK(3,4), NIZ)
C        
C     Assign convenient aliases for end of increment levels, which will 
C     be written to FRWK(4,?) in at the end of this subroutine.
      EQUIVALENCE (FWRK(4,1), MHEOI)
      EQUIVALENCE (FWRK(4,2), SCEOI)
      EQUIVALENCE (FWRK(4,3), EREOI)
      EQUIVALENCE (FWRK(4,4), DeltaZ4)
C
C
C     FWRK(5,?) reserved for period average levels in calling sub, 
C     except for CGMLV, which is actually input variable that gets
C     read instead of computed:
      EQUIVALENCE (CGMLV, FWRK(5,4))
C
C     FWRK(6,?) reserved for period average outflows in calling sub.
C     FWRK(7,?) reserved for beginning of period levels in calling sub.
C     FWRK(8,?) reserved for end of period levels in calling sub.
C     FWRK(9,?) reserved for solution check values in calling sub.
C     DO NOT equivalence this memory in this function!
C
C     Assign convenient aliases for water balance terms, using storage
C     actually allocated to the common workspace.  This is for 
C     legibility and remaining close to GLERL's original code, while
C     still conserving memory.
      EQUIVALENCE (FWRK(14,1), MHWB)
      EQUIVALENCE (FWRK(14,2), SCWB)
      EQUIVALENCE (FWRK(14,3), ERWB)
C
C     Assign convenient aliases for matrix determinants, using storage
C     actually allocated to the common workspace.  This is for legibility
C     and remaining close to GLERL's original code, while still conserving
C     memory.
      EQUIVALENCE (FWRK(10,1), Z1DET)
      EQUIVALENCE (FWRK(10,2), Z2DET)
      EQUIVALENCE (FWRK(10,3), Z3DET)
      EQUIVALENCE (FWRK(10,4), CDET)
C
C     Assign convenient aliases for matrix coeffients, using storage
C     actually allocated to the common workspace.  This is for legibility
C     and remaining close to GLERL's original code, while still conserving
C     memory.
      EQUIVALENCE (FWRK(11,1), C11)
      EQUIVALENCE (FWRK(12,1), C21)
      EQUIVALENCE (FWRK(13,1), C31)
      EQUIVALENCE (FWRK(11,2), C12)
      EQUIVALENCE (FWRK(12,2), C22)
      EQUIVALENCE (FWRK(13,2), C32)
      EQUIVALENCE (FWRK(11,3), C13)
      EQUIVALENCE (FWRK(12,3), C23)
      EQUIVALENCE (FWRK(13,3), C33)
C        
C     Assign convenient aliases for the intermediate outflow terms, using 
C     storage actually allocated to the common workspace.  This is for 
C     legibility and remaining close to GLERL's original code, while still
C     conserving memory.
      EQUIVALENCE (FWRK(32,1), OUTMH)
      EQUIVALENCE (OUT1, FWRK(33,1))
      EQUIVALENCE (FWRK(34,1), OUT2)
      EQUIVALENCE (FWRK(32,2), OUTSC)
      EQUIVALENCE (FWRK(33,2), OUT3)
      EQUIVALENCE (FWRK(34,2), OUT4)
      EQUIVALENCE (FWRK(32,3), OUTER)
      EQUIVALENCE (FWRK(33,3), OUT5)
      EQUIVALENCE (FWRK(34,3), OUT6)
C
C
C 
C  Calls to Function OUTFLOW needed in this subroutine.
C
      OUTMH=OUTFLOW(MHBOI,SCBOI,MHK,MHYM,      MHA,      MHB,MHWT,MHC)
      OUT1= OUTFLOW(MHBOI,SCBOI,MHK,MHYM,MHA - ONE,      MHB,MHWT,MHC)
      OUT2= OUTFLOW(MHBOI,SCBOI,MHK,MHYM,     MHA ,MHB - ONE,MHWT,MHC)
      OUTSC=OUTFLOW(SCBOI,ERBOI,SCK,SCYM,     SCA ,      SCB,SCWT,SCC)
      OUT3= OUTFLOW(SCBOI,ERBOI,SCK,SCYM,SCA - ONE,      SCB,SCWT,SCC)
      OUT4= OUTFLOW(SCBOI,ERBOI,SCK,SCYM,     SCA ,SCB - ONE,SCWT,SCC)
C     For single gage eqn, ERWT = 1, and ERB = 0, so doesn't matter a lick
C     what you use for NIZ.  For double gage eqn, NIZ should be BOI level,
C     which is assumed to be the CGIP level from the previous interval if
C     going through the first increment, or else the level of the current
C     interval if a subsequent increment.
      OUTER=OUTFLOW(ERBOI,NIZ,ERK,ERYM,     ERA ,      ERB,ERWT,ERC) 
      OUT5= OUTFLOW(ERBOI,NIZ,ERK,ERYM,ERA - ONE,      ERB,ERWT,ERC)
      OUT6= OUTFLOW(ERBOI,NIZ,ERK,ERYM,ERA      ,ERB - ONE,ERWT,ERC)
C
C  Calculate matrix coefficients for this time step.
C
C     Replaced variable containing constant surface area with call to
C     function GETAREA().  If ever using variable area, this presumes that 
C     areas calculated for MHZ, SCZ, and ERZ (the BOI lake levels) closely
C     approximate the area corresponding from the mean level.
      C11= HALF*MHA*MHWT*OUT1 + HALF*MHB*OUT2 +GETAREA(2,MHBOI)/DT
      C12 = HALF*MHA*(ONE - MHWT)*OUT1 - HALF*MHB*OUT2  
      C13 = ZERO
      C21 = -(HALF*MHA*MHWT*OUT1+HALF*MHB*OUT2)
      C22 = HALF*SCA*SCWT*OUT3 + HALF*SCB*OUT4 - HALF*MHA*(ONE - 
     &      MHWT) * OUT1 - HALF*MHB*OUT2 + GETAREA(5, SCBOI)/DT 
      C23 = -(-HALF*SCA*(ONE - SCWT)*OUT3 + HALF*SCB*OUT4)
      C31 = ZERO
      C32 = -(HALF*SCA*SCWT*OUT3 + HALF*SCB*OUT4) 
      C33 = HALF*ERA*ERWT*OUT5 + HALF*ERB*OUT6 - HALF*SCA*(ONE-SCWT) *
     &      OUT3 + HALF*SCB*OUT4 + GETAREA(6,ERBOI)/DT 

C  Calculate deltaz4, the change in levels of the Chippewa-
C  Grass Island Pool, which is needed to calculate the Lake Erie water balance 
C  term, ermwb, below when stage discharge relationship is based on more than
C  one gauge (when erb > 0).
C    
C      IF (ERB .GT. 0.) THEN
C         IF (I.EQ.1) THEN
C            DELTAZ4(QM,MON) = NIAQLEV(QM,MON)-NIZ
C         ELSE
C            DELTAZ4(QM,MON) = 0.
C         ENDIF
C      ENDIF
C
C  In case using two-gage eqn, initialize DeltaZ4 as change in level 
C  from previous routing increment.  If first increment, NIZ will have
C  value of CGMLV from previous routing interval.  If not first increment
C  then NIZ will have value of CGMLV from current routing interval, and 
C  DeltaZ4 comes out as zero;  if first (if not using two-gage, will end
C  up as zero, since CGMLV and NIZ are 0.0D):
       DeltaZ4 = CGMLV - NIZ
C
C  WB's are "water balance" summations...NBS, diversions, consumptive 
C  use, inflow, outflow, and ice retardation. These terms are needed
C  to find the matrix solution in the next subroutine.
C
C      MHWB = MHQNBS(QM,MON)+MIQDIV(QM,MON)+HUQDIV(QM,MON)+MIQGDW(QM,MON
C     &  )+HUQGDW(QM,MON)-MIQCON(QM,MON)-HUQCON(QM,MON)+SPQFLO(QM,MON)-
C     &  OUTMH+MHQICE(QM,MON)
C         
       MHWB = MHNTS - OUTMH + MHICE
C             
C      SCWB = SCQNBS(QM,MON)+SCQDIV(QM,MON)+SCQGDW(QM,MON)-
C     &   SCQCON(QM,MON)+OUTMH-MHQICE(QM,MON)-OUTSC+SCQICE(QM,MON)
C         
       SCWB = SCNBS + OUTMH - MHICE - OUTSC + SCICE

C      ERWB = ERQNBS(QM,MON)+ERQDIV(QM,MON)+ERQGDW(QM,MON)-
C     &   ERQCON(QM,MON)+OUTSC-SCQICE(QM,MON)-OUTER+ERQICE(QM,MON)+
C     &   (.5*ERB*OUT6-.5*ERA*(1.-ERWT)*OUT5)*DELTAZ4(QM,MON)
C         
       ERWB = ERNBS + OUTSC - SCICE - OUTER + ERICE +
     &   (.5*ERB*OUT6-.5*ERA*(1.-ERWT)*OUT5)*DeltaZ4
C
C  The follwing code used to be in it's own subroutine called 'MATRIX'.
C  It was combined with the previous code (which used to be called 'COEFF').
C  
C  This subroutine uses the coefficients calculated above to solve for the        
C  unknowns deltaz1, deltaz2, and deltaz3 by way of a matrix.  These deltaz's 
C  define the new end-of-increment levels.  This subroutine
C  was written by Anne H. Clites, NOAA/GLERL, November, 1995.
C
C      MHZ2=0.
C      SCZ2=0.
C      ERZ2=0.

C  cdet is the determinant of the C matrix.
C  z1det, z2det, AND z3det are the determinants of the Michigan-Huron,
C  St. Clair, and Erie matrices, respectively.

      CDET = C11*((C22*C33)-(C23*C32)) - C12*((C21*C33)-(C23*C31))
     &       + C13*((C21*C32)-(C22*C31))

      Z1DET = MHWB*((C22*C33)-(C23*C32)) - C12*((SCWB*C33)-
     &        (C23*ERWB)) + C13*((SCWB*C32)-(C22*ERWB))

      Z2DET = C11*((SCWB*C33)-(C23*ERWB)) - MHWB*((C21*C33)-
     &        (C23*C31)) + C13*((C21*ERWB)-(SCWB*C31))

      Z3DET = C11*((C22*ERWB)-(SCWB*C32)) - C12*((C21*ERWB)-
     &        (SCWB*C31)) + MHWB*((C21*C32)-(C22*C31))

C  The change in level equals the determinant divided by the cdet.

C      DELTAZ(1) = Z1DET/CDET     
C      DELTAZ(2) = Z2DET/CDET   
C      DELTAZ(3) = Z3DET/CDET

C  New end-of-increment levels.

      MHEOI = MHBOI + Z1DET/CDET
      SCEOI = SCBOI + Z2DET/CDET
      EREOI = ERBOI + Z3DET/CDET
C
  620 FORMAT(I3,3X,4A18,A1)
  621 FORMAT(I3,3X,6A15,A1)
C  622 FORMAT(I3,A3,6D15.7,A1)
  623 FORMAT(I3,D21.11,3D18.11,A1)
  624 FORMAT(I3,D18.7,5D15.7,A1)
  625 FORMAT(I3,3X,3D15.7,3D15.7,A1)
  626 FORMAT(I3,A3,6F15.5,A1)
C
C
      IF( MLVerbosity .GE. 9) THEN
         WRITE(CTMP,620)INCNUM,'M-H BOI','STC BOI','ERI BOI',
     &      'CGI BOI',LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,623)INCNUM,MHBOI,SCBOI,ERBOI,NIZ,LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,621)INCNUM,'NBS','ICE','WB','OUT0','OUT1',
     &      'OUT2',LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,626)INCNUM,'M-H',MHNTS,MHICE,MHWB,OUTMH,OUT1,
     &      OUT2,LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,626)INCNUM,'STC',SCNBS,SCICE,SCWB,OUTSC,OUT3,
     &      OUT4,LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,626)INCNUM,'ERI',ERNBS,ERICE,ERWB,OUTER,OUT5,
     &      OUT6,LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,621)INCNUM,'C11','C21','C31','C12','C22','C32',LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,624)INCNUM,C11,C21,C31,C12,C22,C32,LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,621)INCNUM,'C13','C23','C33','Z1DET','Z2DET',
     &      'Z3DET',LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,625)INCNUM,C13,C23,C33,Z1DET,Z2DET,Z3DET,LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,620)INCNUM,'CDET','M-H EOI','STC EOI','ERI EOI',LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,623)INCNUM,CDET,MHEOI,SCEOI,EREOI,LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,620)INCNUM,'    ','M-H AVE','STC AVE','ERI AVE',LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,623)INCNUM,ZERO,(MHBOI+MHEOI)*HALF,
     &      (SCBOI+SCEOI)*HALF, (ERBOI+EREOI)*HALF,  LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,620)INCNUM,'    ','M-H OUT','STC OUT','ERI OUT',LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,623)INCNUM,ZERO,
     &      OUTFLOW( (MHBOI+MHEOI)*HALF, (SCBOI+SCEOI)*HALF, MHK, MHYM,
     &      MHA, MHB, MHWT, MHC) - MHICE,
     &      OUTFLOW( (SCBOI+SCEOI)*HALF, (ERBOI+EREOI)*HALF, SCK, SCYM,
     &      SCA, SCB, SCWT, SCC) - SCICE,
     &      OUTFLOW( (ERBOI+EREOI)*HALF, CGMLV, ERK, ERYM,
     &      ERA, ERB, ERWT, ERC) - ERICE, LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C
C
      END
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C
      SUBROUTINE IterativeRoute( DT, IncNum )
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      DOUBLE PRECISION DT
C                   (Input) The number of seconds per increment.  Intrinsic
C                   part of the calculations, DT stands for "delta t". 
C
      INTEGER       IncNum
C                   (Input) The number of increments within each interval.
C                   Needed here only for logging purposes.
C
C
C
C     Declarations regarding local variables or functions:
C
C      INTEGER       SolvedSC, SolvedMH, SolvedER
      INTEGER       K
C                   These are scratch variables.
C
      DOUBLE PRECISION OUTMH, OUTSC, OUTER
C                   Aliased to values in FWRK for economy of storage
C                   and passing values to/from calling routine.
C
      DOUBLE PRECISION MHBOI, SCBOI, ERBOI
      DOUBLE PRECISION MHNTS, SCNBS, ERNBS
      DOUBLE PRECISION MHICE, SCICE, ERICE
C                   Input values assigned by calling routine.
C                   Aliased to values in FWRK for economy of storage
C                   and passing values to/from calling routine.
C                   MHBOI, etc are beginning of increment levels.
C                   MHICE, etc are retardation values for the period.
C                   MHNTS, SCNBS, ERNBS are summarized inflow/outflow
C                   for each lake (MH includes St. Marys flow, hence NTS).
C
      DOUBLE PRECISION CGMLV
C                   Input value assigned by calling routine.
C                   Mean level for current routing interval, used for 
C                   two-gage equation on Niagara River.  Note that this
C                   storage is actually aliased to memory in FWRK.
C                   
      DOUBLE PRECISION MHIMLV, SCIMLV, ERIMLV
C                   Intermediate working variables used in this subroutine, 
C                   holds the mean level for the increment.  Aliased to 
C                   memory in FWRK.
C
      DOUBLE PRECISION LastMHEOI, LastSCEOI, LastEREOI 
      DOUBLE PRECISION LastMHOF, LastSCOF, LastEROF 
C                   Intermediate working variables used in this subroutine, 
C                   holds the outflow and end-of-increment level computed 
C                   in the previous iteration.  Aliased to memory in FWRK.
C
      DOUBLE PRECISION DiffMH, DiffSC, DiffER
C                   Intermediate working variables used in this subroutine, 
C                   holds the difference in outflow  between iterations.
C                   Aliased to memory in FWRK.
C
      DOUBLE PRECISION MHEOI, SCEOI, EREOI
C                   Output end-of-increment values calculated by this routine.
C                   Aliased to values in FWRK for economy of storage
C                   and passing values to/from calling routine.
C
      DOUBLE PRECISION GETAREA
C                   Function that returns surface area for a lake
C
      DOUBLE PRECISION OUTFLOW
C                   Function that returns connecting channel flow, given 
C                   the parameters for the discharge equation.
C
      DOUBLE PRECISION Round
C                   Function that rounds according the old hand method.
C
      DOUBLE PRECISION ToleranceMH, ToleranceSC, ToleranceER
C                   Threshold at which the difference in outflow between
C                   two iterations becomes close enough, in m3s.
C 
      DOUBLE PRECISION LevelToInteger 
C
      INTEGER       EOILevels2Ago(3)
C
      INTEGER       MaxIterations
C
      LOGICAL       Converge
C                   
C 
C
C     Declarations regarding COMMON block variables:
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                   See description in Block Data Genral
C
      CHARACTER     LLT*1
C                   This is the character that marks the end of a line sent
C                   to the LogIt subroutine.
C 
      CHARACTER     CTMP*255, CWRK*80
C                   These are just buffers used throughout program.
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
      INTEGER       TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &              MissingVal, SecsInUniformMon
      DOUBLE PRECISION MPFT, HALF, ONE, ZERO
C                   See descriptions in Block Data General.  
C
      INTEGER       MidLakeSolution, StMarysSolution, UStLawrSolution,
     &              P77ACompat, BOMRound, 
     &              RoutingInterval(3), RoutingIncrement(3)
C                   See descriptions in Block Data Genral.  Just here
C                   because of other members of MISINT
C
      INTEGER       FlowRound, LevelRound
C                   These values determine how to round the flows and levels 
C                   during the internal routing and regulation computations.
C                   They correspond to powers of ten; for instance, a value
C                   of 3 means round to the thousand, while a value of -2 means
C                   round to the hundredth.  The allowable values are integers
C                   from -3 to 3.
C
C
      DOUBLE PRECISION MHK,MHYM,MHA,MHB,MHWT,MHC
      DOUBLE PRECISION SCK,SCYM,SCA,SCB,SCWT,SCC
      DOUBLE PRECISION ERK,ERYM,ERA,ERB,ERWT,ERC
C                   See descriptions in CHECKPAR and OUTFLOW
C
      DOUBLE PRECISION FWRK
C                   This is a common block used for 64-bit calculations.
      DIMENSION FWRK(225,12)
C                   Dimensioned it this way to provide a 225 element array for
C                   each pool lake, and some scratch.
C       
C
C     Common block statements.

      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /MISINT/   MidLakeSolution, BOMRound, FlowRound, 
     &                  LevelRound, RoutingInterval, RoutingIncrement,
     &                  StMarysSolution, UStLawrSolution, P77ACompat
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /STAGE/    MHK,SCK,ERK,MHYM,SCYM,ERYM,MHA,SCA,ERA,
     &                  MHB,SCB,ERB,MHWT,SCWT,ERWT,MHC,SCC,ERC
      COMMON /FSPACE/   FWRK
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
C        
C         
C     Assign convenient aliases for summary inflow values, which
C     where copied into FRWK(1,?) in calling subroutine.  Remember that
C     the M-H summary inflow already includes St. Marys River contribution.
C     Also recognize that consumptive use and diversions, etc, are already
C     combined into each of these terms, which is a bit of a change from
C     GLERL's original approach.
      EQUIVALENCE (MHNTS, FWRK(1,1))
      EQUIVALENCE (SCNBS, FWRK(1,2))
      EQUIVALENCE (ERNBS, FWRK(1,3))
C       
C     Assign convenient aliases for summary ice/weed retardation values, 
C     which where copied into FRWK(2,?) in calling subroutine.
      EQUIVALENCE (FWRK(2,1), MHICE)
      EQUIVALENCE (FWRK(2,2), SCICE)
      EQUIVALENCE (FWRK(2,3), ERICE)
C        
C     Assign convenient aliases for beginning of increment levels,
C     which where copied into FRWK(3,?) in calling subroutine.
      EQUIVALENCE (FWRK(3,1), MHBOI)
      EQUIVALENCE (FWRK(3,2), SCBOI)
      EQUIVALENCE (FWRK(3,3), ERBOI)
C        
C     Assign convenient aliases for end of increment levels, which will 
C     be written to FRWK(4,?) in at the end of this subroutine.
      EQUIVALENCE (FWRK(4,1), MHEOI)
      EQUIVALENCE (FWRK(4,2), SCEOI)
      EQUIVALENCE (FWRK(4,3), EREOI)
C
C     FWRK(5,?) reserved for period average levels in calling sub, 
C     except for CGMLV, which is actually input variable that gets
C     read instead of computed:
      EQUIVALENCE (CGMLV, FWRK(5,4))
C        
C     FWRK(6,?) reserved for period average outflows in calling subroutine.
C     FWRK(7,?) reserved for beginning of period levels in calling subroutine.
C     FWRK(8,?) reserved for end of period levels in calling subroutine.
C     FWRK(9,?) reserved for solution check values in calling subroutine.
C     DO NOT equivalence this memory in this function!
C
C
C
C     Assign convenient aliases for previous outflows and EOI, which will 
C     are used within this subroutine for intermediate values.  Last??EOI
C     only used for updating EOILevels2Ago(3), but could be used for 
C     convergence test instead of the Last??OF variables.
      EQUIVALENCE (FWRK(11,1), LastMHOF)
      EQUIVALENCE (FWRK(11,2), LastSCOF)
      EQUIVALENCE (FWRK(11,3), LastEROF)
      EQUIVALENCE (FWRK(14,1), LastMHEOI)
      EQUIVALENCE (FWRK(14,2), LastSCEOI)
      EQUIVALENCE (FWRK(14,3), LastEREOI)
      EQUIVALENCE (FWRK(15,4), LevelToInteger) 
C
C     Assign convenient aliases for difference from previous outflows, 
C     which will are used within this subroutine for intermediate values.
      EQUIVALENCE (FWRK(12,1), DiffMH)
      EQUIVALENCE (FWRK(12,2), DiffSC)
      EQUIVALENCE (FWRK(12,3), DiffER)
C
C     Assign convenient aliases for increment levels, which will 
C     are used within this subroutine for intermediate values.
      EQUIVALENCE (FWRK(13,1), MHIMLV)
      EQUIVALENCE (FWRK(13,2), SCIMLV)
      EQUIVALENCE (FWRK(13,3), ERIMLV)
C
C     Assign convenient aliases for the intermediate outflow terms, using 
C     storage actually allocated to the common workspace.
      EQUIVALENCE (FWRK(15,1), OUTMH)
      EQUIVALENCE (FWRK(15,2), OUTSC)
      EQUIVALENCE (FWRK(15,3), OUTER)
C
C
C
C
      DATA MaxIterations /9999/
C
      Converge = .TRUE.
C
C     Consider equations solved if iterate to one FlowRound (ie, if
C     if FlowRound set to 1, must be 10 m3s or less difference between
C     iterations).  If FlowRound == -99, use .1 m3s.
      ToleranceMH = 10.0D0**FlowRound
C
C     Hmmm, experimentally set tolerance to 1/10 of FlowRound, to make
C     sure get as many iterations as existing practice.
      ToleranceMH = 10.0D0**(FlowRound - 1)
C     But no less than 0.1 m3s
      IF( ToleranceMH .LT. 0.1D0 ) ToleranceMH = 0.1D0
      ToleranceSC = ToleranceMH
      ToleranceER = ToleranceSC
C
C
C     Initialize the variables that say how many iterations required for each lake:
C      SolvedMH = 0
C      SolvedSC = 0
C      SolvedER = 0
C
C     Initialize the outflows and EOI levels with bogus numbers
      OUTMH = -999.D0
      OUTSC = OUTMH
      OUTER = OUTSC
      MHEOI = OUTMH
      SCEOI = OUTMH
      EREOI = OUTMH
      EOILevels2Ago(1) = IDNINT(OUTMH)
      EOILevels2Ago(2) = EOILevels2Ago(1)
      EOILevels2Ago(3) = EOILevels2Ago(1)
C
C     Initialize the mean level values with the initial level
      MHIMLV = MHBOI 
      SCIMLV = SCBOI 
      ERIMLV = ERBOI 
C
C     Loop through so many iterations
      DO 1310, K = 1, MaxIterations
C      DO 1310, K = 1, 5
C        Save previous outflow for comparison to tolerance
         LastMHOF = OUTMH
         LastSCOF = OUTSC
         LastEROF = OUTER
C
         LastMHEOI = MHEOI
         LastSCEOI = SCEOI
         LastEREOI = EREOI

C        Approximate the outflows for the iteration
C
C         ===========nuke this section later===============
C         =====vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv=====
C         OUTMH = OUTFLOW( MHIMLV, SCIMLV, MHK, MHYM, MHA,
C     &      MHB, MHWT, MHC) - MHICE
C         OUTSC = OUTFLOW( SCIMLV, ERIMLV, SCK, SCYM, SCA, 
C     &      SCB, SCWT, SCC) - SCICE
C         OUTER = OUTFLOW( ERIMLV, CGMLV, ERK, ERYM, ERA,
C     &      ERB, ERWT, ERC) - ERICE
C         WRITE(CTMP,'(I3,A,I1,3(A,F15.8),A1)')INCNUM, ' Iteration ',K,
C     &       ' OutMH=', OutMH, ' OutSC=', OutSC,  ' OutER=', OutER, LLT
C         CALL LogIt(LOGD,0,0)
C         OUTMH = Round( OUTFLOW( MHIMLV, SCIMLV, MHK, MHYM, MHA,
C     &      MHB, MHWT, MHC) - MHICE, FlowRound ) 
C         OUTSC = Round( OUTFLOW( SCIMLV, ERIMLV, SCK, SCYM, SCA, 
C     &      SCB, SCWT, SCC) - SCICE, FlowRound ) 
C         OUTER = Round( OUTFLOW( ERIMLV, CGMLV, ERK, ERYM, ERA,
C     &      ERB, ERWT, ERC) - ERICE, FlowRound ) 
C         WRITE(CTMP,'(I3,A,I1,2(A,F18.2),A,F18.15,A1)')INCNUM,
C     &      ' Iteration ',K,
C     &      ' MHNTS=',MHNTS,' factor=', 1/(DT /GETAREA(2,MHIMLV)),
C     &      ' delta=',(MHNTS - OUTMH) * DT /GETAREA(2,MHIMLV),LLT
C         CALL LogIt(LOGD,0,0)
C         WRITE(CTMP,'(I3,A,I1,3(A,F18.2),A1)')INCNUM, ' Iteration ',K,
C     &      ' MHNTS=',MHNTS,' GETAREA(2,MHIMLV)=',GETAREA(2,MHIMLV),
C     &      ' factor=', 1/(DT /GETAREA(2,MHIMLV)),LLT
C         CALL LogIt(LOGD,0,0)
C         MHEOI = MHBOI + (MHNTS - OUTMH) * DT /GETAREA(2,MHIMLV)
C         SCEOI = SCBOI + (SCNBS + OUTMH - OUTSC) * DT/GETAREA(5,SCIMLV)
C         EREOI = ERBOI + (ERNBS + OUTSC - OUTER) * DT/GETAREA(6,ERIMLV)
C         WRITE(CTMP,'(I3,A,I1,3(A,F15.8),A1)')INCNUM, ' Iteration ',K,
C     &       ' MHEOI=', MHEOI, ' SCEOI=', SCEOI,  ' EREOI=', EREOI, LLT
C         CALL LogIt(LOGD,0,0)
C         =====^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^====
C         ===========nuke this section later===============
C
         OUTMH = Round( OUTFLOW( MHIMLV, SCIMLV, MHK, MHYM, MHA,
     &      MHB, MHWT, MHC) - MHICE, FlowRound ) 
         OUTSC = Round( OUTFLOW( SCIMLV, ERIMLV, SCK, SCYM, SCA, 
     &      SCB, SCWT, SCC) - SCICE, FlowRound ) 
         OUTER = Round( OUTFLOW( ERIMLV, CGMLV, ERK, ERYM, ERA,
     &      ERB, ERWT, ERC) - ERICE, FlowRound ) 
C
C        Compute EOP levels for the iteration
         MHEOI = Round( MHBOI + (MHNTS - OUTMH) * DT /GETAREA(2,MHIMLV),
     &      LevelRound )
         SCEOI = Round( SCBOI + (SCNBS + OUTMH - OUTSC) * DT /
     &      GETAREA(5,SCIMLV), LevelRound )
         EREOI = Round( ERBOI + (ERNBS + OUTSC - OUTER) * DT /
     &      GETAREA(6,ERIMLV), LevelRound )
C        WRITE(CTMP,11111) INCNUM, K,' ERNBS=', ERNBS, ' OUTSC=',
C    &      OUTSC, ' OUTER=', OUTER, LLT
C        CALL LogIt(LOGD,0,0)
C        WRITE(CTMP,'(I8,I14,A,F15.6,A,F15.2,A1)') INCNUM, K,
C    &      ' DT=', DT, ' Area=', GETAREA(6,ERIMLV), LLT
C        CALL LogIt(LOGD,0,0)
C        Compute mean levels for the iteration
         MHIMLV = (MHEOI + MHBOI) * HALF
         SCIMLV = (SCEOI + SCBOI) * HALF
         ERIMLV = (EREOI + ERBOI) * HALF
C
11111    FORMAT(' Inc ', I3,' Iteration',I4,3(A,F15.8),A1)
         IF( MLVerbosity .GE. 9) THEN
            WRITE(CTMP,11111) INCNUM,
     &         K,' OutMH=',OutMH, ' OutSC=', OutSC,' OutER=', OutER, LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,11111) INCNUM, 
     &         K,' MHEOI=', MHEOI,' SCEOI=', SCEOI, ' EREOI=', EREOI,LLT
C----&------------------------------------------------------------------
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,11111) INCNUM, 
     &         K,' MHIMLV=',MHIMLV,' SCIMLV=',SCIMLV,' ERIMLV=', ERIMLV,
     &         LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
C
C        Compute Diff values
         DiffMH = DABS(OUTMH - LastMHOF)
         DiffSC = DABS(OUTSC - LastSCOF)
         DiffER = DABS(OUTER - LastEROF)
C        the lines with variables SolvedMH, etc, establish that the
C        hardest lake to solve for is lake st. clair (mh and erie
C        converge before it does)
C         IF( SolvedMH .EQ. 0 .AND. DiffMH .LT. ToleranceMH) SolvedMH=K
C         IF( SolvedSC .EQ. 0 .AND. DiffSC .LT. ToleranceSC) SolvedSC=K
C         IF( SolvedER .EQ. 0 .AND. DiffER .LT. ToleranceER) SolvedER=K
C
         IF( DiffMH .LE. ToleranceMH .AND. DiffSC .LE. ToleranceSC
     &      .AND. DiffER .LE. ToleranceER ) GOTO 1350
C
         IF (LevelRound .NE. -99) THEN
C           Check to see if oscillating between two consecutive values
C           of rounded EOI levels.
C    
C 101       FORMAT(' Inc ', I3,' Iteration',I4,3(A,I10),A1)
C           WRITE(CTMP,101) INCNUM, K, ' MH 2 iter ago=', 
C    &         EOILevels2Ago(1),' SC 2 iter ago=', EOILevels2Ago(2),
C    &         ' ER 2 iter ago=', EOILevels2Ago(3), LLT
C           CALL LogIt(LOGD,0,0)
C                                                                                                          
C           WRITE(CTMP,'(3(A,F11.5,I10),A1)')
C    &        ' LastMHEOI= ',LastMHEOI,IDNINT(LastMHEOI*LevelToInteger),
C    &        ' LastSCEOI= ',LastSCEOI,IDNINT(LastSCEOI*LevelToInteger),
C    &        ' LastEREOI= ',LastEREOI,IDNINT(LastEREOI*LevelToInteger),
C    &        LLT
C           CALL LogIt(LOGD,0,0)
C 
            IF(EOILevels2Ago(1) .EQ. IDNINT(MHEOI*LevelToInteger) .AND.
     &         EOILevels2Ago(2) .EQ. IDNINT(SCEOI*LevelToInteger) .AND. 
     &         EOILevels2Ago(3) .EQ. IDNINT(EREOI*LevelToInteger) )THEN
C
               WRITE(CTMP,'(4X,A,I4,A,A1)') 
     &            '. . . solution converges at ', K,
     &            ' iterations, but desired tolerance not met!', LLT
               CALL LogIt(LOGD,0,0)
               GOTO 1350
C
C
            ENDIF            
C            
C           Convert to integers, with ones digit at rounding point, 
C           and update EOILevels2Ago.
            EOILevels2Ago(1) = IDNINT(LastMHEOI*LevelToInteger)
            EOILevels2Ago(2) = IDNINT(LastSCEOI*LevelToInteger)
            EOILevels2Ago(3) = IDNINT(LastEREOI*LevelToInteger)
         ENDIF            
C
 1310    CONTINUE
C     
C----&------------------------------------------------------------------
         WRITE(CTMP,614)MaxIterations, INCNUM, LLT
         CALL LogIt(LOGCD,0,0)
         Converge = .FALSE.
         MLVerbosity = 9
C
 1350 IF( MLVerbosity .GE. 8) THEN
         WRITE(CTMP,605)INCNUM, ' Beginning of',  
     &      MHBOI, SCBOI, ERBOI, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,605)INCNUM, '       End of',  
     &      MHEOI, SCEOI, EREOI, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,605)INCNUM, '         Mean', 
     &      MHIMLV, SCIMLV, ERIMLV, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,606)INCNUM, '  Ave outflow during',
     &      OUTMH, OUTSC, OUTER, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,606)INCNUM, '      Final diff for', 
     &      DiffMH, DiffSC, DiffER, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,612)INCNUM, K, ToleranceMH, LLT
         CALL LogIt(LOGD,0,0)
C         WRITE(CTMP,615)INCNUM, SolvedMH, SolvedSC, SolvedER, LLT
C         CALL LogIt(LOGD,0,0)
      ENDIF
C
C
      IF(Converge) RETURN
C
      CALL ErrMess(1091,'IterativeRoute'//LLT//' ',LOGCDE)
      CALL BuildInfo
      CALL LogIt(LOGCDE,0,0)
      CALL SAVETOD
      CALL ADIOS()
C
  605 FORMAT(I3, A,' increment levels (MH,SC,ER): ',3F14.8,A1)
  606 FORMAT(I3, A,' increment (MH,SC,ER): ',3F14.5, A1)
  612 FORMAT(I3, ' Required ', I3, 
     &   ' iterations to converge within ',F7.3,' m3s.',A1)
  614 FORMAT(' Exceeded max of ', I5,
     &   ' iterations for increment ',I2,A1)
C  615 FORMAT(I3, ' MH needed ', I3, ',  SC needed ', I3, 
C     &   ',  and ER needed ', I3, ' iterations.',A1)
C
C
C
      END
C
C
C
C
C
      SUBROUTINE RouteMiddleLakes(SolveErrors)
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER       SolveErrors       
C                   (Output)  Indicates how many increments the 
C                   simultaneous solution of the routing equations 
C                   failed to check out regarding continuity.
C 
C     Declarations regarding local variables or functions:
C
      INTEGER I
C                   Scratch variables
      DOUBLE PRECISION MHBOP, MHEOP, SCBOP, SCEOP, ERBOP, EREOP
C                   Beginning- and end-of-period levels for each lake.
C                   Note that these storages are actually aliased to
C                   memory in FWRK.
C
      DOUBLE PRECISION MHBOI, MHEOI, SCBOI, SCEOI, ERBOI, EREOI
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
      DOUBLE PRECISION MHForceEOPLev, STCForceEOPLev, ErieForceEOPLev
C                   EOP levels to force, after routing through 
C                   current period.  No effect on current period.  
C                   Simply replace BOP levels when roll over into
C                   next period.  Had to implement it this way (instead
C                   of just forcing BOP level for current period) 
C                   so that Plan 1977A would be aware of forced BOP
C                   levels on the middle lakes.
C                   
      DOUBLE PRECISION MHNTS, SCNBS, ERNBS
      DOUBLE PRECISION MHICE, SCICE, ERICE
C                   Values read from summary input files for net supply
C                   and ice/weed retardation.  Note that these storages
C                   are actually aliased to memory in FWRK.
C
      DOUBLE PRECISION MHMLV, SCMLV, ERMLV, CGMLV
C                   Period-average levels for each lake.  
C                   Note that these storages are actually aliased to
C                   memory in FWRK.
C
      DOUBLE PRECISION USSlip
C                   Level computed at US Slip for this interval
C
      DOUBLE PRECISION StMarysConverted
C                   Flow from St. Marys River for this middle lakes 
C                   routing period.  Note that this storages is 
C                   actually aliased to memory in FWRK.
C
      DOUBLE PRECISION STCRIVFLOW, DETRIVFLOW, NIARIVFLOW, STMRIVFLOW
C                   Flows computed from the period-average level of each
C                   lake. Note that these storages are actually aliased
C                   to memory in FWRK.
C
      DOUBLE PRECISION MLErrorFlow
C                   Maximum error (in terms of river flows) seen during
C                   continuity checks. 
C                   
      DOUBLE PRECISION MHOutI, SCOutI, EROutI
C                   Flows computed from the mean level of the current
C                   increment for each lake.
C                   Note that these storages are actually aliased
C                   to memory in FWRK.
C                   
      DOUBLE PRECISION MHSolutionCheck, SCSolutionCheck, 
     &              ERSolutionCheck, SolutionTest
C                   Variables used to verify that the simultaneous
C                   solution of the  routing equations was accurate.
C                   Note that these storages are actually aliased to
C                   memory in FWRK.
C
C
C      CHARACTER     CHNKST*16
C                   This is the calling subroutine name passed to
C                   ErrMess() in case of error messages.
C
      DOUBLE PRECISION OUTFLOW
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
      DOUBLE PRECISION GetLevelAtUSSlip 
C                   Returns US Slip level corresponding to given 
C                   St Marys flow and MH elevation
C
      CHARACTER*1   EOPType
C                	  Returns 'M' if period ends a month, or 'Q' if period
C                   ends a quarter-month, else returns 'W', or 'D'.
C
C
C      CHARACTER*16  MakeIntervalName
C                   Returns a string with the interval name.  For output.
C 
C     Declarations regarding COMMON block variables:
C
      INTEGER       PeriodEnd(6,3), Day(3), Q(4,3), DaysSinceSaturday
      DOUBLE PRECISION dT(3)
      LOGICAL       NumIncsVary(3) 
      INTEGER       StMarysSolution, UStLawrSolution, P77ACompat
C                   See descriptions in RUNIT()
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
C                   Solution method (1=iterative, 2=matrix).
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
      CHARACTER*3   MonthAbbrev(12)
      INTEGER       MonthLength
      CHARACTER*8   MonthLengthName(2)
C                   See description in Block Data Genral
C
      CHARACTER*8   RCRUNTYP
C                   See description in Block Data Genral
C
      CHARACTER*256 FileNames(100)
C                   See descriptions in Block Data Genral
C     
      CHARACTER*256 RepeatLines(80)
C                   See description in Block Data Genral
C
      DOUBLE PRECISION MHK,MHYM,MHA,MHB,MHWT,MHC
      DOUBLE PRECISION SCK,SCYM,SCA,SCB,SCWT,SCC
      DOUBLE PRECISION ERK,ERYM,ERA,ERB,ERWT,ERC
C                   See descriptions in CHECKPAR
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                   See description in Block Data Genral
C
      CHARACTER     CTMP*255, CWRK*80
C                   These are just buffers used throughout program.
C
      CHARACTER LLT*1
C                   This is the character that marks the end of a line
C                   sent to the LogIt subroutine.
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
      CHARACTER*21  LAKENAME(16), RiverName(16)
C                   See description in Block Data Genral
C
      INTEGER       NLN(16), LkToSlv2LN(8), NRN(16), LkSlv2RtLoopI(8)
C                   See description in Block Data Genral
C
      DOUBLE PRECISION FWRK
C                   This is a common block used for 64-bit calculations.
      DIMENSION FWRK(225,12)
C                   Dimensioned it this way to provide a 225 element
C                   array for each pool.
C       
      INTEGER       TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &              MissingVal, SecsInUniformMon
      DOUBLE PRECISION MPFT, HALF, ONE, ZERO
C                   See descriptions in Block Data General.  
C
      CHARACTER*128 RCOUTDIR
C                   See description in Block Data General.  
C
      CHARACTER*16  RCOUTEXT
C                   See description in Block Data General.  
C
      INTEGER       USSlipEquation
C                   1=ancient, 2=2010
      DOUBLE PRECISION  USSIceWeedAdj
C
C
C  Common block statements.
C       
      COMMON /CLDATE/   dT, Day, Q, PeriodEnd, NumIncsVary,
     &                  DaysSinceSaturday
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /MONAME/   MonthAbbrev
      COMMON /STIME/    STIME,SDATE
      COMMON /TITLES/   RCTITLE, MLHEADER
      COMMON /RCRUNI/   LakesToSolve, 
     &                  IntervalNameLen, MonthLength,
     &                  StartYear, StartMonth, StartDay, 
     &                  EndYear, EndMonth, EndDay
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCRUNC/   RCRUNTYP, IntervalName, MonthLengthName,
     &                  MLSolutionName, SMSolutionName, USLSolutionName
      COMMON /STAGE/    MHK,SCK,ERK,MHYM,SCYM,ERYM,MHA,SCA,ERA,
     &                  MHB,SCB,ERB,MHWT,SCWT,ERWT,MHC,SCC,ERC
      COMMON /WKCHAR/   CTMP,LLT,CWRK
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
      COMMON /OUTNAM/   RCOUTDIR, RCOUTEXT
      COMMON /USSLIP/   USSIceWeedAdj, USSlipEquation
C----&------------------------------------------------------------------
C
C
C
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
C     No NBS applicable to CGIP, so leave FWRK(1,4) empty
C      EQUIVALENCE (SPNTS, FWRK(1,5))
C       
C     Assign convenient aliases for summary ice/weed retardation values, 
C     which will be copied into FRWK(2,?) in this subroutine.
      EQUIVALENCE (FWRK(2,1), MHICE)
      EQUIVALENCE (FWRK(2,2), SCICE)
      EQUIVALENCE (FWRK(2,3), ERICE)
C     No applicable term for CGIP; storage in FWRK(2,4) unused.
C      EQUIVALENCE (FWRK(2,5), SMICE)
C
C     Assign convenient aliases for beginning of increment levels,
C     which will be copied into FRWK(3,?) in this subroutine.
      EQUIVALENCE (FWRK(3,1), MHBOI)
      EQUIVALENCE (FWRK(3,2), SCBOI)
      EQUIVALENCE (FWRK(3,3), ERBOI)
C     No applicable term for CGIP; use the storage for NIZ
      EQUIVALENCE (FWRK(3,4), NIZ)
C      EQUIVALENCE (FWRK(3,5), SPBOI)
C        
C     Assign convenient aliases for end of increment levels, which will 
C     be written to FRWK(4,?) at the end of the routing subroutine.
      EQUIVALENCE (FWRK(4,1), MHEOI)
      EQUIVALENCE (FWRK(4,2), SCEOI)
      EQUIVALENCE (FWRK(4,3), EREOI)
C     No applicable term for CGIP; storage in FWRK(4,4) unused.
C      EQUIVALENCE (FWRK(4,5), SPEOI)
C
C     Assign convenient aliases for the period average levels, using 
C     storage actually allocated to the common workspace.  The levels
C     will be initialized in FRWK(5,?) and computed by this subroutine.
      EQUIVALENCE (MHMLV, FWRK(5,1))
      EQUIVALENCE (SCMLV, FWRK(5,2))
      EQUIVALENCE (ERMLV, FWRK(5,3))
C     CGMLV is actually input variable; gets read instead of computed.
      EQUIVALENCE (CGMLV, FWRK(5,4))
C      EQUIVALENCE (SPMLV, FWRK(5,5))
C
C     Assign convenient aliases for the period final outflow terms, 
C     using storage actually allocated to the common workspace.  The 
C     flows will be copied into FRWK(6,?) in this subroutine.
      EQUIVALENCE (STCRIVFLOW, FWRK(6,1))
      EQUIVALENCE (DETRIVFLOW, FWRK(6,2))
      EQUIVALENCE (NIARIVFLOW, FWRK(6,3))
C     No term for CGIP; use storage in FWRK(6,4) for MLErrorFlow.
      EQUIVALENCE (MLErrorFlow, FWRK(6,4))
      EQUIVALENCE (STMRIVFLOW, FWRK(6,5))
C        
C     Assign convenient aliases for beginning of period levels,
C     which will be copied into FRWK(7,?) in this subroutine.
      EQUIVALENCE (FWRK(7,1), MHBOP)
      EQUIVALENCE (FWRK(7,2), SCBOP)
      EQUIVALENCE (FWRK(7,3), ERBOP)
C     This implementation assumes that the CGIP levels are constant 
C     throughout the routing interval, so that the beginning-of-
C     period, the end-of-period, and period average values are all
C     equivalent.  So no need to specify a BOP variable, and
C     leaving storage in FWRK(7,4) unused.
C     EQUIVALENCE (FWRK(7,4), CGBOP)
C      EQUIVALENCE (FWRK(7,5), SPBOP)
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
C     leaving storage in FWRK(8,4) unused.
C     EQUIVALENCE (FWRK(8,4), CGEOP)
C      EQUIVALENCE (FWRK(8,5), SPEOP)
C
C     Assign convenient aliases for solution check values, which will 
C     be written to FRWK(9,?) after call to a routing subroutine.
      EQUIVALENCE (FWRK(9,1), MHSolutionCheck)
      EQUIVALENCE (FWRK(9,2), SCSolutionCheck)
      EQUIVALENCE (FWRK(9,3), ERSolutionCheck)
C     No applicable term for CGIP, so use the storage for SolutionTest
      EQUIVALENCE (FWRK(9,4), SolutionTest)
C      EQUIVALENCE (FWRK(9,5), SPSolutionCheck)
C
C     FWRK(10,?), FWRK(11,?), FWRK(12,?), FWRK(13,?), and FWRK(14,?)
C     get overwritten each increment in Middle Lakes routing.  Use
C     them for intermediate values if desired.
C
C     Assign convenient aliases for the intermediate outflow terms, using 
C     storage actually allocated to the common workspace.  These 
C     represent outflows per increment.
      EQUIVALENCE (FWRK(15,1), MHOutI)
      EQUIVALENCE (FWRK(15,2), SCOutI)
      EQUIVALENCE (FWRK(15,3), EROutI)
C      EQUIVALENCE (FWRK(15,5), SPOutI)
C
C     StMarysConverted contains the St. Marys Flow for this middle lake 
C     routing period.  
C     XX StMarysConverted may be undefined if not computing St. Marys River flows! 
C     XX would need to update ReadFormat around label 8010 of Control.for in
C     XX order to preserve St. Marys River flow as separate timeseries when not
C     XX solving for St Marys flows
      EQUIVALENCE (FWRK(17,5), StMarysConverted)
C
      EQUIVALENCE (FWRK(31,1), MHForceEOPLev)
      EQUIVALENCE (FWRK(31,2), STCForceEOPLev)
      EQUIVALENCE (FWRK(31,3), ErieForceEOPLev)
C     No applicable term for CGIP; storage in FWRK(31,4) unused.
C      EQUIVALENCE (FWRK(31,4), CGIPForceEOPLev)
C
      SolveErrors = 0
C
C        Initialize BOP levels as the EOP levels from previous period
         MHBOP = MHEOP
         SCBOP = SCEOP
         ERBOP = EREOP
C              
C            
      IF(MLVerbosity .GE. 7 ) THEN
         WRITE(CTMP,'(3X,6(A,F8.1),A1)') ' MHNTS=', MHNTS, 
     &      ' MHICE=', MHICE, ' SCNBS=', SCNBS, ' SCICE=', SCICE,
     &      ' ERNBS=',ERNBS,' ERICE=',ERICE, LLT
         CALL LogIt(LOGD,0,0)
         WRITE(CTMP,'(3X,4(A,F11.4),A1)') ' CGMLV=', CGMLV, 
     &      '  EOP levels to force: MH=', MHForceEOPLev, ' SC=',
     &      STCForceEOPLev, ' ER=', ErieForceEOPLev, LLT
         CALL LogIt(LOGD,0,0)
         IF(MLVerbosity .GE. 8 ) THEN
            WRITE(CTMP,'(4X,A,F13.0,2(A,F12.0), A1)') 
     &         'beginning of interval areas: Mich-Hur=', 
     &         GETAREA(2,MHBOP),        ' St. Clair=', GETAREA(5,SCBOP), 
     &         ' Erie=', GETAREA(6,ERBOP), LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,'(4X,A,I4,A,F12.1,A1)') 
     &         'Middle lakes increments per interval=',
     &         RoutingIncrement(2), ' dT=', dT(2), LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
      ENDIF
C
C        Initialize mean lake level registers for period.  The mean
C        level is computed as the average of BOI and EOI levels for
C        each increment.  The way the math works out, the BOI level 
C        for the first increment and the BOI level for the last
C        increment get weighted by half, while the EOI levels between
C        periods get full weights.  This is the concept that Matt
C        had such difficulty with. Anyway, for example, if 4 increments:
C
C        MLV= (BOI1+EOI1)/2 +(BOI2+EOI2)/2 +(BOI3+EOI3)/2 +(BOI4+EOI4)/2 
C        and BOI2 = EOI1, BOI3 = EOI2, BOI4 = EOI3
C        and BOI1 = BOP, EOI4 = EOP
C           so
C        MLV = 0.5 * BOP + EOI1 + EOI2 + EOI3 + EOI4 - 0.5 * EOP
C
         MHMLV = MHBOP * HALF
         SCMLV = SCBOP * HALF
         ERMLV = ERBOP * HALF
C
C        Initialize EOI values with BOP values; they will be rolled into BOI 
C        values at beginning of loop through increments.
         MHEOI = MHBOP
         SCEOI = SCBOP
         EREOI = ERBOP
C
C        Initialize period-average discharge values with 0.0D.  They will be
C        computed as the mean of the outflows computed for each increment.
C        (The outflow for each increment is calculated using the mean level
C        for the increment).
         STCRIVFLOW = ZERO
         DETRIVFLOW = ZERO
         NIARIVFLOW = ZERO
C
         DO 1222 I=1, RoutingIncrement(2)            
C
C           Initialize BOP levels as EOP levels from previous increment
            MHBOI = MHEOI
            SCBOI = SCEOI
            ERBOI = EREOI
C
            IF( MidLakeSolution .EQ. 2 )THEN 
C              Solve the simultaneous equations using GLERL's MIDLAKES
C              algorithm:
               CALL COEFFMATRIX(dT(2), I)
            ELSE
C              Solve simultaneous equations using iterative procedure
               CALL IterativeRoute(dT(2), I)
            ENDIF
C
C           Determine a mean level and calculate outflow for the
C           increment:
            MHOutI = OUTFLOW((MHBOI + MHEOI)*HALF, (SCBOI + 
     &         SCEOI)*HALF, MHK,MHYM,MHA,MHB,MHWT,MHC) - MHICE
            SCOutI = OUTFLOW((SCBOI + SCEOI)*HALF, (ERBOI + 
     &         EREOI)*HALF, SCK,SCYM,SCA,SCB,SCWT,SCC) - SCICE
            EROutI = OUTFLOW((ERBOI + EREOI)*HALF, CGMLV,
     &         ERK,ERYM,ERA,ERB,ERWT,ERC) - ERICE
C                
C           Accumulate outflow for increment, for later averaging
            STCRIVFLOW = STCRIVFLOW + MHOutI
            DETRIVFLOW = DETRIVFLOW + SCOutI
            NIARIVFLOW = NIARIVFLOW + EROutI
C
C           WRITE(CTMP,'(I3,16X,A3,3(A9,F13.6),A1)')I, '@@@', 'MHMLV=', 
C    &         MHMLV, 'SCMLV=', SCMLV, 'ERMLV=', ERMLV, LLT
C           CALL LogIt(LOGD,0,0)
            IF(MLVerbosity .GE. 5 ) THEN
               WRITE(CTMP, 8067) I, 'BOI','EOI','BOI',
     &            'EOI','BOI','EOI',LLT   
               CALL LogIt(LOGD,0,0)
C----&------------------------------------------------------------------
               WRITE(CTMP, 8063)I,PeriodEnd(1,2), PeriodEnd(2,2), 
     &            PeriodEnd(3,2), PeriodEnd(4,2),      
     &             MHBOI, MHEOI,SCBOI,SCEOI,ERBOI,EREOI,LLT  
               CALL LogIt(LOGD,0,0)
               WRITE(CTMP, '(I3,19X,6A12,A1)') I, 'MH inc mean', 
     &            'St. Clair R.', 'SC inc mean', 'Detroit R.',
     &            'ER inc mean', 'Niagara R.', LLT   
               CALL LogIt(LOGD,0,0)
C----&------------------------------------------------------------------
               WRITE(CTMP,'(I3,19X,6F12.5,A1)')I, (MHBOI + MHEOI)*HALF,
     &            MHOutI, (SCBOI + SCEOI)*HALF, SCOutI, 
     &            (ERBOI + EREOI)*HALF, EROutI, LLT
               CALL LogIt(LOGD,0,1)
            ENDIF
C                   
C           Now check that the simultaneous equations were solved
C           correctly.  Verify that the inflow plus NBS, minus the
C           average rate of change in storage during the period,
C           nearly equals the outflow for the period.  Store the 
C           difference in the SolutionCheck variables for debugging
C           info.  Note that this implemenation uses surface area as
C           BOI value.  Can easily fix this by introducing a new 
C           increment mean level if want to check variable surface 
C           area runs.
            ERSolutionCheck= ROUND(DABS( EROutI - SCOutI -
     &         ERNBS + ((EREOI - ERBOI) * GETAREA(6,EREOI))/dT(2)),
     &         FlowRound)
             SCSolutionCheck= ROUND(DABS( SCOutI - MHOutI -
     &         SCNBS + ((SCEOI - SCBOI) * GETAREA(5,SCEOI))/dT(2)),
     &         FlowRound)
            MHSolutionCheck = ROUND(DABS( MHOutI -    
     &         MHNTS + ((MHEOI - MHBOI) * GETAREA(2,MHEOI))/dT(2)),
     &         FlowRound)
            IF(ERSolutionCheck .GT. SolutionTest .OR. 
     &         SCSolutionCheck .GT. SolutionTest .OR. 
     &         MHSolutionCheck .GT. SolutionTest  ) THEN
               SolveErrors = SolveErrors + 1
               IF(MLVerbosity .GT. 5 ) THEN
                  WRITE(CTMP, 8072) I, PeriodEnd(1,2), 
     &               MonthAbbrev(PeriodEnd(3,2)), PeriodEnd(4,2), LLT
                  CALL LogIt(LOGD,0,0)
                  WRITE(CTMP,8033)PeriodEnd(5,2),I,'MH',MHSolutionCheck,
     &               'SC', SCSolutionCheck,'ER', ERSolutionCheck, LLT
                  CALL LogIt(LOGD,0,0)
               ENDIF
               IF(ERSolutionCheck .GT. MLErrorFlow)
     &            MLErrorFlow = ERSolutionCheck
               IF(SCSolutionCheck .GT. MLErrorFlow)
     &            MLErrorFlow = SCSolutionCheck
               IF(MHSolutionCheck .GT. MLErrorFlow)
     &            MLErrorFlow = MHSolutionCheck
            ELSE
               IF(MLVerbosity .GE. 7 ) THEN
                  WRITE(CTMP, 8071)  I, PeriodEnd(5,2), LLT
                  CALL LogIt(LOGD,0,1)
               ENDIF
            ENDIF
C
C           In case using two-gage Niagara eqn, save current mean level
C           in NIZ for use next time as BOI level:
            NIZ = CGMLV
C
C           Accumulate the EOI levels into the mean level variables
            MHMLV = MHMLV + MHEOI
            SCMLV = SCMLV + SCEOI
 1222       ERMLV = ERMLV + EREOI
C            
C        Ok, have looped through each increment of a period.  The EOP 
C        values are the final EOI levels.
         MHEOP = MHEOI
         SCEOP = SCEOI
         EREOP = EREOI
C           
C        Must adjust the mean level values to subtract half of the final
C        EOI values (was only supposed to add a half inside the
C        loop, but added a whole).  Then divide by number of increments.
C           WRITE(CTMP,'(I3,16X,A3,3(A9,F18.11),A1)')I, '@@@', 'MHMLV=', 
C    &         MHMLV, 'SCMLV=', SCMLV, 'ERMLV=', ERMLV, LLT
C           CALL LogIt(LOGD,0,0)
C           WRITE(CTMP,'(I3,16X,A3,3(A9,F18.11),A1)')I,'@@@','MHEOP/2=', 
C    &         MHEOP * HALF, 'SCEOP/2=', SCEOP * HALF, 'EREOP/2=', 
C    &         EREOP * HALF, LLT
C           CALL LogIt(LOGD,0,0)
         MHMLV = (MHMLV - (MHEOP * HALF)) / DBLE(RoutingIncrement(2))
         SCMLV = (SCMLV - (SCEOP * HALF)) / DBLE(RoutingIncrement(2))
         ERMLV = (ERMLV - (EREOP * HALF)) / DBLE(RoutingIncrement(2))
C            WRITE(CTMP,'(I3,16X,A3,3(A9,F18.11),A1)')I, '@@@', 'MHMLV=', 
C     &         MHMLV, 'SCMLV=', SCMLV, 'ERMLV=', ERMLV, LLT
C            CALL LogIt(LOGD,0,0)
C
C        Now calculate mean outflow for the period.
         STCRIVFLOW = STCRIVFLOW / DBLE(RoutingIncrement(2))
         DETRIVFLOW = DETRIVFLOW / DBLE(RoutingIncrement(2))
         NIARIVFLOW = NIARIVFLOW / DBLE(RoutingIncrement(2))
C
C
C        When using U.S. Slip Equation 1, only default Ice & Weed 
C        adjustments are allowed
         IF (USSlipEquation .EQ. 1) THEN
            USSlip = GetLevelAtUSSlip (StMarysConverted, MHMLV, 
     &        USSlipEquation, PeriodEnd(3,2),.TRUE.)
         ELSE
            USSlip = GetLevelAtUSSlip (StMarysConverted, MHMLV, 
     &        USSlipEquation, PeriodEnd(3,2),.FALSE.)
         END IF
C
C        Send the results to the midlakes output file
         WRITE(CTMP, 8069) PeriodEnd(5,2), Day(2), PeriodEnd(1,2), 
     &      PeriodEnd(3,2), PeriodEnd(4,2), EOPType(2), 
     &      ROUND(USSlip, -3), ROUND(MHMLV,-3),
     &      ROUND(MHEOP,-3), ROUND(STCRIVFLOW, 0), ROUND(SCMLV,-3), 
     &      ROUND(SCEOP,-3), ROUND(DETRIVFLOW, 0), ROUND(ERMLV,-3), 
     &      ROUND(EREOP,-3), ROUND(NIARIVFLOW, 0), ROUND(CGMLV,-3), LLT
         CALL LogIt(LOGM,0,0)
         IF(MLVerbosity .GE. 3) THEN
             CALL LogIt(LOGD,0,0)
             CALL WriteHighPrecision(2, MHMLV, MHEOP, StCRivFlow)
             CALL WriteHighPrecision(5, SCMLV, SCEOP, DetRivFlow)
             CALL WriteHighPrecision(6, ERMLV, EREOP, NiaRivFlow)

         I = RoutingInterval(2)
         WRITE(CTMP, 5997) IntervalName(I)(1:IntervalNameLen(I)), 
     &      PeriodEnd(1,2), PeriodEnd(3,2), PeriodEnd(4,2), EOPType(2),
     &      USSlip, MHMLV, StMarysConverted, LLT
         CALL LogIt(LOGD,0,1)
C 
 5997    FORMAT('US Slip high-precision mean level, based on mean ',
     &      'Michigan-Huron level and mean St. Marys River flow for ',
     &      A, ' period ending', 2I3, I5, ' (EO', A1,
     &      '): ', 2F17.10,F15.5,A1)
C
         ENDIF
C
C
C     Force the EOP level, if necessary, for use as next period's BOP
C     level:
 8088 FORMAT(' Replacing modeled Lake ', A, ' EOP level (', F10.6,
     &   ') with forced level of ',F12.8, A1)
      IF( MHForceEOPLev .NE. MissingVal )THEN
         IF(MLVerbosity .GE. 3 ) THEN
            WRITE(CTMP,8088) 'Michigan-Huron', MHEOP, MHForceEOPLev, LLT
            CALL LogIt(LOGD,1,0)
         ENDIF
         MHEOP = MHForceEOPLev
      ENDIF
      IF( StCForceEOPLev .NE. MissingVal )THEN
         IF(MLVerbosity .GE. 3 ) THEN
            WRITE(CTMP,8088) 'St. Clair', SCEOP, StCForceEOPLev, LLT
            CALL LogIt(LOGD,1,0)
         ENDIF
         SCEOP = StCForceEOPLev
      ENDIF
      IF( ErieForceEOPLev .NE. MissingVal )THEN
         IF(MLVerbosity .GE. 3 ) THEN
            WRITE(CTMP,8088) 'Erie', EREOP, ErieForceEOPLev, LLT
            CALL LogIt(LOGD,1,0)
         ENDIF
         EREOP = ErieForceEOPLev
      ENDIF
C
C----&------------------------------------------------------------------
C
C
 8033 FORMAT(' Day number: ',I9,' increment:',I3,3(A5,'Check=',F9.1),A1)
 8063 FORMAT(I3,I8,I3,I3,I5,6F12.7,A1)
 8067 FORMAT(I3,2X,'EndDay  Q  M   Yr', 6A12, A1 ) 
 8069 FORMAT(I9,3(I3),I5,1X,A1,F10.4,3(F10.4,F9.4,F8.2),F11.4,A1)
 8071 FORMAT(' Successfully solved increment ', I3, 
     &   ' of routing interval ending on day ', I10, A1)
 8072 FORMAT(' Suspicious middle lake routing for increment',
     &   I4, ' of period ending ',I2,A,I5,A1)
C
      END


C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
C
      DOUBLE PRECISION FUNCTION GetLevelAtUSSlip (StMarysFlow, MHLevel, 
     &  USSlipEqArg, MonthIndex, UseDefaultAdj)
C
C
C     This function returns the water level at US Slip to the nearest 
C     1/10 mm that satisfies the US Slip equation, given the period- 
C     average St.Marys' River flow and period-average Mich-Huron level.
C
C     The IGLS1985 version of the historical equation:
C            Q = 347.3*(USS - 172.30)^1.5 * (USS - MH - 0.06)^0.4
C     The old equation gets used by default, or with this parameter:
C  USSlipEquation: 1
C  
C     In March 2010, added ability to use a new equation developed 
C     by Environment Canada in 2010.  
C     USS = (Q/(6.291*(MH-167.0)^2.816))^2.625 + MH + constant
C                  
C     The new equation is used when parameter set as follows:
C  USSlipEquation: 2
C
C     By default this form uses hard-coded, monthly parameters for
C     the adjustment constant at the end , but these can be
C     overrriden by time-series data using parameter like:
C UssIceWeedAdj: m, monthly, 12, 0.51, 0.52, 0.53, -1.4, 1.5, 1.0, -9999, 0.5, 0.5, 0.5, 0.5, 0.5
C
C     (in example above the value for July is missing so equation falls back on
C     default value of 0.004m)
C
C     **BE CAREFUL USING UssIceWeedAdj DURING LAKE SUPERIOR CALCULATIONS! **
C     When calling this routine from the Lake Superior module, ordinarily
C      
C     Data file values for UssIceWeedAdj only get updated when starting 
C     a new middle lakes routing period, which occurs after the Lake
C     Superior routing loop finishes.  If doing monthly Superior routing,
C     then during time that Superior module is active the values in 
C     UssIceWeedAdj will be for previous month.   Criterion B calculations
C     should use the hardcoded constant values anyway, so values from an
C     external file would only for some special analysis or plan formulation
C     purpose. Calls from the Superior module should set UseDefaultAdj to 
C     .TRUE. to ensure that the parameter input values are ignored.
C  
C     If lakes Michigan-Huron are too low, this function returns -9999
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C
      DOUBLE PRECISION StMarysFlow, MHLevel
      INTEGER MonthIndex
C      1=Jan, 12=Dec
      INTEGER USSlipEqArg
C      1=old eqn,  2=new eqn
      LOGICAL UseDefaultAdj
C      Only used if USSlipEqArg=2.  If true will use the
C      hardcoded defaults regardless of presence of data inputs.
C
C
C
      DOUBLE PRECISION testuss, testflow, lastuss, lastflow, radj
      DOUBLE PRECISION USSD(12)
      INTEGER I
C
      DOUBLE PRECISION USSIceWeedAdj
      INTEGER USSlipEquation

      INTEGER      TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &             MissingVal, SecsInUniformMon
      DOUBLE PRECISION MPFT, HALF, ONE, ZERO
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
      CHARACTER    CTMP*255, CWRK*80, LLT*1
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP, LLT, CWRK
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
C                   See descriptions in Block Data General.  
      COMMON /USSLIP/   USSIceWeedAdj, USSlipEquation
C
C     Initialize local variables      
C     Values for 2010 EC Thessalon Equation:
      DATA USSD /0.209D0, 0.310D0, 0.254D0, 0.062D0, -0.009D0, 
     &   -0.026D0, 0.004D0, 0.012D0, 0.006D0, 0.005D0, -0.012D0, 
     &   0.005D0/

C
      WRITE(CTMP, '(A, 2F9.3, A, I1, A, I3, A, L1, A1 )')
     &   ' In GetLevelAtUSSlip StMarysFlow and  MHLevel = ', 
     &   StMarysFlow, MHLevel, ' equation arg=', USSlipEqArg, 
     &   ' month=', MonthIndex, '  UseDefaultAdj=', UseDefaultAdj, LLT 
      CALL LogIt(LOGD,0,0)

C     XX May want to revisit this and make a message or fine tune this test.
C     XX US Slip must be at least 6 cm higher than MH for equation to 
C     XX arithmetically viable.  
C     XX Need to test low end of new equation!
      IF (MHLevel .LE.(172.30D0)) THEN
         GetLevelAtUSSlip = DBLE(-9999)
         RETURN
      ENDIF 
C     There must be some flow to proceed - arbitrarily say 0.5 m3/s
      IF (StMarysFlow .LE.(0.5D0)) THEN
         GetLevelAtUSSlip = DBLE(-9999)
         RETURN
      ENDIF 
C
      IF ((USSlipEqArg .EQ. 1) .AND. (UseDefaultAdj .EQV. .FALSE.)) THEN
         WRITE(CTMP, '(A, A1 )')
     &     ' Cannot use old US Slip equation with adjustment value',LLT
         CALL LogIt(LOGCDE,0,0)
         WRITE(CTMP, '(A, I2, A, L1, A1 )') ' USSlipEqArg=',USSlipEqArg,
     &     '  UseDefaultAdj=', UseDefaultAdj, LLT
         CALL LogIt(LOGCDE,0,0)
         STOP 'BadUSSEq2'
      ENDIF
C
      IF (USSlipEqArg .NE. USSlipEquation) THEN
C     Messy thing, but function started out only using argument to specify which equation 
C     to use, but brought in COMMON block to access time-series 
C     adjustment value, and now have redundant info.  For now just check 
C     to see if consistent.  If inconsisent bark but move on, maybe  
C     is a reason for different value.
         WRITE(CTMP, '(A, I1, A, I1, A1 )')
     &      ' In GetLevelAtUSSlip equation argument=', USSlipEqArg, 
     &      ' but input file says =', USSlipEquation, LLT 
         CALL LogIt(LOGCDE,0,0)
         STOP 'BadUSSEq3'
      ENDIF 
C
      IF (USSlipEqArg .EQ. 1) THEN
C        Using "old" US SLip calculations
         WRITE(CTMP, '(A, A1)') 
     &      ' In GetLevelAtUSSlip using old equation', LLT
         CALL LogIt(LOGD,0,0)
C        Do a first cut by incrementing provisional US Slip level by 10cm
C        until resulting flow exceeds the given St. Marys flow.
C
C        Initialize provisional US Slip test to minimum valid value, less
C        the first increment:
         testuss = MHLevel + 0.0601D0 - 0.1D0    
C
         DO 22, I = 1, 99999
            testuss = testuss + 0.1D0
            testflow =  347.3D0 * ((testuss - 172.30D0)**1.5D0) * 
     &        ((testuss - MHLevel - 0.06D0)**0.4D0)
   22       IF (testflow .GE. StMarysFlow) GOTO 25
C        This should not be possible
         STOP ' InfLoop1 USS'
C
C        Just using a CONTINUE for convenient jumping out of above loop
   25    CONTINUE
C        Now decrement provisional US Slip level by 0.0001 m until fall below
C        given river flow:
         DO 32, I = 1, 99999
            testuss = testuss - 0.0001D0
            testflow =  347.3D0 * ((testuss - 172.30D0)**1.5D0) * 
     &        ((testuss - MHLevel - 0.06D0)**0.4D0)
   32       IF (testflow .LE. StMarysFlow) GOTO 35
C        This should not be possible
         STOP ' InfLoop2 USS'
C        Just using a CONTINUE for convenient jumping out of above loop
   35    CONTINUE
C        Check to see if level from previous increment gave closer flow  
C        than the increment where result passed below the given flow.
         lastuss = testuss + 0.0001D0
         lastflow =  347.3D0 * ((lastuss - 172.30D0)**1.5D0) * 
     &     ((lastuss - MHLevel - 0.06D0)**0.4D0)
         IF(DABS(StMarysFlow - lastflow) .LT. StMarysFlow - testflow) 
     &     THEN
            GetLevelAtUSSlip = lastuss
         ELSE
            GetLevelAtUSSlip = testuss
         ENDIF
C
      ELSE IF (USSLipEqArg .EQ. 2) THEN
C
C        Use the 2010 EC (Thessalon) methodology
C        Set constant adjustment according to input options

         IF( (UseDefaultAdj .EQV. .FALSE.) .AND. 
     &      (IDNINT(USSIceWeedAdj) .NE. MissingVal) ) THEN
            radj = USSIceWeedAdj
            WRITE(CTMP, '(A, F12.7, A1)') 
     &         ' GetLevelAtUSSlip using new equation with input adj ',
     &         radj, LLT
            CALL LogIt(LOGD,0,0)
         ELSE
            radj = USSD(MonthIndex)
            WRITE(CTMP, '(A, F12.7, A1)') 
     &         ' GetLevelAtUSSlip using new equation with default adj ',
     &         radj, LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
         GetLevelAtUSSlip = (StMarysFlow/(6.291D0*(MHLevel-167.0D0)
     &             **2.816D0))**2.625D0 + MHLevel + radj
C
      ELSE
C        No appropriate US Slip Equation given as argument
         STOP ' Stewie Griffin'      
      END IF
C    
C

C From: Caldwell,Rob [Ontario] [mailto:Rob.Caldwell@ec.gc.ca] 
C Sent: Wednesday, March 03, 2010 5:27 AM (Davis)
C To: McPherson, Matthew M HEC; Fay,David [Ontario]
C Cc: Fan,Yin [Ontario]; Dahl, Travis A LRE; Koschik, John A LRE; Woodruff, Carl L LRE; Noorbakhsh, 
C Nanette L LRE; bill@werick.com
C Subject: RE: SWP and USS equations
C 
C Certainly understandable, Matt!  I guess my druthers would be option C.  To clarify, use the constants on 
C top of p. 9 and Table 3 with Equation 3 for MHU (green dots) and use constants on bottom of p. 9 and 
C Table 4 with Equation 4 for TSL (red).  What I recommend is always to use the latter, even if you must 
C feed it the MHU levels.  I was hoping to get some time soon to look at the impact of using one over the 
C other with MHU values, but it now looks like it'll be about a month before I can get to this.  Again, 
C as I elude to in the text, my feelings are that the latter option will still prove better.  If someone 
C else has some time to kill, feel free to check!
C 
C 
C Rob Caldwell
C Great Lakes - St. Lawrence Regulation Office | Bureau de la régularisation des Grands Lacs et du St-Laurent
C Environment Canada | Environnement Canada
C 111 Water St E | 111, rue Water est
C Cornwall ON  K6H 6S2
C Rob.Caldwell@ec.gc.ca
C Telephone | Téléphone  613-938-5725
C Facsimile | Télécopieur 613-937-1302
C 
C 
C 
C 
C -----Original Message-----
C From: McPherson, Matthew M HEC [mailto:Matthew.McPherson@usace.army.mil] 
C Sent: Tuesday, March 02, 2010 8:41 PM
C To: Caldwell,Rob [Ontario]; Fay,David [Ontario]
C Cc: Fan,Yin [Ontario]; Dahl, Travis A LRE; Koschik, John A LRE; Woodruff, Carl L LRE; Noorbakhsh, Nanette L LRE; bill@werick.com
C Subject: RE: SWP and USS equations
C 
C Sorry Rob, I got a little confused at the end of your new report.  For routing and regulation calculations 
C  am I supposed to: 
C 
C A) use the constants in Eq 3 for lake-wide mean Michigan-Huron levels as given at top of page 9 (green dots)
C B) use the constants Eq 4 for Thessalon as given at bottom of page 9, but feed it Michigan-Huron levels with 
C    the assumption that they are equivalent on a period-average basis (red dots - maybe?)
C C) make it yet another option for the modeler to choose either way
C 
C mmm
C 
C 
C -----Original Message-----
C From: Caldwell,Rob [Ontario] [mailto:Rob.Caldwell@ec.gc.ca]
C Sent: Friday, February 26, 2010 1:04 PM (Davis)
C To: McPherson, Matthew M HEC; Fay,David [Ontario]
C Cc: Fan,Yin [Ontario]; Dahl, Travis A LRE; Koschik, John A LRE; Woodruff, Carl L LRE; Noorbakhsh, Nanette L LRE; bill@werick.com
C Subject: RE: SWP and USS equations
C 
C Good afternoon, all:
C  
C I've finally taken the time to write up my efforts from 2000, 2002, and 2006 in one draft document (attached).  
C The final relationship is Equation 4 (at the end of the first report), and I would strongly recommend it's use, 
C even with lake-wide average MHU levels in place of Thessalon levels, as necessary.
C I drafted this report quite quickly, but it contains all relevant information you might desire.  Comments are 
C always appreciated and encouraged.
C  
C For your reference, I also include my original 2000 report for USS, as well as my final 2008 report for SWP.  
C  
C Thanks,
C 
C Rob Caldwell
C 
C   
C   Todo
C   - travis reconcile and review this version
C   - john will hand-check with spreadsheet
C   - rob will verify implementation about which equation to use 
C     with avg MH levels (ie., use eq 4 assumiong MHmean=Thessalon)
C   - matt will modify supcriterions to optionally use new eq if CritBMethod .NE. 1
C   - matt to refernce Rob's rport and explain implementation in IUGLS modeling description

      END
