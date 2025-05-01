C-----------------------------------------------------------------------
C
      LOGICAL FUNCTION InitRepeatData(TSDataIndex, 
     &   DataInterval, UnitTypeIndex, FloatArray)
C
C     This subroutine reads the lines from RepeatLines character array
C     that pertain to particular set of repeating time-series data.  The
C     element in FileNames given by TSDataIndex contains which 
C     entries in RepeatLines to use.  The first field on the first line
C     must specify the units of the data, with the time interval (daily, 
C     monthly, etc) as the second field.  Then the number of values to 
C     be read from rest of the line, followed by a comma-delimited list
C     fo the values.  Additional lines specify the number of values to
C     read from the line, followed by the list of values.  
C
C     ALL SUCH DATA IS REPEATED ANNUALLY.  Repeating data must be given
C     in blocks of an entire year.  This is largely for the programmer's
C     convenience, but in practice I doubt there is much need to repeat
C     fragments of years.  If you're doing a quarter-monthly routing
C     running from June to August, and you have repeating values for 
C     only that period, fine, just pad the leading and trailing values
C     with zero's or something.  It won't get used anyway.  If you
C     desire to repeat sequences stretching beyond one year, then you
C     should be using a file to specify your values.
C
C     So 12 values must be specified for monthly data, beginning with
C     January.  48 quarter-monthly values must be specified, again
C     beginning with the first quarter-month of the year.  Weekly and
C     daily values have not been addressed yet.  You can specify the 
C     values on one record (256 char limit), or spread them across up 
C     to *about* 20 records.  However, the total number of repeating
C     data records that can be accommodated in a run is limited to the
C     dimension of the RepeatLines array (currently 80).
C     
C     Note that there is a goofy workaround concerning use of CTMP.  Most
C     error messages are written to this buffer starting at element 150,
C     so that ErrMess() doesn't overwrite them, then are copied back so
C     as to start at element 1 after ErrMess is called.  Then LogIt()
C     works just fine.
C
      IMPLICIT NONE
C
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER      TSDataIndex
C                  (Input) File's element number in array FileNames.
C
C      INTEGER      BeginDate(3)
C                  Input. See description in FORMATINPUT()
C
      CHARACTER    DataInterval
C                  (Output) A letter designating the data's time interval.
C                  One of 'Q', 'M', 'W', 'D', or 'C'.
C
      INTEGER      UnitTypeIndex
C                  (Output) Element in array UnitType() of file's units.
C
      DOUBLE PRECISION FloatArray(50)
C                  (Output) This is where to store the values read from
C                  the lines in RepeatLines.                    
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER      I, J, K, L, M, N, H, FS, FE
C                  These are scratch variables.
C
      INTEGER      NumLines
C                  The number of lines in RepeatLines associated with
C                  the repeating time-series.  
C
      INTEGER      TotalValues
C                  The number of values required for an entire year.  
C                  12 for monthly, 48 for qmonthly.
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
C                  Just declared because it's in the COMMON block
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
      DOUBLE PRECISION UnitConvFactor(20)
C                  See descriptions in Block Data Genral
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
      COMMON /FNAMEC/   FileNames, RepeatLines
      COMMON /FNAMEI/   RepeatIndex, NPFILZ
      COMMON /UNITSC/   UnitType
      COMMON /UNITSI/   UnitTypeLen, BaseUnitType, UnitConvFactor
C
      CHNKST='InitRepeatData'//LLT//' '
C
 1333 FORMAT(I3)
      InitRepeatData = .FALSE.
C
      IF(TSDVerbosity .GE. 3) THEN
         WRITE(CTMP,'(A,I2,A,A1)') 
     &      ' Initializing repeating data for TimeSeries number ',
     &      TSDataIndex, ' . . .', LLT
         CALL LogIt(LOGD,1,0)
C        J = 16 + INDEX(FileNames(TSDataIndex)(16:256),LLT) - 1
C          WRITE(*,'(2(A,I3),4A)')'FileNames(', TSDataIndex, ')(1:',J,
C    &      ')={', FileNames(TSDataIndex)(1:J), '}', LLT
         IF(TSDVerbosity .GE. 9) THEN
            WRITE(CTMP,'(A,A,A,A1)') ' Number of elements ',
     &         'in RepeatLines array, and their indexes:',
     &         FileNames(TSDataIndex)
     &         (16:16+INDEX(FileNames(TSDataIndex)(16:256),LLT) - 1),LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
      ENDIF
C
C     Get the number of lines to process.
      READ(FileNames(TSDataIndex)(16:18), 1333) NumLines
C      IF(NumLines .LT. 1) STOP 'NumLines < 1'
C
C     Get the index number to the first line in RepeatLines:
      READ(FileNames(TSDataIndex)(19:21), 1333) I
C      IF(I .LT. 1) STOP 'impossible'
C
C     Ok, now try to read the units.  First get the field bounds.
C     Initialize the current position on RepeatLines(I) in J:
      J = 1
C     Get beginning and end of field in FS and FE:
      IF( .NOT. FIELDBND(RepeatLines(I), 256, J, FS, FE, M ) ) THEN 
         CALL ErrMess(1053, CHNKST, LOGCDE)
C        This is pretty skimpy message, but should never happen, 
C        since a comma is needed in INIT.FOR to get here.
         RETURN
      ENDIF 
C     Save the beginning of the next field, if any, in J:
      J = M + 1
C      WRITE(*,'(A)') ' ##'//RepeatLines(I)(FS:FE)//'##'
C
C     Convert to uppercase
      L = 0
      DO 40, M = FS, FE
         K = ICHAR( RepeatLines(I)(M:M) )
         IF( K .GE. 97 .AND. K .LE. 122) K = K-32
         IF((K .GE. 48 .AND. K .LE. 57) .OR.
     &      (K .GE. 65 .AND. K .LE. 90)) THEN
            L = L + 1
            CWRK(L:L) = CHAR(K)
         ENDIF
   40    CONTINUE
C
C     Now search through the list of units (except #20)
      DO 311, M=1,19
C         WRITE(*,'(1X,I2,A)') M,' Does '//CWRK(1:L)//
C     &      ' .EQ. '//UnitType(M)(1:UnitTypeLen(M))//'?'
  311    IF( CWRK(1:L) .EQ. UnitType(M)(1:UnitTypeLen(M))) GOTO 312
C
      IF(CWRK(1:L) .NE. UnitType(M)(1:UnitTypeLen(M)) )THEN
C        Didn't match anything in UnitType!
         WRITE(CTMP,'(A,A,A1)') 
     &      ' Can not recognize units: ',CWRK(1:L), LLT
         CALL LogIt(LOGCDE,1,0)
         CALL ErrMess(1104, CHNKST, LOGCDE)
         WRITE(CTMP,'(A)') RepeatLines(I)(1:79)//LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C                    
  312 UnitTypeIndex = M
C      WRITE(*,'(A,A,A1)') ' Recognized units ', 
C     &      UnitType(UnitTypeIndex)(1:UnitTypeLen(UnitTypeIndex)), LLT
C
C     Get beginning and end of second field in FS and FE:
      IF( .NOT. FIELDBND(RepeatLines(I), 256, J, FS, FE, M ) ) THEN 
         CALL ErrMess(1105, CHNKST, LOGCDE)
         WRITE(CTMP,'(A)') RepeatLines(I)(1:79)//LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF 
C     Save the beginning of the next field, if any, in J:
      J = M + 1
C
      K = ICHAR( RepeatLines(I)(FS:FS) )
      IF( K .GE. 97 .AND. K .LE. 122) K = K-32
      DataInterval = CHAR(K)
C
      IF(DataInterval .EQ. 'D' )THEN
         CALL ErrMess(1142, CHNKST, LOGCDE)
         WRITE(CTMP,'(A)') RepeatLines(I)(1:79)//LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ELSE IF(DataInterval .EQ. 'W' )THEN
         CALL ErrMess(1151, CHNKST, LOGCDE)
         WRITE(CTMP,'(A)') RepeatLines(I)(1:79)//LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ELSE IF(DataInterval .EQ. 'Q' )THEN
         TotalValues = 48
      ELSE IF(DataInterval .EQ. 'M' )THEN
         TotalValues = 12
      ELSE IF(DataInterval .EQ. 'C' )THEN
         TotalValues = 1
C        Advance to next field
         IF(.NOT. FIELDBND(RepeatLines(I), 256,J, FS, FE, M))GOTO 524
C        Read the number, and if OK, skip to the end of routine
         IF(GETFLOAT( RepeatLines(I)(FS:FE), FE-FS+1, 0, 
     &         FloatArray(1))) GOTO 230
         GOTO 522
      ELSE 
         WRITE(CTMP,'(A,A,A1)') ' Can not recognize interval: ',
     &      RepeatLines(I)(FS:FE), LLT
         CALL LogIt(LOGCDE,1,0)
         CALL ErrMess(1106, CHNKST, LOGCDE)
         WRITE(CTMP,'(A)') RepeatLines(I)(1:79)//LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
C     N will be the number of values stored in FloatArray
      N = 0
C
C     Now loop through remaining fields on first record, then
C     subsequent records, until have TotalValues read into FloatArray.
C     The 233 loop will cycle through the entries in RepeatLines.
      DO 233, K = 1, NumLines
C
C        Enter with current position on first line saved in J.
C        Re-read first entry in RepeatLines for simpler looping,
C        but ignore everything on line prior to J.  Get the index
C        to RepeatLines for the first line again in I:
         READ(FileNames(TSDataIndex)(16+K*3:), 1333) I
C
         IF(TSDVerbosity .GE. 5) THEN
C           Write(*, '(A, I8)')  '  001', I
            WRITE(CTMP,'(A,I2,A,A,A1)') ' RepeatLines(',I,
     &         '): ', RepeatLines(I)(1:235), LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
C
C           Write(*, '(A)')  '  002'
C        Get field bounds of expected value (J=1 on 2nd time through)
         IF( .NOT. FIELDBND(RepeatLines(I), 256, J, FS, FE, M ) ) THEN 
            CALL ErrMess(1107, CHNKST, LOGCDE)
            WRITE(CTMP,'(A)') RepeatLines(I)(1:79)//LLT
            CALL LogIt(LOGCDE,0,0)
            RETURN
         ENDIF 
C        Save the beginning of the next field, if any.
         J = M + 1
C
C        Read the number of values on line into L
         READ(RepeatLines(I)(FS:FE), 1333, ERR=521) L
C        Do some bounds checking
         IF( L .LT. 1 ) GOTO 521
         IF( L + N .GT. TotalValues) GOTO 523
C
         IF(TSDVerbosity .GE. 9) THEN
            WRITE(CTMP,'(A,I3,A,A,A,A1)') ' Reading ', L, 
     &         ' values from string ', RepeatLines(I)(J:220),LLT
            CALL LogIt(LOGD,0,0)
         ENDIF
         DO 243, H = 1, L
C
            IF(.NOT. FIELDBND(RepeatLines(I), 256,J, FS, FE, M))GOTO 524
C           Save the beginning of the next field, if any.
            J = M + 1
C
C           Increment the number of values stored
            N = N + 1  
C
  243       IF(.NOT. GETFLOAT( RepeatLines(I)(FS:FE), FE-FS+1, 0, 
     &         FloatArray(N))) GOTO 522
C           WRITE(*,'(A)') ' ##'//RepeatLines(I)(FS:FE)//'##'
C
C
C     Reset J to first character (of next entry in RepeatLines)
  233 J = 1
C
      IF( N .NE. TotalValues) GOTO 523
C
  230 InitRepeatData = .TRUE.
C
      IF(TSDVerbosity .GE. 3) THEN
         WRITE(CTMP,'(A,A,I2,A,A1)')
     &      ' Succesfully initialized repeating data for TimeSeries',
     &      ' number ', TSDataIndex, ' . . .', LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
      RETURN
C
  521 WRITE(CTMP,'(A,A,A1)') ' Non-integer field:', 
     &   RepeatLines(I)(FS:FS+2), LLT
      CALL LogIt(LOGCDE,1,0)
      CALL ErrMess(1108, CHNKST, LOGCDE)
      WRITE(CTMP,'(A)') RepeatLines(I)(1:79)//LLT
      CALL LogIt(LOGCDE,0,0)
      RETURN
C
  522 WRITE(CTMP,'(A,A,A1)') ' Non-floating point field:', 
     &   RepeatLines(I)(FS:FE), LLT
      CALL LogIt(LOGCDE,1,0)
      CALL ErrMess(1110, CHNKST, LOGCDE)
      WRITE(CTMP,'(A)') RepeatLines(I)(1:79)//LLT
      CALL LogIt(LOGCDE,0,0)
      RETURN
C
  523 WRITE(CTMP,'(A,I3,A,I3,A,A1)') ' Trying to use ', L + N,
     &   ' values, but exactly ',TotalValues,' are required.', LLT
      CALL LogIt(LOGCDE,1,0)
      CALL ErrMess(1109, CHNKST, LOGCDE)
      WRITE(CTMP,'(A)') RepeatLines(I)(1:79)//LLT
      CALL LogIt(LOGCDE,0,0)
      RETURN
C
  524 WRITE(CTMP,'(A,I3,A,I3,A,A1)') ' Trouble reading field ', H,
     &   ' out of ',L,' from the input line listed below.', LLT
      CALL LogIt(LOGCDE,1,0)
      CALL ErrMess(1109, CHNKST, LOGCDE)
      WRITE(CTMP,'(A)') RepeatLines(I)(1:79)//LLT
      CALL LogIt(LOGCDE,0,0)
      END
C
C----&------------------------------------------------------------------
C
C
C----&------------------------------------------------------------------
C
      LOGICAL FUNCTION InitDataFile (FileHandle, TSDataIndex, BeginDate, 
     &  FirstDataLine, DataInterval, UnitTypeIndex)
C
C     This subroutine reads lines from a data file until it recognizes
C     the units and time interval of the data, and positions the 
C     file at the first line of data used in the run.  
C
C     Data files have a much more rigid format than parameter files.
C     In addition to a strict naming convention, the contents must conform 
C     to a standard structure:  
C        1)  All lines beginning with '#' are comments; they can appear
C            anywhere in the file.
C        2)  All lines beginning with a digit ('0' - '9') contain time
C            series data values.  The format depends on the interval:
C               - monthly is I5,12F8
C               - quarter-monthly is I5,I2,12F8
C               - weekly is I5,I3,I3 1F8    (Y,M,D, 1 value)
C               - daily is I5,I3,I2,8F8     (Y,M,QM, 7-8 values)
C        3)  A line containing the string "units:" as the first 
C            non-whitepace characters must appear somewhere
C            in the file.  Once of the recognized unit indentifiers
C            (ie. cfs, m3s, m, etc) must follow the ':'
C        4)  A line containing the string "interval:" as the first 
C            non-whitepace characters must appear somewhere
C            in the file.  The next non-whitespace character after the ':'
C            must be a 'q', 'm', 'w', or 'd'.  This interval designation
C            determines the format used to read the data lines.
C     
      IMPLICIT NONE
C
C
C     Declarations regarding subroutine arguments:
C      
      INTEGER      FileHandle
C                  (Input) The unit specifier for the file to initialize.
C
      INTEGER      TSDataIndex
C                  (Input) File's element number in array FileNames.
C
      INTEGER      BeginDate(3)
C                  Input. See description in FORMATINPUT()
C
      INTEGER      FirstDataLine
C                  (Output) Saves the location of the first line of actual 
C                  time-series data, in case the units and interval 
C                  keywords appear later in file, and have to backspace 
C                  after reading them.  Returned less one, so that line
C                  number reporting in SumInData is accurate.
C
      CHARACTER    DataInterval
C                  (Output) A letter designating the data's time interval.
C                  One of 'q', 'm', 'w', or 'd'.
C
      INTEGER      UnitTypeIndex
C                  (Output) Element in array UnitType() of file's units.
C
C
C
C      
C     Declarations regarding local variables or functions:
C
      INTEGER      I, J, K, L, M, QTR 
C                  These are scratch variables.
C
      CHARACTER    ErrTmp*184
C                  Buffer to hold an error message
C
      INTEGER      TSDataLen
C                  Set to the length of the data file name.  Used for
C                  messages.
C
      CHARACTER    INLINE*128
C                  This is a scratch variable used for hold lines from the
C                  master input file and pass them to MIFLINE().
C
      INTEGER      FirstYr
      CHARACTER    FirstDateStr*6  
C                  The stores the year given on FIRSTDATALINE, and the next
C                  six characters after it.  Used to determine earlier data
C                  in file, and to see if subsequent lines are continous and
C                  sequential.
C
      CHARACTER    CHNKST*16
C                  This is the calling subroutine name passed to ErrMess() in
C                  case of error messages.
C
      INTEGER      FIRSTKC, LASTKC, DateSequence
C                  Functions that return the first and last non-whitespace 
C                  keyboard characters.
C 
      LOGICAL      READQMONTHLY, READMONTHLY, ReadDaily, ReadWeekly
C                  Functions that read lines from datafile in the appropriate 
C                  format.
C
C     INTEGER      ENDSTR
C                  Function that returns the end position of the current word.
C 
      INTEGER      GETQMONNUM 
C                  Determines what quarter-month a date is in
C
      CHARACTER*16 MakeIntervalName
C                  Returns a string with the interval name.  For output.
C 
C
C     Declarations regarding COMMON block variables:
C
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                  See description in Block Data Genral
C
      CHARACTER*256 FileNames(100)
C                  See description in Block Data Genral
C
      CHARACTER*256 RepeatLines(80)
C                  See description in Block Data Genral
C
      CHARACTER    CTMP*255,CWRK*80
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
C                  See descriptions in Block Data Genral
C
C      LOGICAL      LakesToSolve(8)
C                  See description in Block Data Genral
C
C      INTEGER      MonthLength
C                  Declared only because in common block along with run start
C                  and end dates.  Description in Block Data Genral.
C
C      INTEGER      IntervalNameLen(4)
C                  Declared only because in common block along with run start
C                  and end dates.  Description in Block Data Genral.
C
C      INTEGER      StartYear, StartMonth, StartDay, 
C     &             EndYear, EndMonth, EndDay
C                  Starting and ending year, month, day.
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
      DOUBLE PRECISION UnitConvFactor(20)
C                  See description in Block Data Genral
C
      DOUBLE PRECISION FWRK
C                  This is a common block used for 64-bit calculations.
      DIMENSION FWRK(54, 50)
C                  Dimensioned it this way to hold return values for 
C                  qmonthyread(), etc.
C
      INTEGER       TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &              MissingVal, SecsInUniformMon
      DOUBLE PRECISION MPFT, HALF, ONE, ZERO
      CHARACTER*3  MonthAbbrev(12)
C                  See descriptions in Block Data Genral
C
C
C
      COMMON /FSPACE/   FWRK
      COMMON /MONAME/   MonthAbbrev
C      COMMON /RCRUNI/   LakesToSolve, 
C     &                  IntervalNameLen, MonthLength,
C     &                  StartYear, StartMonth, StartDay,
C     &                  EndYear, EndMonth, EndDay
      COMMON /WKCHAR/   CTMP,LLT,CWRK
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /FNAMEC/   FileNames, RepeatLines
      COMMON /UNITSC/   UnitType 
      COMMON /UNITSI/   UnitTypeLen, BaseUnitType, UnitConvFactor
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
C
C
C
      CHNKST = 'InitDataFile'//LLT
C 
C     
C     Initialize return variables to bogus values 
      INITDATAFILE  = .FALSE.
      UnitTypeIndex =  0
      DataInterval  = '?'
C
C     Save number of days in starting quarter-month for later use
      QTR = GetQMonNum (BeginDate(3), BeginDate(2), BeginDate(1))
C      QTR = GetQMonNum (StartYear, StartMonth, StartDay)
C
C     Assume TSDataIndex is a valid index to FileNames(), or else wouldn't
C     be running this routine.  Determine length of filename.
      TSDataLen = INDEX(FileNames(TSDataIndex), LLT) - 1
C
      OPEN(UNIT=FileHandle,FILE=FileNames(TSDataIndex)(1:TSDataLen),
     &   STATUS='OLD', IOSTAT=J) 
      IF ( J .NE. 0 ) THEN
         CALL ErrMess(1003,CHNKST,LOGCDE)
         WRITE(CTMP,9015) J,FileNames(TSDataIndex)(1:TSDataLen),LLT
         CALL LogIt(LOGCDE,1,0)
         RETURN
      ENDIF
      IF(TSDVerbosity .GE. 7) THEN
         WRITE(CTMP,9761) ' Successfully opened file ',
     +      FileNames(TSDataIndex)(1:TSDataLen), '  . . .'//LLT
         CALL LogIt(LOGD,1,0)
C        IF(TSDVerbosity .GE. 7) THEN
C            CALL LogIt(LOGC,1,0)
C        ENDIF
      ENDIF
C           
      FIRSTDATALINE =  0 
C
C 
C     Initiate a loop to read in lines until know the first data 
C     line, units, and interval 
      DO 1231 K= 1, MaxInt
C
         READ(FileHandle, '(A)', END=1234, ERR=1233, IOSTAT=J) INLINE
C
C        Must be no error, so don't need J anymore, use it to
C        store first character:
         J = FIRSTKC(INLINE, 127)
         IF( J .NE. 0 .AND. INLINE(J:J) .NE. '#' ) THEN
C           
C           Line has data, use first char to see if it's what we want
            I = ICHAR( INLINE(J:J) )
C             
            IF(I .GE. 48 .AND. I .LE. 57) THEN
C              It's a number!  Must be a data line, set FIRSTDATALINE 
               IF(FIRSTDATALINE .EQ. 0 ) THEN
                  FIRSTDATALINE = K
                  READ(INLINE(1:11),'(I5,A6)')FirstYr, FirstDateStr
                  IF( FirstYr .LE. 0) THEN
                     WRITE(ErrTmp,'(A,A,A1)') 
     &                   ' Invalid year: ', INLINE(1:5), LLT
                     I = 1077
                     GOTO 975
                  ENDIF
               ENDIF
C               
               IF(DataInterval .EQ. 'Q') THEN
                  READ(INLINE(6:7), '(I2)') L
                  IF( L .LT. 1 .OR. L .GT. 4)THEN
                     WRITE(ErrTmp,'(A,I5,A1)') 
     &                   ' Invalid quarter in column 7: ', L, LLT
                     I = 1076
                     GOTO 975
                  ENDIF
               ENDIF
               IF(DataInterval .EQ. 'W') THEN
                  READ(INLINE(6:11), '(I3,I3)') I, L
                  IF( L .LT. 1 .OR. L .GT. 31 .OR. 
     &               I .LT. 1 .OR. I .GT. 12 )THEN
                     WRITE(ErrTmp,'(A,2I5,A1)') 
     &                   ' Invalid month or date in columns 6-11: ', 
     &                   I, L, LLT
                     I = 1153
                     GOTO 975
                  ENDIF
               ENDIF
               IF(DataInterval .EQ. 'D') THEN
                  READ(INLINE(7:10), '(I2,I2)') I, L
                  IF( L .LT. 1 .OR. L .GT. 4 .OR. 
     &               I .LT. 1 .OR. I .GT. 12 )THEN
                     WRITE(ErrTmp,'(A,2I5,A1)') 
     &                   ' Invalid month or quarter in columns 7-10: ', 
     &                   I, L, LLT
                     I = 1152
                     GOTO 975
                  ENDIF
               ENDIF
C                
            ELSE IF( I .EQ. 85 .OR. I .EQ. 117 .OR. 
     &         I .EQ. 73 .OR. I .EQ. 105) THEN
C              It's a 'U' or 'I'.  Must process line for units or interval
C              Find last character on line:
               L = 0
               DO 40, M = 1, LASTKC(INLINE, 127) - J + 1
                  I = ICHAR( INLINE(J+M-1:J+M-1) )
                  IF( I .GE. 97 .AND. I .LE. 122) THEN
                     I = I-32
                  ENDIF
                  IF( (I .GE. 48 .AND. I .LE. 57) .OR.
     &               (I .GE. 65 .AND. I .LE. 90)) THEN
                     L = L + 1
                     IF( L .GT. 80) STOP 'Bullwinkle!'
                     CWRK(L:L) = CHAR(I)
                  ENDIF
C               WRITE(CTMP,'(A,I2,A,A,A,I3,A,A,A1,A1)')'   char ', M,
C     &            ': ',INLINE(J+M-1:J+M-1),' (',I,') CWRK={', CWRK(1:L), 
C     &            '}', LLT
C               CALL LogIt(LOGCD,0,0)
 40            CONTINUE
C               WRITE(CTMP,'(A,A,A1,A1)')' CWRK={', CWRK(1:L), '}', LLT
C               CALL LogIt(LOGCD,0,0)
               IF( CWRK(1:8) .EQ. 'INTERVAL' ) THEN
                  IF(CWRK(9:9) .EQ. 'D' .OR. CWRK(9:9) .EQ. 'W' .OR.
     &               CWRK(9:9) .EQ. 'Q' .OR. CWRK(9:9) .EQ. 'M' ) THEN
                     DataInterval=CWRK(9:9)
                  ELSE
                     WRITE(ErrTmp,'(A,A,A1)') 
     &                   ' Can not recognize interval ', CWRK(9:9), LLT
                     I = 1069
                     GOTO 975
                  ENDIF
               ELSE IF( CWRK(1:5) .EQ. 'UNITS' ) THEN
C                
C                 Now search through the list of units (except #20)
                  DO 311, I=1,19
                    IF(CWRK(6:L) .EQ. UnitType(I)(1:UnitTypeLen(I)))THEN
C                        WRITE(*,'(6X,A,A,A)')CWRK(6:L),' matches ',
C     &                     UnitType(I)(1:UnitTypeLen(I))
                        GOTO 312
C                     ELSE
C                        WRITE(*,'(6X,A,A,A)')CWRK(6:L),
C     &                     ' does not match ',
C     &                     UnitType(I)(1:UnitTypeLen(I))
                     ENDIF
  311             CONTINUE
C
C                 IF( I .EQ. 16 .AND. CWRK(1:L) .NE. 
C    &               UnitType(I)(1:UnitTypeLen(I)) )THEN
                  IF(CWRK(1:L) .NE. UnitType(I)(1:UnitTypeLen(I)) )THEN
C                    Didn't match anything in UnitType!
                     WRITE(ErrTmp,'(A,A,A1)') 
     &                   ' Can not recognize units ', CWRK(1:L), LLT
                     I = 1041
                     GOTO 975
                  ENDIF
C                    
C                 Good.  UnitType(I) matches the units.
  312             UnitTypeIndex = I
               ELSE
C                 Hmmm, started with 'U' or 'I' but didn't match UNITS or INTERVAL        
                  WRITE(ErrTmp,'(A,A)') ' Mis-typed keyword: ',
     &               CWRK(1:L)//LLT
                  I = 1070
                  GOTO 975
               ENDIF
            ELSE
C              Hmmm, not blank, and didn't start with '#', 'U', or 'I'
               WRITE(ErrTmp,'(A,A)') 
     &         ' Expected #, U, or I as first character, but received '
     &         , INLINE(J:J)
               I = 1071
               GOTO 975
            ENDIF
         ENDIF
C      
 1231 IF( UnitTypeIndex .GT. 0 .AND. FIRSTDATALINE .GT. 0 .AND.
     &   DataInterval .NE. '?') GOTO 1236
C
C     Jump here if error reading a line in 
 1233 WRITE(ErrTmp,'(A,I5,A1)')'IOSTAT=',J,LLT
      I = 1003
      GOTO 975
C
C     Read all the lines, but was missing at least one item!
 1234 IF( UnitTypeIndex .EQ. 0 )THEN
         I = 1072
      ELSE IF( FIRSTDATALINE .EQ. 0 )THEN
         I = 1073
      ELSE IF( DataInterval .EQ. '?') THEN
         I = 1074
      ELSE
         STOP 'How the heck did this happen?'
      ENDIF
      CALL ErrMess(I,CHNKST,LOGCDE)
      WRITE(CTMP,'(3A)')' File:',
     &   FileNames(TSDataIndex)(1:TSDataLen), LLT
      CALL LogIt(LOGCDE,0,0)
      RETURN
C
 1235 WRITE(ErrTmp,'(A,I5,A,I5,A,I5,A,I5,A1)')
     &   ' Backspacing from line ', K, ' to line ', 
     &   FIRSTDATALINE, ', got IOSTAT=', L, ' at line ', J, LLT
      I = 1075
      GOTO 975
C
 1237 I = 1079
      GOTO 975
C
  975 CONTINUE
      CALL ErrMess(I,CHNKST,LOGCDE)
      IF( I .LT. 1072 .OR. I .GT. 1074) THEN      
         CTMP = ErrTmp
         CALL LogIt(LOGCDE,0,0)
         WRITE(CTMP,'(A,I5,A1,A1)')' Problem with line ',K,
     &      ':',LLT
         CALL LogIt(LOGCDE,1,0)
         WRITE(CTMP,'(A,A1)')INLINE(1:LASTKC(INLINE, 127)),LLT
         CALL LogIt(LOGCDE,0,0)
      ENDIF
      RETURN
C
C     Well, looks like eveything, is known, back file pointer up
C     to just before the first data line:
 1236 DO 1239, J = K, FIRSTDATALINE, -1
 1239    BACKSPACE(FileHandle, ERR=1235, IOSTAT=L)
C
C
C     See if starting date precedes data on first data line
C     Use I to indicate file doesn't have old enough data
      I = 0
C
      IF(DataInterval .EQ. 'M' .OR. DataInterval .EQ. 'Q') THEN
         L = 1
         IF( BeginDate(3) .LT. FirstYr) I = 1096
      ELSE IF (DataInterval .EQ. 'W' ) THEN
         L = 6
         READ(FirstDateStr, '(I3,I3)') K, J
C         IF( (DateSequence(J,K,FirstYr) - 6) .GT. DateSequence(StartDay, 
C     &      StartMonth, StartYear) ) I = 1096
         IF( (DateSequence(J,K,FirstYr) - 6) .GT. DateSequence(
     &      BeginDate(1), BeginDate(2), BeginDate(3)) ) I = 1096
      ELSE IF (DataInterval .EQ. 'D' ) THEN
         L = 5
         READ(FirstDateStr, '(I3,I2)') K, J
         IF( FirstYr .GT. BeginDate(3) .OR.
     &      (FirstYr .EQ. BeginDate(3) .AND. K .GT. BeginDate(2)) .OR.
     &      (FirstYr .EQ. BeginDate(3) .AND. K .EQ. BeginDate(2) .AND.
     &      J .GT. GetQMonNum (BeginDate(3), BeginDate(2), BeginDate(1))
     &      ))I=1096
      ELSE
         STOP 'Stimpy'
      ENDIF
      IF(I .NE. 0) THEN
         CALL ErrMess(I,CHNKST,LOGCDE)
         WRITE(CTMP,'(A,A,A1)')
     &      ' Start date for run precedes data in file ',
     &               FileNames(TSDataIndex)(1:TSDataLen), LLT
         CALL LogIt(LOGCDE,0,0)
         WRITE(CTMP,'(A,I2,A3,I5,A,I4,A,A1)') ' Run starts ',
     &      BeginDate(1), MonthAbbrev(BeginDate(2)), BeginDate(3), 
     &      ', but first data for: ', FirstYr, FirstDateStr(1:L), LLT
         CALL LogIt(LOGCDE,0,0)
         RETURN
      ENDIF
C
C
C     Now advance the file pointer until get to line with data 
C     for start of run:
C
      IF( DataInterval .EQ. 'M') THEN
         I = FirstYr - 1
         DO 1330, M = 1, MaxInt
            I = I + 1
 6330       READ(FileHandle, '(A)', ERR=1237) INLINE
            IF( INLINE(1:1) .EQ. '#' ) GOTO 6330
            IF(.NOT. READMONTHLY(INLINE, J, FWRK))THEN
               K = M + FIRSTDATALINE - 1
               WRITE(ErrTmp,'(A,I5,A1)')
     &            ' READMONTHLY() failed, IOSTAT=',J,LLT
               I = 1077
               GOTO 975
            ENDIF
            IF( J .NE. I) THEN
               WRITE(ErrTmp,'(A,I5,A,I5,A)')
     &            ' Expected line for year ',I,', but got ',J,
     &            ' instead!'//LLT
               I = 1078
               K = M + FIRSTDATALINE - 1
               GOTO 975
            ENDIF
C            
C           Jump out of loop when get to first data line
            IF( BeginDate(3) .EQ. J) THEN
C              Save routing interval number in J
               J = 4
               GOTO 1390
            ENDIF
 1330       CONTINUE
C             
      ELSE IF( DataInterval .EQ. 'Q') THEN
         DO 1333, M = 1, MaxInt
 6333       READ(FileHandle, '(A)', ERR=1237) INLINE
            IF( INLINE(1:1) .EQ. '#' ) GOTO 6333
            IF(.NOT. READQMONTHLY(INLINE, J,L, FWRK))THEN
               WRITE(ErrTmp,'(A,I5,A, I5,A1)')
     &            ' READQMONTHLY() returned YEAR=',J,' Qtr=',L,LLT
               I = 1076
               K = M + FIRSTDATALINE - 1
               GOTO 975
            ENDIF
            IF( BeginDate(3) .EQ. J) THEN
               IF( L .EQ. 1 ) THEN
                  J = 3
                  GOTO 1390
C              Break (Ok, jump) out of loop when get to first data line
C               IF( (StartMonth + 2)/3 .EQ. L) THEN
C                 Save routing interval number in J
C                  J = 3
C                  GOTO 1390
C               ENDIF
               ELSE
                  WRITE(ErrTmp,'(A,I5,A, I5,A1)')
     &               ' READQMONTHLY() returned YEAR=',J,' Qtr=',L,LLT
                  I = 1076
                  K = M + FIRSTDATALINE - 1
                  GOTO 975
               ENDIF
            ENDIF
 1333       CONTINUE
C         
 1432 FORMAT(' ReadDaily() returned year=', I5, ' month=', I2,
     &            ' qtr-month=', I1, A1)     
      ELSE IF( DataInterval .EQ. 'D') THEN
         DO 1433, M = 1, MaxInt
 6433       READ(FileHandle, '(A)', ERR=1237) INLINE
            IF( INLINE(1:1) .EQ. '#' ) GOTO 6433
            IF(.NOT. ReadDaily(INLINE, J, K, L, FWRK))THEN
               WRITE(ErrTmp,1432) J, K, L, LLT
               I = 1143
               K = M + FIRSTDATALINE - 1
               GOTO 975
            ENDIF
            IF(BeginDate(3) .EQ. J .AND. BeginDate(2) .EQ. K 
     &         .AND. QTR .EQ. L)THEN
               J = 1
               GOTO 1390
            ENDIF
 1433       CONTINUE
C         
 1532 FORMAT(' ReadWeekly() returned year=', I5, ' month=', I2,
     &            ' date=', I2, A1)     
      ELSE IF( DataInterval .EQ. 'W') THEN
         DO 1533, M = 1, MaxInt
 6434       READ(FileHandle, '(A)', ERR=1237) INLINE
            IF( INLINE(1:1) .EQ. '#' ) GOTO 6434
            IF(.NOT. ReadWeekly(INLINE, J, K, L, FWRK))THEN
               WRITE(ErrTmp,1532) J, K, L, LLT
               I = 1154
               K = M + FIRSTDATALINE - 1
               GOTO 975
            ENDIF
C
C           Jump out when date read in on the line is later than 
C           run start date.
            IF( J .GT. BeginDate(3) .OR.
     &         (J .EQ. BeginDate(3) .AND. K .GT. BeginDate(2)) .OR.
     &         (J .EQ. BeginDate(3) .AND. K .EQ. BeginDate(2) .AND.
     &         L .GE. BeginDate(1)) ) THEN
               J = 2
               GOTO 1390
            ENDIF
 1533       CONTINUE
      ENDIF
C
      STOP 'Data unavailable at start date'
C
C     
 1390 BACKSPACE(FileHandle)
C
C     
      IF(TSDVerbosity .GE. 2) THEN
         CWRK = MakeIntervalName(J,.FALSE.,I)
C
         WRITE(CTMP,'(A,A,A,A,A,I2,A)') '  . . . containing ', 
     +      CWRK(1:I), ' data in ',
     +      UnitType(UnitTypeIndex)(1:UnitTypeLen(UnitTypeIndex)), 
     +      ', beginning on line ', FIRSTDATALINE, '.'//LLT 
         CALL LogIt(LOGD,0,0)
         IF(TSDVerbosity .GE. 9) THEN
            WRITE(CTMP,'(A,A1)')
     +         ' First line of data used in simulation:',LLT
            CALL LogIt(LOGD,0,0)
            WRITE(CTMP,'(A,A1)')INLINE(1:LASTKC(INLINE, 127)),LLT
            CALL LogIt(LOGD,0,0)
C           WRITE(CTMP,'(A,12D18.11,A1)') '  floats:', FWRK,LLT
C           CALL LogIt(LOGD,0,0)
         ENDIF
      ENDIF
C
      INITDATAFILE  = .TRUE.
      FirstDataLine = FirstDataLine - 1
      IF(TSDVerbosity .GE. 2) THEN
         WRITE(CTMP,9761) ' Successfully initialized file ',
     +      FileNames(TSDataIndex)(1:TSDataLen), '  . . .'//LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C      
 9015 FORMAT(' Got IOSTAT=',I5,' trying to open file ',A,'.',A1)
 9761 FORMAT(A,A,A)
C           
      END
C
C----&------------------------------------------------------------------
C----&------------------------------------------------------------------
C
      SUBROUTINE FORMATINPUT()
C
C
C     Directs subroutine SUMINDATA() in reading the various time-series 
C     data files and summarizing input files into temporary files that 
C     the regulation and routing computations use.
C
C     FormatInput() checks to see if the input
C     summary file was specified as a parameter, or must be created.
C     If it needs to be generated, FormatInput() assigns appropriate
C     values to the arrays DFNUM and DFTXT.  
C
C     DFNUM lists the indexes to array FileNames (in COMMON block
C     FNAMEC) that correspond to the files needed to summarize the 
C     time-series input data for a given lake.  The files fall into
C     three groups.  The last files listed in DFNUM are the 
C     inflow/outflow files.  The middle group consists of ice/weed 
C     retardation files.  The first files listed in DFNUM are any other 
C     oddball files applicable to a particular lake or its outlet
C     river.  The list in DFNUM is terminated by an element with a
C     value of 0.  For instance, for Lake Erie, the only
C     oddball file is Chippewa-Grass Island Pool levels.  There are
C     three retardation files permitted: one for ice, one for weed, 
C     and another for combined retardation.  (Note that these data
C     must be in flow units, as is traditional practice.  Never really
C     understood this.  Seems to me that the amount of retardation
C     would be a function of the upstream and downstream levels, 
C     rather than some absolute flow rate, but what do i know.  David
C     says he uses some other mechanism on the St. Lawrence, like a 
C     percentage or something.  Such input for a lake would be 
C     accomodated as another oddball time series.)  There six inflow/
C     outflow files available for Lake Erie: net basin supplies, 
C     Detroit River inflow, Welland Canal outflow, consumptive use,
C     "other" diversion, and "other" supplies.  DFNUM(11) gets assigned
C     a value of zero to terminate the list.
C      
C     SumInData() contains an array named DFHEADING, which has the list
C     of possible column headings to use to describe each input file 
C     (such as 'NBS', 'Welland', or 'CGIP Lv').  After storing the 
C     appropriate index values for each file data file in DFNUM, 
C     FormatInput() stores another index value in array DFTXT, which 
C     simply refers to which element in DFHEADING to use.  
C
C     After preparing the arrays DFNUM and DFTXT, FormatInput() calls
C     SumInData(), using the arrays as arguments.  The next two 
C     arguments specify how many of the files in DFNUM are oddballs  
C     and retards.  The fifth argument to SumInData() indicates which
C     lake is being summarized (actually, SumInData() doesn't care
C     very much, using LNUM only as an index to array LAKENAME for 
C     messages and such).  The sixth argument is the index to array 
C     FileNames (in COMMON block FNAMEC), for the file name to write
C     the inflow/outflow summary to.  The final argument specifies
C     which routing interval (daily, weekly, quarter-monthly, or 
C     monthly) to use for generating the summary of time-series data.
C     Frinstance, if the second element of array RoutingInterval (from
C     COMMON block MISINT) is 3, then the middle lakes routing requires
C     quarter-monthly routing summary files for Lakes Michigan-Huron,
C     St. Clair, and Erie.  If Plan 77A is used, FormatInput() calls 
C     SumInData() again to generate a monthly summary (with potentially
C     different input files?).
C
C  The following paragraphs blather on about the complications of 
C  expressing evaporation and precipitation type data in cfs when
C  doing variable surface-area runs.  Proceed ONLY if you specialize
C  in writing research papers on arcane academic issues with no
C  practical value.
C           
C  This is going to be a problem for variable surface area.  Part of the
C  subsequent code tries to sum up the precip, evap, and runoff components,
C  which are *apparently* read in terms of inches of water surface area.  
C  Then they get converted to cfs (or m3s) by multiplying by the water 
C  surface area and dividing by however many seconds there are in the month.
C  Then the NBS's get used in the water balance equations.  The problem 
C  arises when surface area is not constant.  Then you can't convert the
C  inches of supply read in to m3s until you know the surface area, which
C  means you have to know the water level, which is what you are trying to 
C  figure out.  ;-)   So without making some simplifying assumptions, in
C  order to use rainfall or evaporation input expressed in terms of inches
C  on water surface (the best way to express them), it looks like we have
C  add yet another unknown and equation (or three of them) to the existing
C  simultaneous solution.
C
C  Personally, I like simplifying solutions.  So there are some ideas to
C  discuss.  First (and easiest) way of dealing with this problem, is to
C  assume constant area for the purpose of converting inches to cfs.  So
C  for expediency until the issue gets thrown back into GLERL's lap, I 
C  just used whatever level in MHZ, SCZ, or ERZ was left over from the loop
C  through the previous year to compute a constant area for the upcoming
C  year.  I'll leave it up to the big thinkers to figure out if this 
C  introduces some type of systematic bias (like does using last year's
C  lower level under-estimate the impacts of a lot of rain this year, 
C  artificially dampening a high level, or in reverse, helping back-stop
C  a low level?)
C
C  So the approaches seem to include:
C
C     1)  Read all the components in terms of cfs (or m3s), avoiding the
C         problem of converting inches to cfs in the model.  However, this
C         implicitly assumes a constant area, which only moves the 
C         inaccuracy around (the same amount of rainfall actually is more
C         cfs when the lake is higher).
C
C     2)  Read the runoff component data in as cfs, since the water
C         surface area isn't really intrinsic to that data type.
C         Read the evap and precip data in inches, but convert them
C         to cfs during the routing, approximating the mean lake level
C         with the beginning-of-period level.  This doesn't look too
C         onerous, and is probably worth doing if actually start doing
C         variable surface area runs.  Will have to provide another
C         category of time-series files to add up when doing the input
C         summaries; evap and precip inches will get totaled up, with
C         the sum for each timestep listed along with the retardation
C         and inflow/outflow sums.
C
C     3)  Re-derive water balance equations to include simultaneous
C         solution for water surface area (shudder).  'Nuff said.
C
C  Of course, if the input parameters call for constant area, or not using
C  component nbs, then this all becomes moot.  Also, does it really matter?
C
      IMPLICIT NONE
C
C     Declarations regarding local variables or functions:
C
      INTEGER      I
C                  These are scratch variables.
C
      INTEGER      DFNUM(48)
C                  Array passed to SUMINPUT() containing which index numbers
C                  to FileNames() should be summarized.
C
      INTEGER      DFTXT(48)
C                  Array passed to SUMINPUT() containing which index numbers
C                  to DFHEADINGS() correspond to DFNUM.  Used for labeling the
C                  columns in the summary supply files.
C
      INTEGER      DFSIGN(48)
C                  Array passed to SUMINPUT() containing sign to use for flow
C                  direction.  +1 means flow into lake.  -1 means flow out 
C                  of lake.  These are actually implemented as multipliers
C                  to the data values.  
C
      INTEGER      BeginDate(3)
C                  Array passed to SUMINPUT() containing day, month, and year
C                  that times series should be begin.  Unless using Plan 1958D,
C                  this should be the same as the start date parameters.
C
C
C
C     Declarations regarding COMMON block variables:
C
C     INTEGER      ENDSTR
C                  Function that returns the end position of the current word.
C 
C     INTEGER      LastKC
C                  Function that returns the position of last 'keyboard char'
C 
      CHARACTER*256 FileNames(100)
C                  See description in Block Data Genral
C                  Note: Yes, it really is "Genral", because it dates 
C                  back to a time when variables of longer than 6 
C                  characters confused FORTRAN.  Deal with it, Travis. :P
C
      CHARACTER*256 RepeatLines(80)
C                  See description in Block Data Genral
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
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
      LOGICAL      LakesToSolve(8)
      CHARACTER    CTMP*255, LLT*1, CWRK*80
      LOGICAL      LOGC(8), LOGD(8), LOGS(8), 
     &             LOGCD(8), LOGCE(8), LOGCDE(8), LOGSD(8), 
     &             LOGM(8), LOGP(8), LOGO(8), 
     &             LOGDM(8), LOGDP(8), LOGDO(8), 
     &             LOGSM(8), LOGSP(8), LOGSO(8),
     &             LOGSDP(8), LOGSDM(8), LOGSDO(8),
     &             LOGSMPO(8), LOGSDMPO(8), LOGDMPO(8)
C                  See descriptions in Block Data Genral
C
      INTEGER      StartYear, StartMonth, StartDay, 
     &             EndYear, EndMonth, EndDay,
     &             MonthLength, IntervalNameLen(4)
C                  Declared because in same COMMON as LakesToSolve(8)
      INTEGER      P58DStartDay, P58DStartMonth, P58DStartYear
C                  See Desc in CheckPar
C                  
      COMMON /FSAT/     P58DStartDay, P58DStartMonth, P58DStartYear
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /FNAMEC/   FileNames, RepeatLines
      COMMON /MISINT/   MidLakeSolution, BOMRound, FlowRound, 
     &                  LevelRound, RoutingInterval, RoutingIncrement,
     &                  StMarysSolution, UStLawrSolution, P77ACompat
      COMMON /RCRUNI/   LakesToSolve, 
     &                  IntervalNameLen, MonthLength,
     &                  StartYear, StartMonth, StartDay,
     &                  EndYear, EndMonth, EndDay
      COMMON /LOGCOD/   LOGC, LOGD, LOGS, LOGCD, LOGCE, LOGCDE, LOGSD, 
     &                  LOGM, LOGP, LOGO, LOGDM, LOGDP, LOGDO, 
     &                  LOGSM, LOGSP, LOGSO, LOGSDP, LOGSDM, LOGSDO,
     &                  LOGSMPO, LOGSDMPO, LOGDMPO
      COMMON /WKCHAR/   CTMP,LLT,CWRK
C
C
      WRITE(CTMP,'(A,A)') 
     &   ' Pre-processing time series data . . . ', LLT
      IF( GenVerbosity .GE. 2 ) THEN 
         CALL LogIt(LOGCD,1,0)
      ELSE
         CALL LogIt(LOGCD,1,1)
      ENDIF
C
 3456 FORMAT(' Generating ', A,' time-series data summary . . .', A1)
C
C     These values remain in place until get to Lake Ontario
      BeginDate(1) = StartDay
      BeginDate(2) = StartMonth
      BeginDate(3) = StartYear
C
C                     Lake Superior
      I = INDEX(FileNames(41),LLT) + 1
      IF(LakesToSolve(1) .AND. FileNames(41)(I:I+6) .EQ. 'default') THEN
C        Inflow/outflow summary file not specified, must build it
C
         IF(TSDVerbosity .GE. 2 ) THEN
            WRITE(CTMP,3456) 'Lake Superior', LLT
            CALL LogIt(LOGCD,0,0)
         ENDIF
C
C        Initialize DFNUM with the index numbers to Lakes Superior files in FileNames()
C                         FileNames(40) = SUPGATESDATA
                             DFNUM(1)=40
                             DFTXT(1)=36
                             DFSIGN(1)=1 
C                         FileNames(69) = SUPSIDECHANDATA
                             DFNUM(2)=69
                             DFTXT(2)=31
                             DFSIGN(2)=1 
C                         FileNames(68) = SUPOUTFLOWDATA
                             DFNUM(3)=68
                             DFTXT(3)=16
                             DFSIGN(3)=1 
C                         FileNames(49) = SUPDEVIATIONDATA
                             DFNUM(4)=49
                             DFTXT(4)=22
                             DFSIGN(4)=1 
C                         FileNames(70) = SUPEOPLEVELDATA
                             DFNUM(5)=70
                             DFTXT(5)=32
                             DFSIGN(5)=1 
C                         FileNames(93) = SUPINDEX1DATA
                             DFNUM(6)=93
                             DFTXT(6)=52
                             DFSIGN(6)=1 
C                         FileNames(94) = SUPINDEX2DATA
                             DFNUM(7)=94
                             DFTXT(7)=53
                             DFSIGN(7)=1 
C                         FileNames(95) = SUPINDEX3DATA
                             DFNUM(8)=95
                             DFTXT(8)=54
                             DFSIGN(8)=1
C  ------------------end of oddball files to show in first columns of tsd file-----------
C                         FileNames(50) = STMRIVRETRDATA
                             DFNUM(9)=50
                             DFTXT(9)=13
                             DFSIGN(9)=1 
C                         FileNames(51) = STMRIVICEDATA
                             DFNUM(10)=51
                             DFTXT(10)=14
                             DFSIGN(10)=1 
C                         FileNames(52) = STMRIVWEEDDATA
                             DFNUM(11)=52
                             DFTXT(11)=15
                             DFSIGN(11)=1 
C                         FileNames(42) = SUPNBSDATA
                             DFNUM(12)=42
                             DFTXT(12)=1
                             DFSIGN(12)=1 
C                         FileNames(43) = SUPCONUSEDATA
                             DFNUM(13)=43
                             DFTXT(13)=4
                             DFSIGN(13)= -1
C                         FileNames(44) = LONGLACDIVDATA
                             DFNUM(14)=44
                             DFTXT(14)=11
                             DFSIGN(14)=1 
C                         FileNames(45) = OGOKIDIVDATA 
                             DFNUM(15)=45
                             DFTXT(15)=12
                             DFSIGN(15)=1 
C                         FileNames(46) = COMBLLODIVDATA 
                             DFNUM(16)=46
                             DFTXT(16)=10
                             DFSIGN(16)=1 
C                         FileNames(47) = SUPOTHERDIVDATA
                             DFNUM(17)=47
                             DFTXT(17)=2
                             DFSIGN(17)=1 
C                         FileNames(48) = SUPOTHERSUPDATA
                             DFNUM(18)=48
                             DFTXT(18)=3
                             DFSIGN(18)=1 
C
         DFNUM(19)=0
         CALL SUMINDATA( DFNUM, DFTXT, DFSIGN, 8, 3, 1, 41, 
     &      RoutingInterval(1), BeginDate) 
      ENDIF
C
C
      I = INDEX(FileNames(31),LLT) + 1
      IF(LakesToSolve(2) .AND. FileNames(31)(I:I+6) .EQ. 'default') THEN
C        Inflow/outflow summary file not specified, must build it
C
         IF(TSDVerbosity .GE. 2 ) THEN
            WRITE(CTMP,3456) 'Lakes Michigan-Huron', LLT
            CALL LogIt(LOGCD,0,0)
         ENDIF
C
C        separate tallies for Mich and huron not implemented
C
C        Initialize DFNUM with the index numbers to Lakes Michigan-Huron files in FileNames()
C                         FileNames(71) = MHUEOPLEVELDATA
                             DFNUM(1)=71
                             DFTXT(1)=32
                             DFSIGN(1)=1 
C                         FileNames(26) = USSICEWEEDADJ
                             DFNUM(2)=26
                             DFTXT(2)=55
                             DFSIGN(2)=1 
C                         FileNames(16) = STCRIVRETRDATA
                             DFNUM(3)=16
                             DFTXT(3)=13
                             DFSIGN(3)=1 
C                         FileNames(34) = STCRIVICEDATA
                             DFNUM(4)=34
                             DFTXT(4)=14
                             DFSIGN(4)=1 
C                         FileNames(35) = STCRIVWEEDDATA
                             DFNUM(5)=35
                             DFTXT(5)=15
                             DFSIGN(5)=1 
C                         FileNames(19) = STMRIVFLOWDATA
                             DFNUM(6)=19
                             DFTXT(6)=16
                             DFSIGN(6)=1 
C                         FileNames(1)  = MHUNBSDATA
                             DFNUM(7)=1
                             DFTXT(7)=1
                             DFSIGN(7)=1 
C                         FileNames(22) = CHICAGODIVDATA
                             DFNUM(8)=22
                             DFTXT(8)=9
                             DFSIGN(8)= -1 
C                         FileNames(23) = MHUCONUSEDATA
                             DFNUM(9)=23
                             DFTXT(9)=4
                             DFSIGN(9)= -1 
C                         FileNames(24) = MHUOTHERDIVDATA
                             DFNUM(10)=24
                             DFTXT(10)=2
                             DFSIGN(10)=1 
C                         FileNames(25) = MHUOTHERSUPDATA
                             DFNUM(11)=25
                             DFTXT(11)=3
                             DFSIGN(11)=1 
C                         FileNames(4)  = MICCONUSEDATA
                             DFNUM(12)=4
                             DFTXT(12)=4
                             DFSIGN(12)= -1 
C                         FileNames(12) = MICOTHERDIVDATA
                             DFNUM(13)=12
                             DFTXT(13)=2
                             DFSIGN(13)=1 
C                         FileNames(8)  = MICOTHERSUPDATA
                             DFNUM(14)=8
                             DFTXT(14)=3
                             DFSIGN(14)=1 
C                         FileNames(5)  = HURCONUSEDATA
                             DFNUM(15)=5
                             DFTXT(15)=4
                             DFSIGN(15)= -1 
C                         FileNames(13) = HUROTHERDIVDATA
                             DFNUM(16)=13
                             DFTXT(16)=2
                             DFSIGN(16)=1 
C                         FileNames(9)  = HUROTHERSUPDATA 
                             DFNUM(17)=9
                             DFTXT(17)=3
                             DFSIGN(17)=1 
         DFNUM(18)=0
         CALL SUMINDATA( DFNUM, DFTXT, DFSIGN, 2, 3, 2, 31, 
     &      RoutingInterval(2), BeginDate) 
      ENDIF
C
C
      I = INDEX(FileNames(32),LLT) + 1
      IF(LakesToSolve(3) .AND. FileNames(32)(I:I+6) .EQ. 'default') THEN
C        Inflow/outflow summary file not specified, must build it
C
         IF(TSDVerbosity .GE. 2 ) THEN
            WRITE(CTMP,3456) 'Lake St. Clair', LLT
            CALL LogIt(LOGCD,0,0)
         ENDIF
C
C        Initialize DFNUM with the index numbers to Lake St. Clair files in FileNames()
C                         FileNames(72) = STCEOPLEVELDATA
                             DFNUM(1)=72
                             DFTXT(1)=32
                             DFSIGN(1)=1 
C                         FileNames(17) = DETRIVRETRDATA
                             DFNUM(2)=17
                             DFTXT(2)=13
                             DFSIGN(2)=1 
C                         FileNames(36) = DETRIVICEDATA  
                             DFNUM(3)=36
                             DFTXT(3)=14
                             DFSIGN(3)=1 
C                         FileNames(37) = DETRIVWEEDDATA  
                             DFNUM(4)=37
                             DFTXT(4)=15
                             DFSIGN(4)=1 
C                         FileNames(28) = STCRIVFLOWDATA  
                             DFNUM(5)=28
                             DFTXT(5)=17
                             DFSIGN(5)=1 
C                         FileNames(2)  = STCNBSDATA      
                             DFNUM(6)=2
                             DFTXT(6)=1
                             DFSIGN(6)=1 
C                         FileNames(6)  = STCCONUSEDATA   
                             DFNUM(7)=6
                             DFTXT(7)=4
                             DFSIGN(7)= -1 
C                         FileNames(14) = STCOTHERDIVDATA
                             DFNUM(8)=14
                             DFTXT(8)=2
                             DFSIGN(8)=1 
C                         FileNames(10) = STCOTHERSUPDATA
                             DFNUM(9)=10
                             DFTXT(9)=3
                             DFSIGN(9)=1 
         DFNUM(10)=0
         CALL SUMINDATA( DFNUM, DFTXT, DFSIGN, 1, 3, 5, 32, 
     &      RoutingInterval(2), BeginDate) 
      ENDIF
C
C
C
      I = INDEX(FileNames(33),LLT) + 1
      IF(LakesToSolve(4) .AND. FileNames(33)(I:I+6) .EQ. 'default') THEN
C
C        Inflow/outflow summary file not specified, must build it
C
         IF(TSDVerbosity .GE. 2 ) THEN
            WRITE(CTMP,3456) 'Lake Erie', LLT
            CALL LogIt(LOGCD,0,0)
         ENDIF
C
C        Initialize DFNUM with the index numbers to Lake Erie files in FileNames()
C                         FileNames(20) = CGIPAVGLEVELDATA
                             DFNUM(1)=20
                             DFTXT(1)=20
                             DFSIGN(1)=1 
C                         FileNames(73) = ERIEOPLEVELDATA
                             DFNUM(2)=73
                             DFTXT(2)=32
                             DFSIGN(2)=1 
C                         FileNames(21) = WELLANDCANALDATA
                             DFNUM(3)=21
                             DFTXT(3)=8
                             DFSIGN(3)= 1 
C                         FileNames(18) = NIARIVRETRDATA
                             DFNUM(4)=18
                             DFTXT(4)=13
                             DFSIGN(4)=1 
C                         FileNames(38) = NIARIVICEDATA  
                             DFNUM(5)=38
                             DFTXT(5)=14
                             DFSIGN(5)=1 
C                         FileNames(39) = NIARIVWEEDDATA  
                             DFNUM(6)=39
                             DFTXT(6)=15
                             DFSIGN(6)=1 
C                         FileNames(29) = DETRIVFLOWDATA  
                             DFNUM(7)=29
                             DFTXT(7)=18
                             DFSIGN(7)=1 
C                         FileNames(3)  = ERINBSDATA
                             DFNUM(8)=3
                             DFTXT(8)=1
                             DFSIGN(8)=1 
C                         FileNames(7)  = ERICONUSEDATA
                             DFNUM(9)=7
                             DFTXT(9)=4
                             DFSIGN(9)= -1 
C                         FileNames(11) = ERIOTHERSUPDATA
                             DFNUM(10)=11
                             DFTXT(10)=3
                             DFSIGN(10)=1 
C                         FileNames(15) = ERIOTHERDIVDATA
                             DFNUM(11)=15
                             DFTXT(11)=2
                             DFSIGN(11)=1 
         DFNUM(12)=0
         CALL SUMINDATA( DFNUM, DFTXT, DFSIGN, 3, 3, 6, 33, 
     &      RoutingInterval(2), BeginDate) 
      ENDIF
C
C
C                     Lake Ontario
      I = INDEX(FileNames(57),LLT) + 1
      IF(LakesToSolve(5) .AND. FileNames(57)(I:I+6) .EQ. 'default') THEN
C        Inflow/outflow summary file not specified, must build it
C
         IF(TSDVerbosity .GE. 2 ) THEN
            WRITE(CTMP,3456) 'Lake Ontario', LLT
            CALL LogIt(LOGCD,0,0)
         ENDIF
C
         IF(UStLawrSolution .EQ. 2)THEN
C           Using Plan 1958-D, set to first Sat
            BeginDate(1) = P58DStartDay
            BeginDate(2) = P58DStartMonth
            BeginDate(3) = P58DStartYear
         ENDIF
C
C        Initialize DFNUM with the index numbers to Lakes Ontario files in FileNames()
C                         FileNames(30) = NIARIVFLOWDATA       19
                             DFNUM(1)=30
                             DFTXT(1)=19
                             DFSIGN(1)=1 
C                         FileNames(21) = WELLANDCANALDATA       8
                             DFNUM(2)=21
                             DFTXT(2)=8
                             DFSIGN(2)=1 
C                         FileNames(64) = ONTDEVIATIONDATA      22
                             DFNUM(3)=64
                             DFTXT(3)=22
                             DFSIGN(3)=1 
C                         FileNames(58) = STLRIVFLOWDATA        30
                             DFNUM(4)=58
                             DFTXT(4)=30
                             DFSIGN(4)=1 
C                         FileNames(62) = LSLRETRDATA           24 
                             DFNUM(5)=62
                             DFTXT(5)=24
                             DFSIGN(5)=1 
C                         FileNames(66) = STLOUISFLOWDATA       23
                             DFNUM(6)=66
                             DFTXT(6)=23
                             DFSIGN(6)=1 
C                         FileNames(63) = LSROUGHNESSDATA       25 
                             DFNUM(7)=63
                             DFTXT(7)=25
                             DFSIGN(7)=1 
C                         FileNames(65) = DPRRIVFLOWDATA          26
                             DFNUM(8)=65
                             DFTXT(8)=26
                             DFSIGN(8)=1 
C                         FileNames(77) = MIIRIVFLOWDATA          26
                             DFNUM(9)=77
                             DFTXT(9)=35
                             DFSIGN(9)=1 
C                         FileNames(67) = ONTGATESDATA          21
                             DFNUM(10)=67
                             DFTXT(10)=21
                             DFSIGN(10)=1 
C                         FileNames(74) = ONTEOPLEVELDATA
                             DFNUM(11)=74
                             DFTXT(11)=32
                             DFSIGN(11)=1 
C                         FileNames(75) = ONTFORCEFLOWDATA        3
                             DFNUM(12)=75
                             DFTXT(12)=33
                             DFSIGN(12)=1 
C                         FileNames(76) = LSLONTDIFFDATA        3
                             DFNUM(13)=76
                             DFTXT(13)=34
                             DFSIGN(13)=1 
C                         FileNames(78) = FLOWDPMIDATA
                             DFNUM(14)=78
                             DFTXT(14)=37
                             DFSIGN(14)=1 
C                         FileNames(79) = FLOWSTFRANDATA
                             DFNUM(15)=79
                             DFTXT(15)=38
                             DFSIGN(15)=1 
C                         FileNames(80) = FLOWSTMAURDATA
                             DFNUM(16)=80
                             DFTXT(16)=39
                             DFSIGN(16)=1 
C                         FileNames(81) = FLOWRICHDATA
                             DFNUM(17)=81
                             DFTXT(17)=40
                             DFSIGN(17)=1 
C                         FileNames(82) = ROUGHJETTY1DATA
                             DFNUM(18)=82
                             DFTXT(18)=41
                             DFSIGN(18)=1 
C                         FileNames(83) = ROUGHVARENDATA
                             DFNUM(19)=83
                             DFTXT(19)=42
                             DFSIGN(19)=1 
C                         FileNames(84) = ROUGHSORELDATA
                             DFNUM(20)=84
                             DFTXT(20)=43
                             DFSIGN(20)=1 
C                         FileNames(85) = ROUGHPIERREDATA
                             DFNUM(21)=85
                             DFTXT(21)=44
                             DFSIGN(21)=1 
C                         FileNames(86) = ROUGH3RIVERSDATA
                             DFNUM(22)=86
                             DFTXT(22)=45
                             DFSIGN(22)=1 
C                         FileNames(87) = ROUGHBATISDATA
                             DFNUM(23)=87
                             DFTXT(23)=46
                             DFSIGN(23)=1 
C                         FileNames(88) = TIDALDATA
                             DFNUM(24)=88
                             DFTXT(24)=47
                             DFSIGN(24)=1 
C                         FileNames(89) = ONTTOTALSUPDATA
                             DFNUM(25)=89
                             DFTXT(25)=48
                             DFSIGN(25)=1 
C                         FileNames(90) = OTHER06ROUGHDATA
                             DFNUM(26)=90
                             DFTXT(26)=49
                             DFSIGN(26)=1 
C                         FileNames(91) = OTHER06FLOWDATA
                             DFNUM(27)=91
                             DFTXT(27)=50
                             DFSIGN(27)=1 
C
C                         FileNames(59) = STLRIVPPRETRDATA      27
                             DFNUM(28)=59
                             DFTXT(28)=27
                             DFSIGN(28)=1 
C                         FileNames(60) = STLRIVPPICEDATA       28 
                             DFNUM(29)=60
                             DFTXT(29)=28
                             DFSIGN(29)=1 
C                         FileNames(61) = STLRIVPPWEEDDATA      29
                             DFNUM(30)=61
                             DFTXT(30)=29
                             DFSIGN(30)=1 
C
C                         FileNames(53) = ONTNBSDATA             1
                             DFNUM(31)=53
                             DFTXT(31)=1
                             DFSIGN(31)=1 
C                         FileNames(54) = ONTCONUSEDATA          4
                             DFNUM(32)=54
                             DFTXT(32)=4
                             DFSIGN(32)= -1 
C                         FileNames(55) = ONTOTHERDIVDATA        2
                             DFNUM(33)=55
                             DFTXT(33)=2
                             DFSIGN(33)=1 
C                         FileNames(56) = ONTOTHERSUPDATA        3
                             DFNUM(34)=56
                             DFTXT(34)=3
                             DFSIGN(34)=1 
C
         DFNUM(35)=0
         CALL SUMINDATA( DFNUM, DFTXT, DFSIGN, 27, 3, 7, 57, 
     &      RoutingInterval(3), BeginDate) 
      ENDIF
C
      WRITE(CTMP, '(A)')
     &   ' End of pre-processing module for time-series data'//LLT
      CALL LogIt(LOGD,1,1)
C
C
      END
C
C
C----&------------------------------------------------------------------
C
      SUBROUTINE SUMINDATA(FNUM, TNUM, SNUM, NumOddBallFiles, 
     &   NumRetardFiles, LNUM, SumFileIndex, RoutingInt, BeginDate)
C
C
C     Reads the various time-series data files or repeating data lines
C     and creates summary input files that the regulation and routing
C     computations use.  FORMATINPUT determines what interval of
C     output to use, and fills in the index values for FNUM and TNUM.
C     FORMATINPUT calls SUMINDATA, also supplying LNUM, NumOddBallFiles,
C     and NumRetardFiles.  SNUM provides direction of flow; a value of 1
C     flow into lake, -1 is flow out of lake.
C
C     Frinstance, for Lake St. Clair it sees what computational interval
C     to use for the middle lakes routing (daily, weekly, qm, or monthly).  
C     Then it reads in all the Lake St. Clair water balance inputs (nbs,
C     cons use, other diversions, etc) and sums them up.
C
C     This routine is horribly convoluted and tricky and ineffecient.
C     Probably ought to totally re-do it at some other point.  But not
C     now, because it appears to work, and it's really tough, confusing,
C     programming work that i'm completely burned out on.
C
C     Precip and evap are not allowed to be expressed in units of height
C     (inches) yet.  They must be in flow units.  This is because of the
C     variable surface area implications; properly, the data would be
C     converted to flow units during the routing computations, when an
C     approximate surface area is available.  Since the program doesn't
C     assume a surface area, separate tallies would have to be reported 
C     in the summary inflow/outflow files for flows expressed as cm and m3s.
C     
C
C
      IMPLICIT NONE
C
C     Declarations regarding subroutine arguments:
C
      INTEGER      FNUM(48)
C                  Input. See description in FORMATINPUT()
C
      INTEGER      TNUM(48)
C                  Input. See description in FORMATINPUT()
C
      INTEGER      SNUM(48)
C                  Input. See description in FORMATINPUT()
C
      INTEGER      NumOddBallFiles
C                  Input. The number of the files in FNUM that
C                  are not ice/weed retardation files or inflow/outflow
C                  files.  They are not summed up.  These must 
C                  be the first in the list.
C      
      INTEGER      NumRetardFiles
C                  Input. The number of the files in FNUM that
C                  are ice/weed retardation files.  These must 
C                  be the first files after the oddballs, since they 
C                  are added up separately from the inflow/outflow 
C                  files, which are assumed be listed starting at
C                  NumOddBallFiles + NumRetardFiles + 1. 
C      
      INTEGER      LNUM
C                  Input. The lake being summarized, corresponding
C                  to index in array LAKENAME().
C      
      INTEGER      SumFileIndex
C                  Input.  The index number within FileNames corresponding
C                  to the summary input file of the lake.  31 for M-H, 
C                  32 for St Clair, 33 for Erie.
C      
      INTEGER      RoutingInt
C                  (Input) Routing interval to construct input summary for.
C                  1=daily, 2=weekly, 3=quarter-monthly, 4=monthly.
C
      INTEGER      BeginDate(3)
C                  Input. See description in FORMATINPUT()
C
C
C     Declarations regarding local variables or functions:
C
      INTEGER      F, I, J, K, L, M, N, YR, MON, DAY, QTR
      CHARACTER    INLINE*256
C
C                  Scratch.
C
      CHARACTER    ErrTmp*184
C                  Buffer to hold an error message
C
C      LOGICAL      GOTEOF(21)
C                  Flag for each input file to indicate if there are no more 
C                  lines available in the file, in order to provide more
C                  specific information in case get a READ error.  Set to 
C                  .TRUE. when the READ() returns EOF.  
C
C                  This didn't work out as spiffy as hoped, and is not currently 
C                  implemented.  The unix and pc compilers behaved differently, 
C                  with READ() sometimes acknowledging EOF only if the file was 
C                  actually terminated by an EOF character, and getting confused
C                  by the presence or absence of CR characters.  So couldn't 
C                  always tell if READ() was successful or unsuccessful when and
C                  if the EOF was encountered.  Also, blank lines or comment 
C                  lines added to the mess.  
C 
C                  So the way it works now, if READ() returns EOF, jumps down to
C                  label 955.  This happens prematurely (sometimes!) if the last
C                  data line isn't properly terminated.  So if the program can't
C                  read the last line of some data file, make sure that it ends
C                  with a LF (unix), or BOTH a CR and LF (pc's).  Or just edit
C                  the stupid thing and keeping hitting return at the bottom of
C                  the file until you feel satisfied.  That'll fix it!
C
      INTEGER      IOS
C                  Used to hold the return value from file reads.
C
      INTEGER      FileLine(48)
C                  Holds number of last line read from input data files.
C
      INTEGER      InFileHandle(48)
C                  FORTRAN unit spec for the input data files to summarize.
C
      LOGICAL      RepeatingValues(48)
C                  Flag for each input data source to indicate whether source
C                  is a file (.FALSE.) or repeating data (.TRUE.).
C
      LOGICAL      OpenedFile(48)
C                  Flag for each input data source to indicate whether source
C                  is an opened file (.TRUE.).
C
      LOGICAL      NEEDMORE(48)
C                  Flag for each input file to indicate if need to read another
C                  line.  Set to .TRUE. as use the final value on an input line.
C
      INTEGER      LASTROW(48)
C                  Contains the row number of the last data used from the scratch
C                  space for each input data file.
C                  Gets initialized to what the row number will
C                  be of the previous value used from the scratch space
C                  for each input data file.  For example, consider a 
C                  quarter-monthly routing starting in the third quarter
C                  of November.  If the data file contains monthly data,
C                  data lines read into the scratch space will have 12
C                  values.  LastRow(F) would be initialized to 10 (or 
C                  thereabouts - the incrementing gets kinda tricky, see
C                  comments in the code).
C
C
      INTEGER      MaxRows(48)
C                  The number of data pieces in the scratch array.  Gets set to
C                  12 for monthly input files, or repeating monthly data.  Set
C                  to 48 for quarter-monthly data.
C
      INTEGER      PreviousYr(48), PreviousMon(48), PreviousQtrMon(48)
C                  Contains the for year, quarter, month, and 
C                  quarter-month (as applicable) for the last line read in.
C                  Get initialized to the year, quarter, 
C                  month, and quarter-month (as applicable) for the
C                  line in each input data file prior to the first
C                  line with data used in the run. Allows routine
C                  to check continuity of input data as it reads in 
C                  subsequent lines.
C
      INTEGER      LINE9, LINE0
C                  These variables contain the billions and ones of data lines read.
C
      INTEGER      EQTR
C                  This variable tells the quarter-month of the ending date (1-4).
C
      INTEGER      NewQMDays
C                  In a week that straddles two quarter-months, the number of days 
C                  in the second quarter-month.  Used only for weekly routing.
C
      LOGICAL      NEWQTR
C                  Flag used for weekly routing interval calculations, to 
C                  indicate whether a week ends in a different quarter-month.
C
      INTEGER      EWEEK
C                  This variable tells the day-of-month for the end of the week 
C                  containing the ending day.  Used only for weekly routing.
C
      INTEGER      ENDOFQM(4)
C                  Contains the ending day of each quarter-month.
C
      CHARACTER    CHNKST*16
C                  This is the calling subroutine name passed to ErrMess() in
C                  case of error messages.
C
      CHARACTER    DataFormat(48)
C                  This is 'C', 'W', 'D', 'Q', or 'M'; gets returned from 
C                  InitDataFile or InitRepeatData.
C
      INTEGER      UnitsType(48)
C                  This is the index to UnitType(48) for the units used in the 
C                  data files.
C
      INTEGER      NUMFILES
C                  This is the actual number of files given in FNUM(48).
C
      CHARACTER*8  DFHEADING(64)
C                  Column headings for use in summary file.
C
      INTEGER      DecimalsToPrint(20)
C                  The number of decimal positions used when writing numbers of
C                  a certain data type to a summary file of time-series data.  
C                  Order corresponds to that of the UnitType series of arrays. 
C                     
      CHARACTER    LineFormat*120
      INTEGER      FormatLength
C                  The FORTRAN format specifier to use for writing the data
C                  lines of the times-series summary.  It gets constructed
C                  early in the routine based on the routing interval being
C                  used and the units specific to each input file (ie., data
C                  expressed in meters gets printed as F9.2, while m3s gets
C                  printed according to F9.0).  FormatLength is simply the
C                  number of characters used in LineFormat.  
C                  
      CHARACTER    POUND      
C                  These are message fragments used repeatedly.
C
C     INTEGER      ENDSTR
C                  Function that returns the end position of the current word.
C 
      LOGICAL      InitDataFile
C                  Function that figures out the units and time interval of a
C                  time-series data file, and positions it so that the next line
C                  read is the first line of data.
C
      LOGICAL      InitRepeatData
C                  Function that figures out the units and time interval of a
C                  repeating time-series data sequence, and returns an array
C                  filled with the annual sequence of values.
C
      LOGICAL      READQMONTHLY, READMONTHLY, ReadDaily, ReadWeekly
C                  Functions that read lines from datafile in the appropriate 
C                  format.
C 
      INTEGER      LASTKC
C                  Function that returns the last non-whitespace keyboard characters.
C 
      INTEGER      DateSequence
C                  Assigns a day number to a given date
C
      INTEGER      GETQMONNUM 
C                  Determines what quarter-month a date is in
C
      INTEGER      GetDaysInQMon 
C                  Determines how many days are in a quarter-month
C
      CHARACTER*16 MakeIntervalName
C                  Returns a string with the interval name.  For output.
C 
C
C     Declarations regarding COMMON block variables:
C
      DOUBLE PRECISION FWRK
C                  This is a common block used for 64-bit calculations.
      DIMENSION FWRK(0:49,0:53)
C                  Dimensioned it this way to provide a 50x48 matrix for
C                  flowrate conversions; up to 31 daily values (and a monthly mean)
C                  for each of the up to 48 supply files to read in (columns 1-48); 
C                  Row 0 is for the mean values, while rows 1-49 are working space
C                  for each input file.  Column 0 (FWRK(0, 0) through FWRK(49, 0))
C                  is scratch.
C
      CHARACTER*3  MonthAbbrev(12)
      INTEGER      GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &             SupVerbosity, MLVerbosity, OntVerbosity
C                  See descriptions in Block Data Genral
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
      CHARACTER    CTMP*255,CWRK*80
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
C                  See descriptions in Block Data Genral
C
      CHARACTER*8  RCRUNTYP, IntervalName(4)
C                  See description in Block Data Genral
C
      LOGICAL      LakesToSolve(8)
      INTEGER      StartYear, StartMonth, StartDay
      INTEGER      IntervalNameLen(4)
      INTEGER      MonthLength
      CHARACTER*8  MonthLengthName(2)
C                  Declared only because in common block along with end dates. 
C
      CHARACTER*9  MLSolutionName(2)
C                  Strings set to 'interative' and 'matrix', used for output.
C
      CHARACTER*16 SMSolutionName(4), USLSolutionName(5)
C
      INTEGER      EndYear, EndMonth, EndDay
C                  Ending year, month, day.
C
      CHARACTER*21 LAKENAME(16), RiverName(16)
C                  See description in Block Data Genral
C
      INTEGER      NLN(16), LkToSlv2LN(8), NRN(16), LkSlv2RtLoopI(8)
C                  See descriptions in Block Data Genral
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
      DOUBLE PRECISION UnitConvFactor(20)
C                  See description in Block Data Genral
C
      INTEGER      TEN, FTPMI, INPFT, SecsInDay, MinInt, MaxInt,
     &             MissingVal, SecsInUniformMon
      DOUBLE PRECISION MPFT, HALF, ONE, ZERO
C                  See descriptions in Block Data Genral
C
C
C
C
      COMMON /MONAME/   MonthAbbrev
      COMMON /FSPACE/   FWRK
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
      COMMON /RCPARI/   GenVerbosity, TSDVerbosity, ChkVerbosity, 
     &                  SupVerbosity, MLVerbosity, OntVerbosity
      COMMON /FNAMEC/   FileNames, RepeatLines
      COMMON /FNAMEI/   RepeatIndex, NPFILZ
      COMMON /LNAMEC/   LAKENAME, RiverName
      COMMON /LNAMEI/   NLN, LkToSlv2LN, NRN, LkSlv2RtLoopI
      COMMON /UNITSC/   UnitType
      COMMON /UNITSI/   UnitTypeLen, BaseUnitType, UnitConvFactor
      COMMON /CONSTS/   MPFT, HALF, ONE, ZERO, TEN, FTPMI, INPFT,
     &                  SecsInDay, MinInt, MaxInt, MissingVal,
     &                  SecsInUniformMon
C
C
      DATA POUND /'#'/
C
      DATA (DFHEADING(i), i = 1,55)  / '     NBS',' Oth Div',' Oth Sup',  
     1' Con Use','  Runoff','  Precip','    Evap',' Welland',' Chicago',  
     2' Tot LLO',' LongLac','   Ogoki','  Retard','     Ice','    Weed', 
     3' StMarys',' StClair',' Detroit',' Niagara',' CGIP Lv',' reserve', 
     4' Deviate',' LouFlow',' LSLRetd','    LSLn','  Des Pr',' OntRetd', 
     5'  OntIce',' OntWeed',' OntFlow',' Side Ch','  EOP Lv',' ForcedQ', 
     6' SLODiff','   Mille','   Gates','  Q DPMI',' QStFran','  QStMau', 
     7'   QRich',' RJetty1','  RVaren','  RSorel',' RPierre',' R3River', 
     8' Rbatisc','   Tidal',' Onttnbs',' Other6n',' Other6Q','   empty',
     9' SupIdx1',' SupIdx2',' SupIdx3',' USS Adj'/
C
C     If it's m, cm, m, inches, or feet, want to use 2 decimals.
C     If it's LT10, use six decimals.  If it's 10m3s or tcfs,
C     use 2 decimals.  Otherwise print to use one decimal place.
      DATA DecimalsToPrint / 2, 1, 2, 1, 2, 2, 2, 1, 1, 1, 
     &                       1, 1, 1, 1, 1, 6, 1, 1, 1, 1 /
C
C
      CHNKST='SumInData'//LLT
C
      YR = BeginDate(3)
      MON = BeginDate(2)
      QTR = GETQMONNUM (YR, MON, BeginDate(1))
      CALL GETQTRMONBOUNDS (YR, MON, ENDOFQM)
      EQTR = GETQMONNUM (EndYear, EndMonth, EndDay)
C
C     This loop initializes each of the data files for a lake
      DO 6444, F=1,48
C        See if done yet (calling module set FNUM(F) to 0 at end)
         IF(FNUM(F) .EQ. 0) GOTO 6445
C
C        Init flag to see if need to close a file later on:
         OpenedFile(F) = .FALSE.
C
C        Initialize line counter for each potential file
         FileLine(F) = 0
C
C        Store length of filename in scratch variable K
         K = INDEX(FileNames(FNUM(F)),LLT) - 1
C
C        Initialize the various arrays containing info
C        regarding last values read in, etc., and position
C        data files that netx line read contains data for
C        first routing interval.
         IF(FileNames(FNUM(F))(1:K+1) .EQ. 'Unspecified'//LLT 
     &      .OR. FileNames(FNUM(F))(1:4) .EQ. 'NULL' )THEN
C
            DataFormat(F) = '0'
            UnitsType(F) = 20
C           So must use 0's or missing for the values (NULL input 
C           file).  Just set row 0 for that file to '0' or MissingVal, 
C           and don't use LASTROW, etc. to change it, by setting 
C           DataFormat(F) to '0'.  If DataFormat(F) == '0', then falls 
C           through the elseif's that determine how to calculate a value
C           for FWRK(0,F). So the value that we put in it here (ZERO) 
C           gets used for every summary input line.
C
            IF( F .LE. NumOddBallFiles) THEN
               FWRK(0,F) = MissingVal
            ELSE
               FWRK(0,F) = ZERO
            ENDIF
C
C           That was general rule.  Now to exceptions:
            IF(FNUM(F) .EQ. 21 .OR. FNUM(F) .EQ. 30 .OR. 
     &         FNUM(F) .EQ. 58 .OR. FNUM(F) .EQ. 62 .OR. 
     &         FNUM(F) .EQ. 64 .OR. FNUM(F) .EQ. 65 .OR. FNUM(F) .EQ. 66
     &         .OR. FNUM(F) .EQ. 76 .OR. FNUM(F) .EQ. 77 ) THEN
C              - file 21 (Welland)
C              - file 30 (Niagara)
C              - file 58 (St Lawrence Flow)
C              - file 62 (St Louis Retard)
C              - file 64 (Ontario Deviation Flow)
C              - file 65 (Des Pr flow)
C              - file 66 (St Louis Flow)
C              - file 76 (Lake St Louis/Ontario outflow diff)
C              - file 77 (Mille flow)
               FWRK(0,F) = ZERO        
            ENDIF
C
         ELSE
C           So there is data.  See if repeating data or must 
C           read a file:
C
            IF(FileNames(FNUM(F))(1:15) .EQ. 'Repeating Data'//LLT) THEN
C
               RepeatingValues(F) = .TRUE.
C
C              Run InitRepeatData() to position at first data line, get 
C              data units in UnitsType, the interval type in DataFormat, 
C              the previous (or current) element number in the scratch 
C              array, and various description terms for prior data line.  
C
               IF(.NOT. InitRepeatData( FNUM(F), 
     &            DataFormat(F), UnitsType(F), FWRK(1,F) ))THEN
C
                  WRITE(CTMP,9761) 
     &               ' Unable to process repeating time-series '//
     &               'input data for ', LAKENAME(LNUM)(1:NLN(LNUM)), LLT
                  GOTO 977
               ENDIF
C
C              If DataFormat(F) == 'C', then falls through the elseif's
C              that determine how to calculate a value for FWRK(0,F).
C              So the value that we put in it here (constant) gets used
C              for every summary input line.
               IF(DataFormat(F) .EQ. 'C') FWRK(0,F)= FWRK(1,F) * SNUM(F)
C
C
C               UnitsType(F) = 1
C               FWRK(0,F) = MPFT * TEN * TEN * TEN * TEN 
C               DataFormat(F) = '0'
            ELSE
C
               RepeatingValues(F) = .FALSE.
               OpenedFile(F) = .TRUE.
C
C              Assign file handle from 51-71 for (up to) 21 input files
               InFileHandle(F) = 50 + F
C
C              Run INITDATAFILE() to position at first data line, get 
C              data units in UnitsType, the interval type in DataFormat, 
C              the previous (or current) element number in the scratch 
C              array.
C
               IF(.NOT. InitDataFile(InFileHandle(F),FNUM(F), BeginDate,
     &            FileLine(F), DataFormat(F), UnitsType(F)))THEN
C
                  WRITE(CTMP,9761) 
     &               ' Unable to process time-series input file for ',
     &               LAKENAME(LNUM)(1:NLN(LNUM)), LLT
                  GOTO 977
               ENDIF
            ENDIF
C
            I = 0
C           Check that units make sense - first the oddball columns:
            IF( F .LE. NumOddBallFiles) THEN
C
C              - file 68 (Lake Superior forced outflows)
C              - file 69 (Lake Superior forced sidechannel flow)
C              - file 75 (Lake Ontario forced flow)
C              - file 76 (Lake St Louis/Ontario outflow diff)
C              - file 65 (Des Pr flow)
C              - file 77 (Mille flow)
C              - file 30 (Niagara)
C              - file 21 (Welland)
C              - file 64 (Ontario Deviation Flow)
C              - file 58 (St Lawrence Flow)
C              - file 62 (St Louis Retard)
C              - file 66 (St Louis Flow)

               IF( (FNUM(F) .EQ. 68 .OR. FNUM(F) .EQ. 69 .OR.
     &            FNUM(F) .EQ. 30 .OR. FNUM(F) .EQ. 21 .OR.
     &            FNUM(F) .EQ. 64 .OR. FNUM(F) .EQ. 58 .OR.
     &            FNUM(F) .EQ. 62 .OR. FNUM(F) .EQ. 66 .OR.
     &            FNUM(F) .EQ. 65 .OR. FNUM(F) .EQ. 77 .OR.
     &            FNUM(F) .EQ. 75 .OR. FNUM(F) .EQ. 76 ) .AND. 
     &            UnitsType(F) .GT. 4 ) THEN
C                 Must be flow units.    squawk.
                  I = 1080
C
C              - file 20 (Chippewa-Grass Island Pool levels)
C              - file 70 (Lake Superior forced EOP levels)
C              - file 71 (Lakes Michigan-Huron forced EOP levels)
C              - file 72 (Lakes St. Clair forced EOP levels)
C              - file 73 (Lakes Erie forced EOP levels)
C              - file 74 (Lakes Ontario forced EOP levels)
               ELSE IF( (FNUM(F) .EQ. 20 .OR. 
     &            (FNUM(F) .GE. 70 .AND. FNUM(F) .LE. 74)) .AND. 
     &            (UnitsType(F) .NE. 5 .AND. UnitsType(F) .NE. 6 .AND.    
     &            UnitsType(F) .NE. 7 .AND. UnitsType(F) .NE. 11 .AND.    
     &            UnitsType(F) .NE. 12 ) ) THEN
C                 Must be elev units.    squawk.
                  I = 1080
C              
C              - file 26 (U.S. Slip Ice-Weed Adjustments)
C              - file 63 (Long Sault roughness factors)
C                         FileNames(82) = ROUGHJETTY1DATA
C                         FileNames(83) = ROUGHVARENDATA
C                         FileNames(84) = ROUGHSORELDATA
C                         FileNames(85) = ROUGHPIERREDATA
C                         FileNames(86) = ROUGH3RIVERSDATA
C                         FileNames(87) = ROUGHBATISDATA
C                         FileNames(88) = TIDALDATA
C              Now check that "special" data types have correct units
C              LS Roughness, Gates Open, etc., can be entered with
C              units of 'NA', 'NONE', 'Missing', or 'LT10'.  
               ELSE IF( FNUM(F) .EQ. 63 .OR. FNUM(F) .EQ. 82 .OR. 
     &                  FNUM(F) .EQ. 83 .OR. FNUM(F) .EQ. 84 .OR. 
     &                  FNUM(F) .EQ. 85 .OR. FNUM(F) .EQ. 86 .OR. 
     &                  FNUM(F) .EQ. 87 .OR. FNUM(F) .EQ. 88 
     &                  .OR. FNUM(F) .EQ. 26 ) THEN
C                 Roughness units should be set to unit type 'LT10' so 
C                 that 6 decimals print out.  Missing data is not legitimate 
C                 for this column, and will print out as '*********' anyway.
                  UnitsType(F) = 16
C
C              - file 40 (Superior gates force open)
C              - file 67 (Ontario gates force open)
               ELSE IF( FNUM(F) .EQ. 40 .OR. FNUM(F) .EQ. 67 ) THEN
C                 Some sort of Gates Open data, use 'NONE' and print 
C                 with one decimal point. This leaves room in 9-column
C                 field for mixed occurrences of 16.0, 0.5, or -9999.0.
C                 Or maybe some other dimensionless value that only
C                 needs 1 decimal printed out.
                  UnitsType(F) = 18
               ENDIF
            ELSE IF( UnitsType(F) .GT. 4) THEN
C              Not an oddball column - Must be a flow unit (first four
C              unit types).
               I = 1080
            ENDIF
C
            IF (I .GT. 0) THEN
               WRITE(ErrTmp,6471)F+2,LAKENAME(LNUM)(1:NLN(LNUM)),LLT 
               GOTO 975
            ENDIF
         ENDIF
C
 6444    NEEDMORE(F) = .TRUE.
 6445 NUMFILES = F - 1
 6471 FORMAT(' Check column ', I2, ' of summary input file for ',A, A1)
C
C
C     Initialize various description terms for prior data line.                
C
      DO 6470, F = 1, NumFiles
C        Initialize the variables indicating last row read from scratch 
C        space, (which may be 0 at first), and the continuity registers 
C        (which must reflect the yr/mon/qtr/week/day of the line 
C        previous to the line containing data for the start date.
C              
         IF(DataFormat(F) .EQ. 'M') THEN
C
            MaxRows(F) = 12
C
C           Since the incoming data in 1-yr chunks, only need to refer
C           to year to describe previous line (not used for repeating
C           data).
            PreviousYr(F) = BeginDate(3) - 1
C
C           LASTROW(F) gets incremented in strange and mysterious
C           ways.  If starting in first quarter for quarter-monthly 
C           routing, or day 1 of month for daily routing, will get 
C           incremented too soon; same for weekly routing. 
C
C           So for weekly routing,
            IF(RoutingInt .EQ. 2 .OR. 
C           or for qmonthly routing starting in 2nd, 3rd, or 4th QM,
     &         (RoutingInt .EQ. 3 .AND. GETQMONNUM(BeginDate(3),
     &         BeginDate(2), BeginDate(1)) .NE. 1) .OR.
C           or for daily routing starting after first day of month,
     &         (RoutingInt .EQ. 1 .AND. BeginDate(1) .NE. 1)) THEN
               LASTROW(F) =  BeginDate(2)
            ELSE
               LASTROW(F) =  BeginDate(2) - 1
            ENDIF
C
         ELSE IF(DataFormat(F) .EQ. 'Q') THEN
C
C           Data gets read from files or repeating blocks 48 at a time:
            MaxRows(F) = 48
C
C           Init LASTROW(F) to the quarter-month previous to the one 
C           containing startday/startmonth since LASTROW(F) for daily, 
C           qmonthly, and monthly routing are incremented before use.
            LASTROW(F) = (BeginDate(2) - 1) * 4 + 
     &         GETQMONNUM (BeginDate(3), BeginDate(2), BeginDate(1)) - 1
C
C           XX  This is on the TODO list to be changed:  
C           XX BUT, if averaging qmonthly data into monthly values, want 
C           XX to include all qtrmonths of the month, even if starting 
C           XX in 2, 3, or 4th qmonth.  So knock LASTROW(F) back to 
C           XX start of prior month.
            IF(RoutingInt .EQ. 4) THEN
               LASTROW(F) = LASTROW(F) - MOD(LASTROW(F),4)
C 
C           ALSO BUT, if interpolating between qmonthly values for 
C           weekly routing, set LASTROW(F) to the month actually 
C           containing starting day, since the weekly loop is incremented
C           after using the value of LASTROW(F).
            ELSE IF (RoutingInt .EQ. 2 .OR. 
     &         (RoutingInt .EQ. 1 .AND. BeginDate(1) .NE. I)) THEN 
               LASTROW(F) = LASTROW(F) + 1
            ENDIF
C
C           PreviousQtrYr(F) and PreviousYr(F) are ignored if the data
C           is repeating; used only to determine continuity of data
C           read from files.   
C           New way, assuming all four lines for a qmonth present:
            PreviousYr(F) = BeginDate(3) - 1
C           
         ELSE IF(DataFormat(F) .EQ. 'W') THEN
C
C           Can't do much here because haven't read a line, and don't
C           know date that data week ends
C           MaxRows not used for reading weekly data (1 line at a time)
            MaxRows(F) = MissingVal
            LastRow(F) = 0
C
         ELSE IF(DataFormat(F) .EQ. 'D') THEN
C
C           Since the incoming data in 1-qmonth chunks, need to refer
C           to year, month, and qmonth to describe previous line (not 
C           used for repeating data).
            PreviousYr(F) = BeginDate(3)
            PreviousMon(F) = BeginDate(2)
            PreviousQtrMon(F) = 
     &         GetQMonNum(BeginDate(3), BeginDate(2), BeginDate(1)) - 1
            IF( PreviousQtrMon(F) .EQ. 0) THEN
               PreviousQtrMon(F) = 4
               PreviousMon(F) = BeginDate(2) - 1
               LASTROW(F) = BeginDate(1) - 1
               IF( BeginDate(2) .EQ. 1) THEN
                  PreviousYr(F) = BeginDate(3) - 1
                  PreviousMon(F) = 12
               ENDIF
            ELSE
               LASTROW(F) = BeginDate(1) - 
     &            ENDOFQM(PreviousQtrMon(F)) - 1
            ENDIF
C           MaxRows does not remain constant for daily data format.
C           Some lines have 7 values, some have 8.  Gets reset with
C           each READ.  At this point, put in a dummy value, which
C           tells that no reads have taken place yet.
            MaxRows(F) = MissingVal
C
C
         ELSE IF(DataFormat(F) .EQ. 'C') THEN
            MaxRows(F) = 0 
         ELSE IF(DataFormat(F) .EQ. '0') THEN
            MaxRows(F) = 0 
            CONTINUE
         ELSE 
            STOP '? Doh!'
         ENDIF
 6470    CONTINUE
C
C
      IF( TSDVerbosity .GE. 8) THEN
         WRITE(CTMP,'(A,A,A1)') 
     &      ' Generating FORMAT string for ', 
     &      ' summary input file lines . . . ', LLT
         CALL LogIt(LOGD,0,0)
      ENDIF
C     Must compose a FORMAT specification for writing the 
C     summarized time-series data, one line per routing interval.
C     FORMAT depends on the routing interval and type of data.  
C     The string should end up like:
C
C     Daily:    (I1.1,I9.9,I3,I3,I6,2F12.3,NumFilesF9.?)
C     Weekly:   (I1.1,I9.9,I3,I3,I6,2F12.3,NumFilesF9.?)
C     Qmonthly: (I1.1,I9.9,I3,I3,I6,2F12.3,NumFilesF9.?)
C     Monthly:  (I1.1,I9.9,I6,I6,   2F12.3,NumFilesF9.?)
C           where NumFiles is the number of input data files, and
C                 '?' is 0 for flow data, or 2 for levels data. 
C     I9.9 means nine-char int field, padded with up to 9 leading 0's 
C                          
C     Initialize string, and keep track of its length in FormatLength
      LineFormat='(I1.1,I9.9,'
      IF( RoutingInt .EQ. 4) THEN 
C        Monthly format doesn't need date or qmonth field
         LineFormat(12:)='I6'
         FormatLength = 24
      ELSE
         LineFormat(12:)='I3,I3'
         FormatLength = 27
      ENDIF
      LineFormat(FormatLength-10:)=',I6,2F12.3'
C     
C     Count number of consecutive files with same decimals to 
C     print in I.  Initialize to 0 
      I = 0
C     Store previous number of decimals in J.  Initialize to first file.
C     Make sure this isn't LS Roughness, because of special 
C     DecimalsToPrint for that data.
      J = DecimalsToPrint(UnitsType(1))  
C
 6487 FORMAT( ',', I2, 'F9.', I1 )
C
C     Loop through the files to summarize, writing the format 
C     appropriate for each one (whenever differs from the previously
C     appropriate format):
      DO 6488, F = 1, NumFiles
C        Set L to 2 if units for file are one of these:
         L = DecimalsToPrint(UnitsType(F)) 
C        If same as previous file, increment I:
         IF( L .EQ. J) THEN
            I = I + 1
         ELSE  
C           Changing the number of decimals printed:
C           Write a string segment for the previous number of
C           decimals to print, concatenated on CWRK:
            WRITE(LineFormat(FormatLength:), 6487) I, J
            FormatLength = FormatLength + 7
C           Reset I:
            I = 1
C           Reset J to the number of decimals now called for:
            J = L
         ENDIF
C        IF( TSDVerbosity .GE. 8) THEN
C           WRITE(CTMP,'(A, I3, A, I2, A, A, A1)') ' Column:', F,
C    &         ', decimals:', DecimalsToPrint(UnitsType(F)), 
C    &         '  FORMAT string: ', LineFormat(1:FormatLength), LLT
C           CALL LogIt(LOGD,0,0)
C        ENDIF
 6488    CONTINUE
C                          
C                          
      WRITE(LineFormat(FormatLength:), 6487) I, J
      FormatLength = FormatLength + 7
      LineFormat(FormatLength:FormatLength) = ')'
C
C     Temporarily store length of the summary filename in I for
C     use on the next few lines.
      I = INDEX(FileNames(SumFileIndex),LLT) - 1
C
C     write(*,'(A,I2)') 'tsdverb=', TSDVerbosity
      IF( TSDVerbosity .GE. 7) THEN
         WRITE(CTMP,'(A,A,A,A,A1)') 
     &      ' Opening summary input file for ',
     &      LAKENAME(LNUM)(1:NLN(LNUM)), ': ',
     &   FileNames(SumFileIndex)(1:I), LLT
         CALL LogIt(LOGD,1,1)
      ENDIF
C
C     Initialize the file for summary input
      OPEN(UNIT=49, FILE=FileNames(SumFileIndex)(1:I), IOSTAT=IOS)
      IF ( IOS .NE. 0 ) THEN
         CALL ErrMess(1068,CHNKST,LOGCDE)
         WRITE(CTMP,9014) IOS, FileNames(SumFileIndex)(1:I), LLT
         GOTO 977
      ENDIF
C
C     write an empty comment line to summary file
      WRITE(49, 9750) POUND
C
C     Save capitalized name of routing interval to CWRK
      CWRK = MakeIntervalName( RoutingInt, .TRUE., J)
C
C     Calculate leading spaces to center title in I
      I = (80 - (J + NLN(LNUM) + 22 + 3))/2
C
C     Write the line, starting with comment character, then padding
C     with spaces to center, and finally sending the text:
      CTMP(1:1) = POUND
      DO 6122, K = 2, I
 6122    CTMP(K:K) = ' ' 
C     Write title to summary file:                                     
      WRITE(49, '(A)') CTMP(1:I)//' '//LAKENAME(LNUM)(1:NLN(LNUM))//
     &   ' '//CWRK(1:J)//' '//'Inflow/Outflow Summary'
C
C     write an empty comment line to summary file
      WRITE(49, 9750) POUND
C
      WRITE(49, 9750) '#    Column Number    Heading     Source File'
      WRITE(49, '(A,I2,A,I2)') 
     &   '#          1       Net Inflow     sum of columns ', 
     &   NumRetardFiles + NumOddBallFiles + 3, ' - ', NUMFILES + 2
      WRITE(49, '(A,I2,A,I2)') 
     &   '#          2       Net Retard     sum of columns ',
     &   NumOddBallFiles + 3, ' - ',
     &   NumRetardFiles + NumOddBallFiles + 2
C
      DO 6374, F = 1, NUMFILES
C        Determine length of file description/name
         IF(FileNames(FNUM(F))(1:12) .EQ. 'Unspecified'//LLT)THEN
            K = 11
         ELSE IF(FileNames(FNUM(F))(1:15) .EQ.'Repeating Data'//LLT)THEN
            K = 14
         ELSE
            K = LastKC( FileNames(FNUM(F)), 256)
         ENDIF
C
 6374    WRITE(49,'(A1,I11,8X,A9,5X,A)') POUND, F+2, DFHEADING(TNUM(F)),
     &      FileNames(FNUM(F))(1:K)
C
C     write an empty comment line to summary file
      WRITE(49, 9750) POUND
C
C     Write column headings to summary file:                                     
      WRITE(49, '(A1,18X,3I11,34I9)') POUND, (I, I=1, NUMFILES + 2)
      IF ( RoutingInt .EQ. 4 ) CWRK(1:1) = ' '
      WRITE(49, '(A1,2X,A7,2X,A1,A,48A9)')
     &   POUND,'Period',CWRK(1:1),'  M  Year  Net Inflow  Net Retard', 
     &   (DFHEADING(TNUM(I)),I=1,NUMFILES)
C     WRITE(49, '(A1,30X,A3,9X,A3,21A9)') POUND, 'M3S', 'M3S',
C    &   (UnitType(UnitsType(I))(1:UnitTypeLen(UnitsType(I))),
C    &   I=1,NUMFILES)
      WRITE(49, '(A1,30X,A3,9X,A3,48A9)') POUND, 'M3S', 'M3S',
     &   (UnitType(BaseUnitType(UnitsType(I)) )
     &   (1:UnitTypeLen(BaseUnitType(UnitsType(I)))),
     &    I=1,NumOddBallFiles),
     &   (UnitType(UnitsType(I))(1:UnitTypeLen(UnitsType(I))),
     &   I=NumOddBallFiles + 1,NUMFILES)
C
C         
C     Prepare for the infinite loops (1002, 2002, 3002, and 4002) made 
C     for each of the RoutingInt values.  Inside the loops,  
C     determine value appropriate to RoutingInterval for each file, and
C     write them out to the summary file for each and every routing 
C     period.  Set up a number of counters:
C
C     LINE9 is billions of lines read and LINE0 is ones
      LINE9=0  
      LINE0=0
C
C
C
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    &++++  Write inflow/outflow summary for daily routing interval  +++
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(RoutingInt .EQ. 1) THEN
C                 
         DAY = BeginDate(1)
C
C        Get the number of days in each quarter-month
         CALL GETQTRMONBOUNDS (BeginDate(3), BeginDate(2), ENDOFQM)
C
 1002    DO 1333, LINE0 = 1, 1000000000
C
C           loop through each input data file
            DO 1310, F = 1, NUMFILES
C
C              Process monthly input data
               IF(DataFormat(F) .EQ. 'M') THEN
C                 Process monthly input data into FWRK and then
C                 to daily intervals.             
C
C                 See if monthly values from last read are all used up.
C                 If so, read some more in:
                  IF( NEEDMORE(F) ) THEN
                     IF( .NOT. RepeatingValues(F) )THEN
C                       Read another line from file
 1240                   FileLine(F) = FileLine(F) + 1
                        READ( InFileHandle(F), 9750, END=955, ERR=960,
     &                     IOSTAT=IOS) INLINE
                        IF( INLINE(1:1) .EQ. '#' ) GOTO 1240
                        PreviousYr(F) = PreviousYr(F) + 1
                        IF(.NOT. READMONTHLY(INLINE, J, FWRK(1,F)))THEN
                           CALL ErrMess(1077,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9776) J, IOS, LLT
                           GOTO 965
                        ENDIF
                        IF( J .NE. PreviousYr(F)) THEN
                           CALL ErrMess(1078,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9782) PreviousYr(F), J, LLT
                           GOTO 965
                        ENDIF
                     ENDIF
C                    reset LASTROW (unless first time through)
                     IF(LINE0 .GT. 1) LASTROW(F) = 0
                     NEEDMORE(F) = .FALSE.
                  ENDIF
C
C                 Ok, have data in scratch space, 'cuz LASTROW(F) 
C                 is less than MaxRows(F).
C                 Just repeat monthly value for all days in month
C                 See if still in same month
                  IF( DAY .EQ. 1 ) THEN
C                    new month!  move counter to next row
                     LASTROW(F) = LASTROW(F) + 1
                  ELSE IF( LASTROW(F) .GE. MaxRows(F) .AND. 
     &               DAY .EQ. ENDOFQM(4)) THEN
C                    Last day for using data from this line of file
                     NEEDMORE(F) = .TRUE.
                  ENDIF
C                 Copy the value in this row to row 0, for later output
                  FWRK(0,F) = FWRK(LASTROW(F),F) * SNUM(F) 
C
               ELSE IF(DataFormat(F) .EQ. 'Q') THEN
C                 Process quarter-monthly input data into FWRK and then
C                 to daily intervals.             
C
C                 See if qmonthly values from last read are all used up.
C                 If so, read some more in:
                  IF( NEEDMORE(F) ) THEN
                     IF( .NOT. RepeatingValues(F) )THEN
C                       Read four more lines from file
                        PreviousYr(F) = PreviousYr(F) + 1
                        DO 1253, M=1,4
 1250                      FileLine(F) = FileLine(F) + 1
                           READ( InFileHandle(F), 9750, END=955, 
     &                        ERR=960,IOSTAT=IOS) INLINE
                           IF( INLINE(1:1) .EQ. '#' ) GOTO 1250
                           IF(.NOT. READQMONTHLY(INLINE, J,L,FWRK(1,F)))
     &                        THEN
                              CALL ErrMess(1076,CHNKST//LLT,LOGCDE)
                              WRITE(CTMP,9766)  J, L, IOS, LLT
                              GOTO 965
                           ENDIF
                           IF( J .NE. PreviousYr(F) .OR. L .NE. M) THEN
                              CALL ErrMess(1157,CHNKST//LLT,LOGCDE)
                              WRITE(CTMP,9781) M, PreviousYr(F),L,J,LLT
                              GOTO 965
                           ENDIF
 1253                      CONTINUE
                     ENDIF
C                    reset LASTROW (unless first time through)
                     IF(LINE0 .GT. 1) LASTROW(F) = 0
                     NEEDMORE(F) = .FALSE.
                  ENDIF
C                 Ok, now have data in scratch space, so move counter to next row
                  IF( (QTR .GT. 1 .AND. DAY .EQ. ENDOFQM(QTR-1) + 1) 
     &               .OR. DAY .EQ. 1 ) THEN 
C                    new qtrmonth!  move counter to next row
                     LASTROW(F) = LASTROW(F) + 1
                  ELSE IF( LASTROW(F) .GE. MaxRows(F) .AND. 
     &               DAY .GE. ENDOFQM(4)) THEN
C                    Last qtrmonth for using data from this line of file
                     NEEDMORE(F) = .TRUE.
                  ENDIF
C                 Copy the value in this row to row 0, for later output
                  FWRK(0,F) = FWRK(LASTROW(F),F) * SNUM(F)            
C
               ELSE IF(DataFormat(F) .EQ. 'W') THEN
C                 Process weekly input data into daily timesteps
C                 Weekly repeating data not supported 
                  IF( NEEDMORE(F) ) THEN
C                    Read another line from file
 1256                FileLine(F) = FileLine(F) + 1
                     READ(InFileHandle(F), 9750, END=955, ERR=960,
     &                  IOSTAT=IOS) INLINE
                     IF( INLINE(1:1) .EQ. '#' ) GOTO 1256
C
                     IF(.NOT. ReadWeekly(INLINE,J,I,K,FWRK(1,F))) THEN
                        CALL ErrMess(1154,CHNKST//LLT,LOGCDE)
                        WRITE(CTMP,9777) J, I, K, IOS, LLT
                        GOTO 965
                     ENDIF
C                    Date in data file is for END of week.
C                    Using PreviousQtrMon to hold *DAY* for weekly file, not qm
                     IF(MaxRows(F) .EQ. MissingVal) THEN
C                       First time read a line, establish date
                        PreviousYr(F) = J
                        PreviousMon(F) = I
                        PreviousQtrMon(F) = K
C                       Figure out how many days of run are in first data week
                        MaxRows(F) = 1 + DateSequence(K,I,J) - 
     &                     DateSequence(BeginDate(1), BeginDate(2), 
     &                     BeginDate(3)) 
                     ELSE
                        PreviousQtrMon(F) = PreviousQtrMon(F) + 7
                        CALL NUMDAYS(PreviousYr(F),PreviousMon(F),L,M)
                        IF( PreviousQtrMon(F) .GT. M) THEN
                           PreviousQtrMon(F) = PreviousQtrMon(F) - M
                           PreviousMon(F) = PreviousMon(F) + 1
                           IF( PreviousMon(F) .GT. 12) THEN
                              PreviousMon(F) = 1
                              PreviousYr(F) = PreviousYr(F) + 1
                           ENDIF
                        ENDIF
C                       Test to see if it's the correct line
                        IF(J .NE. PreviousYr(F) .OR. 
     &                     I .NE. PreviousMon(F) .OR. 
     &                     K .NE. PreviousQtrMon(F)) THEN
                           CALL ErrMess(1155,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9784) PreviousYr(F), 
     &                        PreviousMon(F), PreviousQtrMon(F), 
     &                        J, I, K, LLT
                           GOTO 965
                        ENDIF
                        MaxRows(F) = 7
                     ENDIF
                     LastRow(F) = 0
C                    Reset line-reading flag 
                     NEEDMORE(F) = .FALSE.
                  ENDIF
                  LastRow(F) = LastRow(F) + 1
                  IF( LASTROW(F) .GE. MaxRows(F))THEN
                     NEEDMORE(F) = .TRUE.
                  ENDIF
C                 Change sign and leave in row 0, for later output
                  FWRK(0,F) = FWRK(1,F) * SNUM(F)                 
C
               ELSE IF(DataFormat(F) .EQ. 'D') THEN
C                 Process daily input data into daily timesteps
C
C                 See if daily values from last read are all used up.
C                 If so, read some more in:
                  IF( NEEDMORE(F) ) THEN
C                    Repeating data not implemented, but leave this check in
C                    place in case of future expansion 
                     IF( .NOT. RepeatingValues(F) )THEN
C                       Read another line from file
 1260                   FileLine(F) = FileLine(F) + 1
                        READ(InFileHandle(F), 9750, END=955, ERR=960,
     &                     IOSTAT=IOS) INLINE
C                       Read again if comment
                        IF( INLINE(1:1) .EQ. '#' ) GOTO 1260
C                       Increment line ID's
                        PreviousQtrMon(F) = PreviousQtrMon(F) + 1
                        IF( PreviousQtrMon(F) .GT. 4) THEN
                           PreviousQtrMon(F) = 1
                           PreviousMon(F) = PreviousMon(F) + 1
                           IF( PreviousMon(F) .GT. 12) THEN
                              PreviousMon(F) = 1
                              PreviousYr(F) = PreviousYr(F) + 1
                           ENDIF
                        ENDIF
C                       Parse the input line
                        IF(.NOT. ReadDaily(INLINE,J,I,K,FWRK(1,F))) THEN
                           CALL ErrMess(1143,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9775) J, I, K, IOS, LLT
                           GOTO 965
                        ENDIF
C                       Test to see if it's the correct line
                        IF(J .NE. PreviousYr(F) .OR. 
     &                     I .NE. PreviousMon(F) .OR. 
     &                     K .NE. PreviousQtrMon(F)) THEN
                           CALL ErrMess(1144,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9783) PreviousYr(F), 
     &                        PreviousMon(F), PreviousQtrMon(F), 
     &                        J, I, K, LLT
                           GOTO 965
                        ENDIF
C                    ELSE
C                       Read daily repeating values from character 
C                       buffer (character*8(365)) prepared in 
C                       InitRepeatData
                     ENDIF
C                    Reset LASTROW (unless first time through)
                     IF(LINE0 .GT. 1) LASTROW(F) = 0
C                    Reset line-reading flag 
                     NEEDMORE(F) = .FALSE.
C                    Determine number of values in array
                     MaxRows(F) = GetDaysInQMon(J, I, K)
                  ENDIF
C
C                 Ok, have data in scratch space, 'cuz LASTROW(F) 
C                 less than MaxRows(F), or just read it.
C                 Move counter to next row
                  LASTROW(F) = LASTROW(F) + 1
                  IF( LASTROW(F) .GE. MaxRows(F))THEN
                     NEEDMORE(F) = .TRUE.
                  ENDIF
C                 Copy the value in this row to row 0, for later output
                  FWRK(0,F) = FWRK(LASTROW(F),F) * SNUM(F)
C
               ELSE IF(DataFormat(F) .EQ. 'C') THEN
C                 Value should already be in FWRK(0,F), do nothing
               ELSE IF(DataFormat(F) .EQ. '0') THEN
C                 Value of 0.0 should already be in FWRK(0,F), do nothing
               ELSE
                  WRITE(CTMP, '(A,I2,A,A1)') ' Bust: DataFormat(',F,
     &               ') set to ',DataFormat(F)
                  PRINT *, CTMP
                  STOP 'oh oh'
               ENDIF
C         
 1310       CONTINUE    
C           Now have a number for each of (up to) 18 input data files
C           saved in row 0 of FWRK, so write out the values for this day
C         
C           Column 0 used for accumulators; use row 0 to inflow/outflow.
C           Initialize accumulator for inflow/outflow files
            FWRK(0,0) = ZERO
            DO 1308, I=NumRetardFiles + NumOddBallFiles + 1, NUMFILES
 1308          FWRK(0,0) = FWRK(0,I) * UnitConvFactor(UnitsType(I))
     &            + FWRK(0,0)
C
C           Column 0 used for accumulators; devote row 1 to retardation.
C           Initialize accumulator for ice/weed retardation files.
            FWRK(1,0) = ZERO
C           Convert and sum up the ins/outs from each input data file 
            DO 1309, I=NumOddBallFiles+1, NumRetardFiles+NumOddBallFiles
 1309          FWRK(1,0) = FWRK(0,I) * UnitConvFactor(UnitsType(I)) 
     &            + FWRK(1,0)
C
C           Write the sum and components to summary input file
C            WRITE(49, '(I1.1,I9.9,I3,I3,I6,2F12.3,18F9.0)')  
            WRITE(49, LineFormat(1:FormatLength))
     &         LINE9,LINE0,DAY,MON,YR,FWRK(0,0),FWRK(1,0),
     &         (FWRK(0,I) * UnitConvFactor(UnitsType(I)), 
     &            I=1, NumOddBallFiles),
     &         (FWRK(0,I), I=NumOddBallFiles + 1, NUMFILES)
C     &         (FWRK(0,I),I=1,NUMFILES)
C
C           Ok, was final interval, jump out of loop
            IF( EndYear .EQ. YR .AND. EndMonth .EQ. MON 
     &         .AND. EndDay .LE. DAY) GOTO 1990
C
            DAY = DAY + 1
            IF( DAY .GT. ENDOFQM(QTR) ) THEN 
               QTR = QTR + 1
               IF( QTR .GT. 4 ) THEN 
                  DAY = 1  
                  QTR = 1
                  MON = MON + 1
                  IF( MON .GT. 12 ) THEN 
                     MON = 1
                     YR = YR + 1
                  ENDIF
                  CALL GETQTRMONBOUNDS (YR, MON, ENDOFQM)
               ENDIF
            ENDIF
 1333    CONTINUE
C
C        Passed 1000000000 lines, increment LINE9
         LINE9 = LINE9 + 1
         GOTO 1002
C        End of infinite loop through daily summary lines
C
 1990    CONTINUE    
C     End of processing for writing daily summary
C
C
C
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    &+++  Write inflow/outflow summary for weekly routing interval  +++
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      ELSE IF(RoutingInt .EQ. 2) THEN
C        Store day that week starts in DAY
         DAY = BeginDate(1)
C
C        Get the number of days in each quarter-month
         CALL GETQTRMONBOUNDS (BeginDate(3), BeginDate(2), ENDOFQM)
C
C        End of period is six days after start
         EWEEK = DAY + 6
C
C        See if period ends in different quarter-month
         IF( EWEEK .GT. ENDOFQM(QTR) ) THEN
C           End of week is in next quarter-month   
            QTR = QTR + 1
C
C           See if next week is in a different month
            IF( QTR .GT. 4 ) THEN 
               EWEEK = EWEEK - ENDOFQM(4)
               IF( DAY .GT. ENDOFQM(4)) DAY = DAY - ENDOFQM(4)
               QTR = 1
               MON = MON + 1
               IF( MON .GT. 12 ) THEN 
                  MON = 1
                  YR = YR + 1
               ENDIF
               CALL GETQTRMONBOUNDS (YR, MON, ENDOFQM)
            ENDIF
         ENDIF
C
C        Compute days in week belonging to this qmonth
         NewQMDays = EWEEK
         IF( QTR .GT. 1 )NewQMDays = NewQMDays - ENDOFQM(QTR-1)
C
C        Adjust NewQMDays and set NEWQTR if week is entirely within qmonth
         IF( NewQMDays .LT. 7 )THEN
C           Week straddles a qmonth
            NEWQTR = .TRUE.
         ELSE
C           EWEEK is 7 or 8 days after end of qmonth
            NewQMDays = 7
            NEWQTR = .FALSE.
         ENDIF
C
 2002    DO 2333, LINE0 = 1, 1000000000
C              
            IF(TSDVerbosity .GE. 9) THEN
               WRITE(49,'(A1)') '#'
               WRITE(49,'(4(A,I3),I6,A,I2,A,4I3,A,L1)') '# Start=',
     &            DAY,' End=',EWEEK,' QM=',QTR,' Mon=',MON, YR,
     &            '   Days in new QM=', NewQMDays, '  QM End Days=',
     &            (ENDOFQM(J), J=1,4),'   NEWQTR=', NEWQTR
            ENDIF
C           loop through each input data file
            DO 2310, F = 1, NUMFILES
C
C              Process monthly input data
               IF(DataFormat(F) .EQ. 'M') THEN
C                 Process monthly input data into FWRK and then
C                 to daily intervals.             
C
C                 See if monthly values from last read are all used up.
C                 If so, read some more in:
 2239             IF( NEEDMORE(F) ) THEN
                     IF( .NOT. RepeatingValues(F) )THEN
C                    Read another line from file
 2240                   FileLine(F) = FileLine(F) + 1
                        READ( InFileHandle(F), 9750, END=955, ERR=960,
     &                     IOSTAT=IOS) INLINE
C                       Repeat reading until not a comment:
                        IF( INLINE(1:1) .EQ. '#' ) GOTO 2240
C                        IF(TSDVerbosity .GE. 9) THEN
C                           WRITE(49,'(A,I3,A,L1,A,I2,3(A,F12.3))')
C     &                        '#m2w1 Col ',F+2, '  NEEDMORE=', 
C     &                        NEEDMORE(F),'  lastrow=', LASTROW(F), 
C     &                        ' prev=', FWRK(MaxRows(F)+1,F), ' new=',
C     &                        FWRK(LASTROW(F),F),' weighted=',FWRK(0,F)
C                        ENDIF
C                       Try to parse the line according monthly format:
                        IF(.NOT. READMONTHLY(INLINE, J, FWRK(1,F)))THEN
                           CALL ErrMess(1077,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9776) J, IOS, LLT
                           GOTO 965
                        ENDIF
C                       Update counter for year:
                        PreviousYr(F) = PreviousYr(F) + 1
C                       See if values from new line are consecutive:
                        IF( J .NE. PreviousYr(F)) THEN
                           CALL ErrMess(1078,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9782) PreviousYr(F), J, LLT
                           GOTO 965
                        ENDIF
                     ENDIF
C                    Reset LASTROW (unless first time through)
                     IF(LINE0 .EQ. 1 .AND. 
     &                  LASTROW(F) .LE. MaxRows(F) ) THEN
C                       Seed row 0 with stday value, in case first week 
C                       wraps around year:
                        FWRK(0,F) = FWRK(LASTROW(F),F)
C                       LASTROW will be advanced if StartDay = 1, so 
C                       anticipate for first time through:
                        IF(EWEEK .EQ. 7)LASTROW(F) = LASTROW(F) - 1
                     ELSE
                        LASTROW(F) = 0
                     ENDIF
                     NEEDMORE(F) = .FALSE.
                  ENDIF
C
                  IF(TSDVerbosity .GE. 9) THEN
                     WRITE(49,'(A,I3,A,L1,A,I2,3(A,F12.3))')
     &                  '#m2w2 Col ',F+2, '  NEEDMORE=', NEEDMORE(F),
     &                  '  lastrow=', LASTROW(F), ' prev=', 
     &                  FWRK(MaxRows(F)+1,F),' new=',FWRK(LASTROW(F),F), 
     &                  ' weighted=', FWRK(0,F)
                  ENDIF
C                 Ok, have data in scratch space, 'cuz LASTROW(F) 
C                 less than MaxRows(F).
C                 Just repeat monthly value for all weeks in month
C                 See if still in same month
                  IF( EWEEK .LE. 7 ) THEN
C                    At least part of week in next month!  move counter
C                    to next row
                     LASTROW(F) = LASTROW(F) + 1
C                    See if need to read another line now 
                     IF( LASTROW(F) .GT. MaxRows(F) ) THEN
                        NEEDMORE(F) = .TRUE.  
                        GOTO 2239
                     ENDIF
C                    Prorate the value in row 0 for days in week from 
C                    new month and last month (value still in row 0)
                     FWRK(0,F) = (FWRK(0,F) * (7 - EWEEK) + 
     &                  FWRK(LASTROW(F),F) * EWEEK) * SNUM(F) / 7
                  ELSE
                     FWRK(0,F) = FWRK(LASTROW(F),F) * SNUM(F)
                  ENDIF
                  IF(TSDVerbosity .GE. 9) THEN
                     WRITE(49,'(A,I3,A,I2,2(A,F12.3))') '# Column ',F+2,
     &                  '  lastrow=',LASTROW(F), ' new=', 
     &                  FWRK(LASTROW(F),F), ' weighted=',FWRK(0,F)
                  ENDIF
C
               ELSE IF(DataFormat(F) .EQ. 'Q') THEN
C                 Process quarter-monthly input data into FWRK and then
C                 to weekly intervals.             
C
C                  IF(TSDVerbosity .GE. 9) THEN
C                     WRITE(49,'(A,I3,A,L1,A,I2,3(A,F12.3))')'# 49 Col ',
C     &                  F+2, '  NEEDMORE=', NEEDMORE(F), '  lastrow=', 
C     &                  LASTROW(F), ' prev=', FWRK(MaxRows(F)+1,F), ' new=',
C     &                  FWRK(LASTROW(F),F), ' weighted=', FWRK(0,F)
C                  ENDIF
C                 See if qmonthly values from last read are used up.
C                 If so, read some more in:
 2249             IF( NEEDMORE(F) ) THEN
C                    Save final value from last line in MaxRows(F), for 
C                    use in pro-rating a value for a week straddling two 
C                    qmonths:
                     FWRK(MaxRows(F)+1,F) = FWRK(MaxRows(F),F)
                     IF( .NOT. RepeatingValues(F) )THEN
C                       Read four more lines from file
                        PreviousYr(F) = PreviousYr(F) + 1
                        DO 2253, M=1,4
 2250                      FileLine(F) = FileLine(F) + 1
                           READ( InFileHandle(F), 9750, END=955, 
     &                        ERR=960,IOSTAT=IOS) INLINE
                           IF( INLINE(1:1) .EQ. '#' ) GOTO 2250
                           IF(.NOT. READQMONTHLY(INLINE, J,L,FWRK(1,F)))
     &                        THEN
                              CALL ErrMess(1076,CHNKST//LLT,LOGCDE)
                              WRITE(CTMP,9766) J, L, IOS, LLT
                              GOTO 965
                           ENDIF
C                          See if values from new line are consecutive:
                           IF( J .NE. PreviousYr(F) .OR. L .NE. M) THEN
                              CALL ErrMess(1157,CHNKST//LLT,LOGCDE)
                              WRITE(CTMP,9781) M, PreviousYr(F),L,J,LLT
                              GOTO 965
                           ENDIF
 2253                      CONTINUE
                     ENDIF
C                    Reset LASTROW (unless first time through,
C                    because data for start date is later on line)
                     IF(LASTROW(F) .GT. MaxRows(F) .OR. LINE0 .NE. 1)
     &                 LASTROW(F) = 0
C                    Reset flag showing that scratch array is refilled.
                     NEEDMORE(F) = .FALSE.
                  ENDIF
C----&------------------------------------------------------------------
C
C                     WRITE(*,'(A,I3,A,L1,A,I2,3(A,F12.3))')'# NQ Col ',
C     &                  F+2, '  NEEDMORE=', NEEDMORE(F), '  lastrow=', 
C     &                  LASTROW(F), ' prev=', FWRK(MaxRows(F),F), ' new=',
C     &                  FWRK(LASTROW(F),F), ' weighted=', FWRK(0,F)
C
C                 See if still in same qmonth
                  IF( NEWQTR )THEN
C                    At least part of week in next qmonth!
C                    If not first qmonth on input line, save prev qtr value  
C                    in MaxRows(F)+1 (which gets updated just after line 2249
C                    in case of a new qmonth).
                     IF(LASTROW(F) .GT. 0) 
     &                  FWRK(MaxRows(F)+1,F)=FWRK(LASTROW(F),F)
C                    move counter to next row
                     LASTROW(F) = LASTROW(F) + 1
C                    Will need to read another line now if start day in 
C                    final qmonth of a line
                     IF( LASTROW(F) .GT. MaxRows(F) ) THEN
C                        WRITE(*,'(A)') ' Jumping back to 2249'
                        NEEDMORE(F) = .TRUE.  
                        GOTO 2249
                     ENDIF
                  ENDIF
C                     WRITE(*,'(A,I3,A,L1,A,I2,3(A,F12.3))')'# EW Col ',
C     &                  F+2, '  NEEDMORE=', NEEDMORE(F), '  lastrow=', 
C     &                  LASTROW(F), ' prev=', FWRK(MaxRows(F)+1,F), ' new=',
C     &                  FWRK(LASTROW(F),F), ' weighted=', FWRK(0,F)
C                 Prorate the value in row 0 for days in week from new qmonth 
C                 and last qmonth (value still in row 0)
                  FWRK(0,F) = (FWRK(MaxRows(F)+1,F) * (7 - NewQMDays) + 
     &                  FWRK(LASTROW(F),F) * NewQMDays) * SNUM(F) / 7 
C
                  IF(TSDVerbosity .GE. 9) THEN
                     WRITE(49,'(A,I3,A,I2,3(A,F12.3))')'# Column ',
     &                  F+2, '  lastrow=', LASTROW(F), 
     &                  ' prev=', FWRK(MaxRows(F)+1,F), ' new=',
     &                  FWRK(LASTROW(F),F), ' weighted=', FWRK(0,F)
                  ENDIF
C                
               ELSE IF(DataFormat(F) .EQ. 'W') THEN
C                 Process weekly input data into weekly timesteps
C                 Only get one per line, so must read another in
C                 Repeating data not implemented, but leave this check in
C                 place in case of future expansion 
C                 Read another line from file
 2256             FileLine(F) = FileLine(F) + 1
                  READ(InFileHandle(F), 9750, END=955, ERR=960,
     &               IOSTAT=IOS) INLINE
                  IF( INLINE(1:1) .EQ. '#' ) GOTO 2256
C                 Copy value from last week into register for last week
                  FWRK(1,F) = FWRK(2,F)
                  IF(.NOT. ReadWeekly(INLINE,J,I,K,FWRK(2,F))) THEN
                     CALL ErrMess(1154,CHNKST//LLT,LOGCDE)
                     WRITE(CTMP,9777) J, I, K, IOS, LLT
                     GOTO 965
                  ENDIF
C
C                 Date in data file is for END of week.
C                 Using PreviousQtrMon to hold *DAY* for weekly file, not qm
                  IF(MaxRows(F) .EQ. MissingVal) THEN
C                    First time read a line, establish date
                     PreviousYr(F) = J
                     PreviousMon(F) = I
                     PreviousQtrMon(F) = K
C                    How many days of routing week are in first data week?
C                    This offset stays constant and needs no updating
                     MaxRows(F) = 1 + DateSequence(K,I,J) - 
     &                  DateSequence(BeginDate(1), BeginDate(2), 
     &                  BeginDate(3)) 
C                    InitDataFile supposed to verify starting date, but
C                    make a simple check anyway
                     IF(MaxRows(F) .LT. 1 .OR. MaxRows(F) .GT. 7)
     &                   STOP 'Apu'
C                    If file week is not same as routing week, read again
                     IF( MaxRows(F) .NE. 7) GOTO 2256
C                    Call it 0 instead of 7, so weighting week just read in
                     MaxRows(F) = 0  
                  ELSE
                     PreviousQtrMon(F) = PreviousQtrMon(F) + 7
                     CALL NUMDAYS(PreviousYr(F),PreviousMon(F),L,M)
                     IF( PreviousQtrMon(F) .GT. M) THEN
                        PreviousQtrMon(F) = PreviousQtrMon(F) - M
                        PreviousMon(F) = PreviousMon(F) + 1
                        IF( PreviousMon(F) .GT. 12) THEN
                           PreviousMon(F) = 1
                           PreviousYr(F) = PreviousYr(F) + 1
                        ENDIF
                     ENDIF
C                    Test to see if it's the correct line
                     IF(J .NE. PreviousYr(F) .OR. I .NE. PreviousMon(F) 
     &                  .OR. K .NE. PreviousQtrMon(F)) THEN
                        CALL ErrMess(1155,CHNKST//LLT,LOGCDE)
                        WRITE(CTMP,9784) PreviousYr(F), PreviousMon(F),  
     &                     PreviousQtrMon(F), J, I, K, LLT
                        GOTO 965
                     ENDIF
                  ENDIF
                  FWRK(0,F) = SNUM(F) * (MaxRows(F) * FWRK(1,F) + 
     &               (7-MaxRows(F)) * FWRK(2,F)) / 7
C
                  IF( TSDVerbosity .GE. 9) WRITE(49,
     &               '(A,I2,A,I2,A3,I5,A,F12.6,2(A,I1,A,F9.2),A1)')
     &               '# Col ',F+2,' mean for week ending ', EWEEK, 
     &               MonthAbbrev(MON),YR, ': ', FWRK(0,F),', based on ',
     &               MaxRows(F), ' days of ', FWRK(1,F), ' and ', 
     &               7-MaxRows(F), ' days of ', FWRK(2,F), LLT
C
C
               ELSE IF(DataFormat(F) .EQ. 'D') THEN
C                 Process daily input data into weekly timesteps
C
C                 Keep number of non-missing values in M
                  M = 0
C
                  FWRK(0,F) = ZERO
C
C                 Get values for each day in week
                  DO 2261, L = 1, 7  
C                    See if daily values from last read are all used up.
C                    If so, read some more in:
                     IF( NEEDMORE(F) ) THEN
C                       Repeating data not implemented, but leave this 
C                       check in place in case of future expansion 
                        IF( .NOT. RepeatingValues(F) )THEN
C                          Read another line from file
 2260                      FileLine(F) = FileLine(F) + 1
                           READ(InFileHandle(F), 9750, END=955, ERR=960,
     &                        IOSTAT=IOS) INLINE
C                          Read again if comment
                           IF( INLINE(1:1) .EQ. '#' ) GOTO 2260
C                          Increment line ID's
                           PreviousQtrMon(F) = PreviousQtrMon(F) + 1
                           IF( PreviousQtrMon(F) .GT. 4) THEN
                              PreviousQtrMon(F) = 1
                              PreviousMon(F) = PreviousMon(F) + 1
                              IF( PreviousMon(F) .GT. 12) THEN
                                 PreviousMon(F) = 1
                                 PreviousYr(F) = PreviousYr(F) + 1
                              ENDIF
                           ENDIF
C                          Parse the input line
                           IF(.NOT. ReadDaily( INLINE, J, I, K, 
     &                           FWRK(1,F))) THEN
                              CALL ErrMess(1143,CHNKST//LLT,LOGCDE)
                              WRITE(CTMP,9775) J, I, K, IOS, LLT
                              GOTO 965
                           ENDIF
C                          Test to see if it's the correct line
                           IF(J .NE. PreviousYr(F) .OR. 
     &                        I .NE. PreviousMon(F) .OR. 
     &                        K .NE. PreviousQtrMon(F)) THEN
                              CALL ErrMess(1144,CHNKST//LLT,LOGCDE)
                              WRITE(CTMP,9783) PreviousYr(F), 
     &                           PreviousMon(F), PreviousQtrMon(F), 
     &                           J, I, K, LLT
                              GOTO 965
                           ENDIF
C                       ELSE
C                          Read daily repeating values from character 
C                          buffer (character*8(365)) prepared in 
C                          InitRepeatData
                        ENDIF
C                       Reset LASTROW (unless first time through)
                        IF(MaxRows(F) .NE. MissingVal) LASTROW(F) = 0
C                       Reset line-reading flag 
                        NEEDMORE(F) = .FALSE.
C                       Determine number of values in array
                        MaxRows(F) = GetDaysInQMon(J, I, K)
                     ENDIF
C
C                    Ok, have data in scratch space, 'cuz LASTROW(F) 
C                    less than MaxRows(F), or just read it.
C                    Move counter to next row
                     LASTROW(F) = LASTROW(F) + 1
                     IF( LASTROW(F) .GE. MaxRows(F))NEEDMORE(F) = .TRUE.
C
C                    Accumulate the values in row 0
                     IF( IDNINT(FWRK(LASTROW(F),F)) .NE. MissingVal)THEN
                        FWRK(0,F)=FWRK(0,F) + FWRK(LASTROW(F),F)*SNUM(F)
                        M = M + 1   
                     ENDIF
                     IF(TSDVerbosity .GE. 9)WRITE( CWRK((L-1)*9+1:L*9),
     &                   '(F9.3)' ) FWRK(LASTROW(F),F)
 2261                CONTINUE 
C
                  IF( M .EQ. 0 )THEN
C                    All data missing.
                     FWRK(0,F) = MissingVal
C                    
                     IF( F .GT. NumOddBallFiles) THEN
                        CALL ErrMess(1145,CHNKST//LLT,LOGCDE)
                        WRITE(CTMP,9765) EWEEK, MonthAbbrev(MON), YR, F,
     &                     LAKENAME(LNUM)(1:NLN(LNUM)), LLT
                        GOTO 965
                     ENDIF
                  ELSE
                     IF( F .GT. NumOddBallFiles .AND. M .LT. 7) THEN
C                       Not all data present, but some is there
                        CALL ErrMess(1156,CHNKST//LLT,LOGCDE)
                        WRITE(CTMP,9764) 'weekly', EWEEK, 
     &                     MonthAbbrev(MON), YR, F+2, 
     &                     LAKENAME(LNUM)(1:NLN(LNUM)), M, 7, LLT
C                       CALL LogIt(LOGD,0,0)
C                       IF( TSDVerbosity .GE. 5 .AND. F .GT. 
C    &                    NumOddBallFiles) CALL LogIt(LOGC,0,0) 
                        IF( TSDVerbosity .GE. 6) WRITE(49, '(A)') 
     &                     '# '//CTMP(1:INDEX(CTMP,LLT)-1)
                        GOTO 965
                     ENDIF
C                    Convert to mean
                     FWRK(0,F) = FWRK(0,F) / M
                  ENDIF
                  IF( TSDVerbosity .GE. 9) WRITE(49,
     &               '(A,I2,A,I2,A3,I5,A,F12.6,A,I1,A,A,A1)')
     &               '# Col ',F+2,' mean for week ending ', EWEEK, 
     &               MonthAbbrev(MON),YR, ': ', FWRK(0,F),', based on ',
     &               M, ' (of 7) daily values: ', CWRK(1:63),LLT
C
               ELSE IF(DataFormat(F) .EQ. 'C') THEN
C                 Value of 0.0 should already be in FWRK(0,F), do nothing
               ELSE IF(DataFormat(F) .EQ. '0') THEN
C                 Value of 0.0 should already be in FWRK(0,F), do nothing
               ELSE
                  WRITE(CTMP, '(A,I2,A,A1)') ' Bust: DataFormat(',F,
     &               ') set to ',DataFormat(F)
                  PRINT *, CTMP
                  STOP 'oh oh'
               ENDIF
C         
 2310       CONTINUE    
C           Now have a number for each of (up to) 18 input data files saved in row 0 of
C           FWRK, so write out the values for this quartermonth
C         
C           Initialize accumulator for ice/weed retardation files
            FWRK(1,0) = ZERO
C           Convert and sum up the ins/outs from each input data file 
            DO 2309, I=NumOddBallFiles+1, NumRetardFiles+NumOddBallFiles
 2309          FWRK(1,0) = FWRK(0,I) * UnitConvFactor(UnitsType(I))
     &            + FWRK(1,0)
C
C           Initialize accumulator for inflow/outflow files
            FWRK(0,0) = ZERO
            DO 2308, I=NumRetardFiles+NumOddBallFiles+1,NUMFILES
C               WRITE (CTMP, '(A,I2,A,F12.6,A1)')'  I=',I, 
C     &            '  UnitConvFactor=', UnitConvFactor(UnitsType(I)), LLT
C               CALL LogIt(LOGD,0,0)
 2308          FWRK(0,0) = FWRK(0,I) * UnitConvFactor(UnitsType(I)) 
     &            + FWRK(0,0)
C           Write the sum and components to summary input file
C            WRITE(49, '(I1.1,I9.9,I3,I3,I6,2F12.3,18F9.0)')  
            WRITE(49, LineFormat(1:FormatLength))
     &         LINE9,LINE0,EWEEK,MON,YR,FWRK(0,0),FWRK(1,0),
     &         (FWRK(0,I) * UnitConvFactor(UnitsType(I)), 
     &            I=1, NumOddBallFiles),
     &         (FWRK(0,I), I=NumOddBallFiles + 1, NUMFILES)
C     &         (FWRK(0,I),I=1,NUMFILES)
C
C           Ok, was final interval, jump out of loop
            IF( (EndYear .EQ. YR .AND. EndMonth .EQ. MON 
     &         .AND. EndDay .LE. EWEEK) .OR. EndYear .LT. YR .OR.
     &         (EndYear .EQ. YR .AND. EndMonth .LT. MON)) GOTO 2990
C
C           Start next period at day after end of this week
            DAY = EWEEK + 1
            EWEEK = DAY + 6
C           See if next week is in a different quarter-month
            IF( EWEEK .GT. ENDOFQM(QTR) ) THEN
C              End of week is in next quarter-month   
               QTR = QTR + 1
               NEWQTR = .TRUE.
C              See if next week is in a different month
               IF( QTR .GT. 4 ) THEN 
                  EWEEK = EWEEK - ENDOFQM(4)
                  IF( DAY .GT. ENDOFQM(4)) DAY = DAY - ENDOFQM(4)
                  QTR = 1
                  MON = MON + 1
                  IF( MON .GT. 12 ) THEN 
                     MON = 1
                     YR = YR + 1
                  ENDIF
                  CALL GETQTRMONBOUNDS (YR, MON, ENDOFQM)
               ENDIF
            ELSE
               NEWQTR = .FALSE.
            ENDIF
C           Compute days in week belonging to this qmonth
            NewQMDays = EWEEK
            IF( QTR .GT. 1 )NewQMDays = NewQMDays - ENDOFQM(QTR-1)
C           Adjust if week is entirely within qmonth
            IF( NewQMDays .GT. 7 ) NewQMDays = 7
C
 2333    CONTINUE
C
C        Passed 1000000000 lines, increment LINE9
         LINE9 = LINE9 + 1
         GOTO 2002
C        End of infinite loop through weekly summary lines
C
 2990    CONTINUE    
C     End of processing for writing weekly summary
C
C
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    & Write inflow/outflow summary for quarter-monthly routing interval  
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      ELSE IF(RoutingInt .EQ. 3) THEN
C                 
C
 3002    DO 3333, LINE0 =  1, 1000000000
C
C           DAY not otherwised used for qmonthly routing; store
C           sequence number of end date of period in it.
            IF( EndYear .EQ. YR .AND. EndMonth .EQ. MON .AND. 
     &         EQTR .EQ. QTR) THEN
               DAY = DateSequence(EndDay, MON, YR)
            ELSE
               DAY = DateSequence(ENDOFQM(QTR), MON, YR)
            ENDIF
C
C           loop through each input data file
            DO 3310, F = 1, NUMFILES
C
C              Process monthly input data
               IF(DataFormat(F) .EQ. 'M') THEN
C                 Process monthly input data into FWRK and then
C                 to quarter-monthly intervals.             
C
C                 See if monthly values from last read are all used up.
C                 If so, read some more in:
                  IF( NEEDMORE(F) ) THEN
                     IF( .NOT. RepeatingValues(F) )THEN
C                       Read another line from file
 3240                   FileLine(F) = FileLine(F) + 1
                        READ( InFileHandle(F), 9750, END=955, ERR=960,
     &                     IOSTAT=IOS) INLINE
                        IF( INLINE(1:1) .EQ. '#' ) GOTO 3240
                        PreviousYr(F) = PreviousYr(F) + 1
                        IF(.NOT. READMONTHLY(INLINE, J, FWRK(1,F)))THEN
                           CALL ErrMess(1077,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9776) J, IOS, LLT
                           GOTO 965
                        ENDIF
                        IF( J .NE. PreviousYr(F)) THEN
                           CALL ErrMess(1078,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9782) PreviousYr(F), J, LLT
                           GOTO 965
                        ENDIF
                     ENDIF
C                    reset LASTROW (unless first time through)
                     IF(LINE0 .GT. 1) LASTROW(F) = 0
                     NEEDMORE(F) = .FALSE.
                  ENDIF
C
C                 Ok, have data in scratch space, 'cuz LASTROW(F) 
C                 less than MaxRows(F).
C                 Just repeat monthly value for all 4 qtrmonths
C                 See if still in same month
                  IF( QTR .EQ. 1 ) THEN
C                    new month!  move counter to next row
                     LASTROW(F) = LASTROW(F) + 1
                  ELSE IF(LASTROW(F) .GE. MaxRows(F) .AND.QTR .EQ.4)THEN
C                    Last qtrmonth for using data from this line of file
                     NEEDMORE(F) = .TRUE.
                  ENDIF
C                 TODO:  allow interpolation between previous and
C                   next month's values
C                 Copy the value in this row to row 0, for later output
                  FWRK(0,F) = FWRK(LASTROW(F),F) * SNUM(F) 
C
               ELSE IF(DataFormat(F) .EQ. 'Q') THEN
C                 Process quarter-monthly input data into FWRK and then
C                 to quarter-monthly intervals.             
C
C                 See if qmonthly values from last read are all used up.
C                 If so, read some more in:
                  IF( NEEDMORE(F) ) THEN
C----&------------------------------------------------------------------

                     IF( .NOT. RepeatingValues(F) )THEN
C                       Read four more lines from file
                        PreviousYr(F) = PreviousYr(F) + 1
                        DO 3253, M=1,4
 3250                      FileLine(F) = FileLine(F) + 1
                           READ( InFileHandle(F), 9750, END=955, 
     &                        ERR=960,IOSTAT=IOS) INLINE
                           IF( INLINE(1:1) .EQ. '#' ) GOTO 3250
                           IF(.NOT. READQMONTHLY(INLINE, J,L,FWRK(1,F)))
     &                        THEN
                              CALL ErrMess(1076,CHNKST//LLT,LOGCDE)
                              WRITE(CTMP,9766) J, L, IOS, LLT
                              GOTO 965
                           ENDIF
C                          See new lines are consecutive:
                           IF( J .NE. PreviousYr(F) .OR. L .NE. M) THEN
                              CALL ErrMess(1157,CHNKST//LLT,LOGCDE)
                              WRITE(CTMP,9781) M, PreviousYr(F),L,J,LLT
                              GOTO 965
                           ENDIF
 3253                      CONTINUE
                     ENDIF
                     NEEDMORE(F) = .FALSE.
C                    reset LASTROW (unless first time through)
                     IF(LINE0 .GT. 1) LASTROW(F) = 0
                  ENDIF
C                 Ok, now have data in scratch space, so move counter to next row
                  LASTROW(F) = LASTROW(F) + 1
                  IF( LASTROW(F) .GE. MaxRows(F) )THEN
                     NEEDMORE(F) = .TRUE.
                  ENDIF
C                 Copy the value in this row to row 0, for later output
                  FWRK(0,F) = FWRK(LASTROW(F),F) * SNUM(F)
C
               ELSE IF(DataFormat(F) .EQ. 'W') THEN
C                 Process weekly input data into quarter-monthly timesteps
C                 Carry over any leftover days from previous week
                  FWRK(0,F) =  LastRow(F) * FWRK(1,F)
C                 Keep count of weights in FWRK(2,F)
                  FWRK(2,F) = LastRow(F)
C                 
                  IF( TSDVerbosity .GE. 9) WRITE(CTMP,'(I1,A,F9.2)') 
     &               LastRow(F), ' * ', FWRK(1,F)
C                 Only get one per line, so must always read at least one line
                  DO 3257, N = 1, MaxInt
C                    Read another line from file
 3256                FileLine(F) = FileLine(F) + 1
                     READ(InFileHandle(F), 9750, END=955, ERR=960,
     &                  IOSTAT=IOS) INLINE
                     IF( INLINE(1:1) .EQ. '#' ) GOTO 3256
C                    Copy value from last week into register for last week
                     IF(.NOT. ReadWeekly(INLINE,J,I,K,FWRK(1,F))) THEN
                        CALL ErrMess(1154,CHNKST//LLT,LOGCDE)
                        WRITE(CTMP,9777) J, I, K, IOS, LLT
                        GOTO 965
                     ENDIF
C
C                    Date in data file is for END of week.
C                    Using PreviousQtrMon to hold *DAY* for weekly file, not qm
                     IF(MaxRows(F) .EQ. MissingVal) THEN
C                       First time read a line, establish date
                        PreviousYr(F) = J
                        PreviousMon(F) = I
                        PreviousQtrMon(F) = K
C                       Get day number of end of data week
                        MaxRows(F) = DateSequence(K,I,J) 
C                       Get number of days in data week in first qm
                        L= DateSequence(BeginDate(1), BeginDate(2), 
     &                     BeginDate(3))
                        IF( MaxRows(F) .GT. DAY) THEN
                           LastRow(F) = DAY - L + 1
                        ELSE
                           LastRow(F) =  MaxRows(F) - L + 1
                        ENDIF
C                       InitDataFile supposed to verify starting date, but
C                       make a simple check anyway
                        IF(LastRow(F) .LT. 1 .OR. LastRow(F) .GT. 7)
     &                      STOP 'Moe'
                     ELSE
                        PreviousQtrMon(F) = PreviousQtrMon(F) + 7
                        CALL NUMDAYS(PreviousYr(F),PreviousMon(F),L,M)
                        IF( PreviousQtrMon(F) .GT. M) THEN
                           PreviousQtrMon(F) = PreviousQtrMon(F) - M
                           PreviousMon(F) = PreviousMon(F) + 1
                           IF( PreviousMon(F) .GT. 12) THEN
                              PreviousMon(F) = 1
                              PreviousYr(F) = PreviousYr(F) + 1
                           ENDIF
                        ENDIF
C                       Test to see if it's the expected line
                        IF(J .NE. PreviousYr(F) .OR. I .NE. 
     &                     PreviousMon(F) .OR. K .NE. 
     &                     PreviousQtrMon(F)) THEN
                           CALL ErrMess(1155,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9784) PreviousYr(F), 
     &                        PreviousMon(F), PreviousQtrMon(F), 
     &                        J, I, K, LLT
                           GOTO 965
                        ENDIF
C                       Increment sequnce number for end of data week
                        MaxRows(F) = MaxRows(F) + 7
                        IF( MaxRows(F) .GT. DAY ) THEN
C                          End of data week is beyond end of routing qm
C                          Set LastRow to number of days within this qm
                           LastRow(F) = 7 - (MaxRows(F) - DAY)
                           IF(LastRow(F) .LT. 1 .OR. LastRow(F) .GT. 6) 
     &                        STOP 'Smithers'
                        ELSE
                           LastRow(F) = 7
                        ENDIF
                     ENDIF
C                    
                     IF( TSDVerbosity .GE. 9) 
     &                  WRITE(CTMP(14+(N-1)*17:14+N*17),'(A,I1,A,F9.2)')  
     &                  ' + ', LastRow(F), ' * ', FWRK(1,F)
                     FWRK(0,F) = FWRK(0,F) + LastRow(F) * FWRK(1,F)
                     FWRK(2,F) = FWRK(2,F) + LastRow(F)
 3257                IF(MaxRows(F) .GE. DAY)GOTO 3258
C                      
C                 Should always jump out; never can fall through here
                  STOP 'Mr. Burns'
C
 3258             FWRK(0,F) = FWRK(0,F) * SNUM(F) / FWRK(2,F)
C                 Reset LastRow to number of unused days in data week
C                 to be picked up at start of next routing qm
                  LastRow(F) = 7 - LastRow(F)
                  IF( TSDVerbosity .GE. 9) WRITE(49,
     &               '(A,I2,A,I2,A3,I5,A,F12.6,A,A,A,I2,A1)')
     &               '# Col ',F+2,' mean for qm ending ', ENDOFQM(QTR), 
     &               MonthAbbrev(MON),YR, ': ', FWRK(0,F),', based on ('
     &               ,CTMP(1:14+N*17), ') / ', IDNINT(FWRK(2,F)), LLT
C
C
C
               ELSE IF(DataFormat(F) .EQ. 'D') THEN
C                 Process daily input data into quarter-monthly timesteps
C                 This logic is supposed to handle averaging for partial 
C                 at the beginning or end of a run.  So that if your run 
C                 starts on 28Apr, the value for 4th qmonth in April is 
C                 the mean of daily values for 28, 29, and 30 Apr, 
C                 instead of all the days in the qmonth.
C
C                 Repeating data not implemented, but leave this 
C                 check in place in case of future expansion 
                  IF( .NOT. RepeatingValues(F) )THEN
C                    Read 1 line at a time (format= 1 qmonth/line)
 3260                FileLine(F) = FileLine(F) + 1
                     READ(InFileHandle(F), 9750, END=955, ERR=960,
     &                  IOSTAT=IOS) INLINE
C                    Read again if comment
                     IF( INLINE(1:1) .EQ. '#' ) GOTO 3260
C                    Increment line ID's
                     PreviousQtrMon(F) = PreviousQtrMon(F) + 1
                     IF( PreviousQtrMon(F) .GT. 4) THEN
                        PreviousQtrMon(F) = 1
                        PreviousMon(F) = PreviousMon(F) + 1
                        IF( PreviousMon(F) .GT. 12) THEN
                           PreviousMon(F) = 1
                           PreviousYr(F) = PreviousYr(F) + 1
                        ENDIF
                     ENDIF
C                    Parse the input line
                     IF(.NOT. ReadDaily( INLINE, J, I, K, 
     &                  FWRK(1,F))) THEN
                        CALL ErrMess(1143,CHNKST//LLT,LOGCDE)
                        WRITE(CTMP,9775) J, I, K, IOS, LLT
                        GOTO 965
                     ENDIF
C                    Test to see if it's the correct line
                     IF(J .NE. PreviousYr(F) .OR. I .NE. PreviousMon(F)
     &                  .OR. K .NE. PreviousQtrMon(F)) THEN
                        CALL ErrMess(1144,CHNKST//LLT,LOGCDE)
                        WRITE(CTMP,9783) PreviousYr(F), PreviousMon(F), 
     &                     PreviousQtrMon(F), J, I, K, LLT
                        GOTO 965
                     ENDIF
C                 ELSE
C                    Read daily repeating values from character 
C                    buffer (character*8(365)) prepared in 
C                    InitRepeatData
                  ENDIF
C                 Reset LASTROW (unless first time through)
                  IF(MaxRows(F) .NE. MissingVal) LASTROW(F) = 0
C
C                 Determine number of values to use in array
                  IF( EndYear .EQ. YR .AND. EndMonth .EQ. MON 
     &              .AND. EQTR .EQ. QTR) THEN
C                    It's the last qmonth - maybe don't want all days
C                    Get the number of days in each quarter-month
                     CALL GETQTRMONBOUNDS (EndYear, EndMonth, ENDOFQM)
                     IF(EQTR .EQ. 1)THEN
                        MaxRows(F) = EndDay
                     ELSE
                        MaxRows(F) = EndDay - ENDOFQM(EQTR - 1)
                     ENDIF
                  ELSE
C                    Use days up to last day of qmonth in average
                     MaxRows(F) = GetDaysInQMon(J, I, K)
                  ENDIF
C
C                 Keep number of non-missing values in M
                  M = 0
                  FWRK(0,F) = ZERO
                  IF( TSDVerbosity .GE. 9) WRITE( CWRK,'(8I9)')
     &                    (MissingVal, I=1,8)
C
                  DO 3264, I = LASTROW(F) + 1, MaxRows(F)
C
C                    Accumulate the values in row 0
                     IF( IDNINT(FWRK(I,F)) .NE. MissingVal)THEN
                        FWRK(0,F) = FWRK(0,F) + FWRK(I,F) * SNUM(F)
                        M = M + 1   
                     ENDIF
 3264                IF(TSDVerbosity .GE. 9) WRITE( CWRK((I-1)*9+1:I*9),
     &                   '(F9.3)' ) FWRK(I,F)
C
C                 Keep track of the values used in I and EWEEK
C                 I has been incremented an extra time:
                  I = I - 1 - LastRow(F)
C
                  IF( M .EQ. 0 )THEN
C                    All data missing.
                     FWRK(0,F) = MissingVal
C                    
                     IF( F .GT. NumOddBallFiles) THEN
                        CALL ErrMess(1145,CHNKST//LLT,LOGCDE)
                        CALL GETQTRMONBOUNDS (YR, MON, ENDOFQM)
                        WRITE(CTMP,9765) ENDOFQM(QTR), MonthAbbrev(MON),
     &                      YR, F, LAKENAME(LNUM)(1:NLN(LNUM)), LLT
                        GOTO 965
                     ENDIF
                  ELSE
                     IF( F .GT. NumOddBallFiles .AND. M .NE. I) THEN
C                       Not all data present, but some is there
                        CALL GETQTRMONBOUNDS (YR, MON, ENDOFQM)
                        CALL ErrMess(1156,CHNKST//LLT,LOGCDE)
                        WRITE(CTMP,9764) 'qmonthly', ENDOFQM(QTR), 
     &                     MonthAbbrev(MON), YR, F+2, 
     &                     LAKENAME(LNUM)(1:NLN(LNUM)), M, I, LLT
C                       CALL LogIt(LOGD,0,0)
C                       IF( TSDVerbosity .GE. 5 .AND. F .GT. 
C    &                    NumOddBallFiles) CALL LogIt(LOGC,0,0) 
                        IF( TSDVerbosity .GE. 6) WRITE(49, '(A)')
     &                      '# '//CTMP(1:INDEX(CTMP,LLT)-1)
                        GOTO 965
                     ENDIF
C                    Convert to mean
                     FWRK(0,F) = FWRK(0,F) / M
                  ENDIF
                  IF( TSDVerbosity .GE. 9) WRITE(49,
     &               '(A,I2,A,I1,A,A3,I5,A,F12.6,A,I1,A,I1,A,A,A1)')
     &               '# Col ',F+2,' mean for qmonth ', QTR, ' of ' ,
     &                MonthAbbrev(MON), YR, ': ', FWRK(0,F), 
     &               ', based on ', M, ' (of ', I, ') daily values: ', 
     &               CWRK(1:I*9), LLT
C
               ELSE IF(DataFormat(F) .EQ. 'C') THEN
C                 Value of 0.0 should already be in FWRK(0,F), do nothing
               ELSE IF(DataFormat(F) .EQ. '0') THEN
C                 Value of 0.0 should already be in FWRK(0,F), do nothing
C                 WRITE(49, '( A,I2,2(A,F13.6) )')  '#  F=',F,'  FWRK(0,F)=', 
C     &              FWRK(0,F), '  FWRK(1,F)=', FWRK(1,F)
               ELSE
                  WRITE(CTMP, '(A,I2,A,A1)') ' Bust: DataFormat(',F,
     &               ') set to ',DataFormat(F)
                  PRINT *, CTMP
                  STOP 'oh oh'
               ENDIF
C         
 3310       CONTINUE    
C           Now have a number for each of (up to) 18 input data files saved in row 0 of
C           FWRK, so write out the values for this quartermonth
C         
C           Initialize accumulator for ice/weed retardation files
            FWRK(1,0) = ZERO
C           Convert and sum up the ins/outs from each input data file 
            DO 3309, I=NumOddBallFiles+1, NumRetardFiles+NumOddBallFiles
 3309          FWRK(1,0) = FWRK(0,I) * UnitConvFactor(UnitsType(I))
     &            + FWRK(1,0)
C
C           Initialize accumulator for inflow/outflow files
            FWRK(0,0) = ZERO
            DO 3308, I=NumRetardFiles + NumOddBallFiles + 1, NUMFILES
C               WRITE (CTMP, '(A,I2,A,F12.6,A1)')'  I=',I, 
C     &            '  UnitConvFactor=', UnitConvFactor(UnitsType(I)), LLT
C               CALL LogIt(LOGD,0,0)
 3308          FWRK(0,0) = FWRK(0,I) * UnitConvFactor(UnitsType(I))
     &            + FWRK(0,0)
C           Write the sum and components to summary input file
C            WRITE(49, '(I1.1,I9.9,I3,I3,I6,2F12.3,18F9.0)')  
C            WRITE(*, '(A)')LineFormat(1:FormatLength)
C            WRITE(*, LineFormat(1:FormatLength))  
            WRITE(49, LineFormat(1:FormatLength))  
     &         LINE9,LINE0,QTR,MON,YR,FWRK(0,0),FWRK(1,0),
     &         (FWRK(0,I) * UnitConvFactor(UnitsType(I)), 
     &            I=1, NumOddBallFiles),
     &         (FWRK(0,I), I=NumOddBallFiles + 1, NUMFILES)
C     &         (FWRK(0,I),I=1,NUMFILES)
C
C
C            Ok, done, break (ok jump) out of loop
             IF( EndYear .EQ. YR .AND. EndMonth .EQ. MON 
     &          .AND. EQTR .LE. QTR) GOTO 3990
C
            QTR = QTR + 1
            IF( QTR .GT. 4 ) THEN 
               QTR = 1
               MON = MON + 1
               IF( MON .GT. 12 ) THEN 
                  MON = 1
                  YR = YR + 1
               ENDIF
               CALL GETQTRMONBOUNDS (YR, MON, ENDOFQM)
            ENDIF
 3333    CONTINUE
C
C        Passed 1000000000 lines, increment LINE9
         LINE9 = LINE9 + 1
         GOTO 3002
C        End of infinite loop through quarter-monthly summary lines
C
 3990    CONTINUE    
C     End of processing for writing quarter-monthly summary
C
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    &++  Write inflow/outflow summary for monthly routing interval  +++
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      ELSE IF(RoutingInt .EQ. 4) THEN
C
C        QTR, MON, and YR should still be set to start date at this point
C
C
 4002    DO 4333, LINE0 =  1, 1000000000
C
C           DAY not otherwised used for monthly routing; store
C           sequence number of end date of period in it.
            DAY = DateSequence(ENDOFQM(4), MON, YR)
C
C           loop through each input data file
            DO 4310, F = 1, NUMFILES
C
C              Process monthly input data
               IF(DataFormat(F) .EQ. 'M') THEN
C                 Process monthly input data into FWRK and then
C                 to monthly intervals.             
C
C                 See if monthly values from last read are all used up.
C                 If so, read some more in:
                  IF( NEEDMORE(F) ) THEN
                     IF( .NOT. RepeatingValues(F) )THEN
C                       Read another line from file
 4240                   FileLine(F) = FileLine(F) + 1
                        READ(InFileHandle(F), 9750, END=955, ERR=960,
     &                     IOSTAT=IOS) INLINE
                        IF( INLINE(1:1) .EQ. '#' ) GOTO 4240
                        PreviousYr(F) = PreviousYr(F) + 1
                        IF(.NOT. READMONTHLY(INLINE, J, FWRK(1,F)))THEN
                           CALL ErrMess(1077,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9776) J, IOS, LLT
                           GOTO 965
                        ENDIF
                        IF( J .NE. PreviousYr(F)) THEN
                           CALL ErrMess(1078,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9782) PreviousYr(F), J, LLT
                           GOTO 965
                        ENDIF
                     ENDIF
C                    reset LASTROW (unless first time through)
                     IF(LINE0 .GT. 1) LASTROW(F) = 0
                     NEEDMORE(F) = .FALSE.
                  ENDIF
C
C                 Ok, have data in scratch space, 'cuz LASTROW(F) 
C                 less than MaxRows(F).
C                 Move counter to next row
                  LASTROW(F) = LASTROW(F) + 1
                  IF( LASTROW(F) .GE. MaxRows(F))THEN
                     NEEDMORE(F) = .TRUE.
                  ENDIF
C                 XX TODO:  allow interpolation between previous and
C                 XX  next month's values
C                 Copy the value in this row to row 0, for later output
                  FWRK(0,F) = FWRK(LASTROW(F),F) * SNUM(F) 
C
               ELSE IF(DataFormat(F) .EQ. 'Q') THEN
C                 Process quarter-monthly input data into FWRK and then
C                 to monthly intervals.             
C
C                 See if qmonthly values from last read are all used up.
C                 If so, read some more in:
                  IF( NEEDMORE(F) ) THEN
                     IF( .NOT. RepeatingValues(F) )THEN
C                       Read four more lines from file
                        PreviousYr(F) = PreviousYr(F) + 1
                        DO 4253, M=1,4
 4250                      FileLine(F) = FileLine(F) + 1
                           READ( InFileHandle(F), 9750, END=955, 
     &                        ERR=960,IOSTAT=IOS) INLINE
                           IF( INLINE(1:1) .EQ. '#' ) GOTO 4250
                           IF(.NOT. READQMONTHLY(INLINE, J,L,FWRK(1,F)))
     &                        THEN
                              CALL ErrMess(1076,CHNKST//LLT,LOGCDE)
                              WRITE(CTMP,9766) J, L, IOS, LLT
                              GOTO 965
                           ENDIF
C                          See if values from new line are consecutive:
                           IF( J .NE. PreviousYr(F) .OR. L .NE. M) THEN
                              CALL ErrMess(1157,CHNKST//LLT,LOGCDE)
                              WRITE(CTMP,9781) M, PreviousYr(F),L,J,LLT
                              GOTO 965
                           ENDIF
 4253                      CONTINUE
                     ENDIF
                     NEEDMORE(F) = .FALSE.
C                    reset LASTROW (unless first time through)
                     IF(LINE0 .GT. 1) LASTROW(F) = 0
                  ENDIF
                  FWRK(0,F) = ZERO
C
C                 Ok, now have data in scratch space.  Take average of 
C                 next four rows, store in row 0 for later output.
                  DO 4255, J=1,4
C                     WRITE(49,'(A,I1,A,I2,A,F10.1)') '# J=', J, 
C     &                  '  LASTROW(F)=', LASTROW(F), 
C     &                  '  FWRK(LASTROW(F) + J, F)=', 
C     &                  FWRK(LASTROW(F) + J, F)
 4255                FWRK(0,F) = FWRK(LASTROW(F) + J, F) + FWRK(0,F)
C                 to do: see if need to use DBLE(4)
                  FWRK(0,F) = FWRK(0,F) * SNUM(F) / 4
C                  WRITE(49,'(A,F10.1)') '#  FWRK(0,F)=', FWRK(0,F)
C                 Advance counter by four
                  LASTROW(F) = LASTROW(F) + 4
                  IF( LASTROW(F) .GE. MaxRows(F) )THEN
                     NEEDMORE(F) = .TRUE.
                  ENDIF
C
               ELSE IF(DataFormat(F) .EQ. 'W') THEN
C                 Process weekly input data into monthly timesteps
C                 Carry over any leftover days from previous week
                  FWRK(0,F) =  LastRow(F) * FWRK(1,F)
C                 Keep count of weights in FWRK(2,F)
                  FWRK(2,F) = LastRow(F)
C                 
                  IF( TSDVerbosity .GE. 9) WRITE(CTMP,'(I1,A,F9.2)') 
     &               LastRow(F), ' * ', FWRK(1,F)
C                 Only get one per line, so must always read at least one line
                  DO 4257, N = 1, MaxInt
C                    Read another line from file
 4256                FileLine(F) = FileLine(F) + 1
                     READ(InFileHandle(F), 9750, END=955, ERR=960,
     &                  IOSTAT=IOS) INLINE
                     IF( INLINE(1:1) .EQ. '#' ) GOTO 4256
                     IF(.NOT. ReadWeekly(INLINE,J,I,K,FWRK(1,F))) THEN
                        CALL ErrMess(1154,CHNKST//LLT,LOGCDE)
                        WRITE(CTMP,9777) J, I, K, IOS, LLT
                        GOTO 965
                     ENDIF
C
C                    Date in data file is for END of week.
C                    Using PreviousQtrMon to hold *DAY* for weekly file, not qm
                     IF(MaxRows(F) .EQ. MissingVal) THEN
C                       First time read a line, establish date
                        PreviousYr(F) = J
                        PreviousMon(F) = I
                        PreviousQtrMon(F) = K
C                       Get day number of end of data week
                        MaxRows(F) = DateSequence(K,I,J) 
C                       Get number of days in data week in first mon
                        L= DateSequence(BeginDate(1), BeginDate(2), 
     &                     BeginDate(3))
                        IF( MaxRows(F) .GT. DAY) THEN
                           LastRow(F) = DAY - L + 1
                        ELSE
                           LastRow(F) =  MaxRows(F) - L + 1
                        ENDIF
C                       InitDataFile supposed to verify starting date, but
C                       make a simple check anyway
                        IF(LastRow(F) .LT. 1 .OR. LastRow(F) .GT. 7)
     &                      STOP 'Flanders'
                     ELSE
                        PreviousQtrMon(F) = PreviousQtrMon(F) + 7
                        CALL NUMDAYS(PreviousYr(F),PreviousMon(F),L,M)
                        IF( PreviousQtrMon(F) .GT. M) THEN
                           PreviousQtrMon(F) = PreviousQtrMon(F) - M
                           PreviousMon(F) = PreviousMon(F) + 1
                           IF( PreviousMon(F) .GT. 12) THEN
                              PreviousMon(F) = 1
                              PreviousYr(F) = PreviousYr(F) + 1
                           ENDIF
                        ENDIF
C                       Test to see if it's the expected line
                        IF(J .NE. PreviousYr(F) .OR. I .NE. 
     &                     PreviousMon(F) .OR. K .NE. 
     &                     PreviousQtrMon(F)) THEN
                           CALL ErrMess(1155,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9784) PreviousYr(F), 
     &                        PreviousMon(F), PreviousQtrMon(F), 
     &                        J, I, K, LLT
                           GOTO 965
                        ENDIF
C                       Increment sequence number for end of data week
                        MaxRows(F) = MaxRows(F) + 7
                        IF( MaxRows(F) .GT. DAY ) THEN
C                          End of data week is beyond end of routing mon
C                          Set LastRow to number of days within this mon
                           LastRow(F) = 7 - (MaxRows(F) - DAY)
                           IF(LastRow(F) .LT. 1 .OR. LastRow(F) .GT. 6) 
     &                        STOP 'Sideshow Bob'
                        ELSE
                           LastRow(F) = 7
                        ENDIF
                     ENDIF
C                    
                     IF( TSDVerbosity .GE. 9) 
     &                  WRITE(CTMP(14+(N-1)*17:14+N*17),'(A,I1,A,F9.2)')  
     &                  ' + ', LastRow(F), ' * ', FWRK(1,F)
                     FWRK(0,F) = FWRK(0,F) + LastRow(F) * FWRK(1,F)
                     FWRK(2,F) = FWRK(2,F) + LastRow(F)
 4257             IF(MaxRows(F) .GE. DAY)GOTO 4258
C                      
C                 Should always jump out; never can fall through here
                  STOP 'Krusty'
C
 4258             FWRK(0,F) = FWRK(0,F) * SNUM(F) / FWRK(2,F)
C                 Reset LastRow to number of unused days in data week
C                 to be picked up at start of next routing month
                  LastRow(F) = 7 - LastRow(F)
                  IF( TSDVerbosity .GE. 9) WRITE(49,
     &               '(A,I2,A,A3,I5,A,F12.6,A,A,A,I2,A1)')
     &               '# Col ',F+2,' mean for ', 
     &               MonthAbbrev(MON),YR, ': ', FWRK(0,F),', based on ('
     &               ,CTMP(1:14+N*17), ') / ', IDNINT(FWRK(2,F)), LLT
C
C
C
               ELSE IF(DataFormat(F) .EQ. 'D') THEN
C                 Process daily input data into monthly timesteps
C                 This logic is supposed to handle averaging for partial 
C                 at the beginning or end of a run.  So that if your run 
C                 starts on 28Apr, the value for  April is 
C                 the mean of daily values for 28, 29, and 30 Apr, 
C                 instead of all the days in the month.
C
C                 Keep number of non-missing values in M
                  M = 0
                  FWRK(0,F) = ZERO
C
C                 Just recycle EWEEK as a counter for days in the month
C                 that we want to average over.
                  EWEEK = 0  
C
C                 Repeating data not implemented, but leave this 
C                 check in place in case of future expansion 
                  IF( .NOT. RepeatingValues(F) )THEN
C
C                    Determine number of lines to read
                     IF( EndYear .EQ. YR .AND. EndMonth .EQ. MON) THEN 
                        J = EQTR
                     ELSE 
                        J = 4
                     ENDIF
                     IF( TSDVerbosity .GE. 9) WRITE(49,
     &                  '(A,I2,A,A3,I5,A,I1,A,I1,A1)')
     &                  '#    Col ',F+2,' monthly mean for ',
     &                  MonthAbbrev(MON), YR, 
     &                  ' includes qmonths ', QTR, ' thru ', J, LLT
C                    Loop through the qmonth's with data being used
                     DO 4267, L = QTR, J
C
C                       Read 1 line at a time (format= 1 qmonth/line)
 4260                   FileLine(F) = FileLine(F) + 1
                        READ(InFileHandle(F), 9750, END=955, ERR=960,
     &                     IOSTAT=IOS) INLINE
C                       Read again if comment
                        IF( INLINE(1:1) .EQ. '#' ) GOTO 4260
C                       Parse the input line
                        IF(.NOT. ReadDaily( INLINE, J, I, K, 
     &                     FWRK(1,F))) THEN
                           CALL ErrMess(1143,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9775) J, I, K, IOS, LLT
                           GOTO 965
                        ENDIF
C                       Test to see if it's the correct line
                        IF(J .NE. YR .AND. I .NE. MON 
     &                     .AND. K .NE. L) THEN
                           CALL ErrMess(1144,CHNKST//LLT,LOGCDE)
                           WRITE(CTMP,9783)YR,MON,L, J, I, K, LLT 
                           GOTO 965
                        ENDIF
C                       Reset LASTROW (unless first line read)
                        IF(MaxRows(F) .NE. MissingVal) LASTROW(F) = 0
C
C                       Determine number of values to use in array
                        IF( EndYear .EQ. YR .AND. EndMonth .EQ. MON 
     &                     .AND. EQTR .EQ. L) THEN
C                          Last qmonth of run - maybe don't want all days
C                          Get the number of days in each quarter-month
                           CALL GETQTRMONBOUNDS (EndYear, EndMonth, 
     &                        ENDOFQM)
                           IF(EQTR .EQ. 1)THEN
                              MaxRows(F) = EndDay
                           ELSE
                              MaxRows(F) = EndDay - ENDOFQM(EQTR - 1)
                           ENDIF
C                           write(*, '(A,I2,A,I2)') 
C     &                        ' MaxRows(', F, ')=', MaxRows(F)
                        ELSE
C                          Use days up to last day of qmonth in average
                           MaxRows(F) = GetDaysInQMon(J, I, K)
                        ENDIF
C
                        IF( TSDVerbosity .GE. 9) WRITE( CWRK,'(8I9)')
     &                     (MissingVal, I=1,8)
C
                        DO 4264, I = LASTROW(F) + 1, MaxRows(F)
C
C                          Accumulate the values in row 0
                           IF( IDNINT(FWRK(I,F)) .NE. MissingVal) THEN
                              FWRK(0,F)= FWRK(0,F) + FWRK(I,F) * SNUM(F)
                              M = M + 1    
                           ENDIF
 4264                      IF( TSDVerbosity .GE. 9) WRITE( 
     &                        CWRK((I-1)*9+1:I*9),'(F9.3)' ) FWRK(I,F)
C
C                       Keep track of the values used in I and EWEEK
C                       I has been incremented an extra time:
                        I = I - 1 - LastRow(F)
C
                        EWEEK = EWEEK + I
C
 4267                   IF( TSDVerbosity .GE. 9) WRITE(49,
     &                     '(A,I2,A,A3,I5,A,I2,A,I2,A,I2,A,A,A1)')
     &                     '#    Col ',F+2,' monthly mean for ',
     &                     MonthAbbrev(MON), YR, 
     &                     ' will use these daily values (', M, 
     &                     ' of ', EWEEK, ') for month ', L, ': ', 
     &                     CWRK(1:I*9), LLT
C
C  
C                 ELSE
C                    Read daily repeating values from character 
C                    buffer (character*8(365)) prepared in 
C                    InitRepeatData
                  ENDIF
C
                  IF( M .EQ. 0 )THEN
C                    All data missing.
                     FWRK(0,F) = MissingVal
C                    
                     IF( F .GT. NumOddBallFiles) THEN
                        CALL ErrMess(1145,CHNKST//LLT,LOGCDE)
                        CALL GETQTRMONBOUNDS (YR, MON, ENDOFQM)
                        WRITE(CTMP,9765) ENDOFQM(QTR), 
     &                     MonthAbbrev(MON), YR, F, 
     &                     LAKENAME(LNUM)(1:NLN(LNUM)), LLT
                        GOTO 965
                     ENDIF
                  ELSE
                     IF( F .GT. NumOddBallFiles .AND. M .NE. EWEEK) THEN
C                       Not all data present, but some is there
                        CALL ErrMess(1156,CHNKST//LLT,LOGCDE)
                        WRITE(CTMP,9764) 'monthly', ENDOFQM(4), 
     &                     MonthAbbrev(MON), YR, F+2, 
     &                     LAKENAME(LNUM)(1:NLN(LNUM)), M, EWEEK, LLT
C                       CALL LogIt(LOGD,0,0)
C                       IF( TSDVerbosity .GE. 5 .AND. F .GT. 
C    &                     NumOddBallFiles) CALL LogIt(LOGC,0,0) 
                        IF( TSDVerbosity .GE. 6) WRITE(49, '(A)')
     &                     '# '//CTMP(1:INDEX(CTMP,LLT)-1)
                        GOTO 965
                     ENDIF
C                    Convert to mean
                     FWRK(0,F) = FWRK(0,F) / M
                  ENDIF
                  IF( TSDVerbosity .GE. 9) WRITE(49,
     &               '(A,I2,A,A3,I5,A,F12.6,A,I2,A,I2,A,A,A1)')
     &               '# Col ',F+2,' mean for ' ,
     &                MonthAbbrev(MON), YR, ': ', FWRK(0,F), 
     &               ', based on ', M, ' (of ', EWEEK, 
     &               ') daily values: ', CWRK(1:I*9), LLT
C
               ELSE IF(DataFormat(F) .EQ. 'C') THEN
C                 Value of 0.0 should already be in FWRK(0,F), do nothing
               ELSE IF(DataFormat(F) .EQ. '0') THEN
C                 Value of 0.0 should already be in FWRK(0,F), do nothing
               ELSE
                  WRITE(CTMP, '(A,I2,A,A1)') ' Bust: DataFormat(',F,
     &               ') set to ',DataFormat(F)
                  PRINT *, CTMP
                  STOP 'oh oh'
               ENDIF
C         
 4310       CONTINUE    
C           Now have a number for each of (up to) 18 input data files saved in row 0 of
C           FWRK, so write out the values for this month
C         
C           Initialize accumulator for ice/weed retardation files
            FWRK(1,0) = ZERO
C           Convert and sum up the ins/outs from each input data file 
            DO 4309, I=NumOddBallFiles+1, NumRetardFiles+NumOddBallFiles
 4309          FWRK(1,0) = FWRK(0,I) * UnitConvFactor(UnitsType(I))
     &            + FWRK(1,0)
C
C           Initialize accumulator for inflow/outflow files
            FWRK(0,0) = ZERO
            DO 4308, I=NumRetardFiles + NumOddBallFiles + 1, NUMFILES
C               WRITE (CTMP, '(A,I2,A,F12.6,A1)')'  I=',I, 
C     &            '  UnitConvFactor=', UnitConvFactor(UnitsType(I)), LLT
C               CALL LogIt(LOGD,0,0)
 4308          FWRK(0,0) = FWRK(0,I) * UnitConvFactor(UnitsType(I))
     &            + FWRK(0,0)
C           Write the sum and components to summary input file
C           Note that if NumOddBallFiles < 1, will blow up!
            WRITE(49, LineFormat(1:FormatLength))  
     &         LINE9,LINE0,MON,YR,FWRK(0,0),FWRK(1,0),
     &         (FWRK(0,I) * UnitConvFactor(UnitsType(I)), 
     &            I=1, NumOddBallFiles),
     &         (FWRK(0,I), I=NumOddBallFiles + 1, NUMFILES)
C
C           Ok, done, break (ok jump) out of loop
            IF( EndYear .EQ. YR .AND. EndMonth .EQ. MON) GOTO 4990
C                    
            QTR = 1
            MON = MON + 1
            IF( MON .GT. 12 ) THEN 
               MON = 1
               YR = YR + 1
            ENDIF
            CALL GETQTRMONBOUNDS (YR, MON, ENDOFQM)
 4333    CONTINUE
C
C        Passed 1000000000 lines, increment LINE9
         LINE9 = LINE9 + 1
         GOTO 4002
C        End of infinite loop through monthly summary lines
C
 4990    CONTINUE    
C
C
C     End of processing for writing monthly summary
C    &++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
      ELSE
            STOP 'Bad arg to SUMINDATA'
      ENDIF
C
C                                       
      DO 922, F = 1, NUMFILES
C         write(*, '(A,I4,A,I4)')  '  F=', F, '   NumFiles=', NumFiles
  922    IF(OpenedFile(F)) CLOSE(InFileHandle(F))
      CLOSE(49)
C
C
      RETURN
C
C     Oops.  Reached EOF while still trying to read more data
  955 CALL ErrMess(1079,CHNKST//LLT,LOGCDE)
      WRITE(CTMP,9745) ' Ran out of data (passed EOF) for column ',
     &   F+2, ' during period ',LINE0, LLT
      CALL LogIt(LOGCDE,0,0)
      WRITE(CTMP,9015) IOS, FileLine(F),
     &      FileNames(FNUM(F))(1:INDEX(FileNames(FNUM(F)),LLT)),LLT
      GOTO 977
C
C     Oops.  Some problem reading file
  960 CALL ErrMess(1003,CHNKST//LLT,LOGCDE)
      WRITE(CTMP,9015) IOS, FileLine(F),
     &      FileNames(FNUM(F))(1:INDEX(FileNames(FNUM(F)),LLT)),LLT
      GOTO 977
C
  965 CALL LogIt(LOGCDE,1,0)
      WRITE(CTMP,9762) FileLine(F), 
     &      FileNames(FNUM(F))(1:INDEX(FileNames(FNUM(F)),LLT)),LLT
      CALL LogIt(LOGCDE,1,0)
      WRITE(CTMP,9760)INLINE(1:LASTKC(INLINE, 149)),LLT
      GOTO 977
C
C
C
  975 CALL ErrMess(I,CHNKST//LLT,LOGCDE)
      CTMP = ErrTmp
  977 CALL LogIt(LOGCDE,0,1)
      CALL BuildInfo
      CALL LogIt(LOGCDE,0,0)
      CALL SAVETOD
      CALL ADIOS()
C
 9745 FORMAT(A,I2,A,I9,A1)
 9014 FORMAT(' Got IOSTAT=',I7,' when accessing file ', A,'.',A1)
 9015 FORMAT(' Got IOSTAT=',I7,' when accessing line ', I7, ' of file ',
     &   A,'.',A1)
 9750 FORMAT(A)
 9760 FORMAT(A,A1)
 9761 FORMAT(A,A,A)
 9762 FORMAT(' Problem on line ', I7, ' of file ',A, A1)
 9764 FORMAT('  ... mean for ', A, ' routing period ending ', I2, A3, 
     &   I5, ' for column ', I2, ' in ', A, ' estimated from ',I2, 
     &   ' (of ', I2, ') data values.', A1)
 9765 FORMAT('  Missing data for routing period ending ', I2, A3, I5, 
     &   ' for column ', I2, ' in ', A, ' summary input data file.', A1)
 9766 FORMAT(' ReadQMonthly failed, returning Year=', I5,
     &   ' and QMonth=', I2,' (IOSTAT=', I7,')',A1)
 9775 FORMAT(' ReadDaily failed, returning Year=', I5, ' Month=', I3,
     &   ' and Qmonth=', I2, ' (IOSTAT=', I7, ')', A1)
 9776 FORMAT(' ReadMonthly failed, returning Year=', I5, 
     &   ' (IOSTAT=', I7,')',A1)
 9777 FORMAT(' ReadWeekly failed, returning Year=', I5, ' Month=', I3,
     &   ' and Day=', I3, ' (IOSTAT=', I7, ')', A1)
 9781 FORMAT(' Expected line ', I1, ' for year ', I5, ', but got line ',
     &   I1, ' of year ', I5, ' instead!', A1)
 9782 FORMAT(' Expected line for year ', I5, ', but got line for year ',
     &   I5, ' instead!', A1)
 9783 FORMAT(' Expected line for year/month/qmonth ', I5, '/', I2, '/', 
     &   I1, ', but got line for year:', I5, ' month:', I3, ' qmonth:',
     &   I2, ' instead!', A1)
 9784 FORMAT(' Expected line for year/month/day ', I5, '/', I2, '/', 
     &   I2, ', but got line for year:', I5, ' month:', I3, ' day:',
     &   I2, ' instead!', A1)
      END
C
C
C----+------------------------------------------------------------------
C----&------------------------------------------------------------------
C-----------------------------------------------------------------------
