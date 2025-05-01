C-----------------------------------------------------------------------      
      SUBROUTINE BuildInfo()
C
C
C     This subroutine writes one line of info pertaining to the compiler
C     used (or other build info) for use with LogIt().
C
C
      IMPLICIT NONE
      CHARACTER*(*) Builder, Buildate, CompilerName, SrcVer
C
      CHARACTER     CTMP*255, CWRK*80, LLT*1
C                   See descriptions in Block Data Genral
C
      COMMON /WKCHAR/  CTMP,LLT,CWRK
C
C
C     Please edit the next three lines to reflect this build
C
      PARAMETER (Builder= 'Caleb Wegener')
C      PARAMETER (Builder=       'John Koschik'         )
C     PARAMETER (Builder=       'Matt McPherson'         )
C     PARAMETER (Builder=       'Yin Fan'         )
C     PARAMETER (Builder=       'Paul M. Yu'         )
C
      PARAMETER (Buildate=      '14Jun2024'              )
C
      PARAMETER (SrcVer=        '2.02beta'                   )
C
C      PARAMETER (CompilerName=  
C     &   'Sun Fortran 95 8.3 SunOS_sparc 2007-05-03')
C      PARAMETER (CompilerName=      'SPARCompiler FORTRAN 3.0')
C      PARAMETER (CompilerName=      'Microsoft FORTRAN 5.0')
C      PARAMETER (CompilerName=      'Lahey LF90 4.5')
C      PARAMETER (CompilerName=      'Lahey LF95 4.5')
C      PARAMETER (CompilerName=      'Lahey F77L 5.0')
       PARAMETER (CompilerName=      'gfortran/gcc 11.3.0')
C      PARAMETER (CompilerName=      'gfortran/gcc 11.3.0')
C       PARAMETER (CompilerName=      'gfortran/gcc 4.4.0 20110419')
C       PARAMETER (CompilerName=      'gfortran/gcc 4.4.0 20110419')
C       PARAMETER (CompilerName=      'Microsoft PowerStation 4.0')
C
      CTMP = 'CGLRRM version '//SrcVer//' modified & compiled with '//
     &   CompilerName//' on '//Buildate//' by '//Builder//LLT
C
      END     
C
C
C-----------------------------------------------------------------------      
      SUBROUTINE MyDate(DSTRNG)
C
C
C     This subroutine returns the system date as a nine-character string 
C     of the form DDMMMYYYY.  
C
C
      IMPLICIT NONE
      CHARACTER*9 DSTRNG
C     ------------------------
C
C
C     Uncomment this section for any compiler that is
C     Fortran 90 compliant.
C     (Time90, Zone90, Values90 will not be used)
C
      INTEGER       Values90(8), M
      CHARACTER*2   D
      CHARACTER*4   Y
      CHARACTER*10  Date90, Time90, Zone90
      CHARACTER     MonthAbbrev(12)*3
      COMMON /MONAME/ MonthAbbrev
      CALL Date_And_Time(Date90, Time90, Zone90, Values90)
      READ(Date90, 9001) Y, M, D
      WRITE(DSTRNG, 9002) D, MonthAbbrev(M), Y
C
 9001 FORMAT(A4, I2, A2)
 9002 FORMAT(A2, A3, A4)
c
C     ------------------------
C
C
C     Uncomment this section for Lahey Fortran (F77L)
C
C
C     CHARACTER LADATE*8, MonthAbbrev(12)*3
C     INTEGER M,Y
C     COMMON /MONAME/ MonthAbbrev 
C
C     DATE returns 8-character string "MM/DD/YY"
C     CALL DATE(LADATE) 
C
C     READ(LADATE, 9504) M, Y
C9504 FORMAT(I2, 4X, I2) 
C
C     IF( Y .GE. 80 ) THEN
C        Y = Y + 1900
C     ELSE
C        Y = Y + 2000
C     ENDIF 
C
C     WRITE(DSTRNG,'(A2,A3,I4)')LADATE(3:4),MonthAbbrev(M),Y
C     ------------------------
C
C
C
C
C
C
C
C     Uncomment this section for SPARCompiler FORTRAN 3.0
C
C      CHARACTER*24 CTIME, DS24
C      INTEGER I, TIME
C      I = TIME()
C      DS24 = CTIME( I )
C      DSTRNG = DS24(9:10)//DS24(5:7)//DS24(21:24)
C     ------------------------
C
C
C     Uncomment this section for M$ FORTRAN 5.0 or M$ PowerStation
C
C      INTEGER*2 YR, MON, DAY
C      CHARACTER MonthAbbrev(12)*3
C      COMMON /MONAME/ MonthAbbrev
C 
C      CALL GETDAT(YR,MON,DAY)
C      WRITE(DSTRNG,'(I2.2,A3,I4)')DAY,MonthAbbrev(MON),YR
C
C
C
      END


C-----------------------------------------------------------------------      
      SUBROUTINE MyTime(TSTRNG)
C
C     This subroutine returns the system time as a Lahey-style string 
C     (ie, in the form HH:MM:SS).  
C
C
      IMPLICIT NONE
      CHARACTER*8 TSTRNG
C     ------------------------
C
C
C
C     Uncomment this section for any compiler that is
C     Fortran 90 compliant.
C     (Date90, Zone90, Values90 will not be used)
C
      INTEGER       Values90(8)
      CHARACTER*10  Date90, Time90, Zone90
      CALL Date_And_Time(Date90, Time90, Zone90, Values90)
      WRITE(TSTRNG, 9050) Time90(1:2), Time90(3:4), Time90(5:6)
 9050 FORMAT(A2, ':', A2, ':', A2)
C
C     ------------------------
C
C
C     Uncomment this section for Lahey Fortran F77L
C
C     CALL TIME(TSTRNG)
C
C     ------------------------
C
C
C
C
C     Uncomment this section for SPARCompiler FORTRAN 3.0
C
C     CHARACTER*24 CTIME, TS24
C     INTEGER I, TIME
C     I = TIME()
C     TS24 = CTIME( I )
C     TSTRNG = TS24(12:19)
C     ------------------------
C
C
C
C
C     Uncomment this section for M$ FORTRAN 5.0
C
C      INTEGER*2 HOUR, MIN, SEC, HSEC
C      
C      CALL GETTIM(HOUR,MIN,SEC,HSEC)
C      WRITE(TSTRNG,'(I2.2,A,I2.2,A,I2.2)')HOUR,':',MIN,':',SEC
C     ------------------------
C
C
      END


