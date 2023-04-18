      !a set of parameters used in the program
      MODULE nrtype
      
      INTEGER, PARAMETER :: I4B = selected_int_kind(9) !4 
      INTEGER, PARAMETER :: I2B = selected_int_kind(4) !2
      INTEGER, PARAMETER :: I1B = selected_int_kind(2) !1
      INTEGER, PARAMETER :: SP = kind(1.0) !4
      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(10) !8
      INTEGER, PARAMETER :: SPC = kind( (1.0,1.0) ) !4
      INTEGER, PARAMETER :: DPC = kind( (1.0_dp,1.0_dp) ) !8
      INTEGER, PARAMETER :: LGT = kind( .true. ) !4
      REAL(DP),PARAMETER :: PI=3.141592653589793_dp
      REAL(DP),PARAMETER :: TWOPI=2._dp*PI
      
      COMPLEX(DPC), PARAMETER :: im = (0._dp,1._dp)
      TYPE sprs2_sp
      INTEGER(I4B) :: n,len
      REAL(SP), DIMENSION(:), POINTER :: val
      INTEGER(I4B), DIMENSION(:), POINTER :: irow
      INTEGER(I4B), DIMENSION(:), POINTER :: jcol
      END TYPE sprs2_sp
      TYPE sprs2_dp
      INTEGER(I4B) :: n,len
      REAL(DP), DIMENSION(:), POINTER :: val
      INTEGER(I4B), DIMENSION(:), POINTER :: irow
      INTEGER(I4B), DIMENSION(:), POINTER :: jcol
      END TYPE sprs2_dp
      !  Some important Parameters, to convert to a.u.
      !  - AUTOA  = 1. a.u. in Angstrom
      !  - RYTOEV = 1 Ry in Ev
      !  - EVTOJ  = 1 eV in Joule
      !  - AMTOKG = 1 atomic mass unit ("proton mass") in kg
      !  - BOLKEV = Boltzmanns constant in eV/K
      !  - BOLK   = Boltzmanns constant in Joule/K
      !  - HPLANK = Plank constant in Joule*s
      !  - V_c0   = Velocity of light in vacuum
      !  From PHYISCS (Chinese) ISSN 0379-4148 CN 11-1957/O4
      REAL(DP), PARAMETER :: AUTOA=0.5291772083_dp,RYTOEV=13.60569172_dp
      REAL(DP), PARAMETER::EVTOJ=1.602176462E-19_dp,AMTOKG=1.67262158E-27_dp, &
     &                     BOLKEV=8.617342E-5_dp,BOLK=BOLKEV*EVTOJ, HPLANK=6.62606876E-34_dp, &
     &                     V_c0=2.99792458d8
      END MODULE nrtype



      module reader
      use nrtype 
      implicit none 
      contains 

      subroutine parameter_reader(sigma, output_MODOS, E0, E1, nE, Mulliken, project_substrate_DOS, restart_CM_SM, check_norm, output_MOOP, MOOP_project_basisfunc, MOOP_mol_basis, MOOP_project_substrate)
      implicit none
      real(8), intent(out) :: sigma, E0, E1
      integer, intent(out) :: nE, MOOP_project_basisfunc
      logical, intent(out) :: output_MODOS, Mulliken, project_substrate_DOS, restart_CM_SM, check_norm, output_MOOP, MOOP_mol_basis, MOOP_project_substrate
      integer :: IERR, N, I, IDUM 
      real(8) :: FDUM 
      complex(8) :: CDUM 
      character :: CHARAC 
      logical :: lopen, LDUM
      integer, parameter :: iu5 = 15, iu6 = 6
      character(30), parameter :: INFILE='MODOS_control.in' 
      
      lopen = .False. 
      open(iu5,file=INFILE,status='old')


      !------------------------------------------------------ 
      !Gaussian broadening factor (Default: sigma = 0.1d0) 
      !------------------------------------------------------ 
      sigma = 0.1d0
      call rdatab(lopen, INFILE, iu5, 'sigma', '=', '#', ';', 'F', & 
     &     IDUM, sigma, CDUM, LDUM, CHARAC, N, 1, IERR) 
      if (((IERR/=0) .and. (IERR/=3)).or. & 
     &   ((IERR==0).and.(N<1))) then 
         write(iu6,*)'error reading item ''sigma'' from file ',INFILE 
         goto 150
      end if
      if (sigma < 0) then 
         write(iu6,*)'Error: ''sigma < 0''' 
         stop 
      end if 


      !------------------------------------------------------ 
      !Do Mulliken analysis instead of MODOS analysis (Default: Mulliken = .False.) 
      !------------------------------------------------------ 
      Mulliken = .False.
      call rdatab(lopen, INFILE, iu5, 'Mulliken', '=', '#', ';', 'L', & 
     &     IDUM, FDUM, CDUM, Mulliken, CHARAC, N, 1, IERR) 
      if (((IERR/=0) .and. (IERR/=3)).or. & 
     &   ((IERR==0).and.(N<1))) then 
         write(iu6,*)'error reading item ''Mulliken'' from file ',INFILE 
         goto 150
      end if


      !------------------------------------------------------ 
      !Project DOS onto substrate orbitals meanwhile (Default: project_substrate_DOS = .False.) 
      !------------------------------------------------------ 
      project_substrate_DOS = .False.
      call rdatab(lopen, INFILE, iu5, 'project_substrate_DOS', '=', '#', ';', 'L', & 
     &     IDUM, FDUM, CDUM, project_substrate_DOS, CHARAC, N, 1, IERR) 
      if (((IERR/=0) .and. (IERR/=3)).or. & 
     &   ((IERR==0).and.(N<1))) then 
         write(iu6,*)'error reading item ''project_substrate_DOS'' from file ',INFILE 
         goto 150
      end if


      !------------------------------------------------------ 
      !Read or write restart files which store CM and SM (Default: restart_CM_SM = .False.) 
      !------------------------------------------------------ 
      restart_CM_SM = .False.
      call rdatab(lopen, INFILE, iu5, 'restart_CM_SM', '=', '#', ';', 'L', & 
     &     IDUM, FDUM, CDUM, restart_CM_SM, CHARAC, N, 1, IERR) 
      if (((IERR/=0) .and. (IERR/=3)).or. & 
     &   ((IERR==0).and.(N<1))) then 
         write(iu6,*)'error reading item ''restart_CM_SM'' from file ',INFILE 
         goto 150
      end if


      !------------------------------------------------------ 
      !Check the normalization condition of eigenvector and overlap matrix (Default: check_norm = .False.) 
      !------------------------------------------------------ 
      check_norm = .False.
      call rdatab(lopen, INFILE, iu5, 'check_norm', '=', '#', ';', 'L', & 
     &     IDUM, FDUM, CDUM, check_norm, CHARAC, N, 1, IERR) 
      if (((IERR/=0) .and. (IERR/=3)).or. & 
     &   ((IERR==0).and.(N<1))) then 
         write(iu6,*)'error reading item ''check_norm'' from file ',INFILE 
         goto 150
      end if


      !------------------------------------------------------ 
      !Output MODOS (turn off for large systems) (Default: output_MODOS = .True.) 
      !------------------------------------------------------ 
      output_MODOS = .False.
      call rdatab(lopen, INFILE, iu5, 'output_MODOS', '=', '#', ';', 'L', & 
     &     IDUM, FDUM, CDUM, output_MODOS, CHARAC, N, 1, IERR) 
      if (((IERR/=0) .and. (IERR/=3)).or. & 
     &   ((IERR==0).and.(N<1))) then 
         write(iu6,*)'error reading item ''output_MODOS'' from file ',INFILE 
         goto 150
      end if

      
      !------------------------------------------------------ 
      !Output MOOP (turn off for large systems) (Default: output_MOOP = .True.) 
      !------------------------------------------------------ 
      output_MOOP = .False.
      call rdatab(lopen, INFILE, iu5, 'output_MOOP', '=', '#', ';', 'L', & 
     &     IDUM, FDUM, CDUM, output_MOOP, CHARAC, N, 1, IERR) 
      if (((IERR/=0) .and. (IERR/=3)).or. & 
     &   ((IERR==0).and.(N<1))) then 
         write(iu6,*)'error reading item ''output_MOOP'' from file ',INFILE 
         goto 150
      end if
      

      !------------------------------------------------------ 
      !The number of energy points in MODOS (Default: nE = 1001) 
      !------------------------------------------------------ 
      nE = 1001
      call rdatab(lopen, INFILE, iu5, 'nE', '=', '#', ';', 'I', & 
     &     nE, FDUM, CDUM, LDUM, CHARAC, N, 1, IERR) 
      if (((IERR/=0) .and. (IERR/=3)).or. & 
     &   ((IERR==0).and.(N<1))) then 
         write(iu6,*)'error reading item ''nE'' from file ',INFILE 
         goto 150
      end if
      if (nE < 1) then 
         write(iu6,*)'Error: ''nE < 1''' 
         stop
      end if 
      
      
      !------------------------------------------------------ 
      !Project MOOP onto specific molecular eigenvector (Default: MOOP_mol_basis = False)
      !if MOOP_mol_basis = .False. MOOP is calculate in atomic basis
      !------------------------------------------------------
      MOOP_mol_basis = .False.
      call rdatab(lopen, INFILE, iu5, 'MOOP_mol_basis', '=', '#', ';', 'L', & 
     &     IDUM, FDUM, CDUM, MOOP_mol_basis, CHARAC, N, 1, IERR) 
      if (((IERR/=0) .and. (IERR/=3)).or. & 
     &   ((IERR==0).and.(N<1))) then 
         write(iu6,*)'error reading item ''MOOP_mol_basis'' from file ',INFILE 
         goto 150
      end if
      
      
      !------------------------------------------------------ 
      !Project MOOP onto specific substrate basis state (Default: MOOP_project_substrate = False)
      !if MOOP_project_substrate = .False. MOOP is calculate in atomic basis
      !------------------------------------------------------
      MOOP_project_substrate = .False.
      call rdatab(lopen, INFILE, iu5, 'MOOP_project_substrate', '=', '#', ';', 'L', & 
     &     IDUM, FDUM, CDUM, MOOP_project_substrate, CHARAC, N, 1, IERR) 
      if (((IERR/=0) .and. (IERR/=3)).or. & 
     &   ((IERR==0).and.(N<1))) then 
         write(iu6,*)'error reading item ''MOOP_project_substrate'' from file ',INFILE 
         goto 150
      end if

      
      !------------------------------------------------------ 
      !Project MOOP onto specific basisfunction (Default: MOOP_project_basisfunc = 0)
      !if MOOP_project_basisfunc = 0 no projection
      !------------------------------------------------------ 
      MOOP_project_basisfunc = 0
      call rdatab(lopen, INFILE, iu5, 'MOOP_project_basisfunc', '=', '#', ';', 'I', & 
     &     MOOP_project_basisfunc, FDUM, CDUM, LDUM, CHARAC, N, 1, IERR) 
      if (((IERR/=0) .and. (IERR/=3)).or. & 
     &   ((IERR==0).and.(N<1))) then 
         write(iu6,*)'error reading item ''MOOP_project_basisfunc'' from file ',INFILE 
         goto 150
      end if
      
      
      !------------------------------------------------------ 
      !Initial energy in MODOS (eV, relative to E-Fermi) (Default: E0 = -10.0d0) 
      !------------------------------------------------------ 
      E0 = -10.0d0
      call rdatab(lopen, INFILE, iu5, 'E0', '=', '#', ';', 'F', & 
     &     IDUM, E0, CDUM, LDUM, CHARAC, N, 1, IERR) 
      if (((IERR/=0) .and. (IERR/=3)).or. & 
     &   ((IERR==0).and.(N<1))) then 
         write(iu6,*)'error reading item ''E0'' from file ',INFILE 
         goto 150 
      end if 


      !------------------------------------------------------ 
      !Final energy in MODOS (eV, relative to E-Fermi) (Default: E1 = 5.0d0) 
      !------------------------------------------------------ 
      E1 = 5.0d0
      call rdatab(lopen, INFILE, iu5, 'E1', '=', '#', ';', 'F', & 
     &     IDUM, E1, CDUM, LDUM, CHARAC, N, 1, IERR) 
      if (((IERR/=0) .and. (IERR/=3)).or. & 
     &   ((IERR==0).and.(N<1))) then 
         write(iu6,*)'error reading item ''E1'' from file ',INFILE 
         goto 150 
      end if 


      close(iu5)
      
      RETURN 
      
  150 continue
      write(iu6,151) IERR, N
      
  151 format(' Error code was IERR=',I1,' ... . Found N=',I5,' data.')
      STOP
      
      END SUBROUTINE parameter_reader


      end module reader



      FUNCTION LENGTH(STRING)
      USE nrtype
      IMPLICIT REAL(8) (A-H,O-Z)
      ! Returns the position of the last non-blank character in STRING
      CHARACTER*(*)   STRING
      CHARACTER*256   B8
      CHARACTER*128   B7
      CHARACTER*64    B6
      CHARACTER*32    B5
      INTEGER         LENGTH,LEN,L,I
      SAVE            B5,B6,B7,B8
      DATA            B5 /' '/, B6 /' '/, B7 /' '/, B8 /' '/

      L=LEN(STRING)
      ! Very crude 'scan' for 'order of length' (performance!!) ...:
    6 IF (L.GE.256) THEN
         IF (STRING(L-255:L).EQ.B8) THEN
            L=L-256
            GOTO 7
         ELSE IF (STRING(L-127:L).EQ.B7) THEN
            L=L-128
            GOTO 7
         ELSE IF (STRING(L-63:L).EQ.B6) THEN
            L=L-64
            GOTO 7
         ELSE IF (STRING(L-31:L).EQ.B5) THEN
            L=L-32
            GOTO 7
         ENDIF
      ENDIF
    7 IF (L.GE.128) THEN
         IF (STRING(L-127:L).EQ.B7) THEN
            L=L-128
            GOTO 8
         ELSE IF (STRING(L-63:L).EQ.B6) THEN
            L=L-64
            GOTO 8
         ELSE IF (STRING(L-31:L).EQ.B5) THEN
            L=L-32
            GOTO 8
         ENDIF
      ENDIF
    8 IF (L.GE.64) THEN
         IF (STRING(L-63:L).EQ.B6) THEN
            L=L-64
            GOTO 9
         ELSE IF (STRING(L-31:L).EQ.B5) THEN
            L=L-32
            GOTO 9
         ENDIF
      ENDIF
    9 IF (L.GE.32) THEN
         IF (STRING(L-31:L).EQ.B5) L=L-32
      ENDIF
      ! Here we should have reached some point where either L is much
      ! smaller than LEN(STRING)/very small at all or more general where
      ! L is quite close to the result for length (so that all goes very
      ! quick now ... --- function LENGTH should show high performance).
      LENGTH=0
      DO 10 I=L,1,-1
         IF (STRING(I:I).NE.' ') THEN
            LENGTH=I
            GOTO 11
         ENDIF
   10 CONTINUE
   11 RETURN
      END FUNCTION LENGTH



      FUNCTION NWORDS(STRING)
      USE nrtype
      
      IMPLICIT REAL(8) (A-H,O-Z)
      ! Find the number of blank-delimited words within some string:
      CHARACTER*(*) STRING
      LOGICAL BLANK
      INTEGER NWORDS,L,LENGTH,I
      EXTERNAL LENGTH

      L=LENGTH(STRING)
      NWORDS=0
      BLANK=.TRUE.
      DO 10 I=1,L
         IF (BLANK.AND.(STRING(I:I).NE.' ')) THEN
            BLANK=.FALSE.
            NWORDS=NWORDS+1
         ELSE IF ((.NOT.BLANK).AND.(STRING(I:I).EQ.' ')) THEN
            BLANK=.TRUE.
         ENDIF
   10 CONTINUE
      RETURN
      END FUNCTION NWORDS



      FUNCTION NITEMS(STRING,WORK,EXPAND,TYPE)
      USE nrtype
      IMPLICIT REAL(8) (A-H,O-Z)
      ! Extract how many data items had to be read if one wanted to read
      ! numerical data from a STRING ... . --> It is not simply the number
      ! of words (because such constructs like 120*0.2 -- counting as 120
      ! different data in this example -- are also allowed in FORTRAN!!).
      ! If one wishes one can resolve such constructs using EXPAND=.TRUE.
      ! (i.e. STRING will be changed so that no more such constructs occur
      ! by translating it in a long list of single items -- of course here
      ! some restrictions apply: STRING and WORK must be long enough to
      ! hold the result and the single items may not exceed 255 characters
      ! and it should be noted that this acts also implicitly like STRIP,
      ! i.e. leading blanks are deleted and all items are separated by a
      ! single blank after all ... . And additionally one can also specify
      ! 'type verification' (means entering some given data via TYPE being
      !   - 'L'   logical items only
      !   - 'I'   integer numbers only
      !   - 'F'   floating point items only
      !   - 'C'   complex numbers only
      !   - or any other type suppressing the check ...
      ! will result in a validity-check of the items found, breaking the
      ! searching and counting at the first invalid item found ...).
      CHARACTER*(*) STRING,WORK
      CHARACTER*1   CHECK,TYPE
      CHARACTER*255 BUFFER,FORM,NUMBER,DUMMY
      LOGICAL EXPAND
      INTEGER NITEMS,NWORDS,NOCCUR,LS,LW,L0,L,I,N,IC,IMULT,LENGTH,J,LEN
      INTRINSIC LEN
      EXTERNAL NWORDS,NOCCUR,LENGTH

      LS=LEN(STRING)
      LW=LEN(WORK)
      IF (EXPAND) CALL STRIP(STRING,L,'B')
  100 NITEMS=0
      N=NWORDS(STRING)
      DO 400 I=1,N
      ! scan word by word ...
         CALL SUBWRD(STRING,WORK,I,1)
         IC=NOCCUR(WORK,'*',0)
      ! invalid item here ---> stop counting items and say good bye ... !
         IF (IC.GT.1) GOTO 500
         IF (IC.EQ.0) THEN
            CHECK='Y'
            IF ((TYPE.EQ.'L').OR.(TYPE.EQ.'I').OR.(TYPE.EQ.'F').OR. &
     &          (TYPE.EQ.'C')) CALL CHKTYP(WORK,DUMMY,TYPE,CHECK,FORM)
      ! invalid item here ---> stop counting items and say good bye ... !
            IF (CHECK.EQ.'N') GOTO 500
      ! 'unreadable' item ---> stop counting items and say good bye ... !
            IF ((CHECK.EQ.'U').AND.(TYPE.EQ.'I')) GOTO 500
            NITEMS=NITEMS+1
         ENDIF
         IF (IC.EQ.1) THEN
            CALL PARSE(WORK,BUFFER,NUMBER,'*',0)
            CHECK='Y'
            IF ((TYPE.EQ.'L').OR.(TYPE.EQ.'I').OR.(TYPE.EQ.'F').OR. &
     &          (TYPE.EQ.'C')) CALL CHKTYP(NUMBER,DUMMY,TYPE,CHECK,FORM)
      ! invalid item here ---> stop counting items and say good bye ... !
            IF (CHECK.EQ.'N') GOTO 500
      ! 'unreadable' item ---> stop counting items and say good bye ... !
            IF ((CHECK.EQ.'U').AND.(TYPE.EQ.'I')) GOTO 500
            IF (EXPAND) CALL STRIP(NUMBER,L0,'A')
            CALL STRIP(BUFFER,L,'A')
      ! invalid item here ---> stop counting items and say good bye ... !
            IF (L.LT.1) GOTO 500
      ! 'multiplier'  MUST  be a positive integer number (strictly!) ...
            CALL CHKINT(BUFFER,DUMMY,CHECK,FORM)
      ! invalid item here ---> stop counting items and say good bye ... !
            IF (CHECK.NE.'Y') GOTO 500
      ! hopefully should work without error/end condition here ... ?
            DUMMY='('//FORM(1:253)//')'
            CALL STRIP(DUMMY,IMULT,'A')
            READ(BUFFER,DUMMY) IMULT
      ! invalid item here ---> stop counting items and say good bye ... !
            IF (IMULT.LE.0) GOTO 500
            IF (EXPAND) THEN
      ! invalid item here ---> stop counting items and say good bye ... !
               IF (L0.LT.1) GOTO 500
               IF (I.EQ.1) WORK=' '
               IF (I.GT.1) CALL SUBWRD(STRING,WORK,1,I-1)
               L=LENGTH(WORK)
               DO 200 J=1,IMULT
                  IF (I.EQ.1) WORK=NUMBER
                  IF (I.GT.1) WORK=WORK(1:L)//' '//NUMBER
                  L=L+L0+1
                  IF (I.EQ.1) L=L-1
                  IF ((L.GT.LW).OR.(L.GT.LS)) THEN
                     WRITE(*,'(A)') ' '
                     WRITE(*,'(A)') 'Error LEXLIB routine'// &
     &                              ' ''NITEMS'': expansion fails,'// &
     &                              ' insufficient CHARACTER length.'
      ! have to stop expansion ---> stop counting items and say good bye ... !
                     GOTO 500
                  ENDIF
  200          CONTINUE
               DO 300 J=I+1,N
                  CALL SUBWRD(STRING,BUFFER,J,1)
                  CALL STRIP(BUFFER,L0,'A')
                  IF (L0.LT.1) GOTO 500
                  WORK=WORK(1:L)//' '//BUFFER
                  L=L+L0+1
                  IF ((L.GT.LW).OR.(L.GT.LS)) THEN
                     WRITE(*,'(A)') ' '
                     WRITE(*,'(A)') 'Error LEXLIB routine'// &
     &                              ' ''NITEMS'': expansion fails,'// &
     &                              ' insufficient CHARACTER length.'
      ! have to stop expansion ---> stop counting items and say good bye ... !
                     GOTO 500
                  ENDIF
  300          CONTINUE
      ! sorry, do not know a simpler way to get new correct 'counters' here:
               STRING=WORK
               GOTO 100
            ELSE
               NITEMS=NITEMS+IMULT
            ENDIF
         ENDIF
  400 CONTINUE
  500 CONTINUE
      RETURN
      END



      SUBROUTINE SUBWRD(STRING,WORDS,IBEG,INUM)
      USE nrtype
      IMPLICIT REAL(8) (A-H,O-Z)
      ! Extracts specified words out of a string ...
      CHARACTER*(*) STRING,WORDS
      LOGICAL BLANK
      INTEGER NWORDS,L,LEN,LENGTH,I,IBEG,INUM,ISTART,ISTOP
      EXTERNAL LENGTH

      L=LENGTH(STRING)
      ISTART=0
      ISTOP=L
      NWORDS=0
      BLANK=.TRUE.
      DO 10 I=1,L
         IF (BLANK.AND.(STRING(I:I).NE.' ')) THEN
            BLANK=.FALSE.
            NWORDS=NWORDS+1
            IF (NWORDS.EQ.IBEG) ISTART=I
         ELSE IF ((.NOT.BLANK).AND.(STRING(I:I).EQ.' ')) THEN
            BLANK=.TRUE.
            IF (NWORDS.EQ.(IBEG+INUM-1)) ISTOP=I-1
         ENDIF
   10 CONTINUE
      IF ((ISTART.GT.0).AND.(INUM.GT.0)) THEN
         WORDS=STRING(ISTART:ISTOP)
         IF ((ISTOP-ISTART+1).GT.LEN(WORDS)) THEN
            WRITE(*,'(A)') ' '
            WRITE(*,'(A)') 'Warning LEXLIB routine ''SUBWRD'': '// &
     &             'Output string will be truncated! The output is'
            WRITE(*,'(A,I5,A,I5,A)') 'a string of length ', &
     &             ISTOP-ISTART+1, &
     &             ' characters but ''WORDS'' can only hold ', &
     &             LEN(WORDS),' characters.'
            WRITE(*,'(A)') 'Continuing execution ...'
            WRITE(*,'(A)') ' '
         ENDIF
      ELSE
         WORDS=' '
      ENDIF
      RETURN
      END



      SUBROUTINE STRIP(STRING,L,MODE)
      USE nrtype
      IMPLICIT REAL(8) (A-H,O-Z)
      ! Strips off blanks in STRING according to setting of MODE and returns
      ! the position L of the last non-blank character after all operations.
      ! MODE may be set to:
      !   - 'L' remove all leading blanks only
      !   - 'I' remove all blanks inside STRING, let leading blanks untouched
      !   - 'S' merge all multiple blanks within STRING into one single blank
      !         but leave all leading blanks untouched!
      !   - 'B' remove leading blanks and merge all multiple blanks into one
      !   - 'A' remove all (but really  all!) blanks
      ! all other settings lead to output L=0 (returns a blank string)!
      CHARACTER*(*) STRING
      CHARACTER*1   MODE
      INTEGER       L,LENGTH,L0,FIRST,POS,I
      EXTERNAL      LENGTH

      L0=LENGTH(STRING)

      ! Here stripping off all leading blanks:
      IF ((MODE.EQ.'L').OR.(MODE.EQ.'A').OR.(MODE.EQ.'B')) THEN
         FIRST=L0+1
         DO 10 I=1,L0
            IF (STRING(I:I).NE.' ') THEN
               FIRST=I
               GOTO 20
            ENDIF
  10     CONTINUE
  20     CONTINUE
         IF (FIRST.LE.L0) THEN
            STRING=STRING(FIRST:L0)
            L=L0
            IF (MODE.EQ.'L') L=LENGTH(STRING)
         ELSE
            STRING=' '
            L=0
         ENDIF
      END IF
      ! Here stripping off all blanks inside STRING except for leading blanks
      IF ((MODE.EQ.'I').OR.(MODE.EQ.'A')) THEN
         FIRST=L0+1
         DO 30 I=1,L0
            IF (STRING(I:I).NE.' ') THEN
               FIRST=I
               GOTO 40
            ENDIF
  30     CONTINUE
  40     CONTINUE
         IF (FIRST.LE.L0) THEN
            POS=FIRST+1
            DO 50 I=FIRST+1,L0
               IF (STRING(POS:POS).EQ.' ') THEN
                  IF (POS.LT.L0) THEN
                     STRING=STRING(1:POS-1)//STRING(POS+1:L0)
                  ELSE
                     STRING=STRING(1:POS-1)
                  ENDIF
               ELSE
                  POS=POS+1
               ENDIF
   50       CONTINUE
            L=LENGTH(STRING)
         ELSE
            STRING=' '
            L=0
         ENDIF
      ENDIF
      ! Here merging multiple blanks into a single blank (except leading blanks)
      IF ((MODE.EQ.'S').OR.(MODE.EQ.'B')) THEN
         FIRST=L0+1
         DO 60 I=1,L0
            IF (STRING(I:I).NE.' ') THEN
               FIRST=I
               GOTO 70
            ENDIF
  60     CONTINUE
  70     CONTINUE
         IF (FIRST.LE.L0) THEN
            POS=FIRST+1
            DO 80 I=FIRST+1,L0-1
               IF (STRING(POS:POS+1).EQ.'  ') THEN
                  IF (POS.LT.L0) THEN
                     STRING=STRING(1:POS-1)//STRING(POS+1:L0)
                  ELSE
                     STRING=STRING(1:POS-1)
                  ENDIF
               ELSE
                  POS=POS+1
               ENDIF
   80       CONTINUE
            L=LENGTH(STRING)
         ELSE
            STRING=' '
            L=0
         ENDIF
      ENDIF
      ! Invalid mode --> return blank string
      IF ((MODE.NE.'L').AND.(MODE.NE.'I').AND.(MODE.NE.'S') &
     &                 .AND.(MODE.NE.'A').AND.(MODE.NE.'B')) THEN
         STRING=' '
         L=0
      ENDIF

      RETURN
      END



      SUBROUTINE UPPER(STRING)
      USE nrtype
      IMPLICIT REAL(8) (A-H,O-Z)
      ! uppercase all letters in STRING ...
      CHARACTER*(*) STRING
      CHARACTER*1   ALPHAL(26),ALPHAU(26)
      INTEGER       L,I,J,LENGTH
      EXTERNAL      LENGTH
      SAVE          ALPHAL,ALPHAU
      DATA ALPHAL /'a','b','c','d','e','f','g','h','i','j','k','l','m', &
     &             'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      DATA ALPHAU /'A','B','C','D','E','F','G','H','I','J','K','L','M', &
     &             'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/

      L=LENGTH(STRING)
      DO 30 I=1,L
         DO 10 J=1,26
            IF (STRING(I:I).EQ.ALPHAL(J)) THEN
               STRING(I:I)=ALPHAU(J)
               GOTO 20
            ENDIF
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE

      RETURN
      END



      SUBROUTINE LOWER(STRING)
      USE nrtype
      IMPLICIT REAL(8) (A-H,O-Z)
      ! lowercase all letters in STRING ...
      CHARACTER*(*) STRING
      CHARACTER*1   ALPHAL(26),ALPHAU(26)
      INTEGER       L,I,J,LENGTH
      EXTERNAL      LENGTH
      SAVE          ALPHAL,ALPHAU
      DATA ALPHAL /'a','b','c','d','e','f','g','h','i','j','k','l','m', &
     &             'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      DATA ALPHAU /'A','B','C','D','E','F','G','H','I','J','K','L','M', &
     &             'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/

      L=LENGTH(STRING)
      DO 30 I=1,L
         DO 10 J=1,26
            IF (STRING(I:I).EQ.ALPHAU(J)) THEN
               STRING(I:I)=ALPHAL(J)
               GOTO 20
            ENDIF
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE

      RETURN
      END



      SUBROUTINE PARSE(STRING,BEFORE,AFTER,PATTERN,MODE)
      USE nrtype
      IMPLICIT REAL(8) (A-H,O-Z)
      ! Parses STRING like REXX-interpreter with rule BEFORE 'PATTERN' AFTER;
      ! MODE determines which "length definition" shall be used for PATTERN
      ! (>=0: take length from LEN, means treat also trailing blanks or <0:
      ! take length from LENGTH, means ignore all trailing blanks ...).
      CHARACTER*(*) STRING,BEFORE,AFTER,PATTERN
      INTEGER       MODE,LENGTH,LS,LP,I
      EXTERNAL      LENGTH

      BEFORE=STRING
      AFTER=' '
      IF (MODE.GE.0) THEN
         LP=LEN(PATTERN)
      ELSE
      ! additional remark: blank string shall act like "single blank" ...
         LP=MAX(LENGTH(PATTERN),1)
      ENDIF
      LS=LENGTH(STRING)
      I=INDEX(STRING(1:LS),PATTERN(1:LP))
      IF (I.EQ.0) RETURN
      IF ((I.GT.1).AND.(I+LP.LE.LS)) THEN
         BEFORE=STRING(1:I-1)
         AFTER=STRING(I+LP:LS)
      ELSE IF ((I.EQ.1).AND.(I+LP.LE.LS)) THEN
         BEFORE=' '
         AFTER=STRING(I+LP:LS)
      ELSE IF ((I.EQ.1).AND.(I+LP.GT.LS)) THEN
         BEFORE=' '
         AFTER=' '
      ELSE IF (I+LP.GT.LS) THEN
         BEFORE=STRING(1:I-1)
         AFTER=' '
      ENDIF

      RETURN
      END



      FUNCTION NOCCUR(STRING,PATTERN,MODE)
      USE nrtype
      IMPLICIT REAL(8) (A-H,O-Z)
      ! Tells how often PATTERN occurs within STRING!
      ! MODE determines which "length definition" shall be used for PATTERN
      ! (>=0: take length from LEN, means treat also trailing blanks or <0:
      ! take length from LENGTH, means ignore all trailing blanks ...).
      CHARACTER*(*) STRING,PATTERN
      INTEGER       NOCCUR,LENGTH,LS,LP,I,LAST,MODE
      EXTERNAL      LENGTH

      NOCCUR=0
      LAST=1
      LS=LENGTH(STRING)
      IF (MODE.GE.0) THEN
         LP=LEN(PATTERN)
      ELSE
      ! again blank strings should be interpreted as 'single blank':
         LP=MAX(LENGTH(PATTERN),1)
      ENDIF
   10 CONTINUE
      I=INDEX(STRING(LAST:LS),PATTERN(1:LP))
      IF (I.EQ.0) RETURN
      NOCCUR=NOCCUR+1
      LAST=LAST+I+LP-1
      IF ((LAST+LP).GT.LS) RETURN
      GOTO 10

      RETURN
      END



      FUNCTION INDEXN(STRING,PATTERN,NTH)
      USE nrtype
      IMPLICIT REAL(8) (A-H,O-Z)
      ! Get the starting position of PATTERN within STRING for the NTHth
      ! occurence (it is some "generalized INDEX-function" ...):
      CHARACTER*(*) STRING,PATTERN
      INTEGER       INDEXN,NTH,LENGTH,LS,LP,I,LAST,NOCCUR
      EXTERNAL      LENGTH

      INDEXN=0
      NOCCUR=0
      LAST=1
      LS=LENGTH(STRING)
      IF (NTH.GE.0) THEN
         LP=LEN(PATTERN)
      ELSE
         LP=LENGTH(PATTERN)
      ENDIF
   10 CONTINUE
      I=INDEX(STRING(LAST:LS),PATTERN(1:LP))
      IF (I.EQ.0) RETURN
      NOCCUR=NOCCUR+1
      IF (NOCCUR.EQ.ABS(NTH)) INDEXN=LAST+I-1
      LAST=LAST+I+LP-1
      IF ((NOCCUR.EQ.ABS(NTH)).OR.((LAST+LP).GT.LS)) RETURN
      GOTO 10

      RETURN
      END



      SUBROUTINE REPLAC(STRING,OLDPAT,NEWPAT,L,MODE)
      USE nrtype
      IMPLICIT REAL(8) (A-H,O-Z)
      ! Replace some pattern OLDPAT by pattern NEWPAT within STRING according
      ! to the setting of MODE which may take following values:
      !   -  0  global replacement (for all --- but really all --- occurences)
      !   -  >0 replacement for the MODEth occurence only
      !   -  <0 replacement for the first ABS(MODE) occurences
      ! on output L returns the position of the last non-blank character
      ! in STRING after all replacements, on input it controls which 'length'
      ! for OLDPAT/NEWPAT should be taken (L>0: length defined by LEN, means
      ! treat also trailing blanks or L<0: length defined by LENGTH, means all
      ! trailing blanks will be ignored ... . L=0 has the special meaning: ``do
      ! the same as for L>0 for OLDPAT, but the same as L<0 for NEWPAT'' (will
      ! become important for the special case when one wants exact replacement
      ! of OLDPAT including trailing blanks by a 'null string' ---> NEWPAT=' '
      ! and L=0 should be taken if one wishes to do that ...).
      CHARACTER*(*) STRING,OLDPAT,NEWPAT
      INTEGER      L,MODE,LENGTH,LOLD,LNEW,N,N1,N2,NOCCUR,I,INDEXN,L0,J
      INTEGER      IU,NXTFRU
      EXTERNAL     LENGTH,NOCCUR,INDEXN,NXTFRU

      L0=LEN(STRING)
      IF (L.GT.0) THEN
         J=1
         LOLD=LEN(OLDPAT)
         LNEW=LEN(NEWPAT)
      ELSE IF (L.EQ.0) THEN
         J=1
         LOLD=LEN(OLDPAT)
         LNEW=LENGTH(NEWPAT)
      ELSE
         J=-1
         LOLD=LENGTH(OLDPAT)
         LNEW=LENGTH(NEWPAT)
      ENDIF
      ! Do not replace 'null strings' (how to do??)!
      IF (LOLD.EQ.0) RETURN
      L=LENGTH(STRING)
      ! Nothing what could be replaced ...
      IF (L.EQ.0) RETURN
      N=NOCCUR(STRING,OLDPAT(1:LOLD),J)
      IF (N.EQ.0) RETURN
      IF (MODE.LT.0) THEN
         N1=1
         N2=ABS(MODE)
         IF (N2.GT.N) N2=N
      ELSE IF (MODE.GT.0) THEN
         IF (MODE.GT.N) RETURN
         N1=MODE
         N2=MODE
      ELSE
         N1=1
         N2=N
      ENDIF
      DO 10 N=N1,N2
         I=INDEXN(STRING,OLDPAT(1:LOLD),J*N1)
         IF ((LNEW.GT.0).AND.(LNEW.LE.LOLD)) THEN
            IF ((I.GT.1).AND.((I+LOLD).LE.L0)) THEN
               STRING=STRING(1:I-1)//NEWPAT(1:LNEW)//STRING(I+LOLD:L0)
            ELSE IF (I.EQ.1) THEN
               STRING=NEWPAT(1:LNEW)//STRING(1+LOLD:L0)
            ELSE IF ((I+LOLD).GT.L0) THEN
               STRING=STRING(1:I-1)//NEWPAT(1:LNEW)
            ENDIF
            L=L+LNEW-LOLD
         ELSE IF ((LNEW.GT.0).AND.(LNEW.GT.LOLD)) THEN
      ! Here we run into trouble due to the order in which the expression and
      ! the assignment are done: some intermediate partially replaced string
      ! will be used for the further evaluation of the expression!! So we have
      ! to take this into account very carefully (if LNEW>LOLD) ... !! Shit!!!
      ! We do not want to waste space for another temporary character variable
      ! so we do it the ugly way: via external I/O (on a scratch file). Sorry!
            IU=NXTFRU()
            OPEN(IU,STATUS='SCRATCH')
            IF ((I.GT.1).AND.((I+LOLD).LE.L0)) THEN
               WRITE(IU,'(A,A,A)') &
     &                 STRING(1:I-1),NEWPAT(1:LNEW),STRING(I+LOLD:L0)
            ELSE IF (I.EQ.1) THEN
               WRITE(IU,'(A,A)') NEWPAT(1:LNEW),STRING(1+LOLD:L0)
            ELSE IF ((I+LOLD).GT.L0) THEN
               WRITE(IU,'(A,A)') STRING(1:I-1),NEWPAT(1:LNEW)
            ENDIF
            REWIND IU
            READ(IU,'(A)') STRING
            CLOSE(IU)
            L=L+LNEW-LOLD
            IF (L.GT.L0) THEN
               WRITE(*,'(A)') ' '
               WRITE(*,'(A)') 'Warning LEXLIB routine ''REPLAC'': '// &
     &                'Output string will be truncated! New string is'
               WRITE(*,'(A)') 'longer than old string and the length'// &
     &                    ' of variable ''STRING'' is too short! Length'
               WRITE(*,'(A,I6,A,I5,A)') 'of new string is ',L, &
     &                ' characters but ''STRING'' can only hold ',L0, &
     &                ' characters.'
               WRITE(*,'(A)') 'Continuing execution ...'
               WRITE(*,'(A)') ' '
            ENDIF
         ELSE
      ! special code for 'null string' replacement ...
            IF ((I.GT.1).AND.((I+LOLD).LE.L0)) THEN
               STRING=STRING(1:I-1)//STRING(I+LOLD:L0)
            ELSE IF (I.EQ.1) THEN
               STRING=STRING(1+LOLD:L0)
            ELSE IF ((I+LOLD).GT.L0) THEN
               STRING=STRING(1:I-1)
            ENDIF
            L=L-LOLD
         ENDIF
   10 CONTINUE
      L=MIN(L,L0)

      RETURN
      END



      SUBROUTINE CHKTYP(STRING,WORK,MATCH,TYPE,FORM)
      USE nrtype
      IMPLICIT real(8) (A-H,O-Z)
      ! Try to check the validity of the data type of what is contained in
      ! STRING ... . The type to be tested must be given in MATCH and must be
      !   - 'I'  string should contain a valid FORTRAN-Integer
      !   - 'F'  string should contain a valid FORTRAN-Float (any format)
      !   - 'C'  string should contain a FORTRAN-Complex number (any format)
      !   - 'L'  string should contain a valid FORTRAN-Logical (any format)
      !   - 'A'  string contains only alphanumeric characters [0-9,A-Z,a-z]
      !   - 'U'  uppercase (alphanumeric) string
      !   - 'l'  lowercase (alphanumeric) string
      !   - 'H'  string would be a valid hexadecimal number (chars 0-9,A-F)
      !   - 'O'  string would be a valid octal number (characters 0-7 only)
      !   - 'B'  string would be a valid binary number (only 0's and 1's)
      !   - 'N'  string is a 'null string' (empty string)
      !   - any other (invalid) inputs will behave as 'test was negative'!
      ! The result is returned in TYPE ('N' means 'no' = false and all other
      ! means 'yes' = true [usually 'Y' returned, sometimes other values ...]
      ! In FORM a format string is returned being needed to read from STRING
      ! GENERAL WARNING: *all* blanks are ignored --> we test only 'one word'
      CHARACTER*(*) STRING,WORK,FORM
      CHARACTER*1   MATCH,TYPE,CH
      CHARACTER*15  FORM1,FORM2
      CHARACTER*255 PART1,PART2
      LOGICAL       LTEST,LPURE
      INTEGER       LENGTH,LEN,NOCCUR,L,I,J
      INTRINSIC     LEN
      EXTERNAL      LENGTH,NOCCUR

      TYPE='N'
      FORM=' '
      LTEST=.TRUE.
      LPURE=.TRUE.
      WORK=STRING
      CALL STRIP(WORK,L,'A')
      IF (MATCH.EQ.'I') THEN
         CALL CHKINT(STRING,WORK,TYPE,FORM)
         RETURN
      ELSE IF (MATCH.EQ.'F') THEN
         CALL CHKFLT(STRING,WORK,TYPE,FORM)
         RETURN
      ELSE IF (MATCH.EQ.'C') THEN
      ! a complex number must be of type '(' float/int ',' float/int ')'
         IF ((NOCCUR(WORK,',',0).NE.1).OR.(WORK(1:1).NE.'(') &
     &          .OR.(WORK(L:L).NE.')').OR.(L.LE.4)) LTEST=.FALSE.
      ! all seems to be okay until here ...
         IF (LTEST) THEN
            CALL PARSE(WORK(2:L-1),PART1,PART2,',',0)
      ! part1 is float (or integer -- does not matter ...)?
            CALL CHKFLT(PART1,WORK,CH,FORM1)
            IF (CH.EQ.'N') LTEST=.FALSE.
            IF (LTEST) I=LENGTH(FORM1)
         ENDIF
      ! if part1 was okay, then we have still to test part2 ...
         IF (LTEST) THEN
      ! part2 is float (or integer -- does not matter ...)?
            CALL CHKFLT(PART2,WORK,CH,FORM2)
            IF (CH.EQ.'N') LTEST=.FALSE.
            IF (LTEST) J=LENGTH(FORM2)
         ENDIF
         IF (LTEST) THEN
            TYPE='Y'
      ! format for this input ...
            IF ((I.NE.0).AND.(J.NE.0).AND.(LEN(FORM).GE.(I+J+10))) &
     &                 FORM='1X,'//FORM1(1:I)//',1X,'//FORM2(1:J)//',1X'
         ENDIF
         RETURN
      ELSE IF (MATCH.EQ.'L') THEN
      ! that is very easy (according to FORTRAN rules it must start with
      ! 'T', 'F' or '.T', '.F' --- and then may follow what you want:
         IF ((WORK(1:1).EQ.'T').OR.(WORK(1:1).EQ.'F').OR. &
     &       (WORK(1:2).EQ.'.T').OR.(WORK(1:2).EQ.'.F')) TYPE='Y'
         IF ((TYPE.EQ.'Y').AND.(LEN(FORM).GE.6)) THEN
            WRITE(FORM,'(A1,I5)') 'L',L
            CALL STRIP(FORM,I,'A')
         ENDIF
         RETURN
      ELSE IF (MATCH.EQ.'A') THEN
      ! Check for 'empty string' ...
         IF (L.LE.0) THEN
            LTEST=.FALSE.
            GOTO 101
         ENDIF
      ! Case is not relevant ...
         CALL UPPER(WORK)
      ! Well, check ... :
         DO 100 I=1,L
            CH=WORK(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND. &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND. &
     &          (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND. &
     &          (CH.NE.'9').AND.(CH.NE.'A').AND.(CH.NE.'B').AND. &
     &          (CH.NE.'C').AND.(CH.NE.'D').AND.(CH.NE.'E').AND. &
     &          (CH.NE.'F').AND.(CH.NE.'G').AND.(CH.NE.'H').AND. &
     &          (CH.NE.'I').AND.(CH.NE.'J').AND.(CH.NE.'K').AND. &
     &          (CH.NE.'L').AND.(CH.NE.'M').AND.(CH.NE.'N').AND. &
     &          (CH.NE.'O').AND.(CH.NE.'P').AND.(CH.NE.'Q').AND. &
     &          (CH.NE.'R').AND.(CH.NE.'S').AND.(CH.NE.'T').AND. &
     &          (CH.NE.'U').AND.(CH.NE.'V').AND.(CH.NE.'W').AND. &
     &          (CH.NE.'X').AND.(CH.NE.'Y').AND.(CH.NE.'Z')) &
     &                                                    LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 101
            LPURE=LPURE.AND.(CH.NE.'0').AND.(CH.NE.'1').AND. &
     &                      (CH.NE.'2').AND.(CH.NE.'3').AND. &
     &                      (CH.NE.'4').AND.(CH.NE.'5').AND. &
     &                      (CH.NE.'6').AND.(CH.NE.'7').AND. &
     &                      (CH.NE.'8').AND.(CH.NE.'9')
  100    CONTINUE
  101    IF (LTEST) TYPE='Y'
      ! 'pure string' [no chars 0-9] shall return TYPE='P' instead of TYPE='Y'
         IF (LTEST.AND.LPURE) TYPE='P'
      ELSE IF (MATCH.EQ.'U') THEN
      ! Check for 'empty string' ...
         IF (L.LE.0) THEN
            LTEST=.FALSE.
            GOTO 201
         ENDIF
      ! Well, check ... :
         DO 200 I=1,L
            CH=WORK(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND. &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND. &
     &          (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND. &
     &          (CH.NE.'9').AND.(CH.NE.'A').AND.(CH.NE.'B').AND. &
     &          (CH.NE.'C').AND.(CH.NE.'D').AND.(CH.NE.'E').AND. &
     &          (CH.NE.'F').AND.(CH.NE.'G').AND.(CH.NE.'H').AND. &
     &          (CH.NE.'I').AND.(CH.NE.'J').AND.(CH.NE.'K').AND. &
     &          (CH.NE.'L').AND.(CH.NE.'M').AND.(CH.NE.'N').AND. &
     &          (CH.NE.'O').AND.(CH.NE.'P').AND.(CH.NE.'Q').AND. &
     &          (CH.NE.'R').AND.(CH.NE.'S').AND.(CH.NE.'T').AND. &
     &          (CH.NE.'U').AND.(CH.NE.'V').AND.(CH.NE.'W').AND. &
     &          (CH.NE.'X').AND.(CH.NE.'Y').AND.(CH.NE.'Z')) &
     &                                                    LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 201
            LPURE=LPURE.AND.(CH.NE.'0').AND.(CH.NE.'1').AND. &
     &                      (CH.NE.'2').AND.(CH.NE.'3').AND. &
     &                      (CH.NE.'4').AND.(CH.NE.'5').AND. &
     &                      (CH.NE.'6').AND.(CH.NE.'7').AND. &
     &                      (CH.NE.'8').AND.(CH.NE.'9')
  200    CONTINUE
  201    IF (LTEST) TYPE='Y'
      ! 'pure string' [no chars 0-9] shall return TYPE='P' instead of TYPE='Y'
         IF (LTEST.AND.LPURE) TYPE='P'
      ELSE IF (MATCH.EQ.'l') THEN
      ! Check for 'empty string' ...
         IF (L.LE.0) THEN
            LTEST=.FALSE.
            GOTO 301
         ENDIF
      ! Well, check ... :
         DO 300 I=1,L
            CH=WORK(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND. &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND. &
     &          (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND. &
     &          (CH.NE.'9').AND.(CH.NE.'a').AND.(CH.NE.'b').AND. &
     &          (CH.NE.'c').AND.(CH.NE.'d').AND.(CH.NE.'e').AND. &
     &          (CH.NE.'f').AND.(CH.NE.'g').AND.(CH.NE.'h').AND. &
     &          (CH.NE.'i').AND.(CH.NE.'j').AND.(CH.NE.'k').AND. &
     &          (CH.NE.'l').AND.(CH.NE.'m').AND.(CH.NE.'n').AND. &
     &          (CH.NE.'o').AND.(CH.NE.'p').AND.(CH.NE.'q').AND. &
     &          (CH.NE.'r').AND.(CH.NE.'s').AND.(CH.NE.'t').AND. &
     &          (CH.NE.'u').AND.(CH.NE.'v').AND.(CH.NE.'w').AND. &
     &          (CH.NE.'x').AND.(CH.NE.'y').AND.(CH.NE.'z')) &
     &                                                    LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 301
            LPURE=LPURE.AND.(CH.NE.'0').AND.(CH.NE.'1').AND. &
     &                      (CH.NE.'2').AND.(CH.NE.'3').AND. &
     &                      (CH.NE.'4').AND.(CH.NE.'5').AND. &
     &                      (CH.NE.'6').AND.(CH.NE.'7').AND. &
     &                      (CH.NE.'8').AND.(CH.NE.'9')
  300    CONTINUE
  301    IF (LTEST) TYPE='Y'
      ! 'pure string' [no chars 0-9] shall return TYPE='P' instead of TYPE='Y'
         IF (LTEST.AND.LPURE) TYPE='P'
      ELSE IF (MATCH.EQ.'H') THEN
      ! Check for 'empty string' ...
         IF (L.LE.0) THEN
            LTEST=.FALSE.
            GOTO 401
         ENDIF
      ! Case shall not be of interest ...
         CALL UPPER(WORK)
      ! Well, check ... :
         DO 400 I=1,L
            CH=WORK(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND. &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND. &
     &          (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND. &
     &          (CH.NE.'9').AND.(CH.NE.'A').AND.(CH.NE.'B').AND. &
     &          (CH.NE.'C').AND.(CH.NE.'D').AND.(CH.NE.'E').AND. &
     &          (CH.NE.'F')) LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 401
  400    CONTINUE
  401    IF (LTEST) TYPE='Y'
      ELSE IF (MATCH.EQ.'O') THEN
      ! Check for 'empty string' ...
         IF (L.LE.0) THEN
            LTEST=.FALSE.
            GOTO 501
         ENDIF
      ! Well, check ... :
         DO 500 I=1,L
            CH=WORK(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND. &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND. &
     &          (CH.NE.'6').AND.(CH.NE.'7')) LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 501
  500    CONTINUE
  501    IF (LTEST) TYPE='Y'
      ELSE IF (MATCH.EQ.'B') THEN
      ! Check for 'empty string' ...
         IF (L.LE.0) THEN
            LTEST=.FALSE.
            GOTO 601
         ENDIF
      ! Well, check ... :
         DO 600 I=1,L
            CH=WORK(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1')) LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 601
  600    CONTINUE
  601    IF (LTEST) TYPE='Y'
      ELSE IF (MATCH.EQ.'N') THEN
      ! that is too easy ...
         IF (L.LE.0) TYPE='Y'
      ENDIF
      ! all what ends up here were checks for things which could only be read
      ! as strings (A-format for FORTRAN read ...), give here the format ...
      IF ((TYPE.NE.'N').AND.(LEN(FORM).GE.6).AND.(MATCH.NE.'N')) THEN
         WRITE(FORM,'(A1,I5)') 'A',L
         CALL STRIP(FORM,I,'A')
      ENDIF

      RETURN
      END



      SUBROUTINE DATTYP(STRING,WORK,TYPE,FORM)
      USE nrtype
      IMPLICIT real(8) (A-H,O-Z)
      ! Try to find out the data type of what is contained in STRING ... .
      ! Of course we do not want to distinguish too may special cases here!
      ! The result will be returned in TYPE which could take the values:
      !   - 'G'  'general string' containing any arbitrary things ...
      !   - 'A'  alphanumeric string (contains only characters [0-9,A-Z,a-z])
      !   - 'F'  valid floating point number
      !   - 'C'  valid complex number
      !   - 'I'  valid integer number (would also match type 'F')
      !   - 'L'  valid logical value (would also match type 'A')
      !   - 'N'  'null string' (empty string)
      ! Additionally we try to return an appropriate format string in FORM
      ! GENERAL WARNING: *all* blanks are ignored --> we test only 'one word'
      ! Let us start ...
      CHARACTER*(*) STRING,WORK,FORM
      CHARACTER*1   TYPE,MATCH
      INTEGER       LENGTH,L,LEN
      INTRINSIC     LEN
      EXTERNAL      LENGTH

      ! Now check possible type for possible type -- beginning with the most
      ! special possibilities and ending up with more and more general types:
      TYPE='N'
      ! 'null string' ...
      CALL CHKTYP(STRING,WORK,TYPE,MATCH,FORM)
      IF (MATCH.NE.'N') RETURN
      TYPE='L'
      ! logical value ... (string beginning with 'F','T' or '.F','.T')
      CALL CHKTYP(STRING,WORK,TYPE,MATCH,FORM)
      IF (MATCH.NE.'N') RETURN
      TYPE='I'
      ! integer ... (strict form only ...)
      CALL CHKTYP(STRING,WORK,TYPE,MATCH,FORM)
      IF (MATCH.EQ.'Y') RETURN
      TYPE='C'
      ! complex number ...
      CALL CHKTYP(STRING,WORK,TYPE,MATCH,FORM)
      IF (MATCH.NE.'N') RETURN
      TYPE='F'
      ! floating point number ...
      CALL CHKTYP(STRING,WORK,TYPE,MATCH,FORM)
      IF (MATCH.NE.'N') RETURN
      TYPE='A'
      ! alphanumeric string ...
      CALL CHKTYP(STRING,WORK,TYPE,MATCH,FORM)
      IF (MATCH.NE.'N') RETURN
      TYPE='G'
      ! okay, arriving here means 'general string' ...
      IF (LEN(FORM).GE.6) THEN
         WRITE(FORM,'(A1,I5)') 'A',LENGTH(STRING)
         CALL STRIP(FORM,L,'A')
      ENDIF

      RETURN
      END



      SUBROUTINE CHKINT(STRING,WORK,TYPE,FORM)
      USE nrtype
      IMPLICIT real(8) (A-H,O-Z)
      ! Supplementary routine for checking validity of type INTEGER, the
      ! answer is returned in TYPE ('Y' or 'N'), FORM is a format string
      ! for this integer number (if it is one ...) for FORTRAN-reading ...
      CHARACTER*(*) STRING,WORK,FORM
      CHARACTER*1   TYPE,CH
      CHARACTER*16  FTEST,FTEMP
      CHARACTER*255 PART1,PART2
      LOGICAL       LTEST
      INTEGER       L,I,J,K,LENGTH,NOCCUR,INDEXN,LEN
      REAL          RTEST,MAXINT,TINY
      INTRINSIC     LEN
      EXTERNAL      LENGTH,NOCCUR,INDEXN
      ! TINY is the machine tolerance and MAXINT the largest possible integer
      ! number (here for 32-bit, generally 2**(bit-1) - 1) --> maybe customize
      PARAMETER (TINY=1.E-10_dp,MAXINT=2147483647.E0_dp)

      ! default values ...
      TYPE='N'
      FORM=' '
      LTEST=.TRUE.
      WORK=STRING
      CALL STRIP(WORK,L,'A')
      ! First character may be a sign (+/-) or a number ...
      CH=WORK(1:1)
      ! if it is a sign something non-blank MUST follow ...
      IF (((CH.EQ.'+').OR.(CH.EQ.'-')).AND.(L.LE.1)) LTEST=.FALSE.
      IF ((CH.NE.'+').AND.(CH.NE.'-').AND.(CH.NE.'0').AND. &
     &    (CH.NE.'1').AND.(CH.NE.'2').AND.(CH.NE.'3').AND. &
     &    (CH.NE.'4').AND.(CH.NE.'5').AND.(CH.NE.'6').AND. &
     &    (CH.NE.'7').AND.(CH.NE.'8').AND.(CH.NE.'9')) LTEST=.FALSE.
      IF (.NOT.LTEST) GOTO 101
      DO 100 I=2,L
      ! Now only numbers may follow ...
         CH=WORK(I:I)
         IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND. &
     &       (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND. &
     &       (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND. &
     &       (CH.NE.'9')) LTEST=.FALSE.
         IF (.NOT.LTEST) GOTO 101
  100 CONTINUE
  101 CONTINUE
      ! now last check! is it within the allowed range for integer numbers???
      IF (LTEST) THEN
         WRITE(FTEST,'(A1,I5,A2)') 'F',L,'.0'
         FTEMP='('//FTEST(1:14)//')'
         CALL STRIP(FTEMP,I,'A')
         READ(WORK,FTEMP) RTEST
      ! sorry, number too large!
         IF (ABS(RTEST).GT.MAXINT) LTEST=.FALSE.
      ENDIF
      ! test was successful ...
      IF (LTEST) THEN
         TYPE='Y'
         IF (LEN(FORM).GE.6) THEN
            WRITE(FORM,'(A1,I5)') 'I',L
            CALL STRIP(FORM,I,'A')
         ENDIF
      ELSE
      ! maybe (and that is also some special case we might allow) it is a
      ! floating point number having an integer value, check it and if yes
      ! return TYPE='F' or 'E' (for 'F'/'E'-format ...) instead of TYPE='Y'
      ! but: allow only some maximum number MAXINT ("valid integer range")
         CALL CHKFLT(STRING,WORK,CH,FTEST)
         IF (CH.NE.'N') THEN
      ! all right it is a simple floating point number ...
            FTEMP='('//FTEST(1:14)//')'
            CALL STRIP(FTEMP,I,'A')
            READ(STRING,FTEMP) RTEST
            RTEST=ABS(RTEST)
            IF (CH.EQ.'Y') THEN
      ! all right it is a simple floating point number and not too large?
               I=LENGTH(FTEST)
               IF ((FTEST(1:1).EQ.'F').AND.(FTEST(I-1:I).EQ.'.0').AND. &
     &                                         (RTEST.LT.MAXINT)) THEN
      ! it must be Fxxx.0-format, then it is okay and quite simple ...
                  TYPE='F'
                  IF (LEN(FORM).GE.9) THEN
                     WRITE(FORM,'(A1,I5,A3)') 'I',I-1,',1X'
                     CALL STRIP(FORM,I,'A')
                  ENDIF
               ENDIF
            ELSE
      ! E-format ...
               IF (((ABS(RTEST-NINT(RTEST))/RTEST).LT.TINY).AND. &
     &                                         (RTEST.LT.MAXINT)) THEN
      ! bingo? -- it represents some (not too large) integer number ...
                  WORK=STRING
                  CALL STRIP(WORK,L,'A')
                  IF ((CH.EQ.'E').OR.(CH.EQ.'D').OR.(CH.EQ.'Q')) &
     &                               CALL PARSE(WORK,PART1,PART2,CH,0)
                  IF (CH.EQ.'S') THEN
                     I=NOCCUR(WORK,'-',0)
                     J=NOCCUR(WORK,'+',0)
                     K=MAX(INDEXN(WORK,'-',I),INDEXN(WORK,'+',J))
                     PART1=WORK(1:K-1)
                     PART2=WORK(K:L)
                  ENDIF
                  CALL STRIP(PART1,I,'A')
                  CALL STRIP(PART2,J,'A')
      ! now some problem ... --> might get trouble finding some I-format; can
      ! only do it with constructs like yyyyyy(.)E(+/-)0, rest is impossible
      ! without changing string ... (so rest returns FORM=' ' and TYPE='U' !)
                  READ(PART2,'(I6)') K
                  IF (K.EQ.0) THEN
      ! exponent is zero, if it is now a Exxx.0-format it is okay and simple
                     K=LENGTH(FTEST)
                     IF (FTEST(K-1:K).EQ.'.0') THEN
                        TYPE='E'
                        IF (LEN(FORM).GE.13) THEN
                           I=I-NOCCUR(PART1,'.',0)
                           IF (CH.NE.'S') J=J+1
                           J=J+NOCCUR(PART1,'.',0)
                           WRITE(FORM,'(A1,I5,A1,I5,A1)')'I',I,',',J,'X'
                           CALL STRIP(FORM,J,'A')
                        ENDIF
                     ELSE
      ! sorry no I-format available for this number ... (should never occur?)
                        TYPE='U'
                     ENDIF
                  ELSE
      ! sorry no I-format available for this number ...
                     TYPE='U'
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END



      SUBROUTINE CHKFLT(STRING,WORK,TYPE,FORM)
      USE nrtype
      IMPLICIT real(8) (A-H,O-Z)
      ! Supplementary routine for checking validity of floating point number,
      ! the answer is returned in TYPE ('Y' or 'N'), FORM is a format string
      ! for this floating point number (if it is one ...) for a FORTRAN-read
      CHARACTER*(*) STRING,WORK,FORM
      CHARACTER*1   TYPE,CH
      LOGICAL       LTEST,LDOT
      INTEGER       L,I,J,LEN
      INTRINSIC     LEN

      ! default values ...
      TYPE='N'
      FORM=' '
      LTEST=.TRUE.
      WORK=STRING
      CALL STRIP(WORK,L,'A')
      LDOT=.FALSE.
      ! First character may be a sign (+/-) a dot (.) or a number ...
      CH=WORK(1:1)
      ! if it is a sign or a dot something non-blank MUST follow ...
      IF (((CH.EQ.'+').OR.(CH.EQ.'-').OR.(CH.EQ.'.')) &
     &                              .AND.(L.LE.1)) LTEST=.FALSE.
      IF ((CH.NE.'+').AND.(CH.NE.'-').AND.(CH.NE.'0').AND. &
     &    (CH.NE.'1').AND.(CH.NE.'2').AND.(CH.NE.'3').AND. &
     &    (CH.NE.'4').AND.(CH.NE.'5').AND.(CH.NE.'6').AND. &
     &    (CH.NE.'7').AND.(CH.NE.'8').AND.(CH.NE.'9').AND. &
     &    (CH.NE.'.')) LTEST=.FALSE.
      ! only one single dot may be there, so remind if we had some already
      IF (CH.EQ.'.') THEN
         LDOT=.TRUE.
      ! the format would be the following (if the rest is correct ...)
         IF (LEN(FORM).GE.12) THEN
            WRITE(FORM,'(A1,I5,A1,I5)') 'F',L,'.',L-1
            CALL STRIP(FORM,J,'A')
         ENDIF
      ENDIF
      IF (.NOT.LTEST) GOTO 101
      DO 100 I=2,L
      ! Now only numbers may follow and somewhere a dot ...
         CH=WORK(I:I)
         IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND. &
     &       (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND. &
     &       (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND. &
     &       (CH.NE.'9').AND.(CH.NE.'.')) LTEST=.FALSE.
      ! only one single dot may be there and remind if we had some already
         IF ((CH.EQ.'.').AND.LDOT) LTEST=.FALSE.
         IF (CH.EQ.'.') THEN
            LDOT=.TRUE.
      ! the format would be the following (if the rest is correct ...)
            IF (LEN(FORM).GE.12) THEN
               WRITE(FORM,'(A1,I5,A1,I5)') 'F',L,'.',L-I
               CALL STRIP(FORM,J,'A')
            ENDIF
         ENDIF
         IF (.NOT.LTEST) GOTO 101
  100 CONTINUE
      ! okay, LTEST should be true and LDOT too (otherwise it is an INTEGER!)
  101 IF (LTEST.AND.LDOT) TYPE='Y'
      ! if LTEST is true but LDOT is false we got an integer -- of course an
      ! integer is also valid for assignments to real values, so in this case
      ! we will not answer 'N' (no), but also not 'Y' (yes) --> answer is 'I'
      IF (LTEST.AND.(.NOT.LDOT)) THEN
         TYPE='I'
         IF (LEN(FORM).GE.12) THEN
            WRITE(FORM,'(A1,I5,A1,I5)') 'F',L,'.',0
            CALL STRIP(FORM,J,'A')
         ENDIF
      ENDIF
      IF (TYPE.EQ.'N') THEN
      ! if it is not a floating point number in direct F-format, it could be
      ! a number in 'E'-,'D'- or 'Q'-format ... --> check this here, it is of
      ! course also a correct and valid definition of a floating point number
         CALL CHKEDQ(STRING,WORK,'E',CH,FORM)
      ! return type='E' instead of type='Y'
         IF (CH.EQ.'Y') TYPE='E'
      ! or was it this ugly strange format without any 'E'/'D'/'Q'??
         IF (CH.EQ.'S') TYPE='S'
         IF (CH.EQ.'N') THEN
            CALL CHKEDQ(STRING,WORK,'D',CH,FORM)
      ! return type='D' instead of type='Y'
            IF (CH.EQ.'Y') TYPE='D'
            IF (CH.EQ.'N') THEN
               CALL CHKEDQ(STRING,WORK,'Q',CH,FORM)
      ! return type='Q' instead of type='Y'
               IF (CH.EQ.'Y') TYPE='Q'
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END



      SUBROUTINE CHKEDQ(STRING,WORK,MATCH,TYPE,FORM)
      USE nrtype
      IMPLICIT real(8) (A-H,O-Z)
      ! Supplementary routine for checking validity of floating point numbers
      ! in the 'E', 'D' or 'Q'-format (if allowed/available ...) ...; the
      ! answer is returned in TYPE ('Y' or 'N'), FORM is a format string for
      ! this floating point number (if it is one ...) for FORTRAN-reading
      CHARACTER*(*) STRING,WORK,FORM
      CHARACTER*1   TYPE,CH,MATCH
      CHARACTER*255 PART1,PART2
      LOGICAL       LTEST,LDOT,ALLOWQ
      INTEGER       L,I,J,K,LEN,NOCCUR,LENGTH,INDEX,INDEXN
      INTRINSIC     LEN,INDEX,MAX
      EXTERNAL      NOCCUR,LENGTH,INDEXN
      ! Q-format (quadruple precision) is only supported on few machines ...
      PARAMETER (ALLOWQ=.TRUE.)

      ! default values ...
      TYPE='N'
      FORM=' '
      ! is 'Q' allowed/possile on this machine?
      IF ((MATCH.EQ.'Q').AND.(.NOT.ALLOWQ)) RETURN
      ! wrong matching type ...
      IF ((MATCH.NE.'E').AND.(MATCH.NE.'D').AND.(MATCH.NE.'Q')) RETURN
      LTEST=.TRUE.
      LDOT=.FALSE.
      K=0
      WORK=STRING
      CALL STRIP(WORK,L,'A')
      ! a float in E-format might contain one single 'E' or 'e', what is in
      ! front should be a float or an integer, what is after it an integer:
      ! (for D-format or Q-format just replace 'E' by 'D' or 'Q' ...). Of
      ! course there is also another valid (but more strange format without
      ! any 'E','D' or 'Q' -- just a signed exponent with some float in front;
      ! support also this valid format by searching such constructs ...
      CALL UPPER(WORK)
      IF (NOCCUR(WORK,MATCH,0).GT.1) LTEST=.FALSE.
      ! 'usual' format with 'E','D' or 'Q' here ...
      IF (LTEST.AND.(NOCCUR(WORK,MATCH,0).EQ.1)) THEN
      ! parse the result string ...
         CALL PARSE(WORK,PART1,PART2,MATCH,0)
         CALL STRIP(PART1,J,'A')
      ELSE IF (LTEST) THEN
      ! here the more strange (but also valid) format ..., it must at least
      ! contain one '+' or one '-' and at maximum two signs at all (the last
      ! is then relevant ...), so check all this here first ...
         I=NOCCUR(WORK,'-',0)
         J=NOCCUR(WORK,'+',0)
         IF (((I+J).NE.1).AND.((I+J).NE.2)) LTEST=.FALSE.
         IF (LTEST) THEN
      ! still all okay --> get last sign ...
            K=MAX(INDEXN(WORK,'-',I),INDEXN(WORK,'+',J))
      ! can not be possible ...
            IF ((K.LE.1).OR.(K.EQ.L)) LTEST=.FALSE.
            IF (LTEST) THEN
      ! well, seems still all to be okay -- parse the string ...
               PART1=WORK(1:K-1)
               CALL STRIP(PART1,J,'A')
               PART2=WORK(K:L)
            ENDIF
         ENDIF
      ENDIF
      ! if all is okay up to now do the rest ...
      IF (LTEST) THEN
      ! first part should be a floating point number ...
         CH=PART1(1:1)
      ! if it is a sign or a dot something non-blank MUST follow ...
         IF (((CH.EQ.'+').OR.(CH.EQ.'-').OR.(CH.EQ.'.')) &
     &                                 .AND.(L.LE.1)) LTEST=.FALSE.
         IF ((CH.NE.'+').AND.(CH.NE.'-').AND.(CH.NE.'0').AND. &
     &       (CH.NE.'1').AND.(CH.NE.'2').AND.(CH.NE.'3').AND. &
     &       (CH.NE.'4').AND.(CH.NE.'5').AND.(CH.NE.'6').AND. &
     &       (CH.NE.'7').AND.(CH.NE.'8').AND.(CH.NE.'9').AND. &
     &       (CH.NE.'.')) LTEST=.FALSE.
      ! only one single dot may be there, so remind if we had some already
         IF (CH.EQ.'.') LDOT=.TRUE.
         IF (.NOT.LTEST) GOTO 101
         DO 100 I=2,J
      ! Now only numbers may follow and somewhere a dot ...
            CH=PART1(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND. &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND. &
     &          (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND. &
     &          (CH.NE.'9').AND.(CH.NE.'.')) LTEST=.FALSE.
      ! only one single dot may be there and remind if we had some already
            IF ((CH.EQ.'.').AND.LDOT) LTEST=.FALSE.
            IF (CH.EQ.'.') LDOT=.TRUE.
            IF (.NOT.LTEST) GOTO 101
  100    CONTINUE
  101    CONTINUE
      ENDIF
      IF (LTEST) THEN
      ! if all is okay until here check finally the exponent (must be integer)
         CALL STRIP(PART2,J,'A')
      ! First character may be a sign (+/-) or a number ...
         CH=PART2(1:1)
      ! if it is a sign something non-blank MUST follow ...
         IF (((CH.EQ.'+').OR.(CH.EQ.'-')).AND.(L.LE.1)) LTEST=.FALSE.
         IF ((CH.NE.'+').AND.(CH.NE.'-').AND.(CH.NE.'0').AND. &
     &       (CH.NE.'1').AND.(CH.NE.'2').AND.(CH.NE.'3').AND. &
     &       (CH.NE.'4').AND.(CH.NE.'5').AND.(CH.NE.'6').AND. &
     &       (CH.NE.'7').AND.(CH.NE.'8').AND.(CH.NE.'9')) LTEST=.FALSE.
         IF (.NOT.LTEST) GOTO 201
         DO 200 I=2,J
      ! Now only numbers may follow ...
            CH=PART2(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND. &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND. &
     &          (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND. &
     &          (CH.NE.'9')) LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 201
  200    CONTINUE
  201    CONTINUE
      ENDIF
      ! if all was okay get the format (use always E-format ...):
      IF (LTEST) THEN
         TYPE='Y'
      ! 'strange format' shall return TYPE='S' instead of TYPE='Y' ...
         IF (K.NE.0) TYPE='S'
         IF (LEN(FORM).GE.12) THEN
      ! 'field length' in format is L, what is number of significant digits?
            I=0
            IF (LDOT) I=LENGTH(PART1)-INDEX(PART1,'.')
            WRITE(FORM,'(A1,I5,A1,I5)') 'E',L,'.',I
            CALL STRIP(FORM,I,'A')
         ENDIF
      ENDIF

      RETURN
      END



      INTEGER FUNCTION NXTFRU()
      USE nrtype
      IMPLICIT REAL(8) (A-H,O-Z)
      ! Find the next free unit number ...
      LOGICAL OCCUP
      INTEGER I
      NXTFRU=-1
      ! Usually the standard FORTRAN range for unit numbers is 0...99;
      ! on some systems also unit numbers beyond 99 might be allowed:
      ! if this is the case and you want to make use of it change it ...
      DO 100 I=0,99
      ! Units 0,5,6 are usually reserved for stderr,stdin,stdout ...
      ! If your system uses other standard I/O-units change it ... !
         IF ((I.EQ.0).OR.(I.EQ.5).OR.(I.EQ.6)) GOTO 100
         INQUIRE(UNIT=I,OPENED=OCCUP)
         IF (.NOT.OCCUP) THEN
            NXTFRU=I
            RETURN
         END IF
  100 CONTINUE
      RETURN
      END



      !!!LOPEN:  .TRUE. or .FALSE.
      !!!FNAM:   the name of the file
      !!!IU:     input Unit of the file 
      !!!WHAT:   the keyword to be searched
      !!!PCHAR:  '=' or ':' links 'WHAT' and the data list
      !!!CCHAR:  '!' or '#'... the comment flag
      !!!SCHAR:  ';' or 'r' seperates different groups of keywords and data
      !!!TYPE:   'S','I','F','C','L'
      !!!INTRES: the name of integer variable
      !!!FLTRES: the name of real variable
      !!!CMOERS: the name of complex variable
      !!!LOGERS: the name of logical variable
      !!!STR:    the name of character variable
      !!!N:      the number of variables in the data list
      !!!NMAX:   the maximum of N
      !!!IERR:   the error information

      !***********************************************************************
      !                                                                      *
      SUBROUTINE RDATAB(LOPEN,FNAM,IU,WHAT,PCHAR,CCHAR,SCHAR,TYPE, &
     &                     INTRES,FLTRES,CMPRES,LOGRES,STR,N,NMAX,IERR)
      USE nrtype
      IMPLICIT REAL(8) (A-H,O-Z)
      !                                                                      *
      !  This routine extracts data from some file (name is stored in FNAM). *
      !  The format is totally free, i.e. there may be as much blanks as you *
      !  like and the data may be located everywhere. The exact rules are:   *
      !                                                                      *
      !     -- no line/definition may be longer than 500   characters        *
      !     -- no 'keyword' may be longer than 255 characters and the        *
      !        'keyword' is the first what appears in each definition!       *
      !     -- PCHAR seperates the 'keyword' to be searched (WHAT) and the   *
      !        data list (usual recommendation is to use '=' or ':' ...)     *
      !     -- all after the character specified in CCHAR (recommendation    *
      !        is to use either '#' or '!' ...) is treated as a comment      *
      !     -- the character given in SCHAR shall seperate all definitions   *
      !        given in one line (so one may place more than one definition  *
      !        within one single line in order have the possibility to group *
      !        things easy together which have some common meaning/function; *
      !        recommendation is to use ';' or ',' ...)                      *
      !     -- different data within one list must be seperated by blanks    *
      !     -- no intermixing of different data types is allowed in a list!! *
      !     -- in string data all blanks are taken 'as is' (from character   *
      !        PCHAR up to the end of line or the next separator or comment) *
      !     -- else multiple blanks act like one blank within a list of data *
      !        or like 'no blank' if inserted between the 'keyword' and the  *
      !        character PCHAR or between character PCHAR and the data list  *
      !     -- because the characters given in PCHAR, CCHAR, SCHAR have now  *
      !        already some special meaning one needs some extra way to make *
      !        them accessible in string data: if one wants to use them in   *
      !        any input strings one has to precede this character with \    *
      !        (similar to UNIX-rules ...), because \ has also a special     *
      !        meaning \ must be written as \\ within string data ...        *
      !     -- data lists can be continued over several lines (with the only *
      !        restriction that the total input may not be longer than 32767 *
      !        characters ...). Continuations must be marked by some \ at    *
      !        the end of the previous line (but really at the end of the    *
      !        line!! -- so lines to be continued may not contain comments!) *
      !     -- the naming of the 'keywords' shall be case-insensitive (this  *
      !        would mean: 'TOKEN1', 'token1' and 'Token1' define in fact    *
      !        all the same keyword !!) -- internally uppercase will be used *
      !     -- for multiple entries (identical 'keywords') always the first  *
      !        occurence defines what is extracted, all other entries will   *
      !        be ignored!                                                   *
      !                                                                      *
      !  Maybe this routine is not yet perfect, but well usuable ...         *
      !                                                                      *
      !  There may be specified five types (given in variable TYPE) for the  *
      !  data to be extracted:                                               *
      !                                                                      *
      !     -- 'S'    read some string       --> result goes to variable STR *
      !     -- 'I'    read integer data list --> result goes to array INTRES *
      !     -- 'F'    read real    data list --> result goes to array FLTRES *
      !     -- 'C'    read complex data list --> result goes to array CMPRES *
      !     -- 'L'    read logical data list --> result goes to array LOGRES *
      !                                                                      *
      !  for TYPE=I,F,C and L the variable N returns the number of elements  *
      !  found in the data list, if an error occurred (e.g. wrong format) N  *
      !  will be zero (and IERR non-zero). For TYPE=S variable N returns the *
      !  position of the last non-blank character in variable STR ... .      *
      !                                                                      *
      !  IERR returns an error code if any problem occured (normal: IERR=0)  *
      !     -- IERR=1   no free I/O-unit found to open file                  *
      !     -- IERR=2   OPEN error (e.g. file not found, ...)                *
      !     -- IERR=3   'keyword' (WHAT) not found on specified file         *
      !     -- IERR=4   invalid data type (TYPE)                             *
      !     -- IERR=5   error reading/parsing data list (check format!)      *
      !     -- IERR=6   cannot open scratch file for conversion of data      *
      !                                                                      *
      !***********************************************************************

      CHARACTER*(*)   FNAM,WHAT,STR
      CHARACTER*1     PCHAR,CCHAR,SCHAR,TYPE,BS
      CHARACTER*2     HIDDEN
      CHARACTER*255   KEY
      CHARACTER*255   BUFLIN
      CHARACTER*255   WORK,WPARSE
      INTEGER         N,IERR,INTRES(*),I,IU,L,NMAX,IKEY,LKEY,INOSL
      REAL(8)            FLTRES(*)
      COMPLEX(8)         CMPRES(*)
      LOGICAL         LOGRES(*),FOUND,CONT,YINTST,LOPEN
      INTEGER         LENGTH,NWORDS,NXTFRU,IUSCR,LMAX,NREAD,NITEMS
      ! Following parameter should probably be customized if porting this
      ! program from one machine to another machine: Usually in FORTRAN77
      ! it is  not   allowed to use the '*'-format for internal reads or
      ! writes, but we try to use it (if possible!). If your OS/compiler
      ! allows internal read/writes with the '*'-format (like AIX on RS6000)
      ! then set YINTST=.TRUE., otherwise you must use YINTST=.FALSE.)!
      PARAMETER       (YINTST=.TRUE.)
      ! Define 'backslash', some UNIX-systems want \\, some only \ ...
      ! so take the following and it should be highly portable ...
      PARAMETER       (BS='\\')
      EXTERNAL        LENGTH,NWORDS,NXTFRU,NITEMS

      ! Initialise N and IERR:
      N=0
      IERR=0

      ! Invalid data type given: return error code '4'
      IF ((TYPE.NE.'I').AND.(TYPE.NE.'F').AND. &
     &    (TYPE.NE.'C').AND.(TYPE.NE.'S').AND.(TYPE.NE.'L')) THEN
         IERR=4
         RETURN
      ENDIF

      ! Search some free unit where to open the input file FNAM if necessary:
      IF ((IU.LT.0).OR.(IU.GT.99)) IU=NXTFRU()
      ! Hmmm ... . No free unit found -- you should not use so much files!
      IF (IU.LT.0) THEN
         IERR=1
         RETURN
      ENDIF

      ! Try to open the data file ...
      IF (LOPEN) THEN
         WORK=FNAM
         CALL STRIP(WORK,L,'A')
         IF (L.EQ.0) GOTO 10
         OPEN(IU,FILE=WORK(1:L),STATUS='OLD',ERR=10)
      ELSE
         REWIND IU
      ENDIF
      GOTO 20
      ! Hmmm ... . No success? Maybe you should first create the file ...
   10 IERR=2
      IF (LOPEN) CLOSE(IU,ERR=11)
   11 CONTINUE
      RETURN
   20 CONTINUE
      ! Set some flag telling us whether we have found what we searched:
      FOUND=.FALSE.
      ! ... and set some flag telling us something about continuation lines:
      CONT=.FALSE.

      ! First we must somehow get information on the 'keyword' stored in WHAT:
      ! how many (non-blank) characters has the keyword really after stripping
      ! all leading and trailing blanks? Try to find out ...
      KEY=WHAT
      CALL STRIP(KEY,LKEY,'A')
      CALL UPPER(KEY)

      ! Now read the file line by line until we find the correct 'keyword'
      ! (or not ...) -- label 30 actually is an entry to an implicit loop!
   30 CONTINUE
      ! Read until the end of the file or until some error occurs ...
      READ(IU,'(A)',ERR=10000,END=10000) BUFLIN

      ! Continuation line has to be appended to current string!
      IF (CONT) THEN
         WORK=WORK(1:LMAX)//BUFLIN
         BUFLIN=WORK
      ENDIF
      ! Where is the 'end of the input string' ...:
      LMAX=LENGTH(BUFLIN)
      ! Continuation line(s) to be expected?
      INOSL=0
      DO 40 I=LMAX,1,-1
         IF (BUFLIN(I:I).NE.BS) THEN
            INOSL=I
            GOTO 50
         ENDIF
   40 CONTINUE
   50 CONT=(MOD(LMAX-INOSL+2,2).EQ.1)
      ! Read next line until this entry is complete ...:
      IF (CONT) THEN
         LMAX=LMAX-1
         WORK=BUFLIN(1:LMAX)
         GOTO 30
      ENDIF

      ! Forget about leading blanks ...:
      CALL STRIP(BUFLIN,LMAX,'L')
      ! The real end is given by character CCHAR: all behind it is a comment!
      IF (BUFLIN(1:1).EQ.CCHAR) THEN
         BUFLIN=' '
         LMAX=0
         GOTO 70
      ENDIF
      DO 60 I=2,LMAX
         IF ((BUFLIN(I:I).EQ.CCHAR).AND.(BUFLIN(I-1:I-1).NE.BS)) THEN
            BUFLIN=BUFLIN(1:I-1)
            LMAX=I-1
            GOTO 70
         ENDIF
   60 CONTINUE
   70 CONTINUE
      ! Following means: only blanks after removing comment, nothing is left:
      IF (LMAX.EQ.0) GOTO 30

      ! Do we have more than one definition on one line? ---> parse!
      IF (LMAX.EQ.1) THEN
         IF (BUFLIN(1:1).EQ.SCHAR) THEN
            WORK=' '
            BUFLIN=' '
            LMAX=0
            GOTO 90
         ENDIF
      ENDIF
      IF (BUFLIN(1:1).EQ.SCHAR) THEN
         WORK=' '
         BUFLIN=BUFLIN(2:LMAX)
         LMAX=LMAX-1
         GOTO 90
      ENDIF
      DO 80 I=2,LMAX
         IF ((BUFLIN(I:I).EQ.SCHAR).AND.(BUFLIN(I-1:I-1).NE.BS)) THEN
            WORK=BUFLIN(1:I-1)
            IF (I.LT.LMAX) BUFLIN=BUFLIN(I+1:LMAX)
            IF (I.EQ.LMAX) BUFLIN=' '
            LMAX=LMAX-I
            GOTO 90
         ENDIF
   80 CONTINUE
      ! Was apparently the last definition (nomore separators found ...):
      WORK=BUFLIN
      BUFLIN=' '
      LMAX=0
   90 WPARSE=WORK
      ! Okay! Is the 'keyword' there?
      CALL STRIP(WPARSE,L,'A')
      ! Well, some separator at the beginning of the line!? Try the next ...
      IF (L.EQ.0) GOTO 70
      CALL UPPER(WPARSE)
      IKEY=INDEX(WPARSE(1:L),KEY(1:LKEY))
      ! The keyword must be the first what appears in the definition string!!
      ! After stripping  all  blanks it must start at the first character ---
      ! otherwise we have found something else (but not meaningful!!!!!!!!!!):
      IF (IKEY.EQ.1) THEN
      ! Maybe we got it! But only if a parsing character follows immediately!
         IF ((1+LKEY).GT.L) THEN
      ! Impossible: PCHAR would be beyond the end of the string if this holds!
      ! Try the next definition ...
            GOTO 70
         ENDIF
         IF (WPARSE(1+LKEY:1+LKEY).NE.PCHAR) THEN
      ! No! We have found the word (substring ...) but no assignment character
      ! immediately after it --- must be some other entry or 'open comment'
      ! containing the keyword just as a substring; try next definition ...
            GOTO 70
         ENDIF
      ! Okay! Got it!!!!
         BUFLIN=WORK
         LMAX=LENGTH(BUFLIN)
         GOTO 100
      ELSE
      ! Try the next definition ...
         GOTO 70
      ENDIF
  100 FOUND=.TRUE.

      ! Parse off the data list now!
      CALL PARSE(BUFLIN,WPARSE,WORK,PCHAR,0)
      BUFLIN=WORK
      L=MAX(LENGTH(WPARSE),1)
      ! Hmmm ... . Have not yet found the correct parsing character, next try!
      IF (WPARSE(L:L).EQ.BS) GOTO 100

      ! And here we have it: get the data from the list!
      IF (TYPE.EQ.'S') THEN
      ! String -- that is tooooooooo easy ... :
         WORK=BUFLIN
      ! 'translate' all 'special characters':
         N=0
      ! the 'parsing character' ...
         HIDDEN=BS//PCHAR
         CALL REPLAC(WORK,HIDDEN,PCHAR,N,0)
      ! the 'comment character' ...
         HIDDEN=BS//CCHAR
         CALL REPLAC(WORK,HIDDEN,CCHAR,N,0)
      ! the 'separation character' ...
         HIDDEN=BS//SCHAR
         CALL REPLAC(WORK,HIDDEN,SCHAR,N,0)
      ! and the 'backslash' ...
         HIDDEN=BS//BS
         CALL REPLAC(WORK,HIDDEN,BS,N,0)
      ! You do not want more than NMAX characters ??
         NREAD=N
         IF (NMAX.GT.0) NREAD=MIN(N,NMAX)
         STR=WORK(1:NREAD)
      ELSE
      ! List of numbers/logicals (checking also types) --> how much valid??
         N=NITEMS(BUFLIN,WORK,.FALSE.,TYPE)
      ! There seems to be a fatal format error?
         IF (N.LT.1) GOTO 110
      ! Well, maybe we want/cant never read more than NMAX data ...
         NREAD=N
         IF (NMAX.GT.0) NREAD=MIN(N,NMAX)
      ! Now really read the data (no responsibilty for core dumps here due to
      ! insufficient dimensioning of the arrays 'coming from above' ... !!):
         IF (YINTST) THEN
      ! Well, here is the 'critical stuff' with YINTST ... . If set .FALSE. the
      ! compiler should hopefully note "code unreachable" and continue without
      ! complaining about the invalid '*'-format ... . If you still run into
      ! trouble because your compiler tries to translate the following section
      ! although it can never be executed then comment out the following block
      ! in order to get the code working ... (shit FORTRAN77 ... !).
            IF (TYPE.EQ.'I') THEN
               READ(BUFLIN,*,ERR=110,END=110) (INTRES(I),I=1,NREAD)
            ELSE IF (TYPE.EQ.'F') THEN
               READ(BUFLIN,*,ERR=110,END=110) (FLTRES(I),I=1,NREAD)
            ELSE IF (TYPE.EQ.'C') THEN
               READ(BUFLIN,*,ERR=110,END=110) (CMPRES(I),I=1,NREAD)
            ELSE IF (TYPE.EQ.'L') THEN
               READ(BUFLIN,*,ERR=110,END=110) (LOGRES(I),I=1,NREAD)
            ENDIF
         ELSE
            IUSCR=NXTFRU()
            IF (IUSCR.LT.0) THEN
               IERR=6
               N=0
               RETURN
            ENDIF
            OPEN(IUSCR,STATUS='SCRATCH',ERR=130)
            WRITE(IUSCR,'(A)') BUFLIN
            REWIND IUSCR
            IF (TYPE.EQ.'I') THEN
               READ(IUSCR,*,ERR=110,END=110) (INTRES(I),I=1,NREAD)
            ELSE IF (TYPE.EQ.'F') THEN
               READ(IUSCR,*,ERR=110,END=110) (FLTRES(I),I=1,NREAD)
            ELSE IF (TYPE.EQ.'C') THEN
               READ(IUSCR,*,ERR=110,END=110) (CMPRES(I),I=1,NREAD)
            ELSE IF (TYPE.EQ.'L') THEN
               READ(IUSCR,*,ERR=110,END=110) (LOGRES(I),I=1,NREAD)
            ENDIF
            CLOSE(IUSCR,ERR=130)
         ENDIF
         GOTO 120
      ! Hmmm ... . Something was not okay with the format ...
  110    IERR=5
         N=0
         GOTO 10000
  120    CONTINUE
         GOTO 140
  130    IERR=6
         N=0
         CLOSE(IUSCR,ERR=10000)
         GOTO 10000
  140    CONTINUE
      ENDIF

      ! Yepeeh! That is the end folks -- successful or not (who cares ...)!
10000 CONTINUE
      ! Close the input file ...
      IF (LOPEN) CLOSE(IU,ERR=10001)
      ! Hmmm ... . No success? Maybe you should add something in your data
      ! file, or you should check the format ... (??)  --  LOOSER!!
      IF (.NOT.FOUND) IERR=3
10001 CONTINUE
      ! ... and bye bye my honey -- would be nice to meet you again ...
      RETURN

      END

