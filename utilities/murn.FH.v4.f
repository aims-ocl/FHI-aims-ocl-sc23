      PROGRAM MURN
!
!     fitting of a murnaghan equation of state. 
!
!     After years, discovered a typo in the variable initializing
!       B0PRIM that invalidated the initialization completely. That 
!       explains a lot of the problems with this fitting program. Well.
!       The initialization is now correct. : VB, Feb 23, 2017
!     increased max. data dimensions from 50 to 10000 : VB, Feb 22, 2017
!     modified from murn.vb.v3.f : Felix Hanke, April 2, 2009
!     modified again: Volker Blum
!     modified from murn2.f : R.S 19.04.91
!     modified from murn1.f : gps 6.12.90
!     least-squares fit of E vs V by Murnaghan equation of state
!
!  ---- begin documentation ----------------------------------------------------
!
!  Usage description for murnaghan equation of state fit. In principle
!  this program can be run with just three "data_point" commands, if all
!  the other defaults are appropriate. Anything beyond the named defaults
!  should be specified using the keywords described below. These keywords
!  require one additional input for each <...> specified in the
!  documentation. All assumptions made during runtime (units, unit cell
!  volumes) etc are documented in the output and can be changed using the
!  keywords below.
!
!  All the input information has to be written into a file named
!  "murn_fit.in", after which the compiled binary should be called in the
!  same directory as that file. The entire output is contained in an
!  output stream, there are no additional output files.
!
!  ---- core keyword for fit input ---------------------------------------------
!
!  data_point <lc / volume> <energy>
!     specification of an input data point for the fit. Parameters are the
!     lattice constants (or the unit cell volume) and the energy. Care
!     must be taken that these are in the proper units of length and
!     energy as specified either by the remaining input or the default
!     settings (angstrom and eV respectively).
!
!     see also keyword input_volumes_only below.
!
!  ---- optional keywords ------------------------------------------------------
!
!  # - lines beginning with "#" are treated as comments
!
!  length_unit <angstrom / bohr / value >
!     specifies unit of length. Default: angstrom
!     if a real-numbered value is specified then the
!     length unit is that number of angstroms.
!
!  energy_unit <eV / hartree / rydberg / value>
!     specifies unit of energy. Default: eV
!     if a real-numbered "value" is specified then
!     the energy unit is that number of eV.
!
!  input_volumes_only
!     specifies that the input (and output!) is only in terms
!     of volumes, rather than considering lattice constants.
!     Default: lattice-constant based input
!
!  unit_cell_volume  <sc / fcc / bcc / value>
!     specifies the conversion factor between (lattice_constant)^3
!     and the volume of each primitive unit cell
!     Default: sc (= 1)
!     pre-defined values are:
!     sc  (simple cubic)        = 1
!     fcc (face-centered cubic) = 0.25
!     bcc (body-centered cubic) = 0.5
!     all other unit cell shapes can be covered by specifying a real
!     number "value"
!
!  fit_range <amin> <amax> <Npoints>
!     requests optional output of a fitted curve. amin and amax are the
!     starting and ending lattice constants (unit cell volumes) for the
!     fit, and Npoints is the number of points plotted.
!
!  print_input
!     prints the original input in the same format as the other fits
!
!  print_fit
!     prints the direct fit to the original input data for comparison
!
!  init_a0  <value>
!     allows optional specification of an initial guess for the lattice
!     constant (or unit cell volume) in case the fit is unstable.
!     Default: lattice constant / volume with minimum energy from input data
!
!  init_e0 <value>
!     optional specification of a lattice energy in case of an unstable
!     fit. Default: minimum value of the input energies
!
!  init_B0 <value>
!     optional specification of an initial guess for the bulk modulus.
!     Default: 1 MBar
!
!  init_B0prime <value>
!     optional specification of an initial guess for the pressure
!     derivative of the bulk modulus. Default: 4
!
!  ---- example input file -----------------------------------------------------
!
!   # sample input for finding the lattice constant of Mg
!   unit_cell_volume 0.25
!   fit_range        4.15 4.35 50
!   print_input
!   data_point       4.15   -2.15930
!   data_point       4.20   -2.20242
!   data_point       4.25   -2.21995
!   data_point       4.30   -2.21451
!   data_point       4.35   -2.18849
!
!
!  ---- end documentation ------------------------------------------------------
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EDATA(10000),VDATA(10000),X(40),AUX(40),VIN(10000)
      DIMENSION EIN(10000)
      character*80 input_name, keyword, length_unit, energy_unit
      integer      ioerr
      logical      init_a0, init_e0, init_B0, init_B0prime
      logical      init_volume, init_length, init_energy, input_volume
      logical      init_fitrange, print_input, print_fit
      EXTERNAL     FUN

C CONVERSION FACTOR FROM RYDBERG TO EV:
      data CONVYY /13.60569193/
      data BOHR   /0.52917720859/
      COMMON /C/ EDATA,VDATA, VOLMIN,VOLMAX,B0MIN,B0MAX,B0PMIN,B0PMAX,U
     &LA,NNN
C
      IOUT=6
      IN  = 20
      ! default input name is murn_fit.in, but one input can be specified if desired
      if (iargc().ge.1) then
         call getarg(1,input_name)
      else
         input_name = "murn_fit.in"
      end if

      write(*,'(A)')'# '
      write(*,'(A)')'# Least-squares fit of etot(v) by Murnaghan eq.'
      write(*,'(A)')'# ' 

C input parser, everything is keyword driven, output is being printed as 
C written during read
      write(*,'(2A)') '# Parsing control file ', trim(input_name)
      open(unit = 20, file = input_name, status = 'old')
      ! set various keyword flags
      init_a0       = .false.
      init_e0       = .false.
      init_B0       = .false.
      init_B0prime  = .false.
      init_volume   = .false.
      init_length   = .false.
      init_energy   = .false.
      input_volume  = .false.
      init_fitrange = .false.
      print_input   = .false.
      print_fit     = .false.
      ioerr         = 0
      NNN           = 0   ! no input data points, count them as we go along
      ! defaults for starting values of guesses:
      ALATT0  = -1d0
      E0EV    = -1d0
      B0      = -1d0
      B0PRIM  = -1d0
      do while (ioerr.eq.0)
         read(unit = IN, fmt = *, iostat = ioerr) keyword
         if (ioerr.eq.0) then
            if (keyword(1:1).eq.'#') then    ! comment
            else if (keyword.eq.'length_unit') then 
               init_length = .true.
               backspace(IN)
               read(IN,*) keyword, keyword
               if (keyword.eq.'angstrom') then
                  ULA = 1d0
                  length_unit = 'angstrom' 
               else if (keyword.eq.'bohr') then
                  ULA = BOHR
                  length_unit = "bohr"
               else
                  backspace(20)
                  read(IN,*) keyword, ULA
                  write(length_unit,'(F10.6,A)') ULA,'angstrom'
               end if 
               write(*,'(A,F15.6,A)')"# Found unit of length ",ULA,
     +              " Angstrom."
            else if (keyword.eq.'energy_unit') then 
               init_energy = .true.
               backspace(IN) 
               read(IN,*) keyword, keyword
               if (keyword.eq.'eV') then
                  CONV_inp = 1d0
                  write(*,'(A)') "# All energies in units of eV."
                  energy_unit = 'eV'
               else if (keyword.eq.'hartree') then
                  CONV_inp = 2d0*CONVYY
                  write(*,'(A)') "# All energies in units of hartree."
                  energy_unit = 'hartree'
               else if (keyword.eq.'rydberg') then
                  CONV_inp = CONVYY
                  write(*,'(A)') "# All energies in units of rydberg."
                  energy_unit = 'rydberg'
               else
                  backspace(IN) 
                  read(IN,*) keyword, CONV_inp
                  write(*,'(A,F10.6,A)') "# Found energy unit ",
     +                 CONV_inp, " eV."
                  write(energy_unit,'(F10.6,A)')CONV_inp, ' eV'
               end if 
            else if (keyword.eq."input_volumes_only") then
               write(*,'(2A)') "# All lattice information is handled ",
     +              "in terms of VOLUMES ONLY!"
               input_volume = .true.
            else if (keyword.eq."unit_cell_volume") then
               init_volume = .true.
               backspace(IN)
               read(IN,*) keyword, keyword
               if (keyword.eq.'sc') then
                  CONVAV = 1d0
                  write(*,'(A)') "# Using simple cubic lattice volume."
               else if (keyword.eq.'bcc') then
                  CONVAV = 5d0
                  write(*,'(A)') "# Using bcc lattice volumes."
               else if (keyword.eq.'fcc') then
                  CONVAV = 0.25d0
                  write(*,'(A)') "# Using fcc lattice volumes."
               else 
                  read(keyword,*) CONVAV
                  write(*,'(2A,F10.6)') "# using a unit cell volume ",
     +                 "of a^3* ",CONVAV
               end if
            else if (keyword.eq."fit_range") then
               backspace(IN)
               read(IN,*) keyword, alat_min, alat_max, nr_alat
               init_fitrange = .true.
               write(*,'(A)') "# Found range for output fit: "
               write(*,'(A,F10.6)') "# | Start: ", alat_min
               write(*,'(A,F10.6)') "# | End:   ", alat_max
               write(*,'(A,I4)')    "# | number of points: ", nr_alat
               if ((nr_alat.lt.2).or.(alat_min.ge.alat_max)) then
                  write(*,'(2A)') "# * WARNING: something wrong with ",
     +                 "this range! Aborting."
                  stop
               end if 
            else if (keyword.eq.'print_input') then
               print_input = .true.
               write(*,'(2A)') "# will print original input data at ",
     +              "end of calculation. "
            else if (keyword.eq.'print_fit') then
               print_fit = .true.
               write(*,'(2A)') "# will print fitted input at end of ",
     +              "calculation. "
            else if (keyword.eq.'init_a0') then
               init_a0      = .true.
               backspace(IN)
               read(IN,*) keyword, ALATT0
               write(*,'(A,F10.6)') "# initial guess for a: ", ALATT0
            else if (keyword.eq.'init_e0') then
               init_e0      = .true.
               backspace(IN)
               read(IN,*) keyword, E0EV
               write(*,'(A,F10.6)') "# initial guess for E: ", E0EV
            else if (keyword.eq.'init_B0') then
               init_B0      = .true.
               backspace(IN)
               read(IN,*) keyword, B0
               write(*,'(A,F10.6)') "# initial guess for B: ", B0
            else if (keyword.eq.'init_B0prime') then 
               init_B0prime = .true.
               backspace(IN)
               read(IN,*) keyword, B0PRIM
               write(*,'(A,F10.6)')"# initial guess for B0prime: ",
     +              B0PRIM
            else if (keyword.eq.'data_point') then
               ! new fitting point found ...
               backspace(IN)
               NNN = NNN + 1
               ! not very elegant but should work for the time being:
               ! no more than 10000 input points, easy to fix if someone really cares
               if (NNN.gt.10000) then
                  write(*,'(2A)')"# * WARNING: only 10000 input points ",
     +                 "allowed"
                  write(*,'(A)') "# * Aborting."
                  stop
               end if 
               ! read point, knowing that it will have to be changed 
               ! to the proper units later ... 
               read(IN,*) keyword, VIN(NNN), EIN(NNN)
               write(*,'(A,F10.6,F30.8)') "# Input data point: ",
     +              VIN(NNN),EIN(NNN)
            else
               write(*,'(2A)') "# * WARNING: Keyword unknown: ",keyword
               write(*,'(A)')  "# * Aborting. "
               stop
            end if
         end if 
      end do
      close(20)

      ! consistency checks of input parameters 
      if (NNN.lt.3) then   ! need at least three points for fit, quit if not present!
         write(*,'(2A)') "# * WARNING: Do not have enough input ",
     +        "for meaningful fit."
         write(*,'(A)') "# * Aborting."
         stop
      end if
      if (.not.init_length)then
         ULA=1.0
         WRITE(*,'(A)')"# using default unit of length of angstrom."
         length_unit = 'angstrom'
      end if
      if (.not.init_energy) then
         CONV_inp = 1d0
         write(*,'(A)') "# Using default energy unit of eV."
         energy_unit = 'eV'
      end if 
      if (.not.init_volume) then
         CONVAV = 1d0
         if (.not.input_volume) then
            ! input lattice constant, should tell user about implicit assumption made
            write(*,'(A)')"# Defaulting to simple cubic lattice volumes"
         end if 
      end if
      if ((.not.init_fitrange).and.(.not.print_fit)) then
         write(*,'(2A)') "# no requested fit output, will only find ",
     +        "optimal lattice parameters"
      else
         if (init_fitrange) then
            ! fit range needs to be in units of A
            if (.not.input_volume) then
               alat_min = alat_min / ULA
               alat_max = alat_max / ULA
            else 
               alat_min = alat_min / (ULA**3)
               alat_max = alat_max / (ULA**3)
            end if 
         end if
      end if
      
      ! change input parameters into those useful for the fitting procedure
      do I=1, NNN
         if (.not.input_volume) then   ! have lattice constants, require volumes ... 
            VDATA(I) = VIN(I) / ULA
            VDATA(I) = VDATA(I)**3*CONVAV
         else 
            VDATA(I) = VIN(I)*CONVAV
         end if
         EDATA(I)=EIN(I)*CONV_inp
      end do

!     set initial guesses if not specified properly
      if ((E0EV.lt.0).or.(ALATT0.lt.0)) then
        E0MIN=1.D6
        VOLMIN=1.D6
        do I = 1, NNN
           if (EDATA(I).lt.E0MIN) then 
              E0MIN = EDATA(I)
              VOLMIN = VDATA(I)              
           end if
        end do 
      end if

!     energy: starting energy in Rydberg
      if (E0EV.lt.0) then
        E0EV=E0MIN
        E0RY=E0EV/CONVYY
      else
        E0RY=E0EV/CONVYY
      end if 

C lattice const.
      if (ALATT0.lt.0d0) then         
        VOL0=VOLMIN        
        ALATT0=(VOL0/CONVAV)**(1.D0/3.D0)
      else 
         if (input_volume) then 
            VOL0 = ALATT0*CONVAV
         else
            VOL0 = (ALATT0**3.D0)*CONVAV
         end if 
      end if
C Bulk modulus and derivative
      if (B0.lt.0) then
        B0=1.D0
      end if

      if (B0PRIM.lt.0) then
        B0PRIM=4.D0
      end if

      write (IOUT,'(A)') "# "      
      write (IOUT,'(A)') "# Initialising fitting procedure to:"
      write (IOUT,'(A)') "# "

      WRITE(IOUT,1130)
 1130 FORMAT('#',7X,'a0 (bohr)',2X, 'VOL0 (ULA**3)',2X,'B0 (MBAR)',
     +     5X,'B0PRIME',6X,'E0 (EV)',6X,'E0 (RY)'/'#')

      WRITE(IOUT,1140)ALATT0,VOL0,B0,B0PRIM,E0EV,E0RY
 1140 FORMAT('#',5X,F9.4,3F13.5,2F13.5)

 1120 FORMAT('# ',79(1H-)/'#')

      ! start actual computational work
      B0MIN=0.01D0
      B0MAX=10.D0
      B0PMIN=0.01D0
      B0PMAX=10.D0
C
      VOLMIN=0.1D0*VOL0
      VOLMAX=3.5D0*VOL0
      WRITE(IOUT,1210)VOLMIN,B0MIN,B0PMIN,VOLMAX,B0MAX,B0PM
     &AX
 1210 FORMAT('# MIN:',9X,3F13.5/'# MAX:',9X,3F13.5/'#')
C
C     A UNIFORM SHIFT OF ALL THE ENERGIES
C     (WILL NOT INFLUENCE THE B0, B0PRIME, VOL0, ONLY SHIFTS E0):
C
      SHIFT=-E0EV-1.D0
      DO 20 I=1,NNN
   20 EDATA(I)=EDATA(I)+SHIFT
C
C     MURNAGHAN LEAST-SQUARES FIT:
      WRITE(IOUT,1120)
      WRITE(IOUT,1280)
 1280 FORMAT('# FIT BY MURNAGHAN EQUATION EQUATION OF STATE')
      WRITE(IOUT,1120)
      WRITE(IOUT,1190)
 1190 FORMAT('# ITER',1X, 'VOL0 (ULA**3)',2X,
     1   'B0 (MBAR)',5X,'B0PRIME',6X,'E0 (EV)',4X,'SUM OF SQUARES'/
     2    '#',74X,'IERR'/'#')
C
      X(1)=E0EV+SHIFT
      X(2)=B0
      X(3)=B0PRIM
      X(4)=VOL0
      NVAR=4
      LIM=20
C
      DO 70 I = 1, 75
      CALL DFMND(FUN,X,FFF,NVAR,LIM,AUX,IERR)
C
C     THE RESULT OF 20 ITERATION :
      WRITE(IOUT,1250) 20*I,X(4),X(2),X(3),
     1                   X(1)-SHIFT,FFF,IERR
 1250 FORMAT('#',I4,5F13.5/'#',63X,I15)
C
      IF(IERR .EQ. 0) GO TO 80
C
   70 CONTINUE
   80 CONTINUE

C     THE RESULT OF THE LAST ITERATION:
      WRITE(IOUT,1250)20*I,X(4),X(2),X(3),
     1                   X(1)-SHIFT,FFF,IERR
C
      WRITE(IOUT,1120)
      VOL0=X(4)
      B0=X(2)
      B0PRIM=X(3)
      E0EV=X(1)-SHIFT
      E0RY=E0EV/CONVYY
      RMS = fun(X)

! final formatted output, to many sig. figures
! units have to match input units, use conv_inp to achieve that. 
      write(*,'(A)') "# Output of the final result for this fit "
      if (.not.input_volume) then
         alat0=(VOL0/CONVAV)**(1./3.)
         alat0 = alat0 * ULA
         write(*,'(A,F20.10,X,A)') "# lattice constant = ", alat0,
     +        trim(length_unit)
         
      else 
         write(*,'(A,F20.10,3A)') "# unit cell volume = ", VOL0,
     +     " (",trim(length_unit),")^3"
      end if         
      write(*,'(A,F20.10,X,A)') "# minimum energy   = ", E0EV/conv_inp,
     +     trim(energy_unit)
!     bulk modulus is hardwired in MBar, see const convxx in function murng1 below
      write(*,'(A,F20.10,A)') "# bulk modulus     = ",b0, " MBar"  
      write(*,'(A,F20.10)')   "# B0prime          = ", b0prim
      write(*,'(A)') "#" 
      write(*,'(A,E8.2,A4)')  "# RMS fit error    = ", sqrt(RMS)*1000.0,
     +                            " meV"
      write(*,'(A)') "#" 

      ! print completed fit if requested
      if (init_fitrange) then 
         alatt0=ula*alatt0
         WRITE(IOUT,1120)
         write(*,'(A)') "# Plot of fitted curve using murnaghan eqn."
         if (input_volume) then 
            write(*,'(5A)') "# unit cell volume (",trim(length_unit),
     +           ")^3  energy (",trim(energy_unit),')'
         else 
            write(*,'(5A)') "# lattice constant (",trim(length_unit),
     +           ")   energy (",trim(energy_unit),")"
         end if
!     print fitted values of lattice constant/volume vs energy
         do i=1,nr_alat
            ! calculate lattice constant in angstrom or volume in (A)^3
            alatt=alat_min + dble(i-1)*(alat_max-alat_min)/
     +           dble(nr_alat-1)
            if (.not.input_volume) then
               vol   = alatt**3 * convav ! calculate unit cell volume in A^3
               alatt = alatt*ula         ! change lc back into original unit of length
            else 
               vol   = alatt * convav    ! calculate proper unit cell volume in units of A^3
               alatt = alatt*(ula)**3    ! change back into original unit ov volume
            end if            
            ! input units: eV & Angstrom
            call murng1(ula,vol,vol0,b0,b0prim,e0ev,etot)
            write(iout,502) alatt, etot/conv_inp   
         enddo
         write(*,*)
      end if

      ! print input data if requested
      if (print_input) then
         write(*,'(A)') "# printout of the original input data"
         if (input_volume) then 
            write(*,'(5A)') "# unit cell volume (",trim(length_unit),
     +           ")^3  energy (",trim(energy_unit),')'
         else 
            write(*,'(5A)') "# lattice constant (",trim(length_unit),
     +           ")   energy (",trim(energy_unit),")"
         end if
         do i = 1, NNN
            write(iout,502) VIN(i), EIN(i)
         end do
         write(*,*)
      end if

      ! print fit to input data if requested
      if (print_fit) then
         write(*,'(A)') "# printout of the direct fit to input data"
         if (input_volume) then 
            write(*,'(5A)') "# unit cell volume (",trim(length_unit),
     +           ")^3  energy (",trim(energy_unit),')'
         else 
            write(*,'(5A)') "# lattice constant (",trim(length_unit),
     +           ")   energy (",trim(energy_unit),")"
         end if
         do i = 1, NNN
            ! calculate input volume
            vol = VIN(i)
            if (input_volume) then 
               vol = vol *convav / (ula**3)
            else
               vol = (vol/ula)**3 * convav
            end if
            call murng1(ula,vol,vol0,b0,b0prim,e0ev,etot)
            write(iout,502) VIN(i), etot/conv_inp               
         end do
         write(*,*)
      end if


  502 format(f15.5,f25.7)
c
      STOP
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DOUBLE PRECISION FUNCTION FUN(X)
C FUNCTION TO BE MINIMIZED IN THE LEAST-SQUARES FIT (BY SUBROUTINE DFMN
C  D)
C FUNCTION
C CALCULATES THE SUM OF THE  A B S O L U T E  DEVIATIONS
C            (E(THEOR)-E(EXP))**2
C DIVIDED BY THE NUMBER OF EXP. POINTS,
C ASSUMING FOR EQUATION OF STATE THE MURNAGHAN EXPRESSION.
C
C MEANING OF VARIABLES:
C      X(1) .... E0
C      X(2) .... B0
C      X(3) .... B0PRIM
C      X(4) .... VOL0
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EDATA(10000),VDATA(10000),X(40)
      COMMON /C/ EDATA,VDATA, VOLMIN,VOLMAX,B0MIN,B0MAX,
     1                                       B0PMIN,B0PMAX,ULA,NNN
C
C         *         *         *         *         *        *
      E0=X(1)
      B0=X(2)
      B0PRIM=X(3)
      VOL0=X(4)
C
      SUM=0.D0
C     THE SUM OF SQUARES:
      DO 10 I=1,NNN
      VOLACT=VDATA(I)
      CALL MURNG1(ULA,VOLACT,VOL0,B0,B0PRIM,E0,ETOT)
      SUM=SUM+(ETOT-EDATA(I))**2
   10 CONTINUE
      FUN=SUM/DFLOAT(NNN)
      RETURN
      END
C -----------------------------------------------------------------
      SUBROUTINE MURNG1(ULA,VOL,VOL0,B0,B0PRIM,E0,ETOT)
C EVALUATION OF THE MURNAGHAN EXPRESSION FOR ENERGY AS A FUNCTION
C OF VOLUME.
C
C INPUT DATA:
C
C      ULA ..... UNIT OF LENGTH, IN ANGSTROEMS, USED HERE.
C      VOL ..... VOLUME, IN THE ABOVE UNITS OF LENGTH CUBED.
C      VOL0 .... VOLUME AT THE ENERGY MINIMUM, IN THE SAME UNITS.
C      B0 ...... BULK MODULUS, IN UNITS MEGABAR.
C      B0PRIM .. PRESSURE DERIVATIVE OF THE BULK MODULUS.
C                SHOULD BE CLOSE NEITHER TO 0 NOR TO 1.
C      E0 ...... AN ADDITIVE CONSTANT (IN ELECTRONVOLTS), ADDED
C                TO THE ENERGY EXPRESSION.
C                (SEE,  PR B28, p 5484: Fu and Ho)
C
C OUTPUT DATA:
C
C      ETOT .... THE ENERGY, INCLUDING THE E0 CONSTANT, IN ELECTRONVOLTS
C
C IF B0 DIFFERS FROM 0. OR FROM 1. BY LESS THAN 1.d-6, THEN
C ETOT IS SET AT +111111.111111 ELECTRONVOLTS.
C
C      *      *      *      *      *      *      *      *
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     CONVERSION FACTOR FROM ERG TO ELECTRONVOLT:
      PARAMETER( CONVXX = 1.60209D0 )
C     1 ELECTRONVOLT = 1.60209 D-12 ERG
C
C
      IF(DABS(B0PRIM)-1.D0 .LT. 1.d-6   .OR.
     1                             DABS(B0PRIM) .LT. 1.d-6) THEN
      ETOT=+111111.111111D0
      RETURN
      END IF
C
      IF(VOL .LT. 0.D0   .OR.   VOL0 .LT. 0.D0) THEN
      ETOT=+111111.111111D0
      RETURN
      END IF
C
      ETOT = E0 - B0*VOL0/B0PRIM *
     1 (((VOL/VOL0)**(1.D0-B0PRIM)-1.D0)/(1.D0-B0PRIM)-VOL/VOL0+1.D0)
     2 *ULA**3/CONVXX
C
      RETURN
      END
C ---------------------------------------------------------------------
C  --
C --------
      SUBROUTINE DFMND(F,X,Y,N,LIM,AUX,IER)
C
C     ******************************************************************
C     *   MINIMIZATION OF A FUNCTION OF SEVERAL VARIABLES              *
C     *   USING POWELL'S ALGORITHM WITHOUT DERIVATIVES                 *
C     ******************************************************************
C
      DOUBLE PRECISION F,X,Y,AUX,EPS,ETA,TOL,
     1     DA,DB,DC,DM,DQ,FA,FB,FC,FL,FM,FS,HD,HQ,HX
      DIMENSION X(1),AUX(1)
C
C SUBROUTINES REQUIRED: THE EXTERNAL FUNCTION F.
C
C INPUT DATA:
C     F .... THE FUNCTION OF N VARIABLES X(1)...X(N) TO BE MINIMIZED
C     X .... X(1) ... X(N) = STARTING GUESS FOR THE VARIABLES;
C            BEWARE: X HAS TO BE DIMENSIONED IN THE MAIN PROGRAM
C            TO CONSIDERABLY MORE THAN N.
C     N .... NUMBER OF VARIABLES; THE DIMENSION OF X AND AUX IN THE
C            CALLING PROGRAM MUST BE, HOWEVER, MUCH HIGHER:
C            - PERHAPS 10 TIMES?
C     LIM .. MAXIMUM NUMBER OF ITERATIONS
C     AUX .. AUXILIARY ARRAY OF THE SAME DIMENSION AS X.
C OUTPUT DATA:
C     X .... X(1) ... X(N) THE RESULTING MINIMUM
C     Y .... VALUE OF THE FUNCTION AT THE MINIMUM
C     IER .. SOME ERROR INDICATION 
C            IERR=0 MEANS 'CONVERGENCE ACHIEVED'.
C
C      *      *      *      *      *      *      *      *
C
      ISW  =IER
      IER  =0
      IF (N) 1,1,2
    1 IER  =1000
      GOTO 109
    2 IF (LIM) 3,3,4
    3 IER  =2000
      GOTO 109
C
C     SET INITIAL VALUES AND SAVE INITIAL ARGUMENT
C
    4 N1   =N+1
      N2   =N+N
      N3   =N+N2
      N4   =N+N3
      N5   =N*N+N4
      EPS  =1.D-15
      ETA  =N*EPS
      DO 5 K=1,N
         AUX(K)=X(K)
         J    =N3+K
    5    AUX(J)=1.D0
      FS   =F(AUX)
      FB   =FS
      I    =1
      IC   =1
      IT   =0
      M    =N4
      MF   =0
      IS   =1
C
C     START ITERATION CYCLE
C
    6 IT   =IT+1
      FL   =FS
      DO 7 K=N1,N2
    7    AUX(K)=0.D0
      ID   =I
      I    =IC
      IW   =0
C
C     START MINIMIZATION ALONG NEXT DIRECTION
C
    8 NS   =0
      IP   =0
      DB   =0.D0
      IF (IW) 10,9,10
    9 HX   =AUX(I)
   10 IF (IT-1) 11,11,14
   11 IF (IS-1) 14,12,14
   12 DM   =.1D0
      IF (DABS(HX)-1.D0) 38,38,13
   13 DM   =-DM*HX
      GOTO 38
   14 IF (IS-2) 18,15,18
   15 IF (IT-1) 17,16,17
   16 DM   =HQ
      GOTO 38
   17 DM   =DQ
      GOTO 38
C
C     INTERPOLATION USING ESTIMATE OF SECOND DERIVATIVE
C
   18 IF (IW-1) 20,19,20
   19 J    =N2+I
      GOTO 21
   20 J    =N3+I
   21 HD   =AUX(J)
      DC   =1.D-2
      IF (IT-2) 23,23,22
   22 DC   =HQ
   23 DM   =DC
      MK   =1
      GOTO 51
   24 DM   =DC*HD
      IF (DM) 26,25,26
   25 DM   =1.D0
   26 DM   =.5D0*DC-(FM-FB)/DM
      MK   =2
      IF (FM-FB) 27,29,29
   27 FC   =FB
      FB   =FM
      DB   =DC
      IF (DM-DB) 28,67,28
   28 DC   =0.D0
      GOTO 51
   29 IF (DM-DB) 31,30,31
   30 DA   =DC
      FA   =FM
      GOTO 37
   31 FC   =FM
      GOTO 51
C
C     ANALYSE INTERPOLATED FUNCTION VALUE
C
   32 IF (FM-FB) 34,33,33
   33 DA   =DM
      FA   =FM
      GOTO 35
   34 DA   =DB
      FA   =FB
      DB   =DM
      FB   =FM
   35 IF ((DC-DA)/(DB-DA)) 36,36,50
   36 IF (DB) 67,37,67
   37 NS   =1
      DM   =-DC
C
C     LINEAR SEARCH FOR SMALLER FUNCTION VALUES
C     ALONG CURRENT DIRECTION
C
   38 IF (NS-15) 43,43,39
   39 IF (FS-FM) 41,40,41
   40 MF   =N+2
      DB   =0.D0
      GOTO 67
   41 IF (DABS(DM)-1.D6) 43,43,42
   42 IER  =100
      GOTO 67
   43 NS   =NS+1
      MK   =3
      GOTO 51
   44 IF (FM-FB) 45,46,47
   45 DA   =DB
      FA   =FB
      DB   =DM
      FB   =FM
      DM   =DM+DM
      GOTO 38
   46 IF (FS-FB) 47,45,47
   47 IF (NS-1) 48,48,49
   48 DA   =DM
      FA   =FM
      DM   =-DM
      GOTO 38
   49 DC   =DM
      FC   =FM
C
C     REFINE MINIMUM USING QUADRATIC INTERPOLATION
C
   50 HD   =(FC-FB)/(DC-DB)+(FA-FB)/(DB-DA)
      DM   =.5D0*(DA+DC)+(FA-FC)/(HD+HD)
      IP   =IP+1
      MK   =4
C
C     STEP ARGUMENT VECTOR AND CALCULATE FUNCTION VALUE
C
   51 IF (IW-1) 54,52,54
   52 DO 53 K=1,N
         L    =M+K
   53    AUX(K)=X(K)+DM*AUX(L)
      GOTO 55
   54 AUX(I)=HX+DM
   55 FM   =F(AUX)
      GOTO (24,32,44,56),MK
C
C     ANALYSE INTERPOLATED FUNCTION VALUE
C
   56 IF (FM-FB) 61,61,57
   57 IF (IP-3) 58,62,62
   58 IF ((DC-DB)/(DM-DB)) 60,60,59
   59 DC   =DM
      FC   =FM
      GOTO 50
   60 DA   =DM
      FA   =FM
      GOTO 50
   61 DB   =DM
      FB   =FM
C
C     CALCULATE NEW ESTIMATE OF SECOND DERIVATIVE
C     ALONG THE CURRENT DIRECTION
C
   62 HD   =(HD+HD)/(DC-DA)
      IF (IW-1) 64,63,64
   63 J    =N2+I
      GOTO 65
   64 J    =N3+I
   65 AUX(J)=HD
      IF (FB-FS) 67,66,67
   66 DB   =0.D0
C
C     SAVE ARGUMENT VECTOR WITH SMALLEST FUNCTION VALUE FOUND
C
   67 IF (IW-1) 70,68,70
   68 DO 69 K=1,N
         L    =M+K
         J    =N+K
         HD   =DB*AUX(L)
         AUX(J)=AUX(J)+HD
         HD   =X(K)+HD
         AUX(K)=HD
   69    X(K) =HD
      GOTO 71
   70 J    =N+I
      AUX(J)=AUX(J)+DB
      HD   =HX+DB
      AUX(I)=HD
      X(I) =HD
   71 IF (IER-100) 72,108,72
C
C     DETERMINE DIRECTION FOR NEXT LINEAR SEARCH
C
   72 FS   =FB
      IF (I-N) 74,73,73
   73 I    =0
   74 I    =I+1
      IF (IS) 75,75,80
   75 IF (DB) 77,76,77
   76 IF (I-IC) 8,77,8
   77 IC   =I
      IS   =1
      IF (IT-N) 79,79,78
   78 IW   =1
   79 I    =ID
      GOTO 8
   80 M    =M+N
      IF (M-N5) 82,81,81
   81 M    =N4
   82 IF (IS-1) 83,83,94
   83 IF (I-1) 84,84,85
   84 IW   =1
   85 IF (I-ID) 8,86,8
   86 HQ   =0.D0
      DO 87 K=N1,N2
   87    HQ   =HQ+AUX(K)*AUX(K)
      IF (HQ) 90,88,90
   88 IF (MF-N1) 108,108,89
   89 IER  =200
      GOTO 108
   90 DQ   =DSQRT(HQ)
      HQ   =DQ
      IF (HQ-1.D0) 92,92,91
   91 HQ   =1.D0
   92 DO 93 K=N1,N2
         L    =M+K-N
   93    AUX(L)=AUX(K)/DQ
      IS   =2
      GOTO 8
C
C     END OF ITERATION CYCLE
C     TEST FOR TERMINATION OF MINIMIZATION
C
   94 IS   =0
      TOL  =EPS
      IF (DABS(FS)-1.D0) 96,96,95
   95 TOL  =EPS*DABS(FS)
   96 IF (FL-FS-TOL) 100,100,97
   97 IF (IT-LIM) 99,99,98
   98 IER  =10
      GOTO 108
   99 MF   =0
      GOTO 6
  100 IF (MF-N1) 102,102,101
  101 IER  =200
      GOTO 108
  102 MF   =MF+1
      DQ   =0.D0
      DO 105 K=1,N
         J    =N+K
         IF (DABS(AUX(K))-1.D0) 103,103,104
  103    DQ   =DQ+DABS(AUX(J))
         GOTO 105
  104    DQ   =DQ+DABS(AUX(J)/AUX(K))
  105    CONTINUE
      IF (DQ-ETA) 108,108,106
  106 IF (MF-N) 6,6,107
  107 IER  =1
  108 Y    =FB
      IF (IER) 111,111,109
  109 IF (ISW+12345) 110,111,110
C 110 CALL WIER(IER,20212)
  110 CONTINUE
  111 RETURN
      END
