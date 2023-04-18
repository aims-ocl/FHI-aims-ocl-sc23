!****s* FHI-aims/confined_out
!  NAME
!    confined_out
!  SYNOPSIS

subroutine confined_out &
     ( i_species, n_grid_co, n_wave, wave_n, wave_l, &
     r_grid_co, pot, wave, eigenval, descript, &
     dim_species, dim_grid, dim_fn )

!  PURPOSE
!  Subroutine confined_out outputs potential  and wave function data
!  for one individual atom.
!
!  It's a slightly modified version of ion_out and should one day replace it entirely!
!
!  USES

  use mpi_tasks,   only : myid
  use localorb_io, only : use_unit, localorb_info
  implicit none

!  ARGUMENTS

  integer dim_species, dim_grid, dim_fn
  integer i_species
  integer n_grid_co
  integer n_wave(dim_species)
  integer wave_n(dim_species, dim_fn)
  integer wave_l(dim_species, dim_fn)
  real*8  r_grid_co (n_grid_co)
  real*8  pot  (dim_grid, dim_species, dim_fn)
  real*8  wave (dim_grid, dim_species, dim_fn)
  real*8  eigenval (dim_species,dim_fn)
  character*4 descript
       

!  INPUTS
!   o dim_species, dim_grid, dim_fn -- dimensions of the tables
!   o i_species -- ion number in question
!   o n_grid_co --  number of radial grid points
!   o n_wave -- number of wave functions
!   o wave_n -- n index of basis functions
!   o wave_l -- l index of basis functions
!   o r_grid_co -- ????????????
!   o pot -- effective wave-function-defining potential, tabulated on grid
!   o rho -- charge density of occupied free orbitals, tabulated on grid
!   o wave -- array of wave functions, tabulated on grid, for each main quantum number n
!     and angular momentum quantum number 
!   o descript -- file descript where the results are written out.
!
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
!






!  local variables

       character*21 pot_file, wave_file
       character l_char

!  counters

       integer i_grid, i_wave

!  functions

       character l_to_str

!  begin work


       if (myid.eq.0) then
          write(use_unit,*) " Writing confined data for species ", &
               i_species, "."
       end if

          do i_wave = 1, n_wave(i_species), 1

             l_char = l_to_str(wave_l(i_species,i_wave))

!     Potential:
             if (myid.eq.0) then
                if (i_wave.le.9) then
                   write (pot_file,'(A1,I1,A1,I1,A1,I1,A1,A1,A4,A8)') &
                        "C",i_species, "_", i_wave, "_", &
                        wave_n(i_species,i_wave), &
                        l_char, "_", descript, "_pot.dat"
                else
                   write (pot_file,'(A1,I1,A1,I2,A1,I1,A1,A1,A4,A8)') &
                        "C",i_species, "_", i_wave, "_", &
                        wave_n(i_species,i_wave), &
                        l_char, "_", descript, "_pot.dat"
                end if

                write(use_unit,*) " | Potential: ", pot_file

                open(50, FILE=pot_file)
                write (50,*) "# ", n_grid_co
                do i_grid = 1, n_grid_co, 1
                   write(50,*) r_grid_co(i_grid), &
                        pot(i_grid,i_species,i_wave)
                enddo
                close(50)

!       Wave functions:

                if (i_wave.le.9) then
                   if (wave_n(i_species,i_wave).le.9) then
                      write (wave_file, &
                           '(A1,I1,A1,I1,A1,I1,A1,A1,A4,A4)') &
                           "C",i_species, "_", i_wave, "_", &
                           wave_n(i_species,i_wave), &
                           l_char, "_", descript, ".dat"
                   else
                      write (wave_file, &
                           '(A1,I1,A1,I1,A1,I2,A1,A1,A4,A4)') &
                           "C",i_species, "_", i_wave, "_", &
                           wave_n(i_species,i_wave), &
                           l_char, "_", descript, ".dat"
                   end if
                else
                   if (wave_n(i_species,i_wave).le.9) then
                      write (wave_file, &
                           '(A1,I1,A1,I2,A1,I1,A1,A1,A4,A4)') &
                           "C",i_species, "_", i_wave, "_", &
                           wave_n(i_species,i_wave), &
                           l_char, "_", descript, ".dat"
                   else
                      write (wave_file, &
                           '(A1,I1,A1,I2,A1,I2,A1,A1,A4,A4)') &
                           "C",i_species, "_", i_wave, "_", &
                           wave_n(i_species,i_wave), &
                           l_char, "_", descript, ".dat"
                   end if
                end if

                write(use_unit,*) " | Atomic shell ", wave_n(i_species,i_wave), &
                     l_char, ": ", wave_file
                write(use_unit,*) " | Eigenvalue: ", eigenval(i_species,i_wave)

                open (50, FILE=wave_file)
                write (50,*) "# ", n_grid_co
                do i_grid = 1, n_grid_co, 1
                   write(50,*) r_grid_co(i_grid), &
                        wave(i_grid, i_species, i_wave)
                enddo
                close(50)

             end if

          end do

          call localorb_info('',use_unit)

!  that's all folks

      return
    end subroutine confined_out
!******
