!****s* FHI-aims/atom_out
!  NAME
!    atom_out
!  SYNOPSIS

      subroutine atom_out &
        ( name, n_grid_ao, r_grid_ao, pot, rho, n_wave, wave, wave_n, &
           wave_l, eigenval, descript, dim_grid, dim_fns )

!  PURPOSE
!
!  Subroutine atom_out outputs potential, density, and wave function data
!  for one individual atom.
!
!  Notice that the variables here have different names than in the subroutines
!  where they came from!

!  USES

      use mpi_tasks,   only : myid
      use localorb_io, only : use_unit, localorb_info
      implicit none

!  ARGUMENTS

      integer dim_grid, dim_fns
       character*2 name
       integer n_grid_ao
       real*8  r_grid_ao (dim_grid)
       real*8  pot  (dim_grid)
       real*8  rho  (dim_grid)
       integer n_wave
       real*8  wave (dim_grid,dim_fns)
       integer  wave_n (dim_fns)
       integer  wave_l (dim_fns)
       real*8  eigenval (dim_fns)
       character*4 descript

!  INPUTS
!
!  o dim_grid -- number of grid points
!  o dim_fns -- number of fns points
!  o name -- name of atom type printed out
!  o n_grid_ao -- ????????
!  o r_grid_ao -- ????????
!  o pot -- potential
!  o rho -- electron density
!  o n_wave -- number of basis functions
!  o wave -- basis functions
!  o wave_n -- n index of the basis functions
!  o wave_l -- l index of the basis functions
!  o eigenval -- Kohn-Sham eigenvalues
!  o descript -- file name descriptor for output
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

       character*15 pot_file, rho_file, wave_file
       character l_char

!  counters

       integer i_grid, i_fn

!  functions

       character l_to_str

!test
!       write(use_unit,*) wave_n(1), wave_l(1)
!       write(use_unit,*) wave_n(2), wave_l(2)
!       write(use_unit,*) wave_n(3), wave_l(3)
!test end

!  begin work

       if (myid.eq.0) then
          write(use_unit,*) " Writing atomic data for atom type ", name, "."
       end if

!  Potential:

       if (myid.eq.0) then
          if (name(2:2).eq.' ') then
             write (pot_file,'(A1,A1,A4,A8)') &
                  name, "_", descript, "_pot.dat"
          else
             write (pot_file,'(A2,A1,A4,A8)') &
                  name, "_", descript, "_pot.dat"
          end if

          write(use_unit,*) " | Potential: ", pot_file

          open(50, FILE=pot_file)
          write (50,*) "# ", n_grid_ao
          do i_grid = 1, n_grid_ao, 1
             write(50,*) r_grid_ao(i_grid), pot(i_grid)
          enddo
          close(50)

!  Charge density:

          if (name(2:2).eq.' ') then
             write (rho_file,'(A1,A1,A4,A8)') &
                  name, "_", descript, "_rho.dat"
          else
             write (rho_file,'(A2,A1,A4,A8)') &
                  name, "_", descript, "_rho.dat"
          end if

          write(use_unit,*) " | Density: ", rho_file

          open(50, FILE=rho_file)
          write (50,*) "# ", n_grid_ao
          do i_grid = 1, n_grid_ao, 1
             write(50,*) r_grid_ao(i_grid), rho(i_grid)
          enddo
          close(50)
       end if

!  Wave functions:

      do i_fn = 1, n_wave, 1

        l_char = l_to_str(wave_l(i_fn))

        if (myid.eq.0) then
           if (name(2:2).eq.' ') then
              write (wave_file, '(A1,A1,A4,A1,I1,A1,A1,A4)') &
                   name, "_", descript, "_", wave_n(i_fn), "_", &
                   l_char, ".dat"
           else
              write (wave_file, '(A2,A1,A4,A1,I1,A1,A1,A4)') &
                   name, "_", descript, "_", wave_n(i_fn), "_", &
                   l_char, ".dat"
           end if

           write(use_unit,*) " | Atomic shell ", wave_n(i_fn), l_char, ": ", &
                wave_file
           write(use_unit,*) " | Eigenvalue: ", eigenval(i_fn)

           open (50, FILE=wave_file)
           write (50,*) "# ", n_grid_ao
           do i_grid = 1, n_grid_ao, 1
              write(50,*) r_grid_ao(i_grid), &
                   wave(i_grid,i_fn)
           enddo
           close(50)
        end if

      enddo

      call localorb_info('',use_unit)

!  that's all folks

      return
    end subroutine atom_out
!******
