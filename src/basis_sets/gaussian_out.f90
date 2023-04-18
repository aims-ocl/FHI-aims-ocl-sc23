!****s* FHI-aims/gaussian_out
!  NAME
!   gaussian_out
!  SYNOPSIS
      subroutine gaussian_out &
        ( i_species, n_grid_go, n_wave, wave_n, wave_l, &
          r_grid_go, wave, kinetic, descript, &
          dim_species, dim_grid, dim_fn )
!  PURPOSE
!    Subroutine gaussian_out outputs potential  and wave function data
!    for one individual atom.
!  USES
      use mpi_tasks,   only : myid
      use localorb_io, only : use_unit, localorb_info
      implicit none
!  ARGUMENTS
      integer     :: dim_species, dim_grid, dim_fn
      integer     :: i_species
      integer     :: n_grid_go
      integer     :: n_wave (dim_species)
      integer     :: wave_n (dim_species, dim_fn)
      integer     :: wave_l (dim_species, dim_fn)
      real*8      :: r_grid_go (n_grid_go)
      real*8      :: wave   (dim_grid, dim_species, dim_fn)
      real*8      :: kinetic(dim_grid, dim_species, dim_fn)
      character*4 :: descript
!  INPUTS
!   o i_species -- ion number in question
!   o n_grid_go -- number of radial grid points
!   o n_wave -- number of basis fns
!   o wave_l -- angular quantum number
!   o r_grid_go -- radial grid values
!   o pot -- effective wave-function-defining potential, tabulated on grid
!   o rho -- charge density of occupied free orbitals, tabulated on grid
!   o wave -- array of wave functions, tabulated on grid, for each main quantum number n
!          and angular momentum quantum number l
!   o descript -- descriptor for output file
!  OUTPUTS
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
       character*21 :: kin_file, wave_file
       character    :: l_char
       integer      :: i_grid, i_wave
       character    :: l_to_str


       if (myid.eq.0) then
          write(use_unit,*) " Writing Gaussian fns. for species ", &
               i_species,"."
       end if

      do i_wave = 1, n_wave(i_species), 1

        l_char = l_to_str(wave_l(i_species,i_wave))

        if (myid.eq.0) then
           if (i_wave.le.9) then
              write (kin_file,'(A1,I1,A1,I1,A1,I1,A1,A1,A4,A8)') &
                   "G",i_species, "_", i_wave, "_", &
                   wave_n(i_species,i_wave), &
                   l_char, "_", descript, "_kin.dat"
           else
              write (kin_file,'(A1,I1,A1,I2,A1,I1,A1,A1,A4,A8)') &
                   "G",i_species, "_", i_wave, "_", &
                   wave_n(i_species,i_wave), &
                   l_char, "_", descript, "_kin.dat"
           end if

           write(use_unit,*) " | Kinetic: ", kin_file

           open(50, FILE=kin_file)
           write (50,*) "# ", n_grid_go
           do i_grid = 1, n_grid_go, 1
              write(50,*) r_grid_go(i_grid), &
                   kinetic(i_grid,i_species,i_wave)
           enddo
           close(50)

!     Wave functions:

           if (i_wave.le.9) then
              write (wave_file, '(A1,I1,A1,I1,A1,I1,A1,A1,A4,A4)') &
                   "G",i_species, "_", i_wave, "_", &
                   wave_n(i_species,i_wave), &
                   l_char, "_", descript, ".dat"
           else
              write (wave_file, '(A1,I1,A1,I2,A1,I1,A1,A1,A4,A4)') &
                   "G",i_species, "_", i_wave, "_", &
                   wave_n(i_species,i_wave), &
                   l_char, "_", descript, ".dat"
           end if

           write(use_unit,*) " | Gaussian L = ", wave_n(i_species,i_wave), &
                ", l = ", l_char, ": ", wave_file

           open (50, FILE=wave_file)
           write (50,*) "# ", n_grid_go
           do i_grid = 1, n_grid_go, 1
              write(50,*) r_grid_go(i_grid), &
                   wave(i_grid, i_species, i_wave)
           enddo
           close(50)
        end if

      end do

      call localorb_info('',use_unit)

!  that's all folks

      return
    end subroutine gaussian_out
!******
