!****s* FHI-aims/output_rho_multipole
!  NAME
!   output_rho_multipole
!  SYNOPSIS

      subroutine output_rho_multipole

!  PURPOSE
!  Subroutine output_rho_multipole writes the multipole components of the
!  partitioned density [rho-rho_free]
!  at each atom into individual files!
!
!  USES

      use dimensions
      use runtime_choices
      use grids
      use geometry
      use species_data
      use spline
      use synchronize_mpi
      use hartree_potential_storage, only : get_rho_multipole_spl
      use localorb_io, only : use_unit
      implicit none

!  ARGUMENTS

!  INPUTS
!    none
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





!     local variables

      character*34 rho_file

      real*8 i_r
      real*8 pot_at_lm

      real*8, dimension(:,:,:), allocatable :: current_rho_multipole_spl

!     counters

      integer i_atom, i_l, i_m
      integer i_grid, i_radial

      integer index_lm

!     begin work


!     allocate temporary spline array
      if (.not.allocated(current_rho_multipole_spl)) then
         allocate(current_rho_multipole_spl &
              ((l_pot_max+1)**2, n_max_spline, n_max_radial+2))
      end if

      if (myid.eq.0) then
         write(use_unit,'(2X,A,A)') &
              "Writing partitioned Hartree potentials (at,l,m) ", &
              "to files rho_multipole.*.dat ."
      end if

      do i_atom = 1, n_atoms, 1

        index_lm = 0
        do i_l = 0, l_hartree(species(i_atom)),1
          do i_m = -i_l, i_l, 1
            index_lm = index_lm+1

!           following is a long series of conditions to determine the
!           proper format for the file name.
            if (i_m .ge.0) then
              if (i_atom.le.9) then
                if (i_l.le.9) then
                  write (rho_file,'(A17,I1,A4,I1,A1,I1,A4)') &
                  "rho_multipole.at_", i_atom, &
                  ".lm_", i_l, "_", i_m,".dat"
                else if (i_l.le.99) then
                  if (i_m.le.9) then
                    write (rho_file,'(A17,I1,A4,I2,A1,I1,A4)') &
                    "rho_multipole.at_", i_atom, &
                    ".lm_", i_l, "_", i_m, &
                    ".dat"
                  else if (i_m.le.99) then
                    write (rho_file,'(A17,I1,A4,I2,A1,I2,A4)') &
                    "rho_multipole.at_", i_atom, &
                    ".lm_", i_l, "_", i_m, &
                    ".dat"
                  end if
                end if
              else if (i_atom.le.99) then
                if (i_l.le.9) then
                  write (rho_file,'(A17,I2,A4,I1,A1,I1,A4)') &
                  "rho_multipole.at_", i_atom, &
                  ".lm_", i_l, "_", i_m,".dat"
                else if (i_l.le.99) then
                  if (i_m.le.9) then
                    write (rho_file,'(A17,I2,A4,I2,A1,I1,A4)') &
                    "rho_multipole.at_", i_atom, &
                    ".lm_", i_l, "_", i_m, &
                    ".dat"
                  else if (i_m.le.99) then
                    write (rho_file,'(A17,I2,A4,I2,A1,I2,A4)') &
                    "rho_multipole.at_", i_atom, &
                    ".lm_", i_l, "_", i_m, &
                    ".dat"
                  end if
                end if
              end if
            else
              if (i_atom.le.9) then
                if (i_l.le.9) then
                  write (rho_file,'(A17,I1,A4,I1,A1,I2,A4)') &
                  "rho_multipole.at_", i_atom, &
                  ".lm_", i_l, "_", i_m,".dat"
                else if (i_l.le.99) then
                  if (i_m.ge.-9) then
                    write (rho_file,'(A17,I1,A4,I2,A1,I2,A4)') &
                    "rho_multipole.at_", i_atom, &
                    ".lm_", i_l, "_", i_m, &
                    ".dat"
                  else if (i_m.ge.-99) then
                    write (rho_file,'(A17,I1,A4,I2,A1,I3,A4)') &
                    "rho_multipole.at_", i_atom, &
                    ".lm_", i_l, "_", i_m, &
                    ".dat"
                  end if
                end if
              else if (i_atom.le.99) then
                if (i_l.le.9) then
                  write (rho_file,'(A17,I2,A4,I1,A1,I2,A4)') &
                  "rho_multipole.at_", i_atom, &
                  ".lm_", i_l, "_", i_m,".dat"
                else if (i_l.le.99) then
                  if (i_m.ge.-9) then
                    write (rho_file,'(A17,I2,A4,I2,A1,I2,A4)') &
                    "rho_multipole.at_", i_atom, &
                    ".lm_", i_l, "_", i_m, &
                    ".dat"
                  else if (i_m.ge.-99) then
                    write (rho_file,'(A17,I2,A4,I2,A1,I3,A4)') &
                    "rho_multipole.at_", i_atom, &
                    ".lm_", i_l, "_", i_m, &
                    ".dat"
                  end if
                end if
              end if
            end if

            if (myid.eq.0) then
               open (50, file=rho_file)
            end if

            call get_rho_multipole_spl(current_rho_multipole_spl, i_atom)

                 do i_radial = 1, n_radial(species(i_atom)), 1

                   if (myid.eq.0) then
                     write(50,*) r_radial(i_radial,species(i_atom)), &
                        current_rho_multipole_spl(index_lm,1,i_radial+1)
                   end if

                 enddo

                 if (myid.eq.0) then
                    close(50)
                 end if

          enddo
        enddo

      enddo

      if (allocated(current_rho_multipole_spl)) then
         deallocate(current_rho_multipole_spl)
      end if

      end subroutine output_rho_multipole
!------------------------------------------------------------------------
!******	
