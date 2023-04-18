!****s* FHI-aims/output_cube_files
!  NAME
!  output_cube_files
!  SYNOPSIS
   subroutine output_cube_files ()
!  PURPOSE
!  Plot charge density for clusters. This routine can still
!  be improved, see FIXME's below.
!
!  USES
!
      use dimensions
      use runtime_choices
      use grids
      use geometry
      use species_data
      use mpi_utilities
      use localorb_io
      use basis
      use cartesian_ylm
      use constants
      use plot
      use physics
!  INPUTS
!   none
!  OUTPUT
!   Writes charge density to file
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



! TODO: figure out automatically which region to plot
! FIXME: ATM all plots use the same grid  (grid of the first plot)
! FIXME: ATM only canonical cartesian basis is supported


      implicit none

!  local variables
      real*8 inv_bohr_3
      real*8 sqrt_inv_bohr_3

      real*8 coord_current(3)
! FIXME for support for different output grids
      real*8 dist_tab(n_atoms, cube_edge_steps(3,1) )
      real*8 dist_tab_sq(n_atoms, cube_edge_steps(3,1) )
      real*8 i_r(n_atoms, cube_edge_steps(3,1) )
      real*8 dir_tab(3, n_atoms, cube_edge_steps(3,1) )
      real*8 trigonom_tab(4, n_atoms, cube_edge_steps(3,1) )
      real*8 radial_wave(n_basis, cube_edge_steps(3,1) )
      real*8 wave(n_basis, cube_edge_steps(3,1) )

      integer :: n_compute
      integer :: i_basis(n_basis)
      integer :: i_cube

      integer :: n_compute_fns
      integer :: i_basis_fns(n_basis_fns*n_atoms)
      integer :: i_basis_fns_inv(n_basis_fns,n_atoms)
      integer :: i_atom_fns(n_basis_fns*n_atoms)

      integer :: n_compute_atoms
      integer :: atom_index(n_atoms)
      integer :: atom_index_inv(n_atoms)

      integer :: spline_array_start(n_atoms)
      integer :: spline_array_end(n_atoms)

!     other local variables
      integer, dimension(n_spin) :: max_occ_number
      real*8, dimension(n_states,n_spin) :: occ_numbers_sqrt

      integer :: l_ylm_max
      integer :: n_points

      real*8, dimension(n_basis, n_states, n_spin) :: &
        KS_vec_times_occ_sqrt

      real*8, dimension(n_states, n_basis, n_spin) :: &
        KS_ev_compute

      real*8, dimension(:,:,:), allocatable ::  KS_orbital
      real*8, dimension(:,:), allocatable :: local_rho
      real*8, dimension(:,:), allocatable ::  local_orb_up
      real*8, dimension(:,:), allocatable ::  local_orb_down

      integer, dimension(:,:), allocatable :: index_lm
      real*8, dimension(:,:,:), allocatable :: ylm_tab

      real*8, dimension(:,:,:), allocatable :: dir_tab_global
      real*8, dimension(:,:), allocatable :: dist_tab_sq_global

      real*8 offset_coord(3)
      real*8 cube_units(3)

!      character*40 cube_filename 
!NOTE: This was commented out to allow the user setting his own name via 
!other routine. Since this routine is never used in the code anymore
!it might break the code if, for some reason, it is ever included again
!However, this should be VERY unlikly!
      character*40 cube_filename_spin_down

!     counters
      integer :: coord_x, coord_y, coord_z
      integer :: i_l
      integer :: i_m
      integer :: i_state
      integer :: i_point

      integer :: i_full_points

      integer :: i_spin = 1
      integer :: i_index

      inv_bohr_3 = 1.0d0/(bohr**3)
      sqrt_inv_bohr_3 = sqrt(inv_bohr_3)

!     begin work

      if (myid.eq.0) then

! FIXME for support of non-cartesian grids
        cube_units(1)= cube_edge_unit(1,1,1)
        cube_units(2)= cube_edge_unit(2,2,1)
        cube_units(3)= cube_edge_unit(3,3,1)

        offset_coord = real(cube_edge_steps(1:3,1)-1)/2.d0* &
            cube_units

        write(use_unit,'(A)') &
          "------------------------------------------------------------"
        write(use_unit,'(2X,A)') "Writing density files:"
        do i_cube = 1, n_cube, 1
            if (cube_type(i_cube).eq.'total_density') then
             cube_filename ="total_density.cube"
          else if (cube_type(i_cube).eq.'spin_density') then
             cube_filename ="spin_density.cube"
            else if (cube_type(i_cube).eq.'eigenstate' &
                  .AND. spin_treatment .ne. 1) then
             write (unit=cube_filename,fmt='(A11,I5.5,A)') 'eigenstate_' &
                ,cube_index(i_cube),".cube"
          else if (cube_type(i_cube).eq.'eigenstate' &
                  .AND. spin_treatment .eq. 1) then
           write (unit=cube_filename,fmt='(A19,I5.5,A)') &
                  'spin_up_eigenstate_' &
                ,cube_index(i_cube),".cube"
            write (unit=cube_filename_spin_down,fmt='(A21,I5.5,A)') &
                   'spin_down_eigenstate_',cube_index(i_cube),".cube"
            else if (cube_type(i_cube).eq.'eigenstate_density' &
                  .and. spin_treatment .ne. 1) then
             write (unit=cube_filename,fmt='(A19,I5.5,A)') &
                'eigenstate_density_',cube_index(i_cube),".cube"
          else if (cube_type(i_cube).eq.'eigenstate_density' &
                  .and. spin_treatment .eq. 1) then
           write (unit=cube_filename,fmt='(A27,I5.5,A)') &
                  'spin_up_eigenstate_density_' &
                ,cube_index(i_cube),".cube"
            write (unit=cube_filename_spin_down,fmt='(A29,I5.5,A)') &
                   'spin_down_eigenstate_density_' &
                  ,cube_index(i_cube),".cube"
            endif

          if (spin_treatment .ne. 1 .or. &
                  (cube_type(i_cube) .eq. 'total_density' &
                  .or. cube_type(i_cube) .eq. 'spin_density')) then
                  write(use_unit,'(2X,A,1X,A)') &
                   "| Creating cube file   :", cube_filename
                  open(unit= 10+i_cube, file=cube_filename)
                  call write_cube_header(10+i_cube,offset_coord &
                  - cube_origin(1:3,1),cube_edge_steps(1:3,1), cube_units)
          else
            write(use_unit,'(2X,A,1X,A)') &
                   "| Creating cube file   :", cube_filename
                  open(unit= 10+i_cube, file=cube_filename)
                  call write_cube_header(10+i_cube,offset_coord &
                  - cube_origin(1:3,1),cube_edge_steps(1:3,1), cube_units)
            write(use_unit,'(2X,A,1X,A)') &
                   "| Creating cube file   :", cube_filename_spin_down
                  open(unit= 100+i_cube, file=cube_filename_spin_down)
                  call write_cube_header(100+i_cube,offset_coord &
                  - cube_origin(1:3,1),cube_edge_steps(1:3,1), cube_units)
          endif
        enddo

        if (n_k_points.ne.1) then
            write(use_unit,*) "WARNING: cube output does not support more than ",&
            "a single k-point."
        endif

        l_ylm_max = l_wave_max

        allocate( ylm_tab( (l_ylm_max+1)**2, n_atoms, &
         cube_edge_steps(3,1) ) )
        allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max) )

!     initialize index_lm
        i_index = 0
        do i_l = 0, l_ylm_max, 1
            do i_m = -i_l, i_l
             i_index = i_index + 1
             index_lm(i_m, i_l) = i_index
            enddo
        enddo

!     find the maximal occupation number
        do i_spin = 1, n_spin, 1
        ! initialize
            max_occ_number(i_spin) = 0
            do i_state = n_states, 1, -1
             if (dabs(occ_numbers(i_state,i_spin,1)).gt.0.d0) then
              max_occ_number(i_spin) = i_state
              exit
             endif
            enddo
        enddo

        allocate( local_rho(1:cube_edge_steps(3,1),1:2) )
        allocate(local_orb_up(max_occ_number(1), &
            1:cube_edge_steps(3,1)))
        local_orb_up = 0.d0
      if (n_spin .gt. 1) then
              allocate(local_orb_down(max_occ_number(2), &
                        1:cube_edge_steps(3,1)))
              local_orb_down = 0.d0
      endif
        allocate( KS_orbital(n_states,1:cube_edge_steps(3,1),1:2) )
        KS_orbital = 0.d0

!     allocate the sqrt array for occupation numbers and fill it
!     up to max_occ_number

      do i_spin = 1, n_spin, 1
        do i_state = 1, max_occ_number(i_spin), 1

          ! JW: WARNING: I very much doubt that the kweights are right
          ! here.  But as this whole subroutine is not used anyway...

          occ_numbers_sqrt(i_state, i_spin) = &
            sqrt(occ_numbers(i_state,i_spin,1))

          KS_vec_times_occ_sqrt(:,i_state,i_spin) = &
          KS_eigenvector(:,i_state,i_spin,1) * &
          occ_numbers_sqrt(i_state,i_spin)

        enddo
      enddo


! ************************************************
! ************************************************
      i_full_points = 0

      do coord_x =   1,cube_edge_steps(1,1),1
        do coord_y = 1,cube_edge_steps(2,1),1
           n_compute = 0
           i_basis = 0
           i_point = 0

           local_rho = 0.d0

           do coord_z = 1,cube_edge_steps(3,1),1
            i_point = i_point+1
!    generate output grid
              coord_current(1)=cube_units(1) &
                *(coord_x-1)-offset_coord(1)+cube_origin(1,1)
              coord_current(2)=cube_units(2) &
                *(coord_y-1)-offset_coord(2)+cube_origin(2,1)
              coord_current(3)=cube_units(3) &
                *(coord_z-1)-offset_coord(3)+cube_origin(3,1)


!     compute atom-centered coordinates of current integration point as viewed from all atoms
                 call tab_atom_centered_coords_v2 &
                      ( coord_current,dist_tab_sq(1,i_point), &
                      dir_tab(1,1,i_point))

!     determine which basis functions are relevant at current grid point,
!     and tabulate their indices

!     next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
                 call prune_basis_v2 &
                      (dist_tab_sq(1,i_point), n_compute, i_basis)

!          end loop over the z component
           enddo

           n_points = i_point
           if (n_compute.gt.0) then

             ! Determine all radial functions, ylm functions and their derivatives that
             ! are best evaluated strictly locally at each individual grid point.
             i_point = 0
             do coord_z = 1,cube_edge_steps(3,1),1
                   i_point = i_point+1
                   n_compute_atoms = 0
                   n_compute_fns = 0
                   i_basis_fns_inv = 0

                ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                ! Are stored in a compact spline array that can be accessed by spline_vector_waves,
                ! without any copying and without doing any unnecessary operations.
                ! The price is that the interface is no longer explicit in terms of physical
                ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.
               call prune_radial_basis_v2 &
                 ( dist_tab_sq(1,i_point), &
                   n_compute_atoms, atom_index, atom_index_inv, &
                   n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                   i_atom_fns, spline_array_start, spline_array_end )

               ! Tabulate distances, unit vectors, and inverse logarithmic grid units
               ! for all atoms which are actually relevant
               call tab_local_geometry &
               ( dist_tab_sq(1, i_point), n_compute_atoms, atom_index, &
                 dir_tab(1,1,i_point), dist_tab(1,i_point), &
                 i_r(1,i_point) &
               )

               ! Now evaluate radial functions u(r) from the previously stored compressed spline arrays
               call evaluate_radial_functions_v2 &
               (   spline_array_start, spline_array_end, &
                   n_compute_atoms, n_compute_fns, &
                   dist_tab(1,i_point), i_r(1,i_point), &
                   atom_index, i_basis_fns_inv, &
                   basis_wave_ordered, radial_wave(1,i_point), &
                   .false. &
                    )

                 call tab_trigonom_v2 &
                       ( n_compute_atoms, dir_tab(1,1,i_point), &
                       trigonom_tab(1,1,i_point) &
                       )

                 ! tabulate distance and Ylm's w.r.t. other atoms
                 call tab_wave_ylm_v2 &
                       ( n_compute_atoms, atom_index, &
                       trigonom_tab(1,1,i_point), l_shell_max, &
                       l_ylm_max, &
                       ylm_tab(1,1,i_point) )


               ! tabulate total wave function value for each basis function in all cases -
               ! but only now are we sure that we have ylm_tab ...
               call evaluate_waves_v2 &
                   (l_ylm_max, ylm_tab(1,1,i_point), &
                    dist_tab(1,i_point), &
                   index_lm, n_compute, &
                   i_basis, radial_wave(1,i_point), &
                   wave(1,i_point), n_compute_atoms, &
                   atom_index_inv, &
                   i_basis_fns_inv &
                   )
!        end loop over z component
         enddo

!     end if (n_compute.gt.0)
      end if

      do i_cube = 1, n_cube, 1
        if (n_compute.gt.0) then

            do i_spin = 1, n_spin, 1
              if (max_occ_number(i_spin).gt.0) then
                if (cube_type(i_cube).eq.'total_density' .OR. &
                        cube_type(i_cube).eq.'spin_density') then
                    call evaluate_KS_density_v2 &
                        (n_points, wave(1,1), n_compute, &
                        i_basis, KS_vec_times_occ_sqrt(1,1,i_spin), &
                        KS_ev_compute(1,1,i_spin), &
                        max_occ_number(i_spin), &
                        occ_numbers_sqrt(1,i_spin), &
                        KS_orbital(1,1,i_spin), local_rho(1,i_spin) &
                    )
                else
             if (i_spin .eq. 1) then
                 call evaluate_KS_orbital_density &
                    (cube_index(i_cube),n_points, wave(1,1), n_compute, &
                    i_basis, KS_eigenvector(1,1,i_spin,1), &
                    KS_ev_compute(1,1,i_spin), &
                    max_occ_number(i_spin), &
                    local_orb_up(1,i_spin), &
                    local_rho(1,i_spin) &
                    )
             else
                 call evaluate_KS_orbital_density &
                    (cube_index(i_cube),n_points, wave(1,1), n_compute, &
                    i_basis, KS_eigenvector(1,1,i_spin,1), &
                    KS_ev_compute(1,1,i_spin), &
                    max_occ_number(i_spin), &
                    local_orb_down(1,i_spin), &
                    local_rho(1,i_spin) &
                    )
             endif
                endif
              end if
            end do

        else

            if ( (cube_type(i_cube).eq.'total_density') .OR. &
                  (cube_type(i_cube).eq.'spin_density') ) then
              local_rho = 0.d0
            else
              local_orb_up(1,1:n_points) = 0.0
             if (n_spin .gt. 1) local_orb_down(1,1:n_points) = 0.0
            endif

        end if

        ! plot data
        do coord_z = 1,cube_edge_steps(3,1),1
       if (n_spin .eq. 1) then
           if (cube_type(i_cube) .ne. 'eigenstate' ) then
                   write (10+i_cube,fmt='(1E13.5, $)') &
                       local_rho(coord_z,1) * inv_bohr_3
           else
                   write (10+i_cube,fmt='(1E13.5, $)') &
                       local_orb_up(1,coord_z) * sqrt_inv_bohr_3
           endif

       else

          if (cube_type(i_cube) .eq. 'total_density' ) then
                  write (10+i_cube,fmt='(1E13.5, $)') &
                      (local_rho(coord_z,1) + local_rho(coord_z,2)) &
                         * inv_bohr_3
        else if (cube_type(i_cube) .eq. 'spin_density' ) then
                  write (10+i_cube,fmt='(1E13.5, $)') &
                      (local_rho(coord_z,1) - local_rho(coord_z,2)) &
                         * inv_bohr_3
          else if (cube_type(i_cube) .eq. 'eigenstate') then
                  write (10+i_cube,fmt='(1E13.5, $)') &
                      (local_orb_up(1,coord_z)) &
                        * sqrt_inv_bohr_3
              write (100+i_cube,fmt='(1E13.5, $)') &
                      (local_orb_down(1,coord_z)) &
                        * sqrt_inv_bohr_3
        else if (cube_type(i_cube) .eq. 'eigenstate_density') then
                  write (10+i_cube,fmt='(1E13.5, $)') &
                      (local_rho(coord_z,1)) &
                        * inv_bohr_3
              write (100+i_cube,fmt='(1E13.5, $)') &
                      (local_rho(coord_z,2)) &
                        * inv_bohr_3

       end if

      end if

          if(mod(coord_z-1,6).eq.5) &
                  write (10+i_cube,*) ""
        if(mod(coord_z-1,6).eq.5 .AND. (cube_type(i_cube) .eq. &
                   'eigenstate' .or. cube_type(i_cube) .eq. 'eigenstate_density') &
                  .AND. n_spin .gt. 1) then
            write (100+i_cube,*) ""
        endif

       enddo

        enddo

        do i_cube = 1, n_cube, 1
                write (10+i_cube,*) ""
            if ((cube_type(i_cube) .eq.'eigenstate' .or. cube_type(i_cube) &
                        .eq.'eigenstate_density' ).AND. &
                        n_spin .gt. 1) write (100+i_cube,*) ""
        enddo




!     end loop over y component integration loop
      end do

!     end loop over x component
      end do

!     finally, deallocate stuff.
      if (allocated(ylm_tab)) then
        deallocate(ylm_tab)
      end if
      if (allocated(index_lm)) then
        deallocate(index_lm)
      end if
      if (allocated(dir_tab_global)) then
         deallocate(dir_tab_global)
      end if
      if (allocated(dist_tab_sq_global)) then
         deallocate(dist_tab_sq_global)
      end if
      if (allocated(local_rho)) then
         deallocate(local_rho)
      end if
      if (allocated(KS_orbital)) then
         deallocate(KS_orbital)
      end if
      if (allocated(local_orb_up)) then
         deallocate(local_orb_up)
      end if
      if (allocated(local_orb_down)) then
         deallocate(local_orb_down)
      end if

      do i_cube = 1, n_cube, 1
        close(10+i_cube)
      if ((cube_type(i_cube) .eq.'eigenstate' &
              .or. cube_type(i_cube) .eq.'eigenstate_density ').AND. &
               n_spin .gt. 1) close (100+i_cube)
      enddo

!     if master thread
      endif
      end subroutine output_cube_files
!******
!----------------------------------------------------------------------------
!****s* FHI-aims/write_cube_header
!  NAME
!    write_cube_header
!  SYNOPSIS
!
!      subroutine write_cube_header(descriptor,offset,num,scaling)
!
!!  PURPOSE
!!   Writes a header to cube format file.
!!
!!  USES
!!
!      use dimensions
!      use runtime_choices
!      use grids
!      use geometry
!      use species_data
!      use mpi_utilities
!      use constants
!      implicit none
!
!!  ARGUMENTS
!
!      integer  descriptor
!      integer  num(3)
!      real*8   offset(3)
!      real*8   scaling(3)
!
!
!!  INPUTS
!!   o descriptor -- file number where the result are written
!!   o num -- number of grid points
!!   o offset -- offset of the coordinates
!!   o scaling -- what is a scale between grid poitns and coordinates
!!
!!  OUTPUT
!!    none
!!  AUTHOR
!!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!!  SEE ALSO
!!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!!    Computer Physics Communications (2008), submitted.
!!  COPYRIGHT
!!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!!   the terms and conditions of the respective license agreement."
!!  HISTORY
!!    Release version, FHI-aims (2008).
!!  SOURCE
!!
!
!
!
!
!
!
!
!
!!     local
!      integer i_coord
!      integer i_atom
!      real*8   local_offset(3)
!
!      write (descriptor,*) "CUBE FILE written by FHI-AIMS"
!      write (descriptor,*) "*****************************"
!
!      local_offset = -offset
!
!      write (descriptor,fmt='(1X,I4,3F12.6)') &
!       n_atoms, (local_offset(i_coord), i_coord=1,3,1)
!
!      write (descriptor,fmt='(1X,I4,3F12.6)') &
!       num(1), scaling(1), 0.d0, 0.d0
!
!      write (descriptor,fmt='(1X,I4,3F12.6)') &
!       num(2), 0.d0, scaling(2), 0.d0
!
!      write (descriptor,fmt='(1X,I4,3F12.6)') &
!       num(3), 0.d0, 0.d0, scaling(3)
!
!      do i_atom = 1,n_atoms,1
!         write (descriptor,fmt='(1X,I4,4F12.6)') &
!         nint(species_z(species(i_atom))), 0.0d0, &
!         (coords(i_coord,i_atom),i_coord=1,3,1)
!      enddo
!
!      end subroutine write_cube_header
!!!!!!!!!!!!******
