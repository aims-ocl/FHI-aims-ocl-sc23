      module physics_ext

      real*8, dimension(:,:,:,:), allocatable :: &
           local_potential_parts_ext
      real*8, dimension(:,:,:,:), allocatable :: xc_gradient_deriv_ext

      real*8, dimension(:,:,:), allocatable :: free_hartree_superpos_ext
      real*8, dimension(:,:,:), allocatable :: free_rho_superpos_ext
      real*8, dimension(:,:,:,:), allocatable  :: rho_ext
      real*8, dimension(:,:,:,:,:), allocatable :: rho_gradient_ext

      contains

            subroutine allocate_physics_ext()

      use dimensions
      use runtime_choices

      implicit none

!     allocate quantities for integrals in s-c loop

      if (.not.allocated(rho_ext)) then
         allocate (rho_ext(n_spin, n_max_angular, &
              n_max_radial, n_atoms))
      end if
      if (.not.allocated(local_potential_parts_ext)) then
         allocate ( local_potential_parts_ext &
              (n_spin, n_max_angular, n_max_radial, n_atoms) )
      end if
      if (.not.allocated(free_hartree_superpos_ext)) then
         allocate ( free_hartree_superpos_ext(n_max_angular, &
              n_max_radial, n_atoms) )
      end if
      if (.not.allocated(free_rho_superpos_ext)) then
         allocate ( &
              free_rho_superpos_ext(n_max_angular, &
              n_max_radial, n_atoms) )
      end if

      if (use_density_gradient) then
         if (.not.allocated(rho_gradient_ext)) then
            allocate ( rho_gradient_ext &
                 (3, n_spin, n_max_angular, n_max_radial, n_atoms) )
         end if
      else
!     dummy allocation; this array should never be used
         allocate ( rho_gradient_ext &
              (1, 1, 1, 1, 1) )
      end if

      if (use_gga) then
         if (.not.allocated(xc_gradient_deriv_ext)) then
            allocate ( xc_gradient_deriv_ext &
                 (n_spin, n_max_angular, n_max_radial, n_atoms) )
         end if
      else
!     dummy allocation; this array should never be used
         allocate ( xc_gradient_deriv_ext (1,1,1,1) )
      end if

      end subroutine allocate_physics_ext

      subroutine deallocate_physics_ext()

      implicit none


      if (allocated(local_potential_parts_ext)) then
         deallocate( local_potential_parts_ext )
      end if
      if (allocated(xc_gradient_deriv_ext)) then
         deallocate( xc_gradient_deriv_ext )
      end if

      if (allocated(rho_ext)) then
         deallocate(rho_ext)
      end if

      if (allocated(rho_gradient_ext)) then
         deallocate ( rho_gradient_ext )
      end if

      end subroutine deallocate_physics_ext

      end module physics_ext
