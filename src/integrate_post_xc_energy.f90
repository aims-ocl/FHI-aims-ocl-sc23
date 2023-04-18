!---------------------------------------------------------------------------
!****s* FHI-aims/integrate_post_xc_energy
!  NAME
!   integrate_post_xc_energy
!  SYNOPSIS

  subroutine integrate_post_xc_energy &
       ( partition_tab, &
         rho,rho_gradient, &
         kinetic_density, &
         en_post_xc, &
         x_energy_post, &
         c_energy_post)

!  PURPOSE
! Subroutine designed to integrate xc energy density for meta-gga post processing
! kinetic_density is used inside but is not modified
!
!  USES

    use dimensions
    use runtime_choices
    use grids
    use xc
    use mpi_tasks
    use mpi_utilities
    use synchronize_mpi
    implicit none

!  ARGUMENTS
!     input
      real*8, dimension(n_full_points)              :: partition_tab
      real*8, dimension(n_spin,n_full_points)       :: rho
      real*8, dimension(3,n_spin,n_full_points)     :: rho_gradient
      real*8, dimension(n_spin, n_full_points)      :: kinetic_density
!     output
      real*8 :: en_post_xc
      real*8 :: x_energy_post
      real*8 :: c_energy_post

!  INPUTS
!  o partition_tab -- the partition function for real space integration     
!  o rho -- the charge density
!  o rho_gradient -- the gradient of the charge density
!  o kinetic_density -- the kinetic density

!  OUTPUTS
!  o xc_energy -- the exchange-correlation energy    
!  o x_energy -- the exchange energy    
!  o c_energy -- the correlation energy
!
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

  !  local variables


  !     and condensed version of partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)

  real*8 :: en_density_xc
  real*8 :: en_density_x
  real*8 :: en_density_c
!  real*8 :: local_xc_density_deriv(n_spin)
!  real*8 :: local_xc_gradient_deriv(3,n_spin)

  real*8 :: xc_energy_mpi
  real*8, dimension(n_spin) :: pot_xc

  integer :: mpierr

  logical :: xc_undefined

  !  counters

  integer i_my_batch
  integer i_index
  integer i_spin


  integer i_full_points
  real*8 :: coord_current(3)

  !  begin work
  !     initialize

  en_post_xc = 0.d0
  x_energy_post = 0.d0
  c_energy_post = 0.d0

  i_full_points = 0
  !     perform partitioned integration, atom by atom, and point by point
  !     This will be the outermost loop, to save evaluations of the potential.
  !     and the Y_lm functions

  do  i_my_batch = 1, n_my_batches, 1



        do i_index = 1, batches(i_my_batch)%size, 1



           i_full_points = i_full_points + 1


           if(partition_tab(i_full_points) .gt. 0.d0) then


              coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)


              call evaluate_post_xc &
                   ( rho(1:n_spin,i_full_points), &
                   rho_gradient(1:3,1:n_spin, i_full_points), &
                   en_density_xc,  en_density_x, en_density_c, &
!                   local_xc_density_deriv(1:n_spin), &
!                   local_xc_gradient_deriv(1:3,1:n_spin), &
                   coord_current(:), &
                   kinetic_density(1:n_spin, i_full_points) &
                   )

              do i_spin = 1, n_spin, 1
                 en_post_xc = en_post_xc + &
                      rho(i_spin,i_full_points)* &
                      en_density_xc * &
                      partition_tab(i_full_points)

                 x_energy_post = x_energy_post + &
                      rho(i_spin,i_full_points)* &
                      en_density_x * &
                      partition_tab(i_full_points)

                 c_energy_post = c_energy_post + &
                      rho(i_spin,i_full_points)* &
                      en_density_c * &
                      partition_tab(i_full_points)
!

              enddo !i_spin

           endif! partition_tab

        enddo !i_index
        !       end of mpi distribution

     !     end integration loop over batches
  enddo !i_my_batch

  if(use_mpi) then
     call sync_real_number(en_post_xc)
     call sync_real_number(x_energy_post)
     call sync_real_number(c_energy_post)
  endif

  return
end subroutine integrate_post_xc_energy
