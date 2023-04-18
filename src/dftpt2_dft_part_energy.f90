!****s* FHI-aims/dftpt2_dft_part_energy
!  NAME
!   dftpt2_dft_part_energy
!  SYNOPSIS

  subroutine dftpt2_dft_part_energy &
       ( partition_tab, &
         rho,rho_gradient, &
         kinetic_density, &
         xc_energy &
         )

!  PURPOSE
!   This subroutine is intended to calculate the exchange and correlation energy
!   separately of the underlying DFT calculation.
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
      real*8, dimension(n_full_points) :: partition_tab
      real*8, dimension(n_spin,n_full_points) :: rho
      real*8, dimension(3,n_spin,n_full_points) :: rho_gradient
      real*8, dimension(n_spin,n_full_points) :: kinetic_density

      real*8 :: xc_energy
      real*8, dimension(n_spin) :: x_energy
      real*8 :: c_energy
      real*8, dimension(n_spin) :: x_energy_lda
      real*8 :: c_energy_lda

!  INPUTS
!  o partition_tab -- the partition function for real space integration     
!  o rho -- the charge density
!  o rho_gradient -- the gradient of the charge density
!  o kinetic_density -- kinetic-energy density
!  OUTPUTS
!  o xc_energy -- the exchange-correlation energy    
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
  real*8, dimension(n_spin) :: en_density_x
  real*8 :: en_density_c
  real*8 :: local_xc_density_deriv(n_spin)
  real*8 :: local_xc_gradient_deriv(3,n_spin)
  real*8 :: local_xc_tau_deriv(n_spin)

  real*8 :: en_density_ldaxc
  real*8, dimension(n_spin) :: en_density_ldax
  real*8 :: en_density_ldac
  real*8 :: xc_energy_mpi
  real*8, dimension(n_spin) :: pot_xc

  integer :: flag_xc_old

  integer :: mpierr

  logical :: xc_undefined

  !  counters

  integer i_my_batch
  integer i_index
  integer i_spin


  integer i_full_points
  real*8, dimension(3) :: coord_current

  !  begin work

  !      if(myid.eq.0) then
  !        write(use_unit,*)
  !        write(use_unit,'(2X,A)')
  !     +    "Integrating the exchange and correlation energy ... "
  !      endif

  !     initialize

  ! switch the flag_xc for doubly-hybrid density functional
  ! evaluation
  flag_xc_old = flag_xc
  flag_xc     = flag_dftpt2_dft_part

  xc_energy = 0.d0

  i_full_points = 0
  !     perform partitioned integration, atom by atom, and point by point
  !     This will be the outermost loop, to save evaluations of the potential.
  !     and the Y_lm functions

  do  i_my_batch = 1, n_my_batches, 1

        do i_index = 1, batches(i_my_batch)%size, 1
           
           coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:) !SAG
           if(n_periodic > 0)then
              call map_to_center_cell(coord_current(1:3) )
           end if

           i_full_points = i_full_points + 1
           if(partition_tab(i_full_points) .gt. 0.d0) then

              call evaluate_xc &
                   ( rho(1:n_spin,i_full_points), &
                   rho_gradient(1:3,1:n_spin, i_full_points), &
                   kinetic_density(1:n_spin,i_full_points), &
                   en_density_xc, &
                   en_density_x, en_density_c,  &
                   local_xc_density_deriv(1:n_spin), &
                   local_xc_gradient_deriv(1:3,1:n_spin), &
                   local_xc_tau_deriv(1:n_spin), &
                   .false., coord_current )  !SAG
              do i_spin = 1, n_spin, 1
                 xc_energy = xc_energy + &
                      rho(i_spin,i_full_points)* &
                      en_density_xc * &
                      partition_tab(i_full_points)
              enddo
           endif
        enddo
        !       end of mpi distribution
     !endif
     !test end
     !     end integration loop over batches
enddo


  if(use_mpi) then
     call sync_real_number(xc_energy)
  endif

  flag_xc = flag_xc_old

  return
end subroutine dftpt2_dft_part_energy

