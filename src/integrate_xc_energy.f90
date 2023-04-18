!****s* FHI-aims/integrate_xc_energy
!  NAME
!   integrate_xc_energy
!  SYNOPSIS

  subroutine integrate_xc_energy &
       ( partition_tab, &
         rho,rho_gradient, &
         kinetic_density, &
         xc_energy, &
         x_energy, &
         c_energy, &
         x_energy_lda, &
         c_energy_lda &
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
    use pseudodata
    use energy_density

    implicit none

!  ARGUMENTS
      real*8, dimension(n_full_points)          :: partition_tab
      real*8, dimension(n_spin,n_full_points)   :: rho
      real*8, dimension(3,n_spin,n_full_points) :: rho_gradient
      real*8, dimension(n_spin,n_full_points)   :: kinetic_density

      real*8 :: xc_energy
      real*8, dimension(n_spin) :: x_energy
      real*8 :: c_energy
      real*8, dimension(n_spin) :: x_energy_lda
      real*8 :: c_energy_lda

!  INPUTS
!  o partition_tab -- the partition function for real space integration     
!  o rho -- the charge density
!  o rho_gradient -- the gradient of the charge density
!  o kinetic_density -- the gradient of the orbitals

!  OUTPUTS
!  o xc_energy -- the exchange-correlation energy    
!  o x_energy -- the exchange energy    
!  o c_energy -- the correlation energy
!  o x_energy_lda -- the LDA exchange energy    
!  o c_energy_lda -- the LDA correlation energy
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

  real*8 :: en_density_xc
  real*8, dimension(n_spin) :: en_density_x
  real*8 :: en_density_c
  real*8 :: local_xc_density_deriv(n_spin)
  real*8 :: local_xc_gradient_deriv(3,n_spin)
  real*8 :: local_xc_tau_deriv(n_spin)

  real*8 :: en_density_ldaxc
  real*8, dimension(n_spin) :: en_density_ldax
  real*8 :: en_density_ldac
  real*8, dimension(n_spin) :: pot_xc


  logical :: xc_undefined

  !  counters

  integer i_my_batch
  integer i_index
  integer i_spin


  integer i_full_points
  real*8, dimension(3) :: coord_current

  real*8, dimension(:,:), allocatable  :: rho_inc_partialcore
  real*8, dimension(:,:,:), allocatable  :: rho_gradient_inc_partialcore

  !  begin work

  !  initialize
  xc_energy = 0.d0
  x_energy = 0.d0
  c_energy = 0.d0
  x_energy_lda =0.d0
  c_energy_lda =0.d0

  if (use_embedding_pp.and.use_nonlinear_core) then
    allocate(rho_inc_partialcore(n_spin,n_full_points))
    allocate(rho_gradient_inc_partialcore(3,n_spin,n_full_points))
    do i_full_points = 1, n_full_points
      do i_spin = 1,n_spin
         rho_inc_partialcore(i_spin,i_full_points) = &
            rho(i_spin,i_full_points) + partial_core_rho(i_full_points)
         rho_gradient_inc_partialcore(1:3,i_spin,i_full_points) = &
            rho_gradient(1:3,i_spin,i_full_points)
         if (use_density_gradient) then
            rho_gradient_inc_partialcore(1:3,i_spin,i_full_points) = &
               rho_gradient_inc_partialcore(1:3,i_spin,i_full_points) + &
               partial_core_rho_grad(1:3,i_full_points)
         endif
      enddo
    enddo
  endif

  !CC: Allocate array for XC-CM-Energy-Density
  if (flag_energy_density) then
    if (.not. allocated(ed_xc_energy_density) ) then
      allocate(ed_xc_energy_density(n_full_points))
    end if
    ed_xc_energy_density(:) = 0.0d0
  end if

  !     perform partitioned integration, atom by atom, and point by point
  i_full_points = 0
  do  i_my_batch = 1, n_my_batches, 1

        do i_index = 1, batches(i_my_batch)%size, 1
           
           coord_current(:) = &
              batches(i_my_batch) % points(i_index) % coords(:) !SAG
           if(n_periodic > 0)then
              call map_to_center_cell(coord_current(1:3) )
           end if

           i_full_points = i_full_points + 1
           if(partition_tab(i_full_points) .gt. 0.d0) then

              if(use_embedding_pp.and.use_nonlinear_core) then

                 call evaluate_xc &
                      ( rho_inc_partialcore(1:n_spin,i_full_points), &
                      rho_gradient_inc_partialcore(1:3,1:n_spin, i_full_points), &
                      kinetic_density(1:n_spin, i_full_points), &
                      en_density_xc, &
                      en_density_x, en_density_c,  &
                      local_xc_density_deriv(1), &
                      local_xc_gradient_deriv(1,1), &
                      local_xc_tau_deriv(1), &
                      .false., coord_current )

              else

                 call evaluate_xc &
                      ( rho(1:n_spin,i_full_points), &
                      rho_gradient(1:3,1:n_spin, i_full_points), &
                      kinetic_density(1:n_spin, i_full_points), &
                      en_density_xc, &
                      en_density_x, en_density_c,  &
                      local_xc_density_deriv(1), &
                      local_xc_gradient_deriv(1,1), &
                      local_xc_tau_deriv(1), &
                      .false., coord_current )

              endif
              
              xc_undefined = .true.

              do i_spin = 1, n_spin, 1
                 ! This is true if both densities are below zero
                 xc_undefined = xc_undefined .and. &
                      (rho(i_spin,i_full_points).le.0.d0)
              enddo

              do i_spin = 1, n_spin, 1
                 ! This is true if one of the densities is less than zero
                 xc_undefined = xc_undefined .or. &
                      (rho(i_spin,i_full_points).lt.0.d0)
              enddo

              if (xc_undefined) then

                 en_density_ldax = 0.d0
                 en_density_ldac = 0.d0
                 en_density_ldaxc = 0.d0

              else

                 call xcpot_pw91_lda &
                      ( rho(1:n_spin,i_full_points), &
                      pot_xc, en_density_ldaxc, &
                      en_density_ldax, en_density_ldac &
                      )

              endif

              do i_spin = 1, n_spin, 1
                 xc_energy = xc_energy + &
                      rho(i_spin,i_full_points)* &
                      en_density_xc * &
                      partition_tab(i_full_points)

                 ! AJL, Apr 2017. Unpleasant hack to get the correct output
                 ! for meta-GGAs. Taken from the post-processing implementation
                 if (use_meta_gga) then
                   x_energy(1) = x_energy(1)+ &
                        rho(i_spin,i_full_points)* &
                        en_density_x(1) * &
                        partition_tab(i_full_points)
                 else
                   x_energy(i_spin) = x_energy(i_spin)+ &
                        rho(i_spin,i_full_points)* &
                        en_density_x(i_spin) * &
                        partition_tab(i_full_points)
                 endif

                 c_energy = c_energy + &
                      rho(i_spin,i_full_points)* &
                      en_density_c * &
                      partition_tab(i_full_points)

                 x_energy_lda(i_spin) = x_energy_lda(i_spin) + &
                      rho(i_spin,i_full_points)* &
                      en_density_ldax(i_spin) * &
                      partition_tab(i_full_points)

                 c_energy_lda = c_energy_lda + &
                      rho(i_spin,i_full_points)* &
                      en_density_ldac * &
                      partition_tab(i_full_points)
              enddo

              !CC: Save local xc_energy_density in ed_array
              if (flag_energy_density) then
                 do i_spin = 1, n_spin, 1
                   ed_xc_energy_density(i_full_points) = &
                        ed_xc_energy_density(i_full_points) &
                      + rho(i_spin,i_full_points) * en_density_xc


                 enddo
              end if

           endif

        enddo !       end of mpi distribution
  enddo !     end integration loop over batches

  if(use_mpi) then
     call sync_real_number(xc_energy)
     call sync_vector(x_energy,n_spin)
     call sync_real_number(c_energy)
     call sync_vector(x_energy_lda,n_spin)
     call sync_real_number(c_energy_lda)
  endif

  if (allocated(rho_inc_partialcore)) deallocate(rho_inc_partialcore)
  if (allocated(rho_gradient_inc_partialcore)) &
     deallocate(rho_gradient_inc_partialcore)

  return
end subroutine integrate_xc_energy

!----------------------------------------------------------------------
!************
!  NAME
!   integrate_crpa
!  SYNOPSIS

  subroutine integrate_crpa &
       ( partition_tab, &
         rho,rho_gradient, &
         c_energy_ldarpa &
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
      real*8, dimension(n_full_points) :: &
          partition_tab

      real*8, dimension(n_spin,n_full_points) :: &
          rho

      real*8, dimension(3,n_spin,n_full_points) :: &
           rho_gradient

      real*8 :: c_energy_ldarpa

!  INPUTS
!  o partition_tab -- the partition function for real space integration     
!  o rho -- the charge density
!  o rho_gradient -- the gradient of the charge density

!  OUTPUTS
!  o c_energy_ldarpa -- the RPA contribution to the LDA correlation energy
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

  real*8 :: en_density_ldac_rpa
  real*8 :: dup, ddn, vcup, vcdn

  !  counters

  integer i_my_batch
  integer i_index
  integer i_spin


  integer i_full_points

  !  begin work

  !     initialize
  c_energy_ldarpa =0.d0

  !     perform partitioned integration, atom by atom, and point by point
  i_full_points = 0
  do  i_my_batch = 1, n_my_batches, 1

        do i_index = 1, batches(i_my_batch)%size, 1

           i_full_points = i_full_points + 1
           if(partition_tab(i_full_points) .gt. 0.d0) then

              if(n_spin.eq.1) then
                 dup = 0.5d0*rho(1,i_full_points)
                 ddn = dup
              elseif(n_spin.eq.2) then
                 dup = rho(1,i_full_points)
                 ddn = rho(2,i_full_points)
              endif

              call crpapw( dup,ddn,vcup,vcdn, &
                   en_density_ldac_rpa &
                   )

              do i_spin = 1, n_spin, 1

                      c_energy_ldarpa = c_energy_ldarpa + &
                      rho(i_spin,i_full_points)* &
                      en_density_ldac_rpa * &
                      partition_tab(i_full_points)
              enddo

           endif

        enddo !       end of mpi distribution     
  enddo !     end integration loop over batches


  if(use_mpi) then
     call sync_real_number(c_energy_ldarpa)
  endif

  return
end subroutine integrate_crpa

!----------------------------------------------------------------------
!************
