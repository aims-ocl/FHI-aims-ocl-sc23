!****s* FHI-aims/evaluate_nonlocal_energy
!  NAME
!   evaluate_nonlocal_energy
!  SYNOPSIS

subroutine evaluate_nlcc_correction_energy(partition_tab, rho, rho_gradient, kinetic_density)

!  PURPOSE
!     Subroutine evaluate_nonlocal_engery
!     calculates energy term contributed by a nonlocal embedding potential.
!
!  USES

  use dimensions
  use pseudodata
  use runtime_choices, only: real_eigenvectors
  use mpi_tasks, only: aims_stop
  use synchronize_mpi
  use grids, only: batches

  implicit none

!  ARGUMENTS



!  INPUTS
!  o KS_coefficents -- coefficient matrix of KS_states
!  o occ_numbers -- occupation numbers of KS_states

!  
!  OUTPUT
!  o en_nonlocal -- value of energy contributed by nonlocal embedding potential
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

  real*8, dimension(n_full_points)             :: partition_tab
  real*8, dimension(n_spin, n_full_points)     :: rho
  real*8, dimension(3, n_spin, n_full_points)     :: rho_gradient
  real*8, dimension(n_spin, n_full_points)     :: kinetic_density

  real*8 :: rho_inc_partialcore, local_xc_derivs(n_spin)
  real*8 :: rho_gradient_inc_partialcore(3)
  real*8 :: en_density_xc, en_density_x(n_spin), en_density_c
  real*8 :: local_xc_gradient_deriv(3,n_spin)
  real*8 :: local_xc_tau_deriv(n_spin)

  !     counters
  integer :: i_full_points, i_index, i_my_batch, i_spin




   character(*), parameter :: func = 'evaluate_nlcc_correction_energy'

  en_xc_nlcc = 0.d0

  i_full_points = 0


  do i_my_batch = 1, n_my_batches, 1

        do i_index = 1, batches(i_my_batch)%size, 1

           i_full_points = i_full_points + 1
           if (partition_tab(i_full_points).gt.0.d0) then


              do i_spin = 1, n_spin, 1
                rho_inc_partialcore = 0.d0


                       
                 rho_inc_partialcore = &
                     rho(i_spin,i_full_points) + partial_core_rho(i_full_points) !0.5

                 if(use_density_gradient) then 
                      rho_gradient_inc_partialcore(1:3) = &
                      rho_gradient(1:3,i_spin,i_full_points) + & 
                      partial_core_rho_grad(1:3,i_full_points)
                 endif

                 call evaluate_xc  &
                      ( rho_inc_partialcore,   &
                        rho_gradient_inc_partialcore(1:3),  &
                        kinetic_density, &
                        en_density_xc,   &
                        en_density_x, en_density_c,  &
                        local_xc_derivs(i_spin),  &
                        local_xc_gradient_deriv(1:3,i_spin), &
                        local_xc_tau_deriv(i_spin), &
                        .false. &
                      )

                 en_xc_nlcc = en_xc_nlcc + partition_tab(i_full_points)* &
                       en_density_xc* partial_core_rho(i_full_points)


          enddo
        endif  
     enddo
  enddo


  call sync_real_number(en_xc_nlcc)


end subroutine evaluate_nlcc_correction_energy

!---------------------------------------------------------------------
!******	
