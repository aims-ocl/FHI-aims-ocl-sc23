!****s* FHI-aims/evaluate_dielectric_constant
!  NAME
!    evaluate_dielectric_constant
!  SYNOPSIS

subroutine evaluate_dielectric_constant( & 
           occ_numbers, Omega_MO, first_order_U_complex, & 
           dielectric_constant)


!  PURPOSE
!    calculate dielectric_constants.

!  dielectric_constants(alpha,beta)^infinite 
!   = delta(alpha,beta) +  4*pi/V * d P_alpha/ dE_beta 
!
!                           4pi                   <i|grad_alpha|j>    <j|H^1_beta|i>
!   = delta(alpha,beta) +  ----- * 4*(sum_(i,j,k)---------------- *------------------   
!                            V                      Eii-Ejj            Ejj-Eii
!  USES

  use constants, only: pi
  use dimensions
  use pbc_lists , only : k_weights
  use mpi_tasks
  use synchronize_mpi
  use localorb_io
  use geometry, only: cell_volume

!  ARGUMENTS

  implicit none
 
  real*8, dimension(n_states, n_spin,n_k_points), intent(IN) :: occ_numbers 
  complex*16, dimension(n_states, n_states, n_k_points_task,3), intent(IN) :: Omega_MO
  complex*16, dimension(n_states, n_states, n_k_points_task,3), intent(IN) :: first_order_U_complex
 
  real*8, dimension(3,3), intent(out) :: dielectric_constant

!  INPUTS
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  dielectric_constant

  integer :: i_state, j_state,i_coord, j_coord, i_spin
  integer :: n_occ_states(n_spin)
  integer :: i_k_task, i_k_point
  integer :: i_basis
  character*1000 :: info_str

  complex*16, dimension(3,3) :: dP_dE
 

  dP_dE(1:3,1:3) = (0.0d0,0.0d0)


   i_k_task = 0
   do i_k_point = 1,n_k_points, 1
   if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
        i_k_task = i_k_task + 1

       do i_spin = 1, n_spin, 1
            n_occ_states(i_spin) = 0
            do i_state = n_states, 1, -1
            !if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.0.d0) then
            ! The following line implies the system does not have any fractional occupation number
            if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.1.e-6) then
               n_occ_states(i_spin) = i_state
            exit
            endif
            enddo
       enddo

       do i_state = 1,n_states
          if(dabs(occ_numbers(i_state,1,i_k_point)).gt.1.e-6.and.dabs(occ_numbers(i_state,1,i_k_point)-2.0d0).gt.1.e-6) then
             write(use_unit,*) 'Warning for fractional occ_numbers:',i_k_point,i_state, occ_numbers(i_state,1,i_k_point)
          endif
       enddo

      do i_coord = 1, 3
      do j_coord = 1, 3

         do i_state = 1, n_occ_states(1)
         !do j_state = 1, n_states
         do j_state = n_occ_states(1)+1, n_states
 
!-------This is for non_frac version:    2 for occ_number, 2 for exchange.
!        dP_dE(i_coord,j_coord) = dP_dE(i_coord,j_coord) & 
!          - k_weights(i_k_point) * 4.0d0 * Omega_MO(i_state,j_state,  i_k_task, i_coord) &  
!          * first_order_U_complex(j_state,i_state, i_k_task, j_coord)

!-------This is the more general version, if there is fractional occ_numbers, this version may
!-------have little difference with the above one, but within 1e-4. (shanghui note)  
        dP_dE(i_coord,j_coord) = dP_dE(i_coord,j_coord) &
          - k_weights(i_k_point) * 2.0d0 * occ_numbers(i_state,1,i_k_point) &
          * Omega_MO(i_state,j_state,  i_k_task, i_coord) &
          * first_order_U_complex(j_state,i_state, i_k_task, j_coord)
          !* (-dconjg(first_order_U_complex(i_state,j_state, i_k_task, j_coord)))


          ! if(i_coord.eq.1.and.j_coord.eq.2) then
          !    write(use_unit,*) 'i,j,k_weights:', i_state,j_state, k_weights(i_k_point)
          !    write(use_unit,*) 'Omega_MO(i,j):', Omega_MO(i_state,j_state,  i_k_task, i_coord)  
          !    write(use_unit,*) 'U1(j,i)      :', first_order_U_complex(j_state,i_state, i_k_task, j_coord) 
          !    write(use_unit,*) 'total:        ',  - k_weights(i_k_point) * 4.0d0 * Omega_MO(i_state,j_state,  i_k_task, i_coord) &
          !    * first_order_U_complex(j_state,i_state, i_k_task, j_coord)
          !    write(use_unit,*) 'total sum:    ', dP_dE(i_coord,j_coord)
          ! endif

         enddo  ! i_state
         enddo  ! j_state
 
      enddo     ! i_coord 
      enddo     ! j_coord 
         
   endif ! i_k_task   
   enddo ! i_k_point


  !-------shanghui begin parallel------
   call  sync_vector_complex(dP_dE, 9)
  !-------shanghui end parallel------

   
   write(info_str,'(A)') 'dP_dE (Bohr^3) at every cycle:--->'
   call localorb_info(info_str, use_unit,'(A)', OL_norm)
   write(info_str,*) dP_dE(1,1:3)
   call localorb_info(info_str, use_unit,'(A)', OL_norm)
   write(info_str,*) dP_dE(2,1:3)
   call localorb_info(info_str, use_unit,'(A)', OL_norm)
   write(info_str,*) dP_dE(3,1:3)
   call localorb_info(info_str, use_unit,'(A)', OL_norm)


   write (info_str,'(A)') 'DFPT polarizability (Bohr^3)        xx        yy        zz        xy        xz        yz'
   call localorb_info(info_str, use_unit,'(A)', OL_norm)
   write (info_str,'(2X,A,1X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3)') &
   '| Polarizability:--->          ', real(dP_dE(1,1)), real(dP_dE(2,2)), &
   & real(dP_dE(3,3)), real(dP_dE(1,2)), real(dP_dE(1,3)), real(dP_dE(2,3))
   call localorb_info(info_str, use_unit,'(A)', OL_norm)


   do i_coord = 1, 3 
   do j_coord = 1, 3

      if(i_coord .eq. j_coord ) then 
        dielectric_constant(i_coord,j_coord) = 1 + dble(4*pi*dP_dE(i_coord,j_coord)/(cell_volume))
      else 
        dielectric_constant(i_coord,j_coord) =  dble(4*pi*dP_dE(i_coord,j_coord)/(cell_volume))
      endif

   enddo     ! i_coord 
   enddo     ! j_coord 

   !write(info_str,'(A)') 'DFPT for dielectric_constant at every cycle:--->'
   !call localorb_info(info_str, use_unit,'(A)', OL_norm)
   !write(info_str,*) dielectric_constant(1,1:3)
   !call localorb_info(info_str, use_unit,'(A)', OL_norm)
   !write(info_str,*) dielectric_constant(2,1:3)
   !call localorb_info(info_str, use_unit,'(A)', OL_norm)
   !write(info_str,*) dielectric_constant(3,1:3)
   !call localorb_info(info_str, use_unit,'(A)', OL_norm)


end subroutine evaluate_dielectric_constant
!******
