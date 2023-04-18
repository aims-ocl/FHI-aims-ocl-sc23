!****s* FHI-aims/evaluate_first_order_U_dielectric
!  NAME
!    evaluate_first_order_U_dielectric
!  SYNOPSIS

subroutine evaluate_first_order_U_dielectric( & 
           Omega_MO, &
           first_order_H_complex,   &
           first_order_U_complex)
               

!  PURPOSE
!    calculate the first-order U with its occupied-unoccupied part.

!  Uia(1)=(-C(0)H(1)C(0))ia/(Eii(1)-Eaa(1))
!  shanghui,2012.04.27

! note: KS_eigenvalue now is only for nspin=1 and ikpoint=1 
! shanghui, 2012.05.30

! shanghui now extended it to pbc case:evaluate_first_order_U_p0.f90
!  USES

  use dimensions
  use runtime_choices, only : real_eigenvectors
  use mpi_tasks
  use physics, only : occ_numbers, KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, chemical_potential 

  use force_occupation

!  ARGUMENTS

  implicit none
  

  complex*16, dimension( n_states, n_states, n_k_points_task), intent(IN) :: Omega_MO
  complex*16, dimension( n_basis, n_basis, n_k_points_task), intent(IN) :: first_order_H_complex
  complex*16, dimension( n_states, n_states, n_k_points_task), intent(OUT) :: first_order_U_complex

!  INPUTS
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  first_order_U_complex 


  integer :: i_state, j_state, i_spin
  integer :: i_k_task, i_k_point
  integer :: max_occ_number(n_spin)
  integer :: n_occ_states

  complex*16, allocatable     ::  temp_first_order(:,:)
  complex*16, allocatable     ::  temp_1(:,:),temp_2(:,:)
  complex*16, allocatable     ::  temp_H(:,:)
  complex*16, allocatable     ::  temp_eigenvector(:,:)

  real*8 :: theta_focc, theta_fvirt, theta_virtocc, lim ! integral of the smearing function
  logical :: warning=.false.
  complex*16     ::  zero, one

  zero = (0.d0,0.d0)
  one = (1.d0,0.d0)


  allocate(temp_first_order(n_basis,n_basis))
  allocate(temp_1(n_states,n_basis))
  allocate(temp_2(n_states,n_states))
  allocate(temp_H(n_states,n_states))
  allocate(temp_eigenvector(n_basis,n_basis))

  i_k_task = 0
  do i_k_point = 1,n_k_points, 1
  if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
        i_k_task = i_k_task + 1

   !-------initialize------------
   first_order_U_complex(:,:,i_k_task) = (0.0d0,0.0d0)

       do i_spin = 1, n_spin, 1
          max_occ_number(i_spin) = 0
          do i_state = n_states, 1, -1
             if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.1.e-6) then
             max_occ_number(i_spin) = i_state
             !print*, 'MAX_OCC', i_state, i_k_point, occ_numbers(i_state,i_spin,i_k_point),KS_eigenvalue(i_state,1,i_k_point)

             !theta_focc = 0.5*(1.+tanh((chemical_potential-KS_eigenvalue(i_state,1,i_k_point))/(2.*DFPT_width)) )
             !print*, 'THETA', theta_focc,occ_numbers(i_state,i_spin,i_k_point)
             exit
             endif
          enddo
       enddo

       n_occ_states=max_occ_number(1) 
  ! shanghui need to change this in the future for n_spin

  ! Nath: check if there are some fractional occupation numbers and print a warning
       if (DFPT_width .eq. 0.0) then
         do i_spin = 1, n_spin, 1
            do i_state = 1, n_occ_states, 1
               if ( (2.0-dabs(occ_numbers(i_state,i_spin,i_k_point))) .gt. 1.e-3) then
                 warning=.true.
               exit
               endif
            enddo
         enddo
       endif

   if(real_eigenvectors) then
          temp_eigenvector(:,:)=          &
          cmplx(KS_eigenvector(:,:,1,i_k_task))
   else
          temp_eigenvector(:,:)=          &
          !dconjg(KS_eigenvector_complex(:,:,1,i_k_task))
          KS_eigenvector_complex(:,:,1,i_k_task)
   endif

     

 !---------(1) Ct*H(1)*C ----------------------------------------
    temp_first_order(1:n_basis,1:n_basis)=          &
    first_order_H_complex(1:n_basis,1:n_basis, i_k_task)

    temp_1(1:n_states,1:n_basis)=(0.0d0,0.0d0)
 
    CALL zgemm("C","N",n_states,n_basis,n_basis,&
               one, temp_eigenvector, n_basis,temp_first_order, n_basis,&
               zero,temp_1,n_states)

    temp_H(1:n_states,1:n_states)=(0.0d0,0.0d0)

    CALL zgemm("N","N",n_states,n_states,n_basis,&
                one,temp_1, n_states,temp_eigenvector, n_basis,&
                zero,temp_H, n_states)


    ! Nath to Hui: I don't think this comment is needed anymore, right ? I leave it here as a memory :D
    !if(i_k_point.eq.1) then
    !write(use_unit,*) ' ' ! only when I add this output the complie -O3 at thnec could get right. 
    !endif 

   !Uia : we add H^(1) (electric field part) here : 
   do i_state=1,n_occ_states
      !do j_state=1, n_states
      do j_state=n_occ_states+1, n_states
      !if(i_state.ne.j_state) then
         !if (abs(KS_eigenvalue(i_state,1,i_k_point)-KS_eigenvalue(j_state,1,i_k_point)).lt. 1.e-2) then 
         ! print*, "DIVERGENCE", KS_eigenvalue(i_state,1,i_k_point)-KS_eigenvalue(j_state,1,i_k_point), i_state, j_state, i_k_point
         !end if

!----------------------we do not need this term, we calculate here just for debug--------

     ! Nath: if no smearing is required, use usual expression (finite gap, no metal)
     if (DFPT_width .eq. 0.0) then

       first_order_U_complex( i_state,j_state, i_k_task)=         &
       -Omega_MO(i_state,j_state,i_k_task)/(KS_eigenvalue(i_state,1,i_k_point)-KS_eigenvalue(j_state,1,i_k_point)) & 
       -temp_H(i_state,j_state)/(KS_eigenvalue(i_state,1,i_k_point)-KS_eigenvalue(j_state,1,i_k_point))

     else ! Extension to metals (fractional occupation numbers)
     ! See 'Lattice dynamics of metals from density-functional perturbation theory', Stefano de Gironcoli,
     ! Phys. Rev.B51, 6773(R), (doi.org/10.1103/PhysRevB.51.6773)

       theta_focc = 0.5*(1.+tanh((chemical_potential-KS_eigenvalue(i_state,1,i_k_point))/(2.*DFPT_width)) )
       theta_fvirt = 0.5*(1.+tanh((chemical_potential-KS_eigenvalue(j_state,1,i_k_point))/(2.*DFPT_width)) )
       theta_virtocc = 0.5*(1.+tanh((KS_eigenvalue(j_state,1,i_k_point)-KS_eigenvalue(i_state,1,i_k_point))/(2.*DFPT_width)) )

       if (abs(KS_eigenvalue(j_state,1,i_k_point)-KS_eigenvalue(i_state,1,i_k_point)) .gt. 5.e-2 ) then
 
         !print*, "DIFF", 2*theta_virtocc*(theta_focc-theta_fvirt),i_state, j_state, i_k_point
 
         first_order_U_complex( i_state,j_state, i_k_task)=         &
         -Omega_MO(i_state,j_state,i_k_task)&
         *theta_virtocc*(theta_focc-theta_fvirt)/(KS_eigenvalue(i_state,1,i_k_point)-KS_eigenvalue(j_state,1,i_k_point)) & 
         -temp_H(i_state,j_state)&
         *theta_virtocc*(theta_focc-theta_fvirt)/(KS_eigenvalue(i_state,1,i_k_point)-KS_eigenvalue(j_state,1,i_k_point))
 
       else ! Switch to the limit when the difference between 2 eigenvalues is too small
         lim = -0.5*DFPT_width*1./(1+cosh((chemical_potential-KS_eigenvalue(i_state,1,i_k_point))/DFPT_width))
        
         first_order_U_complex( i_state,j_state, i_k_task)=         &
         -Omega_MO(i_state,j_state,i_k_task)*lim&
         -temp_H(i_state,j_state)*lim
       endif

     endif ! Extension to metals

!----------------in fact, we only neet U(unocc,occ)=U(a,i) like this:----------

     ! Nath: if no smearing is required, use usual expression (finite gap, no metal)
     if (DFPT_width .eq. 0.0) then

       first_order_U_complex( j_state,i_state,i_k_task)= & 
       -Omega_MO(j_state,i_state,i_k_task)/(KS_eigenvalue(j_state,1,i_k_point)-KS_eigenvalue(i_state,1,i_k_point)) & 
       -temp_H(j_state,i_state)/(KS_eigenvalue(j_state,1,i_k_point)-KS_eigenvalue(i_state,1,i_k_point))

     else ! Extension to metals (fractional occupation numbers)

       if (abs(KS_eigenvalue(j_state,1,i_k_point)-KS_eigenvalue(i_state,1,i_k_point)) .gt. 5.e-2 ) then

        first_order_U_complex( j_state,i_state, i_k_task)=         &
        -Omega_MO(j_state,i_state,i_k_task)&
        !*2*theta_virtocc*(theta_focc-theta_fvirt)/(KS_eigenvalue(j_state,1,i_k_point)-KS_eigenvalue(i_state,1,i_k_point)) & 
        *theta_virtocc*(theta_focc-theta_fvirt)/(KS_eigenvalue(j_state,1,i_k_point)-KS_eigenvalue(i_state,1,i_k_point)) & 
        -temp_H(j_state,i_state)&
        !*2*theta_virtocc*(theta_focc-theta_fvirt)/(KS_eigenvalue(j_state,1,i_k_point)-KS_eigenvalue(i_state,1,i_k_point))
        *theta_virtocc*(theta_focc-theta_fvirt)/(KS_eigenvalue(j_state,1,i_k_point)-KS_eigenvalue(i_state,1,i_k_point))

       else ! Switch to the limit when the difference between 2 eigenvalues is too small
         lim = -0.5*DFPT_width*1./(1+cosh((chemical_potential-KS_eigenvalue(i_state,1,i_k_point))/DFPT_width))
         
         first_order_U_complex( j_state,i_state, i_k_task)=         &
         -Omega_MO(j_state,i_state,i_k_task)*lim&
         -temp_H(j_state,i_state)*lim
       endif

     endif ! Extension to metals

      !endif ! i.ne.j
      enddo
   enddo

  endif ! i_k_task   
  enddo ! i_k_point

  if (myid.eq.0) then
   if (warning) then
    write(use_unit,*) 'WARNING: This system has some fractional occupation numbers. This may result in absurd values for the polarizability. Consider using the keyword DFPT_width in the control.in file. For example DFPT_width 0.01 is a sensible setting.'
   end if
  end if


  deallocate(temp_first_order)
  deallocate(temp_1)
  deallocate(temp_2)
  deallocate(temp_H)
  deallocate(temp_eigenvector)

end subroutine evaluate_first_order_U_dielectric
!*****
