!****s* FHI-aims/evaluate_first_zero_order_DM_EDM_p0
!  NAME
!    evaluate_first_zero_order_DM_EDM_p0
!  SYNOPSIS

subroutine evaluate_first_zero_order_DM_EDM_p0(first_order_S_complex,  &
           first_order_H_complex,  &
           KS_eigenvector_complex, KS_eigenvalue, occ_numbers,  &
           first_order_U_complex,density_matrix,first_order_density_matrix, & 
           energy_density_matrix,first_order_energy_density_matrix)
               

!  PURPOSE
!  (1) calculate zero and first-order DM for PBC case 
!  DM(uR1, vR2)(0)=sum{k} sum{i} Cui*(k) Cvi(k) exp(-ik(R2-R1))
!  DM(uR1, vR2)(1)=sum{k} sum{i} [Cui*(k)(1) Cvi(k)+ Cui*(k) Cvi(k)(1)] exp(-ik(R2-R1))
!  shanghui,2012.10.19

!  (2) calculate zero and first-order EDM for PBC case for phonon
!  EDM()(0)= 
!  EDM()(1)=
!  shanghui,2013.11.28 


!  (3) paralle over k
!  n_k_points : occ_numbers, KS_eigenvalue,k_phase
!  n_k_task   :  KS_eigenvector_complex 


!  USES

  use dimensions
  use pbc_lists
  use synchronize_mpi

!  ARGUMENTS

  implicit none
  
  complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task),  intent(IN) :: KS_eigenvector_complex
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers

  complex*16, dimension(3,n_atoms, n_basis, n_basis,n_k_points_task),intent(IN) :: first_order_S_complex
  complex*16, dimension(3,n_atoms, n_basis, n_basis,n_k_points_task),intent(IN) :: first_order_H_complex
  complex*16, dimension(3,n_atoms, n_basis,n_basis,n_k_points_task), intent(IN) :: first_order_U_complex
  real*8, dimension(n_Cbasis,n_Cbasis),intent(INOUT) :: density_matrix
  real*8, dimension(3, n_atoms, n_Cbasis,n_Cbasis),intent(INOUT) :: & 
                                              first_order_density_matrix

  real*8, dimension(n_Cbasis,n_Cbasis),intent(INOUT) :: energy_density_matrix
  real*8, dimension(3, n_atoms, n_Cbasis,n_Cbasis),intent(INOUT) :: & 
                                              first_order_energy_density_matrix

!  INPUTS
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o first_order_S_complex -- first order S at every K  
!   o first_order_H_complex -- first order H at every K
!   o first_order_U_complex -- first order U at every K

!
!  OUTPUT
!   o density_matrix
!   o energy_density_matrix
!   o first_order_density_matrix
!   o first_order_energy_density_matrix


  integer :: i_basis,j_basis ,i_coord, i_atom,i_spin
  integer :: i_state,j_state,i_k_point,i_k_task
  integer, dimension(n_spin, n_k_points_task)   :: max_occ_number

  complex*16 , dimension(n_basis, n_basis) ::temp_1 
  complex*16 , dimension(n_basis, n_basis) ::temp_H1,temp_S1
  complex*16 , dimension(n_basis, n_basis) ::temp_eigenvector
  complex*16 , dimension(n_basis, n_basis) ::Ck_Ck,Ck_Ck_Ei,Ck_S1_Ck,Ck_H1_Ck
  complex*16 , dimension(n_basis, n_basis) ::ij_term,ia_term,ij_term_Ei,ia_term_Ei



!=======================zero order DM and EDM================
   density_matrix(1:n_Cbasis,1:n_Cbasis)=0.0d0
   energy_density_matrix(1:n_Cbasis,1:n_Cbasis)=0.0d0
 
    i_k_task = 0 
 
    do i_k_point = 1, n_k_points,1
   if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
       i_k_task = i_k_task + 1

!-------------- Find the max_occ_number
        do i_spin = 1, n_spin, 1
            max_occ_number(i_spin,i_k_task) = 0
            do i_state = n_states, 1, -1
             if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.0.d0) then
              max_occ_number(i_spin,i_k_task) = i_state
              exit
             endif
            enddo
        enddo


     Ck_Ck(1:n_basis,1:n_basis)=(0.0d0,0.0d0) 
     Ck_Ck_Ei(1:n_basis,1:n_basis)=(0.0d0,0.0d0) 

     do i_basis=1,n_basis
       do j_basis=1,n_basis
        do i_state=1, max_occ_number(1,i_k_task)
         
            Ck_Ck(i_basis,j_basis)=Ck_Ck(i_basis,j_basis)       &
          + dconjg(KS_eigenvector_complex(i_basis,i_state,1,i_k_task)) & 
          * occ_numbers(i_state,1,i_k_point)                            &
          * KS_eigenvector_complex(j_basis,i_state,1,i_k_task )              

            Ck_Ck_Ei(i_basis,j_basis)=Ck_Ck_Ei(i_basis,j_basis)       &
          + dconjg(KS_eigenvector_complex(i_basis,i_state,1,i_k_task)) & 
          * occ_numbers(i_state,1,i_k_point)                            &
          * KS_eigenvector_complex(j_basis,i_state,1,i_k_task )     & 
          * KS_eigenvalue(i_state,1,i_k_point )         

        enddo
       enddo
     enddo

     do i_basis=1, n_Cbasis
        do j_basis=1, n_Cbasis
   density_matrix(i_basis,j_basis)= &
   density_matrix(i_basis,j_basis)  &
          + dble( Ck_Ck(Cbasis_to_basis(i_basis),Cbasis_to_basis(j_basis))  & 
          * dconjg(k_phase(center_to_cell(Cbasis_to_center(j_basis)), i_k_point)) & 
          * k_phase(center_to_cell(Cbasis_to_center(i_basis)), i_k_point) )   & 
          * k_weights(i_k_point)
          ! exp(-k(R2-R1)) 

   energy_density_matrix(i_basis,j_basis)= &
   energy_density_matrix(i_basis,j_basis)  &
          + dble( Ck_Ck_Ei(Cbasis_to_basis(i_basis),Cbasis_to_basis(j_basis))  & 
          * dconjg(k_phase(center_to_cell(Cbasis_to_center(j_basis)), i_k_point)) & 
          * k_phase(center_to_cell(Cbasis_to_center(i_basis)), i_k_point) )   & 
          * k_weights(i_k_point)
          ! exp(-k(R2-R1)) 
        enddo
     enddo 

   endif
   enddo ! i_k_point


  call sync_density_matrix(density_matrix) 
  call sync_density_matrix(energy_density_matrix) 


!=======================first order DM and EDM================
  first_order_density_matrix(1:3,1:n_atoms,1:n_Cbasis,1:n_Cbasis)=0.0d0
  first_order_energy_density_matrix(1:3,1:n_atoms,1:n_Cbasis,1:n_Cbasis)=0.0d0

  do i_coord=1,3 
  do i_atom=1,n_atoms

!--------------------the ij_term-------------------------------
  !temp2=-CtS(1)C----->
 
  i_k_task = 0
 
  do i_k_point = 1, n_k_points,1
  if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
       i_k_task = i_k_task + 1

    do i_basis=1,n_basis
    temp_S1(i_basis,1:n_basis)=          &
    -first_order_S_complex(i_coord, i_atom, i_basis,1:n_basis, i_k_task)

    temp_H1(i_basis,1:n_basis)=          &
    first_order_H_complex(i_coord, i_atom, i_basis,1:n_basis, i_k_task)

    temp_eigenvector(i_basis,1:n_basis)=          &
    KS_eigenvector_complex(i_basis,1:n_basis,1,i_k_task)
    enddo

!----------for Ck_S1_Ck--------------------
    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("C","N",n_basis,n_basis,n_basis,&
               1.0d0, temp_eigenvector, n_basis,temp_S1, n_basis,&
               0.0d0,temp_1,n_basis)

    Ck_S1_Ck(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                1.0d0,temp_1, n_basis,temp_eigenvector, n_basis,&
                0.0d0,Ck_S1_Ck, n_basis)

!----------for Ck_H1_Ck--------------------
    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("C","N",n_basis,n_basis,n_basis,&
               1.0d0, temp_eigenvector, n_basis,temp_H1, n_basis,&
               0.0d0,temp_1,n_basis)

    Ck_H1_Ck(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                1.0d0,temp_1, n_basis,temp_eigenvector, n_basis,&
                0.0d0,Ck_H1_Ck, n_basis)



     do i_basis=1,n_basis
     do j_basis=1,n_basis

        ij_term(i_basis,j_basis)=(0.0d0,0.0d0) 
        ij_term_Ei(i_basis,j_basis)=(0.0d0,0.0d0) 

        do i_state=1, max_occ_number(1,i_k_task)
        do j_state=1, max_occ_number(1,i_k_task)
         
            ij_term(i_basis,j_basis)=ij_term(i_basis,j_basis)       &
          + dconjg(KS_eigenvector_complex(i_basis,j_state,1,i_k_task)) & 
          * occ_numbers(i_state,1,i_k_point)                            &
          * Ck_S1_Ck(i_state,j_state)                                     &
          * KS_eigenvector_complex(j_basis,i_state,1,i_k_task )              

            ij_term_Ei(i_basis,j_basis)=ij_term_Ei(i_basis,j_basis)       &
          + dconjg(KS_eigenvector_complex(i_basis,j_state,1,i_k_task)) & 
          * occ_numbers(i_state,1,i_k_point)                            &
          * ( Ck_H1_Ck(i_state,j_state)+ & 
              Ck_S1_Ck(i_state,j_state)  & 
              *( KS_eigenvalue(i_state,1,i_k_point)  &    
                +KS_eigenvalue(j_state,1,i_k_point)) ) &  
          * KS_eigenvector_complex(j_basis,i_state,1,i_k_task )              

        enddo
        enddo
   
     enddo
     enddo

!-----------------------end ij_term----------------------------

  
!-----------------------the ia_term-------------------------------
     do i_basis=1,n_basis
        do j_basis=1,n_basis

         ia_term(i_basis,j_basis)=(0.0d0,0.0d0)
         ia_term_Ei(i_basis,j_basis)=(0.0d0,0.0d0)

         do i_state=1, max_occ_number(1,i_k_task)
         do j_state=max_occ_number(1,i_k_task)+1,n_states

           ia_term(i_basis,j_basis)= ia_term(i_basis,j_basis) &
           + dconjg( KS_eigenvector_complex(i_basis,j_state,1,i_k_task)    &
           * first_order_U_complex(i_coord, i_atom, j_state,i_state,i_k_task) )      &
           * occ_numbers(i_state,1,i_k_point)                               &
           * KS_eigenvector_complex(j_basis,i_state,1,i_k_task )           &
           + dconjg( KS_eigenvector_complex(i_basis,i_state,1,i_k_task))   &
           * occ_numbers(i_state,1,i_k_point)                               &
           * KS_eigenvector_complex(j_basis,j_state,1,i_k_task )           &
           * first_order_U_complex(i_coord, i_atom, j_state,i_state,i_k_task)
              
           ia_term_Ei(i_basis,j_basis)= ia_term_Ei(i_basis,j_basis) &
           + dconjg( KS_eigenvector_complex(i_basis,j_state,1,i_k_task)    &
           * first_order_U_complex(i_coord, i_atom, j_state,i_state,i_k_task) )   &
           * occ_numbers(i_state,1,i_k_point)*KS_eigenvalue(i_state,1,i_k_point)  &
           * KS_eigenvector_complex(j_basis,i_state,1,i_k_task )           &
           + dconjg( KS_eigenvector_complex(i_basis,i_state,1,i_k_task))   &
           * occ_numbers(i_state,1,i_k_point)*KS_eigenvalue(i_state,1,i_k_point)  &
           * KS_eigenvector_complex(j_basis,j_state,1,i_k_task )           &
           * first_order_U_complex(i_coord, i_atom, j_state,i_state,i_k_task)              
         enddo
         enddo 
    
         enddo
     enddo
!-----------------------end the ia_term----------------------------

     do i_basis=1, n_Cbasis
        do j_basis=1, n_Cbasis
   first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)= &
   first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)  &
          + dble( (ij_term(Cbasis_to_basis(i_basis),Cbasis_to_basis(j_basis))  & 
                  +ia_term(Cbasis_to_basis(i_basis),Cbasis_to_basis(j_basis))) & 
          * dconjg(k_phase(center_to_cell(Cbasis_to_center(j_basis)), i_k_point)) & 
          * k_phase(center_to_cell(Cbasis_to_center(i_basis)), i_k_point) )   & 
          * k_weights(i_k_point) 


   first_order_energy_density_matrix(i_coord,i_atom,i_basis,j_basis)= &
   first_order_energy_density_matrix(i_coord,i_atom,i_basis,j_basis)  &
          + dble( (ij_term_Ei(Cbasis_to_basis(i_basis),Cbasis_to_basis(j_basis))  & 
                  +ia_term_Ei(Cbasis_to_basis(i_basis),Cbasis_to_basis(j_basis))) & 
          * dconjg(k_phase(center_to_cell(Cbasis_to_center(j_basis)), i_k_point)) & 
          * k_phase(center_to_cell(Cbasis_to_center(i_basis)), i_k_point) )   & 
          * k_weights(i_k_point) 

        enddo
     enddo 
 
    endif
    enddo ! i_k_point

   enddo ! i_atom
   enddo ! i_coord

  !-------shanghui begin parallel------
   call sync_frist_order_density_matrix(first_order_density_matrix)
   call sync_frist_order_density_matrix(first_order_energy_density_matrix)
  !-------shanghui end parallel------

end subroutine evaluate_first_zero_order_DM_EDM_p0
!******
