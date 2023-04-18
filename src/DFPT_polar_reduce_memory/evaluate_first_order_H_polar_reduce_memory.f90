!****s* FHI-aims/evaluate_first_order_H_polar_reduce_memory
!  NAME
!    evaluate_first_order_H_polar_reduce_memory
!  SYNOPSIS

subroutine evaluate_first_order_H_polar_reduce_memory & 
          (first_order_H,n_points, &
           partition_tab, grid_coord,       &
           H_times_psi, n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave, & 
           first_order_rho,v_hartree_gradient,dVxc_drho,& 
           vsigma, v2rho2, v2rhosigma, v2sigma2, & 
           gradient_rho, first_order_gradient_rho, n_matrix_size)

!  PURPOSE
!    calculate the first-order Hamiltion matrix elements.
!    five terms:
!   (1) <X0| -r                 |X0>   ----------
!   (2) <X0| Int{rho(1)/|r-r'|} |X0>     Vscf(1)
!   (3) <X0| dVxc/drho * rho(1) |X0>   ----------

!  shanghui,2012.03.05
!  shanghui,2012.04.23
!  shanghui,2012.05.29 : to complie 
!  shanghui, 2013.12.09 for polarizability
 
!  USES

  use dimensions
  use species_data
  use geometry
  use runtime_choices
  use basis 
  use localorb_io, only: use_unit
  use pbc_lists, only: index_hamiltonian, column_index_hamiltonian 

  ! wyj
  use load_balancing
  use grids
  use mpi_tasks

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) :: partition_tab
  real*8, dimension(n_points), intent(in) :: grid_coord

  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham,n_points,n_spin), intent(in) :: H_times_psi

  real*8, dimension( n_spin, n_points), intent(in) :: first_order_rho
  real*8, dimension( n_points), intent(in) :: v_hartree_gradient

  !-------LDA----------------
  real*8, dimension(3,n_points), intent(in) :: dVxc_drho

  !-------GGA----------------
  real*8, dimension(3,n_points), intent(in) :: v2rho2  ! = dVxc_drho
  real*8, dimension(3,n_points), intent(in) :: vsigma  
  real*8, dimension(6,n_points), intent(in) :: v2rhosigma
  real*8, dimension(6,n_points), intent(in) :: v2sigma2 
  real*8, dimension(3 ,n_spin, n_points), intent(in) :: gradient_rho
  real*8, dimension(3, n_spin, n_points), intent(in) :: first_order_gradient_rho
 
  !real*8, dimension(n_hamiltonian_matrix_size, n_spin), intent(inout) :: first_order_H
  !real*8, dimension(*, n_spin), intent(inout) :: first_order_H
  !real*8, intent(inout) :: first_order_H_arg(*)
  ! wyj
  real*8, dimension(n_matrix_size, n_spin), intent(inout) :: first_order_H
  integer n_matrix_size
!  INPUTS
!  o  rho -- electron density
!  o  wave --
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!
!  OUTPUT
!  first_order_H 


  integer :: i_point, i_basis,j_basis,i_compute,j_compute,i_spin,i_coords2
  real*8 :: point_term
  integer :: i_start, i_end, i_place, i_index_real

  real*8 :: basis_basis, gradient_basis_basis(3),  & 
            gradient_rho_gradient_basis_basis(n_spin), & 
            first_order_gradient_rho_gradient_basis_basis(3), & 
            first_order_sigama(3) ,&
            first_order_sigama2(n_points,n_spin), prefactor1(n_points), prefactor2(n_points), prefactor3(n_points)
            

  real*8, external :: ddot

  real*8 :: contract(n_points,n_compute_c), wave_t(n_points,n_max_compute_ham),&
            contract_grad(n_points,n_compute_c), &
            first_order_H_dense(n_compute_c, n_compute_c, n_spin)

  real*8, allocatable     :: tmp_H(:,:), tmp_H2(:,:)


  ! wyj add
  integer :: i, j ,i_off, n_bp
  integer, allocatable :: ins_idx(:)
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

  ! wyj
  n_bp = use_batch_permutation
  if (use_batch_permutation > 0) then
      batches_work => batch_perm(n_bp)%batches
      allocate(ins_idx(batch_perm(n_bp)%n_basis_local))
  endif

  if(.not.use_gga) then ! LDA case 
!------------------------------------LDA------------------------------------------------------
  if(n_spin.eq.1) then 
     i_spin = 1 

    first_order_H_dense=0.d0

    wave_t=transpose(wave) ! Exchange dimensions


    do i_compute=1,n_compute_c
      do i_point=1,n_points
          contract(i_point,i_compute)=partition_tab(i_point)*wave_t(i_point,i_compute)*&
          (-grid_coord(i_point)+v_hartree_gradient(i_point) &
           +dVxc_drho(i_spin,i_point)*first_order_rho(i_spin,i_point))
      enddo
    enddo

   !! wyj: TODO debug
   !print *, myid, 'H_wave_t=', wave_t(1:3,1)
   !wave_t(:,:) = 1.0d0
   !contract(:,:) = 1.0d0

    call dgemm("T","N",n_compute_c,n_compute_c,n_points,&
               1.d0,contract,n_points,&
               wave_t,n_points,0.d0,first_order_H_dense(:,:,i_spin),n_compute_c)


!---------------begin dense to sparse matrix for cluster case--------------------

   !-------(1) begin the slow version-------------------
   !do i_compute=1,n_compute_c,1
   !   i_basis = i_basis_index(i_compute)
   !   i_start = index_hamiltonian(1,1, i_basis)
   !   i_end   = index_hamiltonian(2,1, i_basis)

   !! shanghui note: here j_compute run over all the n_compute_c, which is
   !! different from the update_full_matrix_p0.f90. If use the full n_compute_c,
   !! then we do not use if-else  inside the i_place loop.  
   !do j_compute=1,n_compute_c,1   
   !   j_basis = i_basis_index(j_compute)
   !     
   !   do i_place = i_start, i_end, 1 

   !      if(column_index_hamiltonian(i_place) == j_basis) then
   !         first_order_H(i_place,i_spin) =    &
   !         first_order_H(i_place,i_spin) + first_order_H_dense(i_compute,j_compute,i_spin) 
   !      endif
   !   
   !   enddo

   !enddo
   !enddo
   !-------(1) end the slow version-------------------
 
   !! wyj: TODO debug
   !print *, myid, 'H_dense=', first_order_H_dense(1:3,1,1)
   !first_order_H_dense(:,:,:) = 1.0d0
    
   ! wyj: THIS PART are correct[PASS]
   if (.not. (use_local_index .and. use_load_balancing)) then 
   !-------(2) begin the faster versin---------------
   !note: the same as the update_full_matrix_p0.f90
     do i_compute = 1, n_compute_c, 1

        i_start =  index_hamiltonian(1,1, i_basis_index(i_compute))
        i_end   =  index_hamiltonian(2,1, i_basis_index(i_compute))

        do j_compute = 1,i_compute,1
           
           j_basis = i_basis_index(j_compute)
           
           place: do i_place = i_start, i_end, 1
              
              if( column_index_hamiltonian( i_place) == j_basis)then
                    
                 if (j_compute.le.i_compute) then
                    first_order_H(i_place, i_spin) = first_order_H(i_place, i_spin) + & 
                    first_order_H_dense(j_compute, i_compute, i_spin)
                 else
                    first_order_H(i_place, i_spin) = first_order_H(i_place, i_spin) + &
                    first_order_H_dense(i_compute, j_compute, i_spin)
                 end if

                 i_index_real = i_place
                 exit place 

              else if(column_index_hamiltonian( i_place) > j_basis)then
                 i_index_real = i_place
                 exit place  

              end if
           end do place
           i_start = i_index_real
        end do
        end do
        !-------(2) end the faster versin---------------

    else
        !print *, myid, 'begin evaluate_first_order_H 0'
        ! wyj: TODO, this part will affect the DM change!
        if(use_batch_permutation > 0) then
            ! If use_batch_permutation > 0 is set, the local hamiltonian is always stored in full form for the local basis functions
            ! Get position of basis functions of current batch within local hamiltonian
            do i=1,n_compute_c
                ins_idx(i) = batch_perm(n_bp)%i_basis_glb_to_loc(i_basis_index(i))
            enddo

            ! Insert hamiltonian_shell of current batch
            do i=1,n_compute_c
            !i_off = (i_spin-1)*ld_hamiltonian + (ins_idx(i)*(ins_idx(i)-1))/2
                i_off = (ins_idx(i)*(ins_idx(i)-1))/2
                do j=1,i ! n_compute_c
                !hamiltonian(ins_idx(j)+i_off) = hamiltonian(ins_idx(j)+i_off) + hamiltonian_shell(j+(i-1)*n_compute_c)
                first_order_H(ins_idx(j)+i_off, i_spin) = first_order_H(ins_idx(j)+i_off, i_spin) + first_order_H_dense(j,i,i_spin)
                if (1 == 2) then
                    print *,"i = ", i, " j = ", j, " ShellIdx = ", j+(i-1)*n_compute_c, " MatIdx = ", ins_idx(j)+i_off, &
                        " Shell = ",first_order_H_dense(j,i,i_spin), " Matrix = ",  first_order_H(ins_idx(j)+i_off, i_spin)
                end if
                enddo
            enddo
        endif
    endif


!---------------end dense to sparse matrix for cluster case--------------------


  else !n_spin=2
  
   first_order_H_dense=0.d0

   do i_compute=1,n_compute_c,1
   do j_compute=1,n_compute_c,1
  
   do i_point = 1, n_points, 1
  
      point_term = partition_tab(i_point) * wave(i_compute,i_point) * wave(j_compute,i_point) 
  
           ! i_spin = 1
           ! up= dK2/(drho_u drho_u) drho_u/dE + dK2/(drho_u drho_d) drho_d/dE 
           first_order_H_dense(i_compute,j_compute,1) =    &
           first_order_H_dense(i_compute,j_compute,1) +    & 
           point_term*(-grid_coord( i_point))   +    & !(1) 
           point_term*v_hartree_gradient(i_point)    +    & !(2)
           point_term*dVxc_drho(1,i_point)*first_order_rho(1,i_point) + & !(3)
           point_term*dVxc_drho(2,i_point)*first_order_rho(2,i_point)

           ! i_spin = 2
           ! down= dK2/(drho_d drho_d) drho_d/dE + dK2/(drho_u drho_d) drho_u/dE 
           first_order_H_dense(i_compute,j_compute,2) =    &
           first_order_H_dense(i_compute,j_compute,2) +    & 
           point_term*(-grid_coord(i_point))   +    & !(1) 
           point_term*v_hartree_gradient(i_point)    +    & !(2)
           point_term*dVxc_drho(3,i_point)*first_order_rho(2,i_point) + & !(3)
           point_term*dVxc_drho(2,i_point)*first_order_rho(1,i_point)
  
  
   end do ! n_points
  
   enddo
   enddo


!---------------begin dense to sparse matrix for cluster case--------------------
    do i_spin = 1, n_spin
    do i_compute=1,n_compute_c,1
       i_basis = i_basis_index(i_compute)
       i_start = index_hamiltonian(1,1, i_basis)
       i_end   = index_hamiltonian(2,1, i_basis)

    do j_compute=1,n_compute_c,1
       j_basis = i_basis_index(j_compute)
         
       do i_place = i_start, i_end, 1 

          if(column_index_hamiltonian(i_place) == j_basis) then
             first_order_H(i_place,i_spin) =    &
             first_order_H(i_place,i_spin) + first_order_H_dense(i_compute,j_compute,i_spin) 
          endif
       
       enddo

    enddo
    enddo
    enddo ! i_spin
!---------------end dense to sparse matrix for cluster case--------------------


  endif !n_spin

  else if(use_gga) then! GGA case 
!------------------------------------GGA------------------------------------------------------
 
  if(n_spin.eq.1) then 
     i_spin = 1

! 28.06.19 Nath: new structure for faster evaluations

  wave_t=transpose(wave) ! Exchange dimensions

    allocate(tmp_H(n_compute_c, n_compute_c)) 
    allocate(tmp_H2(n_compute_c, n_compute_c))

    tmp_H = 0.d0
    tmp_H2 = 0.d0
    contract_grad = 0.0

  do i_point = 1, n_points
    ! Nath: This use of ddot is very slow here...must depend on the position of the argument which is multiplied (here in 2nd position)
    !first_order_sigama(i_spin) = 2.0d0*ddot(3, first_order_gradient_rho(1:3,i_spin,i_point), 1, &
    !                                   gradient_rho(1:3,i_spin,i_point),1)
    first_order_sigama2(i_point,i_spin) = 0.0
    do i_coords2=1,3
      first_order_sigama2(i_point,i_spin) = first_order_sigama2(i_point,i_spin) &
       +2.0d0*first_order_gradient_rho(i_coords2,i_spin,i_point)*gradient_rho(i_coords2,i_spin,i_point)
    enddo

    ! Put together i_point-dependent factors
    prefactor1(i_point) = &
    partition_tab(i_point) *&
    (- grid_coord(i_point) + v_hartree_gradient(i_point) &
     + v2rho2(1,i_point)*first_order_rho(i_spin, i_point) &
     + v2rhosigma(1,i_point)*first_order_sigama2(i_point,i_spin) )

    prefactor2(i_point) = partition_tab(i_point)*(2.0d0*v2rhosigma(1,i_point)*first_order_rho(i_spin,i_point) &
                          + 2.0d0*v2sigma2(1,i_point)*first_order_sigama2(i_point,i_spin))

    prefactor3(i_point) = partition_tab(i_point)*2.0d0*vsigma(1,i_point)

    do i_compute = 1, n_compute_c

      contract(i_point,i_compute) = prefactor1(i_point)*wave(i_compute,i_point)

      !----------(3.2) gradient_basis_basis term---------------
      do i_coords2 = 1, 3
        !----------(3.2.1) gradient_rho term--------------------
        contract_grad(i_point,i_compute) = contract_grad(i_point,i_compute)  &
        + prefactor2(i_point)*gradient_basis_wave(i_compute,i_coords2,i_point)*gradient_rho(i_coords2,i_spin,i_point) &
        !----------(3.2.2) first_order_gradient_rho term--------------------
        + prefactor3(i_point)*gradient_basis_wave(i_compute,i_coords2,i_point) &
          *first_order_gradient_rho(i_coords2,i_spin,i_point)
      enddo


    enddo ! i_compute
  enddo ! i_point

  call dgemm("T","N",n_compute_c,n_compute_c,n_points,&
             1.d0,contract,n_points,&
             wave_t,n_points,0.d0,tmp_H,n_compute_c)

  call dgemm("T","N",n_compute_c,n_compute_c,n_points,&
             1.d0,contract_grad,n_points,&
             wave_t,n_points,0.d0,tmp_H2,n_compute_c)

  tmp_H = tmp_H + tmp_H2 + transpose(tmp_H2)

  ! "Reconstruct" matrix
!---------------begin dense to sparse matrix for cluster case--------------------
    do i_compute=1,n_compute_c,1
       i_basis = i_basis_index(i_compute)
       i_start = index_hamiltonian(1,1, i_basis)
       i_end   = index_hamiltonian(2,1, i_basis)

    do j_compute=1,n_compute_c,1
       j_basis = i_basis_index(j_compute)
         
       do i_place = i_start, i_end, 1 

          if(column_index_hamiltonian(i_place) == j_basis) then
             first_order_H(i_place,i_spin) =    &
             first_order_H(i_place,i_spin) + tmp_H(i_compute,j_compute) 
          endif
       
       enddo

    enddo
    enddo
!---------------end dense to sparse matrix for cluster case--------------------



  deallocate(tmp_H) 
  deallocate(tmp_H2) 

! ENd new structure



  else  !n_spin = 2 for GGA


   first_order_H_dense=0.d0

   do i_compute=1,n_compute_c,1
   do j_compute=1,n_compute_c,1
  
  
   do i_point = 1, n_points, 1
  
         basis_basis = partition_tab(i_point) * wave(i_compute,i_point) * wave(j_compute,i_point) 
         gradient_basis_basis(1:3) = partition_tab(i_point) * ( & 
                                     gradient_basis_wave(i_compute,1:3,i_point)*wave(j_compute,i_point) + &
                                     gradient_basis_wave(j_compute,1:3,i_point)*wave(i_compute,i_point) )
 
         do i_spin = 1,n_spin 
         gradient_rho_gradient_basis_basis(i_spin) = ddot(3, gradient_rho(1:3,i_spin,i_point), 1, & 
                                                          gradient_basis_basis(1:3),1)
         enddo  

 
           do i_spin = 1, n_spin 
              first_order_gradient_rho_gradient_basis_basis(i_spin) = & 
                                               ddot(3, first_order_gradient_rho(1:3,i_spin,i_point), 1, &
                                                       gradient_basis_basis(1:3),1)
           enddo 
 
           first_order_sigama(1) = 2.0d0*ddot(3, first_order_gradient_rho(1:3,1,i_point), 1, &
                                              gradient_rho(1:3,1,i_point),1)
           first_order_sigama(2) = ddot(3, first_order_gradient_rho(1:3,1,i_point), 1, &
                                              gradient_rho(1:3,2,i_point),1) + & 
                                   ddot(3, first_order_gradient_rho(1:3,2,i_point), 1, & 
                                              gradient_rho(1:3,1,i_point),1) 
           first_order_sigama(3) = 2.0d0*ddot(3, first_order_gradient_rho(1:3,2,i_point), 1, &
                                              gradient_rho(1:3,2,i_point),1)
 
         !===================i_spin=1================================
           first_order_H_dense(i_compute,j_compute,1) =    &
           first_order_H_dense(i_compute,j_compute,1) +    & 
           basis_basis*(-grid_coord(i_point))   +    & !(1) 
           basis_basis*v_hartree_gradient(i_point)    +    & !(2)
           !----------(3.1) bassis_basis term-----------------------
           v2rho2(1,i_point)*basis_basis*first_order_rho(1, i_point) + &
           v2rho2(2,i_point)*basis_basis*first_order_rho(2, i_point) + &
           v2rhosigma(1,i_point)*basis_basis*first_order_sigama(1) + &
           v2rhosigma(2,i_point)*basis_basis*first_order_sigama(2) + &
           v2rhosigma(3,i_point)*basis_basis*first_order_sigama(3) + &
           !----------(3.2) gradient_basis_basis term---------------
           !----------(3.2.1) gradient_rho term--------------------
          ( 2.0d0*v2rhosigma(1,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_rho(1,i_point) + &
                  v2rhosigma(2,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_rho(1,i_point)) + &
          ( 2.0d0*v2rhosigma(4,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_rho(2,i_point) + &
                  v2rhosigma(5,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_rho(2,i_point)) + &
          ( 2.0d0*v2sigma2(1,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_sigama(1) + & 
                  v2sigma2(2,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_sigama(1)) + &
          ( 2.0d0*v2sigma2(2,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_sigama(2) + &   
                  v2sigma2(4,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_sigama(2)) + &
          ( 2.0d0*v2sigma2(3,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_sigama(3) + &
                  v2sigma2(5,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_sigama(3)) + &
          !----------(3.2.2) first_order_gradient_rho term--------------------
          ( 2.0d0*vsigma(1,i_point)*first_order_gradient_rho_gradient_basis_basis(1) + &
                  vsigma(2,i_point)*first_order_gradient_rho_gradient_basis_basis(2)) 

         !===================i_spin=2================================
           first_order_H_dense(i_compute,j_compute,2) =    &
           first_order_H_dense(i_compute,j_compute,2) +    & 
           basis_basis*(-grid_coord(i_point))   +    & !(1) 
           basis_basis*v_hartree_gradient(i_point)    +    & !(2)
           !----------(3.1) bassis_basis term-----------------------
           v2rho2(2,i_point)*basis_basis*first_order_rho(1, i_point) + &
           v2rho2(3,i_point)*basis_basis*first_order_rho(2, i_point) + &
           v2rhosigma(4,i_point)*basis_basis*first_order_sigama(1) + &
           v2rhosigma(5,i_point)*basis_basis*first_order_sigama(2) + &
           v2rhosigma(6,i_point)*basis_basis*first_order_sigama(3) + &
          !----------(3.2) gradient_basis_basis term---------------
          !----------(3.2.1) gradient_rho term--------------------
          ( 2.0d0*v2rhosigma(3,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_rho(1,i_point) + &
                  v2rhosigma(2,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_rho(1,i_point)) + &
          ( 2.0d0*v2rhosigma(6,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_rho(2,i_point) + &
                  v2rhosigma(5,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_rho(2,i_point)) + &
          ( 2.0d0*v2sigma2(3,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_sigama(1) + &       
                  v2sigma2(2,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_sigama(1)) + &
          ( 2.0d0*v2sigma2(5,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_sigama(2) + &       
                  v2sigma2(4,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_sigama(2)) + &
          ( 2.0d0*v2sigma2(6,i_point)*gradient_rho_gradient_basis_basis(2)*first_order_sigama(3) + &  
                  v2sigma2(5,i_point)*gradient_rho_gradient_basis_basis(1)*first_order_sigama(3)) + &
          !----------(3.2.2) first_order_gradient_rho term--------------------
          ( 2.0d0*vsigma(3,i_point)*first_order_gradient_rho_gradient_basis_basis(2) + &
                  vsigma(2,i_point)*first_order_gradient_rho_gradient_basis_basis(1)) 

         
  
   end do ! n_points
  
   enddo
   enddo



!---------------begin dense to sparse matrix for cluster case--------------------
    do i_spin = 1, n_spin
    do i_compute=1,n_compute_c,1
       i_basis = i_basis_index(i_compute)
       i_start = index_hamiltonian(1,1, i_basis)
       i_end   = index_hamiltonian(2,1, i_basis)

    do j_compute=1,n_compute_c,1
       j_basis = i_basis_index(j_compute)
         
       do i_place = i_start, i_end, 1 

          if(column_index_hamiltonian(i_place) == j_basis) then
             first_order_H(i_place,i_spin) =    &
             first_order_H(i_place,i_spin) + first_order_H_dense(i_compute,j_compute,i_spin) 
          endif
       
       enddo

    enddo
    enddo
    enddo ! i_spin
!---------------end dense to sparse matrix for cluster case--------------------




  endif !n_spin 

  else 

  
    write(use_unit,'(2X,A)') "DFPT_polar_reduce_memory is only support LDA, GGA at present." 
    stop

  end if ! LDA, GGA 

  ! wyj
  if (allocated(ins_idx)) deallocate(ins_idx) 
end subroutine evaluate_first_order_H_polar_reduce_memory
!******
