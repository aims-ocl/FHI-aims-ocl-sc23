!****s* FHI-aims/evaluate_dipole_polarisability
!  NAME
!   evaluate_dipole_polarisability
!  SYNOPSIS

      subroutine evaluate_dipole_polarisability &
           ( n_occ,n_unocc,n_homo,n_first, &
             n_KS_states,occ_numbers, omega_n, &
             KS_eigenvalue, ovlp_3KS, polar_freq, &
             dipole_mom, mp2_dipole_polar, &
             rpa_dipole_polar &
            )


!  PURPOSE
!  evaluate the dipole polarisability from the dipole moments of
!  KS state pairs.
!
!
!  USES

      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi

      implicit none

!  ARGUMENTS

      integer :: n_occ
      integer :: n_unocc
      integer :: n_homo(n_spin)
      integer :: n_first(n_spin)
      integer :: n_KS_states

      real*8  occ_numbers(n_states,n_spin)
      real*8  omega_n
      real*8  KS_eigenvalue(n_states,n_spin)
      real*8  ovlp_3KS(n_loc_prodbas, n_states,n_KS_states,n_spin)
      real*8  polar_freq(n_basbas, n_loc_prodbas)
      real*8  dipole_mom(n_occ,n_unocc,n_spin,3)

      real*8  mp2_dipole_polar
      real*8  rpa_dipole_polar
      real*8  mp2_dipole_polar_matr(3,3)

!  INPUTS
!  o n_occ -- number of occupied states, n_occ = max(n_homo)
!  o n_unocc -- number of unoccupied states, n_unocc = n_states - min(n_lumo)
!  o n_homo -- the HOMO level, i.e., the number of occupied state 
!  o n_first -- the first non-fully occupied eigenstate, differs from n_lumo
!    in case of fractional occupation
!  o n_KS_states -- the number of single-particle states to be considered in
!    GW calculation, works here provided n_KS_states >= n_homo
!  o occ_number -- occupation number for each state and each spin
!  o KS_eigenvalue -- real array,
!          the eigenvalues of the single-particle calculation. For DFT calculation,
!          this is the KS eigenvalue, but for HF calculation, this is then the HF
!          eigenvalue
!  o ovlp_3KS -- real array
!          this is the transformed 3-cener overlap integration. Now two
!          orbitals of them are KS ones, and one is the auxiliary basis.
!          Note: for parallel calculations, the auxiliary basis are
!          distribuated among the different processors.
!  o polar_freq -- the calculated non-interacting polarisability in
!          terms of the auxiliary basis. What is actaully stored is 
!          v^(1/2) chi_0 v^(1/2)
!
!  o dipole_mom -- dipole moment associated with products of KS (occupied
!          and unoccupied) state pairs 
!  OUTPUTS 
!  o mp2_dipole_polar -- dipole poalrisability at the mp2 level evaluated at current
!          (imaginary )frequency.
!  o rpa_dipole_polar -- dipole poalrisability at the rpa level evaluated at current
!          (imaginary )frequency.
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

      integer ::  n_order
      parameter (n_order = 5)
      real*8 ::  epsilon_zero
!      parameter (epsilon_zero = 5.d0)
      real*8  zeta
      real*8, dimension(:,:), allocatable :: tmp_dipole_prodbas
      real*8, dimension(:,:), allocatable :: dielec_matr
      real*8, dimension(:,:), allocatable :: polar_order_n
      real*8, dimension(:,:), allocatable :: polar_order_2n
      real*8, dimension(:,:), allocatable :: delta_dielec_matr
      real*8, dimension(:,:), allocatable :: aux_dipole_prodbas
      real*8  term1, term2

      integer n_first_min

!     counters


      integer :: j_state
      integer :: k_state

      integer :: i_prodbas
      integer :: i_prodbas_1
      integer :: i_basbas

      integer :: i_index
      integer :: i_order
      integer :: i_spin
      integer :: i_dim

      integer :: i_task

!     begin work

!      n_first_min=n_states
!      do  i_spin = 1, n_spin
!       n_first_min = min(n_first_min, n_first(i_spin))
!      enddo
!
      allocate (tmp_dipole_prodbas(n_loc_prodbas,3),stat=i_index)
      call check_allocation(i_index, 'tmp_dipole_prodbas            ')

      allocate (aux_dipole_prodbas(n_basbas,3),stat=i_index)
      call check_allocation(i_index, 'aux_dipole_prodbas            ')

      allocate (dielec_matr(n_basbas,n_loc_prodbas),stat=i_index)
      call check_allocation(i_index, 'dielec_matr                   ')

      allocate (delta_dielec_matr(n_basbas,n_loc_prodbas),stat=i_index)
      call check_allocation(i_index, 'delta_dielec_matr             ')


      allocate (polar_order_n(n_basbas,n_loc_prodbas),stat=i_index)
      call check_allocation(i_index, 'polar_order_n                 ')

      allocate (polar_order_2n(n_basbas,n_loc_prodbas),stat=i_index)
      call check_allocation(i_index, 'polar_order_2n            ')


      mp2_dipole_polar=0.d0
      mp2_dipole_polar_matr(:,:) = 0.d0
      tmp_dipole_prodbas(:,:) = 0.d0

      i_index =0
      do i_dim = 1, 3, 1
!      do i_dim_1 = 1, 3, 1
       do i_spin = 1, n_spin
        do j_state = 1, n_homo(i_spin), 1

          do k_state =  n_first(i_spin), n_states, 1

            zeta = 2.d0*(KS_eigenvalue(j_state,i_spin) - &
                        KS_eigenvalue(k_state,i_spin))/ &
                       (( KS_eigenvalue(j_state,i_spin) - &
                        KS_eigenvalue(k_state,i_spin))**2 &
                       +omega_n*omega_n) * &
                      occ_numbers(j_state,i_spin)* &
                      (2.d0/n_spin -occ_numbers(k_state,i_spin))* &
                      dble(n_spin*n_spin)/4.d0

             i_index = k_state - n_states + n_unocc

             mp2_dipole_polar_matr(i_dim,i_dim) =  &
                   mp2_dipole_polar_matr(i_dim,i_dim) + &        
                   zeta*2.d0/dble(n_spin)* &
                   dipole_mom(j_state,i_index,i_spin,i_dim)* &
                   dipole_mom(j_state,i_index,i_spin,i_dim) 
                    

             do i_prodbas = 1, n_loc_prodbas
               tmp_dipole_prodbas(i_prodbas,i_dim) = &
                   tmp_dipole_prodbas(i_prodbas,i_dim) + &
                   zeta*2.d0/dble(n_spin)* &
                   dipole_mom(j_state,i_index,i_spin,i_dim)* &
                   ovlp_3KS(i_prodbas,k_state,j_state,i_spin) 
             enddo
           
!  end of k_state
         enddo

!  end of j_state
        enddo
       enddo
!      enddo
      mp2_dipole_polar = mp2_dipole_polar + mp2_dipole_polar_matr(i_dim,i_dim) 
      enddo
      mp2_dipole_polar = mp2_dipole_polar/3.d0

! estimate the ratio of the first term with respect to second term
      call dgemm('N','N',n_basbas,3,n_loc_prodbas,1.d0,&
                  polar_freq,n_basbas, &
                  tmp_dipole_prodbas,n_loc_prodbas,0.d0,&
                  aux_dipole_prodbas,n_basbas)

      call sync_matrix(aux_dipole_prodbas,n_basbas,3)

      term1=0.d0
      term2=0.d0
      i_task = myid+1
      do i_dim = 1, 3, 1
       do i_prodbas=1, n_loc_prodbas, 1

          i_basbas = map_prodbas(i_prodbas,i_task) 
          if(i_basbas.gt.0) then
            term2=term2+ &
                 tmp_dipole_prodbas(i_prodbas,i_dim) * &
                 aux_dipole_prodbas(i_basbas,i_dim)
            term1 = term1 + tmp_dipole_prodbas(i_prodbas,i_dim)* &
             tmp_dipole_prodbas(i_prodbas,i_dim)
         endif
        enddo
      enddo 
      call sync_real_number(term1)
      call sync_real_number(term2)

      term1=term1/3.d0
      term2=term2/3.d0

      epsilon_zero = term2/term1
!      epsilon_zero = -2.d0
!      if(myid.eq.0) then
!       write(use_unit,'(A,f20.12)') "ratio", epsilon_zero
!      endif
! calculate mp2_polarisability using a Taylor-type exansion
! 1/(1-x) = (1/2)*( 1+ (1+x)/2 +(1+x)^2/4 + (1+x)^3/8+ ...)
! here x = v^(1/2) * chi_0 * v^(1/2)

      i_task=myid+1
      do i_prodbas_1=1, n_loc_prodbas, 1 
        i_basbas=map_prodbas(i_prodbas_1,i_task)
        if(i_basbas .le.0) then
          polar_freq(:,i_prodbas_1) = 0.d0
        endif
      enddo

      dielec_matr(:,:) = 0.d0
      i_task = myid+1
      do i_prodbas = 1, n_loc_prodbas, 1

        i_basbas = map_prodbas(i_prodbas, i_task)
        if(i_basbas.gt.0) then
          dielec_matr(i_basbas,i_prodbas) =  1.d0

          polar_freq(i_basbas,i_prodbas) =  &
             polar_freq(i_basbas,i_prodbas) - epsilon_zero 
        endif
      enddo

      polar_freq(:,:) = polar_freq(:,:)/(1.d0-epsilon_zero)
      dielec_matr(:,:) = dielec_matr(:,:) +polar_freq(:,:)
      polar_order_n(:,:) = polar_freq(:,:)

!      aux_polar_freq(:,:) = polar_freq(:,:)
      do i_order = 1, n_order, 1
!         do i_task = 1, n_tasks, 1
! distribute polar_freq into different processors, for parallel calc.
!          do i_prodbas_1 = 1, n_loc_prodbas, 1
!             do i_prodbas = 1, n_loc_prodbas, 1
!               i_basbas = map_prodbas(i_prodbas, i_task)
!               if(i_basbas.gt.0) then
!                 aux_polar_freq(i_prodbas,i_prodbas_1) = &
!                  polar_freq(i_basbas,i_prodbas_1)
!               endif
!             enddo
!          enddo
!
!          call dgemm('N', 'T', n_basbas, n_loc_prodbas, n_loc_prodbas, 1.d0, &
!                 polar_order_n, n_basbas, aux_polar_freq, n_loc_prodbas, 0.d0, &
!                 polar_order_2n, n_basbas)
!
!          call sync_matrix(polar_order_2n,n_basbas,n_loc_prodbas)
!
!          if(myid.eq.i_task-1) then 
!           polar_updated(:,:) = polar_order_2n(:,:)
!          endif

!   end of loop over i_tasks
!        enddo
!
!        polar_order_n(:,:) = polar_updated(:,:)
!        dielec_matr(:,:) = dielec_matr(:,:) +  polar_order_n(:,:)

        call auxiliary_matrix_multi &
             (n_basbas, n_loc_prodbas, map_prodbas, &
               polar_order_n, polar_order_n, polar_order_2n &
             )

        call auxiliary_matrix_multi &
             ( n_basbas, n_loc_prodbas, map_prodbas, &
               dielec_matr, polar_order_2n, delta_dielec_matr &
             )

         dielec_matr(:,:) = dielec_matr(:,:) + delta_dielec_matr(:,:)
         polar_order_n(:,:) = polar_order_2n(:,:)     

!   end of loop over i_order
      enddo
      dielec_matr(:,:) = dielec_matr(:,:)/(1.d0-epsilon_zero)


      call dgemm('N','N',n_basbas,3,n_loc_prodbas,1.d0,&
                  dielec_matr,n_basbas, &
                  tmp_dipole_prodbas,n_loc_prodbas,0.d0,&
                  aux_dipole_prodbas,n_basbas)

      call sync_matrix(aux_dipole_prodbas,n_basbas,3)

      rpa_dipole_polar = 0.d0
      i_task = myid+1
      do i_dim = 1, 3, 1
        do i_prodbas = 1, n_loc_prodbas, 1

         i_basbas = map_prodbas(i_prodbas, i_task)

          if(i_basbas.gt.0) then
           rpa_dipole_polar = rpa_dipole_polar + &
              tmp_dipole_prodbas(i_prodbas,i_dim) * &
              aux_dipole_prodbas(i_basbas,i_dim) 
!            write(use_unit,'(2I5,3f16.8)')i_dim, i_basbas, rpa_dipole_polar,& 
!                 tmp_dipole_prodbas(i_prodbas,i_dim), &
!                 aux_dipole_prodbas(i_basbas,i_dim)
          endif

        enddo
      enddo
      call sync_real_number(rpa_dipole_polar)

      rpa_dipole_polar = rpa_dipole_polar/3.d0 + mp2_dipole_polar

!      if(myid.eq.0) then
!       write(use_unit,'(5X,f16.8)') rpa_dipole_polar
!      endif

      deallocate(tmp_dipole_prodbas)
      deallocate(dielec_matr)
      deallocate(delta_dielec_matr)
      deallocate(polar_order_n)
      deallocate(polar_order_2n)
      deallocate(aux_dipole_prodbas)
      return
      end subroutine evaluate_dipole_polarisability
!---------------------------------------------------------------------
!**********
      subroutine evaluate_dipole_polarisability_2 &
           ( n_occ,n_unocc,n_homo,n_first, &
             n_KS_states,occ_numbers, omega_n, &
             KS_eigenvalue, ovlp_3KS, polar_freq, &
             dipole_mom, mp2_dipole_polar, &
             rpa_dipole_polar &
            )


!  PURPOSE
!  evaluate the dipole polarisability from the dipole moments of
!  KS state pairs.
!
!
!  USES

      use dimensions
      use prodbas
      use mpi_tasks
      use synchronize_mpi

      implicit none

!  ARGUMENTS

      integer :: n_occ
      integer :: n_unocc
      integer :: n_homo(n_spin)
      integer :: n_first(n_spin)
      integer :: n_KS_states

      real*8  occ_numbers(n_states,n_spin)
      real*8  omega_n
      real*8  KS_eigenvalue(n_states,n_spin)
      real*8  ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
      real*8  polar_freq(max_row_2d, max_col_2d)
      real*8  dipole_mom(n_occ,n_unocc,n_spin,3)

      real*8  mp2_dipole_polar
      real*8  rpa_dipole_polar
      real*8  mp2_dipole_polar_matr(3,3)

!  INPUTS
!  o n_occ -- number of occupied states, n_occ = max(n_homo)
!  o n_unocc -- number of unoccupied states, n_unocc = n_states - min(n_lumo)
!  o n_homo -- the HOMO level, i.e., the number of occupied state
!  o n_first -- the first non-fully occupied eigenstate, differs from n_lumo
!    in case of fractional occupation
!  o n_KS_states -- the number of single-particle states to be considered in
!    GW calculation, works here provided n_KS_states >= n_homo
!  o occ_number -- occupation number for each state and each spin
!  o KS_eigenvalue -- real array,
!          the eigenvalues of the single-particle calculation. For DFT calculation,
!          this is the KS eigenvalue, but for HF calculation, this is then the HF
!          eigenvalue
!  o ovlp_3KS -- real array
!          this is the transformed 3-cener overlap integration. Now two
!          orbitals of them are KS ones, and one is the auxiliary basis.
!          Note: for parallel calculations, the auxiliary basis are
!          distribuated among the different processors.
!  o polar_freq -- the calculated non-interacting polarisability in
!          terms of the auxiliary basis. What is actaully stored is
!          v^(1/2) chi_0 v^(1/2)
!
!  o dipole_mom -- dipole moment associated with products of KS (occupied
!          and unoccupied) state pairs
!  OUTPUTS
!  o mp2_dipole_polar -- dipole poalrisability at the mp2 level evaluated at current
!          (imaginary )frequency.
!  o rpa_dipole_polar -- dipole poalrisability at the rpa level evaluated at current
!          (imaginary )frequency.
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

      integer ::  n_order
      parameter (n_order = 5)
      real*8 ::  epsilon_zero
!      parameter (epsilon_zero = 5.d0)
      real*8  zeta
      real*8, dimension(:,:), allocatable :: tmp_dipole_prodbas
      real*8, dimension(:,:), allocatable :: tmp_dipole_prodbas_r
      real*8, dimension(:,:), allocatable :: tmp_dipole_prodbas_c
      real*8, dimension(:,:), allocatable :: dielec_matr
      real*8, dimension(:,:), allocatable :: polar_order_n
      real*8, dimension(:,:), allocatable :: polar_order_2n
      real*8, dimension(:,:), allocatable :: delta_dielec_matr
      real*8, dimension(:,:), allocatable :: aux_dipole_prodbas
      real*8  term1, term2

      integer n_first_min

!     counters


      integer :: j_state
      integer :: k_state

      integer :: i_prodbas
      integer :: i_prodbas_1
      integer :: i_basbas

      integer :: i_index
      integer :: i_order
      integer :: i_spin
      integer :: i_dim

      integer :: i_task
      integer :: icnt_r, icnt_c

!     begin work

      allocate (tmp_dipole_prodbas(n_basbas,3),stat=i_index)
      call check_allocation(i_index, 'tmp_dipole_prodbas            ')

      allocate (tmp_dipole_prodbas_r(max_row_2d,3),stat=i_index)
      call check_allocation(i_index, 'tmp_dipole_prodbas_r          ')

      allocate (tmp_dipole_prodbas_c(max_col_2d,3),stat=i_index)
      call check_allocation(i_index, 'tmp_dipole_prodbas_c          ')

      allocate (aux_dipole_prodbas(max_row_2d,3),stat=i_index)
      call check_allocation(i_index, 'aux_dipole_prodbas            ')

      allocate (dielec_matr(max_row_2d,max_col_2d),stat=i_index)
      call check_allocation(i_index, 'dielec_matr                   ')

      allocate (delta_dielec_matr(max_row_2d,max_col_2d),stat=i_index)
      call check_allocation(i_index, 'delta_dielec_matr             ')


      allocate (polar_order_n(max_row_2d,max_col_2d),stat=i_index)
      call check_allocation(i_index, 'polar_order_n                 ')

      allocate (polar_order_2n(max_row_2d,max_col_2d),stat=i_index)
      call check_allocation(i_index, 'polar_order_2n            ')


      mp2_dipole_polar=0.d0
      mp2_dipole_polar_matr(:,:) = 0.d0
      tmp_dipole_prodbas(:,:) = 0.d0

      i_index =0
      do i_dim = 1, 3, 1
       do i_spin = 1, n_spin
        do j_state = 1, n_homo(i_spin), 1

          do k_state =  n_first(i_spin), n_states, 1

            zeta = 2.d0*(KS_eigenvalue(j_state,i_spin) - &
                        KS_eigenvalue(k_state,i_spin))/ &
                       (( KS_eigenvalue(j_state,i_spin) - &
                        KS_eigenvalue(k_state,i_spin))**2 &
                       +omega_n*omega_n) * &
                      occ_numbers(j_state,i_spin)* &
                      (2.d0/n_spin -occ_numbers(k_state,i_spin))* &
                      dble(n_spin*n_spin)/4.d0

             i_index = k_state - n_states + n_unocc

             mp2_dipole_polar_matr(i_dim,i_dim) =  &
                   mp2_dipole_polar_matr(i_dim,i_dim) + &
                   zeta*2.d0/dble(n_spin)* &
                   dipole_mom(j_state,i_index,i_spin,i_dim)* &
                   dipole_mom(j_state,i_index,i_spin,i_dim)


             if(own_dim1_o3ks(k_state)==myp1_o3KS .and. &
                own_dim2_o3ks(j_state)==myp2_o3KS) then
               do i_prodbas = 1, n_basbas
                 tmp_dipole_prodbas(i_prodbas,i_dim) = &
                     tmp_dipole_prodbas(i_prodbas,i_dim) + &
                     zeta*2.d0/dble(n_spin)* &
                     dipole_mom(j_state,i_index,i_spin,i_dim)* &
                     ovlp_3KS(i_prodbas,loc_dim1_o3ks(k_state),loc_dim2_o3ks(j_state),i_spin)
               enddo
             endif

!  end of k_state
         enddo

!  end of j_state
        enddo
       enddo
!      enddo
      mp2_dipole_polar = mp2_dipole_polar + mp2_dipole_polar_matr(i_dim,i_dim)
      enddo
      mp2_dipole_polar = mp2_dipole_polar/3.d0
      call sync_matrix(tmp_dipole_prodbas,n_basbas,3)

! estimate the ratio of the first term with respect to second term

      term1=0.d0
      term2=0.d0

      term1 = dot_product(tmp_dipole_prodbas(:,1),tmp_dipole_prodbas(:,1)) &
            + dot_product(tmp_dipole_prodbas(:,2),tmp_dipole_prodbas(:,2)) &
            + dot_product(tmp_dipole_prodbas(:,3),tmp_dipole_prodbas(:,3))

      ! Gather tmp_dipole_prodbas in a distribution that fits to the rows
      ! of the current process and in one which fits to the columns

      icnt_r = 0
      icnt_c = 0
      do i_prodbas = 1, n_basbas
         if(MOD((i_prodbas-1)/nb_aux_2d,nprow_aux_2d)==myprow_aux_2d) then
           icnt_r = icnt_r+1
           tmp_dipole_prodbas_r(icnt_r,:) = tmp_dipole_prodbas(i_prodbas,:)
         endif
         if(MOD((i_prodbas-1)/nb_aux_2d,npcol_aux_2d)==mypcol_aux_2d) then
           icnt_c = icnt_c+1
           tmp_dipole_prodbas_c(icnt_c,:) = tmp_dipole_prodbas(i_prodbas,:)
         endif
      enddo

      call dgemm('N','N',max_row_2d,3,max_col_2d,1.d0,&
                  polar_freq,ubound(polar_freq,1), &
                  tmp_dipole_prodbas_c,ubound(tmp_dipole_prodbas_c,1),0.d0,&
                  aux_dipole_prodbas,ubound(aux_dipole_prodbas,1))

      term2 = dot_product(tmp_dipole_prodbas_r(:,1),aux_dipole_prodbas(:,1)) &
            + dot_product(tmp_dipole_prodbas_r(:,2),aux_dipole_prodbas(:,2)) &
            + dot_product(tmp_dipole_prodbas_r(:,3),aux_dipole_prodbas(:,3))

      call sync_real_number(term2)

      term1=term1/3.d0
      term2=term2/3.d0

      epsilon_zero = term2/term1

! calculate mp2_polarisability using a Taylor-type exansion
! 1/(1-x) = (1/2)*( 1+ (1+x)/2 +(1+x)^2/4 + (1+x)^3/8+ ...)
! here x = v^(1/2) * chi_0 * v^(1/2)

      ! set dielec_matr to unit matrix

      call PDLASET( 'Full', n_basbas, n_basbas, 0.d0, 1.d0, dielec_matr, 1, 1, aux_sc_desc_2d )

      ! Subtract epsilon_zero from diagonal of polar_freq

      polar_freq = polar_freq - epsilon_zero*dielec_matr

      polar_freq(:,:) = polar_freq(:,:)/(1.d0-epsilon_zero)
      dielec_matr(:,:) = dielec_matr(:,:) +polar_freq(:,:)
      polar_order_n(:,:) = polar_freq(:,:)

      do i_order = 1, n_order, 1

         call pdgemm('N','N', n_basbas, n_basbas, n_basbas, 1.d0, &
                     polar_order_n, 1, 1, aux_sc_desc_2d,         &
                     polar_order_n, 1, 1, aux_sc_desc_2d, 0.d0,   &
                     polar_order_2n, 1, 1, aux_sc_desc_2d)

         call pdgemm('N','N', n_basbas, n_basbas, n_basbas, 1.d0, &
                     dielec_matr, 1, 1, aux_sc_desc_2d,           &
                     polar_order_2n, 1, 1, aux_sc_desc_2d, 0.d0,  &
                     delta_dielec_matr, 1, 1, aux_sc_desc_2d)

         dielec_matr(:,:) = dielec_matr(:,:) + delta_dielec_matr(:,:)
         polar_order_n(:,:) = polar_order_2n(:,:)

      enddo
      dielec_matr(:,:) = dielec_matr(:,:)/(1.d0-epsilon_zero)

      call dgemm('N','N',max_row_2d,3,max_col_2d,1.d0,&
                  dielec_matr,ubound(dielec_matr,1), &
                  tmp_dipole_prodbas_c,ubound(tmp_dipole_prodbas_c,1),0.d0,&
                  aux_dipole_prodbas,ubound(aux_dipole_prodbas,1))

      rpa_dipole_polar = dot_product(tmp_dipole_prodbas_r(:,1),aux_dipole_prodbas(:,1)) &
                       + dot_product(tmp_dipole_prodbas_r(:,2),aux_dipole_prodbas(:,2)) &
                       + dot_product(tmp_dipole_prodbas_r(:,3),aux_dipole_prodbas(:,3))

      call sync_real_number(rpa_dipole_polar)

      rpa_dipole_polar = rpa_dipole_polar/3.d0 + mp2_dipole_polar

      deallocate(tmp_dipole_prodbas)
      deallocate(dielec_matr)
      deallocate(delta_dielec_matr)
      deallocate(polar_order_n)
      deallocate(polar_order_2n)
      deallocate(aux_dipole_prodbas)
      return
      end subroutine evaluate_dipole_polarisability_2
!---------------------------------------------------------------------
!**********
