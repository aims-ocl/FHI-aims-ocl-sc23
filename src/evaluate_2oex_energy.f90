!****s* FHI-aims/evaluate_2oex_energy
!  NAME
!   evaluate_2oex_energy
!  SYNOPSIS

      subroutine evaluate_2oex_energy &
       (n_low_state,n_KS_states,n_homo, &
        KS_eigenvalue,ovlp_3KS,e_2oex )

!  PURPOSE
!  Calculate the second order exchange energy.
!
!  Note: the usual MP2 energy contain the direct (Hartree-like) part
!  and exchange part. This is the exchange part.

!  USES

      use constants
      use dimensions
      use prodbas
      use mpi_tasks
      use mpi_utilities
      use synchronize_mpi
      use runtime_choices
      use localorb_io, only: use_unit

      implicit none

      integer :: n_homo(n_spin)
      integer :: n_low_state
      integer :: n_KS_states
      real*8  :: KS_eigenvalue(n_states,n_spin,n_k_points)
      real*8  :: ovlp_3KS(n_loc_prodbas,n_states,n_KS_states,n_spin)
      real*8  :: e_2oex

! INPUTS
! o  n_low_state -- the first KS state to be included in frozen-core
!                   calculation
! o  n_homo -- integer array, the HOMO level for each spin channel
! o  n_KS_states -- integer number,
! o          the highest KS/HF eigenstate, the same as the n_high_state.
! o          In the present case, n_high_state >= n_homo should do fine.
! o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
! o  ovlp_3KS -- real array
!            this is the transformed 3-cener overlap integration. Now two orbitals of
!            them are KS ones, and one is the auxiliary basis.
!            Note: for parallel calculations, the auxiliary basis are distribuated
!            among the different processors.
! OUTPUT
! o  e_2oex -- real number, the calculated second-order exchange energy
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

      real*8 :: aux

      real*8, dimension(:,:), allocatable :: aux_eri

      real*8 :: E_mp2_a, E_mp2_b
      real*8 :: e_deno


      integer ::  n_lumo(n_spin)
      integer ::  n_lumo_min

      integer ::  n_shell_aux(n_species)

!     Accuracy
!      real*8 :: nu=1.d-12
!      complex*16 :: inu=(0.d0,1.d-20)



!     Integer counters
      integer i_spin
      integer i_spin2
      integer i_state
      integer i_state2
      integer i_virt
      integer i_virt2
      integer i_basis_1
      integer i_species
      integer i_atom
      integer term
      integer i_empty(n_species)
      integer i_task
      integer i_index
      integer n_unoccupied
      integer i_start_b

!     Error variable
      integer mpierr


      if (myid.eq.0) then
          write(use_unit,'(A)') " --------------------------------------------------------------------------" 
          write(use_unit,'(A)') " | Evaluating 2nd order exchange energy ..."
          write(use_unit,'(A)') 
      endif

      n_lumo(:) = n_homo(:) + 1
   
      n_lumo_min = min(n_lumo(1),n_lumo(n_spin))

      n_unoccupied = n_states-n_lumo_min+1

      e_2oex = 0.d0

      if (.not.allocated(aux_eri)) then
          allocate(aux_eri(n_lumo_min:n_states,n_lumo_min:n_states))
      endif

! -------------------------------------------------------------------
!    NON-Spin polarized systems
! -------------------------------------------------------------------

     if (n_spin.eq.1) then

          i_index = 0
          do i_state=n_low_state, n_homo(n_spin),1

             do i_state2=i_state,n_homo(n_spin),1
!             do i_state2=1,n_homo(n_spin),1

              call dgemm('T', 'N', n_unoccupied, n_unoccupied, &
                      n_loc_prodbas, 1.0d0, &
                      ovlp_3KS(:,n_lumo_min:n_states,i_state,n_spin), &
                      n_loc_prodbas, &
                      ovlp_3KS(:,n_lumo_min:n_states,i_state2, &
                      n_spin), n_loc_prodbas, 0.d0, &
                      aux_eri(n_lumo_min:n_states, &
                              n_lumo_min:n_states), &
                      n_unoccupied &
                     )

              call sync_matrix( aux_eri(n_lumo_min:n_states, &
                                n_lumo_min:n_states),n_unoccupied, &
                                n_unoccupied &
                               )

              i_index = i_index + 1
! MPI task distribution
              if(myid.eq.MOD(i_index,n_tasks)) then

               do i_virt=n_lumo(n_spin), n_states,1
                do i_virt2=n_lumo(n_spin),n_states,1
!                    if (i_state.eq.i_state2.or.i_virt.eq.i_virt2) &
!                       cycle

                   E_mp2_a=aux_eri(i_virt,i_virt2)
                   E_mp2_b=aux_eri(i_virt2,i_virt)

                   E_mp2_a= - E_mp2_b * E_mp2_a


                  if (i_state.ne.i_state2) then
                     E_mp2_a=E_mp2_a*2.d0
                  endif


                  e_deno=  (KS_eigenvalue(i_state,n_spin, 1)) &
                    +  (KS_eigenvalue(i_state2,n_spin,1) ) &
                    -  (KS_eigenvalue(i_virt,n_spin,1)) &
                    -  (KS_eigenvalue(i_virt2,n_spin,1))


                  if (abs(e_deno).lt.1e-6)then

                      write(use_unit,'(10X,A)') &
                           "****************************************"
                      write(use_unit,'(10X,2A)') "| Warning :", &
                        " too close to degeneracy"
                      write(use_unit,'(10X,A)') &
                           "****************************************"
                   endif

                   e_2oex = e_2oex + E_mp2_a/e_deno


!   close b  virtual state 2
                enddo
!   close a virtual state
               enddo

!   end of MPI distribution
              endif
!   close j state2
              enddo

!   close i state
        enddo

      else

! -------------------------------------------------------------------
!     Spin polarized systems
! -------------------------------------------------------------------

        i_index = 0
        do i_spin=1,n_spin,1

             i_spin2 = i_spin
             do i_state=n_low_state, n_homo(i_spin),1

                do i_state2=n_low_state,n_homo(i_spin2) ,1

                  call dgemm('T', 'N', n_unoccupied, n_unoccupied, &
                         n_loc_prodbas, 1.0d0, &
                         ovlp_3KS(:,n_lumo_min:n_states,i_state, &
                                       i_spin), n_loc_prodbas, &
                         ovlp_3KS(:,n_lumo_min:n_states,i_state2, &
                                  i_spin2), n_loc_prodbas, 0.d0, &
                         aux_eri(n_lumo_min:n_states, &
                                 n_lumo_min:n_states), &
                         n_unoccupied &
                        )

                  call sync_matrix( aux_eri(n_lumo_min:n_states, &
                                  n_lumo_min:n_states),n_unoccupied, &
                                  n_unoccupied &
                                 )

               i_index = i_index + 1
! MPI task distribution
               if(myid.eq.MOD(i_index,n_tasks)) then

               do i_virt=n_lumo(i_spin), n_states,1
                do i_virt2=n_lumo(i_spin2),n_states,1
!                    if (i_state.eq.i_state2.or.i_virt.eq.i_virt2) cycle

                        E_mp2_a=aux_eri(i_virt,i_virt2)
                        E_mp2_b=aux_eri(i_virt2,i_virt)

                        E_mp2_a=-E_mp2_b*E_mp2_a


                        e_deno=  KS_eigenvalue(i_state,i_spin, 1) &
                           +  KS_eigenvalue(i_state2,i_spin2,1)

                        e_deno=  e_deno - &
                                KS_eigenvalue(i_virt,i_spin,1)

                        e_deno=  e_deno - &
                                KS_eigenvalue(i_virt2,i_spin2,1)




                       if (abs(e_deno).lt.1e-6)then

                          write(use_unit,'(10X,A)') &
                           "****************************************"
                          write(use_unit,'(10X,2A)') "| Warning :", &
                              "too close to degeneracy"
                          write(use_unit,'(10X,A)') &
                           "****************************************"
                       endif


                       e_2oex= E_mp2_a/e_deno+ e_2oex


!   close b  virtual state 2
                enddo
!   close a virtual state
               enddo
!   end of MPI distribution
              endif
!   close j state2
              enddo

!   close i state
         enddo

       enddo
       e_2oex = 0.5d0* e_2oex

      endif

      call sync_real_number(e_2oex) 

      if (allocated(aux_eri)) then
        deallocate (aux_eri)
      endif

      if(myid.eq.0) then

        write(use_unit,'(2X,A,f19.8,A,f19.8,A)') &
         " Second order exchange energy : ",  &
          e_2oex, '  Ha. ',  e_2oex*hartree,  ' eV '
        write(use_unit,*)
      endif

      return

      end subroutine evaluate_2oex_energy


!******
!---------------------------------------------------------------------------------------------------
!****s* FHI-aims/evaluate_2oex_energy
!  NAME
!   evaluate_2oex_energy
!  SYNOPSIS

      subroutine evaluate_2oex_energy_2 &
       (n_low_state,n_KS_states,n_homo, &
        KS_eigenvalue,ovlp_3KS,e_2oex )

!  PURPOSE
!  Calculate the second order exchange energy.
!
!  Note: the usual MP2 energy contain the direct (Hartree-like) part
!  and exchange part. This is the exchange part.

!  USES

      use constants
      use dimensions
      use prodbas
      use mpi_tasks
      use mpi_utilities
      use synchronize_mpi
      use runtime_choices
      use localorb_io, only: use_unit

      implicit none

      integer :: n_homo(n_spin)
      integer :: n_low_state
      integer :: n_KS_states
      real*8  :: KS_eigenvalue(n_states,n_spin,n_k_points)
      real*8  :: ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin)
      real*8  :: e_2oex

! INPUTS
! o  n_low_state -- the first KS state to be included in frozen-core
!                   calculation
! o  n_homo -- integer array, the HOMO level for each spin channel
! o  n_KS_states -- integer number,
! o          the highest KS/HF eigenstate, the same as the n_high_state.
! o          In the present case, n_high_state >= n_homo should do fine.
! o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
! o  ovlp_3KS -- real array
!            this is the transformed 3-cener overlap integration. Now two orbitals of
!            them are KS ones, and one is the auxiliary basis.
!            Note: for parallel calculations, dimensions 2 and 3 are distribuated
!            among the different processors.
! OUTPUT
! o  e_2oex -- real number, the calculated second-order exchange energy
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

      real*8, allocatable :: tmp_o3KS(:,:)
      real*8, allocatable :: tmp_eri_full(:,:,:,:)
      real*8, allocatable :: tmp_eri_tran(:,:,:)

      real*8 :: E_mp2_a, E_mp2_b
      real*8 :: e_deno


      integer ::  n_lumo(n_spin)
      integer ::  n_lumo_min
      integer ::  n_homo_max
      integer ::  min_lumo_loc
      integer ::  max_homo_loc
      integer ::  max_local_states1

!     Integer counters
      integer i_spin
      integer i_state
      integer i_state2
      integer i_virt
      integer i_virt2
      integer n_p1, n_p2
      integer i_virt_loc, i_virt2_loc, i_virt2_own
      integer i_state2_loc

!     Error variable
      integer mpierr
real*8 ttt0


      if (myid.eq.0) then
          write(use_unit,'(A)') "---------------------------------------------"
          write(use_unit,'(A)') " | Evaluating 2nd-order exchange energy ..."
          write(use_unit,'(A)') 
      endif

      n_lumo(:) = n_homo(:) + 1
      n_homo_max = maxval(n_homo(:))
      n_lumo_min = minval(n_lumo(:))

      e_2oex = 0.d0

      min_lumo_loc = loc_dim1_o3ks(n_lumo_min) ! overall minimal local value corresponding to n_lumo_min
      max_homo_loc = loc_dim2_o3ks(n_homo_max) ! overall maximal local value corresponding to n_homo_max

      max_local_states1 = ndim1_o3KS - min_lumo_loc + 1

      allocate(tmp_o3KS(n_basbas,min_lumo_loc:ndim1_o3KS))
      allocate(tmp_eri_full(min_lumo_loc:ndim1_o3KS,min_lumo_loc:ndim1_o3KS,0:np1_o3KS-1,max_homo_loc))
      allocate(tmp_eri_tran(min_lumo_loc:ndim1_o3KS,min_lumo_loc:ndim1_o3KS,0:np1_o3KS-1))
ttt0 = MPI_Wtime()

      do i_spin = 1, n_spin, 1

        do i_state = n_low_state, n_KS_states, 1

          n_p2 = own_dim2_o3ks(i_state)

          do n_p1 = 0, np1_o3KS-1

            ! Proc (n_p1,n_p2) broadcasts its local part of
            ! ovlp_3KS_full(:,n_lumo_min:n_states,i_state,i_spin)

            if(n_p1==myp1_o3KS .and. n_p2==myp2_o3KS) &
              tmp_o3KS(:,:) = ovlp_3KS(:,min_lumo_loc:ndim1_o3KS,loc_dim2_o3ks(i_state),i_spin)

            call MPI_Bcast(tmp_o3KS,n_basbas*max_local_states1,MPI_REAL8,global_id(n_p1,n_p2),mpi_comm_global,mpierr)

            ! Multiply what we got with our local part of ovlp_3KS_full(:,n_lumo_min:n_states,1:n_homo(i_spin),i_spin)

            do i_state2 = 1, n_homo(i_spin), 1

              if(own_dim2_o3ks(i_state2) /= myp2_o3KS) cycle
              i_state2_loc = loc_dim2_o3ks(i_state2)

              call dgemm('T', 'N', max_local_states1, max_local_states1, n_basbas, 1.d0, &
                         tmp_o3KS, ubound(tmp_o3KS,1), &
                         ovlp_3KS(1,min_lumo_loc,i_state2_loc,i_spin), ubound(ovlp_3KS,1), &
                         0.d0, tmp_eri_full(min_lumo_loc,min_lumo_loc,n_p1,i_state2_loc),max_local_states1)
            enddo

          enddo


          do i_state2=i_state,n_homo(i_spin),1

            ! The original code:
            !
            !call dgemm('T', 'N', n_unoccupied, n_unoccupied, &
            !       n_loc_prodbas, 1.0d0, &
            !       ovlp_3KS(:,n_lumo_min:n_states,i_state,n_spin), &
            !       n_loc_prodbas, &
            !       ovlp_3KS(:,n_lumo_min:n_states,i_state2,n_spin), &
            !       n_loc_prodbas,0.d0, &
            !       aux_eri(n_lumo_min:n_states,n_lumo_min:n_states), &
            !       n_unoccupied &
            !      )

            if(own_dim2_o3ks(i_state2) /= myp2_o3KS) cycle
            i_state2_loc = loc_dim2_o3ks(i_state2)

            ! Transpose the strip we have in tmp_eri_full:

            call MPI_Alltoall(tmp_eri_full(min_lumo_loc,min_lumo_loc,0,i_state2_loc), &
                              max_local_states1*max_local_states1,MPI_REAL8, &
                              tmp_eri_tran,max_local_states1*max_local_states1,MPI_REAL8, &
                              mpi_comm_rows_aux_2d,mpierr)

            do i_virt=n_lumo(i_spin),n_states,1

              if(own_dim1_o3ks(i_virt) /= myp1_o3ks) cycle

              do i_virt2=n_lumo(i_spin),n_states,1

                i_virt_loc  = loc_dim1_o3ks(i_virt)
                i_virt2_own = own_dim1_o3ks(i_virt2)
                i_virt2_loc = loc_dim1_o3ks(i_virt2)

                E_mp2_a=tmp_eri_full(i_virt2_loc,i_virt_loc,i_virt2_own,i_state2_loc)
                E_mp2_b=tmp_eri_tran(i_virt_loc,i_virt2_loc,i_virt2_own)
                E_mp2_a= - E_mp2_b * E_mp2_a

                if (i_state.ne.i_state2) E_mp2_a=E_mp2_a*2.d0

                e_deno = KS_eigenvalue(i_state ,i_spin,1) &
                       + KS_eigenvalue(i_state2,i_spin,1) &
                       - KS_eigenvalue(i_virt ,i_spin,1)  &
                       - KS_eigenvalue(i_virt2,i_spin,1)


                if (abs(e_deno).lt.1e-6)then
                  write(use_unit,'(10X,A)') "****************************************"
                  write(use_unit,'(10X,A)') "| Warning : too close to degeneracy"
                  write(use_unit,'(10X,A)') "****************************************"
                endif

                e_2oex = e_2oex + E_mp2_a/e_deno

              enddo ! i_virt2
            enddo ! i_virt
          enddo ! i_state2
! if(myid==0) print *,i_state,'B',MPI_Wtime()-ttt0
        enddo ! i_state
      enddo ! i_spin

      deallocate(tmp_o3KS)
      deallocate(tmp_eri_full)
      deallocate(tmp_eri_tran)

      e_2oex = e_2oex / n_spin

      call sync_real_number(e_2oex) 

      if(myid.eq.0) then

        write(use_unit,'(2X,A,f19.8,A,f19.8,A)') &
         " Second order exchange energy : ",  &
          e_2oex, '  Ha. ',  e_2oex*hartree,  ' eV '
        write(use_unit,*)
      endif

      return

      end subroutine evaluate_2oex_energy_2


!******
