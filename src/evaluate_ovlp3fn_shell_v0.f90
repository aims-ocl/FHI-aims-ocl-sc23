!****s* FHI-aims/evaluate_ovlp3fn_shell_v0
!  NAME
!   evaluate_ovlp3fn_shell_v0
!  SYNOPSIS

      subroutine evaluate_ovlp3fn_shell_v0 &
           ( n_points,i_task,i_atom,partition_3atoms, &
             n_prod_compute, i_prodbas, &
             n_prodbas_current, i_prodbas_current, &
             n_prodbas_other, i_prodbas_other, &
             n_compute, i_basis, &
             n_compute_current, i_basis_current, &
             n_compact, n_compact_current, &
             compact_pairs, compact_pairs_current, &
             wave, wave_prod, &
             ovlp_3fn &
           )

!  PURPOSE
!  Subroutine evaluate_ovlp3basis_shell evaluates the three basis  overlap
!  integral contribution of a radial shell (or a batch), and adds it to the overall
!  3-center overlap matrix.
!
!  USES
      use dimensions
      use basis
      use prodbas
      use mpi_tasks

      implicit none

!  ARGUEMNTS

      integer :: n_points
      integer :: i_task
      integer :: i_atom


      integer n_compute
      integer i_basis(n_compute)
      integer n_compute_current
      integer i_basis_current(n_compute_current)

      integer n_prod_compute
      integer n_prodbas_current
      integer n_prodbas_other
      integer i_prodbas(n_prod_compute)
      integer i_prodbas_current(n_prodbas_current)
      integer i_prodbas_other(n_prodbas_other)

      integer n_compact(n_basis)
      integer compact_pairs(n_basis,n_basis)
      integer n_compact_current(n_basis)
      integer compact_pairs_current(n_max_basis_atom,n_basis)

      real*8 partition_3atoms(n_points,n_atoms*(n_atoms+1)/2)
      real*8 wave(n_compute,n_points)
      real*8 wave_prod(n_loc_prodbas,n_points)

      real*8 ovlp_3fn(n_basis_pairs, n_loc_prodbas)

!  INPUTS
!  o n_points -- number of relevant points in the current integration shell or batch   
!  o i_task -- the current task number
!  o i_atom -- the current atom under consideration
!  o n_compute -- the number of non-zero basis function
!  o i_basis -- array, specifies these non-zero basis function
!  o n_compute_current -- the non-zero basis funciton belonging to the current atom
!  o i_basis_current -- specifies the non-zero basis funciton belonging to the 
!           current atom
!  o n_prod_compute -- the number of non-zero auxiliary basis function
!  o n_prodbas_current -- the number of non-zero auxiliary basis function assoicated
!             with the current atom
!  o n_prodbas_other -- the number of non-zero auxiliary basis funciton not associated
!             with the current atom
!  o i_prodbas -- specifies all the non-zero auxiliary basis fucntion
!  o i_prodbas_current -- specifies non-zero auxiliary basis function assoicated
!             with the current atom
!  o i_prodbas_other -- specifies non-zero auxiliary basis funciton not associated
!             with the current atom
!  o n_compact -- array, for each basis in the batch, the number of basis  
!           functions in the same batch which has zero overlap with the one under 
!           consideration
!  o compact_pairs -- two-dimentional array, for each basis in the batch, 
!           list the basis functions in the same batch which has zero overlap with 
!           the one under consideration
!  o n_compact_current -- array, for each basis in the batch, the number of basis  
!           functions in the same batch that belong to the current atom and has
!           zero overlap with the one under consideration
!  o compact_pairs_current -- two-dimentional array, for each basis in the batch, 
!           list the basis functions in the same batch which has zero overlap with 
!           the one under consideration and belong to the curren atom
!  o partition_3atoms -- the value of the 3-center partition function
!  o wave -- the values of the basis functions in the batch
!  o wave_prod -- the values of the auxiliary basis functions in the batch
!  OUTPUTS
!  o ovlp_3fn -- 3-center integral over two regular NAO basis and one auxiliary basis
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

!     auxiliary matrices for Level 3 Blas matrix multiplications

! ??? Need we store partition * wave in an extra matrix? ??? need BLAS descriptions to go on ...
! ??? Does it matter that n_points (the summation index in the matrix multiplication) is
! ??? the outer one, not the inner one, in wave?
!     Yes, we need, but as long as we are inside the cache, no harm should follow...


      real*8, allocatable :: psi_times_psi(:,:)
      real*8, allocatable :: chi_times_psi(:,:)

!VB: if these are supposed to be allocatable, they have to be allocated outside ONCE as
!    work arrays, not thousands of times here!

      real*8, allocatable, dimension(:,:) :: aux_ovlp3fn_matrix
      real*8, allocatable, dimension(:,:) :: aux_ovlp3fn_matrix_1

      real*8, allocatable :: wave_prod_compute(:,:)
      real*8, allocatable :: wave_compute(:,:)
      real*8, allocatable :: wave_compute_current(:,:)

      real*8, allocatable :: partition(:,:)

!     counters

      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_basis_3
      integer :: i_basbas
      integer :: i_index

      integer :: i_atom_1
      integer :: i_atom_2
      integer :: i_atom_index

      integer :: i_compute
      integer :: i_compute_1
      integer :: i_compute_2
      integer :: i_compact
      integer :: i_prodbas_2
      integer :: i_offset
      integer :: i_point
      integer :: info

      character(*), parameter :: func = 'evaluate_ovlp3fn_shell_v0'


!     begin work

      allocate(psi_times_psi(n_points,n_compute), stat=info)
      call check_allocation(info, 'psi_times_psi', func)
      allocate(chi_times_psi(n_points,n_prodbas_other), stat=info)
      call check_allocation(info, 'chi_times_psi', func)
      allocate(wave_prod_compute(n_points,n_prod_compute), stat=info)
      call check_allocation(info, 'wave_prod_compute', func)
      allocate(wave_compute(n_points,n_compute), stat=info)
      call check_allocation(info, 'wave_compute', func)
      allocate(wave_compute_current(n_points,n_compute_current), stat=info)
      call check_allocation(info, 'wave_compute_current', func)
      allocate(partition(n_points,n_atoms*(n_atoms+1)/2), stat=info)
      call check_allocation(info, 'partition', func)

      allocate(aux_ovlp3fn_matrix(n_compute, n_prodbas_current), stat=info)
      call check_allocation(info, 'aux_ovlp3fn_matrix', func)
      allocate(aux_ovlp3fn_matrix_1(n_prodbas_other, n_compute_current), stat=info)
      call check_allocation(info, 'aux_ovlp3fn_matrix_1', func)

!  First run: here we multiply the product basis on the present atom
!  with all the possible pairs of regular basis functions.

! not needed - full array filled anyway       wave_compute(:,:) = 0.d0
       do i_compute = 1, n_compute, 1
          wave_compute(1:n_points, i_compute) = &
           wave(i_compute, 1:n_points)
       enddo

       if(n_prodbas_current .gt.0 ) then

! not needed - all relevant parts of array are overwritten        wave_prod_compute(:,:) = 0.d0
        do i_compute = 1, n_prodbas_current, 1
          wave_prod_compute(1:n_points, i_compute) = &
           wave_prod(i_prodbas_current(i_compute), 1:n_points)
        enddo


        i_index = 0
        do i_compute =1, n_compute, 1

         i_basis_1 = i_basis(i_compute)
         i_atom_1 = basis_atom(i_basis_1)
         do i_compact = 1, n_compact(i_compute), 1

           i_compute_2 = compact_pairs(i_compact,i_compute)
           i_atom_2=basis_atom(i_basis(i_compute_2))

!           if(i_atom_1.gt.i_atom_2) then
!             i_atom_index = (i_atom_1-1)*i_atom_1/2 + i_atom_2
!           else
!             i_atom_index = (i_atom_2-1)*i_atom_2/2 + i_atom_1
!           endif

!  "atom_pairs" is defined in prodbas.f and obtained in condense_basis_pairs.
!  this is to eliminate the expensive "if" here.

           i_atom_index = atom_pairs(i_atom_2,i_atom_1)

           do i_point = 1, n_points
             psi_times_psi (i_point,i_compact) = &
                  wave_compute(i_point, i_compute) * &
                  wave_compute(i_point, i_compute_2) * &
                  partition_3atoms(i_point, i_atom_index)
           enddo
          enddo

!not needed, done by dgemm!          aux_ovlp3fn_matrix(:,:) = 0.d0
             call dgemm('T', 'N', n_compact(i_compute), &
                      n_prodbas_current, n_points, 1.0d0, &
                      psi_times_psi, n_points, &
                      wave_prod_compute, n_points, 0.d0, &
                      aux_ovlp3fn_matrix, n_compute &
                     )

             do i_compact = 1, n_compact(i_compute), 1

                 i_compute_2 = compact_pairs(i_compact,i_compute)
                 i_basis_2 = i_basis(i_compute_2)
                 i_index = basis_nghr(i_basis_2,i_basis_1)

              if(i_index.gt.0) then
          do i_compute_1 = 1, n_prodbas_current, 1

             i_basis_3 = i_prodbas(i_prodbas_current(i_compute_1))

                 ovlp_3fn(i_index, i_basis_3) = &
                      ovlp_3fn(i_index, i_basis_3) + &
                      aux_ovlp3fn_matrix(i_compact, i_compute_1)
             enddo
               endif
          enddo
       enddo
!   end of if n_prodbas_current
      endif

!  Second run: here we multiply one regualr basis on the present atom
!  with all the possible pairs of one arbitrary regular basis and one product basis
!  on other atoms.
      if (n_compute_current.gt.0.and.n_prodbas_other.gt.0) then

!        wave_compute(:,:) = 0.d0
!        do i_compute = 1, n_compute, 1
!           wave_compute(1:n_points,i_compute) =
!     +        wave(i_compute, 1:n_points)
!        enddo

!VB: should not be needed!        wave_compute_current(:,:) = 0.d0

!        do i_compute_1 = 1, n_compute_current, 1
!           i_compute = i_basis_current(i_compute_1)
!           wave_compute_current(1:n_points,i_compute_1) =
!     +      wave_compute(1:n_points,i_compute)
!        enddo

! not needed - all relevant parts of array are overwritten        wave_prod_compute(:,:) = 0.d0
        do i_compute = 1, n_prodbas_other, 1
          wave_prod_compute(1:n_points, i_compute) = &
           wave_prod(i_prodbas_other(i_compute), 1:n_points)
        enddo


        do i_compute = 1, n_compute, 1

           i_basis_1=i_basis(i_compute)
           i_atom_1 = basis_atom(i_basis_1)

! VB: Should this really be initialized here????
!     there's that if statement - I'll leave it untouched.
          chi_times_psi(:,:) = 0.d0
          do i_compute_1 = 1, n_prodbas_other, 1

            i_prodbas_2 = i_prodbas_other(i_compute_1)
            i_basis_3 = i_prodbas(i_prodbas_2)
            i_basbas = map_prodbas(i_basis_3,i_task)
            if(i_basbas.gt.0) then

               i_atom_2 = basbas_atom(i_basbas)
               i_atom_index = atom_pairs(i_atom_2, i_atom_1)

               do i_point = 1, n_points
                 chi_times_psi (i_point, i_compute_1) = &
                    wave_prod_compute(i_point,i_compute_1) * &
                    wave_compute(i_point,i_compute) * &
                    partition_3atoms(i_point,i_atom_index)
                enddo
             endif
           enddo

            do i_compact = 1, n_compact_current(i_compute), 1

               i_compute_2 = compact_pairs_current(i_compact,i_compute)
                 wave_compute_current(1:n_points,i_compact) = &
                   wave_compute(1:n_points, i_compute_2)
            enddo

!VB: not needed, dgemm takes care of that!            aux_ovlp3fn_matrix_1(:,:) = 0.d0

!VB: BUG!!! n_compact_current(i_compact) should be
!           n_compact_current(i_compute) !!!

            call dgemm('T', 'N', n_prodbas_other, &
                      n_compact_current(i_compute), n_points, &
                      1.0d0, chi_times_psi, n_points, &
                      wave_compute_current, n_points, 0.0d0, &
                      aux_ovlp3fn_matrix_1,n_prodbas_other &
                     )

            do i_compact = 1, n_compact_current(i_compute), 1

              i_compute_2 = compact_pairs_current(i_compact,i_compute)
              i_basis_2 = i_basis(i_compute_2)
              i_atom_2 = basis_atom(i_basis_2)

              if( i_atom_2.eq.i_atom_1 .and. &
                         i_basis_2.ne.i_basis_1 ) then
                aux_ovlp3fn_matrix_1(:, i_compact) = &
                0.5d0*  aux_ovlp3fn_matrix_1(:, i_compact)
              endif

              i_index = basis_nghr(i_basis_2,i_basis_1)
              if(i_index .gt. 0) then
                do i_compute_1 = 1, n_prodbas_other, 1

                 i_prodbas_2 = i_prodbas_other(i_compute_1)
                 i_basis_3 = i_prodbas(i_prodbas_2)

                 ovlp_3fn(i_index, i_basis_3) = &
                      ovlp_3fn(i_index, i_basis_3) + &
                      aux_ovlp3fn_matrix_1(i_compute_1, i_compact)
!     end of loop over i_compute_1
               enddo
!     end of if i_index >0
             endif
!     end of loop over i_compact
           enddo

!     end of loop over i_compute
        enddo
!     end of if n_compute_current > 0
      endif

      deallocate(psi_times_psi)
      deallocate(chi_times_psi)
      deallocate(wave_prod_compute)
      deallocate(wave_compute)
      deallocate(wave_compute_current)
      deallocate(partition)

      deallocate(aux_ovlp3fn_matrix)
      deallocate(aux_ovlp3fn_matrix_1)

      return
      end subroutine evaluate_ovlp3fn_shell_v0
!---------------------------------------------------------------------
!******
