module plus_u

! PURPOSE
!
! This module contains all DFT+U related routines.
!
!   M.J.Han, T.Ozaki, and J.Yu, Phys. Rev. B 73, 045110 (2006).

! NOTES:
! If multiple sets of radial basis functions are present e.g. 
! atomic 3d and tier1 3d basis functions. The user can build up
! the hubbard projectors as a linear combination of these basis functions.
! A different possibility would be to treat each set as a seperate shell, however,
! this requires to switch of the on site orthogonalization of the basis function
! in order to get physicall meaningful hubbard projectors. 

!  USES
   use physics
   use pbc_lists   

   implicit none

   ! variable definitions

   logical, public   :: plus_u_petukhov_mixing_defined         ! if true double counting is determined self consistently
   logical, public   :: plus_u_ramping_defined 
   logical, public   :: plus_u_matrix_release_defined
   logical, public   :: plus_u_eigenvalues_defined             ! calculates the eigenvalues of the occupation matrix
   logical, public   :: plus_u_occupation_matrix_control_read
   logical, public   :: plus_u_occupation_matrix_control_write 
   logical, public   :: plus_u_hydros_defined                  ! use hubbard correction also for tier1, tier2, ... functions
   logical, public   :: plus_u_matrix_error_defined            ! flag for calculating the idempotence error of the occupation matrix
   logical, public   :: occ_mat_file_exists                    ! flag to read in the occupation matrix from a fail
   real*8,  public   :: plus_u_matrix_release                  ! threshold for matrix release
   real*8,  public   :: plus_u_ramping_accuracy                ! total energy threshold for DFT+U ramping
   real*8,  public   :: plus_u_petukhov_mixing                 ! 0=>AMF, 1=>FLL, otherwise it will be determined self consistently
   real*8,  public   :: plus_u_energy_correction               ! E_U - \sum_i f_i < Psi_i | v_U | Psi_i >

   public allocate_plus_u, plus_u_init_idx, deallocate_plus_u, occ_numbers_plus_u, add_plus_u_to_hamiltonian, &
          plus_u_energy_correction_term, plus_u_write_occ_mat, plus_u_read_occ_mat, plus_u_matrix_error, plus_u_check_occupation, &
          plus_u_eigenvalues, occ_numbers_plus_u_mulliken, add_plus_u_mulliken_to_hamiltonian, &
          occ_numbers_plus_u_full, add_plus_u_full_to_hamiltonian, plus_u_ramping_work

   private

   integer                 :: l_max                         ! angular quantum number of plus_U shell with highst l

   integer                 :: num_shell_max                 ! number of plus_U shells of the species with the most plus_u shells

   integer                 :: n_plusUorb                    ! number of relevant shells within unit cell
   integer, allocatable    :: plusUorb_index(:)             ! (n_plusUorb)
   real*8,  allocatable    :: occ_plus_u(:,:,:,:)           ! (n_plusUorb,2*l_max+1,2*l_max+1,n_spin)
   logical                 :: occ_numbers_initialized
   logical                 :: ramping_initialized
   real*8,  allocatable    :: hubbard_u(:)                  ! (n_plusUorb)
   real*8,  allocatable    :: hubbard_u_target(:)           ! (n_plusUorb) target values for DFT+U ramping
   real*8,  allocatable    :: hubbard_u_increment(:)

   integer, allocatable    :: num_l_shell(:)                ! number of basisfkt shells used for to
   integer, allocatable    :: p_max(:)                      ! include hydrogenic basis fkt.

!#endif
   real*8,  allocatable    :: av_occ_array(:,:)             ! average occupation number of occupation matrix 
   real*8,  allocatable    :: AMF_FLL_array(:)


   real*8,  allocatable    :: hubc(:,:,:)                   ! predefined coefficients for hubbard projectors

   integer, allocatable    :: n_idx(:,:,:,:,:)              ! arrays with sparse indecies
   integer, allocatable    :: ridx(:,:,:,:,:,:),cidx(:,:,:,:,:,:)

   integer                 :: ramping_counter

contains

!-----------------------------------------------------------------------------------------
!------------------------------HELPER-FUNCTIONS-------------------------------------------
!-----------------------------------------------------------------------------------------


   ! only the upper triangular matrix is stored, this 
   ! functions accounts for that
   real*8 function overlap_matrix_w(i,j)
      integer,intent(in) :: i,j
      if(i<=j)then
         overlap_matrix_w = overlap_matrix(i+j*(j-1)/2)
      else
         overlap_matrix_w = overlap_matrix(j+i*(i-1)/2)
      endif
   end function overlap_matrix_w


   ! input:  matrix indecies
   ! output: corresponding sparse index 
   integer function get_idx_sparse(i_cell,ix,jx)
      integer,intent(in) :: i_cell,ix,jx
      integer :: i,j,mn,md,mx
      i=ix;j=jx
      if(j<i)then
         i=jx;j=ix
      endif
      get_idx_sparse = 0
      if (index_hamiltonian(1, i_cell, j) == 0)return
      mn=index_hamiltonian(1, i_cell, j)
      mx=index_hamiltonian(2, i_cell, j)
      if(column_index_hamiltonian(mn) >= i)then
         if(column_index_hamiltonian(mn) == i)get_idx_sparse = mn
         return
      endif
      if(column_index_hamiltonian(mx) <= i)then
         if(column_index_hamiltonian(mx) == i)get_idx_sparse = mx
         return
      endif
      do
         if(mx-mn <= 1)return
         md=(mn+mx)/2
         if(column_index_hamiltonian(md) == i)then
            get_idx_sparse = md
            return
         elseif(column_index_hamiltonian(md) < i)then
            mn=md
         else
            mx=md
         endif
      enddo
   end function

!-----------------------------------------------------------------------------------------
!--------------------------------ROUTINES-FOR-BOOKKEEPING---------------------------------
!-----------------------------------------------------------------------------------------


   subroutine allocate_plus_u

   !  PURPOSE
   !
   !  allcotes the occupation matrix for each shell and atom.
   !  all book keeping is happening here!

   !  USES

      use dimensions
      use constants, only: hartree
      use basis
      use pbc_lists
      use species_data
      use localorb_io, only: use_unit
      use geometry, only: species
      !use mpi_tasks, only: myid

      implicit none

   !  ARGUMENTS

   !  INPUT

   !  OUTPUT

   !  AUTHOR
   !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
   !  HISTORY
   !    Release version, FHI-aims (2017).

   !  SOURCE

   ! local variables

      integer :: i_basis,i_species,i_orb,i_pass,i_shell, i_fn, k
      integer :: prev_center, prev_shell, prev_fn
      integer :: ialloc

      integer :: counter, i_basis_2, i_fn_2
      integer :: prev_center_2, prev_shell_2, prev_fn_2

      if(allocated(plusUorb_index))  deallocate(plusUorb_index,stat=ialloc) ! if here does not hurt
      if(allocated(occ_plus_u))      deallocate(occ_plus_u,stat=ialloc)     ! important for reinizialisation 
      if(allocated(hubbard_u))       deallocate(hubbard_u,stat=ialloc)      ! of scf
      if(allocated(num_l_shell))     deallocate(num_l_shell,stat=ialloc)
      if(allocated(p_max))           deallocate(p_max,stat=ialloc)
      if(allocated(hubc))            deallocate(hubc,stat=ialloc)

      !-------------------------- labelling radial basis functions ------------------------
      ! Equivalent functions on different atoms are considered different


      do i_pass = 1,2
         ! pass 1: determine number of shells within unit cell
         ! allocate
         ! pass 2: fill allocated arrays

         i_orb = 0 ; prev_center = 0 ; prev_fn = 0
         do i_basis = 1, n_basis
            i_fn = basis_fn(i_basis)   ! i_fn = radial basis function number to i_basis		
            if(plus_u_hydros_defined .eqv. .false.) then
               if(basisfn_type(i_fn) /= 'atomic') cycle
               i_species = species(basis_atom(i_basis))
                  do i_shell = 1, n_max_shells_plus_u
                    if(basisfn_n(i_fn) /= plus_u_n(i_species,i_shell)) cycle
                    if(basis_l(i_basis) /= plus_u_l(i_species,i_shell)) cycle 

                    if(Cbasis_to_center(i_basis) /= prev_center.or.i_shell /= prev_shell) then 
                       ! only the first basis fkt can pass here all remaining basis fkt do not pass this if statement
                       i_orb = i_orb + 1
                       if(i_pass==2)then
                          plusUorb_index(i_orb) = i_basis! point to m = -l entry of basis functions of type (n,l)
                          hubbard_u(i_orb) = plus_u_value(i_species,i_shell) / hartree 

                          ! sets the target value for DFT+U ramping
                          if(plus_u_ramping_defined .eqv. .true.) then
                             hubbard_u_target(i_orb) = hubbard_u(i_orb)
                             hubbard_u_increment(i_orb) = plus_u_ramping_increment(i_species) / hartree
                             ! setting the value to 0 for U ramping
                             hubbard_u(i_orb) = 0.0d0
                          endif


                       endif
                    else if (basis_fn(i_basis) /= prev_fn) then
                       write(use_unit,*)'* ERROR: DFT+U found more then one atomic shell on atom'
                       stop
                    endif
                  prev_center = Cbasis_to_center(i_basis)
                  prev_shell = i_shell
                  prev_fn = basis_fn(i_basis)
               enddo
            else
               if(basisfn_type(i_fn) /= 'atomic' .and.  basisfn_type(i_fn) /= 'hydro') cycle
               i_species = species(basis_atom(i_basis))
                  do i_shell = 1, n_max_shells_plus_u
                    if(basisfn_n(i_fn) /= plus_u_n(i_species,i_shell)) cycle
                    if(basis_l(i_basis) /= plus_u_l(i_species,i_shell)) cycle 

                    if(Cbasis_to_center(i_basis) /= prev_center.or.i_shell /= prev_shell) then 
                       i_orb = i_orb + 1
                       if(i_pass==2)then
                          plusUorb_index(i_orb) = i_basis! point to m = -l entry of basis functions of type (n,l)
                          hubbard_u(i_orb) = plus_u_value(i_species,i_shell) / hartree 

                          ! sets the target value for DFT+U ramping
                          if(plus_u_ramping_defined .eqv. .true.) then
                             hubbard_u_target(i_orb) = hubbard_u(i_orb)
                             hubbard_u_increment(i_orb) = plus_u_ramping_increment(i_species) / hartree
                             ! setting the value to 0 for U ramping
                             hubbard_u(i_orb) =0.0d0
                          endif

                          ! now we have the first basis funktion of the shell, now check how many basis functions
                          ! are in the localized sub space

                          counter = 0

                          do i_basis_2 = i_basis, (atom2basis_off(basis_atom(i_basis)) + sp2n_basis_sp(i_species))
                             i_fn_2 = basis_fn(i_basis_2)
                             if(basisfn_type(i_fn_2) /= 'atomic' .and. basisfn_type(i_fn_2) /= 'hydro') cycle
                             if(basisfn_n(i_fn_2) /= plus_u_n(i_species,i_shell)) cycle
                             if(basis_l(i_basis_2) /= plus_u_l(i_species,i_shell)) cycle

                             counter = counter + 1

                          enddo
                          num_l_shell(i_orb) = counter/(2*basis_l(i_basis)+1)
                       endif
                    ! in principle this should never happen
                    !else if (basis_fn(i_basis) /= prev_fn) then
                    !   write(use_unit,*)'* ERROR: DFT+U found more then one atomic shell on atom'
                    !   stop
                    endif
                  prev_center = Cbasis_to_center(i_basis)
                  prev_shell = i_shell
                  prev_fn = basis_fn(i_basis)
                  enddo
            endif ! hydros defined
         enddo
         
         if(i_pass==2)exit
         n_plusUorb = i_orb
         allocate(plusUorb_index(n_plusUorb))
         allocate(hubbard_u(n_plusUorb))
         allocate(num_l_shell(n_plusUorb))
         allocate(p_max(n_plusUorb))
         if(plus_u_ramping_defined .eqv. .true.) then
           allocate(hubbard_u_target(n_plusUorb))
           allocate(hubbard_u_increment(n_plusUorb))
           ramping_initialized = .true.
         endif

      enddo

      ! if not plus_u_use_hydros is specified there is only one possible set of (n,l) reference basis fkt.
      if(plus_u_hydros_defined .eqv. .false.) then
         num_l_shell(:) = 1
      endif  

      p_max(:) = num_l_shell(:)

      l_max = maxval(plus_u_l(:,:))
      num_shell_max = maxval(num_l_shell(:))

      allocate(occ_plus_u(n_plusUorb, 1:(2*l_max+1), 1:(2*l_max+1), n_spin))

      occ_numbers_initialized = .false.

      !allocate linear coefficents for plusU orbital description
      allocate(hubc(n_plusUorb,1:(2*l_max+1),1:num_shell_max))

      hubc(:,:,:) = 0.0d0

      if(plus_u_hydros_defined .eqv. .false.) then
      ! we only have the atomic basis functions as projectors
      ! there is no linear combination of basis functions
      ! hence, the hubbard coefficents are set to 1

        hubc(:,:,:) = 1.0d0

      else

        do i_orb = 1, n_plusUorb
           i_species = species(basis_atom(plusUorb_index(i_orb)))
           do k = 1, num_shell_max
              hubc(i_orb,1:(2*l_max+1),k) = plus_u_hubc(i_species,k)! * 1.0d0

              ! In case we do not have explicit hubbard coefficients
              ! plus_u_hubc contains only 1.0s. In that case every projector
              ! contributes equally to the correlated subspace

              if (plus_u_hubc(i_species,1) .eq. 1.0d0 .and. plus_u_hubc(i_species,2) .eq. 1.0d0) then
                  ! setting up correct normalization
                  ! instead of a primitive normalization we should do here some fitting routine to get
                  ! the expansion coefficents
                  hubc(i_orb,1:(2*l_max+1),k) = hubc(i_orb,1:(2*l_max+1),k)/sqrt(1.0d0*num_shell_max)
                  !if (myid .eq. 0) then
                  !   write(use_unit,'(1X,A)') &
                  !     '* WARNING: Hubbard projector Normalization is switched off, look at plus_u.f90'
                  !endif


              endif

           enddo
        enddo

      endif


   end subroutine allocate_plus_u

   subroutine plus_u_ramping_work(etot,etot_prev)

   !  PURPOSE
   !
   !  Performes the DFT+U ramping. Useful for systems with bad convergence.
   !
   !  We start with a U value of 0. If the convergence criteria (defined in
   !  plus_u_ramping_accuracy) is fullfilled, the U value will be increased by 
   !  <plus_u_ramping_increment> until the target U value is reached.
   !  This is only done once. This means plus_u_ramping will be switched off
   !  if the target U value is reached.
   !  DFT+U ramping is only active for the initial wavefunction. All following
   !  scf-restarts are unaffected!

   !  USES
      use constants, only: hartree
      use mpi_tasks, only: aims_stop_coll

      implicit none

   !  ARGUMENTS

   !  INPUT

   !  OUTPUT

   !  AUTHOR
   !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
   !  HISTORY
   !    Release version, FHI-aims (2017).

   !  SOURCE

   ! local variables

      real*8, intent(in) :: etot
      real*8, intent(in) :: etot_prev
      integer            :: i_orb

      character(*), parameter :: func = 'plus_u'


      if(plus_u_ramping_defined .eqv. .false.) then
         return
      endif

      ! set initial hubbard_u to 0
      ! Do this only once
      if(ramping_initialized .eqv. .true.) then
         !This doesn't hurt
         hubbard_u = 0.0d0
         ! switching off
         ramping_initialized = .false.
         ramping_counter = 0
         do i_orb = 1, n_plusUorb
            if(hubbard_u_increment(i_orb) .eq. 0.0d0) then
               call aims_stop_coll('Some of the ramping values are zero &
                       & , please define a ramping value for each species.', func)
            endif
         enddo
      endif

      ! increase U value if target accuracy was reached, only.
      if((abs(etot-etot_prev))*hartree .le. plus_u_ramping_accuracy) then

         do i_orb = 1, n_plusUorb
            hubbard_u(i_orb) = hubbard_u(i_orb) + hubbard_u_increment(i_orb)
            if(hubbard_u(i_orb) .gt. hubbard_u_target(i_orb)) then
              hubbard_u(i_orb) = hubbard_u_target(i_orb)
              ramping_counter = 0
            endif
         end do

         ! Final check if all target values have been reached.
         ! If reached switch off plus_u_ramping.
         ! This means plus_u_ramping is only active for the
         ! initial wavefunctions.

         if (ramping_counter .eq. n_plusUorb) then
             plus_u_ramping_defined = .false.
         endif

      else

         return

      endif

   end subroutine plus_u_ramping_work


   subroutine plus_u_init_idx

   !  PURPOSE
   !
   !  Precomputes all sparse indecies for later routines.
   !  This is only needed for mulliken charges 
   !  if periodic boundery condistions are applied.
   !  For cluster calculations pre computer indecies are not 
   !  needed.

   !  USES

      use dimensions
      use runtime_choices
      use basis
      use pbc_lists
      use species_data
      use constants
      use mpi_tasks, only: myid, check_allocation
      use localorb_io, only: use_unit

      implicit none

   !  ARGUMENTS

   !  INPUT

   !  OUTPUT

   !  AUTHOR
   !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
   !  HISTORY
   !    Release version, FHI-aims (2017).

   !  SOURCE

   ! local variables

      integer :: i_orb,i_pass,pk,pl
      integer :: i,j,mi,mj,l_orb,i_cell,k,rtmp,ctmp,idx
      integer :: ialloc

      deallocate(n_idx,stat=ialloc)! counts number of spar indecies
      deallocate(cidx,stat=ialloc) ! stores sparse indecies (i,k)
      deallocate(ridx,stat=ialloc) ! stores sparse indecies (k,j)


      if(packed_matrix_format==PM_index)then

         if (use_local_index) then
            write(use_unit,*)'* ERROR: DFT+U not yet implemented with use_local_index'
            write(use_unit,*)'* See implementation notes in plus_u.f90'
            stop
         endif

         allocate(n_idx(n_plusUorb, 1:(2*l_max+1),1:num_shell_max, 1:(2*l_max+1),1:num_shell_max),stat=ialloc)
         call check_allocation(ialloc,'n_idx','allocate_plus_u')

         do i_pass = 1,2
            ! pass 1: determine the dimension of all arrays
            ! allocate
            ! pass 2: fill allocated arrays

            n_idx(:,:,:,:,:) = 0

            ! In principle this is pretty easy. We performe all looping in the initialization of the scf
            ! and store the computed sparse indecies in ridx and cidx.

            ! Both, Hamiltonian and occupation matrix can be calculated via matrix multiplication.
            ! In both cases the loops are identical.

            do i_orb = 1, n_plusUorb
               l_orb = basis_l(plusUorb_index(i_orb))
               ! In general there are more than 1 basis fkts which correspond to one m,
               ! therefore, also loop over pk and pl.
               do pk = 1, num_l_shell(i_orb)
                  do pl = 1, num_l_shell(i_orb)
                     do mj = 1,(2*l_orb+1)
                        j = plusUorb_index(i_orb) + (mj-1) + (pl-1)*(2*l_orb+1)
                        do mi = 1,(2*l_orb+1)
                           i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                           ! for mulliken charges we consider explicitly the overlap to neighbouring cells
                           do i_cell = 1, n_cells_in_hamiltonian - 1
                              do k = 1, n_basis
                                 ! somewhat inefficient, but simpler to code
                                 rtmp = get_idx_sparse(i_cell,i,k)
                                 ctmp = get_idx_sparse(i_cell,k,j)
                                 if(rtmp/=0.and.ctmp/=0)then
                                    n_idx(i_orb,mi,pk,mj,pl) = n_idx(i_orb,mi,pk,mj,pl)+1
                                    if(i_pass==2)then
                                       ridx(i_orb,mi,pk,mj,pl,n_idx(i_orb,mi,pk,mj,pl)) = rtmp
                                       cidx(i_orb,mi,pk,mj,pl,n_idx(i_orb,mi,pk,mj,pl)) = ctmp
                                    endif
                                 endif
                              enddo ! k
                           enddo ! i_cell
                        enddo ! mi
                     enddo ! mj
                  enddo
              enddo
            enddo ! i_orb
            if(i_pass==2)exit
            allocate(ridx(n_plusUorb, 1:(2*l_max+1),1:num_shell_max, 1:(2*l_max+1),1:num_shell_max, maxval(n_idx(:,:,:,:,:))),stat=ialloc)
            call check_allocation(ialloc,'ridx','allocate_plus_u')
            allocate(cidx(n_plusUorb, 1:(2*l_max+1),1:num_shell_max, 1:(2*l_max+1),1:num_shell_max, maxval(n_idx(:,:,:,:,:))),stat=ialloc)
            call check_allocation(ialloc,'cidx','allocate_plus_u')
         enddo ! i_pass
      endif
     
   end subroutine plus_u_init_idx


   subroutine deallocate_plus_u

   !  USES

   !  ARGUMENTS

   !  INPUT

   !  OUTPUT

   !  AUTHOR
   !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
   !  HISTORY
   !    Release version, FHI-aims (2011).

   !  SOURCE
      integer :: ialloc

      deallocate(plusUorb_index,stat=ialloc)
      deallocate(occ_plus_u,stat=ialloc)
      deallocate(hubbard_u,stat=ialloc)
      deallocate(hubbard_u_target,stat=ialloc)
      deallocate(hubbard_u_increment,stat=ialloc)
      deallocate(AMF_FLL_array,stat=ialloc)
      deallocate(av_occ_array,stat=ialloc)

      deallocate(num_l_shell,stat=ialloc)
      deallocate(p_max,stat=ialloc)
      deallocate(hubc,stat=ialloc)

      if(allocated(ridx)) deallocate(ridx,stat=ialloc)
      if(allocated(cidx)) deallocate(cidx,stat=ialloc)
      if(allocated(n_idx)) deallocate(n_idx,stat=ialloc)

   end subroutine deallocate_plus_u

!-----------------------------------------------------------------------------------------
!------------------------------DFT+U-OCCUPATION-MATRIX------------------------------------
!-----------------------------------------------------------------------------------------


   subroutine occ_numbers_plus_u(density_matrix_sparse, density_matrix, i_spin)

   !  PURPOSE
   !
   !  Calculates the occupation matrix.
   !
   !  on-site: n_mm' = D_mm' (with D being the density matrix)

   !  USES

      use dimensions
      use runtime_choices
      use mpi_tasks
      use physics
      use pbc_lists
      use basis
      use localorb_io

      implicit none

   !  ARGUMENTS

   !  INPUT
   
   !  OUTPUT

   !  AUTHOR
   !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
   !  HISTORY
   !    Release version, FHI-aims (2017).
   !  SOURCE

      real*8,   intent(in), dimension(n_hamiltonian_matrix_size) :: density_matrix_sparse
      real*8,   intent(inout), dimension(:,:)                       :: density_matrix
      integer,  intent(in)                                       :: i_spin
      character*20                                               :: fmt

   !  Local variables

      integer        :: i_cell
      integer        :: pk, pl
      integer        :: i, j, mi, mj, i_orb, l_orb
      integer        :: ialloc, index_sparse
      
      if (.not. use_plus_u) return

      ! BEGIN WORK

      ! If plus_u_matrix_release eq. true, the occupation matrix is calculated in the 
      ! first scf iteration  and is kept fixed until the specified convergence
      ! criteria (plus_u_matrix_relase) is reached.
      if (plus_u_matrix_release_defined .and. occ_numbers_initialized) return

      occ_plus_u(:,:,:,i_spin) = 0.d0

      select case (packed_matrix_format)

      ! For the use of on-site charges only, the occupation matrix is just the corresponding 
      ! density matrix block.
      case (PM_none)

            do i_orb = 1, n_plusUorb
               l_orb = basis_l(plusUorb_index(i_orb))
               do pk = 1, num_l_shell(i_orb)
                  do pl = 1, num_l_shell(i_orb)
                     do mj = 1,(2*l_orb+1)
                        j = plusUorb_index(i_orb) + (mj-1) + (pl-1)*(2*l_orb+1) 
                        do mi = 1,(2*l_orb+1)
                           i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)

                           occ_plus_u(i_orb, mi, mj, i_spin) = occ_plus_u(i_orb, mi, mj, i_spin) + &
                                                                hubc(i_orb,mi,pk) * hubc(i_orb,mj,pl)  * density_matrix(i,j)

                        enddo
                     enddo               
                  enddo
               enddo
            enddo


      case (PM_index)

            do i_orb = 1, n_plusUorb
               l_orb = basis_l(plusUorb_index(i_orb))
               do pk = 1, num_l_shell(i_orb)
                  do pl = 1, num_l_shell(i_orb)
                     do mi = 1,(2*l_orb+1)
                        i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                        do mj = 1,(2*l_orb+1)
                           j = plusUorb_index(i_orb) + (mj-1) + (pl-1)*(2*l_orb+1)

                           index_sparse = get_idx_sparse(1,i,j)


                           if (index_sparse /= 0) then
                              occ_plus_u(i_orb, mi, mj, i_spin) = occ_plus_u(i_orb, mi, mj, i_spin) + hubc(i_orb,mi,pk) &
                                                                * hubc(i_orb,mj,pl) * &
                                                                  density_matrix_sparse(index_sparse)
                           endif

                        enddo
                     enddo
                  enddo     
               enddo
            enddo! i_orb      

      end select

      ! At the initialization of the scf, there is no density-matrix.
      if (n_spin .eq. 2) then
         if (i_spin .eq. 2) then
            occ_numbers_initialized = .true.
         endif
      else
         occ_numbers_initialized = .true.
      endif
 
   end subroutine occ_numbers_plus_u


   subroutine occ_numbers_plus_u_mulliken(density_matrix_sparse, density_matrix, i_spin)

   !  PURPOSE
   !
   !  Calculates the occupation matrix using mulliken charges.
   !
   !  dual: n_mm' = 0.5*(S_mi*D_im' + D_mj*S_jm') (with D being the density matrix)

   !  USES

      use dimensions
      use runtime_choices
      use mpi_tasks
      use physics
      use pbc_lists
      use basis
      use localorb_io

      implicit none

   !  ARGUMENTS

   !  INPUT

   !  OUTPUT

   !  AUTHOR
   !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
   !  HISTORY
   !    Release version, FHI-aims (2017).
   !  SOURCE

      real*8,   intent(in), dimension(n_hamiltonian_matrix_size) :: density_matrix_sparse
      real*8,   intent(inout), dimension(:,:)                    :: density_matrix
      integer,  intent(in)                                       :: i_spin
      character*20                                               :: fmt

   !  Local variables

      integer        :: i_cell
      integer        :: pk, pl, k
      integer        :: i, j, mi, mj, i_orb, l_orb
      integer        :: ialloc

      if (.not. use_plus_u) return

      ! BEGIN WORK

      if (plus_u_matrix_release_defined .and. occ_numbers_initialized) return

      occ_plus_u(:,:,:,i_spin) = 0.d0

      select case (packed_matrix_format)

      case (PM_none)

            do i_orb = 1, n_plusUorb
               l_orb = basis_l(plusUorb_index(i_orb))
               do pk = 1, num_l_shell(i_orb)
                  do pl = 1, num_l_shell(i_orb)
                     do mj = 1,(2*l_orb+1)
                        j = plusUorb_index(i_orb) + (mj-1) + (pl-1)*(2*l_orb+1) 
                        do mi = 1,(2*l_orb+1)
                           i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                           do k = 1, n_centers_basis_T

                              occ_plus_u(i_orb, mi, mj, i_spin) = occ_plus_u(i_orb, mi, mj, i_spin) + &
                                                                0.5 * hubc(i_orb,mi,pk) * hubc(i_orb,mj,pl) &
                                                                    * (density_matrix(i,k) * overlap_matrix_w(k,j) &
                                                                    +  overlap_matrix_w(i,k) * density_matrix(k,j))

                           enddo
                        enddo
                     enddo               
                  enddo
               enddo
            enddo

      case (PM_index)

            do i_orb = 1, n_plusUorb
               l_orb = basis_l(plusUorb_index(i_orb))
               do pk = 1, num_l_shell(i_orb)
                  do pl = 1, num_l_shell(i_orb)
                     do mj = 1,(2*l_orb+1)
                        j = plusUorb_index(i_orb) + (mj-1) + (pl-1)*(2*l_orb+1) 
                        do mi = 1,(2*l_orb+1)
                           i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                           do k = 1, n_idx(i_orb,mi,pk,mj,pl)

                              occ_plus_u(i_orb, mi, mj, i_spin) = occ_plus_u(i_orb, mi, mj, i_spin) + &
                                                                0.5 * hubc(i_orb,mi,pk) * hubc(i_orb,mj,pl) &
                                                                    * (density_matrix_sparse(ridx(i_orb,mi,pk,mj,pl,k)) &
                                                                    *  overlap_matrix(cidx(i_orb,mi,pk,mj,pl,k)) &
                                                                    +  overlap_matrix(ridx(i_orb,mi,pk,mj,pl,k)) &
                                                                    *  density_matrix_sparse(cidx(i_orb,mi,pk,mj,pl,k)))

                           enddo
                        enddo
                     enddo               
                  enddo
               enddo
            enddo

      end select

      ! At the initialization of the scf, there is no density-matrix.
      if (n_spin .eq. 2) then
         if (i_spin .eq. 2) then
            occ_numbers_initialized = .true.
         endif
      else
         occ_numbers_initialized = .true.
      endif
 

   end subroutine occ_numbers_plus_u_mulliken


   subroutine occ_numbers_plus_u_full(density_matrix_sparse, density_matrix, i_spin)

   !  PURPOSE
   !
   !  Calculates the occupation matrix using the full representation.
   !  
   !  !!!WARNING!!!: this is only for testing. Matrix indecies have to be precalculated
   !                 as it is done for the mulliken case. However, for cluster
   !                 calculations no rewrite is needed.
   !
   !  full: n_mm' = S_mi*D_ij*S_jm' (with D being the density matrix)

   !  USES

      use dimensions
      use runtime_choices
      use mpi_tasks
      use physics
      use pbc_lists
      use basis
      use localorb_io

      implicit none

   !  ARGUMENTS

   !  INPUT

   !  OUTPUT

   !  AUTHOR
   !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
   !  HISTORY
   !    Release version, FHI-aims (2017).
   !  SOURCE

      real*8,   intent(in), dimension(n_hamiltonian_matrix_size) :: density_matrix_sparse
      real*8,   intent(inout), dimension(:,:)                       :: density_matrix
      integer,  intent(in)                                       :: i_spin
      character*20                                               :: fmt

   !  Local variables

      integer        :: i_cell, index_sparse, index_sparse2, index_sparse3
      integer        :: pk, pl, k, l
      integer        :: i, j, mi, mj, i_orb, l_orb
      integer        :: ialloc

      if (.not. use_plus_u) return

      ! BEGIN WORK

      if (plus_u_matrix_release_defined .and. occ_numbers_initialized) return

      occ_plus_u(:,:,:,i_spin) = 0.d0

      select case (packed_matrix_format)

      case (PM_none)

            do i_orb = 1, n_plusUorb
               l_orb = basis_l(plusUorb_index(i_orb))
               do pk = 1, num_l_shell(i_orb)
                  do pl = 1, num_l_shell(i_orb)
                     do mj = 1,(2*l_orb+1)
                        j = plusUorb_index(i_orb) + (mj-1) + (pl-1)*(2*l_orb+1) 
                        do mi = 1,(2*l_orb+1)
                           i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                           do k = 1, n_centers_basis_T
                              do l = 1, n_centers_basis_T

                                 occ_plus_u(i_orb, mi, mj, i_spin) = occ_plus_u(i_orb, mi, mj, i_spin) &
                                                                   + hubc(i_orb,mi,pk) * hubc(i_orb,mj,pl) &
                                                                   * overlap_matrix_w(i,k) * density_matrix(k,l) &
                                                                   * overlap_matrix_w(l,j)
                                                                    

                              enddo
                           enddo
                        enddo
                     enddo               
                  enddo
               enddo
            enddo

      case (PM_index)

            do i_orb = 1, n_plusUorb
               l_orb = basis_l(plusUorb_index(i_orb))
               do i_cell = 1, n_cells_in_hamiltonian - 1
                  do pk = 1, num_l_shell(i_orb)
                     do pl = 1, num_l_shell(i_orb)
                        do mj = 1,(2*l_orb+1)
                           j = plusUorb_index(i_orb) + (mj-1) + (pl-1)*(2*l_orb+1) 
                           do mi = 1,(2*l_orb+1)
                              i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                              do k = 1, n_basis
                                 do l = 1, n_basis ! density matrix is stored as full matrix

                                    ! This is ugly, I know.
                                    ! However, this is just for testing and only single point calculations are 
                                    ! available for this occupation matrix definition. In the future, the matrix
                                    ! indecies should be precalculated as it is done for the mulliken charge method.

                                    index_sparse  = get_idx_sparse(i_cell,i,k)
                                    index_sparse2 = get_idx_sparse(i_cell,k,l) 
                                    index_sparse3 = get_idx_sparse(i_cell,l,j)


                                    if (index_sparse /= 0 .and. index_sparse2 /= 0 .and. index_sparse3 /= 0) then

                                       occ_plus_u(i_orb, mi, mj, i_spin) = occ_plus_u(i_orb, mi, mj, i_spin) &
                                                                         + hubc(i_orb,mi,pk) * hubc(i_orb,mj,pl) &
                                                                         * overlap_matrix(index_sparse) &
                                                                         * density_matrix_sparse(index_sparse2) &
                                                                         * overlap_matrix(index_sparse3)

                                    endif

                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo               
                  enddo
               enddo
            enddo

      end select

      ! At the initialization of the scf, there is no density-matrix.
      if (n_spin .eq. 2) then
         if (i_spin .eq. 2) then
            occ_numbers_initialized = .true.
         endif
      else
         occ_numbers_initialized = .true.
      endif
 

   end subroutine occ_numbers_plus_u_full

!-----------------------------------------------------------------------------------------
!------------------------ALL-HAMILTONIAN-RELATED-ROUTINES---------------------------------
!-----------------------------------------------------------------------------------------


   subroutine add_plus_u_to_hamiltonian(hamiltonian)

      !  PURPOSE
      !
      !  Adds the DFT+U Hamiltonian correction to the usual DFT Hamiltonian
      !
      !  on-site: h_U = v_mm' = -U ( n_mm' - 0.5)


      !  USES

      use dimensions
      use runtime_choices
      use basis
      use pbc_lists
      use species_data
      use constants
      use mpi_tasks, only: myid, check_allocation
      use localorb_io, only: use_unit, OL_norm, localorb_info

      implicit none

      ! ARGUMENTS

      real*8,intent(inout)  :: hamiltonian( n_hamiltonian_matrix_size, n_spin )

      ! INPUT

      ! OUTPUT

      ! AUTHOR
      !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
      ! HISTORY
      !    Release version, FHI-aims (2017).
      ! SOURCE

      integer :: i_spin, i, j, k, mi, mj, i_orb, l_orb, l, pk, pl
      integer :: ialloc, index_sparse, i_cell
      real*8 :: av_occ(n_spin)
      real*8 :: AMF_FLL_mixing_numer, AMF_FLL_mixing_denom, AMF_FLL_mixing_factor 
      character*20 :: fmt
      character*200  :: info_str

      if(.not. use_plus_u)then
         if (myid .eq. 0) then
            write(use_unit,'(1X,A)') &
            '* ERROR: Internal error in DFT+U - this should never happen!'
         endif
         stop
      endif

      if(use_local_index)then
         write(use_unit,*)'* ERROR: DFT+U not yet implemented with use_local_index'
         write(use_unit,*)'* See implementation notes in plus_u.f90'
         stop
      endif

      if(.not.occ_numbers_initialized)return  ! skip this when occupation numbers are not yet initialized


      if (.not. allocated(AMF_FLL_array)) then ! save values for later routines
         allocate(AMF_FLL_array(n_plusUorb),stat=ialloc)
         call check_allocation(ialloc,'AMF_FLL_array','add_plus_u_to_hamiltonian')
      endif
       
      if (.not. allocated(av_occ_array)) then  ! save values for later routines
         allocate(av_occ_array(n_spin,n_plusUorb),stat=ialloc)
         call check_allocation(ialloc,'av_occ_array','add_plus_u_to_hamiltonian')
      endif

      if(output_level /= 'MD_light') then
            if(myid == 0)then
               write(use_unit,*)" ----------------------------------------------------------------"
               write(use_unit,*)" Adding Hubbard correction to Hamiltonian"
               write(use_unit,*)" ----------------------------------------------------------------" 
            endif
            write(info_str,'(2X,A)') "!--Definition for DFT+U occupation matrix is:    'on-site'" 
            call localorb_info ( info_str,use_unit,'(A)', OL_norm  ) 
      endif
 
      do i_orb = 1, n_plusUorb
         l_orb = basis_l(plusUorb_index(i_orb))

         av_occ(:) = 0.d0

         do i_spin = 1,n_spin
            do mi = 1,(2*l_orb+1)

               av_occ(i_spin) = av_occ(i_spin) + occ_plus_u(i_orb, mi, mi, i_spin)

            enddo
            av_occ(i_spin) = av_occ(i_spin) / (2*l_orb+1)
            av_occ_array(i_spin,i_orb) = av_occ(i_spin)
         enddo

         if(plus_u_petukhov_mixing_defined)then
            AMF_FLL_mixing_factor = plus_u_petukhov_mixing
            AMF_FLL_array(i_orb)  = plus_u_petukhov_mixing
         else
            ! Determine Petukhov mixing factor:
            !   mixing_factor == 0.0  =>  AMF - 'around mean field'
            !   mixing_factor == 1.0  =>  FLL - 'fully localized limit
            !   0.0 < mixing_factor < 1.0 linear interpolation

            AMF_FLL_mixing_numer = 0.d0
            AMF_FLL_mixing_denom = 0.d0
            do i_spin = 1,n_spin
              
               do mi = 1,(2*l_orb+1)
                  do mj = 1,(2*l_orb+1)
                     
                     AMF_FLL_mixing_numer = AMF_FLL_mixing_numer &
                                   + occ_plus_u(i_orb, mi, mj, i_spin) * occ_plus_u(i_orb, mj, mi, i_spin)

                  enddo

                     AMF_FLL_mixing_numer = AMF_FLL_mixing_numer &
                           - 2./p_max(i_orb)*occ_plus_u(i_orb, mi, mi, i_spin) * av_occ(i_spin)

                     AMF_FLL_mixing_numer = AMF_FLL_mixing_numer + av_occ(i_spin) ** 2

               enddo
                  AMF_FLL_mixing_denom = AMF_FLL_mixing_denom &
                     + (2*l_orb+1) * av_occ(i_spin) * (1.d0-av_occ(i_spin))
            enddo
               AMF_FLL_mixing_factor = AMF_FLL_mixing_numer / AMF_FLL_mixing_denom
               AMF_FLL_array(i_orb) = AMF_FLL_mixing_factor
         endif

         if(output_level /= 'MD_light') then

            if(myid .eq. 0)then

               write(use_unit,'(X,A,I2,6X,A,F6.2,X,A)')" &
                   &correlated subspace ",i_orb, '|  U = ',hubbard_u(i_orb)*hartree, 'eV'
               write(use_unit,'(X,A,F8.4,F8.4)')"  avg. occ. hubbard projectors:",av_occ(:)
               write(use_unit,'(X,A,F8.4)')"  petukhov mixing factor      :",     AMF_FLL_mixing_factor
               write(use_unit,*)" ----------------------------------------------------------------"

               do i_spin = 1,n_spin         
                  write(use_unit,*)"  occupation matrix (subspace #",i_orb,", spin ",i_spin,")"
                  write(fmt,'(A,I2,A)')" (",l_orb*2+1,"f10.5)"
                  write(use_unit,fmt)occ_plus_u(i_orb, 1:(2*l_orb+1), 1:(2*l_orb+1), i_spin)
                  write(use_unit,*)" ----------------------------------------------------------------"
               enddo

            endif
         endif

         select case (packed_matrix_format)

         case (PM_none)

            do i_spin = 1,n_spin
               do pk = 1,num_l_shell(i_orb)
                  do pl = 1,num_l_shell(i_orb)
                     do mi = 1,(2*l_orb+1)
                        i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                        do mj = 1,(2*l_orb+1)
                           j = plusUorb_index(i_orb) + (mj-1) + (pl-1)*(2*l_orb+1)
 

                           call hamiltonian_add(i,j,i_spin, &
                              - hubbard_u(i_orb)*hubc(i_orb,mi,pk)*hubc(i_orb,mj,pl)*occ_plus_u(i_orb, mi, mj, i_spin))

                        enddo
                     enddo
                  enddo
               enddo
            enddo
            ! diagonal term
            do i_spin = 1,n_spin
               do pk = 1,num_l_shell(i_orb)
                  do pl = 1,num_l_shell(i_orb)
                     do mi = 1,(2*l_orb+1)
                        i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                        j = plusUorb_index(i_orb) + (mi-1) + (pl-1)*(2*l_orb+1)


                        call hamiltonian_add(i,j,i_spin, &
                           ((1.d0-AMF_FLL_mixing_factor)*av_occ(i_spin) + AMF_FLL_mixing_factor*0.5) &
                           * hubbard_u(i_orb) * hubc(i_orb,mi,pk) * hubc(i_orb,mi,pl))

                     enddo
                  enddo
               enddo
            enddo

         ! periodic systems
         case (PM_index)
         
            do i_spin = 1,n_spin
               do pk = 1,num_l_shell(i_orb)
                  do pl = 1,num_l_shell(i_orb)
                     do mi = 1,(2*l_orb+1)
                        i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                        do mj = mi,(2*l_orb+1)
                           j = plusUorb_index(i_orb) + (mj-1) + (pl-1)*(2*l_orb+1)

                           index_sparse = get_idx_sparse(1,i,j)

                           ! this doesn't hurt, matrices are small in that case.
                           if (index_sparse /= 0) then
                              hamiltonian(index_sparse,i_spin) = hamiltonian(index_sparse,i_spin) &
                                    - hubbard_u(i_orb) &
                                    * occ_plus_u(i_orb, mi, mj, i_spin) * hubc(i_orb,mi,pk) * hubc(i_orb,mj,pl)
                           endif

                        enddo
                     enddo
                  enddo
               enddo
             enddo !i_spin

             ! diagonal term
             do i_spin = 1,n_spin
               do pk = 1,num_l_shell(i_orb)
                  do pl = 1,num_l_shell(i_orb)
                     do mi = 1,(2*l_orb+1)
                        i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                        j = plusUorb_index(i_orb) + (mi-1) + (pl-1)*(2*l_orb+1)

                        index_sparse = get_idx_sparse(1,i,j)

                        if (index_sparse /= 0) then
                           hamiltonian(index_sparse,i_spin) = hamiltonian(index_sparse,i_spin) &
                              + ((1.d0-AMF_FLL_mixing_factor)*av_occ(i_spin) + AMF_FLL_mixing_factor*0.5) &
                              * hubbard_u(i_orb) * hubc(i_orb,mi,pk) * hubc(i_orb,mi,pl)
                        endif

                     enddo
                  enddo
               enddo
            enddo
            
         endselect

      enddo !i_orb


   contains

      subroutine hamiltonian_add(i,j,s,value)
         ! this has to be an inner subroutine of add_plus_u_to_hamiltonian
         ! so that it can act on 'hamiltonian' argument of the outer routine
         integer,intent(in) :: i,j,s
         real*8,intent(in) :: value
         if(i<=j)then
            hamiltonian(i+j*(j-1)/2,s) = hamiltonian(i+j*(j-1)/2,s) + value
         endif
      end subroutine hamiltonian_add

   end subroutine add_plus_u_to_hamiltonian


   subroutine add_plus_u_mulliken_to_hamiltonian(hamiltonian)

      !  PURPOSE
      !
      !  Adds the DFT+U Hamiltonian correction to the usual DFT Hamiltonian
      !  This Hamiltonian correction belongs to mulliken charges method.
      !
      !  dual: h_U = 0.5*(v*S + S*v)
      !  S: overlap matrix
      !  v: effective_potential v


      !  USES

      use dimensions
      use runtime_choices
      use basis
      use pbc_lists
      use species_data
      use constants
      use mpi_tasks, only: myid, check_allocation
      use localorb_io, only: use_unit, OL_norm, localorb_info

      implicit none

      ! ARGUMENTS

      real*8,intent(inout)  :: hamiltonian( n_hamiltonian_matrix_size, n_spin )

      ! INPUT

      ! OUTPUT

      ! AUTHOR
      !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
      ! HISTORY
      !    Release version, FHI-aims (2017).
      ! SOURCE

      integer :: i_spin, i, j, k, mi, mj, i_orb, l_orb, l, pk, pl
      integer :: ialloc, index_sparse, i_cell
      real*8 :: av_occ(n_spin)
      real*8 :: AMF_FLL_mixing_numer, AMF_FLL_mixing_denom, AMF_FLL_mixing_factor 
      character*20 :: fmt
      character*200  :: info_str

      if(.not. use_plus_u)then
         if (myid .eq. 0) then
            write(use_unit,'(1X,A)') &
            '* ERROR: Internal error in DFT+U - this should never happen!'
         endif
         stop
      endif

      if(use_local_index)then
         write(use_unit,*)'* ERROR: DFT+U not yet implemented with use_local_index'
         write(use_unit,*)'* See implementation notes in plus_u.f90'
         stop
      endif

      if(.not.occ_numbers_initialized)return  ! skip this when occupation numbers are not yet initialized


      if (.not. allocated(AMF_FLL_array)) then ! save values for later routines
         allocate(AMF_FLL_array(n_plusUorb),stat=ialloc)
         call check_allocation(ialloc,'AMF_FLL_array','add_plus_u_to_hamiltonian')
      endif
       
      if (.not. allocated(av_occ_array)) then  ! save values for later routines
         allocate(av_occ_array(n_spin,n_plusUorb),stat=ialloc)
         call check_allocation(ialloc,'av_occ_array','add_plus_u_to_hamiltonian')
      endif

      if(output_level /= 'MD_light') then
            if(myid == 0)then
               write(use_unit,*)" ----------------------------------------------------------------"
               write(use_unit,*)" Adding Hubbard correction to Hamiltonian"
               write(use_unit,*)" ----------------------------------------------------------------" 
            endif
            write(info_str,'(2X,A)') "!--Definition for DFT+U occupation matrix is:    'dual'" 
            call localorb_info ( info_str,use_unit,'(A)', OL_norm  ) 
      endif
 
      do i_orb = 1, n_plusUorb
         l_orb = basis_l(plusUorb_index(i_orb))

         av_occ(:) = 0.d0

         do i_spin = 1,n_spin
            do mi = 1,(2*l_orb+1)

               av_occ(i_spin) = av_occ(i_spin) + occ_plus_u(i_orb, mi, mi, i_spin)

            enddo
            av_occ(i_spin) = av_occ(i_spin) / (2*l_orb+1)
            av_occ_array(i_spin,i_orb) = av_occ(i_spin)
         enddo

         if(plus_u_petukhov_mixing_defined)then
            AMF_FLL_mixing_factor = plus_u_petukhov_mixing
            AMF_FLL_array(i_orb)  = plus_u_petukhov_mixing
         else
            ! Determine Petukhov mixing factor:
            !   mixing_factor == 0.0  =>  AMF - 'around mean field'
            !   mixing_factor == 1.0  =>  FLL - 'fully localized limit
            !   0.0 < mixing_factor < 1.0 linear interpolation

            AMF_FLL_mixing_numer = 0.d0
            AMF_FLL_mixing_denom = 0.d0
            do i_spin = 1,n_spin
              
               do mi = 1,(2*l_orb+1)
                  do mj = 1,(2*l_orb+1)
                     
                     AMF_FLL_mixing_numer = AMF_FLL_mixing_numer &
                                   + occ_plus_u(i_orb, mi, mj, i_spin) * occ_plus_u(i_orb, mj, mi, i_spin)

                  enddo

                     AMF_FLL_mixing_numer = AMF_FLL_mixing_numer &
                           - 2./p_max(i_orb)*occ_plus_u(i_orb, mi, mi, i_spin) * av_occ(i_spin)

                     AMF_FLL_mixing_numer = AMF_FLL_mixing_numer + av_occ(i_spin) ** 2

               enddo
                  AMF_FLL_mixing_denom = AMF_FLL_mixing_denom &
                     + (2*l_orb+1) * av_occ(i_spin) * (1.d0-av_occ(i_spin))
            enddo
               AMF_FLL_mixing_factor = AMF_FLL_mixing_numer / AMF_FLL_mixing_denom
               AMF_FLL_array(i_orb) = AMF_FLL_mixing_factor
         endif

         if(output_level /= 'MD_light') then

            if(myid .eq. 0)then
               write(use_unit,'(X,A,I2,6X,A,F6.2,X,A)')" &
                   & correlated subspace ",i_orb, '|  U = ',hubbard_u(i_orb)*hartree, 'eV'
               write(use_unit,'(X,A,F8.4,F8.4)')"  avg. occ. hubbard projectors:",av_occ(:)
               write(use_unit,'(X,A,F8.4)')"  petukhov mixing factor      :",     AMF_FLL_mixing_factor
               write(use_unit,*)" ----------------------------------------------------------------"

               do i_spin = 1,n_spin         
                  write(use_unit,*)"  occupation matrix (subspace #",i_orb,", spin ",i_spin,")"
                  write(fmt,'(A,I2,A)')" (",l_orb*2+1,"f10.5)"
                  write(use_unit,fmt)occ_plus_u(i_orb, 1:(2*l_orb+1), 1:(2*l_orb+1), i_spin)
                  write(use_unit,*)" ----------------------------------------------------------------"
               enddo

            endif
         endif

         select case (packed_matrix_format)

         case (PM_none)

            do i_spin = 1,n_spin
               do pk = 1,num_l_shell(i_orb)
                  do pl = 1,num_l_shell(i_orb)
                     do mi = 1,(2*l_orb+1)
                        i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                        do k = 1, n_basis
                           do mj = 1,(2*l_orb+1)
                              j = plusUorb_index(i_orb) + (mj-1) + (pl-1)*(2*l_orb+1)
 
                              call hamiltonian_add(i,k,i_spin, &
                               - 0.5 * hubbard_u(i_orb) &
                                     * hubc(i_orb,mi,pk)*hubc(i_orb,mj,pl) &
                                     * occ_plus_u(i_orb, mi, mj, i_spin) * overlap_matrix_w(j,k))

                              call hamiltonian_add(k,i,i_spin, &
                               - 0.5 * hubbard_u(i_orb) &
                                     * hubc(i_orb,mi,pk) * hubc(i_orb,mj,pl) &
                                     * overlap_matrix_w(k,j) * occ_plus_u(i_orb, mj, mi, i_spin))

                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            ! diagonal term
            do i_spin = 1,n_spin
               do pk = 1,num_l_shell(i_orb)
                  do pl = 1,num_l_shell(i_orb)
                     do mi = 1,(2*l_orb+1)
                        i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                        j = plusUorb_index(i_orb) + (mi-1) + (pl-1)*(2*l_orb+1)
                        do k = 1, n_basis

                           call hamiltonian_add(i,k,i_spin, &
                              0.5 * ((1.d0-AMF_FLL_mixing_factor)*av_occ(i_spin) + AMF_FLL_mixing_factor*0.5) &
                                  * hubbard_u(i_orb) &
                                  * overlap_matrix_w(i,k) * hubc(i_orb,mi,pk) * hubc(i_orb,mi,pl) &
                                )

                           call hamiltonian_add(k,i,i_spin, &
                              0.5 * ((1.d0-AMF_FLL_mixing_factor)*av_occ(i_spin) + AMF_FLL_mixing_factor*0.5) &
                                  * hubbard_u(i_orb) &
                                  * overlap_matrix_w(k,i) * hubc(i_orb,mi,pk) * hubc(i_orb,mi,pl) &
                                )
                       
                        enddo
                     enddo
                  enddo
               enddo
            enddo

         ! periodic systems
         case (PM_index)
         
            do i_spin = 1,n_spin
               do pk = 1,num_l_shell(i_orb)
                  do pl = 1,num_l_shell(i_orb)
                     do mi = 1,(2*l_orb+1)
                        i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                        do mj = mi,(2*l_orb+1)
                           j = plusUorb_index(i_orb) + (mj-1) + (pl-1)*(2*l_orb+1)
                           do k = 1, n_idx(i_orb,mi,pk,mj,pl)

                              ! ok, this looks strange.
                              ! However, we are using precomputed indecies.
                              ! See, plus_u_init_idx for explaination.

                              hamiltonian(ridx(i_orb,mi,pk,mj,pl,k),i_spin) = hamiltonian(ridx(i_orb,mi,pk,mj,pl,k),i_spin) &
                                                                   - 0.5 * hubbard_u(i_orb) &
                                                                   * occ_plus_u(i_orb, mi, mj, i_spin) &
                                                                   * overlap_matrix(ridx(i_orb,mj,pl,mi,pk,k)) &
                                                                   * hubc(i_orb,mi,pk) * hubc(i_orb,mj,pl)

                              hamiltonian(cidx(i_orb,mj,pl,mi,pk,k),i_spin) = hamiltonian(cidx(i_orb,mj,pl,mi,pk,k),i_spin) &
                                                                   - 0.5 * hubbard_u(i_orb) &
                                                                   * overlap_matrix(cidx(i_orb,mi,pk,mj,pl,k)) &
                                                                   * occ_plus_u(i_orb, mj, mi, i_spin) &
                                                                   * hubc(i_orb,mi,pk) * hubc(i_orb,mj,pl)

                           enddo
                        enddo
                     enddo
                  enddo
               enddo
             enddo !i_spin

             ! diagonal term
             do i_spin = 1,n_spin
               do pk = 1,num_l_shell(i_orb)
                  do pl = 1,num_l_shell(i_orb)
                     do mi = 1,(2*l_orb+1)
                        i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                        j = plusUorb_index(i_orb) + (mi-1) + (pl-1)*(2*l_orb+1)
                        do k = 1, n_idx(i_orb,mi,pk,mi,pl)

                           hamiltonian(ridx(i_orb,mi,pk,mi,pl,k),i_spin) = hamiltonian(ridx(i_orb,mi,pk,mi,pl,k),i_spin) &
                              + 0.5 * ((1.d0-AMF_FLL_mixing_factor)*av_occ(i_spin) + AMF_FLL_mixing_factor*0.5) &
                              * hubbard_u(i_orb) &
                              * overlap_matrix(ridx(i_orb,mi,pk,mi,pl,k)) &
                              * hubc(i_orb,mi,pk) * hubc(i_orb,mi,pl)

                           hamiltonian(cidx(i_orb,mi,pk,mi,pl,k),i_spin) = hamiltonian(cidx(i_orb,mi,pk,mi,pl,k),i_spin) &
                              + 0.5 * ((1.d0-AMF_FLL_mixing_factor)*av_occ(i_spin) + AMF_FLL_mixing_factor*0.5) &
                              * hubbard_u(i_orb) &
                              * overlap_matrix(cidx(i_orb,mi,pk,mi,pl,k)) &
                              * hubc(i_orb,mi,pk) * hubc(i_orb,mi,pl)

                        enddo
                     enddo
                  enddo
               enddo
            enddo
            
         endselect

      enddo !i_orb


   contains

      subroutine hamiltonian_add(i,j,s,value)
         ! this has to be an inner subroutine of add_plus_u_to_hamiltonian
         ! so that it can act on 'hamiltonian' argument of the outer routine
         integer,intent(in) :: i,j,s
         real*8,intent(in) :: value
         if(i<=j)then
            hamiltonian(i+j*(j-1)/2,s) = hamiltonian(i+j*(j-1)/2,s) + value
         endif
      end subroutine hamiltonian_add

   end subroutine add_plus_u_mulliken_to_hamiltonian

   subroutine add_plus_u_full_to_hamiltonian(hamiltonian)

      !  PURPOSE
      !
      !  Adds the DFT+U Hamiltonian correction to the usual DFT Hamiltonian
      !  This Hamiltonian correction belongs to the full representation of
      !  the DFT+U occupation matrix.
      !
      !  !!!WARNING!!!: this routine is for tests only. The matrix indecies
      !                 have to be precomputed as it is done for the mulliken
      !                 charge DFT+U version. This means, for large periodic
      !                 systems, the calculations can take very long.
      !                 However, for cluster calculations
      !                 no rewrite is needed.
      !
      !  full: h_U = S_im*v_mm'*S_mj'
      !  S: overlap matrix
      !  v: effective_potential v


      !  USES

      use dimensions
      use runtime_choices
      use basis
      use pbc_lists
      use species_data
      use constants
      use mpi_tasks, only: myid,check_allocation
      use localorb_io, only: use_unit, OL_norm, localorb_info

      implicit none

      ! ARGUMENTS

      real*8,intent(inout)  :: hamiltonian( n_hamiltonian_matrix_size, n_spin )

      ! INPUT

      ! OUTPUT

      ! AUTHOR
      !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
      ! HISTORY
      !    Release version, FHI-aims (2017).
      ! SOURCE

      integer :: i_spin, i, j, k, mi, mj, i_orb, l_orb, l, pk, pl
      integer :: ialloc, index_sparse, index_sparse2, index_sparse3, i_cell
      real*8 :: av_occ(n_spin)
      real*8 :: AMF_FLL_mixing_numer, AMF_FLL_mixing_denom, AMF_FLL_mixing_factor 
      character*20 :: fmt
      character*200  :: info_str

      if(.not. use_plus_u)then
         if (myid .eq. 0) then
            write(use_unit,'(1X,A)') &
            '* ERROR: Internal error in DFT+U - this should never happen!'
         endif
         stop
      endif

      if(use_local_index)then
         write(use_unit,*)'* ERROR: DFT+U not yet implemented with use_local_index'
         write(use_unit,*)'* See implementation notes in plus_u.f90'
         stop
      endif

      if(.not.occ_numbers_initialized)return  ! skip this when occupation numbers are not yet initialized


      if (.not. allocated(AMF_FLL_array)) then ! save values for later routines
         allocate(AMF_FLL_array(n_plusUorb),stat=ialloc)
         call check_allocation(ialloc,'AMF_FLL_array','add_plus_u_to_hamiltonian')
      endif
       
      if (.not. allocated(av_occ_array)) then  ! save values for later routines
         allocate(av_occ_array(n_spin,n_plusUorb),stat=ialloc)
         call check_allocation(ialloc,'av_occ_array','add_plus_u_to_hamiltonian')
      endif

      if(output_level /= 'MD_light') then
            if(myid == 0)then
               write(use_unit,*)" ----------------------------------------------------------------"
               write(use_unit,*)" Adding Hubbard correction to Hamiltonian"
               write(use_unit,*)" ----------------------------------------------------------------" 
            endif
            write(info_str,'(2X,A)') "!--Definition for DFT+U occupation matrix is:    'full'" 
            call localorb_info ( info_str,use_unit,'(A)', OL_norm  ) 
      endif
 
      do i_orb = 1, n_plusUorb
         l_orb = basis_l(plusUorb_index(i_orb))

         av_occ(:) = 0.d0

         do i_spin = 1,n_spin
            do mi = 1,(2*l_orb+1)

               av_occ(i_spin) = av_occ(i_spin) + occ_plus_u(i_orb, mi, mi, i_spin)

            enddo
            av_occ(i_spin) = av_occ(i_spin) / (2*l_orb+1)
            av_occ_array(i_spin,i_orb) = av_occ(i_spin)
         enddo

         if(plus_u_petukhov_mixing_defined)then
            AMF_FLL_mixing_factor = plus_u_petukhov_mixing
            AMF_FLL_array(i_orb)  = plus_u_petukhov_mixing
         else
            ! Determine Petukhov mixing factor:
            !   mixing_factor == 0.0  =>  AMF - 'around mean field'
            !   mixing_factor == 1.0  =>  FLL - 'fully localized limit
            !   0.0 < mixing_factor < 1.0 linear interpolation

            AMF_FLL_mixing_numer = 0.d0
            AMF_FLL_mixing_denom = 0.d0
            do i_spin = 1,n_spin
              
               do mi = 1,(2*l_orb+1)
                  do mj = 1,(2*l_orb+1)
                     
                     AMF_FLL_mixing_numer = AMF_FLL_mixing_numer &
                                   + occ_plus_u(i_orb, mi, mj, i_spin) * occ_plus_u(i_orb, mj, mi, i_spin)

                  enddo

                     AMF_FLL_mixing_numer = AMF_FLL_mixing_numer &
                           - 2./p_max(i_orb)*occ_plus_u(i_orb, mi, mi, i_spin) * av_occ(i_spin)

                     AMF_FLL_mixing_numer = AMF_FLL_mixing_numer + av_occ(i_spin) ** 2

               enddo
                  AMF_FLL_mixing_denom = AMF_FLL_mixing_denom &
                     + (2*l_orb+1) * av_occ(i_spin) * (1.d0-av_occ(i_spin))
            enddo
               AMF_FLL_mixing_factor = AMF_FLL_mixing_numer / AMF_FLL_mixing_denom
               AMF_FLL_array(i_orb) = AMF_FLL_mixing_factor
         endif

         if(output_level /= 'MD_light') then

            if(myid .eq. 0)then

               write(use_unit,'(X,A,I2,6X,A,F6.2,X,A)')" & 
                   & correlated subspace ",i_orb, '|  U = ',hubbard_u(i_orb)*hartree, 'eV'
               write(use_unit,'(X,A,F8.4,F8.4)')"  avg. occ. hubbard projectors:",av_occ(:)
               write(use_unit,'(X,A,F8.4)')"  petukhov mixing factor      :",     AMF_FLL_mixing_factor
               write(use_unit,*)" ----------------------------------------------------------------"

               do i_spin = 1, n_spin         
                  write(use_unit,*)"  occupation matrix (subspace #",i_orb,", spin ",i_spin,")"
                  write(fmt,'(A,I2,A)')" (",l_orb*2+1,"f10.5)"
                  write(use_unit,fmt)occ_plus_u(i_orb, 1:(2*l_orb+1), 1:(2*l_orb+1), i_spin)
                  write(use_unit,*)" ----------------------------------------------------------------"
               enddo

            endif
         endif

         select case (packed_matrix_format)

         case (PM_none)

            do i_spin = 1,n_spin
               do pk = 1,num_l_shell(i_orb)
                  do pl = 1,num_l_shell(i_orb)
                     do mi = 1,(2*l_orb+1)
                        i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                        do k = 1, n_basis
                           do l = 1, n_basis
                              do mj = 1,(2*l_orb+1)
                                 j = plusUorb_index(i_orb) + (mj-1) + (pl-1)*(2*l_orb+1)

                                 call hamiltonian_add(k,l,i_spin, &
                                      - hubbard_u(i_orb) &
                                      * hubc(i_orb,mi,pk)*hubc(i_orb,mj,pl) &
                                      * overlap_matrix_w(k,j) &
                                      * occ_plus_u(i_orb, mj, mi, i_spin) * overlap_matrix_w(i,l) )

                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            ! diagonal term
            do i_spin = 1,n_spin
               do pk = 1,num_l_shell(i_orb)
                  do pl = 1,num_l_shell(i_orb)
                     do mi = 1,(2*l_orb+1)
                        i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                        j = plusUorb_index(i_orb) + (mi-1) + (pl-1)*(2*l_orb+1)
                        do k = 1, n_basis
                           do l = 1, n_basis

                                call hamiltonian_add(k,l,i_spin, &
                                     ((1.d0-AMF_FLL_mixing_factor)*av_occ(i_spin) + AMF_FLL_mixing_factor*0.5) &
                                     * hubbard_u(i_orb) &
                                     * overlap_matrix_w(k,j) * overlap_matrix_w(i,l) &
                                     * hubc(i_orb,mi,pk) * hubc(i_orb,mi,pl) )

                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo

         ! periodic systems
         case (PM_index)

             do i_spin = 1, n_spin
               do i_cell = 1, n_cells_in_hamiltonian - 1
                  do pk = 1,num_l_shell(i_orb)
                     do pl = 1,num_l_shell(i_orb)
                        do mi = 1,(2*l_orb+1)
                           i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                           do k = 1, n_basis
                              do l = k, n_basis ! only upper triangular matrix is stored
                                 do mj = 1,(2*l_orb+1)
                                    j = plusUorb_index(i_orb) + (mj-1) + (pl-1)*(2*l_orb+1)

                                    ! This is ugly, I know.
                                    ! However, this is just for testing and only single point calculations are 
                                    ! available for this occupation matrix definition. In the future, the matrix
                                    ! indecies should be precalculated as it is done for the mulliken charge method.


                                    index_sparse  = get_idx_sparse(i_cell,k,l)
                                    index_sparse2 = get_idx_sparse(i_cell,k,j)
                                    index_sparse3 = get_idx_sparse(i_cell,i,l)

                                    if (index_sparse /= 0 .and. index_sparse2 /= 0 .and. index_sparse3 /= 0) then

                                       hamiltonian(index_sparse,i_spin) = hamiltonian(index_sparse,i_spin) &
                                            - hubbard_u(i_orb) &
                                            * hubc(i_orb,mi,pk)*hubc(i_orb,mj,pl) &
                                            * overlap_matrix(index_sparse2) &
                                            * occ_plus_u(i_orb, mj, mi, i_spin) * overlap_matrix(index_sparse3) 

                                    endif

                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            ! diagonal term
            do i_spin = 1, n_spin
               do i_cell = 1, n_cells_in_hamiltonian - 1
                  do pk = 1,num_l_shell(i_orb)
                     do pl = 1,num_l_shell(i_orb)
                        do mi = 1,(2*l_orb+1)
                           i = plusUorb_index(i_orb) + (mi-1) + (pk-1)*(2*l_orb+1)
                           j = plusUorb_index(i_orb) + (mi-1) + (pl-1)*(2*l_orb+1)
                           do k = 1, n_basis
                              do l = k, n_basis ! only upper triangular matrix is stored

                                 index_sparse  = get_idx_sparse(i_cell,k,l)
                                 index_sparse2 = get_idx_sparse(i_cell,k,j)
                                 index_sparse3 = get_idx_sparse(i_cell,i,l)

                                 if (index_sparse /= 0 .and. index_sparse2 /= 0 .and. index_sparse3 /= 0) then

                                    hamiltonian(index_sparse,i_spin) = hamiltonian(index_sparse,i_spin) &
                                           + ((1.d0-AMF_FLL_mixing_factor)*av_occ(i_spin) + AMF_FLL_mixing_factor*0.5) &
                                           * hubbard_u(i_orb) &
                                           * overlap_matrix(index_sparse2) * overlap_matrix(index_sparse3) &
                                           * hubc(i_orb,mi,pk) * hubc(i_orb,mi,pl) 

                                endif

                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo

           
         endselect

      enddo !i_orb


   contains

      subroutine hamiltonian_add(i,j,s,value)
         ! this has to be an inner subroutine of add_plus_u_to_hamiltonian
         ! so that it can act on 'hamiltonian' argument of the outer routine
         integer,intent(in) :: i,j,s
         real*8,intent(in) :: value
         if(i<=j)then
            hamiltonian(i+j*(j-1)/2,s) = hamiltonian(i+j*(j-1)/2,s) + value
         endif
      end subroutine hamiltonian_add

   end subroutine add_plus_u_full_to_hamiltonian


!------------------------------------------------------------------------------------------
!-----------------START OF ENERGY CORRECTION-----------------------------------------------
!------------------------------------------------------------------------------------------

   subroutine plus_u_energy_correction_term

      !  PURPOSE
      !
      !  calculates the DFT+U energy correction according to
      !  E_U - \sum_i f_i < Psi_i | v_U | Psi_i >
      !  The second expression is basically the sum of eigenvalues of the hubbard_operator.
      !  FHI-AIMS uses the Harris functional to calculate the total energy, so we have to 
      !  substract this part in order to recover the correct DFT+U energy description 
      !  which is known from literature (E_{DFT+U} = E_{DFT} + E_{U}).
      !
      !  In general, the hubbard correction results in a penalty to the total energy,
      !  therefore, the DFT+U energy correction E_U should always be postitive.
      
      !  USES

      use basis, only: basis_l
      use dimensions, only: n_spin
      use mpi_tasks, only: myid!,aims_stop_coll
      use localorb_io, only: use_unit

      implicit none

      !  ARGUMENTS

      !  INPUTS

      !  OUTPUTS

      !  AUTHOR
      !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
      !  HISTORY
      !    Release version, FHI-aims (2017).
      !  SOURCE

      ! local variables

      integer  :: i_spin, mi, mj, i_orb, l_orb     
      real*8   :: expec_v_U
      real*8   :: E_U
      !character(*), parameter :: func = 'plus_u'


      if (occ_numbers_initialized .eqv. .false.) then ! energy correction not calculated in the initialization of scf.
         plus_u_energy_correction = 0.0d0            ! energy correction is only calculated when occ_plus_u is available.
      else

         expec_v_U = 0.0d0       ! f_i < Psi_i | v_U | Psi_i >
         E_U = 0.0d0             ! E_U

         do i_orb = 1, n_plusUorb
            l_orb = basis_l(plusUorb_index(i_orb))

            do i_spin = 1,n_spin
               do mi = 1,(2*l_orb+1)
                  do mj = 1,(2*l_orb+1)

                     expec_v_U = expec_v_U - hubbard_u(i_orb) &
                                           * occ_plus_u(i_orb, mi, mj, i_spin) &
                                           * occ_plus_u(i_orb, mj, mi, i_spin)

                     E_U =  E_U - 0.5 * hubbard_u(i_orb) &
                                      * occ_plus_u(i_orb, mi, mj, i_spin) &
                                      * occ_plus_u(i_orb, mj, mi, i_spin)

                  enddo

                    expec_v_U = expec_v_U + hubbard_u(i_orb) &
                                          * ((1-AMF_FLL_array(i_orb))*av_occ_array(i_spin,i_orb) &
                                          + AMF_FLL_array(i_orb) * 0.5) &
                                          * occ_plus_u(i_orb, mi, mi, i_spin)

                    E_U = E_U - 0.5 * hubbard_u(i_orb) &
                                    * ( - 2*occ_plus_u(i_orb, mi, mi, i_spin) * av_occ_array(i_spin,i_orb) &
                                    + av_occ_array(i_spin,i_orb) ** 2 )

               enddo

                  E_U = E_U - 0.5 * hubbard_u(i_orb) &
                                  *  ( - (2*l_orb+1) * AMF_FLL_array(i_orb) &
                                  * av_occ_array(i_spin,i_orb) * (1. - av_occ_array(i_spin,i_orb)) )

            enddo
         enddo !end do i_orb

         ! plus_u_energy_correction is not the E_U term known from literatur. Here, it includes also
         ! the sum of eigenvalues of the hubbard_operator.
         ! Hence, the listed energy in the output is not zero, even we hit the FLL or AMF limit!!!
         ! In general, the DFT+U correction energy (E_U) should always be positive.

         if (E_U .lt. 0.0d0) then
             if(myid .eq. 0)then
                write(use_unit,*)" WARNING: DFT+U ENERGY (E_U) IS NEGATIVE!!!"
             endif
         
            !call aims_stop_coll('DFT+U energy (E_U) is negative', func)
         endif


         plus_u_energy_correction = E_U - expec_v_U
   
      endif !occ_numbers_initialized

   end subroutine plus_u_energy_correction_term

!------------------------------------------------------------------------------------------
!-----------------START OF MATRIX CONTROL--------------------------------------------------
!------------------------------------------------------------------------------------------

   subroutine plus_u_write_occ_mat

      !  PURPOSE
      !
      !  Writes the self-consistent occupation matrix to a file.
      !  The file name is fixed to matrix_control.txt
      
      !  USES
      use basis, only: basis_l
      use mpi_tasks, only: myid
      use localorb_io, only: use_unit, localorb_info
      use dimensions, only: n_spin
      
      implicit none

      !  ARGUMENTS

      !  INPUTS

      !  OUTPUTS

      !  AUTHOR
      !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
      !  HISTORY
      !    Release version, FHI-aims (2017).
      !  SOURCE

      ! local variables

      integer  :: i_spin, mi, mj, i_orb, l_orb
      character*50  :: fmt
      character*150 :: info_str

      if (myid.ne.0) then
         return
      endif

      ! Write message to console.
      write(info_str,'(2X,2A)') 'Writing plus_u occupation matrix for occupation matrix control '
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(2X,2A)') ''
      call localorb_info(info_str,use_unit,'(A)')

      ! Open the restart file for writing.
      open(file = 'occupation_matrix_control.txt', unit = 7, status = 'unknown', form = 'formatted')

      do i_orb = 1, n_plusUorb
         l_orb = basis_l(plusUorb_index(i_orb))
         do i_spin = 1, n_spin
            write(7,*)"  occupation matrix (subspace #",i_orb,", spin ",i_spin,")"
            write(fmt,'(A,I2,A)')" (",l_orb*2+1,"f10.5)"
            write(7,fmt)occ_plus_u(i_orb, 1:(2*l_orb+1), 1:(2*l_orb+1), i_spin)
         enddo
      enddo

      close(7)

   end subroutine plus_u_write_occ_mat

   subroutine plus_u_read_occ_mat

      !  PURPOSE
      !
      !  Reads the self-consistent occupation matrix from a file.
      !  The file name is fixed to occupation_matrix_control.txt
      
      !  USES

      use runtime_choices
      use synchronize_mpi
      use basis, only: basis_l
      use localorb_io, only: use_unit, localorb_info
      use dimensions, only: n_spin
      
      implicit none

      !  ARGUMENTS

      !  INPUTS

      !  OUTPUTS

      !  AUTHOR
      !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
      !  HISTORY
      !    Release version, FHI-aims (2017).
      !  SOURCE

      ! local variables

      integer  :: i_spin, mi, mj, i_orb, l_orb, dim_occ
      character*150 :: info_str
      integer :: mpierr

      if (occ_numbers_initialized .eqv. .false.) then

         if (myid.eq.0) then

            write(info_str,'(2X,A)') ''
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,2A)') 'Reading plus_u occupation matrix from file: occupation_matrix_control.txt'
            call localorb_info(info_str,use_unit,'(A)')

            open(file = 'occupation_matrix_control.txt', unit = 7, status = 'old', form = 'formatted')

            do i_orb = 1, n_plusUorb
               l_orb = basis_l(plusUorb_index(i_orb))
               do i_spin = 1, n_spin
                  read(7,*)
                  read(7,*)occ_plus_u(i_orb, 1:(2*l_orb+1), 1:(2*l_orb+1), i_spin)
               enddo
            enddo


            close(7)

         else
            occ_plus_u = 0.0d0
         endif

      ! MPI broadcast to all other threads
         
         dim_occ = (2*l_max+1)
         p_max = num_shell_max

         call sync_plus_u_occupation_matrix(occ_plus_u, n_plusUorb, dim_occ, dim_occ, n_spin)
         occ_numbers_initialized = .true.

      else
         return
      endif

   end subroutine plus_u_read_occ_mat

!------------------------------------------------------------------------------------------
!-----------------calculate Idempotence Error----------------------------------------------
!------------------------------------------------------------------------------------------

   subroutine plus_u_matrix_error

      !  PURPOSE
      !
      !  calculates the idempotence error of occupation matrix for each
      !  shell which is treated by DFT+U
      !  Error is given by: Tr(n-n*n)
      
      !  In principle, this error tells us how good a projector is for 
      !  a certain system. Almost perfect projector functions should 
      !  give a idempotence error around 0.
      
      !  USES

      use basis, only: basis_l
      use mpi_tasks, only: myid
      use localorb_io, only: use_unit
      use dimensions, only: n_spin
      
      implicit none

      !  ARGUMENTS

      !  INPUTS

      !  OUTPUTS

      !  AUTHOR
      !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
      !  HISTORY
      !    Release version, FHI-aims (2017).
      !  SOURCE

      ! local variables

      integer  :: i_spin, mi, mj, i_orb, l_orb
      real*8   :: trace_1 ! Tr(n)
      real*8   :: trace_2 ! Tr(n*n)
      real*8, allocatable :: plus_u_idempotence_error(:,:)

      allocate(plus_u_idempotence_error(n_plusUorb,n_spin))

      plus_u_idempotence_error = 0.0d0

      trace_1 = 0.0d0
      trace_2 = 0.0d0

      do i_orb = 1, n_plusUorb
         l_orb = basis_l(plusUorb_index(i_orb))
         do i_spin = 1,n_spin
            do mi = 1, (2*l_orb+1)

               trace_1 = trace_1 + occ_plus_u(i_orb, mi, mi, i_spin)

               do mj = 1, (2*l_orb+1)

                  trace_2 = trace_2 + occ_plus_u(i_orb, mi, mj, i_spin) * occ_plus_u(i_orb, mj, mi, i_spin)

               enddo
            enddo

            plus_u_idempotence_error(i_orb,i_spin) = trace_1 - trace_2

            trace_1 = 0.0d0
            trace_2 = 0.0d0

         enddo

      enddo

      if(myid .eq. 0)then
         write(use_unit,*)" -----------------------------------------------------------"
         write(use_unit,*)" ----calculating idempotence error for occupation matrix-----"
         write(use_unit,*)" -----------------------------------------------------------"
         do i_orb = 1, n_plusUorb
            do i_spin = 1,n_spin         
               write(use_unit,*)"  occupation matrix (subspace #",i_orb,", spin ",i_spin,")"
               write(use_unit,*)"  plus_u_idempotence_error :", plus_u_idempotence_error(i_orb,i_spin)
               write(use_unit,*)" -----------------------------------------------------------"
            enddo
         enddo

      endif

      deallocate(plus_u_idempotence_error)


   end subroutine plus_u_matrix_error

!------------------------------------------------------------------------------------------
!-----------------check DFT+U occupation numbers-------------------------------------------
!------------------------------------------------------------------------------------------

   subroutine plus_u_check_occupation

      !  PURPOSE
      !
      !  Performes a check on the occupation numbers of the occupation matrix.
      !  If the occupation number of one hubbard projector is >> 1,
      !  the calculation is stopped because of physical nonsens.
      !  
      !  Reasons could be:
      !
      !   o) bad hubbard projectors for this material
      !   o) the U value is to high
      !   o) bad initial guess for the density  
      !
      !  How to fix it?:
      !
      !   o) define a physical meaningful hubbard projector via the
      !      hubbard_coefficient keyword
      !   o) start from a LDA/GGA calculation and increase the U value
      !      slowly
      !   o) use the matrix_control approach
      !   o) use the DFT+U ramping method
      !   o) use combinations of the above mentioned strategies

      !  USES

      use basis, only: basis_l
      use mpi_tasks, only: myid, aims_stop_coll 
      use localorb_io, only: use_unit
      use dimensions, only: n_spin

      implicit none

      !  ARGUMENTS

      !  INPUTS

      !  OUTPUTS

      !  AUTHOR
      !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
      !  HISTORY
      !    Release version, FHI-aims (2017).
      !  SOURCE

      ! local variables

      integer  :: i_spin, mi, mj, i_orb, l_orb
      integer  :: counter
      logical  :: bad_occupation
      character(*), parameter :: func = 'plus_u'

      if(.not.occ_numbers_initialized)return

      ! loop through the diagonal elements of each occupation matrix
      ! and check the entries

      bad_occupation = .false.

      do i_orb = 1, n_plusUorb
         l_orb = basis_l(plusUorb_index(i_orb))
         do i_spin = 1, n_spin
            do mi = 1, (2*l_orb+1)

               ! In principle there could be any number greater 1 ( or 2 in none spin polarized calculations)
               ! however, sometimes the occupation matrix recovers itself during the run.
               ! In the case of an occupation of > 4.0, the occupation matrix is so messed up that everything
               ! will collapse anyway.
 
               if(occ_plus_u(i_orb, mi, mi, i_spin) .gt. 4.0) then
                  bad_occupation = .true.
               endif

            enddo
         enddo
      enddo

      if(bad_occupation .eqv. .true.) then

         if(myid .eq. 0)then
            write(use_unit,*)"#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#"
            write(use_unit,*)"    ERROR: Bad DFT+U occupation matrix"
            write(use_unit,*)"---------------------------------------------"
            write(use_unit,*)"    some of the matrix elements are > 4"     
            write(use_unit,*)"#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#"
            write(use_unit,*)" "
            write(use_unit,*)" Reasons could be:"
            write(use_unit,*)" "
            write(use_unit,*)"   o) bad hubbard projectors for this material"
            write(use_unit,*)"   o) the U value is to high"
            write(use_unit,*)"   o) bad initial guess for the density"
            write(use_unit,*)" "
            write(use_unit,*)" How to fix it?:"
            write(use_unit,*)" "
            write(use_unit,*)"   o) define a more physical meaningful hubbard projector via the"
            write(use_unit,*)"      hubbard_coefficient keyword"
            write(use_unit,*)"   o) start from a LDA/GGA calculation and increase the U value"
            write(use_unit,*)"      slowly"
            write(use_unit,*)"   o) use the matrix_control approach"
            write(use_unit,*)"   o) DFT+U ramping could also help"
            write(use_unit,*)"   o) use combinations of the above mentioned strategies"
            write(use_unit,*)"   o) our observation is, that small basis sets performe better &
                                    &with DFT+U"
            write(use_unit,*)" "
            write(use_unit,*)" otherwise"
            write(use_unit,*)" "
            !write(use_unit,*)" MAY THE FORCE BE WITH +U"
            write(use_unit,*)" LET THE HATE FLOW THROUGH +U"

         endif

         call aims_stop_coll('Bad occupation matrix in DFT+U', func)

      endif
         

   end subroutine plus_u_check_occupation
 
 !------------------------------------------------------------------------------------------
 !----------------------Eigenvalues of occupation numers------------------------------------
 !------------------------------------------------------------------------------------------
   

   subroutine plus_u_eigenvalues

      !  PURPOSE
      !
      !  Calculates the eigenvalues of the DFT+U occupation matrix.
      !
      !  In principle, the diagonal elements of the occupation matrix
      !  are pretty good approximations for the occupation numbers of
      !  the defined correlated subspace. However, in general, the 
      !  eigenvalues are the actual occupation numbers.
      !
      !  Eigenvalues are written to standard output, only.
      !
      !  Eigenvalues are calculated after scf convergence.
      !  

      !  USES

      use basis, only: basis_l
      use mpi_tasks, only: myid, aims_stop_coll
      use localorb_io, only: use_unit
      use dimensions, only: n_spin

      implicit none

      !  ARGUMENTS

      !  INPUTS

      !  OUTPUTS

      !  AUTHOR
      !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
      !  HISTORY
      !    Release version, FHI-aims (2017).
      !  SOURCE

      ! local variables

      integer                 :: N, LDVL, LDVR, LDA
      integer                 :: INFO, LWORK
      integer                 :: LWMAX
      parameter                  ( LWMAX = 1000 )
      integer                 :: i_orb, l_orb, i_spin
      character(*), parameter :: func = 'plus_u'
      character*20            :: fmt

      real*8, allocatable     :: VL(:,:), VR(:,:), WR(:), WI(:), A(:,:)
      real*8                  :: WORK(LWMAX)

      if(myid .eq. 0)then
         write(use_unit,*)" -----------------------------------------------------------"
         write(use_unit,*)" ----calculating eigenvalues of DFT+U occupation matrix-----"
         write(use_unit,*)" -----------------------------------------------------------"
      endif

      do i_orb = 1, n_plusUorb
         l_orb = basis_l(plusUorb_index(i_orb))

         N    = 2*l_orb + 1
         LDA  = N
         LDVL = N
         LDVR = N

         ! in principle each occupation matrix can have a different dimension
         ! allocating working arrays
         allocate(A(LDA,N))
         allocate(VL(LDVL,N))
         allocate(VR(LDVR,N))
         allocate(WR(N))
         allocate(WI(N))



         do i_spin = 1, n_spin

            WR = 0.0d0
            WI = 0.0d0
            VL = 0.0d0
            VR = 0.0d0

            A(1:LDA,1:N) = occ_plus_u(i_orb, 1:LDA, 1:N, i_spin)

            ! first: get correct length of WORK
            LWORK = -1
            call DGEEV('N', 'N', N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
            LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

            ! second: calculate eigenvalues
            call DGEEV('N', 'N', N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)


            if( INFO .gt. 0) then

               call aims_stop_coll('Failed to compute eigenvalues of DFT+U occupation matrix', func)
               
            endif

            if(myid .eq. 0)then
               write(use_unit,*)"  eigenvalues (subspace #",i_orb,", spin ",i_spin,")"
               write(fmt,'(A,I2,A)')" (",l_orb*2+1,"f10.5)"
               write(use_unit,fmt) WR
               write(use_unit,*)" -----------------------------------------------------------"
            endif

         enddo

         deallocate(A)
         deallocate(VL)
         deallocate(VR)
         deallocate(WR)
         deallocate(WI)

      enddo

   end subroutine plus_u_eigenvalues

  
end module plus_u
