!****s* FHI-aims/evaluate_coulomb_matr_shell_v0
!  NAME
!   evaluate_coulomb_matr_shell_v0
!  SYNOPSIS

      subroutine evaluate_coulomb_matr_shell_v0 &
           ( n_points, i_task, i_atom, &
             partition, &
             n_compute, i_basis, &
             n_compute_current, i_basis_current, &
             n_prod_compute, i_prodbas, &
             n_prodbas_current, i_prodbas_current, &
             v_times_psi, &
             wave, coulomb_matr &
           )

!  PURPOSE
!  Subroutine evaluate_coulomb_matr_shell evaluates the kinetic energy matrix
!  element between between two basis functions
!  for a radial shell, and adds it to the overall
!  kinetic energy matrix.
!
!  USES

      use dimensions
      use basis
      use prodbas
      use constants

      implicit none

!  ARGUMENTS

      integer :: n_points
      integer :: i_task
      integer :: i_atom

      integer ::  n_compute
      integer ::  n_compute_current
      integer ::  n_prod_compute
      integer ::  n_prodbas_current
      integer ::  i_basis(n_compute)
      integer ::  i_basis_current(n_compute_current)
      integer ::  i_prodbas(n_prod_compute)
      integer ::  i_prodbas_current(n_prodbas_current)

      real*8 partition(n_atoms,n_points)
      real*8 wave(n_basbas,n_points)
      real*8 v_times_psi(n_loc_prodbas,n_points)

      real*8 coulomb_matr( n_basbas,n_loc_prodbas )

!  INPUTS
!  o n_points -- number of relevant points in the current integration shell or batch
!  o i_task -- the current task number
!  o  n_compute -- the number of non-zero basis function
!  o  n_compute_current -- the non-zero basis funciton belonging to the current atom
!  o  n_prod_compute -- the number of non-zero auxiliary basis function
!  o  n_prodbas_current -- the number of non-zero auxiliary basis function assoicated
!             with the current atom
!  o  i_basis -- array, specifies these non-zero basis function
!  o  i_basis_current -- specifies the non-zero basis funciton belonging to the
!            current atom
!  o  i_prodbas -- specifies all the non-zero auxiliary basis fucntion
!  o  i_prodbas_current -- specifies non-zero auxiliary basis function assoicated
!             with the current atom
!  o  partition -- the value the two-center partition function
!  o  wave -- the value of auxiliary basis functions at points in the current batch
!  o  v_times_psi -- the value of the v times the auxiliary basis function
!  OUTPUTS
!  o  coulomb_matr -- the Coulomb interaciton matrix within the auxiliary basis 
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

      integer ::  n_compute_other
      integer ::  i_basis_other(n_compute)

!      real*8 :: aux_coulomb_matrix
!     +       (n_compute_current,n_prod_compute)


!      real*8 :: aux_coulomb_matrix_1
!     +       (n_compute, n_prodbas_current)

      real*8, allocatable, dimension(:,:) ::   aux_coulomb_matrix
      real*8, allocatable, dimension(:,:) ::   aux_coulomb_matrix_1
      real*8, allocatable, dimension(:,:) ::   wave_compute
      real*8, allocatable, dimension(:,:) ::   v_waves_compute
!      real*8 wave_compute(n_compute,n_points)
!      real*8 v_waves_compute(n_prod_compute,n_points)

!     counters

      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_basbas

      integer :: i_atom_1

      integer :: i_compute
      integer :: i_compute_1
      integer :: i_point

!     begin work

!     Now Integrate:
!     We only integrate the upper triangle -- all matrices must be Hermitian.

!test
!      write(use_unit,*) "psi_times_psi(1): "
!      do i_point = 1, n_points, 1
!        if (i_basis(1).eq.1) then
!          write(use_unit,*) i_point, psi_times_psi(1,i_point),
!     +                 psi_times_psi(n_basis,i_point)
!        end if
!      enddo
!test end

!     First run: we multiply the product basis on the current atom
!     with all the "exchange potential" of the product basis.

      if(.not.allocated(aux_coulomb_matrix)) then
        allocate(aux_coulomb_matrix(n_compute_current,n_prod_compute))
      endif
      if(.not.allocated(aux_coulomb_matrix_1)) then
        allocate(aux_coulomb_matrix_1(n_compute, n_prodbas_current))
      endif
      if(.not.allocated(wave_compute)) then
        allocate(wave_compute(n_compute,n_points))
      endif
      if(.not.allocated(v_waves_compute)) then
        allocate(v_waves_compute(n_prod_compute,n_points))
      endif

      wave_compute(:,:) = 0.d0
      do i_compute = 1, n_compute_current, 1
         i_basis_1 = i_basis_current(i_compute)
         wave_compute(i_compute, :) =  wave(i_basis_1, :)
      enddo

      v_waves_compute(:,:) = 0.d0
      do i_compute = 1, n_prod_compute

         i_basis_1 = i_prodbas(i_compute)
         i_basbas = map_prodbas(i_basis_1, i_task)
         if (i_basbas .gt. 0) then
           i_atom_1 = basbas_atom(i_basbas)
           do i_point = 1, n_points
             v_waves_compute(i_compute,i_point) &
              = v_times_psi(i_compute, i_point) * &
                partition(i_atom_1,i_point)
           enddo
         endif
       enddo

!      allocate(aux_coulomb_matrix(n_compute_current, n_prod_compute))
!      allocate(aux_coulomb_matrix_1(n_compute, n_prodbas_current))

!      aux_coulomb_matrix = 0.0d0

!     compute wave*v*wave and add this to coulomb matrix
      call dgemm('N', 'T', n_compute_current, n_prod_compute, &
           n_points, 1.0d0, &
           wave_compute, n_compute, v_waves_compute, &
           n_prod_compute, 0.0d0, aux_coulomb_matrix, &
           n_compute_current )
!       do i_compute =1, n_compute
!         do i_compute_1 =1, n_compute
!           do i_point =1, n_points
!            aux_coulomb_matrix (i_compute, i_compute_1) =
!     +      aux_coulomb_matrix (i_compute, i_compute_1) +
!     +        wave_compute(i_compute,i_point)*
!     +        v_times_psi(i_compute_1,i_point)
!           enddo
!         enddo
!       enddo
!      now add the aux. coulomb matrix to the actual coulomb matrix
      do i_compute =1, n_prod_compute, 1

         i_basis_1=i_prodbas(i_compute)
         do i_compute_1 = 1, n_compute_current, 1

              i_basis_2=i_basis(i_basis_current(i_compute_1))

              coulomb_matr(i_basis_2,i_basis_1) = &
                coulomb_matr(i_basis_2,i_basis_1) &
               + aux_coulomb_matrix(i_compute_1, i_compute)

!     end of i_compute_1
           enddo
!     end of i_compute
       enddo

!     Second run: we multiply the product basis on the current atom
!     with all the "exchange potential" of the product basis.

      n_compute_other = 0
      i_basis_other = 0
      wave_compute(:,:)  = 0.d0
      do i_compute = 1, n_compute, 1

         i_basis_1 = i_basis(i_compute)
         i_atom_1 = basbas_atom(i_basis_1)

         if(i_atom_1 .ne. i_atom) then

            n_compute_other = n_compute_other + 1
            i_basis_other(n_compute_other) = i_basis_1

            do i_point = 1, n_points, 1
              wave_compute(n_compute_other, i_point) = &
               wave(i_compute, i_point) * &
               partition(i_atom_1, i_point)
            enddo
         endif
      enddo

      v_waves_compute(:,:) = 0.d0
      do i_compute = 1, n_prodbas_current
         i_basis_1 = i_prodbas_current(i_compute)
         v_waves_compute(i_compute,:) &
             = v_times_psi(i_basis_1, :)
      enddo

!      aux_coulomb_matrix_1 = 0.d0
!     compute wave*v*wave and add this to coulomb matrix
      call dgemm('N', 'T', n_compute_other, n_prodbas_current, &
           n_points, 1.0d0, &
           wave_compute, n_compute, v_waves_compute, &
           n_prod_compute, 0.0d0, aux_coulomb_matrix_1, &
           n_compute )

!      now add the aux. coulomb matrix to the actual coulomb matrix
      do i_compute =1, n_prodbas_current, 1

         i_basis_1=i_prodbas(i_prodbas_current(i_compute))
         do i_compute_1 = 1, n_compute_other, 1

              i_basis_2=i_basis_other(i_compute_1)

              coulomb_matr(i_basis_2,i_basis_1) = &
                coulomb_matr(i_basis_2,i_basis_1) &
               + aux_coulomb_matrix_1(i_compute_1, i_compute)

!     end of i_compute_1
           enddo
!     end of i_compute
      enddo

      if(allocated(aux_coulomb_matrix)) then
        deallocate(aux_coulomb_matrix)
      endif
      if(allocated(aux_coulomb_matrix_1)) then
        deallocate(aux_coulomb_matrix_1)
      endif
      if(allocated(wave_compute)) then
        deallocate(wave_compute)
      endif
      if(allocated(v_waves_compute)) then
        deallocate(v_waves_compute)
      endif

      end subroutine evaluate_coulomb_matr_shell_v0
!---------------------------------------------------------------------
!******
