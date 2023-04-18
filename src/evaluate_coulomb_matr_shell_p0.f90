!****s* FHI-aims/evaluate_coulomb_matr_shell_p0
!  NAME
!   evaluate_coulomb_matr_shell_p0
!  SYNOPSIS

      subroutine evaluate_coulomb_matr_shell_p0 &
           ( n_points, i_task, i_atom, &
             partition, &
             n_compute, i_prodbas, &
             n_compute_current, &
             n_loc_compute, i_loc_prodbas, &
             n_loc_compute_current, &
             v_times_waves, &
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
      use mpi_tasks
      use constants

      implicit none

!  ARGUMENTS

      integer :: n_points
      integer :: i_task
      integer :: i_atom

      integer ::  n_compute
      integer ::  n_compute_current
      integer ::  n_loc_compute
      integer ::  n_loc_compute_current
      integer ::  i_prodbas(n_compute)
      integer ::  i_loc_prodbas(n_loc_compute)

      real*8 partition(n_atoms,2,n_points)
      real*8 wave(n_compute,n_points)
      real*8 v_times_waves(n_loc_compute,n_points)

      real*8 coulomb_matr( n_basbas,n_loc_prodbas )

!  INPUTS
!  o n_points -- number of relevant points in the current integration shell or batch
!  o i_task -- the current task number
!  o  n_compute -- the number of non-zero basis function
!  o  n_compute_current -- the non-zero basis funciton belonging to the current atom
!  o  n_loc_compute -- the number of non-zero auxiliary basis function
!  o  n_loc_compute_current -- the number of non-zero auxiliary basis function assoicated
!             with the current atom
!  o  i_prodbas -- array, specifies these non-zero basis function
!            current atom
!  o  i_loc_prodbas -- specifies all the non-zero auxiliary basis fucntion
!  o  i_loc_prodbas_current -- specifies non-zero auxiliary basis function assoicated
!             with the current atom
!  o  partition -- the value the two-center partition function
!  o  wave -- the value of auxiliary basis functions at points in the current batch
!  o  v_times_waves -- the value of the v times the auxiliary basis function
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
      integer ::  n_loc_compute_other
      integer ::  i_prodbas_other(n_compute)

!      real*8 :: aux_coulomb_matrix
!     +       (n_compute_current,n_loc_compute)


!      real*8 :: aux_coulomb_matrix_1
!     +       (n_compute, n_loc_compute_current)

      real*8, allocatable, dimension(:,:) ::   aux_coulomb_matrix
      real*8, allocatable, dimension(:,:) ::   aux_coulomb_matrix_1
      real*8, allocatable, dimension(:,:) ::   wave_compute
      real*8, allocatable, dimension(:,:) ::   v_waves_compute
!      real*8 wave_compute(n_compute,n_points)
!      real*8 v_waves_compute(n_loc_compute,n_points)

!     counters

      integer :: i_prodbas_1
      integer :: i_prodbas_2
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
!        if (i_prodbas(1).eq.1) then
!          write(use_unit,*) i_point, psi_times_psi(1,i_point),
!     +                 psi_times_psi(n_basis,i_point)
!        end if
!      enddo
!test end

!     First run: we multiply the product basis on the current atom
!     with all the "exchange potential" of the product basis.


      n_compute_other = n_compute - n_compute_current
      n_loc_compute_other = n_loc_compute - n_loc_compute_current

      if(.not.allocated(aux_coulomb_matrix)) then
        allocate(aux_coulomb_matrix(n_compute_current,n_loc_compute_other))
      endif
      if(.not.allocated(aux_coulomb_matrix_1)) then
        allocate(aux_coulomb_matrix_1(n_compute_other, n_loc_compute_current))
      endif
      if(.not.allocated(wave_compute)) then
        allocate(wave_compute(n_compute_other,n_points))
      endif
      if(.not.allocated(v_waves_compute)) then
        allocate(v_waves_compute(n_loc_compute_other,n_points))
      endif

!      wave_compute(:,:) = 0.d0
!      wave_compute(1:i_compute, :) =  wave(i_prodbas_1, :)

      v_waves_compute(:,:) = 0.d0
      do i_compute = 1, n_loc_compute_other, 1

         i_compute_1 = i_compute + n_loc_compute_current  
         i_prodbas_1 = i_loc_prodbas(i_compute_1)
         i_basbas = map_prodbas(i_prodbas_1, i_task)
         if (i_basbas .gt. 0) then
           i_atom_1 = basbas_atom(i_basbas)
           do i_point = 1, n_points
             v_waves_compute(i_compute,i_point) &
              = v_times_waves(i_compute_1, i_point) * &
                partition(i_atom_1,1,i_point)
           enddo
         endif
       enddo

      aux_coulomb_matrix = 0.0d0
!     compute wave*v*wave and add this to coulomb matrix
      call dgemm('N', 'T', n_compute_current, n_loc_compute_other, &
        n_points, 1.0d0, &
        wave, n_compute, v_waves_compute, &
        n_loc_compute_other, 0.0d0, aux_coulomb_matrix, &
        n_compute_current )
!        do i_compute=1, n_loc_compute_other,1
!         do i_compute_1=1, n_compute_current,1
!           do i_point=1,n_points,1
!             aux_coulomb_matrix(i_compute_1, i_compute) = &
!              aux_coulomb_matrix(i_compute_1, i_compute) + &
!              wave(i_compute_1,i_point)*v_waves_compute(i_compute,i_point)
!            if(myid.eq.0) then
!             write(use_unit,'(3I4,4f12.8)') i_compute,i_compute_1,i_point, &
!              wave(i_compute_1,i_point),v_waves_compute(i_compute,i_point),&
!              v_times_waves(i_compute_1+n_loc_compute_current,i_point), &
!              aux_coulomb_matrix(i_compute_1,i_compute)
!            endif
!           enddo
!          enddo
!        enddo
!     now add the aux. coulomb matrix to the actual coulomb matrix
      do i_compute =1, n_loc_compute_other, 1

         i_prodbas_1=i_loc_prodbas(i_compute+n_loc_compute_current)
         do i_compute_1 = 1, n_compute_current, 1

              i_prodbas_2=i_prodbas(i_compute_1)

              coulomb_matr(i_prodbas_2,i_prodbas_1) = &
                coulomb_matr(i_prodbas_2,i_prodbas_1) &
               + aux_coulomb_matrix(i_compute_1, i_compute)

!     end of i_compute_1
           enddo
!     end of i_compute
       enddo

!     Second run: we multiply all the product basis with all the 
!    "exchange potential" (v_times_wave) of the product basis at current atom

!      n_compute_other = 0
!      i_prodbas_other = 0
      wave_compute(:,:)  = 0.d0
      do i_compute = 1, n_compute_other, 1

         i_prodbas_1 = i_prodbas(i_compute+n_compute_current)
         i_atom_1 = basbas_atom(i_prodbas_1)

         do i_point = 1, n_points, 1
           wave_compute(i_compute, i_point) = &
           wave(i_compute+n_compute_current, i_point) * &
           partition(i_atom_1,2,i_point)
         enddo
      enddo

!      v_waves_compute(:,:) = 0.d0
!      do i_compute = 1, n_loc_compute_current
!         i_prodbas_1 = i_loc_prodbas_current(i_compute)
!         v_waves_compute(i_compute,:) &
!             = v_times_waves(i_prodbas_1, :)
!      enddo

!      aux_coulomb_matrix_1 = 0.d0
!     compute wave*v*wave and add this to coulomb matrix
      if(n_loc_compute_current.gt.0) then
        call dgemm('N', 'T', n_compute_other, n_loc_compute_current, &
         n_points, 1.0d0, &
         wave_compute, n_compute_other, v_times_waves, &
         n_loc_compute, 0.0d0, aux_coulomb_matrix_1, &
         n_compute_other )

!      now add the aux. coulomb matrix to the actual coulomb matrix
         do i_compute =1, n_loc_compute_current, 1

           i_prodbas_1=i_loc_prodbas(i_compute)
           do i_compute_1 = 1, n_compute_other, 1

              i_prodbas_2=i_prodbas(i_compute_1+n_compute_current)

              coulomb_matr(i_prodbas_2,i_prodbas_1) = &
                coulomb_matr(i_prodbas_2,i_prodbas_1) &
               + aux_coulomb_matrix_1(i_compute_1, i_compute)

!              coulomb_matr(i_prodbas_1,i_prodbas_2) = &
!                coulomb_matr(i_prodbas_1,i_prodbas_2) &
!               + aux_coulomb_matrix_1(i_compute_1, i_compute)

!     end of i_compute_1
          enddo
!     end of i_compute
        enddo
!     end of if n_loc_compute_current >0
      endif 
!      write(use_unit,'(4f18.10)') coulomb_matr(1,2), coulomb_matr(2,1), &
!         coulomb_matr(1,6),  coulomb_matr(6,1)

      if(allocated(aux_coulomb_matrix)) then
        deallocate(aux_coulomb_matrix)
      endif
      if(allocated(aux_coulomb_matrix)) then
        deallocate(aux_coulomb_matrix_1)
      endif
      if(allocated(wave_compute)) then
        deallocate(wave_compute)
      endif
      if(allocated(v_waves_compute)) then
        deallocate(v_waves_compute)
      endif

      end subroutine evaluate_coulomb_matr_shell_p0
!---------------------------------------------------------------------
!******
