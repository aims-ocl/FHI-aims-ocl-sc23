!****s*  FHI-aims/evaluate_xc_matr_shell
!  NAME
!    evaluate_xc_matr_shell
!  SYNOPSIS

      subroutine evaluate_xc_matr_shell &
           ( n_points, partition, &
             n_compute, i_basis, xc_times_psi, &
             wave, xc_matr &
           )

!  PURPOSE
!  Subroutine evaluate_xc_matr_shell evaluates the kinetic energy matrix
!  element between between two basis functions
!  for a radial shell, and adds it to the overall
!  kinetic energy matrix.
!
!  USES

      use dimensions

      implicit none

!  ARGUMENTS

      integer :: n_points

      integer n_compute
      integer i_basis(n_compute)

      real*8 partition(n_points)
      real*8 wave(n_basis,n_points)
      real*8 xc_times_psi(n_basis,n_points)

      real*8 xc_matr( n_basis,n_basis )

!  INPUTS
!  o n_points -- number of relevant points in the current integration shell 
!  o n_compute -- number of non-zero basis functions in the current integration shell
!  o i_basis -- array, specifies the non-zero basis functions  
!  o partition -- the values of the partition function (for integration) at each 
!         spatial grid point
!  o wave -- value of the basis functions
!  o xc_times_psi -- the local exchange-correlation potential times basis functions
!
!  OUTPUTS
!  o xc_matr -- the matrix of the exchange-correlation potential within basis functions
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

      real*8 aux_xc_matrix &
             (n_compute,n_compute)

      real*8 wave_compute(n_compute,n_points)

!     counters

      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_index

      integer :: i_atom

      integer :: i_compute
      integer :: i_compute_1
      integer :: i_compute_2
      integer :: i_offset
      integer :: i_index_real
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

!     First we take only the waves that are associated with atom i

      do i_point = 1, n_points, 1

         wave_compute(1:n_compute, i_point) = &
              wave(1:n_compute, i_point)*partition(i_point)

      enddo

!     now, integrate, initialisation of aux. overlap matrix is needed
!     since we add to it the values in the second stage below

      aux_xc_matrix = 0.0d0

!     compute wave*(psi*psi) and add this to aux. overlap matrix
      call dgemm('N', 'T', n_compute, n_compute, &
        n_points, 1.0d0, &
        wave_compute, n_compute, xc_times_psi, &
        n_basis, 0.0d0, aux_xc_matrix, &
        n_compute )
!       do i_compute =1, n_compute
!         do i_compute_1 =1, n_compute
!           do i_point =1, n_points
!            aux_xc_matrix (i_compute, i_compute_1) =
!     +      aux_xc_matrix (i_compute, i_compute_1) +
!     +        wave_compute(i_compute,i_point)*
!     +        xc_times_psi(i_compute_1,i_point)
!     +        wave(i_compute_1,i_point)
!           enddo
!         enddo
!       enddo
!      now add the aux. overlap matrix to the actual overlap matrix
!  this requires translating between the actually computed matrix elements and the
!  full overlap matrix ...
      if (use_density_gradient) then

        do i_compute = 1, n_compute
            do i_compute_1 = 1, i_compute -1

             aux_xc_matrix(i_compute, i_compute_1) = &
             0.5d0 * ( aux_xc_matrix(i_compute, i_compute_1) + &
                       aux_xc_matrix(i_compute_1, i_compute) )


             aux_xc_matrix(i_compute_1, i_compute) = &
              aux_xc_matrix(i_compute, i_compute_1)

            enddo
         enddo

      endif

      i_index = 0
      i_index_real = 0

      do i_compute =1, n_compute, 1

         i_basis_1=i_basis(i_compute)
         i_index=  (i_basis(i_compute)-1)*i_basis(i_compute)/2

         do i_compute_1 = 1, n_compute, 1

              i_basis_2=i_basis(i_compute_1)

              i_index_real =  i_index + &
                   i_basis (i_compute_1)

              xc_matr(i_basis_1,i_basis_2) = &
                xc_matr(i_basis_1,i_basis_2) &
               + aux_xc_matrix(i_compute, i_compute_1)

!     end of i_compute_1
           enddo

!     end of i_compute
       enddo
      end subroutine evaluate_xc_matr_shell
!---------------------------------------------------------------------
!***************
