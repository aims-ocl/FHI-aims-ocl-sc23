!****s* FHI-aims/evaluate_ovlp_shell_p0
!  NAME
!   evaluate_ovlp_shell_p0
!  SYNOPSIS

      subroutine evaluate_ovlp_shell_p0 &
           ( n_points, partition,  &
           n_compute_a, n_compute_c, wave, matrix_shell, n_basis_list &
           )

!  PURPOSE
!  The subroutine evaluates the overlap integral contribution
!  of several integration points, and adds it to the overall overlap matrix.
!
!  USES

!      use dimensions
      implicit none

!  ARGUMENTS

      integer :: n_points
      integer :: n_basis_list

      real*8 :: partition(n_points)

      integer :: n_compute_a
      integer :: n_compute_c

      real*8 :: wave( n_basis_list,n_points)
      real*8 :: matrix_shell( n_compute_c, n_compute_a)


!  INPUTS
!  o n_points -- number of grid points
!  o partition -- values of partition function
!  o n_compute_a -- number of relevan basis functions
!  o n_compute_c -- number of relevan basis functions
!  o wave -- values of basis functions
!  o n_basis_list -- list of relevan basis functions
! 
!  OUTPUT
!  o matrix_shell -- results of the overlap multiplications
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
!
!  local variables

!     auxiliary matrices for Level 3 Blas matrix multiplications

! ??? Need we store partition * wave in an extra matrix?
! ??? Need BLAS descriptions to go on ...
! ??? Does it matter that n_points (the summation index in the matrix multiplication) is
! ??? the outer one, not the inner one, in wave?
!     Yes, we need, but as long as we are inside the cache, no harm should follow...

      real*8, dimension( n_compute_c,n_points) ::  wave_compute_c 
      real*8, dimension( n_compute_a,n_points) ::  wave_compute_a 

!     counters
      integer :: i_point
!     retired counters
!      integer :: i_basis_1
!      integer :: i_basis_2
!      integer :: i_index

!      integer :: i_compute
!      integer :: i_compute_1
!      integer :: i_compute_2
!      integer :: i_offset
!      integer :: i_index_real

!     first, allocate
      
!     begin work

!     Now Integrate:
!     We only integrate the upper triangle -- all matrices must be Hermitian.

!     We must use the packed format for all matrices, due to speed. 
!     This means we must 
!     go through each matrix column | phi_j > in the outer loop, 
!     and through each line < phi_i | in the inner loop, where
!     column j = i_basis_2
!     line   i = i_basis_1

!     NEC_CBq
      if (n_compute_c.eq.n_compute_a) then

        do i_point = 1, n_points, 1
           wave_compute_c(1:n_compute_c, i_point) =  &
                sqrt(partition(i_point))* &
                wave(1:n_compute_c, i_point)
        end do

!        call dgemm('N', 'T', n_compute_c, n_compute_c,  &
!             n_points, 1.0d0,  &
!             wave_compute_c, n_compute_c, wave_compute_c,  &
!             n_compute_c, 0.0d0, matrix_shell, n_compute_c )
! 
! Lapack Routine DSYRK performs the symmetric rank k operations:
! matrix_shell := 1.0d0*wave_compute_c*wave_compute_c^T + 0.0d0*matrix_shell

         call dsyrk('U','N', n_compute_c,n_points, 1.0d0, &
                     wave_compute_c, n_compute_c,  0.0d0, matrix_shell, n_compute_c )
      else

        do i_point = 1, n_points, 1
           wave_compute_c(1:n_compute_c, i_point) =  &
                sqrt(partition(i_point))* &
                wave(1:n_compute_c, i_point)

           wave_compute_a(1:n_compute_a,i_point) =  wave_compute_c(1:n_compute_a, i_point)

        enddo
! Lapack Routine DGEMM performs the matrix-matrix operation:
! matrix_shell := 1.0d0*wave_compute_c*wave_compute_a^T + 0.0d0*matrix_shell

        call dgemm('N', 'T', n_compute_c, n_compute_a,  &
          n_points, 1.0d0,  &
          wave_compute_c, n_compute_c, wave_compute_a,  &
          n_compute_a, 0.0d0, matrix_shell, n_compute_c )

      end if


    end subroutine evaluate_ovlp_shell_p0
!---------------------------------------------------------------------
!****** 
