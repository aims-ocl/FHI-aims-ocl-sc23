!****s* FHI-aims/evaluate_zora_vector_p1
!  NAME
!   evaluate_zora_vector_p1
!  SYNOPSIS

      subroutine evaluate_zora_vector_p1 & 
           ( zora_operator,  partition_tab, gradient_basis_wave, n_compute, &
           zora_vector1, zora_vector2, n_basis_list, t_zora &
           )

!  PURPOSE
!  Subroutine add_zora_matrix adds the ZORA difference term to the
!  non-relativistic Hamiltonian matrix.
!
!  The idea, here, is to do only the minimum number of operations
!  inside both loops i_basis_1 and i_compute_2 ... these operations
!  lead to a quadratic dependence of the integrations on n_basis,
!  and should be avoided as much as possible.
!
!  USES

      use dimensions
      use runtime_choices
      implicit none

!  ARGUMENTS

      real*8  :: zora_operator
      real*8  :: partition_tab
      integer :: n_compute
      integer :: n_basis_list
      real*8  :: gradient_basis_wave(n_compute, 3)
      real*8  :: zora_vector1( n_basis_list, 3 )
      real*8  :: zora_vector2( n_basis_list, 3 )
      logical :: t_zora

!  INPUTS
!    o  zora_operator -- prefactor for zora
!    o  partition_tab -- values of partition function
!    o  n_compute -- number of non-zero basis functions
!    o  n_basis_list -- total number of basis functions including periodic mirror images
!    o  gradient_basis_wave -- gradient of basis functions
!
!  OUTPUT
!   o  zora_vector1 -- first zora vector
!   o  zora_vector2 -- second zora vector
!   o  t_zora -- is this grid point relevant for relativistic calculations?
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

!     counters

      integer :: i_coord
      integer :: i_compute_1


!     begin work

      do i_compute_1 = 1, n_compute, 1

         do i_coord = 1,3,1
            
            zora_vector1(i_compute_1,i_coord) = &
                 zora_operator* &
                 gradient_basis_wave(i_compute_1, i_coord)
            
            zora_vector2(i_compute_1,i_coord) = &
                 partition_tab*gradient_basis_wave(i_compute_1, i_coord)
            
         end do
      enddo

      t_zora = .false.
      do i_coord = 1,3,1
         if ( ( maxval(abs(zora_vector1(1:n_compute,i_coord))) * & 
              maxval(abs(zora_vector2(1:n_compute,i_coord))) ) .gt. zora_threshold ) then
            t_zora = .true.
         end if
      end do
      
         
    end subroutine evaluate_zora_vector_p1
!---------------------------------------------------------------------
!******
