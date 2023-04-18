!****s* FHI-aims/evaluate_waves_p2
!  NAME
!   evaluate_waves_p2
!  SYNOPSIS

      subroutine evaluate_waves_p2 &
      ( n_compute, n_compute_atoms, n_compute_fns, &
        l_ylm_max, ylm_tab, one_over_dist_tab,  &
        radial_wave, wave, &
        rad_index, wave_index, l_index, l_count, fn_atom, &
        n_zero_compute, zero_index_point )

!  PURPOSE
!     Prepares only those wave function components which are later needed
!     to evaluate the Hamiltonian or density.
!     There is previous versions of this routines: v0 -> v1 -> p0 -> p2
!
!  USES
      implicit none

!  ARGUMENTS


      integer, intent(IN) :: n_compute_atoms
      integer, intent(IN) :: n_compute_fns
      integer, intent(IN) :: n_compute
      integer, intent(IN) :: l_ylm_max
      real*8, intent(IN) :: ylm_tab ((l_ylm_max+1)**2, n_compute_atoms )
      real*8, intent(IN) :: one_over_dist_tab(n_compute_atoms)
      real*8, intent(IN) :: radial_wave(n_compute_fns)
      integer, intent(IN) :: rad_index(n_compute_atoms)
      integer, intent(IN) :: wave_index(n_compute_fns)
      integer, intent(IN) :: l_index(n_compute_fns)
      integer, intent(IN) :: l_count(n_compute_fns)
      integer, intent(IN) :: fn_atom(n_compute_fns)
      integer, intent(IN) :: n_zero_compute
      integer, intent(IN) :: zero_index_point(n_compute)
      real*8, intent(OUT)  :: wave(n_compute)

!  INPUTS
!    o n_compute_atoms -- number of relevant atoms
!    o n_compute_fns -- number of non-zero basis fns
!    o n_compute -- number of non-zero basis functions 
!    o n_basis_list -- total number of basis functions in whole grid
!    o l_ylm_max -- maximum of l
!    o ylm_tab -- spherical harmonic functions
!    o one_over_dist_tab -- 1/r
!    indices for basis functions that are nonzero at current point
!     o rad_index
!     o wave_index
!     o l_count
!     o fn_atom
!    indices for known zero basis functions at current point
!     o n_zero_compute
!     o zero_index_point
!
!  OUTPUT
!   o wave -- total basis functions. wave is defined for only one point here, but may be stored on the outside
!         of the subroutine
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

      real*8  :: aux_radial (n_compute_fns)

!     counters

      integer :: i_compute_point

      integer :: i_compute
      integer :: i_compute_fn

      integer :: i_compute_atom

      integer :: index_start
      integer :: index_end

!     begin work

!     do work on radial functions first
    
      

      index_start = 1
      do i_compute_atom = 1, n_compute_atoms, 1

        index_end = rad_index(i_compute_atom)

        aux_radial   ( index_start:index_end ) = &
          radial_wave( index_start:index_end ) * &
          one_over_dist_tab( i_compute_atom )

        index_start = index_end+1

      enddo

!      write(use_unit,*) 
!      write(use_unit,*) aux_radial(1:n_compute_fns)
!      write(use_unit,*)

!     tabulate total wave function value for each basis function

      ! first, the nonzero functions
      do i_compute_fn = 1, n_compute_fns, 1

            call mul_vec ( &
                 wave(wave_index(i_compute_fn)), l_count(i_compute_fn)+1, &
                 ylm_tab(l_index(i_compute_fn),fn_atom(i_compute_fn)), &
                 aux_radial(i_compute_fn) &
                 )

!        wave ( wave_index(i_compute_fn):wave_index(i_compute_fn)+l_count(i_compute_fn) ) = &
!          ylm_tab(l_index(i_compute_fn):   l_index(i_compute_fn)+l_count(i_compute_fn) ,fn_atom(i_compute_fn)) * &
!          aux_radial(i_compute_fn)

      enddo

      ! then, the zero functions
      do i_compute_point = 1, n_zero_compute, 1
        i_compute = zero_index_point(i_compute_point)

        wave(i_compute) = 0.0d0

      enddo


    end subroutine evaluate_waves_p2

!******
!---------------------------------------------------------------------
!****s* FHI-aims/mul_vec
!  NAME
!    mul_vec
!  SYNOPSIS

    subroutine mul_vec &
    ( wave, n_mul, ylm, factor &
    )

!  PURPOSE
!  write an explicitly vectorizable multiplication to avoid an index mess
!
!  ARGUMENTS

    implicit none
    integer :: n_mul
    real*8 :: wave(1:n_mul)
    real*8 :: ylm(1:n_mul)
    real*8 :: factor

!  INPUTS
!    o n_mul -- dimensions of the vectors
!    o ylm -- sperical harmonics
!    o factor -- multiplication factor
! 
!  OUTPUT
!    o wave -- basis functions, results from the multiplication
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





    ! begin work

    wave(1:n_mul) = ylm(1:n_mul) * factor

!    write(use_unit,*) ylm
!    write(use_unit,*) factor

    end subroutine mul_vec
!
!******
!---------------------------------------------------------------------
