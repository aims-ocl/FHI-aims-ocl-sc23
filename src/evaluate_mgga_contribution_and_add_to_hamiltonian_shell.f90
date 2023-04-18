
!****** FHI-aims/evaluate_mgga_contribution_and_add_to_hamiltonian_shell
!  NAME
!    evaluate_mgga_contribution_and_add_to_hamiltonian_shell
!  SYNOPSIS
!    Calculates mgga contribution to hamiltonian. Updated to work also with
!    the stress vector routines as implemented in forces_densmat

      subroutine evaluate_mgga_contribution_and_add_to_hamiltonian_shell  &
           ( n_compute_1, n_compute_2, n_points, &
             left_side_of_mgga_dot_product, &
             gradient_basis_wave_store, &
             hamiltonian_shell )
 
!  PURPOSE
!  
!  USES

!      use dimensions
      implicit none

!  ARGUMENTS

      integer :: n_compute_1
      integer :: n_compute_2
      integer :: n_points
      real*8  :: left_side_of_mgga_dot_product(n_compute_1,3*n_points)
      real*8  :: gradient_basis_wave_store(n_compute_2,3*n_points)
      real*8  :: hamiltonian_shell(n_compute_1,n_compute_2) 

!  INPUTS
!   o n_compute -- number of nonzero basis functions in this grid batch.
!   o n_points  -- number of points in this grid patch
!   o left_side_of_dot_product -- Left side of the meta-GGA dot product
!   o gradient_basis_wave_store -- gradient of the basis functions for all points 
!  
!  OUTPUT
!   o hamiltonian_shell -- Updated hamiltonian shell with meta-GGA contributions.
!  
!  AUTHOR
!    AJL, UCL. 2015. Contact: a.logsdail@ucl.ac.uk
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

!     begin work

      ! Implementation using dgemm, calculates dot product in one step
      call dgemm('N','T', n_compute_1, n_compute_2, 3*n_points, 1.0d0, left_side_of_mgga_dot_product, n_compute_1, &
                 gradient_basis_wave_store, n_compute_2, 1.0d0, hamiltonian_shell, n_compute_1)

      end subroutine evaluate_mgga_contribution_and_add_to_hamiltonian_shell 
!---------------------------------------------------------------------

!****** FHI-aims/evaluate_mgga_left_side_of_dot_product
!  NAME
!    evaluate_mgga_left_side_of_dot_product
!  SYNOPSIS

      subroutine evaluate_mgga_left_side_of_dot_product  &
           ( n_compute, n_points, i_point, &
             partition, xc_tau_deriv, &
             gradient_basis_wave, &
             left_side_of_dot_product)

!  PURPOSE
!
!  USES

      implicit none

!  ARGUMENTS

      integer :: n_compute
      integer :: n_points
      integer :: i_point
      real*8  :: partition
      real*8  :: xc_tau_deriv
      real*8  :: gradient_basis_wave(n_compute,3)
      real*8  :: left_side_of_dot_product(n_compute,3*n_points)

!  INPUTS
!   o n_compute -- number of nonzero basis functions in this grid batch.
!   o n_points  -- number of points in this grid patch
!   o i_point   -- current point in this grid patch
!   o partition -- current partition grid weighting
!   o xc_tau_deriv -- Derivative of XC with respect to tau at current point   
!   o gradient_basis_wave -- gradient of the basis functions
!
!  OUTPUT
!   o left_side_of_dot_product -- Left side of the meta-GGA dot product
!
!  AUTHOR
!    AJL, UCL. 2015. Contact: a.logsdail@ucl.ac.uk
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

!     local variables

      integer :: i_coord

!     begin work

      do i_coord = 1, 3, 1
         left_side_of_dot_product(1:n_compute,((i_point-1)*3)+i_coord) = &
            gradient_basis_wave(1:n_compute, i_coord) * partition * xc_tau_deriv
      enddo

    end subroutine evaluate_mgga_left_side_of_dot_product
!---------------------------------------------------------------------

!****** FHI-aims/store_gradient_basis_wave
!  NAME
!    store_gradient_basis_wave
!  SYNOPSIS

      subroutine store_gradient_basis_wave &
           ( n_compute, n_points, i_point, &
           gradient_basis_wave, &
           gradient_basis_wave_store )

!  PURPOSE
!  Store gradient basis wave for meta-GGA calculations
!
!  USES

      implicit none


!  ARGUMENTS

      integer :: n_compute
      integer :: n_points
      integer :: i_point
      real*8  :: gradient_basis_wave(n_compute,3)
      real*8  :: gradient_basis_wave_store(n_compute,3*n_points)

!  INPUTS
!   o n_compute -- number of nonzero basis functions in this grid batch.
!   o n_points  -- number of points in this grid patch
!   o i_point   -- current point in this grid patch
!   o gradient_basis_wave -- gradient of the basis functions
!
!  OUTPUT
!   o gradient_basis_wave_store -- updated array storing all the basis function gradients
!
!  AUTHOR
!    AJL, UCL. 2015. Contact: a.logsdail@ucl.ac.uk
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
!  Local variables

      integer :: i_coord

!     begin work

      do i_coord=1, 3, 1
         gradient_basis_wave_store(1:n_compute, ((i_point-1)*3)+i_coord) = &
            gradient_basis_wave(1:n_compute, i_coord)
      enddo

      end subroutine store_gradient_basis_wave
!---------------------------------------------------------------------
