!****h* FHI-aims/inner_product
!  NAME
!    inner_product
!  SYNOPSIS
MODULE inner_product
!  PURPOSE
!    This module provides the inner product functions for (overlap) matrix
!    based inner products.
!  USES
  USE runtime_choices

  implicit none
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

  INTERFACE S_inner_product
     MODULE PROCEDURE S_inner_product_real, S_inner_product_complex
  END INTERFACE

  PUBLIC :: S_inner_product

  real*8, external :: ddot
  complex*16, external :: zdotc

!******
CONTAINS
!****s* FHI-aims/S_inner_product_real
!  NAME
!    S_inner_product_real
!  SYNOPSIS
  REAL*8 FUNCTION S_inner_product_real(vector1, vector2, matrix, length) RESULT(value)
!  PURPOSE
!    Real-valued (overlap) matrix based inner product.
!  USES
!  ARGUMENTS 
    integer :: length
    real*8, dimension(length) :: vector1, vector2
    real*8, dimension(length*(length+1)/2) :: matrix
!  INPUTS
!    o length -- length of the arrays
!    o vector1 -- first vector of the inner product
!    o vector2 -- second vector of the inner product
!    o matrix -- matrix defining the inner product in packed storage
!  OUTPUT
!    o value -- value of the inner product
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    ! locals
    real*8, dimension(length) :: tempvec

    call dspmv('U', length, 1.0d0, matrix, vector2, 1, 0.0d0, tempvec, 1)
    value = ddot( length, vector1, 1, tempvec, 1)   

  END FUNCTION S_inner_product_real
!******
!****s* FHI-aims/S_inner_product_complex
!  NAME
!    S_inner_product_complex
!  SYNOPSIS
  REAL*8 FUNCTION S_inner_product_complex(vector1, vector2, matrix, length) RESULT(value)
!  PURPOSE
!    Real-valued (overlap) matrix based inner product with a complex matrix
!  USES
!  ARGUMENTS
    integer :: length
    complex*16, dimension(length) :: vector1, vector2
    complex*16, dimension(length*(length+1)/2) :: matrix
!  INPUTS
!    o length -- length of the arrays
!    o vector1 -- first vector of the inner product
!    o vector2 -- second vector of the inner product
!    o matrix -- matrix defining the inner product in packed storage
!  OUTPUT
!    o value -- real part of the value of the inner product
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    ! locals
    complex*16, dimension(length) :: tempvec
    complex*16 :: value_complex
    integer :: i_index

    call zhpmv('U', length, (1.0d0, 0.0d0), matrix, vector2, 1, (0.0d0, 0.0d0), tempvec, 1)
    
    value_complex = zdotc( length, vector1, 1, tempvec, 1)
    value = dble(value_complex) 

  END FUNCTION S_inner_product_complex
!******
end MODULE inner_product
