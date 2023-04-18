!****h* FHI-aims/mpe_constants
!  NAME
!    mpe_constants
!  SYNOPSIS

module mpe_constants

!  PURPOSE
!    This module contains all constants for the MPE continuum solvation module.
!    All values are collected in a (private) derived type and then exposed via
!    a constant public instance.
!  USES
   use types, only: dp
   implicit none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180 (2009), 2175-2196.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   private

   type :: T_MPE_Constants
      ! basis types
      integer :: BASIS_INVAL = -1
      integer :: BASIS_UNDEF = 0
      integer :: BASIS_REG = 1
      integer :: BASIS_IRR = 2
      ! XML logging
      integer :: XML_NOLOG = 0
      integer :: XML_BASIC = 10
      integer :: XML_MEDIUM = 50
      integer :: XML_DETAILED = 100
      ! SLE factorization types
      integer :: FACTZN_INVAL = -1
      integer :: FACTZN_UNDEF = 0
      integer :: FACTZN_QR = 1
      integer :: FACTZN_SVD = 2
      integer :: FACTZN_QRpSVD = 3
      ! non-electrostatic model
      integer :: NONEL_INVAL = -1
      integer :: NONEL_UNDEF = 0
      integer :: NONEL_linOV = 1
   end type
   type(T_MPE_Constants), parameter, public :: MPE_CONST = T_MPE_Constants()

   type :: T_ISC_Constants
      ! cavity types
      integer :: CAVITY_INVAL = -1
      integer :: CAVITY_UNDEF = 0
      integer :: CAVITY_OvlpSph = 1
      integer :: CAVITY_RhoFree = 2
      integer :: CAVITY_RhoMPStat = 3
      integer :: CAVITY_RhoMPDyn = 4
      ! constraint dynamics parameters
      integer :: CDYN_EXIT_NOTCONV = -1
      integer :: CDYN_EXIT_CONV = 0
      integer :: CDYN_EXIT_NOSTEP = 1
      ! sampling types
      integer :: CAV_SAMP_Inval = -1
      integer :: CAV_SAMP_Undef = 0
      integer :: CAV_SAMP_EqSph = 1
      integer :: CAV_SAMP_Even = 2
   end type
   type(T_ISC_Constants), parameter, public :: ISC_CONST = T_ISC_Constants()

end module mpe_constants

