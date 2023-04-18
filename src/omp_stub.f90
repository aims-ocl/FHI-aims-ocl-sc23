!****s* FHI-aims/OMP_GET_NUM_THREADS
!  NAME
!   OMP_GET_NUM_THREADS
!  SYNOPSIS
integer function OMP_GET_NUM_THREADS( )
!  PURPOSE
!    This is a omp stub that is only compiled into the code if 
!    the OpenMP support is not available. Its only purpose is to make the
!    (in that case unused) omp parts of the code compile. 
!    please refer to the OpenMP  to infer about the purpose
!    and use of the actual routine corresponding to this stub!
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
!    Release version, FHI-aims (2014).
!  SOURCE

  use localorb_io
  
  implicit none
  
  character*100 :: info_str
  integer :: num_thr

  num_thr = 0

  write(info_str,'(1X,A,A)') '* You have called a OMP stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str, use_unit, '(A)')
  stop
  
end function OMP_GET_NUM_THREADS
!****s* FHI-aims/OMP_GET_WTIME
!  NAME
!   OMP_GET_WTIME
!  SYNOPSIS
DOUBLE PRECISION function OMP_GET_WTIME( )
!  PURPOSE
!    This is a OMP stub that is only compiled into the code if 
!    the OpenMP support is not available. Its only purpose is to make the
!    (in that case unused) omp parts of the code compile. 
!    please refer to the OpenMP  to infer about the purpose
!    and use of the actual routine corresponding to this stub!
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
!    Release version, FHI-aims (2014).
!  SOURCE

  use localorb_io
  
  implicit none
  
  character*100 :: info_str
  integer :: num_thr

  num_thr = 0

  write(info_str,'(1X,A,A)') '* You have called a OMP stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str, use_unit, '(A)')
  stop
  
end function OMP_GET_WTIME
