!****s* FHI-aims/BLACS_Gridinit
!  NAME
!   BLACS_Gridinit
!  SYNOPSIS
subroutine BLACS_Gridinit( comm, uplo, nprow, npcol )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE
  use localorb_io

  implicit none

  character*100 :: info_str
  integer :: comm
  character :: uplo
  integer :: nprow, npcol

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine BLACS_Gridinit
!******

!****s* FHI-aims/BLACS_GET
!  NAME
!   BLACS_GET
!  SYNOPSIS
subroutine BLACS_GET( int1, int2, my_blacs_ctxt_aux)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE
  use localorb_io
  implicit none
  integer:: int1, int2, my_blacs_ctxt_aux
  character*100 :: info_str


  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine BLACS_GET
!******

!****s* FHI-aims/blacs_freebuff
!  NAME
!   blacs_freebuff
!  SYNOPSIS
subroutine blacs_freebuff(icontxt, wait)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
!    and use of the actual routine corresponding to this stub!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE
  use localorb_io
  implicit none

  integer :: icontxt, wait
  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine blacs_freebuff
!******



!****s* FHI-aims/pdgetrf
!  NAME
!   pdgetrf
!  SYNOPSIS
subroutine pdgetrf( n_basbas, n_basbas2, v_times_polar, int1, int2, &
                           sc_desc, ipiv_scal, info)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE
  use localorb_io
  implicit none
  integer:: int1, int2, n_basbas, n_basbas2,  ipiv_scal, sc_desc, info
  real*8:: v_times_polar
  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdgetrf
!******

!****s* FHI-aims/BLACS_Pinfo
!  NAME
!   BLACS_Pinfo
!  SYNOPSIS
subroutine BLACS_Pinfo( mypnum, nprocs )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE
  implicit none
  integer :: mypnum, nprocs
  mypnum = 0
  nprocs = -5
end subroutine BLACS_Pinfo
!******

!****s* FHI-aims/BLACS_Gridinfo
!  NAME
!   BLACS_Gridinfo
!  SYNOPSIS
subroutine BLACS_Gridinfo( comm, uplo, nprow, npcol, myprow, mypcol )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  integer :: comm
  character :: uplo
  integer :: nprow, npcol
  integer :: myprow, mypcol

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine BLACS_Gridinfo
!******

!****s* FHI-aims/BLACS_Gridexit
!  NAME
!   BLACS_Gridexit
!  SYNOPSIS
subroutine BLACS_Gridexit( comm )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  integer :: comm
  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine BLACS_Gridexit
!******

!----------------------------------------------------------------------------
!****s* FHI-aims/BLACS_Exit
!  NAME
!   BLACS_Exit
!  SYNOPSIS
subroutine BLACS_Exit( value )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  integer :: value
  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine BLACS_Exit
!******
!----------------------------------------------------------------------------
!****s* FHI-aims/info2gl
!  NAME
!   info2gl
!  SYNOPSIS
subroutine infog2l( row, col, desc, nprow, npcol, myprow, mypcol, &
     localrow, localcol, prow, pcol )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  integer, dimension(9) :: desc
  integer :: row, col, nprow, npcol, myprow, mypcol, localrow, localcol
  integer :: prow, pcol

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine infog2l
!******

!****s* FHI-aims/infog1l
!  NAME
!   infog1l
!  SYNOPSIS
subroutine infog1l( gindx, nb, nprocs, myproc, isrcproc, &
     lindx, rocsrc )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
!    and use of the actual routine corresponding to this stub!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  integer :: gindx, nb, nprocs, myproc, isrcproc, lindx, rocsrc

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine infog1l
!******

!****s* FHI-aims/iceil
!  NAME
!   iceil
!  SYNOPSIS
integer function iceil( inum, idenom )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
!    and use of the actual routine corresponding to this stub!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  integer :: inum, idenom

  iceil = 0

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end function iceil
!******

!****s* FHI-aims/pilaenv
!  NAME
!   pilaenv
!  SYNOPSIS
integer function pilaenv(ictxt, prec)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
!    and use of the actual routine corresponding to this stub!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE
  use localorb_io
  implicit none

  integer :: ictxt
  character :: prec
  character*100 :: info_str

  pilaenv = 0

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end function pilaenv
!******

!****s* FHI-aims/pjlaenv
!  NAME
!   pjlaenv
!  SYNOPSIS
integer function pjlaenv( comm, i, string, string2, j1, j2, j3, j4 )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  character :: string, string2
  integer :: comm
  integer :: i, j1, j2, j3, j4

  pjlaenv = 0

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end function pjlaenv
!******

!****s* FHI-aims/descinit
!  NAME
!   descinit
!  SYNOPSIS
subroutine descinit( desc, n, m, nb, mb, rsrc, csrc, comm, mxld, info )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  integer, dimension(9) :: desc
  integer :: n, m, nb, mb, rsrc, csrc
  integer :: comm, mxld, info

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine descinit
!******

!****s* FHI-aims/numroc
!  NAME
!   numroc
!  SYNOPSIS
integer function numroc( n, nb, mycol, csrc, npcol )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  integer :: n, nb, mycol, csrc, npcol

  numroc = 0

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end function numroc
!******

!****s* FHI-aims/CHK1MAT
!  NAME
!   CHK1MAT
!  SYNOPSIS
SUBROUTINE CHK1MAT( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA, &
     DESCAPOS0, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  INTEGER :: DESCAPOS0, IA, INFO, JA, MA, MAPOS0, NA, NAPOS0
  integer, dimension(9) :: desca

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end SUBROUTINE CHK1MAT
!******

!****s* FHI-aims/PCHK1MAT
!  NAME
!   PCHK1MAT
!  SYNOPSIS
SUBROUTINE PCHK1MAT( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA, &
     DESCAPOS0, NEXTRA, EX, EXPOS, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  INTEGER :: DESCAPOS0, IA, INFO, JA, MA, MAPOS0, NA, NAPOS0, NEXTRA
  INTEGER :: EX( NEXTRA ), EXPOS( NEXTRA )

  integer, dimension(9) :: desca

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

END SUBROUTINE PCHK1MAT
!******

!****s* FHI-aims/PCHK2MAT
!  NAME
!   PCHK2MAT
!  SYNOPSIS
SUBROUTINE PCHK2MAT( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA, &
     DESCAPOS0, MB, MBPOS0, NB, NBPOS0, IB, JB, &
     DESCB, DESCBPOS0, NEXTRA, EX, EXPOS, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  INTEGER :: DESCAPOS0, DESCBPOS0, IA, IB, INFO, JA, JB, MA
  INTEGER :: MAPOS0, MB, MBPOS0, NA, NAPOS0, NB, NBPOS0
  INTEGER :: NEXTRA

  INTEGER :: DESCA(9), DESCB( 8 ), EX( NEXTRA ), EXPOS( NEXTRA )

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

END SUBROUTINE PCHK2MAT
!******

!****s* FHI-aims/INDXL2G
!  NAME
!   INDXL2G
!  SYNOPSIS
INTEGER FUNCTION INDXL2G( INDXLOC, NB, IPROC, ISRCPROC, NPROCS )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  INTEGER :: INDXLOC, IPROC, ISRCPROC, NB, NPROCS

  character*100 :: info_str

  INDXL2G = 0

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

END FUNCTION INDXL2G
!******

!****s* FHI-aims/INDXG2L
!  NAME
!   INDXG2L
!  SYNOPSIS
INTEGER FUNCTION INDXG2L( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  INTEGER :: INDXGLOB, IPROC, ISRCPROC, NB, NPROCS

  character*100 :: info_str

  indxg2l = 0

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end FUNCTION INDXG2L
!******

!****s* FHI-aims/INDXG2P
!  NAME
!   INDXG2P
!  SYNOPSIS
INTEGER FUNCTION INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  INTEGER :: INDXGLOB, IPROC, ISRCPROC, NB, NPROCS

  character*100 :: info_str

  indxg2p = 0

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end FUNCTION INDXG2P
!******


!****s* FHI-aims/PXERBLA
!  NAME
!   PXERBLA
!  SYNOPSIS
SUBROUTINE PXERBLA( ICTXT, SRNAME, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  INTEGER :: ICTXT, INFO
  CHARACTER :: SRNAME

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

END SUBROUTINE PXERBLA
!******

!------------------------------------------------------------------------------
!****s* FHI-aims/pdelget
!  NAME
!   pdelget
!  SYNOPSIS
subroutine pdelget( scope, top, value, array, row, col, desc )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  real*8 :: array
  integer, dimension(9) :: desc
  integer :: row, col
  real*8 :: value
  character :: scope, top

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdelget
!******

!****s* FHI-aims/pdelset
!  NAME
!   pdelset
!  SYNOPSIS
subroutine pdelset( array, row, col, desc, value )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  real*8 :: array
  integer, dimension(9) :: desc
  integer :: row, col
  real*8 :: value

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdelset
!******


!****s* FHI-aims/pdscal
!  NAME
!   pdscal
!  SYNOPSIS
subroutine pdscal( n, alpha, array, row, col, desc, m )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  real*8 :: array
  integer, dimension(9) :: desc
  integer :: row, col, m, n
  real*8 :: alpha

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdscal
!******

!****s* FHI-aims/pdrscl
!  NAME
!   pdrscl
!  SYNOPSIS
subroutine pdrscl( n, alpha, array, row, col, desc, m )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  real*8 :: array
  integer, dimension(9) :: desc
  integer :: row, col, m, n
  real*8 :: alpha

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdrscl
!******

!****s* FHI-aims/pdsyevx
!  NAME
!   pdsyevx
!  SYNOPSIS
subroutine pdsyevx( jobz, range, uplo, n, a, ia, ja, desca, &
     vl, vu, il, iu, abstol, m, nz, w, orfac, z, &
     iz, jz, descz, work, lwork, iwork, liwork, ifail, iclustr, &
     gap, info )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  real*8 :: a, w, z
  real*8 :: gap, work
  integer :: iwork, ifail, iclustr
  real*8 :: vl, vu, abstol, orfac, info
  integer, dimension(9) :: desca, descz
  integer :: n, ia, ja, il, iu, m, nz, iz, jz, lwork, liwork
  character :: jobz, range, uplo

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdsyevx
!******

!****s* FHI-aims/pdsgvx
!  NAME
!   pdsgvx
!  SYNOPSIS
subroutine pdsygvx( jobz, range, uplo, n, a, ia, ja, desca, &
     b, ib, jb, descb, vl, vu, il, iu, abstol, m, nz, w, orfac, z, &
     iz, jz, descz, work, lwork, iwork, liwork, ifail, iclustr, &
     gap, info )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  real*8 :: a, b, w, z
  real*8 :: gap, work
  integer :: iwork, ifail, iclustr
  real*8 :: vl, vu, abstol, orfac, info
  integer, dimension(9) :: desca, descb, descz
  integer :: n, ia, ja, ib, jb, il, iu, m, nz, iz, jz, lwork, liwork
  character :: jobz, range, uplo

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdsygvx
!******

!****s* FHI-aims/pdsyrk
!  NAME
!   pdsyrk
!  SYNOPSIS
subroutine pdsyrk( uplo, tran, n, k, alpha, array, row, col, desc, beta, &
     array2, row2, col2, desc2)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  character :: uplo, tran
  real*8 :: array, array2
  integer, dimension(9) :: desc, desc2
  integer :: row, col, row2, col2, n, k
  real*8 :: alpha, beta

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdsyrk
!******

!****s* FHI-aims/pdtran
!  NAME
!   pdtran
!  SYNOPSIS
subroutine pdtran( n, m, alpha, array, row, col, desc, beta, &
     array2, row2, col2, desc2)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  real*8 :: array, array2
  integer, dimension(9) :: desc, desc2
  integer :: row, col, row2, col2, m, n
  real*8 :: alpha, beta

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdtran
!******

!****s* FHI-aims/pdsyngst
!  NAME
!   pdsyngst
!  SYNOPSIS
subroutine pdsyngst(IBTYPE, UPLO, N, A, IA, JA, DESCA, B, IB, JB, &
     DESCB, SCALE, WORK, LWORK, INFO)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  integer :: IA, IB, IBTYPE, INFO, JA, JB, LWORK, N
  character :: UPLO
  real*8 :: SCALE, A, B, WORK
  integer, dimension(9) :: desca, descb

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdsyngst
!******

!****s* FHI-aims/PDOTRF
!  NAME
!   PDPOTRF
!  SYNOPSIS
SUBROUTINE PDPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: UPLO
  INTEGER :: IA, INFO, JA, N

  integer, dimension(9) :: desca
  real*8 :: A

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end SUBROUTINE PDPOTRF
!******

!****s* FHI-aims/pdlacpy
!  NAME
!   pdlacpy
!  SYNOPSIS
subroutine pdlacpy( uplo, m, n, a, ia, ja, desca, b, ib, jb, descb )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
!    and use of the actual routine corresponding to this stub!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

  use localorb_io

  implicit none

  character :: uplo
  integer :: m, n, ia, ja, desca(*), ib, jb, descb(*)
  real*8 :: a(*), b(*)

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdlacpy
!******

!****s* FHI-aims/pdlatra
!  NAME
!   pdlatra
!  SYNOPSIS
real*8 function pdlatra( n, a, ia, ja, desca )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
!    and use of the actual routine corresponding to this stub!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

  use localorb_io

  implicit none

  integer :: n, ia, ja, desca(*)
  real*8 :: a(*)

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end function pdlatra
!******

!****s* FHI-aims/PDLAMCH
!  NAME
!   PDLAMCH
!  SYNOPSIS
real*8 FUNCTION PDLAMCH( ICTXT, CMACH )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character :: CMACH
  integer :: ICTXT

  character*100 :: info_str

  pdlamch = 0.0d0

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end FUNCTION PDLAMCH
!******

!****s* FHI-aims/PSYEVD
!  NAME
!   PSYEVD
!  SYNOPSIS
SUBROUTINE PDSYEVD( JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, &
     DESCZ, WORK, LWORK, IWORK, LIWORK, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: JOBZ, UPLO
  INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LWORK, N, IWORK

  integer, dimension(9) :: desca, descz
  real*8 :: A, W, WORK, Z
  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end SUBROUTINE PDSYEVD
!******

!****s* FHI-aims/PDLANSY
!  NAME
!   PDLANSY
!  SYNOPSIS
real*8 FUNCTION PDLANSY( NORM, UPLO, N, A, IA, JA, &
     DESCA, WORK )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: NORM, UPLO
  INTEGER :: IA, JA, N

  integer, dimension(9) :: desca

  real*8 :: A( * ), WORK( * )
  character*100 :: info_str

  pdlansy = 0.0d0

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end FUNCTION PDLANSY
!******

!****s* FHI-aims/PDSTEC
!  NAME
!   PDSTEC
!  SYNOPSIS
SUBROUTINE PDSTEDC( COMPZ, N, D, E, Q, IQ, JQ, DESCQ, WORK, LWORK, &
     IWORK, LIWORK, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: COMPZ
  INTEGER :: INFO, IQ, JQ, LIWORK, LWORK, N, IWORK

  integer, dimension(9) :: descq

  real*8 :: D, E, Q, WORK

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end SUBROUTINE PDSTEDC
!******

!****s* FHI-aims/PDLASET
!  NAME
!   PDLASET
!  SYNOPSIS
SUBROUTINE PDLASET( UPLO, M, N, ALPHA, BETA, A, IA, JA, DESCA )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: UPLO
  INTEGER :: IA, JA, M, N
  real*8 :: ALPHA, BETA, A

  integer, dimension(9) :: desca

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

END SUBROUTINE PDLASET
!******

!****s* FHI-aims/PDLARED1d
!  NAME
!   PDLARED1D
!  SYNOPSIS
SUBROUTINE PDLARED1D( N, IA, JA, DESC, BYCOL, BYALL, WORK, LWORK )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  INTEGER :: IA, JA, LWORK, N

  real*8 :: BYALL, BYCOL
  real*8 :: WORK( LWORK )

  integer, dimension(9) :: desc

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end SUBROUTINE PDLARED1D
!******

!****s* FHI-aims/PDORMTR
!  NAME
!   PDORMTR
!  SYNOPSIS
SUBROUTINE PDORMTR( SIDE, UPLO, TRANS, M, N, A, IA, JA, DESCA, &
     TAU, C, IC, JC, DESCC, WORK, LWORK, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: SIDE, TRANS, UPLO
  INTEGER :: IA, IC, INFO, JA, JC, LWORK, M, N
  integer, dimension(9) :: desca, descc
  real*8 :: A, C, TAU, WORK

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end SUBROUTINE PDORMTR
!******

!****s* FHI-aims/PDSYTRD
!  NAME
!   PDSYTRD
!  SYNOPSIS
SUBROUTINE PDSYTRD( UPLO, N, A, IA, JA, DESCA, D, E, TAU, WORK, &
     LWORK, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: UPLO
  INTEGER :: IA, INFO, JA, LWORK, N

  integer, dimension(9) :: desca

  real*8 :: D, E
  real*8 :: A, TAU, WORK

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end SUBROUTINE PDSYTRD
!******

!****s* FHI-aims/pdtrsm
!  NAME
!   pdtrsm
!  SYNOPSIS
subroutine pdtrsm( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, &
     A, IA, JA, DESCA, B, IB, JB, DESCB )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character :: DIAG, SIDE, TRANS, UPLO
  integer :: IA, IB, JA, JB, M, N
  real*8 :: ALPHA, A, B

  integer, dimension(9) :: desca, descb

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdtrsm
!******

!****s* FHI-aims/pdtrtrs
!  NAME
!   pdtrtrs
!  SYNOPSIS
subroutine pdtrtrs( uplo, trans, diag, n, nrhs, a, ia, ja, desca, &
               b, ib, jb, descb, info )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
!    and use of the actual routine corresponding to this stub!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

  use localorb_io

  implicit none

  character :: uplo, trans, diag
  integer :: n, nrhs, ia, ja, desca(*), ib, jb, descb(*), info
  real*8 :: a(*), b(*)

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdtrtrs
!******

!****s* FHI-aims/pdsymv
!  NAME
!   pdsymv
!  SYNOPSIS
subroutine pdsymv( UPLO, N, ALPHA, A, IA, JA, DESCA, X, IX, JX, DESCX, &
     INCX, BETA, Y, IY, JY, DESCY, INCY )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character ::  UPLO
  integer :: IA, INCX, INCY, IX, IY, JA, JX, JY, N
  real*8 :: ALPHA, BETA, A, X, Y

  integer, dimension(9) :: desca, descx, descy

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdsymv
!******

!****s* FHI-aims/pddot
!  NAME
!   pddot
!  SYNOPSIS
subroutine pddot( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  integer :: INCX, INCY, IX, IY, JX, JY, N
  real*8 :: DOT
  real*8 :: X, Y

  integer, dimension(9) :: descx, descy

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pddot
!******

!****s* FHI-aims/PDLASCL
!  NAME
!   PDLASCL
!  SYNOPSIS
SUBROUTINE PDLASCL( TYPE, CFROM, CTO, M, N, A, IA, JA, DESCA, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: TYPE
  INTEGER :: IA, INFO, JA, M, N
  real*8 :: CFROM, CTO

  real*8 :: A

  integer, dimension(9) :: desca

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end SUBROUTINE PDLASCL
!******

!****s* FHI-aims/pdgemm
!  NAME
!   pdgemm
!  SYNOPSIS
subroutine pdgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, &
     B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  character :: TRANSA, TRANSB
  integer :: M, N, K, IA, JA, IB, JB, IC, JC
  real*8 :: ALPHA, BETA, A, B , C

  integer, dimension(9) :: DESCA, DESCB, DESCC

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdgemm
!******


!****s* FHI-aims/pdgemv
!  NAME
!   pdgemv
!  SYNOPSIS
subroutine pdgemv( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, &
     B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  character :: TRANSA, TRANSB
  integer :: M, N, K, IA, JA, IB, JB, IC, JC
  real*8 :: ALPHA, BETA, A, B , C

  integer, dimension(9) :: DESCA, DESCB, DESCC

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdgemv
!******

!****s* FHI-aims/pzgemv
!  NAME
!   pdzgemv
!  SYNOPSIS
subroutine pzgemv( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, &
     B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  character :: TRANSA, TRANSB
  integer :: M, N, K, IA, JA, IB, JB, IC, JC
  real*8 :: ALPHA, BETA, A, B , C

  integer, dimension(9) :: DESCA, DESCB, DESCC

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzgemv
!******




!****s* FHI-aims/PDSYMM
!  NAME
!   PDSYMM
!  SYNOPSIS
subroutine PDSYMM( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, &
     B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

!('L','U',n_basis,n_states,1.d0,ham(1,1,i_spin),1,1,sc_desc, &
!                            eigenvec(1,1,i_spin),1,1,sc_desc,0.d0,tmp,1,1,sc_desc)

!subroutine pdgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, &
!     B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )

  use localorb_io

  character :: TRANSA, TRANSB
  integer :: M, N, K, IA, JA, IB, JB, IC, JC
  real*8 :: ALPHA, BETA, A, B , C

  integer, dimension(9) :: DESCA, DESCB, DESCC

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine PDSYMM
!******

!****s* FHI-aims/PDAXPY
!  NAME
!   PDAXPY
!  SYNOPSIS
subroutine PDAXPY( N, ALPHA, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  integer :: N, IX, JX, INCX, IY, JY, INCY
  real*8 :: ALPHA, X, Y

  integer, dimension(9) :: DESCX, DESCY

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine PDAXPY
!******

!****s* FHI-aims/PDASUM
!  NAME
!   PDASUM
!  SYNOPSIS
subroutine PDASUM( N, ASUM, X, IX, JX, DESCX, INCX )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io
  implicit none

  integer :: N, IX, JX, INCX
  real*8 :: ASUM, X(*)

  integer, dimension(9) :: DESCX

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine PDASUM
!******

!****s* FHI-aims/PDNRM2
!  NAME
!   PDNRM2
!  SYNOPSIS
subroutine PDNRM2( N, NORM2, X, IX, JX, DESCX, INCX )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io
  implicit none

  integer :: N, IX, JX, INCX
  real*8 :: NORM2, X(*)

  integer, dimension(9) :: DESCX

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine PDNRM2
!******
!****s* FHI-aims/PDCOPY
!  NAME
!   PDCOPY
!  SYNOPSIS
subroutine PDCOPY( N, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  integer :: N, IX, JX, INCX, IY, JY, INCY
  real*8 :: X, Y

  integer, dimension(9) :: DESCX, DESCY

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine PDCOPY
!******
!****s* FHI-aims/PZCOPY
!  NAME
!   PZCOPY
!  SYNOPSIS
subroutine PZCOPY( N, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  integer :: N, IX, JX, INCX, IY, JY, INCY
  complex*16 :: X, Y

  integer, dimension(9) :: DESCX, DESCY

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine PZCOPY
!******
!--------------------------------------------------------------------------------
!****s* FHI-aims/pzelget
!  NAME
!   pzelget
!  SYNOPSIS
subroutine pzelget( scope, top, value, array, row, col, desc )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  complex*16 :: array
  integer, dimension(9) :: desc
  integer :: row, col
  complex*16 :: value
  character :: scope, top

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzelget
!******

!****s* FHI-aims/pzelset
!  NAME
!   pzelset
!  SYNOPSIS
subroutine pzelset( array, row, col, desc, value )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  complex*16 :: array
  integer, dimension(9) :: desc
  integer :: row, col
  complex*16 :: value

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzelset
!******

!****s* FHI-aims/pzhegvx
!  NAME
!   pzhegvx
!  SYNOPSIS
subroutine pzhegvx( ibtype, jobz, range, uplo, n, a, ia, ja, desca, &
     b, ib, jb, descb, vl, vu, il, iu, abstol, m, nz, w, orfac, z, &
     iz, jz, descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, &
     gap, info )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  complex*16 :: a, b, w, z
  complex*16 :: work
  real*8 :: gap, rwork
  integer :: ibtype, iwork, ifail, iclustr
  real*8 :: vl, vu, abstol, orfac, info
  integer, dimension(9) :: desca, descb, descz
  integer :: n, ia, ja, ib, jb, il, iu, m, nz, iz, jz, lwork, lrwork, liwork
  character :: jobz, range, uplo

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzhegvx
!******

!****s* FHI-aims/pzherk
!  NAME
!   pzherk
!  SYNOPSIS
subroutine pzherk( uplo, tran, n, k, alpha, array, row, col, desc, beta, &
     array2, row2, col2, desc2)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  character :: uplo, tran
  complex*16 :: array, array2
  integer, dimension(9) :: desc, desc2
  integer :: row, col, row2, col2, n, k
  complex*16 :: alpha, beta

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzherk
!******

!****s* FHI-aims/pzscal
!  NAME
!   pzscal
!  SYNOPSIS
subroutine pzscal( n, alpha, array, row, col, desc, m )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  complex*16 :: array
  integer, dimension(9) :: desc
  integer :: row, col, m, n
  complex*16 :: alpha

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzscal
!******

!****s* FHI-aims/pztranc
!  NAME
!   pztranc
!  SYNOPSIS
subroutine pztranc( n, m, alpha, array, row, col, desc, beta, &
     array2, row2, col2, desc2)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  complex*16 :: array, array2
  integer, dimension(9) :: desc, desc2
  integer :: row, col, row2, col2, m, n
  complex*16 :: alpha, beta

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pztranc
!******

!****s* FHI-aims/pzhengst
!  NAME
!   pzhengst
!  SYNOPSIS
subroutine pzhengst(IBTYPE, UPLO, N, A, IA, JA, DESCA, B, IB, JB, &
     DESCB, SCALE, WORK, LWORK, INFO)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  integer :: IA, IB, IBTYPE, INFO, JA, JB, LWORK, N
  character :: UPLO
  real*8 :: SCALE
  complex*16 :: A, B, WORK
  integer, dimension(9) :: desca, descb

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzhengst
!******

!****s* FHI-aims/PZUNMTR
!  NAME
!   PZUNMTR
!  SYNOPSIS
SUBROUTINE PZUNMTR( SIDE, UPLO, TRANS, M, N, A, IA, JA, DESCA, &
     TAU, C, IC, JC, DESCC, WORK, LWORK, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: SIDE, TRANS, UPLO
  INTEGER :: IA, IC, INFO, JA, JC, LWORK, M, N
  integer, dimension(9) :: desca, descc
  COMPLEX*16 :: A, C, TAU, WORK

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end SUBROUTINE PZUNMTR
!******

!****s* FHI-aims/PZPOTRF
!  NAME
!   PZPOTRF
!  SYNOPSIS
SUBROUTINE PZPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: UPLO
  INTEGER :: IA, INFO, JA, N

  integer, dimension(9) :: desca
  complex*16 :: A

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end SUBROUTINE PZPOTRF
!******

!****s* FHI-aims/PZLASET
!  NAME
!   PZLASET
!  SYNOPSIS
SUBROUTINE PZLASET( UPLO, M, N, ALPHA, BETA, A, IA, JA, DESCA )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: UPLO
  INTEGER :: IA, JA, M, N
  complex*16 :: ALPHA, BETA, A

  integer, dimension(9) :: desca

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

END SUBROUTINE PZLASET
!******

!****s* FHI-aims/PSLANHE
!  NAME
!   PZLANHE
!  SYNOPSIS
real*8 FUNCTION PZLANHE( NORM, UPLO, N, A, IA, JA, DESCA, WORK )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: NORM, UPLO
  INTEGER :: IA, JA, N

  integer, dimension(9) :: desca

  real*8 :: WORK
  COMPLEX*16 :: A

  character*100 :: info_str

  pzlanhe = 0.0d0

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

END FUNCTION PZLANHE
!******

!****s* FHI-aims/PZHETRD
!  NAME
!   PZHETRD
!  SYNOPSIS
SUBROUTINE PZHETRD( UPLO, N, A, IA, JA, DESCA, D, E, TAU, WORK, &
     LWORK, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: UPLO
  INTEGER :: IA, INFO, JA, LWORK, N

  integer, dimension(9) :: desca

  real*8 :: D, E
  COMPLEX*16 :: A, TAU, WORK

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end SUBROUTINE PZHETRD
!******

!****s* FHI-aims/pztrsm
!  NAME
!   pztrsm
!  SYNOPSIS
subroutine pztrsm( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, &
     A, IA, JA, DESCA, B, IB, JB, DESCB )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character :: DIAG, SIDE, TRANS, UPLO
  integer :: IA, IB, JA, JB, M, N
  complex*16 :: ALPHA, A, B

  integer, dimension(9) :: desca, descb

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pztrsm
!******

!****s* FHI-aims/pzhemv
!  NAME
!   pzhemv
!  SYNOPSIS
subroutine pzhemv( UPLO, N, ALPHA, A, IA, JA, DESCA, X, IX, JX, DESCX, &
     INCX, BETA, Y, IY, JY, DESCY, INCY )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character ::  UPLO
  integer :: IA, INCX, INCY, IX, IY, JA, JX, JY, N
  complex*16 :: ALPHA, BETA, A, X, Y

  integer, dimension(9) :: desca, descx, descy

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzhemv
!******

!****s* FHI-aims/pzdotc
!  NAME
!   pzdotc
!  SYNOPSIS
subroutine pzdotc( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  integer :: INCX, INCY, IX, IY, JX, JY, N
  complex*16 :: DOT
  complex*16 :: X, Y

  integer, dimension(9) :: descx, descy

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzdotc
!******

!****s* FHI-aims/PLZASCL
!  NAME
!   PSLASCL
!  SYNOPSIS
SUBROUTINE PZLASCL( TYPE, CFROM, CTO, M, N, A, IA, JA, DESCA, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  CHARACTER :: TYPE
  INTEGER :: IA, INFO, JA, M, N
  real*8 :: CFROM, CTO

  COMPLEX*16 :: A

  integer, dimension(9) :: desca

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end SUBROUTINE PZLASCL
!******

!****s* FHI-aims/pzgemm
!  NAME
!   pzgemm
!  SYNOPSIS
subroutine pzgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, &
     B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  character :: TRANSA, TRANSB
  integer :: M, N, K, IA, JA, IB, JB, IC, JC
  complex*16 :: ALPHA, BETA, A, B , C

  integer, dimension(9) :: DESCA, DESCB, DESCC

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzgemm
!******

!****s* FHI-aims/PZHEMM
!  NAME
!   PZHEMM
!  SYNOPSIS
subroutine PZHEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, &
     B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

!('L','U',n_basis,n_states,(1.d0,0.d0),ham_complex(1,1,i_spin),1,1,sc_desc, &
!                            eigenvec_complex(1,1,i_spin),1,1,sc_desc,(0.d0,0.d0),tmp_complex,1,1,sc_desc)

!subroutine pzgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, &
!     B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )

  use localorb_io

  character :: TRANSA, TRANSB
  integer :: M, N, K, IA, JA, IB, JB, IC, JC
  complex*16 :: ALPHA, BETA, A, B , C

  integer, dimension(9) :: DESCA, DESCB, DESCC

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine PZHEMM
!******

!****s* FHI-aims/PZHEEVD
!  NAME
!   PSHEEVD
!  SYNOPSIS
SUBROUTINE PZHEEVD( JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ, &
     DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, &
     LIWORK, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  CHARACTER :: JOBZ, UPLO
  INTEGER :: IA, INFO, IZ, JA, JZ, LIWORK, LRWORK, LWORK, N

  INTEGER, dimension(9) :: DESCA, DESCZ
  INTEGER :: IWORK
  real*8 :: RWORK, W
  COMPLEX*16 :: A, WORK, Z

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

END SUBROUTINE PZHEEVD
!******

!****s* FHI-aims/PDGEMR2D
!  NAME
!   PDGEMR2D
!  SYNOPSIS
SUBROUTINE PDGEMR2D( M, N, A, IA, JA, ADESC, B, IB, JB, BDESC, CTXT)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  INTEGER :: M, N, IA, JA, IB, JB, CTXT
  INTEGER, dimension(9) :: ADESC, BDESC
  REAL*8 :: A, B
  CHARACTER*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

END SUBROUTINE PDGEMR2D
!******

!****s* FHI-aims/PZGEMR2D
!  NAME
!   PZGEMR2D
!  SYNOPSIS
SUBROUTINE PZGEMR2D( M, N, A, IA, JA, ADESC, B, IB, JB, BDESC, CTXT)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  INTEGER :: M, N, IA, JA, IB, JB, CTXT
  INTEGER, dimension(9) :: ADESC, BDESC
  COMPLEX*16 :: A, B
  CHARACTER*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

END SUBROUTINE PZGEMR2D
!******

!****s* FHI-aims/PDGEQRF
!  NAME
!   PDGEQRF
!  SYNOPSIS
SUBROUTINE PDGEQRF( M, N, A, IA, JA, DESCA, TAU, WORK, LWORK, INFO)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  INTEGER :: M, N, IA, JA, LWORK, INFO
  INTEGER :: DESCA(*)
  REAL*8 :: A(*), TAU(*), WORK(*)
  CHARACTER*100 :: info_str

  write(info_str,'(1X,A,A)') &
      "* You have called the ScaLAPACK stub routine PDGEQRF. ", &
      "Please check control.in and the libraries you linked."
  call localorb_info(info_str,use_unit,'(A)')
  stop

END SUBROUTINE PDGEQRF
!******
!****s* FHI-aims/PDORMQR
!  NAME
!   PDORMQR
!  SYNOPSIS
SUBROUTINE PDORMQR(SIDE, TRANS, M, N, K, A, IA, JA, DESCA, TAU, C, IC, JC, DESCC, WORK, LWORK, INFO)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  INTEGER :: M, N, K, IA, JA, IC, JC, LWORK, INFO
  INTEGER :: DESCA(*), DESCC(*)
  REAL*8 :: A(*), C(*), TAU(*), WORK(*)
  CHARACTER*100 :: info_str

  write(info_str,'(1X,A,A)') &
      "* You have called the ScaLAPACK stub routine PDORMQR. ", &
      "Please check control.in and the libraries you linked."
  call localorb_info(info_str,use_unit,'(A)')
  stop

END SUBROUTINE PDORMQR
!******
!****s* FHI-aims/pdgels
!  NAME
!   pdgels
!  SYNOPSIS
subroutine pdgels( trans, m, n, nrhs, a, ia, ja, desca, &
               b, ib, jb, descb, work, lwork, info )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
!    and use of the actual routine corresponding to this stub!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

  use localorb_io

  implicit none

  character :: trans
  integer :: m, n, nrhs, ia, ja, desca(*), ib, jb, descb(*), lwork, info
  real*8 :: a(*), b(*), work(*)

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdgels
!******

!****s* FHI-aims/pdgesvd
!  NAME
!    pdgesvd
!  SYNOPSIS

subroutine pdgesvd( JOBU, JOBVT, M, N, &
&                   A, IA, JA, DESCA, &
&                   S, &
&                   U, IU,  JU,  DESCU,  &
&                   VT,  IVT, JVT, DESCVT, &
&                   WORK, LWORK, INFO )

  !  PURPOSE
  !
  !  USES

  use mpi_tasks
  implicit none

  !  ARGUMENTS

  character :: JOBU, JOBVT
  integer :: IA, INFO, IU, IVT, JA, JU, JVT, LWORK, M, N
  integer :: DESCA( * ), DESCU( * ), DESCVT( * )
  double precision :: A( * ), S( * ), U( * ), VT( * ), WORK( * )

  !  INPUTS
  !    dummy
  !  OUTPUTS
  !    dummy
  !  SOURCE

  character(*), parameter :: func = 'pdgesvd'

  call aims_stop('* A ScaLAPACK stub routine has been called.', func)

end subroutine pdgesvd
!******
!****s* FHI-aims/PDPOTRI
!  NAME
!   PDPOTRI
!  SYNOPSIS
SUBROUTINE PDPOTRI( UPLO, N, A, IA, JA, ADESC, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  CHARACTER :: UPLO
  INTEGER :: N, IA, JA, INFO
  INTEGER, dimension(9) :: ADESC
  REAL*8 :: A
  CHARACTER*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

END SUBROUTINE PDPOTRI
!******
!****s* FHI-aims/PZGETRF
!  NAME
!   PZGETRF
!  SYNOPSIS
SUBROUTINE PZGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  INTEGER :: M, N, IA, JA, INFO
  INTEGER DESCA (*), IPIV(*)
  COMPLEX*16 A (*)
  CHARACTER*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

END SUBROUTINE PZGETRF
!******
!****s* FHI-aims/PZGETRS
!  NAME
!   PZGETRS
!  SYNOPSIS
SUBROUTINE PZGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB, INFO )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      COMPLEX*16         A( * ), B( * )

  CHARACTER*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

END SUBROUTINE PZGETRS
!******

!****s* FHI-aims/pdsyev
!  NAME
!   pdsyev
!  SYNOPSIS
subroutine pdsyev( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, info )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  real*8 :: a, w, z
  real*8 :: gap, work
  integer :: iwork, ifail, iclustr
  real*8 :: vl, vu, abstol, orfac, info
  integer, dimension(9) :: desca, descz
  integer :: n, ia, ja, il, iu, m, nz, iz, jz, lwork, liwork
  character :: jobz, range, uplo

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdsyev
!******

subroutine pzheev(jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, rwork, lrwork, info)
  use localorb_io

  implicit none

  complex*16 :: a, z, work
  real*8 :: w, rwork
  integer, dimension(9) :: desca, descz
  integer :: n, ia, ja, iz, jz, lwork, lrwork, info
  character :: jobz, uplo

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzheev

!****s* FHI-aims/pzgesv
!  NAME
!   pzgesv
!  SYNOPSIS
subroutine pzgesv( n, nrhs, a, ia, ja, desca, ipiv, b, ib, jb, descb, info )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  integer       :: ia, ib, info, ja, jb, n, nrhs
  integer       :: desca( * ), descb( * ), ipiv( * )
  complex*16    :: a( * ), b( * )

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzgesv
!******


!****s* FHI-aims/pzlatra
!  NAME
!   pzlatra
!  SYNOPSIS
subroutine pzlatra( N, A, IA, JA, DESCA )
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE

  use localorb_io

  implicit none

  character*100 :: info_str
  integer       :: ia, ja, n
  integer       :: desca( * )
  complex*16    :: a( * )

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzlatra

!******
! VVG for mbd routine
subroutine pdgetri(int1,rel1,int2,int3,int4,int5,rel2,int6,int7,int8,info)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE
  use localorb_io
  implicit none
  integer:: int1,int2,int3,int4,int5,int6,int7,int8,info
  real*8:: rel1, rel2
  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdgetri

subroutine pzgetri(int1,cmplx1,int2,int3,int4,int5,cmplx2,int6,int7,int8,info)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE
  use localorb_io
  implicit none
  integer:: int1,int2,int3,int4,int5,int6,int7,int8,info
  complex*16 :: cmplx1, cmplx2
  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzgetri


subroutine pztranu(m, n, alpha, a, ia, ja, desc_a, beta, c, ic, jc, desc_c)
!  PURPOSE
!    This is a scalapack stub that is only compiled into the code if
!    the real scalapack is not available. Its only purpose is to make the
!    (in that case unused) scalapack parts of the code compile.
!    please refer to the scalapack/blacs manuals to infer about the purpose
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
!    Release version, FHI-aims (2008).
!  SOURCE
  use localorb_io
  implicit none

  integer                   :: m
  integer                   :: n
  complex*16                :: alpha
  complex*16, dimension(:)  :: a
  integer                   :: ia
  integer                   :: ja
  integer,    dimension(9)  :: desc_a
  complex*16                :: beta
  complex*16, dimension(:)  :: c
  integer                   :: ic
  integer                   :: jc
  integer,    dimension(9)  :: desc_c

  character*100 :: info_str


  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ', &
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pztranu

subroutine pdgeadd(TRANS,M,N,ALPHA,A,IA,JA,DESCA,BETA,C,IC,JC,DESCC)
  use localorb_io

  character :: TRANS
  integer :: M,N,IA,JA,IC,JC
  real*8 :: ALPHA,BETA,A,C
  integer, dimension(9) :: DESCA,DESCC

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ',&
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdgeadd

subroutine pzgeadd(TRANS,M,N,ALPHA,A,IA,JA,DESCA,BETA,C,IC,JC,DESCC)
  use localorb_io

  character :: TRANS
  integer :: M,N,IA,JA,IC,JC
  complex*16 :: ALPHA,BETA,A,C
  integer, dimension(9) :: DESCA,DESCC

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ',&
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzgeadd

subroutine pdsygst(IBTYPE,UPLO,N,A,IA,JA,DESCA,B,IB,JB,DESCB,SCA,INFO)
  use localorb_io

  character :: UPLO
  integer :: IA,IB,IBTYPE,INFO,JA,JB,N
  real*8 :: SCA
  integer, dimension(9) :: DESCA,DESCB
  real*8 :: A,B

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ',&
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdsygst

subroutine pzhegst(IBTYPE,UPLO,N,A,IA,JA,DESCA,B,IB,JB,DESCB,SCA,INFO)
  use localorb_io

  character :: UPLO
  integer :: IA,IB,IBTYPE,INFO,JA,JB,N
  real*8 :: SCA
  integer, dimension(9) :: DESCA,DESCB
  complex*16 :: A,B

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ',&
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pzhegst

subroutine pdtrmm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,IA,JA,DESCA,B,IB,JB,DESCB)
  use localorb_io

  character :: DIAG,SIDE,TRANS,UPLO
  integer :: IA,IB,JA,JB,M,N
  real*8 :: ALPHA,A,B
  integer, dimension(9) :: DESCA,DESCB

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ',&
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pdtrmm

subroutine pztrmm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,IA,JA,DESCA,B,IB,JB,DESCB)
  use localorb_io

  character :: DIAG,SIDE,TRANS,UPLO
  integer :: IA,IB,JA,JB,M,N
  complex*16 :: ALPHA,A,B
  integer, dimension(9) :: DESCA,DESCB

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ',&
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine pztrmm

subroutine blacs_pcoord(ICONTXT,PNUM,PROW,PCOL)
  use localorb_io

  integer :: ICONTXT,PNUM,PROW,PCOL

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ',&
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine blacs_pcoord

subroutine dgsum2d(ICONTXT,SCOPE,TOP,M,N,A,LDA,RDEST,CDEST)
  use localorb_io

  integer :: ICONTXT,M,N,LDA,RDEST,CDEST
  character :: SCOPE,TOP
  real*8 :: A

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ',&
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine dgsum2d

subroutine zgsum2d(ICONTXT,SCOPE,TOP,M,N,A,LDA,RDEST,CDEST)
  use localorb_io

  integer :: ICONTXT,M,N,LDA,RDEST,CDEST
  character :: SCOPE,TOP
  complex*16 :: A

  character*100 :: info_str

  write(info_str,'(1X,A,A)') '* You have called a ScaLAPACK stub routine. ',&
       'Please check control.in and the libraries you linked.'
  call localorb_info(info_str,use_unit,'(A)')
  stop

end subroutine zgsum2d
