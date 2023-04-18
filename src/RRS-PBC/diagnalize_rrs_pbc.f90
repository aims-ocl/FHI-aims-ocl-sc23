!
! Diagnalization coding for complex matrix equation 2009-04-21
!
        subroutine diagnalize_rrs_pbc(PBC_Dim,faex,olex,oe,IOP,k)
!    IOP  :: 0 for lower triangular matrix
!         :: 1 for upper triangular matrix
        use localorb_io, only: use_unit
        implicit none
!       IO variables
        Integer      :: PBC_Dim,IOP,k
        Complex*16   :: faex(PBC_Dim,PBC_Dim),olex(PBC_Dim,PBC_Dim)
        Real*8       :: oe(PBC_Dim)
!       Temp variables
        Integer :: LWORK
        Complex*16, allocatable :: WORK(:)
        Real*8, allocatable :: RWORK(:)

        LWORK=PBC_Dim*10
        allocate(WORK(LWORK),RWORK(LWORK))
        If (IOP.eq.0) then
!  Now diagnalizing for the low-triangle matrix
          call ZHEGV(1,'V','L',PBC_Dim,faex,PBC_Dim,olex,&
                      PBC_Dim,oe,WORK,LWORK,RWORK,k)
        ElseIf (IOP.eq.1) then
!  Now solve the eigenvalues from the upper-triangle matrix
          call ZHEGV(1,'V','U',PBC_Dim,faex,PBC_Dim,olex,&
                      PBC_Dim,oe,WORK,LWORK,RWORK,k)
        EndIf
        If (k.ne.0) then
          write(use_unit,*) "Error for Diagnalization",k
          !stop
        EndIf
        End subroutine

