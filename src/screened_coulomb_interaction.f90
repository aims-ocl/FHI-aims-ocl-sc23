!****s* FHI-aims/screened_coulomb_interaction
!  NAME
!   screened_coulomb_interaction
!  SYNOPSIS

      subroutine screened_coulomb_interaction &
          ( polar_freq, &
            screened_coulomb &
          )

!  PURPOSE
!  Subroutine screened_coulomb_interaction  evaluates the screened
!  Coulomb interaction W= 1/(V^(-1) - P).
!
!  USES

      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use runtime_choices
      implicit none

!  ARGUMENTS

      real*8  polar_freq(n_basbas, n_loc_prodbas)
      real*8  screened_coulomb(n_basbas, n_loc_prodbas)
!aux

!  INPUT
!  o polar_freq -- the polarisability within the auxiliary basis functions
!
!  OUTPUT
!  o screened_coulomb -- the screened Coulomb interaction matrix with the
!    auxiliary basis
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

!      real*8  temp_screened_coulomb(n_basbas, n_loc_prodbas)

!      real*8  wk(n_basbas*n_basbas)

!      integer info
!      integer ipiv(n_basbas)

!     counters


      integer :: i_basbas
      integer :: j_basbas
      integer :: i_task
      integer :: i_index

      integer info
      integer, dimension(:),  allocatable :: ipiv
      real*8, dimension(:),   allocatable :: wk
      real*8, dimension(:,:), allocatable :: temp_coulomb
!      real*8, dimension(1) :: wk

!     begin work

!     This is actually the inverse of the screened coulomb interaciton,
!     only the upper triangle is needed.

      i_task = myid + 1
      do j_basbas = 1, n_loc_prodbas, 1

        i_index = map_prodbas(j_basbas, i_task)
        do i_basbas = 1, n_basbas, 1

          if(i_basbas.eq.i_index) then
            screened_coulomb(i_basbas,j_basbas)= &
               1.d0 - &
               polar_freq(i_basbas,j_basbas)*2.d0/dble(n_spin)
          else
            screened_coulomb(i_basbas,j_basbas)= &
               - polar_freq(i_basbas,j_basbas)*2.d0/dble(n_spin)
          endif

!        if(myid.eq.0) then
!         write(STDERR,'(2I6,3f20.10)') i_basbas, j_basbas,
!     +      temp_screened_coulomb(i_basbas,j_basbas),
!     +      coulomb_matr(i_basbas,j_basbas),
!     +      polar_freq(i_basbas,j_basbas)
!         endif
        enddo
      enddo

      if (.not. (use_scgw .or. use_scgw0.or. use_dmft_gw))then
        if(use_scalapack) then
           call power_auxmat_scalapack(screened_coulomb,-1.d0,'')
        else
           call power_auxmat_lapack(screened_coulomb,-1.d0,'')
        endif
      elseif((use_scgw .or. use_scgw0.or.use_dmft_gw))then
!    Now invert the inverted screened interaction
        allocate(ipiv(n_basbas))
        allocate(temp_coulomb(n_basbas,n_basbas))
        allocate(wk(n_basbas*n_basbas))
        temp_coulomb (:,:) = 0.d0

        call gather_auxmat(temp_coulomb, screened_coulomb, n_basbas)
!----------------------------------------------------------------------
!        call  dsytrf ('U',n_basbas,temp_coulomb, &
!                n_basbas,ipiv,wk,n_basbas*n_basbas,info)
!
!        if (info.ne.0) then
!            write(STDERR,'(A,A,I5)') "factorization of the screened Coulomb ",&
!                       "interaction matrix fail"
!!              stop
!        else
!              call  dsytri ('U',n_basbas,temp_coulomb,&
!                   n_basbas,ipiv,wk1,info)
!              if(info.ne.0) then
!                write(STDERR,'(A,A,I5)') "inversion of the screened Coulomb ",&
!                  "matrix  fails "
!                stop
!              endif
!        endif
!-------------------------------------------------------------------------
        call dgetrf( n_basbas, n_basbas, temp_coulomb , &
                     n_basbas, ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then

            write(STDERR,*) " * Error info = ", info
            stop
          endif
        endif
! this invert the matrix in the LU factorized form
        call dgetri(n_basbas, temp_coulomb, n_basbas, &
                 ipiv, wk, n_basbas, info)

       if (info.ne.0) then
         if(myid.eq.0)then
            write(STDERR,*) " * Failure of matrix inversion "
            write(STDERR,*) " * Error info = ", info
            stop
          endif
        endif
!----------------------------------------------------------------------


        call scatter_auxmat(temp_coulomb , screened_coulomb, n_basbas)

        deallocate(ipiv)
        deallocate(temp_coulomb)
        deallocate(wk)



      endif 


      return
      end subroutine screened_coulomb_interaction
!---------------------------------------------------------------------
!****s* FHI-aims/screened_coulomb_interaction
!  NAME
!   screened_coulomb_interaction
!  SYNOPSIS

      subroutine screened_coulomb_interaction_2 &
          ( screened_coulomb )

!  PURPOSE
!  Subroutine screened_coulomb_interaction  evaluates the screened
!  Coulomb interaction W= 1/(V^(-1) - P).
!
!  USES

      use dimensions
      use prodbas
      use constants
      use mpi_tasks
      use runtime_choices
      use synchronize_mpi
      implicit none

!  ARGUMENTS

      real*8  screened_coulomb(n_basbas, n_basbas)
!aux

!  INPUT
!  o polar_freq -- the polarisability within the auxiliary basis functions
!
!  OUTPUT
!  o screened_coulomb -- the screened Coulomb interaction matrix with the
!    auxiliary basis
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

!      real*8  wk(n_basbas*n_basbas)

!      integer info
!      integer ipiv(n_basbas)

!     counters


      integer :: j_basbas
      integer :: j_basis_1
      integer :: i_task
      integer :: i_index

      integer info
      integer, dimension(:),  allocatable :: ipiv
      real*8, dimension(:),   allocatable :: wk
      real*8, dimension(:,:),   allocatable :: screened_coulomb_loc
!      real*8, dimension(1) :: wk

!     begin work

!     This is actually the inverse of the screened coulomb interaciton,
!     only the upper triangle is needed.

      screened_coulomb(:,:)= -screened_coulomb(:,:)*2.d0/dble(n_spin)
      do j_basbas = 1, n_basbas, 1

            screened_coulomb(j_basbas,j_basbas)= &
               1.d0 + screened_coulomb(j_basbas,j_basbas)
      enddo

      if (.not. (use_scgw .or. use_scgw0 .or. use_dmft_gw))then
        if(use_scalapack) then
           if (.not. allocated(screened_coulomb_loc)) then
             allocate(screened_coulomb_loc(n_basbas, n_loc_prodbas))
           endif

! ugly, only a temporary solution for the sake of time
           screened_coulomb_loc(:,:) = 0.d0
           do j_basis_1 = 1, n_loc_prodbas, 1
              j_basbas = map_prodbas(j_basis_1, myid+1)
              if(j_basbas .gt. 0 ) then
                screened_coulomb_loc(:, j_basis_1) = screened_coulomb(:,j_basbas)
              endif
           enddo

           call power_auxmat_scalapack(screened_coulomb_loc,-1.d0,'')

           screened_coulomb(:,:) = 0.d0
           do j_basis_1 = 1, n_loc_prodbas, 1
              j_basbas = map_prodbas(j_basis_1, myid+1)
              if(j_basbas == 0) cycle
              screened_coulomb(:, j_basbas) = screened_coulomb_loc(:,j_basis_1)
           enddo
           call sync_matrix(screened_coulomb, n_basbas, n_basbas)

           if(allocated(screened_coulomb_loc)) then
              deallocate(screened_coulomb_loc)
           endif
        else


          call power_genmat_lapack(n_basbas, screened_coulomb, -1.d0, &
                                safe_minimum, prodbas_threshold, '')


        endif
      elseif((use_scgw .or. use_scgw0 .or. use_dmft_gw))then

!    Now invert the inverted screened interaction
        allocate(ipiv(n_basbas))
        allocate(wk(n_basbas*n_basbas))

!        call gather_auxmat(temp_coulomb, screened_coulomb, n_basbas)
!----------------------------------------------------------------------
!        call  dsytrf ('U',n_basbas,temp_coulomb, &
!                n_basbas,ipiv,wk,n_basbas*n_basbas,info)
!
!        if (info.ne.0) then
!            write(STDERR,'(A,A,I5)') "factorization of the screened Coulomb ",&
!                       "interaction matrix fail"
!!              stop
!        else
!              call  dsytri ('U',n_basbas,temp_coulomb,&
!                   n_basbas,ipiv,wk1,info)
!              if(info.ne.0) then
!                write(STDERR,'(A,A,I5)') "inversion of the screened Coulomb ",&
!                  "matrix  fails "
!                stop
!              endif
!        endif
!-------------------------------------------------------------------------
        call dgetrf( n_basbas, n_basbas, screened_coulomb , &
                     n_basbas, ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then

            write(STDERR,*) " * Error info = ", info
            stop
          endif
        endif
! this invert the matrix in the LU factorized form
        call dgetri(n_basbas, screened_coulomb, n_basbas, &
                 ipiv, wk, n_basbas, info)

       if (info.ne.0) then
         if(myid.eq.0)then
            write(STDERR,*) " * Failure of matrix inversion "
            write(STDERR,*) " * Error info = ", info
            stop
          endif
        endif
!----------------------------------------------------------------------


!        call scatter_auxmat(temp_coulomb , screened_coulomb, n_basbas)

        deallocate(ipiv)
        deallocate(wk)



      endif 


      return
      end subroutine screened_coulomb_interaction_2
!---------------------------------------------------------------------
!*****
