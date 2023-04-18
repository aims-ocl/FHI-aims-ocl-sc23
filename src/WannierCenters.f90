!****s* FHI-aims/WannierCenters
!  NAME
!    WannierCenters
!  SYNOPSIS
MODULE WannierCenters
!*  PURPOSE
!*    This module contains various subroutines used for computing Wannier centers
!*    of charge. Loosely based on a limited and undocumented implementation by 
!*    Carlos Mera Acosta for computing Z2 based on the FHI-aims 2016 release.
!*  USES
  implicit none
!*  AUTHOR
!*    Christian Carbogno (FHI, carbogno@fhi-berlin.mpg.de)
!*  NOTES
!*    TODOs: 
!*       - NR: Test Hybrids
!*       - NR: Check edge cases (spin collinear, no SOC) carefully
!*       - NR: Add regression tests
!*       - NR: Add documentation and improve control.in directives
!*       - ZKY: Check for use_unit 
!*       - ZKY: LAPCKify            
!*       - ZKY: ScalaLAPCKify          
!*       - ZKY: Better Parallelization over k
!*       - ZKY: Only USE necessary stuff in out_plot_band.f90
!*       - Improve Output
!*       - Is it useful to store and retain D?
!*  HISTORY
!*    August 2019 - first release
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!*    the terms and conditions of the respective license agreement."
!*  SOURCE
SAVE
  complex*16, dimension(:,:,:,:), allocatable,public :: KS_eigenvector_for_WCC
  complex*16, dimension(:,:)    , allocatable,public :: overlap_for_WCC
  complex*16, dimension(:,:), allocatable,public     :: WCC_imag
  integer,public                                     :: WCC_nr_of_bands
  integer,public                                     :: WCC_highest_occ
  logical,public                                     :: WCC_calc = .false.

PRIVATE

  public :: WCC_get_dimensions
  public :: WCC_get_paths
  public :: WCC_convert_and_store_KS_evec
  public :: WCC_convert_and_store_KS_evec_no_SOC
  public :: WCC_check_occ_numbers
  public :: compute_WCC
  public :: ouput_WCC_to_file

CONTAINS

  ! CC: Get the primary dimensions for the arrays to be allocated
  subroutine WCC_get_dimensions( WCC_plane_index, WCC_nr_of_bands, n_k_points, n_k_points_task )
    use localorb_io
    use mpi_tasks
    use dimensions, only: Z2_n_plot, Z2_n_k_points
    implicit none
    integer,intent(in)  :: WCC_plane_index
    integer,intent(out) :: WCC_nr_of_bands, n_k_points, n_k_points_task
    !local:
    integer :: i_k_point
    character*300 :: info_str

    if ( (WCC_plane_index .eq. 1) .or. (WCC_plane_index .eq. 2) ) then
      WCC_nr_of_bands = Z2_n_plot
      write(info_str,'(2X,A,I5)') "Number of k-points for which the Wannier centers of charge are calculated for: ", WCC_nr_of_bands
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(2X,A)') "  Number of k-points used at each of this k-points to obtain the Wannier center of charge:"
      call localorb_info(info_str,use_unit,'(A)')
      call WCC_get_n_k_points(n_k_points)
      if ( Z2_n_k_points .lt. 0 ) then
        write(info_str,'(2X,A,I4)') "   Using default number of commensurate k-points per Wannier center: ", n_k_points
        call localorb_info(info_str)
        write(info_str,'(2X,A)') "   *** WARNING: This might not be a quantitatively safe default."
        call localorb_info(info_str)
      else
        write(info_str,'(2X,A,I4)') "   Default number of commensurate k-points per Wannier center would be: ", n_k_points
        call localorb_info(info_str)
        if (Z2_n_k_points < n_k_points ) then
          write(info_str,'(2X,A,I4,A,A)') "   *** WARNING: Requesting less k-points (",Z2_n_k_points,") than the default.",& 
                                   & " This is definitely not safe." 
          call aims_stop(info_str)
        else
          write(info_str,'(2X,A,I4,A,A)') "   Using user-requested number of k-points per Wannier center:", Z2_n_k_points
          call localorb_info(info_str)
          n_k_points = Z2_n_k_points
        end if
      end if
      n_k_points_task = 0
      do i_k_point = 1, n_k_points, 1
         if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points )then
            n_k_points_task = n_k_points_task + 1
         end if
      end do
    else
      call aims_stop('WCC_plane_index not defined')
    end if
  end subroutine WCC_get_dimensions

  ! CC: Construct plane-dependent path through BZ
  subroutine WCC_get_paths( i_plane, n_bands, n_k_points, band_begin, band_end, band_points )
    use localorb_io
    use mpi_tasks, only: aims_stop
    IMPLICIT NONE
    integer,intent(in) :: i_plane, n_bands, n_k_points
    real*8, dimension(n_bands,3), intent(out) :: band_begin, band_end
    integer,dimension(n_bands),   intent(out) :: band_points
    integer :: i_point
    character*300 :: info_str

    if ( (i_plane .lt. 1) .or. (i_plane .gt. 2) ) then
      call aims_stop('WCC_plane_index not defined')
    end if

    band_points(:) = n_k_points
    if(i_plane.eq.1) then
      band_begin(:,3) =  0.0d0
      band_end(:,3)   =  0.0d0
    elseif(i_plane.eq.2) then
      band_begin(:,3) =  0.5
      band_end(:,3)   =  0.5
    endif
    band_begin(:,1) =  0.5d0
    band_end(:,1)   = -0.5d0
    do i_point = 1, n_bands
      band_begin(i_point,2) =   0.0d0 + dble(i_point-1)*(0.5d0/dble(n_bands-1))
      band_end(i_point,2)   =   0.0d0 + dble(i_point-1)*(0.5d0/dble(n_bands-1))
    enddo
    call localorb_info('')
    write(info_str,'(2X,A,F7.3,A,F7.3,A,I4,A)') "Wannier center evolution is sampled along the k_y direction from k_y=", band_begin(1,2), " to ky=", band_end(n_bands,2), " in ",n_bands, " steps."
    call localorb_info(info_str)
    write(info_str,'(2X,A,F7.3,A,F7.3,A,F7.3,A,F7.3,A)') "- For each k_y, we sample the path between (",band_begin(1,1),", k_y, ",band_begin(1,3),") and (",band_end(1,1),", k_y, ",band_end(1,3),")"
    call localorb_info(info_str)
  end subroutine WCC_get_paths

  ! CC: Store the KS eigenvectors (SOC case)
  subroutine WCC_convert_and_store_KS_evec( i_k_point, &
                 n_basis, n_states, n_spin, n_k_points_task, KS_eigenvector_complex , &
                 n_states_soc, eigenvec_soc_wf_basis, &
                 KS_eigenvector_for_WCC )
    use soc_utilities, only: convert_wf_basis_to_compute_basis
    implicit none
    integer, intent(in) :: i_k_point
    integer, intent(in) :: n_basis, n_states, n_spin, n_k_points_task, n_states_soc
    complex*16,dimension(n_basis,n_states,n_spin,n_k_points_task), intent(in)    :: KS_eigenvector_complex
    complex*16,dimension(n_states_soc,n_states_soc),               intent(in)    :: eigenvec_soc_wf_basis
    complex*16,dimension(n_basis,n_states_soc,2,n_k_points_task),  intent(inout) :: KS_eigenvector_for_WCC
    !local:
    real*8,dimension(1,1,1,1) :: dummy_real_KS
    complex*16,dimension(2*n_basis,2*n_states, 1) :: tmp_KS
    integer :: i_spin, basis_offset

    ! From n_states_soc --> n_basis
    call convert_wf_basis_to_compute_basis( &
       n_basis,      n_states,           dummy_real_KS, KS_eigenvector_complex(:,:,:,i_k_point), &
       n_states_soc, n_states_soc,       eigenvec_soc_wf_basis, &
       2*n_basis ,   2*n_states,         tmp_KS )

    ! From implicit to explicit spin
    do i_spin = 1, 2, 1
      if (i_spin .eq. 1) then
        ! Spin-up components are requested
        basis_offset = 0
      else
        ! Spin-dn components are requested
        basis_offset = n_basis
      end if

      KS_eigenvector_for_WCC(:,:,i_spin,i_k_point) = tmp_KS(basis_offset+1:basis_offset+n_basis,:,1)
    end do
  end subroutine WCC_convert_and_store_KS_evec

  ! CC: Store the KS eigenvectors (non-SOC) case
  subroutine WCC_convert_and_store_KS_evec_no_SOC( &
             n_basis, n_states, n_spin, n_k_points, n_k_points_task, occ_numbers, &
             KS_eigenvalue, KS_eigenvector_complex , KS_eigenvector_for_WCC, WCC_highest_occ )
    use dimensions_soc, only: n_states_sr, sr_state_start, n_states_soc,&
        n_basis_soc, n_basis_soc_coll, n_basis_soc_ncoll, n_saved_states_soc, &
        soc_saved_state_start, n_core_states_omit_from_soc
    use mpi_tasks
    implicit none
    integer, intent(in) :: n_basis, n_states, n_spin, n_k_points, n_k_points_task, WCC_highest_occ
    real*8,dimension(n_states,n_spin,n_k_points), intent(in)                  :: occ_numbers, KS_eigenvalue
    complex*16,dimension(n_basis,n_states,n_spin,n_k_points_task), intent(in) :: KS_eigenvector_complex
    complex*16,dimension(n_basis,2*n_states,2,n_k_points_task),intent(out)    :: KS_eigenvector_for_WCC
    ! local
    integer :: i_state,i_k_point, i_k
    real*8,dimension(2*n_states,1,n_k_points)  :: occ_numbers_soc, KS_eigenvalue_soc
    complex*16,dimension(2*n_states,2*n_states) :: eigenvec_soc_wf_basis, eigenvec_soc_wf_basis_def
    integer :: n_states_sr_tmp, sr_state_start_tmp, n_states_soc_tmp, n_basis_soc_tmp, &
             & n_basis_soc_coll_tmp, n_basis_soc_ncoll_tmp, n_saved_states_soc_tmp, &
             & soc_saved_state_start_tmp, n_core_states_omit_from_soc_tmp

    ! CC: Ugly: Just mimic SOC behaviour and exploits its functions
    ! Dummy setup of eigenvec_soc_wf_basis:
    eigenvec_soc_wf_basis_def(:,:) = 0.0d0
    do i_state=1,n_states,1
      eigenvec_soc_wf_basis_def( i_state,          (2*i_state)-1 ) = (1.0d0,0.0d0)
      eigenvec_soc_wf_basis_def( n_states+i_state, (2*i_state)   ) = (1.0d0,0.0d0)
    end do

    ! CC: Store original settings as safety measure
    n_states_sr_tmp                  =   n_states_sr
    sr_state_start_tmp               =   sr_state_start
    n_states_soc_tmp                 =   n_states_soc
    n_basis_soc_tmp                  =   n_basis_soc
    n_basis_soc_coll_tmp             =   n_basis_soc_coll
    n_basis_soc_ncoll_tmp            =   n_basis_soc_ncoll
    n_saved_states_soc_tmp           =   n_saved_states_soc
    soc_saved_state_start_tmp        =   soc_saved_state_start
    n_core_states_omit_from_soc_tmp  =   n_core_states_omit_from_soc
    ! CC: Create KS in SOC notation:
    n_states_sr                  = n_states
    sr_state_start               = 1
    n_states_soc                 = 2*n_states
    n_basis_soc                  = 2*n_basis
    n_basis_soc_coll             = 2*n_basis
    n_basis_soc_ncoll            = 0
    n_saved_states_soc           = 2*n_states
    soc_saved_state_start        = 1
    n_core_states_omit_from_soc  = 0


    i_k = 0
    do i_k_point = 1, n_k_points

      ! CC: i_spin 1 and 2 are obviously different
      if (n_spin .eq. 1) then
        do i_state=1,n_states,1
          KS_eigenvalue_soc( (2*i_state)-1 ,1, i_k_point ) = KS_eigenvalue(i_state,1,i_k_point)
          KS_eigenvalue_soc( (2*i_state)   ,1, i_k_point ) = KS_eigenvalue(i_state,1,i_k_point)
          occ_numbers_soc( (2*i_state)-1 ,1, i_k_point )   = occ_numbers(i_state,1,i_k_point)/2
          occ_numbers_soc( (2*i_state)   ,1, i_k_point )   = occ_numbers(i_state,1,i_k_point)/2
        end do
      else
        do i_state=1,n_states,1
          KS_eigenvalue_soc( (2*i_state)-1 ,1, i_k_point ) = KS_eigenvalue(i_state,1,i_k_point)
          KS_eigenvalue_soc( (2*i_state)   ,1, i_k_point ) = KS_eigenvalue(i_state,2,i_k_point)
          occ_numbers_soc( (2*i_state)-1 ,1, i_k_point )   = occ_numbers(i_state,1,i_k_point)
          occ_numbers_soc( (2*i_state)   ,1, i_k_point )   = occ_numbers(i_state,2,i_k_point)
        end do
      end if
      eigenvec_soc_wf_basis(:,:) = eigenvec_soc_wf_basis_def(:,:)

      ! CC: For i_spin 2, sort eigenvalues again in the unified pic, not
      !      separatly by spin channel, since KS_eigenvalues are now not
      !      necessarily ascending for the SOC 2*n_states notation:
      if (n_spin .eq. 2) then
        call WCC_bubblesort( 2*n_states, KS_eigenvalue_soc(:,1,i_k_point), &
               occ_numbers_soc(:,1,i_k_point), eigenvec_soc_wf_basis )
      end if

      ! CC: Now finally convert the evecs
      if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points )then
        i_k = i_k + 1
        call WCC_convert_and_store_KS_evec( i_k, &
            n_basis, n_states, n_spin, n_k_points_task, KS_eigenvector_complex , &
            2*n_states, eigenvec_soc_wf_basis, KS_eigenvector_for_WCC )
      end if

    end do

    ! CC: Check sorted occ. numbers:
    call WCC_check_occ_numbers( 2*n_states, 1,  n_k_points, occ_numbers_soc, WCC_highest_occ )

    ! CC: Revert to original settings:
    n_states_sr                  =   n_states_sr_tmp
    sr_state_start               =   sr_state_start_tmp
    n_states_soc                 =   n_states_soc_tmp
    n_basis_soc                  =   n_basis_soc_tmp
    n_basis_soc_coll             =   n_basis_soc_coll_tmp
    n_basis_soc_ncoll            =   n_basis_soc_ncoll_tmp
    n_saved_states_soc           =   n_saved_states_soc_tmp
    soc_saved_state_start        =   soc_saved_state_start_tmp
    n_core_states_omit_from_soc  =   n_core_states_omit_from_soc_tmp
  end subroutine WCC_convert_and_store_KS_evec_no_SOC

  ! CC: Check that all N-lowest states are indeed occupied along the path
  subroutine WCC_check_occ_numbers( n_states, n_spin,  n_k_points, occ_numbers, WCC_highest_occ )
    use localorb_io
    use mpi_tasks, only: aims_stop
    implicit none
    integer, intent(in)  :: n_states, n_spin,  n_k_points, WCC_highest_occ
    real*8, dimension(n_states, n_spin,  n_k_points) :: occ_numbers
    !local
    integer       :: this_k_point, i_spin
    character*200 :: info_str
    real*8        :: lowest_N_bands

    do  this_k_point = 1,  n_k_points
      do i_spin=1,n_spin
        lowest_N_bands = sum(occ_numbers(1:WCC_highest_occ,i_spin,this_k_point))
        if ( ( abs(lowest_N_bands - WCC_highest_occ) .gt. 1.0d-4) ) then
          write(info_str,'(2X,A)')            " *** WARNING: Number of occupied states not constant across Brillouin zone:"
          call localorb_info(info_str)
          write(info_str,'(2X,A,I6,A,F12.4)') " ***          ", WCC_highest_occ, " != ", lowest_N_bands
          call localorb_info(info_str)
          if ( ( abs(lowest_N_bands - WCC_highest_occ) .gt. 1.0d-2) ) then
             call aims_stop(" *** Do we have a metal here? Abort! ")
          end if
        end if
      end do
    end do
  end subroutine WCC_check_occ_numbers

  ! CC: Main Routine for computing the Wannier centers
  subroutine compute_WCC(i_band,n_states,n_basis,n_k_points,n_k_points_task, &
    & WCC_highest_occ,WCC_n_plot,KS_eigenvector, overlap_matrix_w_complex, Wcc_ima)
    use localorb_io
    use mpi_tasks
    use synchronize_mpi_basic
    implicit none
    integer, intent(in)    :: i_band, n_states, n_basis, n_k_points, & 
                            & n_k_points_task, WCC_highest_occ, WCC_n_plot
    complex*16,intent(in),dimension(n_basis,2*n_states,2,n_k_points_task)    :: KS_eigenvector
    complex*16,intent(in),dimension(n_basis*(n_basis+1)/2,n_k_points_task)   :: overlap_matrix_w_complex
    complex*16, dimension(WCC_highest_occ,WCC_n_plot), intent(inout)         :: Wcc_ima
    ! Local
    complex*16, dimension(n_basis, 2*n_states,2) :: KS_eigenvector_reference
    complex*16, dimension(2*n_states,2*n_states) :: WCC_shell_kkprime
    complex*16, dimension(2*n_states,2*n_states) :: D
    integer :: i_state, j_state, i_k, i_k_point,i_basis, info
    complex*16 :: prod
    character*300 :: info_str
    complex*16, dimension(2*n_states,2*n_states) :: sum_shell
    complex*16, dimension(2*n_states,2*n_states) :: Dold

    WCC_shell_kkprime(:,:) = (0.d0,0.d0)
    KS_eigenvector_reference(:,:,:) = (0.0d0,0.0d0)

    ! CC: D is initialized as 0.0
    !     except for diagonal elements 
    !     with D(i,i)=1.0 for i < WCC_highest_occ
    D(:,:) = (0.0d0,0.0d0)
    do i_state = 1, WCC_highest_occ
      D(i_state,i_state) = (1.0d0,0.d0)
    enddo

    ! CC: Send 
    !   KS_eigenvector(:,:,:,1) --> KS_eigenvector_reference
    !    We thus start with kp=1 as a reference
    if(myid == MOD(1,n_tasks) .and. myid <= n_k_points) then
      KS_eigenvector_reference(:,:,:) = KS_eigenvector(:,:,:,1)
    else
      KS_eigenvector_reference(:,:,:) = (0.0d0,0.0d0)
    endif
    call sync_vector_complex( KS_eigenvector_reference, 2*n_states*2*n_basis )

    ! CC: Loop over i_k_point, n_k_points_task --> FIXME clean i_k notation
    !  Iteration | i_k_point     | i_k              | i_k_ref
    !  -----------------------------------------------------------------------
    !      1     | n_k_points-1  | n_k_points_task    | 1
    !      2     | n_k_points-2  | n_k_points_task-1  | n_k_points_task
    !      3     | n_k_points-3  | n_k_points_task-2  | n_k_points_task-1
    Dold = D
    i_k = n_k_points_task
    sum_shell(:,:) = (0.0d0, 0.0d0)
    do i_k_point = n_k_points, 1, -1

      ! CC: Since the reference KS_eigenvectors at i_k_ref is already present on all nodes,
      !      here only the ONE node that has the KS_eigenvector at i_k in memory
      !      does the work, all others set WCC_shell_kkprime == 0. Then, the result is ``synced''
      !      via an MPI_SUM
      WCC_shell_kkprime(:,:) = (0.0d0,0.0d0)
      if (myid == MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then

         ! CC: Compute entry (i_state,j_state) of WCC_shell_kkprime by contracting KS_vec at i_k and KS_vec at i_k_ref
         ! SHELL(l,m) = \sum_i,j [ VEC^*(i,l,1,i_k)*VECREF(j,m,1,i_kref) + VEC^*(i,l,2,i_k)*VECREF(j,m,2,i_kref) ] * OVL(index)
         call compute_shell_kkprime( n_states, n_basis, n_k_points_task, i_k, KS_eigenvector, KS_eigenvector_reference, overlap_matrix_w_complex, WCC_shell_kkprime )

         ! CC: The old Eigenvec i_k now becomes the reference for the next iteration
         KS_eigenvector_reference(:,:,:) = KS_eigenvector(:,:,:,i_k)
         i_k = i_k - 1
      else
         WCC_shell_kkprime(:,:) = (0.d0,0.d0)
         KS_eigenvector_reference(:,:,:) = (0.0d0,0.0d0)
      endif

      ! CC: Synchronize the KS eigenvectors and shells
      call sync_vector_complex( KS_eigenvector_reference, 2*n_states*2*n_basis )
      call sync_vector_complex( WCC_shell_kkprime, 2*n_states*2*n_states )
      sum_shell = sum_shell + WCC_shell_kkprime

      ! CC: Now compute D(:,:): D(i,j) = sum_k WCC_shell_kkprime(i,k)*old_D(k,j)
      call update_D( n_states, WCC_highest_occ, WCC_shell_kkprime, D)
    enddo

    ! CC: Diagonalize the part of D pertaining to occupied states and sort the
    !     eigenvalues so to obtain Wcc at each k_point i_band
    if(myid == 0) then
      call from_D_to_Wcc( D, WCC_highest_occ, n_states,  Wcc_ima(:,i_band) )
    else
      Wcc_ima(:,i_band) = (0.d0,0.d0)
    end if
  end subroutine compute_WCC

  ! CC: Output of WCC for plane i_plane:
  subroutine ouput_WCC_to_file(i_plane, WCC_highest_occ, WCC_n_plot, Wcc_ima )
    IMPLICIT NONE
    integer,intent(in)                                          :: i_plane, WCC_highest_occ, WCC_n_plot
    complex*16,intent(in),dimension(WCC_highest_occ, WCC_n_plot)  :: Wcc_ima
    integer    :: i_state, i_band
    character*50 :: file_name

    write(file_name, '(A3,I1,A4)') 'WCC',i_plane,'.dat'
    open(88, file=file_name)
    do i_state = 1, WCC_highest_occ, 1
      do i_band = 1, WCC_n_plot, 1
        write(88,*) i_band, dimag(cdlog(Wcc_ima(i_state,i_band)))
      enddo
      write(88,*) '                                                          '
    enddo
    close(88)
  end subroutine ouput_WCC_to_file


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! CC: Here start the private routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! CC: Diagonalize a Complex, Non-Symmetric matrix
  subroutine diagD(n,H,Eigenvalue)
   implicit none
   integer, intent(in)                       :: n
   complex*16, intent(in), dimension(n,n)    :: H
   complex*16, dimension(n,n)                :: cualL
   complex*16, dimension(n,n)                :: cualR
   complex*16, intent(out), dimension(n)     :: Eigenvalue
   integer                                   :: information
   double precision, dimension(n*n)          :: rwork
   complex*16, dimension(n*n*n)              :: work
   integer                                 :: lwork

   lwork = n*n*n
   rwork = 0.d0
   work  = (0.d0,0.d0)

   call zgeev('N','N',n,H,n,Eigenvalue,cualL,1,cualR,1,work,lwork,rwork,information)
   if(information.ne.0) then
     write(*,*) "Falha ao calcular os autovalores."
     stop
   endif
   return
   endsubroutine diagD


  ! CC: Diagonalize the part of D pertaining to occupied states and sort the
  !     eigenvalues so to obtain Wcc at one specific k_point
  subroutine from_D_to_Wcc( D, WCC_highest_occ, n_states, Wcc )
    use localorb_io
    IMPLICIT NONE
    integer,intent(in)                                     :: WCC_highest_occ,n_states
    complex*16,dimension(2*n_states,2*n_states),intent(in) :: D
    complex*16,dimension(WCC_highest_occ),intent(out)       :: Wcc
    !local
    complex*16,dimension(WCC_highest_occ,WCC_highest_occ) :: D_occ
    complex*16  :: dummy_comp
    integer                                   :: i_state,j_state
    character*300 :: info_str

    ! CC: Store D in D_occ: Dimensionality reduction:
    !  D(2*n_states, 2*n_states) -->  D_occ(WCC_highest_occ,WCC_highest_occ)
    D_occ(:,:) = (0.0d0,0.0d0)
    do i_state = 1, WCC_highest_occ,1
      do j_state = 1,WCC_highest_occ,1
        D_occ(i_state,j_state) = D(i_state,j_state)
      enddo
    enddo

    ! CC: Diagonalize D_occ to obtain Wcc
    call diagD(WCC_highest_occ, D_occ, Wcc)
    ! CC: Sort the eigenvalues
    do i_state = 1,WCC_highest_occ,1
      do j_state = i_state + 1,WCC_highest_occ,1
        if (dimag(cdlog(Wcc(i_state))) > dimag(cdlog(Wcc(j_state)))) then
          dummy_comp = Wcc(i_state)
          Wcc(i_state) = Wcc(j_state)
          Wcc(j_state) = dummy_comp
        end if
      end do
    end do

  end subroutine from_D_to_Wcc

  ! CC: Now compute D(:,:): D(i,j) = sum_k WCC_shell_kkprime(i,k)*old_D(k,j)
  subroutine update_D( n_states, WCC_highest_occ, WCC_shell_kkprime, D)
    IMPLICIT NONE
    integer,intent(in)                                        :: WCC_highest_occ,n_states
    complex*16,dimension(2*n_states,2*n_states),intent(in)    :: WCC_shell_kkprime
    complex*16,dimension(2*n_states,2*n_states),intent(inout) :: D
    !local:
    integer :: i_state, j_state, k_state
    complex*16,dimension(2*n_states,2*n_states)               :: Dold

    Dold(:,:) = D(:,:)
    D(:,:) = (0.d0, 0.d0)
    do i_state = 1, n_states*2  , 1
      do j_state = 1, n_states*2  , 1
        do k_state = 1, WCC_highest_occ, 1
          D(i_state,j_state) =  D(i_state,j_state) + &
          WCC_shell_kkprime(i_state,k_state)*Dold(k_state,j_state)
        enddo
      enddo
    enddo
  end subroutine update_D


  ! CC: Compute entry (i_state,j_state) of WCC_shell_kkprime by contracting KS_vec at i_k and KS_vec at i_k_ref
  subroutine compute_shell_kkprime(n_states, n_basis, n_k_points_task, i_k, KS_eigenvector, KS_eigenvector_reference, overlap_matrix_w_complex, shell )
    implicit none
    integer, intent(in) :: n_states, n_basis, n_k_points_task, i_k
    complex*16,intent(in),dimension(n_basis,2*n_states,2,n_k_points_task)    :: KS_eigenvector
    complex*16,intent(in),dimension(n_basis, 2*n_states,2)                   :: KS_eigenvector_reference
    complex*16,intent(in),dimension(n_basis*(n_basis+1)/2,n_k_points_task)   :: overlap_matrix_w_complex
    complex*16,intent(out),dimension(2*n_states,2*n_states)                  :: shell
    ! Local vars
    integer :: i_basis, j_basis
    integer :: i_state, j_state

    ! FIXME: LAPACK-ify
    shell(:,:) = (0.0d0,0.0d0)
    do i_state = 1, 2*n_states, 1
      do j_state = 1, 2*n_states, 1
        do i_basis = 1, n_basis, 1
          do j_basis = 1, n_basis, 1
            if(j_basis.ge.i_basis) then
              shell(i_state,j_state) = shell(i_state,j_state) + &
               ((conjg(KS_eigenvector(i_basis,i_state,1,i_k))&
                      * KS_eigenvector_reference(j_basis,j_state,1))&
              + (conjg(KS_eigenvector(i_basis,i_state,2,i_k)) &
                      * KS_eigenvector_reference(j_basis,j_state,2)))&
              * overlap_matrix_w_complex((i_basis + (((j_basis)*(j_basis-1))/2)),i_k)
            elseif(j_basis.lt.i_basis) then
              shell(i_state,j_state) = shell(i_state,j_state) + &
               ((conjg(KS_eigenvector(i_basis,i_state,1,i_k))&
                      * KS_eigenvector_reference(j_basis,j_state,1))&
              + (conjg(KS_eigenvector(i_basis,i_state,2,i_k)) &
                      * KS_eigenvector_reference(j_basis,j_state,2)))&
              * conjg(overlap_matrix_w_complex((j_basis + (((i_basis)*(i_basis-1))/2)),i_k))
            endif
          enddo
        enddo
      enddo
    enddo
  end subroutine compute_shell_kkprime

  ! CC: get n-kpoints -- Check relation to cell_index
  subroutine WCC_get_n_k_points(n_k_points)
    use pbc_lists, only: cell_index,n_cells
    IMPLICIT NONE
    integer, intent(out) :: n_k_points
    integer :: i,i_cell_n
    i = 0
    do i_cell_n = 1, n_cells
       if(cell_index(i_cell_n,1).ne.0.and.cell_index(i_cell_n,2).eq.0.and.cell_index(i_cell_n,3).eq.0)  then
          i = i + 1
       endif
    end do
    n_k_points = i + 1
  end subroutine WCC_get_n_k_points

  ! CC: Sort KS_eval and re-arrange occ/matrix accordingly
  subroutine WCC_bubblesort( n_states_soc, KS_eval, occ,  matrix )
    implicit none
    integer,intent(in)                                            :: n_states_soc
    real*8,dimension(n_states_soc),intent(inout)                  :: KS_eval, occ
    complex*16,dimension(n_states_soc,n_states_soc),intent(inout) :: matrix
    !local
    integer                            :: i_state,j_state
    real*8,dimension(n_states_soc)     :: KS_eval_tmp, occ_tmp
    real*8                             :: dummy_eval, dummy_occ
    complex*16,dimension(n_states_soc,n_states_soc) :: matrix_tmp
    complex*16,dimension(n_states_soc) :: dummy_matrix

    KS_eval_tmp = KS_eval
    occ_tmp     = occ
    matrix_tmp  = matrix

    do i_state = 1, n_states_soc, 1
      do j_state = i_state + 1, n_states_soc, 1
        if (KS_eval_tmp(i_state) > KS_eval_tmp(j_state)) then
          dummy_eval      = KS_eval_tmp(i_state)
          dummy_occ       = occ_tmp(i_state)
          dummy_matrix(:) = matrix_tmp(:,i_state)

          KS_eval_tmp(i_state)  = KS_eval_tmp(j_state)
          occ_tmp(i_state)      = occ_tmp(j_state)
          matrix_tmp(:,i_state) = matrix_tmp(:,j_state)

          KS_eval_tmp(j_state)  = dummy_eval
          occ_tmp(j_state)      = dummy_occ
          matrix_tmp(:,j_state) = dummy_matrix(:)
        end if
      end do
    end do

    KS_eval = KS_eval_tmp
    occ     = occ_tmp
    matrix  = matrix_tmp

  end subroutine WCC_bubblesort

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! CC: Here start the debug routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! CC: Below follow some debug, input, and output routines that were extremly
  !     useful when debugging, but that are neither well tested nor cleaned up.
  !     For bug-hunting in this module, especially for parallel bugs, they might
  !     be helpful, though:

  !! DBG: ! CC: Output D for plane i_plane:
  !! DBG: subroutine ouput_D_to_file(i_plane, n_states, Z2_n_plot, Dstore, Z2_highest_occ )
  !! DBG:   IMPLICIT NONE
  !! DBG:   integer,intent(in)                                                  :: i_plane, n_states, Z2_n_plot, Z2_highest_occ
  !! DBG:   complex*16,intent(in),dimension(2*n_states, 2*n_states, Z2_n_plot)  :: Dstore
  !! DBG:   integer    :: i_state, i_band, j_state
  !! DBG:   character*99 :: file_name

  !! DBG:   write(file_name, '(A,I1,A4)') 'Dstore_diag',i_plane,'.dat'
  !! DBG:   open(88, file=trim(file_name))
  !! DBG:   write(file_name, '(A,I1,A4)') 'Dstore_offdiag',i_plane,'.dat'
  !! DBG:   open(89, file=trim(file_name))
  !! DBG:   do i_state = 1, Z2_highest_occ
  !! DBG:     do i_band = 1, Z2_n_plot, 1
  !! DBG:       write(88,'(2I8,2F20.3)') i_state, i_band, Dstore(i_state,i_state,i_band)
  !! DBG:     enddo
  !! DBG:     write(88,*) '                                                          '
  !! DBG:     do j_state = 1, Z2_highest_occ
  !! DBG:       do i_band = 1, Z2_n_plot, 1
  !! DBG:         write(89,'(3I8,2F20.3)') i_state,j_state,i_band, Dstore(i_state,j_state,i_band)
  !! DBG:       enddo
  !! DBG:       write(89,*) '                                                          '
  !! DBG:     enddo
  !! DBG:   enddo
  !! DBG:   close(88)
  !! DBG:   close(89)

  ! DBG: end subroutine ouput_D_to_file

  !! DBG subroutine Z2_print_ovl_on_tasks(n_basis, overlap_matrix_w_complex)
  !! DBG   use localorb_io
  !! DBG   use mpi_tasks
  !! DBG   use synchronize_mpi_basic
  !! DBG   implicit none
  !! DBG   integer, intent(in)                                                      :: n_basis
  !! DBG   complex*16,intent(in),dimension(n_basis*(n_basis+1)/2)                   :: overlap_matrix_w_complex
  !! DBG   ! local
  !! DBG   integer :: i_task,i_state
  !! DBG   complex*16,dimension(n_basis*(n_basis+1)/2)                              :: ovlp_tmp
  !! DBG   character*300 :: info_str

  !! DBG   do i_task=0,n_tasks-1,16
  !! DBG     if (myid == i_task) then
  !! DBG       ovlp_tmp = overlap_matrix_w_complex
  !! DBG     else
  !! DBG       ovlp_tmp = (0.0d0,0.0d0)
  !! DBG     end if
  !! DBG     call sync_vector_complex( ovlp_tmp , n_basis*(n_basis+1)/2 )

  !! DBG     do i_state=1,10
  !! DBG       write(info_str,'(2X,A,2I6,2F20.5)') " Z2: OVLP on TASK", i_task, i_state, ovlp_tmp(i_state)
  !! DBG       call localorb_info(info_str)
  !! DBG     end do

  !! DBG   end do
  !! DBG end subroutine Z2_print_ovl_on_tasks

  !! DBG subroutine Z2_print_KS_on_tasks_v2(n_basis, n_states, n_k_points_task, KS_eigenvector_soc_perturbed )
  !! DBG   use localorb_io
  !! DBG   use mpi_tasks
  !! DBG   use synchronize_mpi_basic
  !! DBG   implicit none
  !! DBG   integer, intent(in)                                                      :: n_basis, n_states, n_k_points_task
  !! DBG   complex*16,intent(in),dimension(n_basis,2*n_states,2,n_k_points_task)    :: KS_eigenvector_soc_perturbed
  !! DBG   ! local
  !! DBG   integer :: i_task,i_state,i_basis,i_k,i_spin, i_k_point
  !! DBG   complex*16,dimension(:,:,:),allocatable               :: KS_tmp,KS_tmp_ref
  !! DBG   character*300 :: info_str
  !! DBG   integer,dimension(0:n_tasks-1) :: n_basis_arr, n_states_arr, n_k_points_task_arr
  !! DBG   complex*16   :: prod,semphas
  !! DBG   ! DBG:
  !! DBG   integer :: nt_dbg, n_k_points


  !! DBG   ! Get max and min dimensions first
  !! DBG   write(info_str,'(2X,A)') " Z2: Z2_print_KS_on_tasks HERES"
  !! DBG   call localorb_info(info_str)
  !! DBG   write(info_str,'(2X,A,2I4)') " Z2: myid, n_tasks", myid, n_tasks
  !! DBG   call localorb_allinfo(info_str)
  !! DBG   call mpi_barrier(mpi_comm_global,nt_dbg)

  !! DBG   n_basis_arr(:) = 0
  !! DBG   n_states_arr(:) = 0
  !! DBG   n_k_points_task_arr(:) = 0
  !! DBG   n_basis_arr(myid) = n_basis
  !! DBG   n_states_arr(myid) = n_states
  !! DBG   n_k_points_task_arr(myid) = n_k_points_task
  !! DBG    write(info_str,'(2X,A,3I6)') " Z2: Z2_print_KS_on_tasks n_basis_arr bef sync",myid,n_tasks,SIZE(n_basis_arr)
  !! DBG    call localorb_allinfo(info_str)
  !! DBG   call sync_integer_vector( n_basis_arr, n_tasks )
  !! DBG    write(info_str,'(2X,A)') " Z2: Z2_print_KS_on_tasks n_basis_arr"
  !! DBG    call localorb_info(info_str)
  !! DBG   call sync_integer_vector( n_states_arr, n_tasks )
  !! DBG    write(info_str,'(2X,A)') " Z2: Z2_print_KS_on_tasks n_tasks"
  !! DBG    call localorb_info(info_str)
  !! DBG   call sync_integer_vector( n_k_points_task_arr, n_tasks )
  !! DBG    write(info_str,'(2X,A,3I6)') " Z2: MAX_DIMS", maxval(n_basis_arr,1), maxval(n_states_arr,1), maxval(n_k_points_task_arr,1)
  !! DBG    call localorb_info(info_str)
  !! DBG    write(info_str,'(2X,A,3I6)') " Z2: MIN_DIMS", minval(n_basis_arr,1), minval(n_states_arr,1), minval(n_k_points_task_arr,1)
  !! DBG    call localorb_info(info_str)
  !! DBG   allocate( KS_tmp( maxval(n_basis_arr,1), maxval(n_states_arr,1 ), 2 ))
  !! DBG   KS_tmp(:,:,:) = (0.0d0,0.0d0)
  !! DBG    do i_task=0,n_tasks-1
  !! DBG      write(info_str,'(2X,A,3I6)') " Z2: DIMS", n_basis_arr(i_task), n_states_arr(i_task), n_k_points_task_arr(i_task)
  !! DBG      call localorb_info(info_str)
  !! DBG    end do

  !! DBG   n_k_points = sum( n_k_points_task_arr(:) )
  !! DBG   write(info_str,'(2X,A,3I6)') " Z2: We have n_k_points : ", n_k_points
  !! DBG   call localorb_info(info_str)
  !! DBG   write(info_str,'(2X,A,3I6)') " Z2: We have n_tasks    : ", n_tasks
  !! DBG   call localorb_info(info_str)

  !! DBG   i_k = 1
  !! DBG   do i_k_point = 1, n_k_points, 1
  !! DBG     if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
  !! DBG       KS_tmp(:,:,:) = KS_eigenvector_soc_perturbed(:,:,:,i_k)
  !! DBG       i_k = i_k + 1
  !! DBG     else
  !! DBG       KS_tmp(:,:,:) = (0.0d0,0.0d0)
  !! DBG     end if
  !! DBG     call sync_vector_complex( KS_tmp , maxval(n_basis_arr,1)*maxval(n_states_arr,1 )*2 )

  !! DBG     do i_state=1,maxval(n_states_arr,1)
  !! DBG       prod = (0.0d0,0.0d0)
  !! DBG       do i_spin=1,2
  !! DBG         do i_basis=1,maxval(n_basis_arr,1)
  !! DBG           prod = prod + CONJG( KS_tmp(i_basis,i_state,i_spin) ) * KS_tmp(i_basis,i_state,i_spin)
  !! DBG           !if ( Z2_amp_of_complex_nr( KS_tmp(i_basis,i_state,i_spin) ) .gt. 1.0d-3) then
  !! DBG           !  write(info_str,'(2X,A,4I6,2F20.5)') " Z2: KS SOC PER on TASK", i_k_point, i_state, i_spin, i_basis, KS_tmp(i_basis,i_state,i_spin)
  !! DBG           !  call localorb_info(info_str)
  !! DBG           !end if
  !! DBG         end do
  !! DBG       end do
  !! DBG       write(info_str,'(2X,A,3I6,2F20.5)') " Z2: KS PROD", i_k_point, i_state, i_spin, prod
  !! DBG       call localorb_info(info_str)
  !! DBG     end do
  !! DBG   end do

  !! DBG   i_k = 1
  !! DBG   do i_k_point = 1, n_k_points, 1
  !! DBG     if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
  !! DBG       write(info_str,'(A,I3.3,A)') "kp_",i_k_point,".dat"
  !! DBG       open(1313,FILE=trim(info_str))
  !! DBG       write(1313,'(2X,A,4I6,2F20.5)') " Z2: ", i_k, myid
  !! DBG       do i_state=1,maxval(n_states_arr,1)
  !! DBG         do i_spin=1,2
  !! DBG           prod = (0.0d0,0.0d0)
  !! DBG           do i_basis=1,maxval(n_basis_arr,1)
  !! DBG             prod = prod + CONJG( KS_eigenvector_soc_perturbed(i_basis,i_state,i_spin,i_k)) * KS_eigenvector_soc_perturbed(i_basis,i_state,i_spin,i_k)
  !! DBG           end do
  !! DBG           write(1313,'(2X,A,3I6,2F20.5)') " Z2: KS PROD", i_k_point, i_state, i_spin, prod
  !! DBG           ! do i_basis=1,maxval(n_basis_arr,1)
  !! DBG           !   if ( Z2_amp_of_complex_nr( KS_eigenvector_soc_perturbed(i_basis,i_state,i_spin,i_k) ) .gt. 1.0d-3) then
  !! DBG           !     write(1313,'(2X,A,4I6,4F20.5)') " Z2: KS SOC PER on TASK", i_k_point, i_state, i_spin, i_basis, &
  !! DBG           !       & KS_eigenvector_soc_perturbed(i_basis,i_state,i_spin,i_k), i
  !! DBG           !       & Z2_amp_of_complex_nr( KS_eigenvector_soc_perturbed(i_basis,i_state,i_spin,i_k) ),
  !! DBG           !       & Z2_phase_of_complex_nr( KS_eigenvector_soc_perturbed(i_basis,i_state,i_spin,i_k) )
  !! DBG           !   end if
  !! DBG           ! end do
  !! DBG         end do
  !! DBG       end do
  !! DBG       close(1313)
  !! DBG       i_k = i_k + 1
  !! DBG     end if

  !! DBG   end do

  !! DBG   deallocate( KS_tmp )

  !! DBG end subroutine Z2_print_KS_on_tasks_v2

  !! DBG subroutine Z2_print_KS_on_tasks(n_basis, n_states, n_k_points_task, KS_eigenvector_soc_perturbed )
  !! DBG   use localorb_io
  !! DBG   use mpi_tasks
  !! DBG   use synchronize_mpi_basic
  !! DBG   implicit none
  !! DBG   integer, intent(in)                                                      :: n_basis, n_states, n_k_points_task
  !! DBG   complex*16,intent(in),dimension(n_basis,2*n_states,2,n_k_points_task)    :: KS_eigenvector_soc_perturbed
  !! DBG   ! local
  !! DBG   integer :: i_task,i_state,i_basis,i_k,i_spin
  !! DBG   complex*16,dimension(:),allocatable               :: KS_tmp,KS_tmp_ref
  !! DBG   character*300 :: info_str
  !! DBG   integer,dimension(0:n_tasks-1) :: n_basis_arr, n_states_arr, n_k_points_task_arr
  !! DBG   complex*16   :: prod,semphas
  !! DBG   ! DBG:
  !! DBG   integer :: nt_dbg


  !! DBG   ! Get max and min dimensions first
  !! DBG   ! write(info_str,'(2X,A)') " Z2: Z2_print_KS_on_tasks HERES"
  !! DBG   ! call localorb_info(info_str)
  !! DBG   ! write(info_str,'(2X,A,2I4)') " Z2: myid, n_tasks", myid, n_tasks
  !! DBG   ! call localorb_allinfo(info_str)
  !! DBG   ! call mpi_barrier(mpi_comm_global,nt_dbg)

  !! DBG   n_basis_arr(:) = 0
  !! DBG   n_states_arr(:) = 0
  !! DBG   n_k_points_task_arr(:) = 0
  !! DBG   n_basis_arr(myid) = n_basis
  !! DBG   n_states_arr(myid) = n_states
  !! DBG   n_k_points_task_arr(myid) = n_k_points_task
  !! DBG   ! write(info_str,'(2X,A,3I6)') " Z2: Z2_print_KS_on_tasks n_basis_arr bef sync",myid,n_tasks,SIZE(n_basis_arr)
  !! DBG   ! call localorb_allinfo(info_str)
  !! DBG   call sync_integer_vector( n_basis_arr, n_tasks )
  !! DBG   ! write(info_str,'(2X,A)') " Z2: Z2_print_KS_on_tasks n_basis_arr"
  !! DBG   ! call localorb_info(info_str)
  !! DBG   call sync_integer_vector( n_states_arr, n_tasks )
  !! DBG   ! write(info_str,'(2X,A)') " Z2: Z2_print_KS_on_tasks n_tasks"
  !! DBG   ! call localorb_info(info_str)
  !! DBG   call sync_integer_vector( n_k_points_task_arr, n_tasks )
  !! DBG   ! write(info_str,'(2X,A,3I6)') " Z2: MAX_DIMS", maxval(n_basis_arr,1), maxval(n_states_arr,1), maxval(n_k_points_task_arr,1)
  !! DBG   ! call localorb_info(info_str)
  !! DBG   ! write(info_str,'(2X,A,3I6)') " Z2: MIN_DIMS", minval(n_basis_arr,1), minval(n_states_arr,1), minval(n_k_points_task_arr,1)
  !! DBG   ! call localorb_info(info_str)
  !! DBG   allocate( KS_tmp( maxval(n_basis_arr,1) ) )
  !! DBG   allocate( KS_tmp_ref( maxval(n_basis_arr,1) ) )
  !! DBG   KS_tmp(:) = (0.0d0,0.0d0)
  !! DBG   KS_tmp_ref(:) = (0.0d0,0.0d0)
  !! DBG   ! do i_task=0,n_tasks-1
  !! DBG   !   write(info_str,'(2X,A,3I6)') " Z2: DIMS", n_basis_arr(i_task), n_states_arr(i_task), n_k_points_task_arr(i_task)
  !! DBG   !   call localorb_info(info_str)
  !! DBG   ! end do


  !! DBG   do i_task=0,n_tasks-1
  !! DBG     do i_k=1,n_k_points_task_arr(i_task)
  !! DBG       write(info_str,'(2X,A,2I6)') " Z2: LOOP TASK / i_k", i_task, i_k
  !! DBG       call localorb_info(info_str)
  !! DBG       do i_state=1,maxval(n_states_arr,1)
  !! DBG         do i_spin=1,2
  !! DBG           KS_tmp = (0.0d0,0.0d0)
  !! DBG           if (myid == i_task) then
  !! DBG             KS_tmp(1:n_basis_arr(i_task)) = KS_eigenvector_soc_perturbed(1:n_basis_arr(i_task),i_state,i_spin,i_k)
  !! DBG           else
  !! DBG             KS_tmp = (0.0d0,0.0d0)
  !! DBG           end if
  !! DBG           call sync_vector_complex( KS_tmp , maxval(n_basis_arr,1) )

  !! DBG           ! if ( (i_state .eq. 97) .or. (i_state .eq. 98) ) then
  !! DBG           do i_basis=1,maxval(n_basis_arr,1)
  !! DBG             write(info_str,'(2X,A,5I6,4F20.5)') " Z2: KS SOC PER on TASK", i_task, i_basis, i_state, i_spin, i_k, Z2_amp_of_complex_nr(KS_tmp(i_basis)), Z2_phase_of_complex_nr(KS_tmp(i_basis)), KS_tmp(i_basis)
  !! DBG             call localorb_info(info_str)
  !! DBG           end do
  !! DBG           ! end if
  !! DBG           prod    = (0.0d0,0.0d0)
  !! DBG           semphas = (0.0d0,0.0d0)
  !! DBG           do i_basis=1,maxval(n_basis_arr,1)
  !! DBG             prod = prod + CONJG( KS_tmp(i_basis) ) * KS_tmp(i_basis)
  !! DBG             semphas = semphas + ( KS_tmp(i_basis) ) * KS_tmp(i_basis)
  !! DBG           end do
  !! DBG           semphas = semphas / prod
  !! DBG           ! write(info_str,'(2X,A,4I6,6F12.3)') " Z2: KS SOC DOT on TASK", i_task, i_state, i_spin, i_k, prod, semphas , Z2_amp_of_complex_nr( semphas ), Z2_phase_of_complex_nr( semphas )
  !! DBG           ! call localorb_info(info_str)
  !! DBG           ! prod = (0.0d0,0.0d0)
  !! DBG           ! do i_basis=1,maxval(n_basis_arr,1)
  !! DBG           !   prod = prod + CONJG( KS_tmp(i_basis) ) * KS_tmp_ref(i_basis)
  !! DBG           ! end do
  !! DBG           ! write(info_str,'(2X,A,4I6,2F20.5)') " Z2: KS SOC DOTREF on TASK", i_task, i_state, i_spin, i_k, prod
  !! DBG           ! call localorb_info(info_str)
  !! DBG           ! KS_tmp_ref(:) = KS_tmp

  !! DBG         end do
  !! DBG       end do
  !! DBG     end do
  !! DBG   end do
  !! DBG   deallocate( KS_tmp )
  !! DBG   deallocate( KS_tmp_ref )

  !! DBG end subroutine Z2_print_KS_on_tasks

  !! DBG subroutine Z2_print_KS_wf_basis_on_tasks( n_dim, KS_wf )
  !! DBG   use localorb_io
  !! DBG   use mpi_tasks
  !! DBG   use synchronize_mpi_basic
  !! DBG   implicit none
  !! DBG   integer, intent(in)                             :: n_dim
  !! DBG   complex*16,intent(in),dimension(n_dim,n_dim)    :: KS_wf
  !! DBG   ! local
  !! DBG   integer :: i_task,i_state,j_state,r_state
  !! DBG   complex*16,dimension(:,:),allocatable               :: KS_tmp
  !! DBG   character*300 :: info_str
  !! DBG   integer,dimension(0:n_tasks-1) :: n_dim_arr
  !! DBG   complex*16   :: prod

  !! DBG   write(info_str,'(2X,A)') " Z2: CHECKING ORTHONORMALITY OF KS_WF_BASIS"
  !! DBG   call localorb_info(info_str)

  !! DBG   ! Get max and min dimensions first
  !! DBG   n_dim_arr(:) = 0
  !! DBG   n_dim_arr(myid) = n_dim
  !! DBG   call sync_integer_vector( n_dim_arr, n_tasks )
  !! DBG   write(info_str,'(2X,A,I6)') " Z2: MAX_DIMS", maxval(n_dim_arr,1)
  !! DBG   call localorb_info(info_str)
  !! DBG   write(info_str,'(2X,A,I6)') " Z2: MIN_DIMS", minval(n_dim_arr,1)
  !! DBG   call localorb_info(info_str)
  !! DBG   allocate( KS_tmp( maxval(n_dim_arr,1) , maxval(n_dim_arr,1)  ) )
  !! DBG   KS_tmp(:,:) = (0.0d0,0.0d0)
  !! DBG   do i_task=0,n_tasks-1
  !! DBG     write(info_str,'(2X,A,I6)') " Z2: DIMS", n_dim_arr(i_task)
  !! DBG     call localorb_info(info_str)
  !! DBG   end do


  !! DBG   do i_task=0,n_tasks-1
  !! DBG     write(info_str,'(2X,A,I6)') " Z2: LOOP", i_task
  !! DBG     call localorb_info(info_str)
  !! DBG     KS_tmp(:,:) = (0.0d0,0.0d0)
  !! DBG     if (myid == i_task) then
  !! DBG       KS_tmp(1:n_dim_arr(i_task),1:n_dim_arr(i_task)) = KS_wf( 1:n_dim_arr(i_task),1:n_dim_arr(i_task) )
  !! DBG     else
  !! DBG       KS_tmp(:,:) = (0.0d0,0.0d0)
  !! DBG     end if
  !! DBG     call sync_matrix_complex( KS_tmp , maxval(n_dim_arr,1) , maxval(n_dim_arr,1) )

  !! DBG     !j_state=97
  !! DBG     !i_state=49
  !! DBG     !write(info_str,'(2X,A,3I6,4F20.5)') " Z2: KS_WF_BASIS SOC WF PER on TASK", i_task, i_state, j_state, KS_tmp(i_state,j_state), Z2_amp_of_complex_nr( KS_tmp(i_state,j_state) ), Z2_phase_of_complex_nr( KS_tmp(i_state,j_state) )
  !! DBG     !call localorb_info(info_str)
  !! DBG     !i_state=113
  !! DBG     !write(info_str,'(2X,A,3I6,4F20.5)') " Z2: KS_WF_BASIS SOC WF PER on TASK", i_task, i_state, j_state, KS_tmp(i_state,j_state), Z2_amp_of_complex_nr( KS_tmp(i_state,j_state) ), Z2_phase_of_complex_nr( KS_tmp(i_state,j_state) )
  !! DBG     !call localorb_info(info_str)
  !! DBG     !j_state=98
  !! DBG     !i_state=49
  !! DBG     !write(info_str,'(2X,A,3I6,4F20.5)') " Z2: KS_WF_BASIS SOC WF PER on TASK", i_task, i_state, j_state, KS_tmp(i_state,j_state), Z2_amp_of_complex_nr( KS_tmp(i_state,j_state) ), Z2_phase_of_complex_nr( KS_tmp(i_state,j_state) )
  !! DBG     !call localorb_info(info_str)
  !! DBG     !i_state=113
  !! DBG     !write(info_str,'(2X,A,3I6,4F20.5)') " Z2: KS_WF_BASIS SOC WF PER on TASK", i_task, i_state, j_state, KS_tmp(i_state,j_state), Z2_amp_of_complex_nr( KS_tmp(i_state,j_state) ), Z2_phase_of_complex_nr( KS_tmp(i_state,j_state) )
  !! DBG     !call localorb_info(info_str)



  !! DBG     !do j_state=97,98
  !! DBG     !  do i_state=1,n_dim
  !! DBG     !    write(info_str,'(2X,A,3I6,4F20.3)') " Z2: KS_WF_BASIS SOC WF PER on TASK", i_task, i_state, j_state, KS_tmp(i_state,j_state), Z2_amp_of_complex_nr( KS_tmp(i_state,j_state) ), Z2_phase_of_complex_nr( KS_tmp(i_state,j_state) )
  !! DBG     !    call localorb_info(info_str)
  !! DBG     !  end do
  !! DBG     !  ! write(info_str,'(2X,A,3I6,4F20.5)') " Z2: KS_WF_BASIS SOC WF PER on TASK", i_task, i_state, i_state, Z2_amp_of_complex_nr( KS_tmp(i_state,i_state) ), Z2_phase_of_complex_nr( KS_tmp(i_state,i_state) ), KS_tmp(i_state,i_state)
  !! DBG     !  ! call localorb_info(info_str)
  !! DBG     !end do
  !! DBG     do i_state=1,n_dim
  !! DBG       do j_state=1,n_dim
  !! DBG         write(info_str,'(2X,A,3I6,4F20.5)') " Z2: KS_WF_BASIS SOC WF PER on TASK", i_task, i_state, j_state, KS_tmp(i_state,j_state) !, Z2_amp_of_complex_nr( KS_tmp(i_state,j_state) ), Z2_phase_of_complex_nr( KS_tmp(i_state,j_state) )
  !! DBG         call localorb_info(info_str)
  !! DBG       end do
  !! DBG       ! write(info_str,'(2X,A,3I6,4F20.5)') " Z2: KS_WF_BASIS SOC WF PER on TASK", i_task, i_state, i_state, Z2_amp_of_complex_nr( KS_tmp(i_state,i_state) ), Z2_phase_of_complex_nr( KS_tmp(i_state,i_state) ), KS_tmp(i_state,i_state)
  !! DBG       ! call localorb_info(info_str)
  !! DBG     end do


  !! DBG     ! do i_state=1,maxval(n_dim_arr,1)
  !! DBG     !   prod = (0.0d0,0.0d0)
  !! DBG     !   do j_state=1,maxval(n_dim_arr,1)
  !! DBG     !     prod = prod + CONJG( KS_tmp(i_state,j_state) ) * KS_tmp(i_state,j_state)
  !! DBG     !   end do
  !! DBG     !   if ( ( abs(real(prod-1.0d0)) .gt. 1e-8) .or.  ( abs(aimag(prod)) .gt. 1e-8) ) then
  !! DBG     !     write(info_str,'(2X,A,2I6,2F20.5)') " Z2: KS_WF_BASIS SOC DOT on TASK", i_task, i_state, prod
  !! DBG     !     call localorb_info(info_str)
  !! DBG     !   end if
  !! DBG     !   do r_state=i_state+1,maxval(n_dim_arr,1)
  !! DBG     !     prod = (0.0d0,0.0d0)
  !! DBG     !     do j_state=1,maxval(n_dim_arr,1)
  !! DBG     !       prod = prod + CONJG( KS_tmp(i_state,j_state) ) * KS_tmp(r_state,j_state)
  !! DBG     !     end do
  !! DBG     !     if ( ( abs(real(prod)) .gt. 1e-8) .or.  ( abs(aimag(prod)) .gt. 1e-8) ) then
  !! DBG     !       write(info_str,'(2X,A,3I6,2F20.5)') " Z2: KS_WF_BASIS SOC DOT I/J on TASK", i_task, i_state,r_state, prod
  !! DBG     !       call localorb_info(info_str)
  !! DBG     !     end if
  !! DBG     !   end do
  !! DBG     ! end do
  !! DBG   end do
  !! DBG   deallocate( KS_tmp )

  !! DBG end subroutine Z2_print_KS_wf_basis_on_tasks

  !! DBG subroutine Z2_bubble_sort_ev( n_dim,  KS_eval, n_basis, KS_evec )
  !! DBG   implicit none
  !! DBG   integer,intent(in) :: n_dim, n_basis
  !! DBG   real*8,dimension(n_dim),intent(inout) :: KS_eval
  !! DBG   complex*16,dimension(n_basis,n_dim,2),intent(inout),optional :: KS_evec
  !! DBG   !local
  !! DBG   integer                     :: i_state,j_state
  !! DBG   real*8,dimension(n_dim)     :: KS_eval_tmp
  !! DBG   real*8                      :: dummy
  !! DBG   complex*16,dimension(:,:,:),allocatable :: KS_evec_tmp
  !! DBG   complex*16,dimension(:,:),allocatable   :: dummy_cmplx

  !! DBG   KS_eval_tmp = KS_eval
  !! DBG   if( present( KS_evec ) ) then
  !! DBG     allocate( KS_evec_tmp(n_basis,n_dim,2) )
  !! DBG     allocate( dummy_cmplx(n_basis,2) )
  !! DBG     KS_evec_tmp = KS_evec
  !! DBG     dummy_cmplx(:,:) = (0.0d0,0.0d0)
  !! DBG   else
  !! DBG     allocate( KS_evec_tmp(1,1,1) )
  !! DBG     allocate( dummy_cmplx(1,1) )
  !! DBG     KS_evec_tmp(:,:,:) = (0.0d0,0.0d0)
  !! DBG     dummy_cmplx(:,:) = (0.0d0,0.0d0)
  !! DBG   end if

  !! DBG   do i_state = 1, n_dim, 1
  !! DBG     do j_state = i_state + 1, n_dim, 1
  !! DBG       if (KS_eval_tmp(i_state) > KS_eval_tmp(j_state)) then
  !! DBG         dummy = KS_eval_tmp(i_state)
  !! DBG         dummy_cmplx = KS_evec_tmp(:,i_state,:)

  !! DBG         KS_eval_tmp(i_state) = KS_eval_tmp(j_state)
  !! DBG         KS_evec_tmp(:,i_state,:) = KS_evec_tmp(:,j_state,:)

  !! DBG         KS_eval_tmp(j_state) = dummy
  !! DBG         KS_evec_tmp(:,j_state,:) = dummy_cmplx
  !! DBG       end if
  !! DBG     end do
  !! DBG   end do

  !! DBG   KS_eval = KS_eval_tmp
  !! DBG   if( present( KS_evec ) ) then
  !! DBG     KS_evec = KS_evec_tmp
  !! DBG   end if
  !! DBG   deallocate( KS_evec_tmp )
  !! DBG   deallocate( dummy_cmplx )

  !! DBG end subroutine Z2_bubble_sort_ev

  !! DBG function Z2_amp_of_complex_nr( z )
  !! DBG   implicit none
  !! DBG   complex*16 :: z
  !! DBG   real*8     :: Z2_amp_of_complex_nr
  !! DBG 
  !! DBG   Z2_amp_of_complex_nr = sqrt(dble( CONJG(z) * z ))
  !! DBG 
  !! DBG end function Z2_amp_of_complex_nr
  !! DBG 
  !! DBG function Z2_phase_of_complex_nr( z )
  !! DBG   use constants
  !! DBG   implicit none
  !! DBG   complex*16 :: z
  !! DBG   real*8     :: Z2_phase_of_complex_nr
  !! DBG   real*8     :: amp
  !! DBG 
  !! DBG   amp = Z2_amp_of_complex_nr(z)
  !! DBG   if (amp .gt. 1e-8) then
  !! DBG     Z2_phase_of_complex_nr = datan2( dimag(z), real( z ) ) / pi *180.0d0
  !! DBG   else
  !! DBG     Z2_phase_of_complex_nr = 0.0d0
  !! DBG   end if
  !! DBG 
  !! DBG 
  !! DBG end function Z2_phase_of_complex_nr
  !! DBG 
  !! DBG subroutine Z2_write_matrix_ndim_sq_to_file( string,  n_dim, KS_wf, i_band )
  !! DBG   use localorb_io
  !! DBG   use mpi_tasks
  !! DBG   use synchronize_mpi_basic
  !! DBG   implicit none
  !! DBG   character(len=*),intent(in)                     :: string
  !! DBG   integer, intent(in)                             :: n_dim,i_band
  !! DBG   complex*16,intent(in),dimension(n_dim,n_dim)    :: KS_wf
  !! DBG   ! local
  !! DBG   integer :: i_task,i_state,j_state,r_state
  !! DBG   complex*16,dimension(:,:),allocatable               :: KS_tmp
  !! DBG   character*300 :: info_str
  !! DBG   integer,dimension(0:n_tasks-1) :: n_dim_arr
  !! DBG   complex*16   :: prod
  !! DBG 
  !! DBG   ! Get max and min dimensions first
  !! DBG   n_dim_arr(:) = 0
  !! DBG   n_dim_arr(myid) = n_dim
  !! DBG   call sync_integer_vector( n_dim_arr, n_tasks )
  !! DBG   write(info_str,'(2X,A,I6)') " Z2: MAX_DIMS", maxval(n_dim_arr,1)
  !! DBG   call localorb_info(info_str)
  !! DBG   write(info_str,'(2X,A,I6)') " Z2: MIN_DIMS", minval(n_dim_arr,1)
  !! DBG   call localorb_info(info_str)
  !! DBG   allocate( KS_tmp( maxval(n_dim_arr,1) , maxval(n_dim_arr,1)  ) )
  !! DBG   KS_tmp(:,:) = (0.0d0,0.0d0)
  !! DBG   do i_task=0,n_tasks-1
  !! DBG     write(info_str,'(2X,A,I6)') " Z2: DIMS", n_dim_arr(i_task)
  !! DBG     call localorb_info(info_str)
  !! DBG   end do
  !! DBG 
  !! DBG   do i_task=0,n_tasks-1
  !! DBG     write(info_str,'(2X,A,I6)') " Z2: LOOP", i_task
  !! DBG     call localorb_info(info_str)
  !! DBG     KS_tmp(:,:) = (0.0d0,0.0d0)
  !! DBG     if (myid == i_task) then
  !! DBG       KS_tmp(1:n_dim_arr(i_task),1:n_dim_arr(i_task)) = KS_wf( 1:n_dim_arr(i_task),1:n_dim_arr(i_task) )
  !! DBG     else
  !! DBG       KS_tmp(:,:) = (0.0d0,0.0d0)
  !! DBG     end if
  !! DBG     call sync_matrix_complex( KS_tmp , maxval(n_dim_arr,1) , maxval(n_dim_arr,1) )
  !! DBG 
  !! DBG     write(info_str,'(A,A,I2.2,A,I4.4,A)') trim(string),"_ktask_",i_task,'_i_band_',i_band,'.dat'
  !! DBG     open(1313,FILE=trim(info_str),STATUS='replace',ACTION='write')
  !! DBG     do i_state=1,n_dim
  !! DBG       do j_state=1,n_dim
  !! DBG         write(1313,'(3I6,2EN22.12)') i_task, i_state, j_state, KS_tmp(i_state,j_state)
  !! DBG       end do
  !! DBG     end do
  !! DBG     close(1313)
  !! DBG 
  !! DBG   end do
  !! DBG 
  !! DBG   deallocate( KS_tmp )
  !! DBG 
  !! DBG end subroutine Z2_write_matrix_ndim_sq_to_file
  !! DBG 
  !! DBG subroutine Z2_read_matrix_ndim_sq_to_file( string, n_dim, KS_wf, i_band )
  !! DBG   use localorb_io
  !! DBG   use mpi_tasks
  !! DBG   use synchronize_mpi_basic
  !! DBG   implicit none
  !! DBG   character(len=*),intent(in)                      :: string
  !! DBG   integer, intent(in)                              :: n_dim,i_band
  !! DBG   complex*16,intent(out),dimension(n_dim,n_dim)    :: KS_wf
  !! DBG   ! local
  !! DBG   integer :: i_task,i_state,j_state,i_task_tmp,i_state_tmp,j_state_tmp
  !! DBG   complex*16,dimension(:,:),allocatable               :: KS_tmp
  !! DBG   character*300 :: info_str
  !! DBG   integer,dimension(0:n_tasks-1) :: n_dim_arr
  !! DBG   complex*16   :: prod
  !! DBG 
  !! DBG   ! Get max and min dimensions first
  !! DBG   n_dim_arr(:) = 0
  !! DBG   n_dim_arr(myid) = n_dim
  !! DBG   call sync_integer_vector( n_dim_arr, n_tasks )
  !! DBG   write(info_str,'(2X,A,I6)') " Z2: MAX_DIMS", maxval(n_dim_arr,1)
  !! DBG   call localorb_info(info_str)
  !! DBG   write(info_str,'(2X,A,I6)') " Z2: MIN_DIMS", minval(n_dim_arr,1)
  !! DBG   call localorb_info(info_str)
  !! DBG   allocate( KS_tmp( maxval(n_dim_arr,1) , maxval(n_dim_arr,1)  ) )
  !! DBG   KS_tmp(:,:) = (0.0d0,0.0d0)
  !! DBG   KS_wf(:,:)  = (0.0d0,0.0d0)
  !! DBG   do i_task=0,n_tasks-1
  !! DBG     write(info_str,'(2X,A,I6)') " Z2: DIMS", n_dim_arr(i_task)
  !! DBG     call localorb_info(info_str)
  !! DBG   end do
  !! DBG 
  !! DBG   do i_task=0,n_tasks-1
  !! DBG     write(info_str,'(2X,A,I6)') " Z2: LOOP", i_task
  !! DBG     call localorb_info(info_str)
  !! DBG 
  !! DBG     if (myid == i_task) then
  !! DBG       write(info_str,'(A,A,I2.2,A,I4.4,A)') trim(string),"_ktask_",i_task,'_i_band_',i_band,'.dat'
  !! DBG       open(1313,FILE=trim(info_str),STATUS='OLD',ACTION='read')
  !! DBG       KS_tmp(:,:) = (0.0d0,0.0d0)
  !! DBG       do i_state=1,n_dim_arr(i_task)
  !! DBG         do j_state=1,n_dim_arr(i_task)
  !! DBG           read(1313,'(3I6,2EN22.12)') i_task_tmp, i_state_tmp, j_state_tmp, KS_tmp(i_state,j_state)
  !! DBG           if ( ( i_task_tmp .ne. i_task) .or. (i_state_tmp .ne. i_state) .or. (j_state_tmp .ne. j_state) ) then
  !! DBG             call aims_stop("SHIT")
  !! DBG           end if
  !! DBG         end do
  !! DBG       end do
  !! DBG       close(1313)
  !! DBG       KS_wf( 1:n_dim_arr(i_task),1:n_dim_arr(i_task) ) = KS_tmp(1:n_dim_arr(i_task),1:n_dim_arr(i_task))
  !! DBG     end if
  !! DBG 
  !! DBG   end do
  !! DBG 
  !! DBG   deallocate( KS_tmp )
  !! DBG 
  !! DBG end subroutine Z2_read_matrix_ndim_sq_to_file
  !! DBG 
  !! DBG subroutine Z2_check_orthogonality( n_dim, KS_wf )
  !! DBG   use localorb_io
  !! DBG   use mpi_tasks
  !! DBG   use synchronize_mpi_basic
  !! DBG   implicit none
  !! DBG   integer, intent(in)                             :: n_dim
  !! DBG   complex*16,intent(in),dimension(n_dim,n_dim)    :: KS_wf
  !! DBG   ! local
  !! DBG   integer :: i_task,i_state,j_state,k_state
  !! DBG   complex*16,dimension(:,:),allocatable               :: KS_tmp
  !! DBG   character*300 :: info_str
  !! DBG   integer,dimension(0:n_tasks-1) :: n_dim_arr
  !! DBG   complex*16   :: prod
  !! DBG 
  !! DBG   ! Get max and min dimensions first
  !! DBG   n_dim_arr(:) = 0
  !! DBG   n_dim_arr(myid) = n_dim
  !! DBG   call sync_integer_vector( n_dim_arr, n_tasks )
  !! DBG   write(info_str,'(2X,A,I6)') " Z2: MAX_DIMS", maxval(n_dim_arr,1)
  !! DBG   call localorb_info(info_str)
  !! DBG   write(info_str,'(2X,A,I6)') " Z2: MIN_DIMS", minval(n_dim_arr,1)
  !! DBG   call localorb_info(info_str)
  !! DBG   allocate( KS_tmp( maxval(n_dim_arr,1) , maxval(n_dim_arr,1)  ) )
  !! DBG   KS_tmp(:,:) = (0.0d0,0.0d0)
  !! DBG   do i_task=0,n_tasks-1
  !! DBG     write(info_str,'(2X,A,I6)') " Z2: DIMS", n_dim_arr(i_task)
  !! DBG     call localorb_info(info_str)
  !! DBG   end do
  !! DBG 
  !! DBG   do i_task=0,n_tasks-1
  !! DBG     write(info_str,'(2X,A,I6)') " Z2: LOOP", i_task
  !! DBG     call localorb_info(info_str)
  !! DBG     KS_tmp(:,:) = (0.0d0,0.0d0)
  !! DBG     if (myid == i_task) then
  !! DBG       KS_tmp(1:n_dim_arr(i_task),1:n_dim_arr(i_task)) = KS_wf( 1:n_dim_arr(i_task),1:n_dim_arr(i_task) )
  !! DBG     else
  !! DBG       KS_tmp(:,:) = (0.0d0,0.0d0)
  !! DBG     end if
  !! DBG     call sync_matrix_complex( KS_tmp , maxval(n_dim_arr,1) , maxval(n_dim_arr,1) )
  !! DBG 
  !! DBG     do i_state=1,n_dim_arr(i_task)
  !! DBG       prod = (0.0d0,0.0d0)
  !! DBG       do j_state=1,n_dim_arr(i_task)
  !! DBG         prod = prod + ( conjg( KS_tmp(j_state,i_state) ) * KS_tmp(j_state,i_state) )
  !! DBG       end do
  !! DBG       if (  ( abs( DBLE(prod) - 1.0d0 ) .gt. 1.0d-12 ) .or.  ( abs( DIMAG(prod) ) .gt. 1.0d-12 ) ) then
  !! DBG         write(info_str,'(2X,A)') " Z2: **** STATE NOT NORMED"
  !! DBG         call localorb_info(info_str)
  !! DBG         write(info_str,'(2X,A,I6,2F20.8)') " Z2: Norm of state i:", i_state, prod
  !! DBG         call localorb_info(info_str)
  !! DBG       end if
  !! DBG     end do
  !! DBG 
  !! DBG     do i_state=1,n_dim_arr(i_task)
  !! DBG       do j_state=1,n_dim_arr(i_task)
  !! DBG         if (j_state .eq. i_state) cycle
  !! DBG         prod = (0.0d0,0.0d0)
  !! DBG         do k_state=1,n_dim_arr(i_task)
  !! DBG           prod = prod + ( conjg( KS_tmp(k_state,i_state) ) * KS_tmp(k_state,j_state) )
  !! DBG         end do
  !! DBG         if (  ( abs( DBLE(prod) ) .gt. 1.0d-12 ) .or.  ( abs( DIMAG(prod) ) .gt. 1.0d-12 ) ) then
  !! DBG           write(info_str,'(2X,A)') " Z2: **** STATES NOT ORTHOGONAL"
  !! DBG           call localorb_info(info_str)
  !! DBG           write(info_str,'(2X,A,2I6,2F20.8)') " Z2: Ovlp  of states i j:", i_state,j_state, prod
  !! DBG           call localorb_info(info_str)
  !! DBG         end if
  !! DBG       end do
  !! DBG     end do
  !! DBG 
  !! DBG   end do
  !! DBG 
  !! DBG   deallocate( KS_tmp )
  !! DBG 
  !! DBG end subroutine Z2_check_orthogonality
  !! DBG 
  !! DBG subroutine Z2_check_orthogonality_compare( n_dim, KS_wf, KS_wf_read )
  !! DBG   use localorb_io
  !! DBG   use mpi_tasks
  !! DBG   use synchronize_mpi_basic
  !! DBG   implicit none
  !! DBG   integer, intent(in)                             :: n_dim
  !! DBG   complex*16,intent(in),dimension(n_dim,n_dim)    :: KS_wf,KS_wf_read
  !! DBG   ! local
  !! DBG   integer :: i_task,i_state,j_state,k_state
  !! DBG   complex*16,dimension(:,:),allocatable               :: KS_tmp, KS_tmp_read
  !! DBG   character*300 :: info_str
  !! DBG   integer,dimension(0:n_tasks-1) :: n_dim_arr
  !! DBG   complex*16   :: prod
  !! DBG 
  !! DBG   ! Get max and min dimensions first
  !! DBG   n_dim_arr(:) = 0
  !! DBG   n_dim_arr(myid) = n_dim
  !! DBG   call sync_integer_vector( n_dim_arr, n_tasks )
  !! DBG   write(info_str,'(2X,A,I6)') " Z2: MAX_DIMS", maxval(n_dim_arr,1)
  !! DBG   call localorb_info(info_str)
  !! DBG   write(info_str,'(2X,A,I6)') " Z2: MIN_DIMS", minval(n_dim_arr,1)
  !! DBG   call localorb_info(info_str)
  !! DBG   allocate( KS_tmp( maxval(n_dim_arr,1) , maxval(n_dim_arr,1)  ) )
  !! DBG   allocate( KS_tmp_read( maxval(n_dim_arr,1) , maxval(n_dim_arr,1)  ) )
  !! DBG   KS_tmp(:,:) = (0.0d0,0.0d0)
  !! DBG   do i_task=0,n_tasks-1
  !! DBG     write(info_str,'(2X,A,I6)') " Z2: DIMS", n_dim_arr(i_task)
  !! DBG     call localorb_info(info_str)
  !! DBG   end do
  !! DBG 
  !! DBG   do i_task=0,n_tasks-1
  !! DBG     write(info_str,'(2X,A,I6)') " Z2: LOOP", i_task
  !! DBG     call localorb_info(info_str)
  !! DBG     KS_tmp(:,:) = (0.0d0,0.0d0)
  !! DBG     KS_tmp_read(:,:) = (0.0d0,0.0d0)
  !! DBG     if (myid == i_task) then
  !! DBG       KS_tmp     (1:n_dim_arr(i_task),1:n_dim_arr(i_task)) = KS_wf     ( 1:n_dim_arr(i_task),1:n_dim_arr(i_task) )
  !! DBG       KS_tmp_read(1:n_dim_arr(i_task),1:n_dim_arr(i_task)) = KS_wf_read( 1:n_dim_arr(i_task),1:n_dim_arr(i_task) )
  !! DBG     else
  !! DBG       KS_tmp     (:,:) = (0.0d0,0.0d0)
  !! DBG       KS_tmp_read(:,:) = (0.0d0,0.0d0)
  !! DBG     end if
  !! DBG     call sync_matrix_complex( KS_tmp      , maxval(n_dim_arr,1) , maxval(n_dim_arr,1) )
  !! DBG     call sync_matrix_complex( KS_tmp_read , maxval(n_dim_arr,1) , maxval(n_dim_arr,1) )
  !! DBG 
  !! DBG     do i_state=1,n_dim_arr(i_task)
  !! DBG       prod = (0.0d0,0.0d0)
  !! DBG       do j_state=1,n_dim_arr(i_task)
  !! DBG         prod = prod + ( conjg( KS_tmp(j_state,i_state) ) * KS_tmp_read(j_state,i_state) )
  !! DBG       end do
  !! DBG       if (  ( abs( DBLE(prod) - 1.0d0 ) .gt. 0.1 ) .or.  ( abs( DIMAG(prod) ) .gt. 0.1 ) ) then
  !! DBG         write(info_str,'(2X,A)') " Z2: **** STATE NOT NORMED"
  !! DBG         call localorb_info(info_str)
  !! DBG         write(info_str,'(2X,A,I6,4F20.8)') " Z2: Norm of state i:", i_state, prod, Z2_amp_of_complex_nr( prod ), Z2_phase_of_complex_nr( prod )
  !! DBG         call localorb_info(info_str)
  !! DBG       end if
  !! DBG     end do
  !! DBG 
  !! DBG     do i_state=1,n_dim_arr(i_task)
  !! DBG       do j_state=1,n_dim_arr(i_task)
  !! DBG         if (j_state .eq. i_state) cycle
  !! DBG         prod = (0.0d0,0.0d0)
  !! DBG         do k_state=1,n_dim_arr(i_task)
  !! DBG           prod = prod + ( conjg( KS_tmp(k_state,i_state) ) * KS_tmp_read(k_state,j_state) )
  !! DBG         end do
  !! DBG         if (  ( abs( DBLE(prod) ) .gt. 0.1 ) .or.  ( abs( DIMAG(prod) ) .gt. 0.1 ) ) then
  !! DBG           write(info_str,'(2X,A)') " Z2: **** STATES NOT ORTHOGONAL"
  !! DBG           call localorb_info(info_str)
  !! DBG           write(info_str,'(2X,A,2I6,4F20.8)') " Z2: Ovlp  of states i j:", i_state,j_state, prod, Z2_amp_of_complex_nr( prod ), Z2_phase_of_complex_nr( prod )
  !! DBG           call localorb_info(info_str)
  !! DBG         end if
  !! DBG       end do
  !! DBG     end do
  !! DBG 
  !! DBG   end do
  !! DBG 
  !! DBG   deallocate( KS_tmp )
  !! DBG   deallocate( KS_tmp_read )
  !! DBG 
  !! DBG end subroutine Z2_check_orthogonality_compare
  
  
  
  !! DBG !CONVERT FROM NEW INDEX TO OLD:
  !! DBG ! CC: NO LONGER NEED
  !! DBG subroutine Z2_convert_old_to_new_index( n_basis_soc_coll, n_saved_states_soc, KS_WITHOUT_SPIN , n_basis, KS_WITH_SPIN )
  !! DBG   implicit none
  !! DBG   integer, intent(in) :: n_basis_soc_coll, n_saved_states_soc, n_basis
  !! DBG   complex*16,dimension(n_basis_soc_coll,n_saved_states_soc,1),intent(in) :: KS_WITHOUT_SPIN
  !! DBG   complex*16,dimension(n_basis,n_saved_states_soc,2),intent(out)         :: KS_WITH_SPIN
  !! DBG   ! Local:
  !! DBG   integer :: i_spin, basis_offset, i_state, i_basis_1
  !! DBG   do i_spin = 1, 2, 1
  !! DBG     if (i_spin .eq. 1) then
  !! DBG       ! Spin-up components are requested
  !! DBG       basis_offset = 0
  !! DBG     else
  !! DBG       ! Spin-dn components are requested
  !! DBG       basis_offset = n_basis_soc_coll/2
  !! DBG     end if
  !! DBG 
  !! DBG     do i_state = 1, n_saved_states_soc, 1
  !! DBG       do i_basis_1 = 1, n_basis_soc_coll/2
  !! DBG         KS_WITH_SPIN(i_basis_1, i_state, i_spin) = KS_WITHOUT_SPIN(basis_offset+i_basis_1, i_state,1)
  !! DBG       end do
  !! DBG     end do
  !! DBG   end do
  !! DBG end subroutine Z2_convert_old_to_new_index

END MODULE WannierCenters
