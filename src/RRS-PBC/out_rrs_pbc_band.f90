!****s* FHI-aims/RRS-PBC/out_rrs_pbc_band()
!  NAME
!    out_rrs_pbc_band
!  SYNOPSIS

    subroutine out_rrs_pbc_band()

!  PURPOSE
!  This routine implement the actual RRS-PBC calculation, and can only be called
!  after parse_rrs_pbc()
!
!  USES

      use localorb_io
      use constants
      use dimensions
      use physics
      use geometry
      use numerical_utilities
      use basis
      use runtime_choices
      use mpi_tasks
      implicit none

!  ARGUMENTS

!  INPUTS
!    none
!  OUTPUT
!    none
!  COPYRIGHT
!   
!   
!   
!  HISTORY
!   
!  SOURCE
!

    ! imported variables

    ! local variables
    character*132 :: rrs_pbc_info, file_name, info_str
    integer :: i, j, o, p
    integer :: i_state
    integer :: info
    real*8, dimension(3)   :: tmp_k_vec
    real*8, dimension(3)   :: tmp_d_vec
    logical,save  :: rrs_pbc_plot_band_allocated = .false.

    real*8  :: HVBM, LCBM
    HVBM = -100.0d0
    LCBM = 100.0d0

    write(info_str,'()') 
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(4X,A)') "|------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(4X,A)') "| Writing the requested band structure output:"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(4X,A)') "|------------------------------------------------------------------"
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'()') 
    call localorb_info(info_str,use_unit,'(A)')
    do i = 1, rrs_pbc_n_plot_band, 1
        ! Print out several useful information about this k point
        if(i < 10) then
           write(file_name, '(A7,I1,A4)') 'band100', i,'.out'
        else if(i < 100)then
           write(file_name, '(A6,I2,A4)') 'band10', i,'.out'
        else if(i < 1000)then
           write(file_name, '(A5,I3,A4)') 'band1', i,'.out'
        else
           write(use_unit,*) 'Band output error: automatic file name does not work with more than 999 bands!'
           stop
        end if
        open(88,file=file_name)
        if(n_spin > 1) then
            if(i < 10) then
               write(file_name, '(A7,I1,A4)') 'band200', i,'.out'
            else if(i < 100)then
               write(file_name, '(A6,I2,A4)') 'band20', i,'.out'
            else if(i < 1000)then
               write(file_name, '(A5,I3,A4)') 'band2', i,'.out'
            else
               write(use_unit,*) 'Band output error: automatic file name does not work with more than 999 bands!'
               stop
            end if
            open(89,file=file_name)
        endif

        write(info_str,'()') 
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'(4X,A,I4,A,I4,A)') "|Treating all ",rrs_pbc_n_points_in_band(i),&
            " k-points in band plot segment #", i, ":"
        call localorb_info(info_str,use_unit,'(A)')
        write(info_str,'()') 
        call localorb_info(info_str,use_unit,'(A)')

        do j = 1, rrs_pbc_n_points_in_band(i), 1
            tmp_k_vec(:) = rrs_pbc_band_begin(i,:) &
                + real(j-1)/real(rrs_pbc_n_points_in_band(i)-1) &
                  *( rrs_pbc_band_end(i,:) - rrs_pbc_band_begin(i,:))

            !write(info_str,'(4X,A,3(F8.4))') '| RRS-PBC :: getting info. of the k point of (',&
            !    (tmp_k_vec(p), p=1,3),')'
            !!call localorb_info(info_str)
            !write(use_unit,*) info_str

            call build_rrs_pbc_k_matrices(tmp_k_vec)
            call calculate_rrs_pbc_k_matrices(0)

            write(88,'(I4,2X,3F15.7)',ADVANCE='NO') j, tmp_k_vec(1), tmp_k_vec(2), tmp_k_vec(3)
            do  i_state = 1,  rrs_pbc_n_center_basis, 1
               write(88,'(F12.5,F15.5)',ADVANCE='NO') rrs_pbc_occ_num(i_state,1), &! *n_k_points, &
                    rrs_pbc_band_info(i_state,1)* hartree
            end do
            write(88,'()') 

            if(n_spin ==2)then
               write(89,'(I4,2X,3F15.7)',ADVANCE='NO') j, tmp_k_vec(1), tmp_k_vec(2), tmp_k_vec(3)
               do  i_state = 1,  rrs_pbc_n_center_basis, 1
                  write(88,'(F12.5,F15.5)',ADVANCE='NO') rrs_pbc_occ_num(i_state,2), &! *n_k_points, &
                       rrs_pbc_band_info(i_state,2)* hartree
               end do
               write(89,'()') 
            end if
        enddo
        close(88)
        if(n_spin==2) close(89)
    enddo

    end subroutine

