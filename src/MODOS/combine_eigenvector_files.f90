!combine eigenvector files 'eigenvec.out***' into one file 'eigenvec.out'
  program combine_eigenvector_files
  implicit none
  character*40 :: filename, filename_new
  integer :: IRECL, iline
  integer :: n_basis, n_states, n_spin, n_k_points, n_atoms
  integer :: n_basis2, n_spin2, n_k_points2
  integer, allocatable :: basis_atom(:)

  real*8, allocatable :: KS_eigenvalue(:,:,:)
  complex*16, allocatable :: KS_eigenvector(:,:,:,:)
  integer :: i_state, i_spin, i_k_point
  integer :: i, itmp(5), fileID
  real(8) :: rtmp
  character*3 c_num
  integer :: LEN, itmp2

  INQUIRE (IOLENGTH=LEN) itmp2

  filename = 'eigenvec.out'
  write(filename_new,'(A)') trim(adjustl(filename))//'001'

! read some preminary data (n_basis, n_states, n_spin, n_k_points, n_atoms) from the first line
  IRECL = 5*LEN
  open (UNIT=50, FILE=filename_new, STATUS='OLD', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='READ', RECL=IRECL)
  read(50, REC=1) n_basis, n_states, n_spin, n_k_points, n_atoms
  close(50)

  !output the information
  write(*,"('n_basis    = ', I6)") n_basis
  write(*,"('n_states   = ', I6)") n_states
  write(*,"('n_spin     = ', I6)") n_spin
  write(*,"('n_k_points = ', I6)") n_k_points
  write(*,"('n_atoms    = ', I6)") n_atoms


  allocate(basis_atom(n_basis))
  allocate(KS_eigenvalue(n_states, n_spin, n_k_points))
  allocate(KS_eigenvector(n_basis, n_states, n_spin, n_k_points))

  do i_k_point = 1, n_k_points

     !construct file name:
     if (i_k_point < 10) then
        write(c_num,'("00",I1)') i_k_point
     else if (i_k_point < 100) then
        write(c_num,'("0",I2)') i_k_point
     else if (i_k_point < 1000) then
        write(c_num,'(I3)') i_k_point
     else
        write(*,*)"Error: too many files to output!"
        stop
     end if
     write(filename_new,'(A)') trim(adjustl(filename))//c_num

     !each line (except line 1 and 2) contains: 1 real number and n_basis complex numbers
     IRECL = (2 + n_basis*4)*LEN
     fileID = 50 + i_k_point
     open (UNIT=fileID, FILE=filename_new, STATUS='Old', FORM='UNFORMATTED', ACCESS='DIRECT', ACTION='READ', RECL=IRECL)
           
     !output the two heading lines
     iline = 1
     read(fileID, REC=iline) itmp(1:5)
     if (itmp(1)/=n_basis .or.  itmp(2)/=n_states .or.  itmp(3)/=n_spin .or.  itmp(4)/=n_k_points .or.  itmp(5)/=n_atoms) then
        write(*,"('Error: the files have different headlines.')")
        write(*,"(5I6)") itmp(1:5)
        stop
     end if
     iline = iline + 1
     read(fileID, REC=iline) basis_atom(1:n_basis)
     
     !read the bulk data: eigenvector(:, :, :, i_k_point)
     do i_spin = 1, n_spin
        do i_state = 1, n_states
           iline = iline + 1
           read(fileID, REC=iline) KS_eigenvalue(i_state, i_spin, i_k_point), KS_eigenvector(1:n_basis, i_state, i_spin, i_k_point)
        end do
     end do
     
     close(fileID)

  end do


  write(*,*)
  write(*,"('Successfully read eigenvectors!')")

  !output all the data into: 'eigenvec.out'
  filename = 'eigenvec.out'
  IRECL = (2 + n_basis*4)*LEN
  open (UNIT=50, FILE=filename, STATUS='REPLACE', FORM='UNFORMATTED', ACCESS='DIRECT', RECL=IRECL)
  
  iline = 1
  write(50, REC=iline) n_basis, n_states, n_spin, n_k_points, n_atoms
  iline = iline + 1
  write(50, REC=iline) basis_atom(1:n_basis)
  
  !write the bulk data
  do i_k_point = 1, n_k_points
     do i_spin = 1, n_spin, 1
        do i_state = 1, n_states, 1
           iline = iline + 1
           write(50, REC=iline) KS_eigenvalue(i_state, i_spin, i_k_point), KS_eigenvector(1:n_basis, i_state, i_spin, i_k_point)
        end do
     end do
  end do

  close(50)
  write(*,*)
  write(*,"('Successfully combine the eigenvector files in to one file: ', A)") filename


  stop
  end program combine_eigenvector_files



  
