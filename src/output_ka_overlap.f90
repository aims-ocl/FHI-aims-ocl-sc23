!****s* FHI-aims/output_ka_overlap
!  NAME
!    output_ka_overlap
!  SYNOPSIS

subroutine output_ka_overlap(overlap_matrix)

!  PURPOSE
!  Subroutine output_ka_overlap writes out all basis function properties, 
!  and the (real-valued) overlap matrix in a format compatible with
!  "aitranss" module (Karlsruhe transport code).
!
!  USES

   use basis
   use dimensions
   use localorb_io
   use runtime_choices, only: out_matrices_elsi

   implicit none

!  ARGUMENTS
!  imported variables

   real*8 overlap_matrix(n_basis*(n_basis+1)/2)

!  INPUTS
!   o overlap_matrix -- overlap matrix
!  OUTPUT
!    none
!  AUTHOR
!    Alexej Bagrets: Institute of Nanotechnology (INT) 
!                    Karlsruhe Institute of Technology (KIT)
!  HISTORY
!    Development version, FHI-aims (2012).
!  SOURCE
!
!  local variables

   real*8 :: output_element(n_basis)
   character*18 :: file_name(n_spin)
   character*128 :: info_str

!  counters
   integer :: i_basis
   integer :: i_basis_1
   integer :: i_basis_2
   integer :: i_index
   integer :: i_spin
   integer :: i_fn
   integer :: n
   integer :: nn
   integer :: j
   integer :: mincol
   integer :: maxcol

!  external files
   integer :: basisfile = 50
   integer :: omatfile = 51

!  begin work

   write (info_str,'(2X,A)') &
      'Writing basis function properties.'
   call localorb_info(info_str,use_unit,'(A)')
   write (info_str,'(2X,I5,A)') n_basis,' basis functions used.'
   call localorb_info(info_str,use_unit,'(A)')

   open(basisfile,file='basis-indices.out')

   write(info_str,'(a)') '$basis.indices.aims' 
   call localorb_info(info_str,basisfile,'(A)')
   write(info_str,'(a,i6)') '# amount of basis functions used: ',n_basis 
   call localorb_info(info_str,basisfile,'(A)')
   write(info_str,'(A5,1X,A8,1X,A3,1X,A3,1X,A3,1X,A3)') &
      '# fn.', '  type  ', 'at.', 'n', 'l', 'm'
   call localorb_info(info_str,basisfile,'(A)')

   do i_basis = 1,n_basis
      i_fn = basis_fn(i_basis)
      write(info_str,'(I5,1X,A8,1X,I3,1X,I3,1X,I3,1X,I3)') &
         i_basis,basisfn_type(i_fn),basis_atom(i_basis),basisfn_n(i_fn),&
         basis_l(i_basis),basis_m(i_basis)
      call localorb_info(info_str,basisfile,'(A)')
   enddo

   write(info_str,'(a)') '$end' 
   call localorb_info(info_str,basisfile,'(A)')

   close(basisfile)

   if(.not. out_matrices_elsi) then
      open(omatfile,file='omat.aims')

!     write(info_str,fmt='(a,i6)') '$omat.aims   format(4d26.20)   nsaos=',n_basis
      write(info_str,fmt='(a,i6)') '$omat.aims   format(4d20.14)   nsaos=',n_basis
      call localorb_info(info_str,omatfile,'(a)')

      do n = 1, n_basis
         write(info_str,fmt='(a,i5)') '#column ',n
         call localorb_info(info_str,omatfile,'(a)')

         nn = n*(n-1)/2
         maxcol = 0
         do while (maxcol<n)
            mincol = maxcol+1
            maxcol = min(maxcol+4,n)
!           write(info_str,fmt='(4d26.20)') (overlap_matrix(j + nn),j=mincol,maxcol)
            write(info_str,fmt='(4d20.14)') (overlap_matrix(j + nn),j=mincol,maxcol)
            call localorb_info(info_str,omatfile,'(a)')
         enddo
      enddo

      write(info_str,fmt='(a)') '$end'
      call localorb_info(info_str,omatfile,'(a)')
      close(omatfile)
   endif

end subroutine output_ka_overlap
      
!----------------------------------------------------------------------
!******	
