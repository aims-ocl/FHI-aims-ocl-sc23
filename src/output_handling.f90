!****h* FHI-aims/output_handling
!  NAME
!   Routines for handling the output
!  SYNOPSIS

module output_handling

   use localorb_io,                  only:use_unit
   use mpi_tasks,                    only: aims_stop

   implicit none
   private

   public :: open_file, close_file
 
contains

! **************************************************************************************************
!> brief open file
!  o filename -- name of file
!  o front_name -- first part of the name
!  o middle_name -- middle part of the name
!  o end_name -- end part of the name
!  o unit_number -- unit number connected to file name 
! **************************************************************************************************
  subroutine open_file(filename, file_action, front_name, middle_name, end_name, extension,&
                       unit_number, open_failed)
 
     character(len=40), intent(out)                           :: filename
     character(len=*), intent(in)                             :: file_action
     character(len=*), intent(in)                             :: front_name
     character(len=*), intent(in), optional                   :: middle_name
     character(len=*), intent(in), optional                   :: end_name
     character(len=*), intent(in), optional                   :: extension

     integer, intent(out)                                     :: unit_number
     logical, intent(out), optional                           :: open_failed 
 
     character(*), parameter :: func = 'open_file'
 
     character(len=40)                                        :: my_middle_name, &
                                                                 my_end_name, &
                                                                 my_extension,&
                                                                 my_file_action
     logical                                                  :: is_open, exists
     integer                                                  :: istat

     my_middle_name = ""
     my_end_name = ""
     my_extension = 'dat'
     my_file_action = trim(adjustl(file_action))
     if(present(open_failed)) open_failed = .false.
     if(present(middle_name)) then
       if(front_name /= "") then
         my_middle_name = "_"//trim(adjustl(middle_name))
       else
         my_middle_name = middle_name
       endif 
     endif
                 
     if(present(end_name)) then
       if(middle_name /= "" .or. front_name /= "") then
         my_end_name = "_"//trim(adjustl(end_name))
       else
         my_end_name = end_name
       endif 
     endif

     if(present(extension)) then
       my_extension = extension
     endif

     filename =  trim(adjustl(front_name)) // trim(adjustl(my_middle_name)) // &
                 trim(adjustl(my_end_name)) // trim('.')// trim(adjustl(my_extension))
     inquire (FILE=trim(filename), opened=is_open, iostat=istat)
     if(is_open) call aims_stop('try to re-open the file '// func)
     if(istat /= 0) call aims_stop('An error occurred inquiring the file'//func) 
 
     unit_number = get_unit_number()
     if(my_file_action == 'read') then
       inquire(file=trim(adjustl(filename)), exist=exists)
       if(.not.exists) then 
         if(present(open_failed)) open_failed = .true.
         return
       endif
       open(unit=unit_number,file=filename,status='OLD',action=my_file_action,iostat=istat)
     elseif(my_file_action == 'write') then
       open(unit=unit_number,file=filename,status='REPLACE',action=my_file_action,iostat=istat)
     else
       open(unit=unit_number,file=filename,status='UNKNOWN',action=my_file_action,iostat=istat)
     endif 
     if(istat /= 0) call aims_stop('An error occurred opening the file '//func) 
     
  end subroutine open_file

! **************************************************************************************************
!> brief close file
!  o unit_number -- unit number to close 
! **************************************************************************************************
  subroutine close_file(unit_number)
 
     integer, intent(in)                                      :: unit_number
 
     character(len=*), parameter :: func = 'close_file'
 
     logical                                                  :: is_open, exists
     integer                                                  :: istat
            
     inquire (unit=unit_number, exist=exists, opened=is_open, iostat=istat)
     if(istat /= 0) call aims_stop('An error occurred inquiring the file '//func) 
 
     if (.not.exists) then
        call aims_stop('Unit number cannot be closed '//func)
     endif
     if (unit_number == use_unit) then
        call aims_stop('Attempt to close the default input unit number '//func)
     endif     
   
     if (is_open) then
       close(unit=unit_number,iostat=istat)
       if(istat /= 0) call aims_stop('An error occurred closing the file '//func)
     else
       call aims_stop('Attempt to close a file that is not opened '//func)
     endif
 
  end subroutine close_file

! **************************************************************************************************
!> brief find unit number that is not in use
! **************************************************************************************************
  function get_unit_number() result(unit_number)
  
    integer :: unit_number
    logical :: file_open
  
    unit_number = 9
    file_open = .TRUE.
    do while (file_open)
      unit_number = unit_number + 1
      inquire (unit_number, opened=file_open)
    enddo
  
  end function get_unit_number

end module output_handling
