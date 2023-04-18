!****h* FHI-aims/xml_write
! PURPOSE
!   Routines for writing output in XML format
! AUTHOR
!   Jan Hermann
! CREATION DATE
!   2015-05-13
! EXAMPLE
!   This shows how can the module be used.
!
!       use xml_write, only: &
!           xml_open_file, xml_close_file, &
!           xml_open, xml_close, xml_elem
!
!       real(8) :: energy, potential(:), kernel(:, :)
!       real(8) :: electron_repulsion_integrals(:, :, :, :)
!
!       call xml_open_file('my_file.xml', 'results')
!       call xml_open('exact_XC_functional')
!       call xml_elem('energy', energy)
!       call xml_elem('potential', potential)
!       call xml_elem('kernel', kernel_ts)
!       call xml_elem('eris', electron_repulsion_integrals)
!       call xml_close() ! closes the <exact_XC_functional> tag
!       call xml_close_file()
!
! COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e.V. Please note 
!   that any use of the "FHI-aims-Software" is subject to the terms and 
!   conditions of the respective license agreement.
!******
module xml_write

implicit none

private

public :: &
    xml_open_file, xml_close_file, xml_open, xml_close, xml_elem, &
    xml_attr, tostr

integer, parameter :: &
    root_name_length = 100, &
    tag_name_length = 100, &
    attrs_buff_str_length = 100, &
    elem_queue_length = 10, &
    buff_str_length = 200

integer, parameter :: &
    unit_min = 10, &
    unit_max = 100

!****t* xml_write/Xml_file_t
! PURPOSE
!   Keeps information about an XML file.
!******
type, public :: Xml_file_t
    logical :: muted = .false.
    integer :: io_unit = -1
    integer :: indent_level = 0
    integer :: tab_spaces = 4
    character(len=root_name_length) :: root_name
    character(len=tag_name_length) :: opened_tags(elem_queue_length)
end type Xml_file_t

!****t* xml_write/Xml_attr_t
! PURPOSE
!   Represents an XML attribute.
!******
type, public :: Xml_attr_t
    private

    character(len=attrs_buff_str_length) :: key
    character(len=attrs_buff_str_length) :: val
end type Xml_attr_t

!****v* xml_write/default_file
! PURPOSE
!   A default xml file that is used by the module functions is no file is 
!   explicitly specified.
!******
type(Xml_file_t) :: default_file

!****f* xml_write/tostr
! PURPOSE
!   Converts various types into string.
! INPUTS
!   * tostr (result) string -- The converted string.
!   * x (in) {real(8),integer} -- Input to be converted.
!******
interface tostr
    module procedure tostr_int_
    module procedure tostr_dble_
end interface

!****f* xml_write/xml_attr
! PURPOSE
!   Constructs XML attributes from different types.
! INPUTS
!   * xml_attr (result) Xml_attr_t -- The XML attribute.
!   * key (in) string -- The attribute keyword.
!   * val (in) {string,real(8),integer} -- The attribute value.
!******
interface xml_attr
    module procedure xml_attr_str_
    module procedure xml_attr_dble_
    module procedure xml_attr_int_
    module procedure xml_attr_vector_str_
    module procedure xml_attr_vector_int_
end interface

!****f* xml_write/xml_elem
! PURPOSE
!   A generic function that can write various types of values
!   into the XML file.
! 
!   The format of the array is (example)
! 
!       <tag type="dble" size="3">
!           1.0 2.0 3.0
!       </tag>
!
!   If the name or unit argument are given, the format will be
!  
!       <tag>
!           <name>NAME</name>
!           <unit>UNIT</unit>
!           <value type=...>
!               ...
!           </value>
!       </tag>
!
! INPUTS
!   * tag (in) string -- Name of the element being written.
!   * val (in) ? -- Can be a string, 0-to-7-rank double or 0-to-1-rank 
!     integer.
!   * name (in) optional string -- Optional name element for the object.
!   * unit (in) optional string -- Optional physical unit specification.
!   * file (inout) optional Xml_file_t -- Object representing a XML file. If not 
!     given, a default private module object will be used.
!******
interface xml_elem
    module procedure xml_elem_str_
    module procedure xml_scalar_dble_
    module procedure xml_scalar_int_
    module procedure xml_vector_dble_
    module procedure xml_vector_int_
    module procedure xml_matrix_dble_
    module procedure xml_3d_dble_
    module procedure xml_4d_dble_
    module procedure xml_5d_dble_
    module procedure xml_6d_dble_
    module procedure xml_7d_dble_
end interface

contains


!****f* xml_write/xml_open_file
! PURPOSE
!   Initializes an object representing a XML file, opens the associated file, 
!   writes the XML version.
! INPUTS
!   * filename (in) string -- Name of the XML file that will be used.
!   * root (in) optional string -- Name of the XML root element. If not 
!     given, "root" will be used.
!   * file (out) optional Xml_file_t -- Object representing a XML file. If not 
!     given, a default private module object will be used.
!******
subroutine xml_open_file(filename, root, file, tab_spaces)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: root
    type(Xml_file_t), intent(out), optional :: file
    integer, intent(in), optional :: tab_spaces

    type(Xml_file_t) :: xml_file

    if (present(root)) then
        xml_file%root_name = root
    else
        xml_file%root_name = "root"
    endif
    open(unit=newunit(xml_file%io_unit), file=filename)
    if (present(tab_spaces)) then
        xml_file%tab_spaces = tab_spaces
    endif
    write (xml_file%io_unit, "(a)") '<?xml version="1.0"?>'
    call xml_open(xml_file%root_name, file=xml_file)
    if (present(file)) then
        file = xml_file
    else
        default_file = xml_file
    endif
end subroutine xml_open_file


!****f* xml_write/xml_close_file
! PURPOSE
!   Finalizes and closes a XML file.
! INPUTS
!   * file (inout) optional Xml_file_t -- Object representing the XML file. If 
!     not given, a default private module object will be used.
!******
subroutine xml_close_file(file)
    type(Xml_file_t), intent(inout), optional :: file

    type(Xml_file_t) :: xml_file

    if (present(file)) then
        xml_file = file
    else
        xml_file = default_file
    endif
    if (xml_file%io_unit < unit_min) return
    call xml_close(file=xml_file)
    if (xml_file%io_unit >= 0) then
        if (xml_file%indent_level > 0) then
            write (xml_file%io_unit, "(a)") &
                "**error: there are unclosed tags left"
        endif
        close(xml_file%io_unit)
        xml_file%io_unit = -1
    endif
    if (present(file)) then
        file = xml_file
    else
        default_file = xml_file
    endif
end subroutine xml_close_file


!****f* xml_write/xml_close
! PURPOSE
!   Opens a new element and leaves it open.
! INPUTS
!   * tag (in) string -- Name of the tag being open.
!   * attrs(:) (in) optional Xml_attr_t -- XML attributes of the element.
!   * file (inout) optional Xml_file_t -- Object representing the XML file.  If 
!     not given, a default private module object will be used.
!******
subroutine xml_open(tag, attrs, file)
    character(len=*), intent(in) :: tag
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    type(Xml_file_t), intent(inout), optional :: file

    type(Xml_file_t) :: xml_file
    integer :: io_unit
    integer :: i_attr

    if (present(file)) then
        xml_file = file
    else
        xml_file = default_file
    endif
    if (xml_file%io_unit < unit_min .or. xml_file%muted) return
    io_unit = xml_file%io_unit
    call xml_indent(xml_file)
    write (io_unit, "(a,a)", advance="no") "<", trim(tag)
    if (present(attrs)) then
        do i_attr = 1, size(attrs)
            write (io_unit, "(1x,a,a,a,a)", advance="no") &
                trim(attrs(i_attr)%key), '="', trim(attrs(i_attr)%val), '"'
        enddo
    endif
    write (io_unit, "(a)") ">"
    xml_file%indent_level = xml_file%indent_level + 1
    xml_file%opened_tags(xml_file%indent_level) = tag
    if (present(file)) then
        file = xml_file
    else
        default_file = xml_file
    endif
end subroutine xml_open


!****f* xml_write/xml_close
! PURPOSE
!   Closes last element which was opened and left open.
! INPUTS
!   * file (inout) optional Xml_file_t -- Object representing the XML file.  If 
!     not given, a default private module object will be used.
!******
subroutine xml_close(file)
    type(Xml_file_t), intent(inout), optional :: file

    type(Xml_file_t) :: xml_file
    integer :: io_unit

    if (present(file)) then
        xml_file = file
    else
        xml_file = default_file
    endif
    if (xml_file%io_unit < unit_min .or. xml_file%muted) return
    io_unit = xml_file%io_unit
    if (xml_file%indent_level < 1) then
        write (io_unit, "(a)") "**Error: there is not tag left to close"
        close (xml_file%io_unit)
        xml_file%io_unit = -1
    else
        xml_file%indent_level = xml_file%indent_level - 1
        call xml_indent(xml_file)
        write (io_unit, "(a,a,a)") &
            "</", trim(xml_file%opened_tags(xml_file%indent_level+1)), ">"
    endif
    if (present(file)) then
        file = xml_file
    else
        default_file = xml_file
    endif
end subroutine xml_close


function xml_attr_str_(key, val) result(attr)
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: val
    type(Xml_attr_t) :: attr

    attr = Xml_attr_t(key, val)
end function xml_attr_str_


function xml_attr_dble_(key, val) result(attr)
    character(len=*), intent(in) :: key
    real(8), intent(in) :: val
    type(Xml_attr_t) :: attr

    attr = Xml_attr_t(key, tostr(val))
end function xml_attr_dble_


function xml_attr_int_(key, val) result(attr)
    character(len=*), intent(in) :: key
    integer, intent(in) :: val
    type(Xml_attr_t) :: attr

    attr = Xml_attr_t(key, tostr(val))
end function xml_attr_int_


function xml_attr_vector_str_(key, val) result(attr)
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: val(:)
    type(Xml_attr_t) :: attr

    character(len=buff_str_length) :: format, converted
    integer :: i

    write(format,'(a,i0,a)') '(', size(val), '(a,", "))'
    write(converted, trim(format)) (trim(val(i)), i=1,size(val))
    attr = Xml_attr_t(key, "("//trim(converted)//")")
end function xml_attr_vector_str_


function xml_attr_vector_int_(key, val) result(attr)
    character(len=*), intent(in) :: key
    integer, intent(in) :: val(:)
    type(Xml_attr_t) :: attr

    attr = xml_attr(key, tostr(val))
end function xml_attr_vector_int_


subroutine xml_elem_str_(tag, val, attrs, file)
    character(len=*), intent(in) :: tag
    character(len=*), intent(in) :: val
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    type(Xml_file_t), intent(inout), optional :: file

    type(Xml_file_t) :: xml_file
    integer :: io_unit
    integer :: i_attr

    if (present(file)) then
        xml_file = file
    else
        xml_file = default_file
    endif
    if (xml_file%io_unit < unit_min .or. xml_file%muted) return
    io_unit = xml_file%io_unit
    call xml_indent(xml_file)
    write (io_unit, "(a,a)", advance="no") "<", trim(tag)
    if (present(attrs)) then
        do i_attr = 1, size(attrs)
            write (io_unit, "(1x,a,a,a,a)", advance="no") &
                trim(attrs(i_attr)%key), '="', trim(attrs(i_attr)%val), '"'
        enddo
    endif
    if (val == "") then
        write (io_unit, "(1x,a)") "/>"
    else
        write (io_unit, "(a,a,a,a,a)") ">", trim(val), "</", trim(tag), ">"
    endif
    if (present(file)) then
        file = xml_file
    else
        default_file = xml_file
    endif
end subroutine xml_elem_str_


!****f* xml_write/xml_vector
! PURPOSE
!   Routine for writing 1-rank vectors. See `xml_elem` for details.
! INPUTS
!   * tag (in) string
!   * arr_dble(:) (in) optional real(8) -- Double vector to be written.
!   * arr_int(:) (in) optional integer -- Integer vector to be written.
!   * attrs(:) (in) optional Xml_attr_t
!   * name (in) optional string
!   * unit (in) optional string
!   * file (inout) optional Xml_file_t
!******
subroutine xml_vector(tag, arr_dble, arr_int, attrs, name, unit, file)
    character(len=*), intent(in) :: tag
    real(8), intent(in), optional :: arr_dble(:)
    integer, intent(in), optional :: arr_int(:)
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: unit
    type(Xml_file_t), intent(inout), optional :: file

    type(Xml_file_t) :: xml_file
    integer :: n_block, i_arr, i, n_arr
    integer :: io_unit
    character(len=buff_str_length) :: inner_tag
    type(Xml_attr_t) :: type_attr

    if (present(file)) then
        xml_file = file
    else
        xml_file = default_file
    endif
    if (xml_file%io_unit < unit_min .or. xml_file%muted) return
    io_unit = xml_file%io_unit
    if (present(name) .or. present(unit)) then
        inner_tag = 'value'
        call xml_open(tag, attrs, file=xml_file)
        if (present(name)) then
            call xml_elem('name', name, file=xml_file)
        end if
        if (present(unit)) then
            call xml_elem('unit', unit, file=xml_file)
        end if
    else
        inner_tag = tag
    end if
    type_attr%key = 'type'
    if (present(arr_dble)) then
        type_attr%val = 'dble'
        n_block = 3
        n_arr = size(arr_dble)
    elseif (present(arr_int)) then
        type_attr%val = 'int'
        n_block = 10
        n_arr = size(arr_int)
    end if
    if (present(attrs)) then
        call xml_open(inner_tag, attrs, file=xml_file)
    else
        call xml_open( &
            inner_tag, &
            (/ type_attr, xml_attr('size', n_arr) /), &
            file=xml_file)
    end if
    i_arr = 0
    do while (i_arr < n_arr)
        call xml_indent(xml_file)
        do i = 1, min(n_block, n_arr-i_arr)
            if (present(arr_dble)) then
                write (io_unit, '(a," ")', advance='no') &
                    trim(tostr(arr_dble(i_arr+i)))
            elseif (present(arr_int)) then
                write (io_unit, '(a," ")', advance='no') &
                    trim(tostr(arr_int(i_arr+i)))
            end if
        end do
        write (io_unit, '()')
        i_arr = i_arr+n_block
    end do
    call xml_close(xml_file)
    if (present(name) .or. present(unit)) then
        call xml_close(xml_file)
    end if
    if (present(file)) then
        file = xml_file
    else
        default_file = xml_file
    endif
end subroutine xml_vector


subroutine xml_vector_dble_(tag, arr, attrs, name, unit, file)
    character(len=*), intent(in) :: tag
    real(8), intent(in) :: arr(:)
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: unit
    type(Xml_file_t), intent(inout), optional :: file

    call xml_vector(tag, arr_dble=arr, attrs=attrs, name=name, unit=unit, file=file)
end subroutine


subroutine xml_vector_int_(tag, arr, attrs, name, unit, file)
    character(len=*), intent(in) :: tag
    integer, intent(in) :: arr(:)
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: unit
    type(Xml_file_t), intent(inout), optional :: file

    call xml_vector(tag, arr_int=arr, attrs=attrs, name=name, unit=unit, file=file)
end subroutine


subroutine xml_scalar_dble_(tag, val, attrs, name, unit, file)
    character(len=*), intent(in) :: tag
    real(8), intent(in) :: val
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: unit
    type(Xml_file_t), intent(inout), optional :: file

    call xml_elem(tag, (/ val /), attrs, name, unit, file)
end subroutine


subroutine xml_scalar_int_(tag, val, attrs, name, unit, file)
    character(len=*), intent(in) :: tag
    integer, intent(in) :: val
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: unit
    type(Xml_file_t), intent(inout), optional :: file

    call xml_elem(tag, (/ val /), attrs, name, unit, file)
end subroutine


!****f* xml_write/xml_array_dble
! PURPOSE
!   Routine for writing N-rank arrays as vectors of vectors of..., see 
!   `xml_elem` and `xml_vector` for details.
! INPUTS
!   * tag (in) string
!   * arr(*) (in) real(8) -- A general N-rank double array.
!   * dims(:) (in) integer -- Array specifying the dimension of the array. Can 
!     be obtained with `shape`.
!   * attrs(:) (in) optional Xml_attr_t
!   * name (in) optional string
!   * unit (in) optional string
!   * file (inout) optional Xml_file_t
!******
subroutine xml_array_dble(tag, arr, dims, attrs, name, unit, file)
    character(len=*), intent(in) :: tag
    real(8), intent(in) :: arr(*)
    integer, intent(in) :: dims(:)
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: unit
    type(Xml_file_t), intent(inout), optional :: file

    type(Xml_file_t) :: xml_file
    integer :: i_dim, n_dim, i
    integer :: running_index
    integer :: multi_index(size(dims))
    integer :: n_x, i_x
    character(len=attrs_buff_str_length) :: str_buff
    character(len=tag_name_length) :: inner_tag
    type(Xml_attr_t) :: type_attr, size_attr, idx_attr

    if (present(file)) then
        xml_file = file
    else
        xml_file = default_file
    endif
    if (present(name) .or. present(unit)) then
        inner_tag = 'value'
        call xml_open(tag, attrs, file=xml_file)
        if (present(name)) then
            call xml_elem('name', name, file=xml_file)
        end if
        if (present(unit)) then
            call xml_elem('unit', unit, file=xml_file)
        end if
    else
        inner_tag = tag
    end if
    type_attr = Xml_attr_t('type', 'dble')
    size_attr%key = 'size'
    str_buff(1:100) = tostr(dims(1))
    n_dim = size(dims)
    do i_dim = 2, n_dim
        str_buff = trim(str_buff) // " " // tostr(dims(i_dim))
    end do
    size_attr%val = str_buff
    if (.not. present(attrs) .or. present(name) .or. present(unit)) then
        call xml_open( &
            inner_tag, &
            (/ type_attr, size_attr /), &
            file=xml_file)
    else
        call xml_open( &
            inner_tag, &
            (/ (attrs(i), i = 1, size(attrs)), &
            type_attr, size_attr /), &
            file=xml_file)
    end if
    idx_attr%key = 'index'
    running_index = 1
    multi_index(:) = 1
    if (n_dim == 1) then
        n_x = 1
    else
        n_x = n_dim-1
    end if
    do
        do
            str_buff = 'x'
            do i_x = 2, n_x
                str_buff = trim(str_buff) // ' x'
            end do
            do i_dim = n_x+1, n_dim
                str_buff = trim(str_buff) // ' ' // tostr(multi_index(i_dim))
            end do
            idx_attr%val = str_buff
            if (n_x == 1) then
                call xml_elem( &
                    "vector", &
                    arr(running_index:running_index+dims(1)-1), &
                    (/ idx_attr /), &
                    file=xml_file)
                exit
            else
                call xml_open("vector", (/ idx_attr /), file=xml_file)
                n_x = n_x-1
            end if
        end do
        if (n_dim == 1) then
            exit
        end if
        running_index = running_index+dims(1)
        multi_index(2) = multi_index(2)+1
        do i_dim = 2, n_dim-1
            if (multi_index(i_dim) > dims(i_dim)) then
                multi_index(i_dim) = 1
                multi_index(i_dim+1) = multi_index(i_dim+1)+1
                call xml_close(xml_file)
                n_x = n_x+1
            else
                exit
            endif
        enddo
        if (multi_index(n_dim) > dims(n_dim)) then
            exit
        end if
    enddo
    call xml_close(xml_file) ! close inner tag
    if (present(name) .or. present(unit)) then
        call xml_close(xml_file) ! close outer tag
    end if
    if (present(file)) then
        file = xml_file
    else
        default_file = xml_file
    endif
end subroutine xml_array_dble


subroutine xml_matrix_dble_(tag, matrix, attrs, name, unit, file)
    character(len=*), intent(in) :: tag
    real(8), intent(in) :: matrix(:, :)
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: unit
    type(Xml_file_t), intent(inout), optional :: file

    call xml_array_dble(tag, matrix, shape(matrix), attrs, name, unit, file)
end subroutine xml_matrix_dble_


subroutine xml_3d_dble_(tag, arr, attrs, name, unit, file)
    character(len=*), intent(in) :: tag
    real(8), intent(in) :: arr(:, :, :)
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: unit
    type(Xml_file_t), intent(inout), optional :: file

    call xml_array_dble(tag, arr, shape(arr), attrs, name, unit, file)
end subroutine


subroutine xml_4d_dble_(tag, arr, attrs, name, unit, file)
    character(len=*), intent(in) :: tag
    real(8), intent(in) :: arr(:, :, :, :)
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: unit
    type(Xml_file_t), intent(inout), optional :: file

    call xml_array_dble(tag, arr, shape(arr), attrs, name, unit, file)
end subroutine


subroutine xml_5d_dble_(tag, arr, attrs, name, unit, file)
    character(len=*), intent(in) :: tag
    real(8), intent(in) :: arr(:, :, :, :, :)
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: unit
    type(Xml_file_t), intent(inout), optional :: file

    call xml_array_dble(tag, arr, shape(arr), attrs, name, unit, file)
end subroutine


subroutine xml_6d_dble_(tag, arr, attrs, name, unit, file)
    character(len=*), intent(in) :: tag
    real(8), intent(in) :: arr(:, :, :, :, :, :)
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: unit
    type(Xml_file_t), intent(inout), optional :: file

    call xml_array_dble(tag, arr, shape(arr), attrs, name, unit, file)
end subroutine


subroutine xml_7d_dble_(tag, arr, attrs, name, unit, file)
    character(len=*), intent(in) :: tag
    real(8), intent(in) :: arr(:, :, :, :, :, :, :)
    type(Xml_attr_t), intent(in), optional :: attrs(:)
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: unit
    type(Xml_file_t), intent(inout), optional :: file

    call xml_array_dble(tag, arr, shape(arr), attrs, name, unit, file)
end subroutine


character(len=buff_str_length) elemental function tostr_int_(k, format)
    integer, intent(in) :: k
    character(*), intent(in), optional :: format

    if (present(format)) then
        write (tostr_int_, format) k
    else
        write (tostr_int_, "(i20)") k
    end if
    tostr_int_ = adjustl(tostr_int_)
end function tostr_int_


character(len=buff_str_length) elemental function tostr_dble_(x, format)
    double precision, intent(in) :: x
    character(*), intent(in), optional :: format

    if (present(format)) then
        write (tostr_dble_, format) x
    else
        write (tostr_dble_, "(g50.17e3)") x
    end if
    tostr_dble_ = adjustl(tostr_dble_)
end function tostr_dble_


!****f* xml_write/xml_indent
! PURPOSE
!   Writes the current indentation to the file.
! INPUTS
!   * file (in) Xml_file_t -- Object representing a XML file.
!   * add (in) optional integer -- If present, indent by `add` levels more.
!******
subroutine xml_indent(xml_file, add)
    type(Xml_file_t), intent(in) :: xml_file
    integer, intent(in), optional :: add

    character(len=xml_file%tab_spaces) :: tab
    integer :: i_tab
    integer :: indent_level

    tab = " "
    indent_level = xml_file%indent_level
    if (present(add)) then
        indent_level = indent_level+add
    end if
    do i_tab = 1, indent_level
        write (xml_file%io_unit, "(a)", advance="no") tab
    end do
end subroutine xml_indent


integer function newunit(unit) result(i_io)
    integer, intent(out), optional :: unit

    logical :: unit_open

    do i_io = unit_min, unit_max
        inquire (unit=i_io, opened=unit_open)
        if (.not. unit_open) then
            if (present(unit)) then
                unit = i_io
                return
            end if
        endif
    enddo
    i_io = -1
end function newunit


end module xml_write
