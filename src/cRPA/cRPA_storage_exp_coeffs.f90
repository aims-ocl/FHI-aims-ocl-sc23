!Storage of RI-LVL expansion coeffs
MODULE cRPA_storage_exp_coeffs
    USE basis
    USE pbc_lists
    USE prodbas
    USE dimensions
!only debug
  use hartree_fock
  use hartree_fock_p0
!---


    USE cRPA_view
    USE cRPA_parallelism_storage
    USE cRPA_flow_realspace_vector
    USE cRPA_calculation_integration
    USE cRPA_flow_adaptive_grid
    USE cRPA_storage_norms
    IMPLICIT NONE

    type expansion_coefficient
        REAL*8, ALLOCATABLE, DIMENSION(:) :: c_compressed


        REAL*8, ALLOCATABLE, DIMENSION(:) :: column_norm_max_row, &
                                             column_norm_max_column

        INTEGER :: n_basbas
        INTEGER :: last_stored_cell, &
                   dictionary_length, &
                   dictionary_pos, &
                   c_compressed_size, &
                   all_cells_first_filled_row_onsite, &
                   all_cells_last_filled_row_onsite, &
                   all_cells_extend_filled_row_onsite, &
                   all_cells_first_filled_column_onsite, &
                   all_cells_last_filled_column_onsite, &
                   all_cells_extend_filled_column_onsite, &
                   all_cells_first_filled_row_offsite, &
                   all_cells_last_filled_row_offsite, &
                   all_cells_extend_filled_row_offsite, &
                   all_cells_first_filled_column_offsite, &
                   all_cells_last_filled_column_offsite, &
                   all_cells_extend_filled_column_offsite

        INTEGER, ALLOCATABLE, DIMENSION(:) :: cell_start, &
                                              cell_end, &
                                              dictionary_start, &
                                              dictionary_end, &
                                              dictionary
        LOGICAL :: is_stored_locally

    end type expansion_coefficient

    type expansion_coefficients
        type(expansion_coefficient), ALLOCATABLE, DIMENSION(:) :: coeff

        REAL*8, ALLOCATABLE, DIMENSION(:) :: loaded_column_norm_max_row, &
                                             loaded_column_norm_max_column, &
                                             requested_coeff_compressed
        INTEGER, ALLOCATABLE, DIMENSION(:) :: basbas_to_id
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: basbas_local_by_id


        type(realspace_vectors) :: rvecs
        type(basisfunc_distribution) :: basis_dist, &
                                        basis_dist_requested

        INTEGER :: n_rvecs, n_basbas,n_cells, &
                   n_basbas_local_max,n_basbas_local, &
                   n_max_all_cells_filled_rows

        INTEGER :: requested_id

        LOGICAL :: is_exp_coeff_storage_inited

    end type expansion_coefficients

    type(expansion_coefficients) :: exp_coeffs_local, &
                                    exp_coeffs_local_offsite

    INTEGER :: EXP_COEFFICIENT_UNASSIGNED
    PARAMETER (EXP_COEFFICIENT_UNASSIGNED=-1)

    INTEGER :: EXP_COEFFICIENT_NO_REQUESTED
    PARAMETER (EXP_COEFFICIENT_NO_REQUESTED=-2)

    INTEGER :: EXP_COEFFICIENT_LOCAL
    PARAMETER (EXP_COEFFICIENT_LOCAL=1)

    INTEGER :: EXP_COEFFICIENT_REQUESTED
    PARAMETER (EXP_COEFFICIENT_REQUESTED=2)


    REAL*8 :: SPARSITY_LIMIT
    PARAMETER (SPARSITY_LIMIT = 1.d-12)

CONTAINS

!expansion coefficient functions (storage by basbas)
    SUBROUTINE init_exp_coeff_storage(this,n_basbas,n_cells,my_basis_dist)
        type(expansion_coefficients) :: this
        INTEGER, INTENT(IN) :: n_basbas,n_cells
        type(basisfunc_distribution),INTENT(IN) :: my_basis_dist

        CHARACTER(*), PARAMETER :: funcname = 'init_exp_coeff_storage'
        INTEGER :: n_vecs, i_basbas, allocation_info

        if(this%is_exp_coeff_storage_inited) then
            write(use_unit,*) "Storage already inited!"
            stop
        end if

!        CALL write_stdout("Initialize expansion coefficient storage for "&
!                          //num2str(n_basbas)//"exp coeffs")

        n_vecs = my_basis_dist%rvecs%n_vecs
        ALLOCATE(this%coeff(n_basbas), stat=allocation_info)
        CALL check_allocation(allocation_info, 'coeff', funcname)
        this%coeff(1:n_basbas)%n_basbas = EXP_COEFFICIENT_UNASSIGNED


        ALLOCATE(this%loaded_column_norm_max_row(n_vecs), stat=allocation_info)
        CALL check_allocation(allocation_info, 'loaded_column_norm_max_row', funcname)
        ALLOCATE(this%loaded_column_norm_max_column(n_vecs), stat=allocation_info)
        CALL check_allocation(allocation_info, 'loaded_column_norm_max_column', funcname)


        this%basis_dist = my_basis_dist

        this%n_basbas_local_max = EXP_COEFFICIENT_UNASSIGNED
        this%n_max_all_cells_filled_rows = EXP_COEFFICIENT_UNASSIGNED
        this%n_basbas = n_basbas
        this%n_cells = n_cells
        this%n_rvecs = n_vecs
        this%requested_id = EXP_COEFFICIENT_NO_REQUESTED
        CALL clone_realspace_vectors(my_basis_dist%rvecs, this%rvecs)

        this%is_exp_coeff_storage_inited = .TRUE.

        ALLOCATE(this%basbas_to_id(n_basbas), stat=allocation_info)
        CALL check_allocation(allocation_info, 'basbas_to_id', funcname)

        do i_basbas = 1,n_basbas
             this%coeff(i_basbas)%last_stored_cell = 0
             this%coeff(i_basbas)%dictionary_length = 1024
             this%coeff(i_basbas)%dictionary_pos = 0
             this%coeff(i_basbas)%c_compressed_size = 0
             this%coeff(i_basbas)%all_cells_first_filled_row_onsite = 0
             this%coeff(i_basbas)%all_cells_last_filled_row_onsite = 0
             this%coeff(i_basbas)%all_cells_extend_filled_row_onsite = 0
             this%coeff(i_basbas)%all_cells_first_filled_column_onsite = 0
             this%coeff(i_basbas)%all_cells_last_filled_column_onsite = 0
             this%coeff(i_basbas)%all_cells_extend_filled_column_onsite = 0
             this%coeff(i_basbas)%all_cells_first_filled_row_offsite = 0
             this%coeff(i_basbas)%all_cells_last_filled_row_offsite = 0
             this%coeff(i_basbas)%all_cells_extend_filled_row_offsite = 0
             this%coeff(i_basbas)%all_cells_first_filled_column_offsite = 0
             this%coeff(i_basbas)%all_cells_last_filled_column_offsite = 0
             this%coeff(i_basbas)%all_cells_extend_filled_column_offsite = 0


             ALLOCATE(this%coeff(i_basbas)%dictionary(this%coeff(i_basbas)%dictionary_length), stat=allocation_info)
             CALL check_allocation(allocation_info, 'dictionary', funcname)

             ALLOCATE(this%coeff(i_basbas)%cell_start(n_cells+1), stat=allocation_info)
             CALL check_allocation(allocation_info, 'cell_start', funcname)
             ALLOCATE(this%coeff(i_basbas)%cell_end(n_cells+1), stat=allocation_info)
             CALL check_allocation(allocation_info, 'cell_end', funcname)
             ALLOCATE(this%coeff(i_basbas)%dictionary_start(n_cells+1), stat=allocation_info)
             CALL check_allocation(allocation_info, 'dictionary_start', funcname)
             ALLOCATE(this%coeff(i_basbas)%dictionary_end(n_cells+1), stat=allocation_info)
             CALL check_allocation(allocation_info, 'dictionary_end', funcname)

             this%coeff(i_basbas)%is_stored_locally = .FALSE.
        enddo
    END SUBROUTINE init_exp_coeff_storage

    SUBROUTINE set_exp_coeff_requested_basisfunc_distribution(this,dist_requested)
        type(expansion_coefficients) :: this
        type(basisfunc_distribution), INTENT(IN) :: dist_requested

        this%basis_dist_requested = dist_requested

    END SUBROUTINE set_exp_coeff_requested_basisfunc_distribution

    LOGICAL FUNCTION is_exp_coeffs_complete(this)
        type(expansion_coefficients) :: this

        INTEGER :: i_basbas

        do i_basbas = 1, n_basbas
            if(.NOT. is_exp_coeff_stored(this, i_basbas)) then
                is_exp_coeffs_complete = .FALSE.
                return
            endif
        enddo

        is_exp_coeffs_complete = .TRUE.
    END FUNCTION is_exp_coeffs_complete

    LOGICAL FUNCTION is_exp_coeff_stored(this,n_basbas)
        type(expansion_coefficients) :: this
        INTEGER, INTENT(IN) :: n_basbas

        is_exp_coeff_stored = ANY(this%coeff(:)%n_basbas == n_basbas)
    END FUNCTION is_exp_coeff_stored

    INTEGER FUNCTION find_exp_coeff_by_basbas(this,n_basbas)
        type(expansion_coefficients) :: this
        INTEGER, INTENT(IN) :: n_basbas
        INTEGER :: i

        if(.NOT. is_exp_coeff_stored(this,n_basbas)) then
            write(use_unit,*) "Atemp to find basbas", n_basbas, " that is not stored"
            stop
        end if

        do i = 1, ubound(this%coeff,1)
            if(this%coeff(i)%n_basbas == n_basbas) then
                find_exp_coeff_by_basbas = i
                return
            end if
        end do
    END FUNCTION find_exp_coeff_by_basbas

    SUBROUTINE load_exp_column(this,i_basbas, i_vec_inversion,&
                               i_local_or_requested, &
                               n_vecs,column_onsite,column_offsite)
        type(expansion_coefficients) :: this
        INTEGER,INTENT(IN) :: i_basbas
        INTEGER, DIMENSION(this%rvecs%n_vecs),INTENT(IN) :: i_vec_inversion
        INTEGER, INTENT(IN) :: i_local_or_requested
        INTEGER, INTENT(OUT) :: n_vecs
        DOUBLE COMPLEX,DIMENSION(n_basis, &
                                 n_basis, &
                                 this%n_cells), INTENT(OUT) :: column_onsite, &
                                                               column_offsite

        DOUBLE PRECISION :: tmp_real(n_basis,n_basis)
        INTEGER :: i_cell, i_vector,i_vector_inversion
        type(realspace_vector) :: vector_inversion

!CALL write_debug("basbas"//num2str(i_basbas))

        n_vecs = this%rvecs%n_vecs

        do i_cell = 1, this%n_cells+1
            if(i_cell == 1) then
                CALL decompress_cell(this,i_basbas,1, i_local_or_requested,&
                                       tmp_real(:,:))
                column_onsite(:,:,1) = tmp_real(:,:) * DCMPLX(1.d0,0.d0)
            elseif(i_cell > 1 .AND. i_cell <= this%n_cells) then
                   CALL decompress_cell(this,i_basbas,i_cell, i_local_or_requested, &
                                        tmp_real)
                   column_onsite(:,:,i_cell) = tmp_real(:,:) * DCMPLX(1.d0,0.d0)

                   column_offsite(:,:,(i_cell)) = transpose(column_onsite(:,:,i_cell))
            elseif(i_cell == n_cells +1) then
                   CALL decompress_cell(this,i_basbas,n_cells+1, i_local_or_requested, &
                                        tmp_real)
                   column_onsite(:,:,1) = column_onsite(:,:,1) + tmp_real(:,:) * DCMPLX(1.d0,0.d0)
                   column_offsite(:,:,1) = transpose(column_onsite(:,:,1))
            endif
        enddo
column_offsite(:,:,:) = column_onsite(:,:,:)
 !       do i_cell = 1, this%n_cells+1
 !           do i_vector=1,n_basis
 !               write(use_unit,*) i_cell,i_vector,"---",column_offsite(:,i_cell,i_vector)
 !           enddo
 !       enddo

!        do i_vector = 1, n_vecs
!            i_vector_inversion = i_vec_inversion(i_vector)

!           this%loaded_column(:,:,i_vector) = &
!                                                  this%coeff(pos)%c(:,:,i_vector)

!            this%loaded_column_norm_max_row(:) = this%coeff(i_basbas)%column_norm_max_row(:)
!            this%loaded_column_norm_max_column(:) = this%coeff(i_basbas)%column_norm_max_column(:)
!        enddo

    END SUBROUTINE load_exp_column

    SUBROUTINE store_cell_compressed(this, i_basbas, i_cell, data)
        type(expansion_coefficients) :: this
        INTEGER, INTENT(IN) :: i_basbas, &
                               i_cell
        REAL*8, DIMENSION(n_basis,n_basis),INTENT(IN) :: data
        INTEGER :: i,j, i_dropped, i_accepted, i_total_add, i_dictionary_add, &
                   i_first_filled_row, i_last_filled_row, &
                   i_first_filled_column, i_last_filled_column

        if(this%coeff(i_basbas)%last_stored_cell+1/=i_cell) then
            stop "You have to store the exp coeff cells incremental!"
        endif
!write(use_unit,*) "Storing basbas,cell", i_basbas,i_cell
        i_accepted = 0
        i_dropped = 0
        i_total_add = 0
        i_dictionary_add = 0
        i_first_filled_row = 0
        i_last_filled_row = 0
        i_first_filled_column = 0
        i_last_filled_column = 0

        do j=1,n_basis
            do i=1,n_basis
                if(abs(data(i,j)) > SPARSITY_LIMIT) then
                   if(i_dropped > 0) then
                      CALL add_dropped()
                   endif

                   i_accepted = i_accepted + 1
                   i_total_add = i_total_add + 1

                   CALL set_first_filled_row_or_column(i, i_first_filled_row)
                   CALL set_last_filled_row_or_column(i, i_last_filled_row)
                   CALL set_first_filled_row_or_column(j, i_first_filled_column)
                   CALL set_last_filled_row_or_column(j, i_last_filled_column)
                else
                   if(i_accepted > 0 .OR. &
                      (i_accepted == 0 .AND. i==1 .AND. j == 1)) then
                       CALL add_accepted()
                   endif

                   i_dropped = i_dropped + 1
                endif
            enddo
        enddo

        if(i_accepted > 0) CALL add_accepted()
        CALL add_dropped()

        if(i_cell > 1) then
            this%coeff(i_basbas)%cell_start(i_cell) = this%coeff(i_basbas)%cell_end(i_cell-1) + 1
            this%coeff(i_basbas)%cell_end(i_cell) = this%coeff(i_basbas)%cell_start(i_cell) + &
                                                      i_total_add !- 1


           this%coeff(i_basbas)%dictionary_start(i_cell) = this%coeff(i_basbas)%dictionary_end(i_cell-1) + 1
           this%coeff(i_basbas)%dictionary_end(i_cell) = this%coeff(i_basbas)%dictionary_start(i_cell) + &
                                                          i_dictionary_add - 1
        else
            this%coeff(i_basbas)%cell_start(1) = 1
            this%coeff(i_basbas)%cell_end(1) = i_total_add

            this%coeff(i_basbas)%dictionary_start(1) = 1
            this%coeff(i_basbas)%dictionary_end(1) = i_dictionary_add
        endif

        CALL check_c_compress_size(this%coeff(i_basbas)%cell_end(i_cell))
        CALL store_cell_compressed_by_dictionary(this,i_basbas, i_cell, data)

!            write(use_unit,*) i_cell,i_first_filled_row, i_last_filled_row

        CALL set_all_cells_x_row()

        this%coeff(i_basbas)%last_stored_cell = this%coeff(i_basbas)%last_stored_cell+1
        this%coeff(i_basbas)%is_stored_locally = .TRUE.
    CONTAINS
        SUBROUTINE add_accepted()
           CALL add_to_dictionary(i_accepted)
           i_dictionary_add = i_dictionary_add + 1
           i_accepted = 0
        END SUBROUTINE

        SUBROUTINE add_dropped()
           CALL add_to_dictionary(i_dropped)
           i_dictionary_add = i_dictionary_add + 1
            i_dropped = 0
        END SUBROUTINE add_dropped

        SUBROUTINE set_first_filled_row_or_column(i_row_or_column, i_first_filled_row_or_column)
            INTEGER, INTENT(IN) :: i_row_or_column
            INTEGER, INTENT(INOUT) :: i_first_filled_row_or_column

            if(i_first_filled_row_or_column==0 &
               .OR. &
               i_first_filled_row_or_column > i_row_or_column) then
                  i_first_filled_row_or_column = i_row_or_column
            endif
        END SUBROUTINE set_first_filled_row_or_column

        SUBROUTINE set_last_filled_row_or_column(i_row_or_column,i_last_filled_row_or_column)
            INTEGER, INTENT(IN) :: i_row_or_column
            INTEGER, INTENT(INOUT) :: i_last_filled_row_or_column


            if(i_last_filled_row_or_column < i_row_or_column) then
                  i_last_filled_row_or_column = i_row_or_column
            endif
        END SUBROUTINE set_last_filled_row_or_column

        SUBROUTINE set_first_last_extend_all_cells(i_first, i_last, &
                                                   i_all_first, i_all_last, i_all_extend)
            INTEGER, INTENT(IN) :: i_first, i_last
            INTEGER, INTENT(INOUT) :: i_all_first,i_all_last, i_all_extend
            if(i_first/=0) then
                if(i_all_first == 0 &
                   .OR. &
                   i_all_first > i_first) then
                      i_all_first = i_first
                endif
            endif

            if(i_all_last == 0 &
               .OR. &
               i_all_last < i_last) then
                    i_all_last = i_last
            endif

            if( i_all_last == 0 &
               .AND. &
                i_all_first == 0) then
                    i_all_extend = 0
            else
                i_all_extend = i_all_last - i_all_first +1
            endif

        END SUBROUTINE

        SUBROUTINE set_all_cells_x_row()

            if(i_cell <= this%n_cells) then
                CALL set_first_last_extend_all_cells(i_first_filled_row, i_last_filled_row, &
                                                     this%coeff(i_basbas)%all_cells_first_filled_row_onsite,&
                                                     this%coeff(i_basbas)%all_cells_last_filled_row_onsite, &
                                                     this%coeff(i_basbas)%all_cells_extend_filled_row_onsite)

                CALL set_first_last_extend_all_cells(i_first_filled_column, i_last_filled_column, &
                                                     this%coeff(i_basbas)%all_cells_first_filled_column_onsite,&
                                                     this%coeff(i_basbas)%all_cells_last_filled_column_onsite, &
                                                     this%coeff(i_basbas)%all_cells_extend_filled_column_onsite)

            elseif(i_cell == this%n_cells+1) then
                CALL set_first_last_extend_all_cells(i_first_filled_row, i_last_filled_row, &
                                                     this%coeff(i_basbas)%all_cells_first_filled_row_offsite,&
                                                     this%coeff(i_basbas)%all_cells_last_filled_row_offsite, &
                                                     this%coeff(i_basbas)%all_cells_extend_filled_row_offsite)

!write(use_unit,*) "col",myid,i_basbas,i_first_filled_column,i_last_filled_column, this%coeff(i_basbas)%all_cells_first_filled_column_offsite


                CALL set_first_last_extend_all_cells(i_first_filled_column, i_last_filled_column, &
                                                     this%coeff(i_basbas)%all_cells_first_filled_column_offsite,&
                                                     this%coeff(i_basbas)%all_cells_last_filled_column_offsite, &
                                                     this%coeff(i_basbas)%all_cells_extend_filled_column_offsite)


!write(use_unit,*) "colafter",myid,i_basbas,i_first_filled_column,i_last_filled_column, this%coeff(i_basbas)%all_cells_first_filled_column_offsite

            endif

        END SUBROUTINE set_all_cells_x_row

        SUBROUTINE add_to_dictionary(i_value)
            INTEGER, INTENT(IN) :: i_value
            REAL*8,ALLOCATABLE,DIMENSION(:) :: tmp

            CHARACTER(*), PARAMETER :: funcname = 'add_to_dictionary'
            INTEGER :: allocation_info

            this%coeff(i_basbas)%dictionary_pos = this%coeff(i_basbas)%dictionary_pos &
                                                    +1
!if(i_basbas==11) write(use_unit,*) "add value", this%coeff(i_basbas)%dictionary_pos, i_value
            if(this%coeff(i_basbas)%dictionary_pos > this%coeff(i_basbas)%dictionary_length) then
                ALLOCATE(tmp(this%coeff(i_basbas)%dictionary_length), stat=allocation_info)
                CALL check_allocation(allocation_info, 'tmp', funcname)
!write(use_unit,*) "realloc dictionary"

                tmp(:) = this%coeff(i_basbas)%dictionary(:)
                DEALLOCATE(this%coeff(i_basbas)%dictionary)

                this%coeff(i_basbas)%dictionary_length = this%coeff(i_basbas)%dictionary_length &
                                                            + 1024
                ALLOCATE(this%coeff(i_basbas)%dictionary(this%coeff(i_basbas)%dictionary_length), stat=allocation_info)
                CALL check_allocation(allocation_info, 'dictionary', funcname)
                this%coeff(i_basbas)%dictionary(1:size(tmp)) = tmp(:)

                DEALLOCATE(tmp)
            endif

            this%coeff(i_basbas)%dictionary(this%coeff(i_basbas)%dictionary_pos) &
                                                 = i_value

        END SUBROUTINE add_to_dictionary

        SUBROUTINE store_cell_compressed_by_dictionary(this, i_basbas, i_cell, data)
        type(expansion_coefficients) :: this
        INTEGER, INTENT(IN) :: i_basbas, &
                               i_cell
        REAL*8, DIMENSION(n_basis,n_basis),INTENT(IN) :: data

        INTEGER :: i,j, &
                   i_accept, &
                   i_drop, &
                   i_pos, &
                   i_total_pos, &
                   i_dic_entry, &
                   i_run

        i_pos = this%coeff(i_basbas)%cell_start(i_cell)

        i = 0
        j = 1
        do i_dic_entry = this%coeff(i_basbas)%dictionary_start(i_cell), &
                         this%coeff(i_basbas)%dictionary_end(i_cell), 2
            i_accept = this%coeff(i_basbas)%dictionary(i_dic_entry)
            i_drop = this%coeff(i_basbas)%dictionary(i_dic_entry+1)
!write(use_unit,*) "accept",i_accept, "drop",i_drop
            do i_run = 1,i_accept
                i = i+1
                if(i>n_basis) then
                     i = 1
                     j=j+1
                endif

                this%coeff(i_basbas)%c_compressed(i_pos) =  data(i,j)
                i_pos = i_pos + 1
            enddo

            i_total_pos = (j-1)*n_basis + i
            i_total_pos = i_total_pos + i_drop

            j = 1 + (i_total_pos-1) / n_basis
            i = 1 + mod(i_total_pos-1,n_basis)
!write(use_unit,*) "tot,i,j",i_total_pos, i,j
        enddo


        END SUBROUTINE store_cell_compressed_by_dictionary

        SUBROUTINE check_c_compress_size(i_min_size)
            INTEGER, INTENT(IN) :: i_min_size

            CHARACTER(*), PARAMETER :: funcname = 'check_c_compress_size'
            INTEGER :: allocation_info
            REAL*8, ALLOCATABLE, DIMENSION(:) :: tmp
            LOGICAL :: copy_data

            copy_data = this%coeff(i_basbas)%c_compressed_size > 0

            if(this%coeff(i_basbas)%c_compressed_size < i_min_size) then

 !write(use_unit,*) "realloc basbas,cur size,min_size:" ,i_basbas, this%coeff(i_basbas)%c_compressed_size,i_min_size

                if(copy_data) then
                    ALLOCATE(tmp(this%coeff(i_basbas)%c_compressed_size), stat=allocation_info)
                    CALL check_allocation(allocation_info, 'mo_coeff_k', funcname)

                    tmp(:) = this%coeff(i_basbas)%c_compressed(:)
                    DEALLOCATE(this%coeff(i_basbas)%c_compressed)
                endif

                this%coeff(i_basbas)%c_compressed_size = i_min_size &
                                                            + 1024 * 4
                ALLOCATE(this%coeff(i_basbas)%c_compressed(this%coeff(i_basbas)%c_compressed_size), stat=allocation_info)
                CALL check_allocation(allocation_info, 'mo_coeff_k', funcname)

                if(copy_data) then
 !               write(use_unit,*) "copy data"
                    this%coeff(i_basbas)%c_compressed(1:size(tmp)) = tmp(:)
                    DEALLOCATE(tmp)
                endif
            endif


        END SUBROUTINE check_c_compress_size
    END SUBROUTINE store_cell_compressed

    SUBROUTINE decompress_cell(this, i_basbas, i_cell, i_local_or_requested, data)
        type(expansion_coefficients) :: this
        INTEGER, INTENT(IN) :: i_basbas, &
                               i_cell, &
                               i_local_or_requested
        REAL*8, DIMENSION(n_basis,n_basis),INTENT(OUT) :: data

        INTEGER :: i,j, &
                   i_accept, &
                   i_drop, &
                   i_pos, &
                   i_total_pos, &
                   i_dic_entry, &
                   i_run

        data(:,:) = 0.d0

        i_pos = this%coeff(i_basbas)%cell_start(i_cell)

!            write(use_unit,*) i_basbas, i_cell, i_pos
!write(use_unit,*) this%coeff(i_basbas)%dictionary_start(i_cell), &
!                         this%coeff(i_basbas)%dictionary_end(i_cell)

        i = 0
        j = 1
        do i_dic_entry = this%coeff(i_basbas)%dictionary_start(i_cell), &
                         this%coeff(i_basbas)%dictionary_end(i_cell), 2
            i_accept = this%coeff(i_basbas)%dictionary(i_dic_entry)
            i_drop = this%coeff(i_basbas)%dictionary(i_dic_entry+1)

            do i_run = 1,i_accept
                i = i+1
                if(i>n_basis) then
                     i = 1
                     j=j+1
                endif

                if(i_local_or_requested == EXP_COEFFICIENT_LOCAL) then
                    data(i,j) = this%coeff(i_basbas)%c_compressed(i_pos)
                elseif(i_local_or_requested == EXP_COEFFICIENT_REQUESTED) then
                    data(i,j) = this%requested_coeff_compressed(i_pos)
                endif
                i_pos = i_pos + 1
            enddo

            i_total_pos = (j-1)*n_basis + i
            i_total_pos = i_total_pos + i_drop

            j = 1 + (i_total_pos-1) / n_basis
            i = 1 + mod(i_total_pos-1,n_basis)
        enddo

    END SUBROUTINE decompress_cell

    SUBROUTINE drop_exp_coeff(this,n_basbas)
        type(expansion_coefficients) :: this
        INTEGER, INTENT(IN) :: n_basbas
        INTEGER :: pos

        if(.NOT. is_exp_coeff_stored(this,n_basbas)) then
            write(use_unit,*) "Tried to drop basbas:", n_basbas, " that is not stored"
            stop
        end if

        pos = find_exp_coeff_by_basbas(this,n_basbas)

        DEALLOCATE(this%coeff(pos)%column_norm_max_row)
        DEALLOCATE(this%coeff(pos)%column_norm_max_column)

        this%coeff(pos)%n_basbas = EXP_COEFFICIENT_UNASSIGNED

        this%coeff(pos)%last_stored_cell = 0
        this%coeff(pos)%dictionary_length = -1
        this%coeff(pos)%dictionary_pos = -1


        DEALLOCATE(this%coeff(pos)%cell_start)
        DEALLOCATE(this%coeff(pos)%cell_end)
        DEALLOCATE(this%coeff(pos)%dictionary_start)
        DEALLOCATE(this%coeff(pos)%dictionary_end)
        DEALLOCATE(this%coeff(pos)%dictionary)

        if(this%coeff(pos)%is_stored_locally) then
            DEALLOCATE(this%coeff(pos)%c_compressed)
        endif

    END SUBROUTINE drop_exp_coeff

    SUBROUTINE condense_exp_coeffs(exp_coeffs_n_cells, &
                                   exp_coeffs_n_cells_bvk)
       type(expansion_coefficients) :: exp_coeffs_n_cells, &
                                       exp_coeffs_n_cells_bvk
       INTEGER :: i_cell, i_basbas,i_k_point, i_cell_bvk
       REAL*8, DIMENSION(n_basis,n_basis) :: tmp2
       DOUBLE COMPLEX, DIMENSION(n_basis,n_basis) :: tmp,tmp_bvk

       CALL write_stdout("Doing expansion coefficients cell condensation n_cells -> n_cells_bvk")

       CALL init_exp_coeff_storage(exp_coeffs_n_cells_bvk, &
                                   exp_coeffs_n_cells%n_basbas, &
                                   n_cells_bvk, &
                                   exp_coeffs_n_cells%basis_dist)

       do i_basbas = 1, exp_coeffs_n_cells%n_basbas
           if(.NOT. exp_coeffs_n_cells%coeff(i_basbas)%is_stored_locally) CYCLE

           do i_cell_bvk = 1,n_cells_bvk
               tmp_bvk(:,:) = 0.d0

               do i_k_point = 1,n_k_points
                  tmp(:,:) = 0.d0
                  do i_cell = 1,n_cells
                     CALL decompress_cell(exp_coeffs_n_cells, i_basbas, i_cell, &
                                          EXP_COEFFICIENT_LOCAL, tmp2)
                     tmp(:,:) = tmp(:,:) + tmp2(:,:) * k_phase(i_cell,i_k_point)
                  enddo

                  tmp_bvk(:,:) = tmp_bvk(:,:) + tmp(:,:)* conjg(k_phase_exx(i_cell_bvk, i_k_point))* k_weights(i_k_point)
               enddo

!CALL write_debug("Storing basbas"//num2str(i_basbas))
               CALL store_cell_compressed(exp_coeffs_n_cells_bvk, &
                                          i_basbas, &
                                          i_cell_bvk, &
                                          dble(tmp_bvk))
           enddo

           CALL decompress_cell(exp_coeffs_n_cells, i_basbas, &
                                exp_coeffs_n_cells%n_cells+1, &
                                EXP_COEFFICIENT_LOCAL, tmp2)
           CALL store_cell_compressed(exp_coeffs_n_cells_bvk, &
                                      i_basbas, &
                                      exp_coeffs_n_cells_bvk%n_cells+1, &
                                      tmp2)

       enddo

       CALL write_stdout("Finished expansion coefficients cell condensation")

    END SUBROUTINE condense_exp_coeffs

    SUBROUTINE set_norms(this)
        type(expansion_coefficients), INTENT(INOUT) :: this

        INTEGER :: i_vector,i_vector_inversion
        type(realspace_vector) :: vector_inversion
        INTEGER :: i, n_vecs

        do i=1, ubound(this%coeff,1)
            if(this%coeff(i)%n_basbas/= EXP_COEFFICIENT_UNASSIGNED) then

                do i_vector = 1, n_vecs
                    CALL get_point_inversion_of_periodic_vector(this%rvecs%vecs(i_vector), &
                                                                this%rvecs, &
                                                                vector_inversion)

                    CALL find_vector_in_vector_array(vector_inversion, &
                                                     this%rvecs, &
                                                     i_vector_inversion)


!                    CALL get_norm_max_column(ubound(this%coeff(i)%c(:,:,i_vector),1), &
!                                             ubound(this%coeff(i)%c(:,:,i_vector),2), &
!                                             this%coeff(i)%c(:,:,i_vector), &
!                                             this%coeff(i)%column_norm_max_row(i_vector_inversion))
!
!
!                   CALL get_norm_max_row(ubound(this%coeff(i)%c(:,:,i_vector),1), &
!                                         ubound(this%coeff(i)%c(:,:,i_vector),2), &
!                                         this%coeff(i)%c(:,:,i_vector), &
!                                         this%coeff(i)%column_norm_max_column(i_vector_inversion))

                enddo
            endif
        enddo
    END SUBROUTINE set_norms

    REAL*8 FUNCTION chksum_exp_coeffs(this)
        type(expansion_coefficients) :: this

        INTEGER :: i

        chksum_exp_coeffs = 0.d0
        do i=1, ubound(this%coeff,1)
            if(this%coeff(i)%is_stored_locally) then
            chksum_exp_coeffs = chksum_exp_coeffs + &
                                sum(abs(this%coeff(i)%c_compressed(:)))
            endif
        enddo
    END FUNCTION chksum_exp_coeffs

    SUBROUTINE print_exp_coeffs_all(this)
        type(expansion_coefficients) :: this
        INTEGER :: i_basbas

        do i_basbas=1, this%n_basbas
            CALL print_exp_coeff(this,i_basbas)
        enddo
    END SUBROUTINE print_exp_coeffs_all

    SUBROUTINE print_exp_coeff(this,i_basbas)
        type(expansion_coefficients) :: this
        INTEGER, INTENT(IN) :: i_basbas
        INTEGER :: i_pos, i_vec, i_row, i_column

        i_pos = find_exp_coeff_by_basbas(this,i_basbas)

        CALL write_stdout("Basbas, pos"//num2str(i_basbas)//num2str(i_pos))

        do i_vec = 1,this%rvecs%n_vecs
            CALL print_vector(this%rvecs%vecs(i_vec))
!            do i_row = 1, ubound(this%coeff(i_basbas)%c,1)
!                do i_column = 1, ubound(this%coeff(i_basbas)%c,2)
!                    write(use_unit,*) i_row,i_column, this%coeff(i_basbas)%c(i_row,i_column,i_vec)
!                enddo
!            enddo
        enddo
    END SUBROUTINE print_exp_coeff

    SUBROUTINE set_requested_id(this, i_id_master)
        type(expansion_coefficients) :: this
        INTEGER, INTENT(IN) :: i_id_master

        this%requested_id = i_id_master
    END SUBROUTINE set_requested_id

    SUBROUTINE drop_requested_basbas(this)
        type(expansion_coefficients) :: this
        INTEGER :: i_basbas

        if(this%requested_id/=EXP_COEFFICIENT_NO_REQUESTED) then
            do i_basbas = 1,this%n_basbas
                if(this%basbas_to_id(i_basbas)==this%requested_id) then
                    DEALLOCATE(this%coeff(i_basbas)%c_compressed)
                endif
            enddo

            this%requested_id=EXP_COEFFICIENT_NO_REQUESTED
        endif
    END SUBROUTINE drop_requested_basbas

    SUBROUTINE finalize_exp_coeff_storage(this)
        type(expansion_coefficients) :: this
        INTEGER :: i

        if(this%is_exp_coeff_storage_inited) then

            do i = 1, ubound(this%coeff,1)
                if (this%coeff(i)%n_basbas/= EXP_COEFFICIENT_UNASSIGNED) then
                    CALL drop_exp_coeff(this,this%coeff(i)%n_basbas)
                end if
            end do

            DEALLOCATE(this%coeff)
            DEALLOCATE(this%loaded_column_norm_max_row)
            DEALLOCATE(this%loaded_column_norm_max_column)

            DEALLOCATE(this%basbas_to_id)

            if(ALLOCATED(this%requested_coeff_compressed)) then
                DEALLOCATE(this%requested_coeff_compressed)
            endif

            if(ALLOCATED(this%basbas_local_by_id)) then
                DEALLOCATE(this%basbas_local_by_id)
            endif

            CALL free_all_lattice_vectors(this%rvecs)

            this%is_exp_coeff_storage_inited = .FALSE.
        end if
    END SUBROUTINE finalize_exp_coeff_storage

END MODULE cRPA_storage_exp_coeffs
