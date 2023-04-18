!Storage of molecular orbitals C
module cRPA_storage_mo_coeffs
    USE cRPA_view
    USE cRPA_parallelism_storage
    USE cRPA_storage_norms
    IMPLICIT NONE

    type mo_coefficients

        !(n_my_basis,n_states,n_k_points,n_spin)
        DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: c
        REAL*8, ALLOCATABLE, DIMENSION(:,:) :: c_norm_row, &
                                               c_norm_column

        INTEGER :: n_k_points, n_k_points_my, n_basis, n_states, n_spin
        INTEGER,DIMENSION(:), ALLOCATABLE :: k_point_to_id, &
                                             k_point_to_k_point_local
        LOGICAL is_mo_coeffs_storage_inited
    end type

    type(mo_coefficients) :: mo_coeffs_local, &
                             mo_coeffs_requested
    INTEGER :: MO_UNASSIGNED
    PARAMETER(MO_UNASSIGNED = -1)
CONTAINS

    SUBROUTINE init_mo_coefficient_storage(this,n_k_points, &
                                           n_basis, n_states, n_spin)
        type(mo_coefficients), INTENT(INOUT) :: this
        INTEGER,INTENT(IN) :: n_k_points , n_basis
        INTEGER,INTENT(IN) :: n_states, n_spin

        if(this%is_mo_coeffs_storage_inited) then
            write(use_unit,*) "MO coeffcient storage already inited!"
            stop
        end if


!TODO check allocation
        this%n_k_points = n_k_points
        this%n_basis = n_basis
        this%n_states = n_states
        this%n_spin = n_spin

        this%n_k_points_my = MO_UNASSIGNED

        ALLOCATE(this%k_point_to_id(this%n_k_points))
        ALLOCATE(this%k_point_to_k_point_local(this%n_k_points))

        this%k_point_to_k_point_local(:) = MO_UNASSIGNED

        this%is_mo_coeffs_storage_inited = .TRUE.

    END SUBROUTINE init_mo_coefficient_storage

    SUBROUTINE allocate_for_my_k_points(this, n_k_points_my)
        type(mo_coefficients), INTENT(INOUT) :: this
        INTEGER,INTENT(IN) :: n_k_points_my

        this%n_k_points_my = n_k_points_my

        ALLOCATE(this%c(this%n_basis,this%n_states,this%n_spin,this%n_k_points_my))
        ALLOCATE(this%c_norm_row(this%n_k_points_my,this%n_spin))
        ALLOCATE(this%c_norm_column(this%n_k_points_my,this%n_spin))

    END SUBROUTINE allocate_for_my_k_points

    SUBROUTINE set_mo_coeffs_norms(this)
        type(mo_coefficients), INTENT(INOUT) :: this

        INTEGER :: i_k_point, n_my_basis, i_spin

        n_my_basis = this%n_basis

        do i_spin = 1, this%n_spin
            do i_k_point = 1, this%n_k_points_my
                CALL get_norm_max_row(n_my_basis,n_states, abs(this%c(:,:,i_k_point,i_spin)), &
                                      this%c_norm_row(i_k_point,i_spin))

                CALL get_norm_max_column(n_my_basis,n_states, abs(this%c(:,:,i_k_point,i_spin)), &
                                         this%c_norm_column(i_k_point,i_spin))
            enddo
        enddo
    END SUBROUTINE set_mo_coeffs_norms


    REAL*8 FUNCTION chksum_mo_coeffs(this)
        type(mo_coefficients), INTENT(INOUT) :: this
        REAL*8 :: chksum
        INTEGER :: i_spin

        chksum = 0.d0
        do i_spin = 1,this%n_spin
            chksum = chksum + &
                     chksum_3d_matrix(abs(this%c(:,:,:,i_spin)))
        enddo

        chksum_mo_coeffs = chksum
    END FUNCTION chksum_mo_coeffs


    SUBROUTINE print_mo_coeffs(this)
        type(mo_coefficients), INTENT(INOUT) :: this


    END SUBROUTINE print_mo_coeffs

    SUBROUTINE finalize_mo_coeffs_storage(this)
        type(mo_coefficients), INTENT(INOUT) :: this

        if(this%n_k_points_my/=MO_UNASSIGNED) then
            DEALLOCATE(this%c)
            DEALLOCATE(this%c_norm_row)
            DEALLOCATE(this%c_norm_column)
        end if

        DEALLOCATE(this%k_point_to_id)
        DEALLOCATE(this%k_point_to_k_point_local)

        this%is_mo_coeffs_storage_inited = .FALSE.

    END SUBROUTINE finalize_mo_coeffs_storage
end module cRPA_storage_mo_coeffs
