!Calculates continuous and discrete fourier transform
MODULE cRPA_calculation_fouriertrans
    USE cRPA_flow_realspace_vector
    USE cRPA_parallelism_storage
    USE cRPA_calculation_integration
    implicit none

CONTAINS

    SUBROUTINE fouriertransform_x_to_t(int_grid, &
                                       n_elements, f_of_x, &
                                       t, f_of_t)
        type(integration_grid), INTENT(IN) :: int_grid
        INTEGER :: n_elements
        REAL*8,DIMENSION(n_elements,int_grid%n_points), INTENT(IN) :: f_of_x
        REAL*8, INTENT(IN) :: t
        COMPLEX*16,DIMENSION(n_elements), INTENT(OUT) :: f_of_t

        REAL*8 :: pi
        COMPLEX*16,DIMENSION(int_grid%n_points) :: weights_fourier


        pi = acos(-1.d0)

        weights_fourier(:) = DCMPLX(cos(int_grid%points(:)*t),0.d0) &!exp( CMPLX(0.d0,&
                             !           1.d0*int_grid%points(:)*t) ) &
                             * int_grid%weights_complex(:)
        weights_fourier(:) = weights_fourier(:)

!        f_of_t(1) = sum(f_of_x(:,1)*weights_fourier(:))
         f_of_t = MATMUL(f_of_x,weights_fourier)
!         f_of_t(:) = f_of_t(:) / (2.d0*pi)

!        CALL zgemv("T", &
!                   int_grid%n_points,n_elements, &
!                   1.d0, &
!                   DCMPLX(f_of_x(:,1),0.d0), int_grid%n_points, &
!                   weights_fourier, 1, &
!                   0.d0, &
!                   f_of_t,1)
!        write(use_unit,*) t,f_of_x(1,1),f_of_t
    END SUBROUTINE fouriertransform_x_to_t

    SUBROUTINE fouriertransform_t_to_x(int_grid,n_elements, f_of_t,x, f_of_x)
        type(integration_grid), INTENT(IN) :: int_grid
        INTEGER, INTENT(IN) :: n_elements
        COMPLEX*16,DIMENSION(n_elements,int_grid%n_points), INTENT(IN) :: f_of_t
        REAL*8, INTENT(IN) :: x
        COMPLEX*16,DIMENSION(n_elements), INTENT(OUT) :: f_of_x

        COMPLEX*16 :: tmp
        REAL*8 :: pi
        COMPLEX*16,DIMENSION(int_grid%n_points) :: weights_fourier
        INTEGER :: i

        pi = acos(-1.d0)

        f_of_x(:) = 0.d0

        weights_fourier(:) = exp( DCMPLX(0.d0,&
                                        -1.d0*int_grid%points(:)*x) ) &
                                        * int_grid%weights_complex(:)
        do i=1, n_elements
            tmp = sum(f_of_t(i,:)*weights_fourier(:))
            f_of_x(i) = tmp!real(tmp)
        enddo
         !   write(use_unit,*) tmp,f_of_x(j)

         f_of_x(:) = f_of_x(:) / (2.d0 * pi)
    END SUBROUTINE fouriertransform_t_to_x


    SUBROUTINE discrete_fouriertransform_k_to_R(rvecs,n_elements, f_of_k, f_of_R)
        type(realspace_vectors),INTENT(IN) :: rvecs
        INTEGER, INTENT(IN) :: n_elements
        COMPLEX*16, DIMENSION(n_elements,n_k_points),INTENT(IN) :: f_of_k
        COMPLEX*16, DIMENSION(n_elements,n_k_points),INTENT(OUT) :: f_of_R

        REAL*8,DIMENSION(n_k_points) :: tmp_real, tmp_img

        INTEGER :: i,n, i_cell_bvk,i_vec, i_k_point

        do i_vec = 1, rvecs%n_vecs
            do n=1, n_elements
               i_cell_bvk = rvecs%vecs(i_vec)%i_cell_bvk

               tmp_real(:) = REAL(f_of_k(n,:) * conjg(k_phase_exx(i_cell_bvk,:))*k_weights(:),8)
               tmp_img(:) = dimag(f_of_k(n,:) * conjg(k_phase_exx(i_cell_bvk,:))*k_weights(:))

               f_of_R(n,i_cell_bvk) = DCMPLX(sort_and_sum_vector(n_k_points,tmp_real), &
                                            sort_and_sum_vector(n_k_points,tmp_img))
           enddo
       enddo

!       do i_vec = 1, rvecs%n_vecs
!           f_of_R(:,i_vec) = f_of_R(:,i_vec) * k_weights(1)
!       enddo
!write(use_unit,*) k_weights*real(rvecs%n_vecs)
!stop
    END SUBROUTINE discrete_fouriertransform_k_to_R

    SUBROUTINE discrete_real_fouriertransform_R_to_k(rvecs,n_elements,f_of_R_real, f_of_k )
        type(realspace_vectors),INTENT(IN) :: rvecs
        INTEGER, INTENT(IN) :: n_elements
        REAL*8, DIMENSION(n_elements,rvecs%n_vecs),INTENT(IN) :: f_of_R_real
        COMPLEX*16, DIMENSION(n_elements,n_k_points),INTENT(OUT) :: f_of_k

        COMPLEX*16, DIMENSION(n_elements,rvecs%n_vecs) :: f_of_R

        f_of_R(:,:) = f_of_R_real(:,:)

        CALL discrete_fouriertransform_R_to_k(rvecs,n_elements,f_of_R, f_of_k )
    END SUBROUTINE

    SUBROUTINE discrete_fouriertransform_R_to_k(rvecs,n_elements,f_of_R, f_of_k )
        type(realspace_vectors),INTENT(IN) :: rvecs
        INTEGER, INTENT(IN) :: n_elements
        COMPLEX*16, DIMENSION(n_elements,rvecs%n_vecs),INTENT(IN) :: f_of_R
        COMPLEX*16, DIMENSION(n_elements,n_k_points),INTENT(OUT) :: f_of_k

        INTEGER :: i,n,i_cell_bvk,i_vec, i_k_point
        REAL*8,DIMENSION(n_k_points) :: tmp_real, tmp_img


        if(n_k_points /= rvecs%n_vecs) then
            stop "Can not do discrete fouriertransform for unequal number of k and R vecs!"
        endif

        f_of_k(:,:) = DCMPLX(0.d0,0.d0)

        do i_k_point = 1, n_k_points
            do n=1,n_elements
                do i_vec = 1, rvecs%n_vecs
                    i_cell_bvk = rvecs%vecs(i_vec)%i_cell_bvk

                    f_of_k(n,i_k_point) = f_of_k(n,i_k_point) + &
                                        f_of_R(n,i_vec) * &
                                        k_phase_exx(i_cell_bvk,i_k_point)!* k_weights(i_k_point)

                   tmp_real(i_vec) = REAL(f_of_R(n,i_vec) * k_phase_exx(i_cell_bvk,i_k_point),8)
                   tmp_img(i_vec) = dimag(f_of_R(n,i_vec) * k_phase_exx(i_cell_bvk,i_k_point))

!                   write(use_unit,*) tmp_real(i_vec),tmp_img(i_vec)

                enddo

                f_of_k(n,i_k_point) = DCMPLX(sort_and_sum_vector(n_k_points,tmp_real), &
                                            sort_and_sum_vector(n_k_points,tmp_img))
!                write(use_unit,*) n,i_k_point, f_of_k(n,i_k_point)
            enddo
       enddo

!       f_of_k(:,:) = f_of_k(:,:) !* dble(n_k_points)

    END SUBROUTINE discrete_fouriertransform_R_to_k
END MODULE cRPA_calculation_fouriertrans
