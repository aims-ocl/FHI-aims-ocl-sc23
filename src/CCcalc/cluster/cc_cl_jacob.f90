!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_E_corr()

  Use CC_cl

  Implicit None

  Integer :: iii,jjj,aaa,bbb,a_start,a_end,n_a,b_start,b_end,n_b,a_run
  Double precision :: rlt,Nrlt

  CC_E_corr = 0.0D0
  CC_Norm = 0.0D0

  rlt = 0.0D0
  Nrlt = 0.0D0

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  n_b = CC_mem_aa_G(CC_mpi_gid+1)
  b_start = CC_index_aa_G(CC_mpi_gid+1)
  b_end = b_start - 1 + n_b

  !write(myid+80,*) 'w_d'
  !do iii = 1, CC_n_occ
  !  do jjj = 1, CC_n_occ
  !    do a_run = 1, n_a
  !      do bbb = 1, CC_n_vir
  !        aaa = a_start - 1 + a_run
  !        write(80+myid,*) iii,jjj,aaa,bbb,CC_t_d(iii,jjj,a_run,bbb)
  !      end do
  !    end do
  !  end do
  !end do

  !$OMP PARALLEL Default(Shared) Reduction(+:rlt,Nrlt) &
  !$OMP Private(a_run,iii,jjj,aaa,bbb)
  !$OMP DO
  do iii = 1, CC_n_occ
    do a_run = 1, n_a
      do jjj = 1, CC_n_occ
        do bbb = b_start, b_end

          aaa = a_start - 1 + a_run

          rlt = rlt + (2.0D0 * CC_intl_iajb(iii,jjj,a_run,bbb)  &
                             - CC_intl_iajb(jjj,iii,a_run,bbb)) &
                    * (CC_t_d(iii,jjj,a_run,bbb) + CC_t_s(iii,aaa) * CC_t_s(jjj,bbb))

          Nrlt = Nrlt + CC_t_d(iii,jjj,a_run,bbb) ** 2
          Nrlt = Nrlt + 0.5D0 * (CC_t_d(jjj,jjj,a_run,bbb) &
                                 - CC_t_d(iii,jjj,a_run,bbb)) ** 2
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  CC_E_corr = rlt
  CC_Norm = Nrlt

  Call CC_mpi_real_number(CC_E_corr, MPI_COMM_WORLD)
  Call CC_mpi_real_number(CC_Norm, MPI_COMM_WORLD)

  do iii = 1, CC_n_occ
    do aaa = 1, CC_n_vir
      CC_Norm = CC_Norm + CC_t_s(iii,aaa) ** 2
    end do
  end do

  CC_Norm = sqrt(1.0D0 + CC_Norm)

  End Subroutine CC_cl_E_corr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_Jacob()

  Use CC_cl

  Implicit None

  Integer :: iii,jjj,aaa,bbb
  Integer :: ij_start,ij_end,n_ij,ij_run,i_ij
  Integer :: n_a,a_run,a_start,a_end,i_a,n_b,b_run,b_start,b_end
  Double precision :: CC_E_diff_0_s,et1,et2

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  n_ij = CC_mem_ij_D(CC_mpi_did+1)
  ij_start = CC_index_ij_D(CC_mpi_did+1)

  ! For singlet excitation
  !$OMP PARALLEL Default(Shared) &
  !$OMP Private(iii,aaa,et1)
  !$OMP DO
  do iii = 1, CC_n_occ
    do aaa = 1, CC_n_vir
    
      Call CC_cl_ev_diff(iii,aaa,et1)
      CC_w_s(iii,aaa) = - CC_w_s(iii,aaa) / et1

    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! For double excitation
  !$OMP PARALLEL Default(Shared) &
  !$OMP Private(a_run,ij_run,i_ij,jjj,iii,aaa,bbb,et1,et2,CC_E_diff_0_s)
  !$OMP DO
  do a_run = 1, n_a
    do bbb = 1, CC_n_vir
      do ij_run = 1, n_ij

        aaa = a_start - 1 + a_run
        i_ij = ij_start - 1 + ij_run

        Call CC_cl_decode(i_ij,iii,jjj,CC_n_occ,2)

        Call CC_cl_ev_diff(iii,aaa,et1)
        Call CC_cl_ev_diff(jjj,bbb,et2)
        CC_E_diff_0_s = - et1 - et2

        CC_w_d(ij_run,a_run,bbb,1) = CC_w_d(ij_run,a_run,bbb,1) / CC_E_diff_0_s

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  !print*,'w_s jacob'
  !  if (myid.eq.0) then
  !    print*,'w_s'
  !    do iii = 1, CC_n_occ
  !      do aaa = 1, CC_n_vir
  !        print*,CC_w_s(iii,aaa)
  !      end do
  !    end do
  !  end if

  !write(90+myid,*) 'w_d jacob'
  !do a_run = 1, CC_mem_aa_G(CC_mpi_gid+1)
  !  do bbb = 1, CC_n_vir
  !    do ij_run = 1, CC_mem_ij_D(CC_mpi_did+1)
  !      aaa = CC_index_aa_G(CC_mpi_gid+1) - 1 + a_run
  !      Call CC_cl_decode(ij_run,iii,jjj,CC_n_occ,2)
  !      write(90+myid,*) iii,jjj,aaa,bbb,CC_w_d(ij_run,a_run,bbb,1)
  !    end do
  !  end do
  !end do


  End Subroutine CC_cl_Jacob

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

