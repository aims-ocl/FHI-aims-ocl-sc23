  Subroutine CC_cl_DIIS()

  Use CC_cl
 

  Implicit None

  Integer :: iii,jjj

  Call CC_cl_Jacob()

  Call CC_cl_DIIS_save('t',CC_DIIS_m_bgn)

  Call CC_cl_DIIS_Bmat()

  if ((CC_i_scf.ge.1).and.(Mod((CC_i_scf - 1),CC_DIIS_step).eq.0)) then
    Call CC_cl_DIIS_solution()
  end if

  CC_DIIS_ndc = CC_DIIS_ndc + 1

  if (CC_DIIS_ndc.gt.CC_DIIS_n_sv) then
    CC_DIIS_ndc = CC_DIIS_n_sv
  end if

  CC_DIIS_m_bgn = CC_DIIS_m_bgn + 1
  if (CC_DIIS_m_bgn.gt.CC_DIIS_ndc) then
    CC_DIIS_m_bgn = 1
  end if

  Call CC_cl_DIIS_save('r',CC_DIIS_m_bgn)

  End Subroutine CC_cl_DIIS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_DIIS_solution()

  Use CC_cl

  Implicit None

  Integer :: iii,jjj,aaa,bbb
  Integer :: a_run,n_a,n_ij
  Integer :: errnum,kkk,lll
  Integer :: NNN,INFO
  Integer , dimension(:) , allocatable :: IPIV
  Double precision , dimension(:,:) , allocatable :: Rmat,tau

  NNN = CC_DIIS_ndc + 1
  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  n_ij = CC_mem_ij_D(CC_mpi_did+1)

  ! Form R matrix and a vector
  Allocate(Rmat(NNN,NNN),stat=errnum)
  Call check_allocation(errnum, 'Rmat in CC')

  Allocate(tau(NNN,1),stat=errnum)
  Call check_allocation(errnum, 'tau in CC')

  Rmat(1,2:NNN) = -1.0D0
  Rmat(2:NNN,1) = -1.0D0

  tau(2:NNN,1) = 0.0D0

  Rmat(1,1) = 0.0D0
  tau(1,1) = -1.0D0

  do kkk = 2, NNN
    do lll = 2, NNN 
      Rmat(kkk,lll) = CC_DIIS_Bmat(kkk - 1,lll - 1)
    end do
  end do

  Allocate(IPIV(NNN),stat=errnum)
  Call check_allocation(errnum, 'IPIV in CC')

  Call DGESV(NNN,1,Rmat,NNN,IPIV,tau,NNN,INFO)

  Call CC_cl_clean_w()

  do kkk = 1, NNN - 1
 
    CC_w_s = CC_w_s + CC_DIIS_t_s(:,:,kkk,1) * tau(kkk+1,1)

    !$OMP PARALLEL Default(Shared) &
    !$OMP Private(a_run,bbb,iii)
    !$OMP DO
    do a_run = 1, n_a
      do bbb = 1, CC_n_vir
        do iii = 1, n_ij
          CC_w_d(iii,a_run,bbb,1) = CC_w_d(iii,a_run,bbb,1) &
                                  + CC_DIIS_t_d(iii,a_run,bbb,kkk) * tau(kkk+1,1)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end do

  Deallocate(Rmat,tau,IPIV)

  End Subroutine CC_cl_DIIS_solution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_DIIS_Bmat()

  Use CC_cl

  Implicit None

  Integer :: iii,jjj,kkk,lll,nnn,aaa,bbb,a_run
  Integer :: n_i,i_start,i_end,n_a,a_start,a_end,n_ij,ij_run,ij_start,ij_end,i_ij
  Double precision :: Bij,rlt

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  n_i = CC_mem_ii_D(CC_mpi_did+1)
  i_start = CC_index_ii_D(CC_mpi_did+1)
  i_end = i_start - 1 + n_i

  n_ij = CC_mem_ij_D(CC_mpi_did+1)
  ij_start = CC_index_ij_D(CC_mpi_did+1)
  ij_end = ij_start - 1 + n_ij

  do nnn = 1, CC_DIIS_ndc

    rlt = 0.0D0

    !$OMP PARALLEL Default(Shared) Reduction(+:rlt) &
    !$OMP Private(iii,aaa)
    !$OMP DO
    do iii = i_start, i_end
      do aaa = a_start, a_end
        rlt = rlt + CC_DIIS_r_s(iii,aaa,nnn,1) &
                  * CC_DIIS_r_s(iii,aaa,CC_DIIS_m_bgn,1)
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    Bij = rlt
    rlt = 0.0D0

    !$OMP PARALLEL Default(Shared) Reduction(+:rlt) &
    !$OMP Private(ij_run,a_run,i_ij,iii,jjj,aaa,bbb)
    !$OMP DO
    do ij_run = 1, n_ij
      do a_run = 1, n_a
        do bbb = 1, CC_n_vir

          aaa = a_start - 1 + a_run
          i_ij = ij_start - 1 + ij_run
          Call CC_cl_decode(i_ij,iii,jjj,CC_n_occ,2)
          
          if ((iii.ne.jjj).or.(aaa.le.bbb)) then
            rlt = rlt + CC_DIIS_r_d(ij_run,a_run,bbb,nnn) &
                      * CC_DIIS_r_d(ij_run,a_run,bbb,CC_DIIS_m_bgn)
          end if

        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    Bij = Bij + rlt

    Call CC_mpi_real_number(Bij, MPI_COMM_WORLD)

    CC_DIIS_Bmat(nnn,CC_DIIS_m_bgn) = Bij
    CC_DIIS_Bmat(CC_DIIS_m_bgn,nnn) = Bij

  end do 

 ! print*,'DIIS_B'
 ! do iii = 1, CC_DIIS_ndc
 !   print*,(CC_DIIS_Bmat(iii,jjj),jjj=1,CC_DIIS_ndc)
 ! end do

  End Subroutine CC_cl_DIIS_Bmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_DIIS_save(jobtype,i_sv)

  Use CC_cl

  Implicit None

  Character (len = 1) , intent(in) :: jobtype
  Integer , intent(in) :: i_sv

  Integer :: ij_start,ij_end,n_ij,ab_start,ab_end,n_ab
  Integer :: i_start,i_end,n_i,a_start,a_end,n_a

  Integer :: aaa,bbb,iii,jjj,i_run,ab_run,i_ab,a_run,ij_run

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  n_ij = CC_mem_ij_D(CC_mpi_did+1)

  if (jobtype.eq.'t') then

    !write(70+myid,*) 'i_sv',i_sv

    CC_DIIS_t_s(:,:,i_sv,1) = CC_w_s
    CC_DIIS_t_d(:,:,:,i_sv) = CC_w_d(:,:,:,1)

    CC_DIIS_r_s(:,:,i_sv,1) = CC_w_s - CC_DIIS_r_s(:,:,i_sv,1)

    !$OMP PARALLEL Default(Shared) &
    !$OMP Private(aaa,bbb,ij_run)
    !$OMP DO
    do bbb = 1, CC_n_vir
      do aaa = 1, n_a
        do ij_run = 1, n_ij
          CC_DIIS_r_d(ij_run,aaa,bbb,i_sv) = CC_w_d(ij_run,aaa,bbb,1) &
                                           - CC_DIIS_r_d(ij_run,aaa,bbb,i_sv)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  else if (jobtype.eq.'r') then

    CC_DIIS_r_s(:,:,i_sv,1) = CC_w_s
    CC_DIIS_r_d(:,:,:,i_sv) = CC_w_d(:,:,:,1)

  end if

  End Subroutine CC_cl_DIIS_save

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
