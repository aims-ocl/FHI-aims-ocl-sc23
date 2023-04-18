!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium, 
!    consisting of the following organizations:
!
!    - Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG), 
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen , 
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie, 
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaftrn, 
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition, 
!      and  
!    - IBM Deutschland GmbH
!
!
!    More information can be found here:
!    http://elpa.rzg.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the 
!    GNU Lesser General Public License as published by the Free 
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
!
! --------------------------------------------------------------------------------------------------
!
! This file contains the compute intensive kernels for the Householder transformations.
!
! *** Special IBM BlueGene/Q version with QPX intrinsics in Fortran ***
! 
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".
!
! --------------------------------------------------------------------------------------------------

subroutine double_hh_trafo(q, hh, nb, nq, ldq, ldh)

   implicit none

   integer, intent(in) :: nb, nq, ldq, ldh
   real*8, intent(inout) :: q(ldq,*)
   real*8, intent(in) :: hh(ldh,*)

   real*8 s
   integer i

   ! Safety only:

   if(mod(ldq,4) /= 0) STOP 'double_hh_trafo: ldq not divisible by 4!'

   call alignx(32,q)

   ! Calculate dot product of the two Householder vectors

   s = hh(2,2)*1
   do i=3,nb
      s = s+hh(i,2)*hh(i-1,1)
   enddo
   
   do i=1,nq-20,24
      call hh_trafo_kernel_24_bgq(q(i   ,1), hh, nb, ldq, ldh, s)
   enddo
   
   if(nq-i+1 > 16) then
      call hh_trafo_kernel_16_bgq(q(i  ,1), hh, nb, ldq, ldh, s)
      call hh_trafo_kernel_4_bgq(q(i+16,1), hh, nb, ldq, ldh, s)
   else if(nq-i+1 > 12) then
      call hh_trafo_kernel_8_bgq(q(i  ,1), hh, nb, ldq, ldh, s)
      call hh_trafo_kernel_8_bgq(q(i+8,1), hh, nb, ldq, ldh, s)
   else if(nq-i+1 > 8) then
      call hh_trafo_kernel_8_bgq(q(i  ,1), hh, nb, ldq, ldh, s)
      call hh_trafo_kernel_4_bgq(q(i+8,1), hh, nb, ldq, ldh, s)
   else if(nq-i+1 > 4) then
      call hh_trafo_kernel_8_bgq(q(i  ,1), hh, nb, ldq, ldh, s)
   else if(nq-i+1 > 0) then
      call hh_trafo_kernel_4_bgq(q(i  ,1), hh, nb, ldq, ldh, s)
   endif
   
end


! --------------------------------------------------------------------------------------------------
! The following kernels perform the Householder transformation on Q for 24/16/8/4 rows.
! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_kernel_24_bgq(q, hh, nb, ldq, ldh, s)

   implicit none

   include 'mpif.h'
   
   integer, intent(in) :: nb, ldq, ldh
   
   real*8, intent(inout)::q(ldq,*)
   real*8, intent(in)::hh(ldh,*), s
   
   VECTOR(REAL(8))::QPX_x1, QPX_x2, QPX_x3, QPX_x4, QPX_x5, QPX_x6
   VECTOR(REAL(8))::QPX_y1, QPX_y2, QPX_y3, QPX_y4, QPX_y5, QPX_y6
   VECTOR(REAL(8))::QPX_q1, QPX_q2, QPX_q3, QPX_q4, QPX_q5, QPX_q6
   VECTOR(REAL(8))::QPX_h1, QPX_h2, QPX_tau1, QPX_tau2, QPX_s
   integer i

   call alignx(32,q)

   !--- multiply Householder vectors with matrix q ---

   QPX_x1 = VEC_LD(0,q(1,2))
   QPX_x2 = VEC_LD(0,q(5,2))
   QPX_x3 = VEC_LD(0,q(9,2))
   QPX_x4 = VEC_LD(0,q(13,2))
   QPX_x5 = VEC_LD(0,q(17,2))
   QPX_x6 = VEC_LD(0,q(21,2))
   
   QPX_h2 = VEC_SPLATS(hh(2,2))
   QPX_q1 = VEC_LD(0,q(1,1))
   QPX_q2 = VEC_LD(0,q(5,1))
   QPX_q3 = VEC_LD(0,q(9,1))
   QPX_q4 = VEC_LD(0,q(13,1))
   QPX_q5 = VEC_LD(0,q(17,1))
   QPX_q6 = VEC_LD(0,q(21,1))
   QPX_y1 = VEC_MADD(QPX_x1, QPX_h2, QPX_q1)
   QPX_y2 = VEC_MADD(QPX_x2, QPX_h2, QPX_q2)
   QPX_y3 = VEC_MADD(QPX_x3, QPX_h2, QPX_q3)
   QPX_y4 = VEC_MADD(QPX_x4, QPX_h2, QPX_q4)
   QPX_y5 = VEC_MADD(QPX_x5, QPX_h2, QPX_q5)
   QPX_y6 = VEC_MADD(QPX_x6, QPX_h2, QPX_q6)

   do i=3,nb,1

      QPX_q1 = VEC_LD(0,q(1,i))
      QPX_q2 = VEC_LD(0,q(5,i))
      QPX_q3 = VEC_LD(0,q(9,i))
      QPX_q4 = VEC_LD(0,q(13,i))
      QPX_q5 = VEC_LD(0,q(17,i))
      QPX_q6 = VEC_LD(0,q(21,i))
      QPX_h1 = VEC_SPLATS(hh(i-1,1))
      QPX_x1 = VEC_MADD(QPX_q1, QPX_h1, QPX_x1)
      QPX_x2 = VEC_MADD(QPX_q2, QPX_h1, QPX_x2)
      QPX_x3 = VEC_MADD(QPX_q3, QPX_h1, QPX_x3)
      QPX_x4 = VEC_MADD(QPX_q4, QPX_h1, QPX_x4)
      QPX_x5 = VEC_MADD(QPX_q5, QPX_h1, QPX_x5)
      QPX_x6 = VEC_MADD(QPX_q6, QPX_h1, QPX_x6)
      QPX_h2 = VEC_SPLATS(hh(i,2))
      QPX_y1 = VEC_MADD(QPX_q1, QPX_h2, QPX_y1)
      QPX_y2 = VEC_MADD(QPX_q2, QPX_h2, QPX_y2)
      QPX_y3 = VEC_MADD(QPX_q3, QPX_h2, QPX_y3)
      QPX_y4 = VEC_MADD(QPX_q4, QPX_h2, QPX_y4)
      QPX_y5 = VEC_MADD(QPX_q5, QPX_h2, QPX_y5)
      QPX_y6 = VEC_MADD(QPX_q6, QPX_h2, QPX_y6)

   enddo

   QPX_h1 = VEC_SPLATS(hh(nb,1))
   QPX_q1 = VEC_LD(0,q(1,nb+1))
   QPX_q2 = VEC_LD(0,q(5,nb+1))
   QPX_q3 = VEC_LD(0,q(9,nb+1))
   QPX_q4 = VEC_LD(0,q(13,nb+1))
   QPX_q5 = VEC_LD(0,q(17,nb+1))
   QPX_q6 = VEC_LD(0,q(21,nb+1))
   QPX_x1 = VEC_MADD(QPX_q1, QPX_h1, QPX_x1)
   QPX_x2 = VEC_MADD(QPX_q2, QPX_h1, QPX_x2)
   QPX_x3 = VEC_MADD(QPX_q3, QPX_h1, QPX_x3)
   QPX_x4 = VEC_MADD(QPX_q4, QPX_h1, QPX_x4)
   QPX_x5 = VEC_MADD(QPX_q5, QPX_h1, QPX_x5)
   QPX_x6 = VEC_MADD(QPX_q6, QPX_h1, QPX_x6)
   
   !--- multiply T matrix ---
   
   QPX_tau1 = VEC_SPLATS(-hh(1,1))
   QPX_x1 = VEC_MUL(QPX_x1, QPX_tau1)
   QPX_x2 = VEC_MUL(QPX_x2, QPX_tau1)
   QPX_x3 = VEC_MUL(QPX_x3, QPX_tau1)
   QPX_x4 = VEC_MUL(QPX_x4, QPX_tau1)
   QPX_x5 = VEC_MUL(QPX_x5, QPX_tau1)
   QPX_x6 = VEC_MUL(QPX_x6, QPX_tau1)
   QPX_tau2 = VEC_SPLATS(-hh(1,2))
   QPX_s = VEC_SPLATS(-hh(1,2)*s)
   QPX_y1 = VEC_MUL(QPX_y1, QPX_tau2)
   QPX_y2 = VEC_MUL(QPX_y2, QPX_tau2)
   QPX_y3 = VEC_MUL(QPX_y3, QPX_tau2)
   QPX_y4 = VEC_MUL(QPX_y4, QPX_tau2)
   QPX_y5 = VEC_MUL(QPX_y5, QPX_tau2)
   QPX_y6 = VEC_MUL(QPX_y6, QPX_tau2)
   QPX_y1 = VEC_MADD(QPX_x1, QPX_s, QPX_y1)
   QPX_y2 = VEC_MADD(QPX_x2, QPX_s, QPX_y2)
   QPX_y3 = VEC_MADD(QPX_x3, QPX_s, QPX_y3)
   QPX_y4 = VEC_MADD(QPX_x4, QPX_s, QPX_y4)
   QPX_y5 = VEC_MADD(QPX_x5, QPX_s, QPX_y5)
   QPX_y6 = VEC_MADD(QPX_x6, QPX_s, QPX_y6)

   !--- rank-2 update of q ---

   QPX_q1 = VEC_LD(0,q(1,1))
   QPX_q2 = VEC_LD(0,q(5,1))
   QPX_q3 = VEC_LD(0,q(9,1))
   QPX_q4 = VEC_LD(0,q(13,1))
   QPX_q5 = VEC_LD(0,q(17,1))
   QPX_q6 = VEC_LD(0,q(21,1))
   QPX_q1 = VEC_ADD(QPX_q1, QPX_y1)
   QPX_q2 = VEC_ADD(QPX_q2, QPX_y2)
   QPX_q3 = VEC_ADD(QPX_q3, QPX_y3)
   QPX_q4 = VEC_ADD(QPX_q4, QPX_y4)
   QPX_q5 = VEC_ADD(QPX_q5, QPX_y5)
   QPX_q6 = VEC_ADD(QPX_q6, QPX_y6)
   call VEC_ST(QPX_q1, 0, q(1,1))
   call VEC_ST(QPX_q2, 0, q(5,1))
   call VEC_ST(QPX_q3, 0, q(9,1))
   call VEC_ST(QPX_q4, 0, q(13,1))
   call VEC_ST(QPX_q5, 0, q(17,1))
   call VEC_ST(QPX_q6, 0, q(21,1))
   
   QPX_h2 = VEC_SPLATS(hh(2,2))
   QPX_q1 = VEC_LD(0,q(1,2))
   QPX_q2 = VEC_LD(0,q(5,2))
   QPX_q3 = VEC_LD(0,q(9,2))
   QPX_q4 = VEC_LD(0,q(13,2))
   QPX_q5 = VEC_LD(0,q(17,2))
   QPX_q6 = VEC_LD(0,q(21,2))
   QPX_q1 = VEC_MADD(QPX_y1, QPX_h2, QPX_q1)
   QPX_q2 = VEC_MADD(QPX_y2, QPX_h2, QPX_q2)
   QPX_q3 = VEC_MADD(QPX_y3, QPX_h2, QPX_q3)
   QPX_q4 = VEC_MADD(QPX_y4, QPX_h2, QPX_q4)
   QPX_q5 = VEC_MADD(QPX_y5, QPX_h2, QPX_q5)
   QPX_q6 = VEC_MADD(QPX_y6, QPX_h2, QPX_q6)
   QPX_q1 = VEC_ADD(QPX_q1, QPX_x1)
   QPX_q2 = VEC_ADD(QPX_q2, QPX_x2)
   QPX_q3 = VEC_ADD(QPX_q3, QPX_x3)
   QPX_q4 = VEC_ADD(QPX_q4, QPX_x4)
   QPX_q5 = VEC_ADD(QPX_q5, QPX_x5)
   QPX_q6 = VEC_ADD(QPX_q6, QPX_x6)
   call VEC_ST(QPX_q1, 0, q(1,2))
   call VEC_ST(QPX_q2, 0, q(5,2))
   call VEC_ST(QPX_q3, 0, q(9,2))
   call VEC_ST(QPX_q4, 0, q(13,2))
   call VEC_ST(QPX_q5, 0, q(17,2))
   call VEC_ST(QPX_q6, 0, q(21,2))
   
   do i=3,nb,1

      QPX_q1 = VEC_LD(0,q(1,i))
      QPX_q2 = VEC_LD(0,q(5,i))
      QPX_q3 = VEC_LD(0,q(9,i))
      QPX_q4 = VEC_LD(0,q(13,i))
      QPX_q5 = VEC_LD(0,q(17,i))
      QPX_q6 = VEC_LD(0,q(21,i))
      QPX_h1 = VEC_SPLATS(hh(i-1,1))
      QPX_q1 = VEC_MADD(QPX_x1, QPX_h1, QPX_q1)
      QPX_q2 = VEC_MADD(QPX_x2, QPX_h1, QPX_q2)
      QPX_q3 = VEC_MADD(QPX_x3, QPX_h1, QPX_q3)
      QPX_q4 = VEC_MADD(QPX_x4, QPX_h1, QPX_q4)
      QPX_q5 = VEC_MADD(QPX_x5, QPX_h1, QPX_q5)
      QPX_q6 = VEC_MADD(QPX_x6, QPX_h1, QPX_q6)
      QPX_h2 = VEC_SPLATS(hh(i,2))
      QPX_q1 = VEC_MADD(QPX_y1, QPX_h2, QPX_q1)
      QPX_q2 = VEC_MADD(QPX_y2, QPX_h2, QPX_q2)
      QPX_q3 = VEC_MADD(QPX_y3, QPX_h2, QPX_q3)
      QPX_q4 = VEC_MADD(QPX_y4, QPX_h2, QPX_q4)
      QPX_q5 = VEC_MADD(QPX_y5, QPX_h2, QPX_q5)
      QPX_q6 = VEC_MADD(QPX_y6, QPX_h2, QPX_q6)
      
      call VEC_ST(QPX_q1, 0, q(1,i))
      call VEC_ST(QPX_q2, 0, q(5,i))
      call VEC_ST(QPX_q3, 0, q(9,i))
      call VEC_ST(QPX_q4, 0, q(13,i))
      call VEC_ST(QPX_q5, 0, q(17,i))
      call VEC_ST(QPX_q6, 0, q(21,i))

   enddo

   QPX_h1 = VEC_SPLATS(hh(nb,1))
   QPX_q1 = VEC_LD(0,q(1,nb+1))
   QPX_q2 = VEC_LD(0,q(5,nb+1))
   QPX_q3 = VEC_LD(0,q(9,nb+1))
   QPX_q4 = VEC_LD(0,q(13,nb+1))
   QPX_q5 = VEC_LD(0,q(17,nb+1))
   QPX_q6 = VEC_LD(0,q(21,nb+1))
   QPX_q1 = VEC_MADD(QPX_x1, QPX_h1, QPX_q1)
   QPX_q2 = VEC_MADD(QPX_x2, QPX_h1, QPX_q2)
   QPX_q3 = VEC_MADD(QPX_x3, QPX_h1, QPX_q3)
   QPX_q4 = VEC_MADD(QPX_x4, QPX_h1, QPX_q4)
   QPX_q5 = VEC_MADD(QPX_x5, QPX_h1, QPX_q5)
   QPX_q6 = VEC_MADD(QPX_x6, QPX_h1, QPX_q6)
   call VEC_ST(QPX_q1, 0, q(1,nb+1))
   call VEC_ST(QPX_q2, 0, q(5,nb+1))
   call VEC_ST(QPX_q3, 0, q(9,nb+1))
   call VEC_ST(QPX_q4, 0, q(13,nb+1))
   call VEC_ST(QPX_q5, 0, q(17,nb+1))
   call VEC_ST(QPX_q6, 0, q(21,nb+1))

end subroutine

! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_kernel_16_bgq(q, hh, nb, ldq, ldh, s)

   implicit none

   include 'mpif.h'
   
   integer, intent(in) :: nb, ldq, ldh
   
   real*8, intent(inout)::q(ldq,*)
   real*8, intent(in)::hh(ldh,*), s
   
   VECTOR(REAL(8))::QPX_x1, QPX_x2, QPX_x3, QPX_x4
   VECTOR(REAL(8))::QPX_y1, QPX_y2, QPX_y3, QPX_y4
   VECTOR(REAL(8))::QPX_q1, QPX_q2, QPX_q3, QPX_q4
   VECTOR(REAL(8))::QPX_h1, QPX_h2, QPX_tau1, QPX_tau2, QPX_s
   integer i

   call alignx(32,q)

   !--- multiply Householder vectors with matrix q ---

   QPX_x1 = VEC_LD(0,q(1,2))
   QPX_x2 = VEC_LD(0,q(5,2))
   QPX_x3 = VEC_LD(0,q(9,2))
   QPX_x4 = VEC_LD(0,q(13,2))
   
   QPX_h2 = VEC_SPLATS(hh(2,2))
   QPX_q1 = VEC_LD(0,q(1,1))
   QPX_q2 = VEC_LD(0,q(5,1))
   QPX_q3 = VEC_LD(0,q(9,1))
   QPX_q4 = VEC_LD(0,q(13,1))
   QPX_y1 = VEC_MADD(QPX_x1, QPX_h2, QPX_q1)
   QPX_y2 = VEC_MADD(QPX_x2, QPX_h2, QPX_q2)
   QPX_y3 = VEC_MADD(QPX_x3, QPX_h2, QPX_q3)
   QPX_y4 = VEC_MADD(QPX_x4, QPX_h2, QPX_q4)

   do i=3,nb,1

      QPX_q1 = VEC_LD(0,q(1,i))
      QPX_q2 = VEC_LD(0,q(5,i))
      QPX_q3 = VEC_LD(0,q(9,i))
      QPX_q4 = VEC_LD(0,q(13,i))
      QPX_h1 = VEC_SPLATS(hh(i-1,1))
      QPX_x1 = VEC_MADD(QPX_q1, QPX_h1, QPX_x1)
      QPX_x2 = VEC_MADD(QPX_q2, QPX_h1, QPX_x2)
      QPX_x3 = VEC_MADD(QPX_q3, QPX_h1, QPX_x3)
      QPX_x4 = VEC_MADD(QPX_q4, QPX_h1, QPX_x4)
      QPX_h2 = VEC_SPLATS(hh(i,2))
      QPX_y1 = VEC_MADD(QPX_q1, QPX_h2, QPX_y1)
      QPX_y2 = VEC_MADD(QPX_q2, QPX_h2, QPX_y2)
      QPX_y3 = VEC_MADD(QPX_q3, QPX_h2, QPX_y3)
      QPX_y4 = VEC_MADD(QPX_q4, QPX_h2, QPX_y4)

   enddo

   QPX_h1 = VEC_SPLATS(hh(nb,1))
   QPX_q1 = VEC_LD(0,q(1,nb+1))
   QPX_q2 = VEC_LD(0,q(5,nb+1))
   QPX_q3 = VEC_LD(0,q(9,nb+1))
   QPX_q4 = VEC_LD(0,q(13,nb+1))
   QPX_x1 = VEC_MADD(QPX_q1, QPX_h1, QPX_x1)
   QPX_x2 = VEC_MADD(QPX_q2, QPX_h1, QPX_x2)
   QPX_x3 = VEC_MADD(QPX_q3, QPX_h1, QPX_x3)
   QPX_x4 = VEC_MADD(QPX_q4, QPX_h1, QPX_x4)
   
   !--- multiply T matrix ---
   
   QPX_tau1 = VEC_SPLATS(-hh(1,1))
   QPX_x1 = VEC_MUL(QPX_x1, QPX_tau1)
   QPX_x2 = VEC_MUL(QPX_x2, QPX_tau1)
   QPX_x3 = VEC_MUL(QPX_x3, QPX_tau1)
   QPX_x4 = VEC_MUL(QPX_x4, QPX_tau1)
   QPX_tau2 = VEC_SPLATS(-hh(1,2))
   QPX_s = VEC_SPLATS(-hh(1,2)*s)
   QPX_y1 = VEC_MUL(QPX_y1, QPX_tau2)
   QPX_y2 = VEC_MUL(QPX_y2, QPX_tau2)
   QPX_y3 = VEC_MUL(QPX_y3, QPX_tau2)
   QPX_y4 = VEC_MUL(QPX_y4, QPX_tau2)
   QPX_y1 = VEC_MADD(QPX_x1, QPX_s, QPX_y1)
   QPX_y2 = VEC_MADD(QPX_x2, QPX_s, QPX_y2)
   QPX_y3 = VEC_MADD(QPX_x3, QPX_s, QPX_y3)
   QPX_y4 = VEC_MADD(QPX_x4, QPX_s, QPX_y4)

   !--- rank-2 update of q ---

   QPX_q1 = VEC_LD(0,q(1,1))
   QPX_q2 = VEC_LD(0,q(5,1))
   QPX_q3 = VEC_LD(0,q(9,1))
   QPX_q4 = VEC_LD(0,q(13,1))
   QPX_q1 = VEC_ADD(QPX_q1, QPX_y1)
   QPX_q2 = VEC_ADD(QPX_q2, QPX_y2)
   QPX_q3 = VEC_ADD(QPX_q3, QPX_y3)
   QPX_q4 = VEC_ADD(QPX_q4, QPX_y4)
   call VEC_ST(QPX_q1, 0, q(1,1))
   call VEC_ST(QPX_q2, 0, q(5,1))
   call VEC_ST(QPX_q3, 0, q(9,1))
   call VEC_ST(QPX_q4, 0, q(13,1))
   
   QPX_h2 = VEC_SPLATS(hh(2,2))
   QPX_q1 = VEC_LD(0,q(1,2))
   QPX_q2 = VEC_LD(0,q(5,2))
   QPX_q3 = VEC_LD(0,q(9,2))
   QPX_q4 = VEC_LD(0,q(13,2))
   QPX_q1 = VEC_MADD(QPX_y1, QPX_h2, QPX_q1)
   QPX_q2 = VEC_MADD(QPX_y2, QPX_h2, QPX_q2)
   QPX_q3 = VEC_MADD(QPX_y3, QPX_h2, QPX_q3)
   QPX_q4 = VEC_MADD(QPX_y4, QPX_h2, QPX_q4)
   QPX_q1 = VEC_ADD(QPX_q1, QPX_x1)
   QPX_q2 = VEC_ADD(QPX_q2, QPX_x2)
   QPX_q3 = VEC_ADD(QPX_q3, QPX_x3)
   QPX_q4 = VEC_ADD(QPX_q4, QPX_x4)
   call VEC_ST(QPX_q1, 0, q(1,2))
   call VEC_ST(QPX_q2, 0, q(5,2))
   call VEC_ST(QPX_q3, 0, q(9,2))
   call VEC_ST(QPX_q4, 0, q(13,2))
   
   do i=3,nb,1

      QPX_q1 = VEC_LD(0,q(1,i))
      QPX_q2 = VEC_LD(0,q(5,i))
      QPX_q3 = VEC_LD(0,q(9,i))
      QPX_q4 = VEC_LD(0,q(13,i))
      QPX_h1 = VEC_SPLATS(hh(i-1,1))
      QPX_q1 = VEC_MADD(QPX_x1, QPX_h1, QPX_q1)
      QPX_q2 = VEC_MADD(QPX_x2, QPX_h1, QPX_q2)
      QPX_q3 = VEC_MADD(QPX_x3, QPX_h1, QPX_q3)
      QPX_q4 = VEC_MADD(QPX_x4, QPX_h1, QPX_q4)
      QPX_h2 = VEC_SPLATS(hh(i,2))
      QPX_q1 = VEC_MADD(QPX_y1, QPX_h2, QPX_q1)
      QPX_q2 = VEC_MADD(QPX_y2, QPX_h2, QPX_q2)
      QPX_q3 = VEC_MADD(QPX_y3, QPX_h2, QPX_q3)
      QPX_q4 = VEC_MADD(QPX_y4, QPX_h2, QPX_q4)
      
      call VEC_ST(QPX_q1, 0, q(1,i))
      call VEC_ST(QPX_q2, 0, q(5,i))
      call VEC_ST(QPX_q3, 0, q(9,i))
      call VEC_ST(QPX_q4, 0, q(13,i))

   enddo

   QPX_h1 = VEC_SPLATS(hh(nb,1))
   QPX_q1 = VEC_LD(0,q(1,nb+1))
   QPX_q2 = VEC_LD(0,q(5,nb+1))
   QPX_q3 = VEC_LD(0,q(9,nb+1))
   QPX_q4 = VEC_LD(0,q(13,nb+1))
   QPX_q1 = VEC_MADD(QPX_x1, QPX_h1, QPX_q1)
   QPX_q2 = VEC_MADD(QPX_x2, QPX_h1, QPX_q2)
   QPX_q3 = VEC_MADD(QPX_x3, QPX_h1, QPX_q3)
   QPX_q4 = VEC_MADD(QPX_x4, QPX_h1, QPX_q4)
   call VEC_ST(QPX_q1, 0, q(1,nb+1))
   call VEC_ST(QPX_q2, 0, q(5,nb+1))
   call VEC_ST(QPX_q3, 0, q(9,nb+1))
   call VEC_ST(QPX_q4, 0, q(13,nb+1))

end subroutine

! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_kernel_8_bgq(q, hh, nb, ldq, ldh, s)

   implicit none

   include 'mpif.h'
   
   integer, intent(in) :: nb, ldq, ldh
   
   real*8, intent(inout)::q(ldq,*)
   real*8, intent(in)::hh(ldh,*), s
   
   VECTOR(REAL(8))::QPX_x1, QPX_x2, QPX_y1, QPX_y2
   VECTOR(REAL(8))::QPX_q1, QPX_q2
   VECTOR(REAL(8))::QPX_h1, QPX_h2, QPX_tau1, QPX_tau2, QPX_s
   integer i

   call alignx(32,q)

   !--- multiply Householder vectors with matrix q ---

   QPX_x1 = VEC_LD(0,q(1,2))
   QPX_x2 = VEC_LD(0,q(5,2))
   
   QPX_h2 = VEC_SPLATS(hh(2,2))
   QPX_q1 = VEC_LD(0,q(1,1))
   QPX_q2 = VEC_LD(0,q(5,1))
   QPX_y1 = VEC_MADD(QPX_x1, QPX_h2, QPX_q1)
   QPX_y2 = VEC_MADD(QPX_x2, QPX_h2, QPX_q2)

   do i=3,nb,1

      QPX_q1 = VEC_LD(0,q(1,i))
      QPX_q2 = VEC_LD(0,q(5,i))
      QPX_h1 = VEC_SPLATS(hh(i-1,1))
      QPX_x1 = VEC_MADD(QPX_q1, QPX_h1, QPX_x1)
      QPX_x2 = VEC_MADD(QPX_q2, QPX_h1, QPX_x2)
      QPX_h2 = VEC_SPLATS(hh(i,2))
      QPX_y1 = VEC_MADD(QPX_q1, QPX_h2, QPX_y1)
      QPX_y2 = VEC_MADD(QPX_q2, QPX_h2, QPX_y2)

   enddo

   QPX_h1 = VEC_SPLATS(hh(nb,1))
   QPX_q1 = VEC_LD(0,q(1,nb+1))
   QPX_q2 = VEC_LD(0,q(5,nb+1))
   QPX_x1 = VEC_MADD(QPX_q1, QPX_h1, QPX_x1)
   QPX_x2 = VEC_MADD(QPX_q2, QPX_h1, QPX_x2)
   
   !--- multiply T matrix ---
   
   QPX_tau1 = VEC_SPLATS(-hh(1,1))
   QPX_x1 = VEC_MUL(QPX_x1, QPX_tau1)
   QPX_x2 = VEC_MUL(QPX_x2, QPX_tau1)
   QPX_tau2 = VEC_SPLATS(-hh(1,2))
   QPX_s = VEC_SPLATS(-hh(1,2)*s)
   QPX_y1 = VEC_MUL(QPX_y1, QPX_tau2)
   QPX_y2 = VEC_MUL(QPX_y2, QPX_tau2)
   QPX_y1 = VEC_MADD(QPX_x1, QPX_s, QPX_y1)
   QPX_y2 = VEC_MADD(QPX_x2, QPX_s, QPX_y2)

   !--- rank-2 update of q ---

   QPX_q1 = VEC_LD(0,q(1,1))
   QPX_q2 = VEC_LD(0,q(5,1))
   QPX_q1 = VEC_ADD(QPX_q1, QPX_y1)
   QPX_q2 = VEC_ADD(QPX_q2, QPX_y2)
   call VEC_ST(QPX_q1, 0, q(1,1))
   call VEC_ST(QPX_q2, 0, q(5,1))
   
   QPX_h2 = VEC_SPLATS(hh(2,2))
   QPX_q1 = VEC_LD(0,q(1,2))
   QPX_q2 = VEC_LD(0,q(5,2))
   QPX_q1 = VEC_MADD(QPX_y1, QPX_h2, QPX_q1)
   QPX_q2 = VEC_MADD(QPX_y2, QPX_h2, QPX_q2)
   QPX_q1 = VEC_ADD(QPX_q1, QPX_x1)
   QPX_q2 = VEC_ADD(QPX_q2, QPX_x2)
   call VEC_ST(QPX_q1, 0, q(1,2))
   call VEC_ST(QPX_q2, 0, q(5,2))
   
   do i=3,nb,1

      QPX_q1 = VEC_LD(0,q(1,i))
      QPX_q2 = VEC_LD(0,q(5,i))
      QPX_h1 = VEC_SPLATS(hh(i-1,1))
      QPX_q1 = VEC_MADD(QPX_x1, QPX_h1, QPX_q1)
      QPX_q2 = VEC_MADD(QPX_x2, QPX_h1, QPX_q2)
      QPX_h2 = VEC_SPLATS(hh(i,2))
      QPX_q1 = VEC_MADD(QPX_y1, QPX_h2, QPX_q1)
      QPX_q2 = VEC_MADD(QPX_y2, QPX_h2, QPX_q2)
      
      call VEC_ST(QPX_q1, 0, q(1,i))
      call VEC_ST(QPX_q2, 0, q(5,i))

   enddo

   QPX_h1 = VEC_SPLATS(hh(nb,1))
   QPX_q1 = VEC_LD(0,q(1,nb+1))
   QPX_q2 = VEC_LD(0,q(5,nb+1))
   QPX_q1 = VEC_MADD(QPX_x1, QPX_h1, QPX_q1)
   QPX_q2 = VEC_MADD(QPX_x2, QPX_h1, QPX_q2)
   call VEC_ST(QPX_q1, 0, q(1,nb+1))
   call VEC_ST(QPX_q2, 0, q(5,nb+1))

end subroutine

! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_kernel_4_bgq(q, hh, nb, ldq, ldh, s)

   implicit none

   include 'mpif.h'
   
   integer, intent(in) :: nb, ldq, ldh
   
   real*8, intent(inout)::q(ldq,*)
   real*8, intent(in)::hh(ldh,*), s
   
   VECTOR(REAL(8))::QPX_x1, QPX_y1
   VECTOR(REAL(8))::QPX_q1
   VECTOR(REAL(8))::QPX_h1, QPX_h2, QPX_tau1, QPX_tau2, QPX_s
   integer i

   call alignx(32,q)

   !--- multiply Householder vectors with matrix q ---

   QPX_x1 = VEC_LD(0,q(1,2))
   
   QPX_h2 = VEC_SPLATS(hh(2,2))
   QPX_q1 = VEC_LD(0,q(1,1))
   QPX_y1 = VEC_MADD(QPX_x1, QPX_h2, QPX_q1)
   
   do i=3,nb,1
   
      QPX_q1 = VEC_LD(0,q(1,i))
      QPX_h1 = VEC_SPLATS(hh(i-1,1))
      QPX_x1 = VEC_MADD(QPX_q1, QPX_h1, QPX_x1)
      QPX_h2 = VEC_SPLATS(hh(i,2))
      QPX_y1 = VEC_MADD(QPX_q1, QPX_h2, QPX_y1)
   
   enddo
   
   QPX_h1 = VEC_SPLATS(hh(nb,1))
   QPX_q1 = VEC_LD(0,q(1,nb+1))
   QPX_x1 = VEC_MADD(QPX_q1, QPX_h1, QPX_x1)
   
   !--- multiply T matrix ---
   
   QPX_tau1 = VEC_SPLATS(-hh(1,1))
   QPX_x1 = VEC_MUL(QPX_x1, QPX_tau1)
   QPX_tau2 = VEC_SPLATS(-hh(1,2))
   QPX_s = VEC_SPLATS(-hh(1,2)*s)
   QPX_y1 = VEC_MUL(QPX_y1, QPX_tau2)
   QPX_y1 = VEC_MADD(QPX_x1, QPX_s, QPX_y1)

   !--- rank-2 update of q ---

   QPX_q1 = VEC_LD(0,q(1,1))
   QPX_q1 = VEC_ADD(QPX_q1, QPX_y1)
   call VEC_ST(QPX_q1, 0, q(1,1))
   
   QPX_h2 = VEC_SPLATS(hh(2,2))
   QPX_q1 = VEC_LD(0,q(1,2))
   QPX_q1 = VEC_MADD(QPX_y1, QPX_h2, QPX_q1)
   QPX_q1 = VEC_ADD(QPX_q1, QPX_x1)
   call VEC_ST(QPX_q1, 0, q(1,2))
   
   do i=3,nb,1

      QPX_q1 = VEC_LD(0,q(1,i))
      QPX_h1 = VEC_SPLATS(hh(i-1,1))
      QPX_q1 = VEC_MADD(QPX_x1, QPX_h1, QPX_q1)
      QPX_h2 = VEC_SPLATS(hh(i,2))
      QPX_q1 = VEC_MADD(QPX_y1, QPX_h2, QPX_q1)
      
      call VEC_ST(QPX_q1, 0, q(1,i))

   enddo

   QPX_h1 = VEC_SPLATS(hh(nb,1))
   QPX_q1 = VEC_LD(0,q(1,nb+1))
   QPX_q1 = VEC_MADD(QPX_x1, QPX_h1, QPX_q1)
   call VEC_ST(QPX_q1, 0, q(1,nb+1))

end subroutine

! --------------------------------------------------------------------------------------------------

! BL: The following kernels are taken from the generic elpa library. 
!     They are not optimized for BlueGene/Q
!     TODO

subroutine single_hh_trafo_complex(q, hh, nb, nq, ldq)

   implicit none

   integer, intent(in) :: nb, nq, ldq
   complex*16, intent(inout) :: q(ldq,*)
   complex*16, intent(in) :: hh(*)

   integer i

   ! Safety only:

   if(mod(ldq,4) /= 0) STOP 'double_hh_trafo: ldq not divisible by 4!'

   ! Do the Householder transformations

   ! Always a multiple of 4 Q-rows is transformed, even if nq is smaller

   do i=1,nq-8,12
      call hh_trafo_complex_kernel_12(q(i,1),hh, nb, ldq)
   enddo

   ! i > nq-8 now, i.e. at most 8 rows remain

   if(nq-i+1 > 4) then
      call hh_trafo_complex_kernel_8(q(i,1),hh, nb, ldq)
   else if(nq-i+1 > 0) then
      call hh_trafo_complex_kernel_4(q(i,1),hh, nb, ldq)
   endif

end

! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_complex_kernel_12(q, hh, nb, ldq)

   implicit none

   integer, intent(in) :: nb, ldq
   complex*16, intent(inout) :: q(ldq,*)
   complex*16, intent(in) :: hh(*)

   complex*16 x1, x2, x3, x4, x5, x6, x7, x8, x9, xa, xb, xc
   complex*16 h1, tau1
   integer i


   x1 = q(1,1)
   x2 = q(2,1)
   x3 = q(3,1)
   x4 = q(4,1)
   x5 = q(5,1)
   x6 = q(6,1)
   x7 = q(7,1)
   x8 = q(8,1)
   x9 = q(9,1)
   xa = q(10,1)
   xb = q(11,1)
   xc = q(12,1)

!DEC$ VECTOR ALIGNED
   do i=2,nb
      h1 = conjg(hh(i))
      x1 = x1 + q(1,i)*h1
      x2 = x2 + q(2,i)*h1
      x3 = x3 + q(3,i)*h1
      x4 = x4 + q(4,i)*h1
      x5 = x5 + q(5,i)*h1
      x6 = x6 + q(6,i)*h1
      x7 = x7 + q(7,i)*h1
      x8 = x8 + q(8,i)*h1
      x9 = x9 + q(9,i)*h1
      xa = xa + q(10,i)*h1
      xb = xb + q(11,i)*h1
      xc = xc + q(12,i)*h1
   enddo

   tau1 = hh(1)

   h1 = -tau1
   x1 = x1*h1
   x2 = x2*h1
   x3 = x3*h1
   x4 = x4*h1
   x5 = x5*h1
   x6 = x6*h1
   x7 = x7*h1
   x8 = x8*h1
   x9 = x9*h1
   xa = xa*h1
   xb = xb*h1
   xc = xc*h1

   q(1,1) = q(1,1) + x1
   q(2,1) = q(2,1) + x2
   q(3,1) = q(3,1) + x3
   q(4,1) = q(4,1) + x4
   q(5,1) = q(5,1) + x5
   q(6,1) = q(6,1) + x6
   q(7,1) = q(7,1) + x7
   q(8,1) = q(8,1) + x8
   q(9,1) = q(9,1) + x9
   q(10,1) = q(10,1) + xa
   q(11,1) = q(11,1) + xb
   q(12,1) = q(12,1) + xc

!DEC$ VECTOR ALIGNED
   do i=2,nb
      h1 = hh(i)
      q(1,i) = q(1,i) + x1*h1
      q(2,i) = q(2,i) + x2*h1
      q(3,i) = q(3,i) + x3*h1
      q(4,i) = q(4,i) + x4*h1
      q(5,i) = q(5,i) + x5*h1
      q(6,i) = q(6,i) + x6*h1
      q(7,i) = q(7,i) + x7*h1
      q(8,i) = q(8,i) + x8*h1
      q(9,i) = q(9,i) + x9*h1
      q(10,i) = q(10,i) + xa*h1
      q(11,i) = q(11,i) + xb*h1
      q(12,i) = q(12,i) + xc*h1
   enddo

end

! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_complex_kernel_8(q, hh, nb, ldq)

   implicit none

   integer, intent(in) :: nb, ldq
   complex*16, intent(inout) :: q(ldq,*)
   complex*16, intent(in) :: hh(*)

   complex*16 x1, x2, x3, x4, x5, x6, x7, x8
   complex*16 h1, tau1
   integer i


   x1 = q(1,1)
   x2 = q(2,1)
   x3 = q(3,1)
   x4 = q(4,1)
   x5 = q(5,1)
   x6 = q(6,1)
   x7 = q(7,1)
   x8 = q(8,1)

!DEC$ VECTOR ALIGNED
   do i=2,nb
      h1 = conjg(hh(i))
      x1 = x1 + q(1,i)*h1
      x2 = x2 + q(2,i)*h1
      x3 = x3 + q(3,i)*h1
      x4 = x4 + q(4,i)*h1
      x5 = x5 + q(5,i)*h1
      x6 = x6 + q(6,i)*h1
      x7 = x7 + q(7,i)*h1
      x8 = x8 + q(8,i)*h1
   enddo

   tau1 = hh(1)

   h1 = -tau1
   x1 = x1*h1
   x2 = x2*h1
   x3 = x3*h1
   x4 = x4*h1
   x5 = x5*h1
   x6 = x6*h1
   x7 = x7*h1
   x8 = x8*h1

   q(1,1) = q(1,1) + x1
   q(2,1) = q(2,1) + x2
   q(3,1) = q(3,1) + x3
   q(4,1) = q(4,1) + x4
   q(5,1) = q(5,1) + x5
   q(6,1) = q(6,1) + x6
   q(7,1) = q(7,1) + x7
   q(8,1) = q(8,1) + x8

!DEC$ VECTOR ALIGNED
   do i=2,nb
      h1 = hh(i)
      q(1,i) = q(1,i) + x1*h1
      q(2,i) = q(2,i) + x2*h1
      q(3,i) = q(3,i) + x3*h1
      q(4,i) = q(4,i) + x4*h1
      q(5,i) = q(5,i) + x5*h1
      q(6,i) = q(6,i) + x6*h1
      q(7,i) = q(7,i) + x7*h1
      q(8,i) = q(8,i) + x8*h1
   enddo

end

! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_complex_kernel_4(q, hh, nb, ldq)

   implicit none

   integer, intent(in) :: nb, ldq
   complex*16, intent(inout) :: q(ldq,*)
   complex*16, intent(in) :: hh(*)

   complex*16 x1, x2, x3, x4
   complex*16 h1, tau1
   integer i


   x1 = q(1,1)
   x2 = q(2,1)
   x3 = q(3,1)
   x4 = q(4,1)

!DEC$ VECTOR ALIGNED
   do i=2,nb
      h1 = conjg(hh(i))
      x1 = x1 + q(1,i)*h1
      x2 = x2 + q(2,i)*h1
      x3 = x3 + q(3,i)*h1
      x4 = x4 + q(4,i)*h1
   enddo

   tau1 = hh(1)

   h1 = -tau1
   x1 = x1*h1
   x2 = x2*h1
   x3 = x3*h1
   x4 = x4*h1

   q(1,1) = q(1,1) + x1
   q(2,1) = q(2,1) + x2
   q(3,1) = q(3,1) + x3
   q(4,1) = q(4,1) + x4

!DEC$ VECTOR ALIGNED
   do i=2,nb
      h1 = hh(i)
      q(1,i) = q(1,i) + x1*h1
      q(2,i) = q(2,i) + x2*h1
      q(3,i) = q(3,i) + x3*h1
      q(4,i) = q(4,i) + x4*h1
   enddo

end

! --------------------------------------------------------------------------------------------------
