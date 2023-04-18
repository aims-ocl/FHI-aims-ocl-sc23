      !
      ! Three quadrature rules for 13,    25,      and 27 degrees with octohedral symmetry.
      !                          (HX0078) (HX0248) (HX0288)
      ! taken from
      ! Heo, Sangwoo, and Yuan Xu. "Constructing fully symmetric cubature formulae for the sphere." Mathematics of computation 70.233 (2001): 269-279.
      ! 
      ! further optimized by local optimization to reach double precision accuracy.
      ! 
      ! These are used to replace LD0074, LD0230, and LD0266 since they have negative weights in the quadrature and result in some
      ! odd behavior in aims calculation.
      !
      ! Yi Yao 
      ! Nov. 2019
      !
      !
       SUBROUTINE HX0078(X,Y,Z,W,N)
       DOUBLE PRECISION X(  78)
       DOUBLE PRECISION Y(  78)
       DOUBLE PRECISION Z(  78)
       DOUBLE PRECISION W(  78)
       INTEGER N
       DOUBLE PRECISION A,B,V
       N=1
       V=0.01386659210450215D+0
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2866401467665042D+0
       V=0.01305093186259184D+0
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6599050016563905D+0
       V=0.013206423223080691D+0
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5394900987058646D+0
       V=0.011942663554868659D+0
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END

       SUBROUTINE HX0248(X,Y,Z,W,N)
       DOUBLE PRECISION X( 248)
       DOUBLE PRECISION Y( 248)
       DOUBLE PRECISION Z( 248)
       DOUBLE PRECISION W( 248)
       INTEGER N
       DOUBLE PRECISION A,B,V
       N=1
       V=0.004313243133051708D+0
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.11169169091881402D+0
       V=0.003986365505453804D+0
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3150671668212475D+0
       V=0.0036630315485030413D+0
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.45946201454087915D+0
       V=0.004204049921581658D+0
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6607534971560259D+0
       V=0.004269004376172789D+0
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.702154945166355D+0
       V=0.004203472415480687D+0
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5320202557315966D+0
       V=0.0041424831176478D+0
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5196950515085538D+0
       B=0.8223599116864062D+0
       V=0.004090305598510224D+0
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3293373852026761D+0
       B=0.9382009660264478D+0
       V=0.003789950436894593D+0
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END

       SUBROUTINE HX0288(X,Y,Z,W,N)
       DOUBLE PRECISION X( 288)
       DOUBLE PRECISION Y( 288)
       DOUBLE PRECISION Z( 288)
       DOUBLE PRECISION W( 288)
       INTEGER N
       DOUBLE PRECISION A,B,V
       N=1
       A=0.11033297862407622D+0
       V=0.0038938290774785506D+0
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3190754015518861D+0
       V=0.0036062862032021923D+0
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.453117779552791D+0
       V=0.0038085043591521775D+0
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.614431551406119D+0
       V=0.002421634084601314D+0
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.702545464356561D+0
       V=0.004077606558264675D+0
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5301185129074026D+0
       V=0.004062279726900243D+0
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6222834318053699D+0
       B=0.708196634679607D+0
       V=0.0022635166912379273D+0
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5172773855252042D+0
       B=0.826495160268811D+0
       V=0.0037829348342885D+0
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3260968774773577D+0
       B=0.9390544873860379D+0
       V=0.0038518118030072936D+0
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
