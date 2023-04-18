/* 
  This file was generated automatically with /nfs/data-012/marques/software/source/libxc/svn/scripts/maple2c.pl.
  Do not edit this file directly as it can be overwritten!!

  Maple version     : Maple 2016 (X86 64 LINUX)
  Maple source      : ../maple/gga_x_airy.mpl
  Type of functional: work_gga_x
*/

void xc_gga_x_airy_enhance
  (const xc_func_type *p,  xc_gga_work_x_t *r)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t10, t11, t14, t16, t18, t19, t21, t22;
  double t24, t28, t29, t33, t37, t41, t43, t44;
  double t45, t46, t48, t51, t53, t54, t55, t58;
  double t62, t63, t67, t71, t75, t77, t81, t82;
  double t83, t85, t88, t90, t93, t95, t98, t99;
  double t103, t104, t108, t109, t112, t115, t128, t131;
  double t134, t140;


  t1 = M_CBRT6;
  t2 = t1 * t1;
  t3 = 0.31415926535897932385e1 * 0.31415926535897932385e1;
  t4 = cbrt(t3);
  t5 = 0.1e1 / t4;
  t6 = t2 * t5;
  t7 = t6 * r->x;
  t8 = pow(t7, 0.2626712e1);
  t10 = 0.1e1 + 0.13471619689594796103e-3 * t8;
  t11 = pow(t10, -0.657946e0);
  t14 = pow(t7, 0.3217063e1);
  t16 = pow(t7, 0.3223476e1);
  t18 = 0.1e1 - 0.45212413010769857073e-1 * t14 + 0.45402221956620378581e-1 * t16;
  t19 = pow(t7, 0.3473804e1);
  t21 = 0.1e1 + 0.47702180224903349918e-3 * t19;
  t22 = 0.1e1 / t21;
  r->f = 0.60146019220211109872e-4 * t8 * t11 + t18 * t22;

  if(r->order < 1) return;

  t24 = pow(t7, 0.1626712e1);
  t28 = pow(t7, 0.4253424e1);
  t29 = pow(t10, -0.1657946e1);
  t33 = pow(t7, 0.2217063e1);
  t37 = pow(t7, 0.2223476e1);
  t41 = -0.14545118103766630870e0 * t33 * t2 * t5 + 0.14635297282383883147e0 * t37 * t2 * t5;
  t43 = t21 * t21;
  t44 = 0.1e1 / t43;
  t45 = t18 * t44;
  t46 = pow(t7, 0.2473804e1);
  t48 = t46 * t2 * t5;
  r->dfdx = 0.15798627043795916483e-3 * t24 * t11 * t6 - 0.14003268362272376394e-7 * t28 * t29 * t6 + t41 * t22 - 0.16570802447399015656e-2 * t45 * t48;

  if(r->order < 2) return;

  t51 = pow(t7, 0.626712e0);
  t53 = t4 * t4;
  t54 = 0.1e1 / t53;
  t55 = t1 * t54;
  t58 = pow(t7, 0.3253424e1);
  t62 = pow(t7, 0.5880136e1);
  t63 = pow(t10, -0.2657946e1);
  t67 = pow(t7, 0.1217063e1);
  t71 = pow(t7, 0.1223476e1);
  t75 = -0.19348465907094694762e1 * t67 * t1 * t54 + 0.19524739356147472178e1 * t71 * t1 * t54;
  t77 = t41 * t44;
  t81 = 0.1e1 / t43 / t21;
  t82 = t18 * t81;
  t83 = pow(t7, 0.4947608e1);
  t85 = t83 * t1 * t54;
  t88 = pow(t7, 0.1473804e1);
  t90 = t88 * t1 * t54;
  r->d2fdx2 = 0.15419889717400405736e-2 * t51 * t11 * t55 - 0.57806634466158731181e-6 * t58 * t29 * t55 + 0.49292780404469177134e-10 * t62 * t63 * t55 + t75 * t22 - 0.33141604894798031312e-2 * t77 * t48 + 0.32950979250087024843e-4 * t82 * t85 - 0.24595750426551284716e-1 * t45 * t90;

  if(r->order < 3) return;

  t93 = pow(t7, -0.373288e0);
  t95 = 0.1e1 / t3;
  t98 = pow(t7, 0.2253424e1);
  t99 = t98 * t29;
  t103 = pow(t7, 0.4880136e1);
  t104 = t103 * t63;
  t108 = pow(t7, 0.7506848e1);
  t109 = pow(t10, -0.3657946e1);
  t112 = pow(t7, 0.217063e0);
  t115 = pow(t7, 0.223476e0);
  t128 = t43 * t43;
  t131 = pow(t7, 0.7421412e1);
  t134 = pow(t7, 0.3947608e1);
  t140 = pow(t7, 0.473804e0);
  r->d3fdx3 = 0.57982979547428658478e-2 * t93 * t11 * t95 - 0.21825052433571428779e-6 * t99 - 0.11284169515885680230e-4 * t99 * t95 + 0.20617302531200757315e-9 * t104 + 0.17390895155784826161e-8 * t104 * t95 - 0.28184688941798707944e-13 * t108 * t109 + (-0.14128981177371834295e2 * t112 * t95 + 0.14332830005101130802e2 * t115 * t95) * t22 - 0.49712407342197046968e-2 * t75 * t44 * t48 + 0.98852937750261074529e-4 * t41 * t81 * t85 - 0.73787251279653854147e-1 * t77 * t90 - 0.99582866925677770903e-7 * t18 / t128 * t131 + 0.97817117127338738886e-3 * t82 * t134 * t95 + 0.49554730439115758451e-4 * t82 * t134 - 0.21749589216991793772e0 * t45 * t140 * t95;

  if(r->order < 4) return;


}

#define maple2c_order 3
#define maple2c_func  xc_gga_x_airy_enhance
