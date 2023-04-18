/*
*******************************************************************************
*                                                                             *
*                                PLUMED                                       *
*   A Portable Plugin for Free Energy Calculations with Molecular Dynamics    *
*                              VERSION 1.1                                    *
*                                                                             *
*******************************************************************************
*
*  
*  Copyright (c) 2009 The PLUMED team.
*  See http://merlino.mi.infn.it/plumed for more information. 
*
*  This file is part of PLUMED.
*
*  PLUMED is free software: you can redistribute it and/or modify
*  it under the terms of the GNU Lesser General Public License as 
*  published by the Free Software Foundation, either version 3 of 
*  the License, or (at your option) any later version.
*
*  PLUMED is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General
*  Public License along with PLUMED.  
*  If not, see <http://www.gnu.org/licenses/>.
*
*  For more info, see:  http://merlino.mi.infn.it/plumed
*  or subscribe to plumed-users@googlegroups.com
*
*/
#include "metadyn.h"

void PREFIX rmsdtor_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int firstAtom, secondAtom, thirdAtom, fourthAtom, i, index;
  rvec rij, sij, tij, vv;
  real mod_rij, mod_sij, mod_tij;
  real t1, t2, t3, t4, t5, t6, t7, t8, t9;
  real caa, cbb, ccc, cab, cac, cbc, dab, dbc;
  real ac, cost, ddd, Vac, dchi;
  real sign0, s1, eps;

  colvar.ss0[i_c] = 0.;
  index = 0;

  // this CV does a cycle over the torsion angle cv
  for(i=0;i<colvar.type[i_c];i++){
    firstAtom  = colvar.cvatoms[i_c][4*i];
    secondAtom = colvar.cvatoms[i_c][1+4*i];
    thirdAtom  = colvar.cvatoms[i_c][2+4*i];
    fourthAtom = colvar.cvatoms[i_c][3+4*i];

    minimal_image(mtd_data->pos[firstAtom],  mtd_data->pos[secondAtom], &mod_rij, rij);
    minimal_image(mtd_data->pos[secondAtom], mtd_data->pos[thirdAtom],  &mod_sij, sij);
    minimal_image(mtd_data->pos[thirdAtom],  mtd_data->pos[fourthAtom], &mod_tij, tij);

    t1 = rij[0];
    t2 = rij[1];
    t3 = rij[2];
    t4 = sij[0];
    t5 = sij[1];
    t6 = sij[2];
    t7 = tij[0];
    t8 = tij[1];
    t9 = tij[2];

    caa = t1*t1+t2*t2+t3*t3;
    cbb = t4*t4+t5*t5+t6*t6;
    ccc = t7*t7+t8*t8+t9*t9;

    cab = t1*t4+t2*t5+t3*t6;
    cac = t1*t7+t2*t8+t3*t9;
    cbc = t4*t7+t5*t8+t6*t9;

    dab = caa*cbb-cab*cab;
    dbc = cbb*ccc-cbc*cbc;

    ac = -(cab*cbc-cac*cbb)/sqrt(dab*dbc);
    if(fabs(ac)>0.999999999999) ac = 0.999999999999;

    Vac = acos(ac);

    vv[0] = t4*(t3*(-t5*t7 + t4*t8) + t2*(t6*t7 - t4*t9) + t1*(-t6*t8 + t5*t9));
    vv[1] = t5*(t3*(-t5*t7 + t4*t8) + t2*(t6*t7 - t4*t9) + t1*(-t6*t8 + t5*t9));
    vv[2] = t6*(t3*(-t5*t7 + t4*t8) + t2*(t6*t7 - t4*t9) + t1*(-t6*t8 + t5*t9));
    s1 = -vv[0]*t4-vv[1]*t5-vv[2]*t6;

    sign0 = +1.;
    if(s1<0.){
      sign0 = -1.;
      Vac = -Vac;
    }

    cost = sign0/sqrt(1.0-ac*ac);
    ddd = cost/sqrt(dab*dbc);

    dchi = Vac-colvar.map0[i_c][i];
    if(dchi>M_PI) dchi -= 2.*M_PI;
    if(dchi<-M_PI) dchi += 2.*M_PI;

    colvar.ss0[i_c] += 1.0/((real) colvar.type[i_c])*dchi*dchi;

    eps = 1.e-12;
    if(dab>eps && dbc>eps) {
      colvar.myder[i_c][index][0] = +2.0/(real) colvar.type[i_c]*dchi*ddd*(cbc*t4-cbb*t7-(cab*cbc-cac*cbb)/dab*(cbb*t1-cab*t4));
      colvar.myder[i_c][index][1] = +2.0/(real) colvar.type[i_c]*dchi*ddd*(cbc*t5-cbb*t8-(cab*cbc-cac*cbb)/dab*(cbb*t2-cab*t5));
      colvar.myder[i_c][index][2] = +2.0/(real) colvar.type[i_c]*dchi*ddd*(cbc*t6-cbb*t9-(cab*cbc-cac*cbb)/dab*(cbb*t3-cab*t6));

      colvar.myder[i_c][index+1][0] = +2.0/(real) colvar.type[i_c]*dchi*ddd*(cbc*t1-cbc*t4+cab*t7+cbb*t7-2*cac*t4-
             (cab*cbc-cac*cbb)/dbc*(ccc*t4-cbc*t7)-
             (cab*cbc-cac*cbb)/dab*(caa*t4-cbb*t1-cab*t1+cab*t4));
      colvar.myder[i_c][index+1][1] = +2.0/(real) colvar.type[i_c]*dchi*ddd*(cbc*t2-cbc*t5+cab*t8+cbb*t8-2*cac*t5-
             (cab*cbc-cac*cbb)/dbc*(ccc*t5-cbc*t8)-
             (cab*cbc-cac*cbb)/dab*(caa*t5-cbb*t2-cab*t2+cab*t5));
      colvar.myder[i_c][index+1][2] = +2.0/(real) colvar.type[i_c]*dchi*ddd*(cbc*t3-cbc*t6+cab*t9+cbb*t9-2*cac*t6-
             (cab*cbc-cac*cbb)/dbc*(ccc*t6-cbc*t9)-
             (cab*cbc-cac*cbb)/dab*(caa*t6-cbb*t3-cab*t3+cab*t6));

      colvar.myder[i_c][index+2][0] = +2.0/(real) colvar.type[i_c]*dchi*ddd*(-cbc*t1+cab*t4-cab*t7-cbb*t1+2.*cac*t4-
             (cab*cbc-cac*cbb)/dbc*(cbb*t7-ccc*t4-cbc*t4+cbc*t7)-
             (cab*cbc-cac*cbb)/dab*(-caa*t4+cab*t1));
      colvar.myder[i_c][index+2][1] = +2.0/(real) colvar.type[i_c]*dchi*ddd*(-cbc*t2+cab*t5-cab*t8-cbb*t2+2.*cac*t5-
             (cab*cbc-cac*cbb)/dbc*(cbb*t8-ccc*t5-cbc*t5+cbc*t8)-
             (cab*cbc-cac*cbb)/dab*(-caa*t5+cab*t2));
      colvar.myder[i_c][index+2][2] = +2.0/(real) colvar.type[i_c]*dchi*ddd*(-cbc*t3+cab*t6-cab*t9-cbb*t3+2.*cac*t6-
             (cab*cbc-cac*cbb)/dbc*(cbb*t9-ccc*t6-cbc*t6+cbc*t9)-
             (cab*cbc-cac*cbb)/dab*(-caa*t6+cab*t3));

      colvar.myder[i_c][index+3][0] = +2.0/(real) colvar.type[i_c]*dchi*ddd*(-cab*t4+cbb*t1-(cab*cbc-cac*cbb)/dbc*(-cbb*t7+cbc*t4));
      colvar.myder[i_c][index+3][1] = +2.0/(real) colvar.type[i_c]*dchi*ddd*(-cab*t5+cbb*t2-(cab*cbc-cac*cbb)/dbc*(-cbb*t8+cbc*t5));
      colvar.myder[i_c][index+3][2] = +2.0/(real) colvar.type[i_c]*dchi*ddd*(-cab*t6+cbb*t3-(cab*cbc-cac*cbb)/dbc*(-cbb*t9+cbc*t6));
      index = index+4;
    }
  }
}

//---------------------------------------------------------------------------------------------------------

int PREFIX read_rmsdtor(char **word, int count, t_plumed_input *input,int *iline,FILE *fplog)
{
  int iat, atnr1, atnr2, atnr3, atnr4, iw, help;
  double delta = 0.0;
  double psi;
  char string[400];
  help=0;

  iw = seek_word(word,"NDIH");
  if(iw>=0) {sscanf(word[iw+1],"%i", &colvar.type[count]); } else {fprintf(fplog,"|- NEEDED NDIH KEYWORD FOR RMSDTOR\n"); help=1;}
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

  if(help){
          fprintf(fplog, "\n-RMSDTOR CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "RMSDTOR NDIH 2 SIGMA 0.1\n");
          fprintf(fplog,"168 170 172 188 0.4\n");
          fprintf(fplog,"178 180 182 188 1.4\n");
          fprintf(stderr, "PluMed dead with errors: check log file  \n");
          EXIT(); 
  }

  colvar.type_s[count]   = 14;
  colvar.natoms[count]   = 4*colvar.type[count]; 

  snew(colvar.cvatoms[count], colvar.natoms[count]);
  snew(colvar.map0[count], colvar.type[count]);
  snew(colvar.myder[count], colvar.natoms[count]);

  fprintf(fplog, "\n%1i-RMSDTOR: %i DIHEDRAL; ", count+1, colvar.type[count]);
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 


  for(iat=0;iat<colvar.type[count];iat++){
    (*iline)++;
    sscanf(input->words[*iline][0],"%i",&atnr1);
    sscanf(input->words[*iline][1],"%i",&atnr2);
    sscanf(input->words[*iline][2],"%i",&atnr3);
    sscanf(input->words[*iline][3],"%i",&atnr4);
    sscanf(input->words[*iline][4],"%lf",&psi);
    atnr1--;
    atnr2--;
    atnr3--;
    atnr4--;
    psi *= M_PI/180.;
    colvar.cvatoms[count][4*iat] = atnr1;
    colvar.cvatoms[count][4*iat+1] = atnr2;
    colvar.cvatoms[count][4*iat+2] = atnr3;
    colvar.cvatoms[count][4*iat+3] = atnr4;
    colvar.map0[count][iat] = (real) psi;
    fprintf(fplog, "|--DIH %i, ATOMS: %i %i %i %i  R_0: %f\n", iat+1, atnr1+1, atnr2+1, atnr3+1, atnr4+1, psi*180./M_PI);
  }
  fprintf(fplog,"\n");
  return colvar.type[count];
}
