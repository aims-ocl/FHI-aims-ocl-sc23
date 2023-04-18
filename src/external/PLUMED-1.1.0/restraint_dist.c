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
void PREFIX dist_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int i, iat;
  real mod_rij, mass1, mass2;
  rvec rij, sum1, sum2;

  mass1 = mass2 = 0.;
  sum1[0] = sum1[1] = sum1[2] = 0.;
  sum2[0] = sum2[1] = sum2[2] = 0.;


  for(i=0;i<colvar.list[i_c][0];i++) {
    iat = colvar.cvatoms[i_c][i];
    sum1[0] += mtd_data->mass[iat]*mtd_data->pos[iat][0];
    sum1[1] += mtd_data->mass[iat]*mtd_data->pos[iat][1];
    sum1[2] += mtd_data->mass[iat]*mtd_data->pos[iat][2];
    mass1 += mtd_data->mass[iat];
  }
  for(i=colvar.list[i_c][0];i<colvar.natoms[i_c];i++) {
    iat = colvar.cvatoms[i_c][i];
    sum2[0] += mtd_data->mass[iat]*mtd_data->pos[iat][0];
    sum2[1] += mtd_data->mass[iat]*mtd_data->pos[iat][1];
    sum2[2] += mtd_data->mass[iat]*mtd_data->pos[iat][2];
    mass2 += mtd_data->mass[iat];
  }

  sum1[0] /= mass1; sum1[1] /= mass1; sum1[2] /= mass1;
  sum2[0] /= mass2; sum2[1] /= mass2; sum2[2] /= mass2;
  if(colvar.cell_pbc[i_c]){
    minimal_image(sum1, sum2, &mod_rij, rij);
  } else {
    rij[0] = sum1[0]-sum2[0];
    rij[1] = sum1[1]-sum2[1];
    rij[2] = sum1[2]-sum2[2];
    mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
  };

  for(i=0;i<colvar.list[i_c][0];i++) { 
    iat = colvar.cvatoms[i_c][i];
    colvar.myder[i_c][i][0] =  (mtd_data->mass[iat]*rij[0])/(mass1*mod_rij);
    colvar.myder[i_c][i][1] =  (mtd_data->mass[iat]*rij[1])/(mass1*mod_rij);
    colvar.myder[i_c][i][2] =  (mtd_data->mass[iat]*rij[2])/(mass1*mod_rij);
  }
  for(i=colvar.list[i_c][0];i<colvar.natoms[i_c];i++) {
    iat = colvar.cvatoms[i_c][i];
    colvar.myder[i_c][i][0] = -(mtd_data->mass[iat]*rij[0])/(mass2*mod_rij);
    colvar.myder[i_c][i][1] = -(mtd_data->mass[iat]*rij[1])/(mass2*mod_rij);
    colvar.myder[i_c][i][2] = -(mtd_data->mass[iat]*rij[2])/(mass2*mod_rij);
  }

  colvar.ss0[i_c] = mod_rij;

}

// ------------------------------------------------------------------------------------------------

int PREFIX read_dist(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iw, iat,j;
  double delta = 0.0;
  char string[400];
  int help;
  help=0;

  colvar.cell_pbc[count]=1; // default is PBC

  iw=seek_word(word,"NOPBC");
  if(iw>=0) {colvar.cell_pbc[count] = 0;}
  iw=seek_word(word,"PBC");
  if(iw>=0) {colvar.cell_pbc[count] = 1;}

  iw = seek_word(word,"LIST");
  if(iw>=0){   
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
             j=plumed_get_group(word[iw+2],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][1]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR DISTANCE\n"); help=1;}
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);  
             colvar.delta_r[count]  = (real) delta; }
 // else if (logical.do_hills) { fprintf(fplog,"|- NEEDED SIGMA KEYWORD FOR DISTANCE\n"); help=1;}
  if(help){
          fprintf(fplog, "\n-DISTANCE CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "      DISTANCE LIST 15 30 SIGMA 1.0 \n");
          fprintf(fplog, "  \n");
          fprintf(fplog, "or in case of groups    \n");
          fprintf(fplog, "      DISTANCE LIST <g1> <g2> SIGMA 1.0 \n");
          fprintf(fplog, "         g1->    \n");
          fprintf(fplog, "         6 10    \n");
          fprintf(fplog, "         g1<-    \n");
          fprintf(fplog, "                 \n");
          fprintf(fplog, "         g2->    \n");
          fprintf(fplog, "         8 15 21 \n");
          fprintf(fplog, "         g2<-    \n");
          fprintf(fplog, "                 \n");
          fprintf(stderr, "PluMed dead with errors: check log file  \n");
          EXIT();
  }

  colvar.type_s[count]   = 1;

  fprintf(fplog, "\n%1i-DISTANCE: (1st SET: %i ATOMS), (2nd SET: %i ATOMS); ", count+1, colvar.list[count][0], colvar.list[count][1]);
  if(colvar.cell_pbc[count]) fprintf(fplog, " PBC ON");
  else                       fprintf(fplog, " PBC OFF");
  if (logical.do_hills){
 	if (colvar.delta_r[count]>0){
        	 fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
        }
  }
  else fprintf(fplog,"\n");

  iat=0;
  fprintf(fplog,"|- 1st SET MEMBERS: ");
  for(i=0;i<colvar.list[count][0];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 2nd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][1];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");

  snew(colvar.myder[count], colvar.natoms[count]);

  return colvar.natoms[count];
}
