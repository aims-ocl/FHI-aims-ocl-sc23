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

void PREFIX radgyr_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int i, j, firstAtom;
  rvec rcom, rdiff, pos0;
  real totmass, Rg, mod_rij;
  real inertia_tensor[3][3];
  real trace;

  totmass = Rg = trace = rcom[0] = rcom[1] = rcom[2] = 0.;

  /* Initialise inertia tensor matrix */
//  if(colvar.nn[i_c]>1)
//    for(i=0;i<=2;i++)
//      for(j=0;j<=2;j++)
//        inertia_tensor[i][j]=0.0;


  firstAtom = colvar.cvatoms[i_c][0];
  totmass += mtd_data->mass[firstAtom];
  for(i=0;i<3;i++) pos0[i] = mtd_data->pos[firstAtom][i];

  for(i=1;i<colvar.natoms[i_c];i++){
    firstAtom = colvar.cvatoms[i_c][i];
    totmass += mtd_data->mass[firstAtom];
    if(colvar.cell_pbc[i_c]) {
     minimal_image(mtd_data->pos[firstAtom], pos0, &mod_rij, rdiff);
     for(j=0;j<3;j++) rcom[j] += mtd_data->mass[firstAtom]*rdiff[j];
    } else {
     for(j=0;j<3;j++) rcom[j] += mtd_data->mass[firstAtom]*(mtd_data->pos[firstAtom][j]-pos0[j]);
    }
  }

  for(j=0;j<3;j++) rcom[j] = rcom[j]/totmass+pos0[j];

  for(i=0;i<colvar.natoms[i_c];i++){
    firstAtom = colvar.cvatoms[i_c][i];
    if(colvar.cell_pbc[i_c]) {
     minimal_image(mtd_data->pos[firstAtom], rcom, &mod_rij, rdiff);
    } else {
     for(j=0;j<3;j++) rdiff[j] = mtd_data->pos[firstAtom][j]-rcom[j];
     mod_rij  = sqrt(rdiff[0]*rdiff[0]+rdiff[1]*rdiff[1]+rdiff[2]*rdiff[2]);    
    }
 
    if(colvar.nn[i_c]<2){
      Rg += mtd_data->mass[firstAtom]*mod_rij*mod_rij;
      for(j=0;j<3;j++) colvar.myder[i_c][i][j] = rdiff[j]*mtd_data->mass[firstAtom]; 
    }else{
//      inertia_tensor[0][0]+=mtd_data->mass[firstAtom]*(rdiff[2]*rdiff[2]+rdiff[3]*rdiff[3]);
//      inertia_tensor[1][1]+=mtd_data->mass[firstAtom]*(rdiff[1]*rdiff[1]+rdiff[3]*rdiff[3]);
//      inertia_tensor[2][2]+=mtd_data->mass[firstAtom]*(rdiff[1]*rdiff[1]+rdiff[2]*rdiff[2]);
//      inertia_tensor[0][1]+=mtd_data->mass[firstAtom]*(rdiff[0]*rdiff[1]);
//      inertia_tensor[0][2]+=mtd_data->mass[firstAtom]*(rdiff[0]*rdiff[2]);
//      inertia_tensor[1][2]+=mtd_data->mass[firstAtom]*(rdiff[1]*rdiff[2]);
    }
  }

// trace of the inertia tensor
    if(colvar.nn[i_c]==0){
      //colvar.ss0[i_c] = 2*Rg/totmass;
      colvar.ss0[i_c] = 2*Rg;
      for(i=0;i<colvar.natoms[i_c];i++) {
       for(j=0;j<3;j++) colvar.myder[i_c][i][j] *= 4;  
      }
// gyration radius
    }else if(colvar.nn[i_c]==1){
      colvar.ss0[i_c] = sqrt(Rg/totmass);
      for(i=0;i<colvar.natoms[i_c];i++) {
       for(j=0;j<3;j++) colvar.myder[i_c][i][j] /= colvar.ss0[i_c]*totmass; 
      }
    }
}

//------------------------------------------------------------------------------------------------

int PREFIX read_rgyr(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iw , j;
  double delta = 0.0;
  char string[400];
  int help;

  help=0;
  colvar.nn[count] = 0;

  iw=seek_word(word,"RGYR");
  if(iw>=0)colvar.nn[count] = 1;

  iw = seek_word(word,"LIST");
  if(iw>=0){
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR RGYR\n"); help=1;} 

  iw=seek_word(word,"TYPE");
  if(iw>0){
    if(seek_word2(word,"TRACE",iw)>iw)colvar.nn[count] = 0;
    if(seek_word2(word,"RGYR",iw)>iw)colvar.nn[count] = 1;
  }

  iw=seek_word(word,"PBC"); 
  if(iw>0) colvar.cell_pbc[count] = 1;
 
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }

  if(help){
          printf("%i \n",help);
          fprintf(fplog, "\n-INERTIA/RGYR CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "INERTIA/RGYR LIST <g1> SIGMA 0.1 [PBC] [TYPE RGYR/TRACE]\n");
          fprintf(fplog, "         g1->            \n");
          fprintf(fplog, "         6 10 16         \n");
          fprintf(fplog, "         LOOP 20 30 2    \n");
          fprintf(fplog, "         g1<-            \n");
          fprintf(fplog, "                         \n");
          fprintf(stderr, "PluMed dead with errors: check log file  \n");
          EXIT(); 
  }

  colvar.type_s[count]   = 11;

  snew(colvar.myder[count], colvar.natoms[count]);


  if(colvar.nn[count] == 0)
    fprintf(fplog,"%i-TRACE OF THE INERTIA TENSOR; ATOMS INVOLVED: %i; ", count+1, colvar.natoms[count]);

  if(colvar.nn[count] == 1)
    fprintf(fplog,"%i-GYRATION RADIUS; ATOMS INVOLVED: %i; ", count+1, colvar.natoms[count]);

  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  
  if (colvar.cell_pbc[count]) {fprintf(fplog,"|- PBC is ON \n");}
  else {fprintf(fplog,"|- PBC is OFF \n");} 

  fprintf(fplog,"|- SET MEMBERS: ");
  for(i=0;i<colvar.natoms[count];i++){fprintf(fplog," %d ",colvar.cvatoms[count][i]+1);if((i+1)%20==0)fprintf(fplog,"\n               ");}fprintf(fplog,"\n\n");

  return colvar.natoms[count];
}
