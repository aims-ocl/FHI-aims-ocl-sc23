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
#include <assert.h>

void PREFIX read_restraint(struct mtd_data_s *mtd_data)
{
  double uno, due, tre, quattro;
  int i, icv, count, exp, iw, nw, iline;
  FILE *file;
  char metafile[120], **word, tmpmeta[120];

// object containing parsed input
  t_plumed_input input;

// open PluMeD parameters file  
  file = fopen(mtd_data->metaFilename, "r");
  if(file==NULL) {
    if(mtd_data->repl==-1){
      char buf[1024];
      sprintf(buf, "MISSING PLUMED INPUT FILE %s",mtd_data->metaFilename);
      plumed_error(buf);
    }
    else if(mtd_data->repl>-1) {
      strcpy(tmpmeta,mtd_data->metaFilename);
      tmpmeta[strlen(mtd_data->metaFilename) - 4] = '\0';
      sprintf(tmpmeta+strlen(tmpmeta),"%d",mtd_data->repl);
      sprintf(metafile, "%s.dat", tmpmeta);
      file = fopen(metafile, "r");
      if(file==NULL) {
        char buf[1024];
        sprintf(buf, "MISSING PLUMED INPUT FILE %s",mtd_data->metaFilename);
        plumed_error(buf);
      }
    }
  }

// CV counter initialization and calling routine to set the default values for many variables
  count = 0;
  iline = 0;
  read_defaults();

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                               INPUT PARSER
//.............................................................................
// first word must be keyword (like COORD), then the parser
// seeks on the same line for additional input (like SIGMA).
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  fprintf(mtd_data->fplog, " \n");  
  fprintf(mtd_data->fplog, "::::::::::::::::: READING PLUMED INPUT :::::::::::::::::\n");

  plumed_read_input(&input,file,mtd_data->fplog);
// everything is now in the "input" structure, the actual file is not needed anymore
  fclose(file);

  for(iline=0;iline<input.nlines;iline++){

    nw   = input.nwords[iline];
    word = input.words[iline];

// empty line
    if(nw==0) continue;

// explicit comment, to be copied on the log file
    if(!strcmp(word[0],"NOTE") || !strcmp(word[0],"COMMENT")){
      fprintf(mtd_data->fplog, "\nCOMMENT: ");
      for(i=1;i<nw;i++)fprintf(mtd_data->fplog, "%s ", word[i]);
      fprintf(mtd_data->fplog,"\n");
// untested features
    } else if(!strcmp(word[0],"ENABLE_UNTESTED_FEATURES")){
      fprintf(mtd_data->fplog, "|-##################################################################\n");
      fprintf(mtd_data->fplog, "|- ENABLE_UNTESTED_FEATURES\n");
      fprintf(mtd_data->fplog, "|- THIS FLAG ENABLES FEATURES WHICH ARE NOT EXPLAINED IN THE MANUAL\n");
      fprintf(mtd_data->fplog, "|- AND COULD BE BUGGY\n");
      fprintf(mtd_data->fplog, "|- USE IT ONLY IF YOU ARE A PLUMED DEVELOPERS\n");
      fprintf(mtd_data->fplog, "|-##################################################################\n");
      logical.enable_untested_features=1;
// commitment analysis
    } else if(!strcmp(word[0],"COMMITMENT")){
      logical.commitment = 1;
      fprintf(mtd_data->fplog, "|-COMMITOR ANALYSIS: YOU WILL ONLY MONITOR YOUR CVs MICRODYNAMICS\n");
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"NCV")){
          iw++; sscanf(word[iw],"%i",&colvar.commit);
        } else {
          plumed_error("Unknown option for COMMITMENT keyword");
        }
      }
      for(i=0;i<colvar.commit;i++){
        iline++; // this mimicks a fgets
        sscanf(input.words[iline][0],"%lf",&uno);
        sscanf(input.words[iline][1],"%lf",&due);
        sscanf(input.words[iline][2],"%lf",&tre);
        sscanf(input.words[iline][3],"%lf",&quattro);
        colvar.Amin[i] = (real) uno;
        colvar.Amax[i] = (real) due;
        colvar.Bmin[i] = (real) tre;
        colvar.Bmax[i] = (real) quattro;
        fprintf(mtd_data->fplog, "|--CV %i: A min %f, max %f -- B min %f, max %f\n", i, colvar.Amin[i], colvar.Amax[i], colvar.Bmin[i], colvar.Bmax[i]);
      }
      fprintf(mtd_data->fplog, "\n");
// metadynamics
    } else if(!strcmp(word[0],"HILLS")) {
      logical.do_hills = 1;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"HEIGHT")){
          iw++; sscanf(word[iw],"%lf",&uno); hills.wwr = (real) uno*mtd_data->eunit;
        } else if(!strcmp(word[iw],"W_STRIDE")){
          iw++; sscanf(word[iw],"%i",&hills.nt_hills);
        } else if(!strcmp(word[iw],"R_STRIDE")){
          iw++; sscanf(word[iw],"%i",&hills.nr_hills);
        } else if(!strcmp(word[iw],"RESTART")){
          logical.restart_hills = 1;
        } else if(!strcmp(word[iw],"RATE")){
          iw++; sscanf(word[iw],"%lf",&uno); hills.rate = (real) uno;
        } else if(!strcmp(word[iw],"MAX_HEIGHT")){
          iw++; sscanf(word[iw],"%lf",&uno); hills.max_height = (real) uno;
        } else if(!strcmp(word[iw],"MAX_STRIDE")){
          iw++; sscanf(word[iw],"%i",&hills.max_stride);
        } else {
          plumed_error("Unknown option for HILLS keyword");
        }
      };
      if(hills.rate==0.0) hills.rate = hills.wwr/hills.nt_hills/mtd_data->dt;
      if(hills.wwr==0.0) hills.wwr=hills.rate*hills.nt_hills*mtd_data->dt;
      fprintf(mtd_data->fplog,"|-HILLS:\n");
      if(hills.max_height==0.0){
        fprintf(mtd_data->fplog,"|--HEIGHT %f  WRITING STRIDE %i DEPOSITION RATE %f \n",
                hills.wwr/mtd_data->eunit, hills.nt_hills, hills.rate/mtd_data->eunit);
      } else {
        fprintf(mtd_data->fplog,"|--DEPOSITION RATE %f \n",hills.rate/mtd_data->eunit);
        fprintf(mtd_data->fplog,"|--MAXIMUM STRIDE BETWEEN HILLS %i \n",hills.max_stride);
        fprintf(mtd_data->fplog,"|--MAXIMUM HEIGHT               %f \n",hills.max_height);
      }
      if(logical.restart_hills) fprintf(mtd_data->fplog,"|-RESTARTING METADYNAMICS!\n");
      if(hills.nr_hills!=1) fprintf(mtd_data->fplog,"|--READING STRIDE %i\n", hills.nr_hills);
      fprintf(mtd_data->fplog, "\n");
// multiple walkers (shared file)
    } else if(!strcmp(word[0],"MULTIPLE_WALKERS")) {
      nwalkers=2;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"HILLS_DIR")){
          iw++; sscanf(word[iw], "%s", hills.dir);
        } else if(!strcmp(word[iw],"WALKERS")){
          iw++; iw++;
        } else {
          plumed_error("Unknown option for MULTIPLE_WALKERS keyword");
        }
      };
      fprintf(mtd_data->fplog,"|-MULTIPLE WALKERS:\n DIRECTORY FOR HILLS I/O %s\n\n", hills.dir);
    } else if(!strcmp(word[0],"PRINT")) {
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"W_STRIDE")){
          iw++; sscanf(word[iw],"%i", &colvar.nt_print);
        } else {
          plumed_error("Unknown option for PRINT keyword");
        }
      }
      fprintf(mtd_data->fplog,"|-PRINTING ON COLVAR FILE EVERY %i STEPS\n",colvar.nt_print);
      logical.print = 1;
    } else if(!strcmp(word[0],"ALIGN_ATOMS")) {
       int list_found;
       list_found=0;
       fprintf(mtd_data->fplog,"|- ALIGNING ATOMS\n");
       for(iw=1;iw<nw;iw++){
         if(!strcmp(word[iw],"LIST")){
           list_found=1;
           iw++; colvar.align_atoms+=plumed_get_group(word[iw],&colvar.align_list,colvar.align_atoms,&input,mtd_data->fplog);
         } else {
           plumed_error("Unknown option for ALIGN keyword");
         }
       }
       if(!list_found)plumed_error("NEEDED LIST KEYWORD FOR ALIGN_ATOMS\n");
       fprintf(mtd_data->fplog,"|- SET MEMBERS: ");
       for(i=0;i<colvar.align_atoms;i++){
         fprintf(mtd_data->fplog," %d ",colvar.align_list[i]+1);if((i+1)%20==0)fprintf(mtd_data->fplog,"\n               ");
       }fprintf(mtd_data->fplog,"\n\n");
    } else if(!strcmp(word[0],"DEBUG_DERIVATIVES")){
      logical.debug_derivatives=1;
    } else if(!strcmp(word[0],"PARALLEL_HILLS")){
#if ! defined (PLUMED_GROMACS) && ! defined (DL_POLY)
        plumed_error("PARALLEL_HILLS NOT YET IMPLEMENTED IN THIS CODE");
#endif
        for(iw=1;iw<nw;iw++){
          if(!strcmp(word[iw],"ON")){
            logical.parallel_hills=1;
          } else if(!strcmp(word[iw],"OFF")){
            logical.parallel_hills=0;
          } else {
            plumed_error("Unknown flag for keyword PARALLEL_HILLS");
          }
        }
        fprintf(mtd_data->fplog, "|-PARALLEL HILLS ");
        if(logical.parallel_hills) fprintf(mtd_data->fplog, "ON");
        else                       fprintf(mtd_data->fplog, "OFF");
        fprintf(mtd_data->fplog, "\n\n");
    } else if(!strcmp(word[0],"PTMETAD")){
      #if defined (NAMD) || defined (DL_POLY) || defined (AMBER)  || defined(FHIAIMS)
          plumed_error("PTMETAD: NOT YET IMPLEMENTED IN THIS CODE");
      #else
      fprintf(mtd_data->fplog, "|-PARALLEL TEMPERING METADYNAMICS\n");
      logical.remd = 1;
      fprintf(mtd_data->fplog, "|--REPLICA 0 TEMPERATURE = %f\n", mtd_data->rte0);
      fprintf(mtd_data->fplog, "|--REPLICA %i TEMPERATURE = %f\n", mtd_data->repl, mtd_data->rteio);
      if(mtd_data->repl==-1) {
        fprintf(mtd_data->fplog, "\n!!!! mdrun not in replica exchange mode, keyword PTMETAD will not be considered !!!!\n");
        logical.remd = 0;
      }
      fprintf(mtd_data->fplog, "\n");
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"NEIGHBOUR")){
          iw++; sscanf(word[iw],"%i",&colvar.ptmetad_neighbours);
        } else if(!strcmp(word[iw],"SIGMA")){
          iw++; sscanf(word[iw],"%lf",&uno); colvar.ptmetad_sigma = (real) uno;
        } else {
          plumed_error("Unknown flag for keyword PTMETAD");
        }
      }
      if(colvar.ptmetad_neighbours){
        fprintf(mtd_data->fplog, "|--SIGMA = %lf\n",(double) colvar.ptmetad_sigma);
        fprintf(mtd_data->fplog, "|--NEIGHBOUR = %i\n",colvar.ptmetad_neighbours);
      }
      fprintf(mtd_data->fplog, "\n");
      #endif
    } else if(!strcmp(word[0],"BIASXMD")){
      #if defined (NAMD) || defined (DL_POLY) || defined (AMBER)
          plumed_error("|-BIASXMD: NOT YET IMPLEMENTED IN THIS CODE");
      #else
      fprintf(mtd_data->fplog, "|-BIAS EXCHANGE METADYNAMICS\n");
      logical.rpxm = 1;
      logical.remd = 1;
      if(mtd_data->repl==-1) {
        fprintf(mtd_data->fplog, "\n!!!! mdrun not in replica exchange mode, keyword BIASXMD will not be considered !!!!\n");
        logical.rpxm = 0;
        logical.remd = 0;
      }
      fprintf(mtd_data->fplog, "\n");
      #endif
    } else if(!strcmp(word[0],"WELLTEMPERED")){
      logical.welltemp = 1;
      int read_biasfactor = 0;
      int read_cvtemp = 0;
      int read_simtemp = 0;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CVTEMP")){
          iw++; sscanf(word[iw], "%lf", &uno); colvar.wtemp   = (real) uno; read_cvtemp = 1;
        } else if(!strcmp(word[iw],"BIASFACTOR")){
          iw++; sscanf(word[iw], "%lf", &uno); colvar.wfactor = (real) uno; read_biasfactor = 1;
        } else if(!strcmp(word[iw],"SIMTEMP")){
          iw++; sscanf(word[iw], "%lf", &uno); colvar.simtemp = (real) uno; read_simtemp = 1;
        } else {
          plumed_error("Unknown flag for keyword WELLTEMPERED");
        }
      }
      if(!read_simtemp)
        plumed_error("WITH WELLTEMPERED YOU ALWAYS HAVE TO SPECIFY THE \"SIMTEMP \" KEYWORD");
      if(read_biasfactor==read_cvtemp)
        plumed_error("WITH WELLTEMPERED YOU HAVE TO SPECIFY EITHER \"CVTEMP \" OR \"BIASFACTOR \" KEYWORD");
      if(read_cvtemp)     colvar.wfactor = colvar.wtemp / colvar.simtemp;  
      if(read_biasfactor) colvar.wtemp   = colvar.wfactor * colvar.simtemp;

      if (colvar.wfactor<=1.0) {  
        char buf[1024];
        sprintf(buf,"WELLTEMPERED, temperature factor less than or equal to 1.0 ( %f ) \n",colvar.wfactor);
        plumed_error(buf);
      } 

      fprintf(mtd_data->fplog, "|-WELL TEMPERED METADYNAMICS WITH BIAS FACTOR %f (CVTEMP = %f) \n\n", colvar.wfactor,colvar.wtemp);
   
    } else if(!strcmp(word[0],"DEBUG")){
      fprintf(mtd_data->fplog,"|- CV DERIVATIVES DEBUGGING \n"); 
      logical.debug = 1;
    } else if(!strcmp(word[0],"DISTANCE")){
      read_dist(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"MINDIST")){
      read_mindist(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"COORD")){
      read_coord(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"ANGLE")){
      read_angle(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"HBONDS")){
      read_hbonds(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"TORSION")){
      read_torsion(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"RGYR") || !strcmp(word[0],"INERTIA")){
      read_rgyr(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"RMSDTOR")){
      read_rmsdtor(word, count, &input,&iline,mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"DIPOLE")){
      read_dipole(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"DIHCOR")) {
      read_dihcor(word, count, &input,&iline,mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"WATERBRIDGE")) {
      read_waterbridge(word, count, &input, mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"ALPHABETA")) {
      read_alfabeta(word, count,&input,&iline,mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"S_PATH")) {
      colvar.type_s[count]   = 30;
      read_path(word, count, &input,mtd_data->fplog);
      count++;
    } else if(!strcmp(word[0],"Z_PATH")) {
      colvar.type_s[count]   = 31;
      read_path(word, count, &input,mtd_data->fplog);
      count++;   
    } else if(!strcmp(word[0],"TARGETED")) {
      colvar.type_s[count]   = 31;
      read_path(word, count, &input,mtd_data->fplog);
      count++;   
/* Atom position */      
    } else if(!strcmp(word[0],"POSITION")) {
      colvar.type_s[count]   = 32;
      read_position(word, count, &input, mtd_data->fplog);
      count++;   
    } else if(!strcmp(word[0],"ELSTPOT")) {
      colvar.type_s[count]   = 33;
      read_elstpot(word, count, &input, mtd_data->fplog);
      count++;   
    } else if(!strcmp(word[0],"PUCKERING")) {
      colvar.type_s[count]   = 34;
      read_puckering(word, count, &input, mtd_data->fplog);
      count++;   
    } else if(!strcmp(word[0],"UWALL")){
      int read_limit=0;
      int read_kappa=0;
// first we select the proper CV
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv);}
      else{plumed_error("WITH UWALL YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");}
// then we parse the line
      logical.upper[icv-1]=1;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CV")) {
          iw++;  // already read
        } else if(!strcmp(word[iw],"LIMIT")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.upper[icv-1]=(real)uno; read_limit=1;
        } else if(!strcmp(word[iw],"KAPPA")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.sigma[icv-1]=(real)uno*mtd_data->eunit; read_kappa=1;
        } else if(!strcmp(word[iw],"EXP")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.uexp[icv-1]=(real)uno;
        } else if(!strcmp(word[iw],"EPS")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.ueps[icv-1]=(real)uno;
        } else if(!strcmp(word[iw],"OFF")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.uoff[icv-1]=(real)uno;
        } else {
          plumed_error("Unknown flag for keyword UWALL");
        };
      }
      if(!read_limit)
        plumed_error("WITH UWALL YOU ALWAYS HAVE TO SPECIFY THE \"LIMIT\" KEYWORD\n");
      if(!read_kappa)
        plumed_error("WITH UWALL YOU ALWAYS HAVE TO SPECIFY THE \"KAPPA\" KEYWORD\n");
      fprintf(mtd_data->fplog, "|-WALL ON COLVAR %i: UPPER LIMIT = %f, KAPPA = %f, EXPONENT = %i, REDUX = %f, OFFSET = %f \n\n",
             icv, cvw.upper[icv-1], cvw.sigma[icv-1]/mtd_data->eunit, cvw.uexp[icv-1], cvw.ueps[icv-1], cvw.uoff[icv-1]);
     } else if(!strcmp(word[0],"LWALL")){
      int read_limit=0;
      int read_kappa=0;
// first we select the proper CV
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv);}
      else{plumed_error("WITH UWALL YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");}
// then we parse the line
      logical.lower[icv-1]=1;
      for(iw=1;iw<nw;iw++){
        if(!strcmp(word[iw],"CV")) {
          iw++;  // already read
        } else if(!strcmp(word[iw],"LIMIT")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.lower[icv-1]=(real)uno; read_limit=1;
        } else if(!strcmp(word[iw],"KAPPA")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.lsigma[icv-1]=(real)uno*mtd_data->eunit; read_kappa=1;
        } else if(!strcmp(word[iw],"EXP")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.lexp[icv-1]=(real)uno;
        } else if(!strcmp(word[iw],"EPS")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.leps[icv-1]=(real)uno;
        } else if(!strcmp(word[iw],"OFF")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvw.loff[icv-1]=(real)uno;
        } else {
          plumed_error("Unknown flag for keyword LWALL");
        };
      }
      if(!read_limit)
        plumed_error("WITH LWALL YOU ALWAYS HAVE TO SPECIFY THE \"LIMIT\" KEYWORD\n");
      if(!read_kappa)
        plumed_error("WITH LWALL YOU ALWAYS HAVE TO SPECIFY THE \"KAPPA\" KEYWORD\n");
      fprintf(mtd_data->fplog, "|-WALL ON COLVAR %i: LOWER LIMIT = %f, KAPPA = %f, EXPONENT = %i, REDUX = %f, OFFSET = %f \n\n",
             icv, cvw.lower[icv-1], cvw.lsigma[icv-1]/mtd_data->eunit, cvw.lexp[icv-1], cvw.leps[icv-1], cvw.loff[icv-1]);
    } else if(!strcmp(word[0],"STEER")){
      int read_max=0;
      int read_delta=0;
      int read_kappa=0;
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv); cvsteer.impose_start[icv-1] = 0;}
      else {plumed_error("WITH STEER YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");} 
      for(iw=1;iw<nw;iw++){ 
        if(!strcmp(word[iw],"CV")){
          iw++; // already read
        } else if(!strcmp(word[iw],"FROM")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvsteer.start[icv-1]=(real)uno; cvsteer.impose_start[icv-1]=1;
        } else if(!strcmp(word[iw],"TO")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvsteer.max[icv-1]=(real)uno; read_max=1;
        } else if(!strcmp(word[iw],"VEL")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvsteer.delta[icv-1]=(real) fabs(uno); read_delta=1;
        } else if(!strcmp(word[iw],"KAPPA")) {
          iw++; sscanf(word[iw], "%lf", &uno); cvsteer.spring[icv-1]=(real) uno*mtd_data->eunit; read_kappa=1;
        } else {
          plumed_error("Unknown flag for keyword LWALL");
        }
      }
      if(!read_max)  plumed_error("WITH STEER YOU ALWAYS HAVE TO SPECIFY THE \"TO\" KEYWORD\n");
      if(!read_delta)plumed_error("WITH STEER YOU ALWAYS HAVE TO SPECIFY THE \"VEL\" KEYWORD\n");
      if(!read_kappa)plumed_error("WITH STEER YOU ALWAYS HAVE TO SPECIFY THE \"KAPPA\" KEYWORD\n");

      logical.steer[icv-1]  = 1; 
      if(cvsteer.impose_start[icv-1]==0){
        fprintf(mtd_data->fplog, "|-STEERING COLVAR %i TO %f: VELOCITY=%lf cvunit/kstep, SPRING=%lf\n\n", icv,cvsteer.max[icv-1],cvsteer.delta[icv-1],cvsteer.spring[icv-1]/mtd_data->eunit);
      } else { 
        fprintf(mtd_data->fplog, "|-STEERING COLVAR %i FROM %f TO %f: VELOCITY=%lf cvunit/kstep, SPRING=%lf\n\n", icv, cvsteer.start[icv-1],cvsteer.max[icv-1],cvsteer.delta[icv-1],cvsteer.spring[icv-1]/mtd_data->eunit);
      } 

    } else if(!strcmp(word[0],"UMBRELLA")){
      iw = seek_word(word,"CV");
      if(iw>=0){ sscanf(word[iw+1], "%i", &icv); 
        cvsteer.impose_start[icv-1]  = 1;
        cvsteer.delta[icv-1]         = 0.;
      }
      else {plumed_error("WITH UMBRELLA YOU ALWAYS HAVE TO SPECIFY THE \"CV\" KEYWORD\n");} 
      iw = seek_word(word,"KAPPA");
      if(iw>=0){ sscanf(word[iw+1], "%lf", &tre);cvsteer.spring[icv-1] = (real) tre*mtd_data->eunit;}
      else{plumed_error("WITH UMBRELLA YOU ALWAYS HAVE TO SPECIFY THE \"KAPPA\" KEYWORD\n");}
      iw = seek_word(word,"AT");
      if(iw>=0){ sscanf(word[iw+1], "%lf", &quattro);cvsteer.max[icv-1]=cvsteer.start[icv-1]=(real) quattro;}
      else{plumed_error("WITH UMBRELLA YOU ALWAYS HAVE TO SPECIFY THE \"AT\" KEYWORD\n");}   

      logical.steer[icv-1]  = 1; 
      fprintf(mtd_data->fplog, "|-UMBRELLA SAMPLING OF COLVAR %i AT %f: SPRING=%lf\n\n", icv,cvsteer.max[icv-1],cvsteer.spring[icv-1]/mtd_data->eunit);

    } else if(!strcmp(word[0],"NOHILLS")){
      iw = seek_word(word,"CV");
      if(iw>=0) sscanf(word[iw+1], "%i", &icv);
      colvar.on[icv-1] = 0;
      fprintf(mtd_data->fplog, "|-NO HILLS ON COLVAR %i\n", icv);
    } else if(!strcmp(word[0],"UREFLECT")){
      iw = seek_word(word,"CV");
      if(iw>=0) sscanf(word[iw+1], "%i", &icv);
      iw = seek_word(word,"LIMIT");
      if(iw>=0) sscanf(word[iw+1], "%lf", &uno);
      cvw.upper[icv-1] = (real) uno;
      logical.ureflect[icv-1] = 1;
      fprintf(mtd_data->fplog, "|-UPPER REFLECTING WALL ON CV %i, AT %f\n\n", icv, cvw.upper[icv-1]);
    } else if(!strcmp(word[0],"LREFLECT")){
      iw = seek_word(word,"CV");
      if(iw>=0) sscanf(word[iw+1], "%i", &icv);
      iw = seek_word(word,"LIMIT");
      if(iw>=0) sscanf(word[iw+1], "%lf", &uno);
      cvw.lower[icv-1] = (real) uno;
      logical.lreflect[icv-1] = 1;
      fprintf(mtd_data->fplog, "|-LOWER REFLECTING WALL ON CV %i, AT %f\n\n", icv, cvw.lower[icv-1]);
    } else if(!strcmp(word[0],"DEBUG_GRID")){
      logical.debug_grid=1;
    } else if(!strcmp(word[0],"NOSPLINE")){
      logical.donot_spline=1;
      fprintf(mtd_data->fplog, "|- GRID SPLINE TURNED OFF\n");
    } else if(!strcmp(word[0],"GRID")){
      iw = seek_word(word,"CV");
      if(iw>=0) sscanf(word[iw+1], "%i", &icv);
      iw = seek_word(word,"MIN");
      if(iw>=0) sscanf(word[iw+1], "%s", &tmpmeta);
      uno = plumed_atof(tmpmeta); 
      grid.min[grid.ncv] = (real) uno;
      iw = seek_word(word,"MAX");
      if(iw>=0) sscanf(word[iw+1], "%s", &tmpmeta);
      due = plumed_atof(tmpmeta); 
      grid.max[grid.ncv] = (real) due;
      iw = seek_word(word,"NBIN");
      if(iw>=0) sscanf(word[iw+1], "%i", &(grid.bin[grid.ncv]));
      iw = seek_word(word,"PBC"); if(iw>=0) grid.period[grid.ncv] = 1; 
      fprintf(mtd_data->fplog, "|-GRID ACTIVE ON CV %i NBIN %d MIN %f MAX %f \n", icv, grid.bin[grid.ncv],grid.min[grid.ncv],grid.max[grid.ncv]);
      if(grid.period[grid.ncv]) fprintf(mtd_data->fplog, "|-- PERIODIC GRID IS ON\n");
      grid.index[grid.ncv] = icv-1;
      for(i=0;i<grid.ncv;i++) if(grid.index[i]==grid.index[grid.ncv]) plumed_error("GRID is already ACTIVE for this CV");
      grid.ncv += 1;
      logical.do_grid = 1;
    } else if(!strcmp(word[0],"PROJ_GRAD")){ 
       iw = seek_word(word,"CV"); // look for a group
       colvar.pg.nlist=plumed_get_group(word[iw+1],&colvar.pg.list,0,&input,mtd_data->fplog); 
    } else {
      char buf[1024];
      sprintf(buf, "Line %i Unkwown Keyword %s \n", iline+1, word[0]);
      plumed_error(buf);
    }
  }
// clean input parser
  plumed_clear_input(&input);

// set the number of collective variables
  colvar.nconst = count;

// unset the SIGMA<0 CVs for hills
  for(i=0;i<colvar.nconst;i++){
     if(colvar.delta_r[i]<0.){colvar.on[i]=0;}
     if(colvar.on[i]==0){           fprintf(mtd_data->fplog, "|-NO HILLS     ON COLVAR %i\n", i+1);} 
     else               {           fprintf(mtd_data->fplog, "|-HILLS ACTIVE ON COLVAR %i\n", i+1);}
  }

// check the correctenes of the input parsed
  if(colvar.nconst > nconst_max) {
    plumed_error("Too many colvars. Change NCONST_MAX in metadyn.h !!!!!!!!!!!\n");	
  }

// checking for conflicts in directive keywords
  if(!logical.do_hills && !logical.commitment){
    fprintf(mtd_data->fplog, "|-ANALYSIS: YOU WILL ONLY MONITOR YOUR CVs MICRODYNAMICS\n\n");
  }

  if(logical.welltemp && !logical.do_hills)  plumed_error("WELLTEMPERED must be used with HILLS keyword");

  if(logical.commitment && logical.do_hills) plumed_error("KEYWORD 'COMMITMENT' AND 'HILLS' ARE NOT COMPATIBLE");

// in case of parallel or solute tempering rescale hills heigth with temperature
  if(logical.do_hills&&logical.remd&&(!logical.rpxm)) hills.wwr *= mtd_data->rteio/mtd_data->rte0;

// in case of PTMETAD and well tempered set the right simtemp
  if(logical.do_hills&&logical.remd&&(!logical.rpxm)&&logical.welltemp) colvar.simtemp = mtd_data->rteio;

// check for untested features
  if(!logical.enable_untested_features) {
   if(logical.debug_derivatives) plumed_error("DEBUG_DERIVATIVES NOT ENABLED");
   if(colvar.ptmetad_neighbours) plumed_error("NEIGHBOUR HILLS NOT ENABLED");
   if(logical.debug_grid)        plumed_error("DEBUG_GRID NOT ENABLED");
   if(hills.max_height>0.0)      plumed_error("MAX_HEIGHT NOT ENABLED");
  }

// checking if GRID and HILLS active variables are consistent
  if(logical.do_grid) {
   icv = 0;
   for(i=0;i<colvar.nconst;i++) if(colvar.on[i]) icv++;  
   if(icv!=grid.ncv) plumed_error("Inconsistency between GRID and HILLS variables. Please, check !!!!!!!!!!!\n"); 
   for(i=0;i<grid.ncv;i++) if(!colvar.on[grid.index[i]] || grid.index[i]>=colvar.nconst) 
     plumed_error("Inconsistency between GRID and HILLS variables. Please, check !!!!!!!!!!!\n");
// in case initialize grid stuff
   grid_initialize(&grid);
  } 

// check for needed projection
  if(colvar.pg.nlist!=0){
        // make the projection tables
       fprintf(mtd_data->fplog, "|- FOUND PROJ_GRAD KEYWORD: NCV involved %d\n",colvar.pg.nlist);
       int j; 
       fprintf(mtd_data->fplog, "|- WHICH ARE: ");
       for(j=0;j<colvar.pg.nlist;j++){fprintf(mtd_data->fplog, " %d",colvar.pg.list[j]);}
       fprintf(mtd_data->fplog, "\n");
       setup_projections( &(colvar.pg));          
  } 

// printout PLEASE_CITE
  cite_please("bono+09cpc",mtd_data->fplog);
  if(logical.do_hills) cite_please("laio-parr02pnas",mtd_data->fplog);
  if(logical.remd && !logical.rpxm) cite_please("buss+06jacs",mtd_data->fplog);
  if(logical.rpxm) cite_please("pian-laio07jpcb",mtd_data->fplog);
  if(logical.welltemp) cite_please("bard+08prl",mtd_data->fplog);
  int logical_path=0;
  for(i=0;i<colvar.nconst;i++) if(colvar.type_s[i]==30 || colvar.type_s[i]==31) logical_path=1; 
  if(logical_path) cite_please("bran+07jcp",mtd_data->fplog);
  int logical_puckering=0;
  for(i=0;i<colvar.nconst;i++) if(colvar.type_s[i]==34) logical_puckering=1; 
  if(logical_puckering) cite_please("sega+09arxiv",mtd_data->fplog);
  if(nwalkers>1) cite_please("rait+06jpcb",mtd_data->fplog);
  fprintf(mtd_data->fplog,"\n"); 

  disclaimer(mtd_data->fplog);
// flushing output
  fflush(mtd_data->fplog);

}

//-----------------------------------------------------------------------------------------------------------------

void PREFIX read_defaults()
{
  int icv;
 
  colvar.nt_print 		= 10;
  colvar.nconst 		= 0;
  logical.restart_hills 	= 0;
  logical.remd 			= 0;
  logical.rpxm			= 0;
  logical.do_hills 		= 0;
  logical.commitment 		= 0;
  logical.print 		= 0;
  logical.widthadapt            = 0;
  logical.welltemp              = 0;
  logical.debug                 = 0;
  logical.parallel_hills        = 0;
#if defined(PLUMED_GROMACS) || defined(DL_POLY)
  logical.parallel_hills        = 1;
#endif
  logical.debug_derivatives     = 0;
  logical.enable_untested_features = 0;
  logical.do_grid               = 0;
  logical.donot_spline          = 0;
  logical.debug_grid            = 0;
  hills.wwr 			= 0.;
  hills.rate			= 0.;
  hills.max_height              = 0.;
  hills.max_stride              = 0;
  nwalkers 			= 1;
  hills.n_hills			= 0;
  hills.nt_hills                = 999999999; 
  hills.nr_hills                = 1; 
  hills.read                    = 0;
  nsz                           = 0;
  hills.first_read              = 1;
  colvar.ptmetad_neighbours     = 0;
  colvar.ptmetad_sigma          = 0.0;
  colvar.align_atoms            = 0;
  colvar.align_list             = NULL;
  grid.ncv                      = 0;
  grid.nhills                   = 0;
  colvar.pg.list		=NULL;
  colvar.pg.nlist		=0;

  for(icv=0;icv<nconst_max;icv++){
    logical.steer[icv]          = 0; 
    logical.upper[icv] 		= 0;
    logical.lower[icv] 		= 0;
    logical.ureflect[icv]       = 0;
    logical.lreflect[icv]       = 0;
    cvw.sigma[icv] 		= 0.;
    cvw.upper[icv] 		= 0.;
    cvw.lower[icv] 		= 0.;
    cvw.lsigma[icv] 		= 0.;
    cvw.fwall[icv] 		= 0;
    cvw.uexp[icv] 		= 4;
    cvw.lexp[icv] 		= 4;
    cvw.ueps[icv] 		= 1.;
    cvw.leps[icv] 		= 1.;
    cvw.uoff[icv]               = 0.;
    cvw.loff[icv]               = 0.;
    colvar.on[icv] 		= 1;
    colvar.Mss0[icv]            = 0.;
    colvar.M2ss0[icv]           = 0.;
    colvar.type_s[icv]          = 0;
    colvar.logic[icv]           = 0;
    colvar.natoms[icv]          = 0;
    colvar.cell_pbc[icv]        = 0;
    colvar.delta_r[icv]         = -1.; // default synonim of NOHILLS
    grid.min[icv]               = 0.;
    grid.max[icv]               = 0.; 
    grid.lbox[icv]              = 0.;
    grid.dx[icv]                = 0.;
    grid.bin[icv]               = 1;
    grid.minibin[icv]           = 1;
    grid.period[icv]            = 0;
    grid.index[icv]             = 0;
    grid.oldelta[icv]           = 0.;
  }

}

//-----------------------------------------------------------------------------------------------------------------

// seek_word WILL BE REMOVED SOON (as soon as it will be replaced everywhere)

int PREFIX seek_word(char **word, const char *wanted)
{
  int i;

  for (i=0;;i++) {
    if (word[i]==NULL) return -1;
    if (strcmp(word[i],wanted)==0) return i;
  }
  return -1;
}

// Added By Paolo to progrssively seek in the input string
int PREFIX seek_word2(char **word, const char *wanted, int is)
{
  int i;

  for (i=is;;i++) {
    if (word[i]==NULL) return -1;
    if (strcmp(word[i],wanted)==0) return i;
  }
  return -1;
}

void PREFIX cite_please (const char* re, FILE *fplog){


 fprintf(fplog, "\n++++ PLEASE READ AND CITE THE FOLLOWING REFERENCE ++++\n"); 
 if(!strcmp(re,"laio-parr02pnas")){
    fprintf(fplog, "  A. Laio and M. Parrinello\n");
    fprintf(fplog, "  Escaping free energy minima\n");
    fprintf(fplog, "  Proc. Natl. Acad. Sci. USA. 2002 vol. 99 (20) pp. 12562-6\n");
 } else if(!strcmp(re,"bono+09cpc")){
    fprintf(fplog, "  M. Bonomi, D. Branduardi, G. Bussi, C. Camilloni, D. Provasi, P. Raiteri, \n");
    fprintf(fplog, "  D. Donadio, F. Marinelli, F. Pietrucci, R. A. Broglia and M. Parrinello \n");
    fprintf(fplog, "  PLUMED: a portable plugin for free-energy calculations with molecular dynamics\n");
    fprintf(fplog, "  arXiv:0902.0874 [physics.comp-ph] \n");
 } else if(!strcmp(re,"pian-laio07jpcb")){
    fprintf(fplog, "  S. Piana and A. Laio\n");
    fprintf(fplog, "  A Bias-Exchange Approach to Protein Folding \n");
    fprintf(fplog, "  J. Phys. Chem. B. 2007 vol. 111 (17) pp. 4553-9\n");
 } else if(!strcmp(re,"bard+08prl")){
    fprintf(fplog, "  A. Barducci, G. Bussi and M. Parrinello\n");
    fprintf(fplog, "  Well-Tempered Metadynamics: A Smoothly Converging and Tunable Free-Energy Method \n");
    fprintf(fplog, "  Phys. Rev. Lett. 2008 vol. 100 (2) pp. 020603 \n");
 } else if(!strcmp(re,"buss+06jacs")){
    fprintf(fplog, "  G. Bussi, F.L. Gervasio, A. Laio and M. Parrinello \n");
    fprintf(fplog, "  Free-energy landscape for beta hairpin folding from combined parallel tempering and metadynamics\n");
    fprintf(fplog, "  J. Am. Chem. Soc. 2006 vol. 128 (41) pp. 13435-41 \n");
 } else if(!strcmp(re,"bran+07jcp")){
    fprintf(fplog, "  D. Branduardi, F.L. Gervasio and M. Parrinello \n");
    fprintf(fplog, "  From A to B in free energy space\n");
    fprintf(fplog, "  Jour. Chem. Phys. 2007 vol. 126 (5) pp. 054103\n");
 } else if(!strcmp(re,"rait+06jpcb")){
    fprintf(fplog, "  P. Raiteri, A. Laio, F.L. Gervasio, C. Micheletti and M. Parrinello \n");
    fprintf(fplog, "  Efficient Reconstruction of Complex Free Energy Landscapes by Multiple Walkers Metadynamics \n"); 
    fprintf(fplog, "  J. Phys. Chem. B. 2006 vol. 110 (8) pp. 3533-3539 \n");
 } else if(!strcmp(re,"sega+09arxiv")){
    fprintf(fplog, "  M. Sega, E. Autieri and F. Pederiva\n");
    fprintf(fplog, "  On the Calculation of Puckering Free Energy Surfaces \n"); 
    fprintf(fplog, "  arXiv:0904.1473 [physics.chem-ph] \n");
 } else {
    assert(1); // wrong bib name
 }

 fprintf(fplog, "-------- -------- --- Thank You --- -------- --------\n\n"); 
};

void PREFIX disclaimer (FILE *fplog){

 fprintf(fplog,"** PLUMED is free software: you can redistribute it and/or modify \n");
 fprintf(fplog,"** it under the terms of the GNU Lesser General Public License as published by \n");
 fprintf(fplog,"** the Free Software Foundation, either version 3 of the License, or \n");
 fprintf(fplog,"** (at your option) any later version. \n\n");
 fprintf(fplog,"** PLUMED is distributed in the hope that it will be useful,\n");
 fprintf(fplog,"** but WITHOUT ANY WARRANTY; without even the implied warranty of \n");
 fprintf(fplog,"** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n");
 fprintf(fplog,"** GNU Lesser General Public License for more details. \n\n");
 fprintf(fplog,"** You should have received a copy of the GNU Lesser General Public License\n");
 fprintf(fplog,"** along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.  \n\n");
 fprintf(fplog,"** For more info, see:  http://merlino.mi.infn.it/plumed \n");
 fprintf(fplog,"** or subscribe to plumed-users@googlegroups.com \n\n");
#ifdef NAMD
 fprintf(fplog,"                              WARNING!! \n");
 fprintf(fplog,"** Starting from version 2.7b1, NAMD has its own module for collective  \n");
 fprintf(fplog,"** variable-based calculations including metadynamics, adaptive biasing \n");
 fprintf(fplog,"** force method, umbrella sampling and steered molecular dynamics. \n");
 fprintf(fplog,"** Please, have a look at the NAMD manual for more info. \n\n");
#endif
}; 


//.................................
//.. HERE WE HAVE THE NEW PARSER ..
//.................................

// very long lines allowed
#define PLUMED_LINEMAX  50000

// split a line into words
int PREFIX plumed_get_words(char* line,char*** words){
  char* ww;
  int i;
  (*words)=NULL;
  
  ww=strtok(line," \t\n"); if(ww==NULL) return 0;
  srenew(*words,1);
  (*words)[0]=ww;
  for(i=1;ww=strtok(NULL," \t\n");i++){
    srenew(*words,i+1);
    (*words)[i]=ww;
  };
  return i;
};

void PREFIX plumed_error(const char*s){
  fprintf(stderr,"!!!!! PLUMED ERROR: %s\n",s);
  fprintf(stderr,"!!!!! ABORTING RUN \n");
  if(mtd_data.fplog) {
    fprintf(mtd_data.fplog,"!!!!! PLUMED ERROR: %s\n",s);
    fprintf(mtd_data.fplog,"!!!!! ABORTING RUN \n");
    fflush(mtd_data.fplog);
  };
  EXIT();
};

// parse a word:
//   if the word ends with postfix, delete the postfix and return 1
//   otherwise return 0
// example:
//   char* word; word=malloc(100); strcpy(word,"pippo->");
//   plumed_parse_word(word,"<-"); // returns 0 and does not change word
//   plumed_parse_word(word,"->"); // returns 1 and changes word to "pippo"
int plumed_parse_word(char* word,const char* postfix){
  int lword;
  int lpostfix;
  lword=strlen(word);
  lpostfix=strlen(postfix);
  if(lword<lpostfix) return 0;
  if(strcmp(& word[lword-lpostfix],postfix)) return 0;
  word[lword-lpostfix]=0;
  return 1;
};

int PREFIX plumed_atoi(const char* word){
  int n,i;
  char* buf;
  n=-1;
  sscanf(word,"%i%n",&i,&n);
  if(n!=strlen(word)) {
    snew(buf,strlen(word)+100);
    sprintf(buf,"parsing integer %s\n",word);
    plumed_error(buf);
  }
  return i;
};

double PREFIX plumed_atof(const char* word)
{
 double x;
 int n, found=0; 

 n=-1;
 sscanf(word,"%lf%n",&x,&n);
 if(n!=strlen(word)){
 if(strcmp(word,"pi")==0)    { x=M_PI;      found=1; } 
 if(strcmp(word,"+pi")==0)   { x=M_PI;      found=1; } 
 if(strcmp(word,"-pi")==0)   { x=-1.*M_PI;  found=1; }
 if(strcmp(word,"2pi")==0)   { x=M_2PI;     found=1; }
 if(strcmp(word,"+2pi")==0)  { x=M_2PI;     found=1; }
 if(strcmp(word,"-2pi")==0)  { x=-1.*M_2PI; found=1; }
 if(found==0) plumed_error("Special symbol not recognized. Please, read the manual for accepted symbols.");
 }
 return x;
};
/*
 get group takes a word  and looks if it's in <mygroup> format:
 if it's so it reads the parsed input and looks for a group specified as

 mygroup->       
   4567 678 678 567 4567
   67 89 678 34 234
   5 7
 mygroup-<       

  IT returns the number of found field for the group ( in the case above 12 ) 
  it places the member of the group in  a vector vec , starting from position n
  ( so the new positions will be stored in vec[i] i=n,i<n+j,i++  ) and reallocate the 
  vector if necessary    

  if the word is not in <mygroup> format, it interprets it as a single number and adds it to the atoms list
   
*/

int PREFIX plumed_get_group(const char *word,int **atoms,int n,t_plumed_input* input,FILE *log){
  int lword,justoneatom,foundgroup,nadd;
  int* toadd;
  lword=strlen(word);
  foundgroup=0;
// check for group syntax
  if(lword>2) if(word[0]=='<' && word[lword-1]=='>') foundgroup=1;
// if not group, just read the atomx index;
  if(!foundgroup){
    justoneatom=plumed_atoi(word)-1;
    nadd=1;
    toadd=&justoneatom;
  } else {
// if group, search for it on the list
    int igroup;
    int found;
    char* groupname;
    snew(groupname,strlen(word)-1);
    strncpy(groupname,& word[1],strlen(word)-2);
    found=0;
    for(igroup=0;igroup<input->ngroups;igroup++){
      if(!strcmp(groupname,input->groupnames[igroup])){
        found=1;
        break;
      }
    }
    sfree(groupname);
    if(!found) plumed_error("group not found");
    nadd=input->natoms[igroup];
    toadd=input->atoms[igroup];
  }
  srenew((*atoms),n+nadd);
  int i;
  for(i=0;i<nadd;i++) (*atoms)[i+n]=toadd[i];
  return nadd;
  
};

// routine to parse input file
// * remove comments
// * join lines with continuation
// * find and stores the group definitions
// * save everything which is not a group in the array input->words[iline][iword]
//   where iline=0...(input->nlines-1) and iword=0...(input=->nwords[iline)
// * line numbers are preserved to allow a better error reporting
void PREFIX plumed_read_input(t_plumed_input* input,FILE* file,FILE* log){
  char* line;
  int i;
  int iline,iword;
  char** words;
  int nwords;
  char* inside_group;
  int inside_loop,loop_start,loop_end,loop_stride;

// initial values
  input->nlines=0;
  input->nwords=NULL;
  input->words=NULL;
  input->ngroups=0;
  input->groupnames=NULL;
  input->natoms=NULL;
  input->atoms=NULL;

  inside_group=NULL;
  inside_loop=0;

  snew(line,PLUMED_LINEMAX);

  while(fgets(line,PLUMED_LINEMAX,file)){

   iline=input->nlines;
   input->nlines++;
   srenew(input->nwords,input->nlines);
   input->nwords[iline]=0;
   srenew(input->words,input->nlines);
   input->words[iline]=NULL;

// merge lines ending with "backslash" or "ampersand"
   int linelength;
   linelength=strlen(line);
   if(linelength>1) while(line[linelength-2]=='\\' || line[linelength-2]=='&'){
     if(!fgets(&line[linelength-2],PLUMED_LINEMAX-linelength+2,file))
       plumed_error("last line is not ending");
     linelength=strlen(line);
// append an empty line
// in this way the line count corresponds to the file (better for error reporting)
     iline=input->nlines;
     input->nlines++;
     srenew(input->nwords,input->nlines);
     input->nwords[iline]=0;
     srenew(input->words,input->nlines);
     input->words[iline]=NULL;
// AN EXTRA EMPTY WORD IS ADDED TO BE COMPATIBLE WITH OLDER seek_word
       srenew(input->words[iline],1);
       input->words[iline][0]=NULL;
   }

// Remove comments (beginning with sharp or esclamation)
    for(i=0;line[i];i++) if(line[i]=='#' || line[i]=='!') line[i]=0;

// Split into words:
    nwords=plumed_get_words(line,&words);

// Check for ENDMETA or ENDPLUMED
    if(nwords>0) if(!strcmp(words[0],"ENDMETA") || !strcmp(words[0],"ENDPLUMED")) {
      free(words);
      break;
    }

// loop over all the input words
    for(iword=0;iword<nwords;iword++){

// begin group
      if(plumed_parse_word(words[iword],"->")){
        int igroup;
        if(inside_group) plumed_error("nested groups are not allowed");
        if(iword>0) plumed_error("a group cannot begin in the middle of a line");
        igroup=input->ngroups;
        input->ngroups++;
//   store group name
        srenew(input->groupnames,input->ngroups);
        snew(input->groupnames[igroup],strlen(words[iword])+1);
        strcpy(input->groupnames[igroup],words[iword]);
//   initialize its atom list
        srenew(input->natoms,input->ngroups);
        input->natoms[igroup]=0;
        srenew(input->atoms,input->ngroups);
        input->atoms[igroup]=NULL;
        inside_group=input->groupnames[igroup];

        fprintf(log,"|- GROUP FOUND: %s\n",inside_group);

//   check if other groups with the same name have been defined
        for(i=0;i<igroup;i++) if(!strcmp(inside_group,input->groupnames[i]))
          plumed_error("two groups cannot have the same name");

// end group
      } else if(plumed_parse_word(words[iword],"<-")){
        int igroup,iatom;
        igroup=input->ngroups-1;
        if(!inside_group) plumed_error("end group without begin group");
        if(strcmp(inside_group,words[iword])) plumed_error("end group different from begin group");
        if(inside_loop==1 || inside_loop==2) plumed_error("wrong LOOP syntax");
// this is for backward compatibility with "LOOP 1 10", without stride
        if(inside_loop==3){
          for(i=loop_start;i<=loop_end;i++){
            iatom=input->natoms[igroup];
            input->natoms[igroup]++;
            srenew(input->atoms[igroup],input->natoms[igroup]);
            input->atoms[igroup][iatom]=i-1;
          }
        };

// log the list of members
        fprintf(log,"|- GROUP MEMBERS: ");
        for(i=0;i<input->natoms[igroup];i++){
          if((i+1)%20==0) fprintf(log,"\n|-                ");
          fprintf(log," %i",input->atoms[igroup][i]+1);
        }
        fprintf(log,"\n");
        inside_group=NULL;

// if we are within a group definition, add atom of check for loop syntax
      } else if(inside_group) {
        int i,igroup,iatom;
        igroup=input->ngroups-1;
// NOTE: this should be triggered only by a LOOP keyword inside a group definition
//       it should allow for a hypothetical LOOP keyword in a standard directive
        if(!strcmp("LOOP",words[iword])){
          inside_loop=1;
        } else if(inside_loop==1){
          loop_start=plumed_atoi(words[iword]);
          inside_loop=2;
        } else if(inside_loop==2){
          loop_end=plumed_atoi(words[iword]);
          inside_loop=3;
        } else if(inside_loop==3){
          loop_stride=plumed_atoi(words[iword]);
          for(i=loop_start;i<=loop_end;i+=loop_stride){
            iatom=input->natoms[igroup];
            input->natoms[igroup]++;
            srenew(input->atoms[igroup],input->natoms[igroup]);
            input->atoms[igroup][iatom]=i-1;
          }
          inside_loop=0;
        } else {
//   add a single atom to the list
          i=plumed_atoi(words[iword]);
          iatom=input->natoms[igroup];
          input->natoms[igroup]++;
          srenew(input->atoms[igroup],input->natoms[igroup]);
          input->atoms[igroup][iatom]=i-1;
        }
// if we are on a normal line, just copy the word
      } else {
        int iw;
        iw=input->nwords[iline];
        input->nwords[iline]++;
        srenew(input->words[iline],input->nwords[iline]);
// // AN EXTRA EMPTY WORD IS ADDED TO BE COMPATIBLE WITH OLDER seek_word
          srenew(input->words[iline],input->nwords[iline]+1);
          input->words[iline][iw+1]=NULL;
        snew(input->words[iline][iw],strlen(words[iword])+1);
        strcpy(input->words[iline][iw],words[iword]);
      };
    }

// Finally delete word pointer for this line
    free(words);
  };

// This buffer is not needed anymore
  sfree(line);
  
// DEBUG
//  for(iline=0;iline<input->nlines;iline++) for(iword=0;iword<input->nwords[iline];iword++)
//  fprintf(log,"%i %i : '%s'\n",iline,iword,input->words[iline][iword]);
};

// Deallocate memory
void PREFIX plumed_clear_input(t_plumed_input*input){
  int i,j;
  for(i=0;i<input->nlines;i++) for(j=0;j<input->nwords[i];j++) sfree(input->words[i][j]);
  sfree(input->nwords);
  for(i=0;i<input->nlines;i++) sfree(input->words[i]);
  sfree(input->words);
  for(i=0;i<input->ngroups;i++) sfree(input->groupnames[i]);
  for(i=0;i<input->ngroups;i++) sfree(input->atoms[i]);
  sfree(input->groupnames);
  sfree(input->atoms);
  sfree(input->natoms);
  input->nlines=0;
  input->nwords=NULL;
  input->words=NULL;
  input->ngroups=0;
  input->groupnames=NULL;
  input->natoms=NULL;
  input->atoms=NULL;
}

void PREFIX couple2list( struct coupling_ll **first_elem ,int *at1,int nat1,int *at2,int nat2){
     int i,j,k;
     struct coupling_ll *newelem,mycouple; 
     newelem= (struct coupling_ll *)malloc(sizeof( struct coupling_ll)); 

     printf("SUMMARY************************** %d \n",*first_elem);
     // this add the element to the linked list
     (* newelem).nat1=nat1;   
     (* newelem).at1=(int *)malloc(nat1*sizeof(int));   
     for(i=0;i<nat1;i++){ (* newelem).at1[i]=at1[i];}
     (* newelem).nat2=nat2;   
     (* newelem).at2=(int *)malloc(nat2*sizeof(int));   
     for(i=0;i<nat2;i++){ (* newelem).at2[i]=at2[i];}

     // the address of newelem is pointing to first elem
     (* newelem).next_elem= (*first_elem);
     // now the pointer first elem is pointing to the new elem  
     (* first_elem) = newelem; 
     printf("AT1 "); 
     for(i=0;i<(**first_elem).nat1;i++){
          printf(" %d ",(**first_elem).at1[i]);
     } 
     printf("\n");
     printf("AT2 "); 
     for(i=0;i<(**first_elem).nat2;i++){
         printf(" %d ",(**first_elem).at2[i]);
     } 
     printf("\n");
     printf("ENDSUMMARY************************** %d \n",(*first_elem));
};
void PREFIX freecouple2list( struct coupling_ll **first_elem ){
     struct coupling_ll *newelem; 
        // this add the element to the linked list
         while( (* first_elem)!=NULL){
             free((**first_elem).at1);
             free((**first_elem).at2);
             newelem= (* first_elem);
             (* first_elem)=(* newelem).next_elem;
             free(newelem);
         }
};
void PREFIX scancouple( struct coupling_ll *first_elem ){
     int i,j,k;
     struct coupling_ll *ptr; 
     ptr=first_elem;
     printf("SCANCOUPLE*************************\n");
     while(ptr!=NULL){
        printf("NEWCOUPLE************************* %d\n",ptr);
        printf("NAT1 %d ",(*ptr).nat1);  
        //EXIT(); 
        for(j=0;j<(*ptr).nat1;j++){
           k=(*ptr).at1[j];
           printf(" AT %d ",k);
        }
        printf("\n") ;
        if((*ptr).nat2){
            printf("NAT2 %d ",(*ptr).nat2);  
            for(j=0;j<(*ptr).nat2;j++){
               k=(*ptr).at2[j];
                printf(" AT %d ",k);
            }
        }
        printf("\n");
        ptr=ptr->next_elem;
     }   
     printf("ENDSCANCOUPLE*************************\n");
 
};
void PREFIX setup_projections(struct proj_grad_s *proj ){
    int i,j,k,dimension,ncv;
    int ii,jj,kk,ll,iii,jjj;
    int nat1,*at1;
    int nat2,*at2;
    int *skip; 
    ncv=proj->nlist;
    dimension=mtd_data.natoms; // total number of atoms
    skip=(int *)malloc(dimension*sizeof(int)); 
    at1=(int *)malloc(dimension*sizeof(int)); 
    at2=(int *)malloc(dimension*sizeof(int)); 

    proj->matrix=(real **)malloc(colvar.nconst*sizeof(real *));
    for(i=0;i<colvar.nconst;i++){
       proj->matrix[i]=(real *)malloc(colvar.nconst*sizeof(real));
    }
 
    printf("DIMENSION %d\n",dimension);
    proj->couple=(struct el_couple *)malloc((ncv*(ncv-1)/2)*sizeof( struct el_couple )); //one element for couple         
    for(i=0;i<(ncv*(ncv-1)/2);i++)(proj->couple[i]).first_elem=NULL; // set each pointer for the linked list  to null 
    k=0;// progressive for couple counting
    for(iii=0;iii<proj->nlist-1;iii++){
      i=proj->list[iii]; 
      for(jjj=iii+1;jjj<proj->nlist;jjj++){
          j=proj->list[jjj];
          proj->couple[k].cv1=i;
          proj->couple[k].cv2=j;
          for(ii=0;ii<dimension;ii++){skip[ii]=0;} 
          for(ii=0;ii<colvar.natoms[i];ii++){
               nat1=0;nat2=0;
               kk=colvar.cvatoms[i][ii]; // get the index of atom 
               if(skip[ii]){
                   printf("SKIPPING ATOM %d \n",ii); 
               }else{  
                   // find within the  other set  
                   for(jj=0;jj<colvar.natoms[j];jj++){ // find within the other cv
                           ll=colvar.cvatoms[j][jj]; 
          //                 printf("ATOM1 %d ATOM2 %d \n",kk,ll);
                           if(kk==ll){// found common indexes
                              at2[nat2]=jj;            
                              nat2++;
          //                    printf("ATOM1=ATOM2 \n");
                           }
                   } 
                   if(nat2){ // makes sense only when there is at least one nat2
                        // find within the  same set ( for additive cv ) and exclude computation 
                        for(jj=ii;jj<colvar.natoms[i];jj++){ // find within the other cv
                            ll=colvar.cvatoms[i][jj]; 
                            if(kk==ll){// found common indexes
                               at1[nat1]=jj;            
                               nat1++;
                               skip[jj]=1;
                            }
                        }
                        // transfer it to the right vectors of the linked list 
                   }  
                   if(nat2){
			couple2list( &((proj->couple[k]).first_elem) ,at1,nat1,at2,nat2);
		   }	
               }
          } 
          // freecouple2list(&((proj.couple[k]).first_elem) );
          scancouple((proj->couple[k]).first_elem);
          // increment couple counter

          k++;
      }  
    }

     proj->ncouples=k;
    // printf("NCOUPLES FOUND %d\n",proj.ncouples);
     // DIAGONAL CONTRIBUTION
     printf("DIAGONAL CONTRIBUTION \n"); 
     proj->diagonal=(struct el_diagonal *)malloc(ncv*sizeof( struct el_diagonal )); //one element for couple         
     for(i=0;i<ncv;i++)(proj->diagonal[i]).first_elem=NULL; // set each pointer for the linked list  to null 
     for(i=0;i<ncv;i++){
               for(ii=0;ii<dimension;ii++)skip[ii]=0; 
               // look for all the atoms which  have the same index
               for(ii=0;ii<colvar.natoms[i];ii++){
                   nat1=0;nat2=0;
                   kk=colvar.cvatoms[i][ii]; // get the index of atom 
                   if(skip[ii]){
                      printf("SKIPPING ATOM %d\n",ii); 
                   }else{  
                   // find within the  same set ( for additive cv ) and exclude computation 
                      for(jj=ii;jj<colvar.natoms[i];jj++){ // find within the other cv
                         ll=colvar.cvatoms[i][jj]; 
                         if(kk==ll){// found common indexes
                            at1[nat1]=jj;            
                            nat1++;
                            skip[jj]=1;
                         }
                      }
                   }
                   // transfer it to the right vectors of the linked list 
                   if(nat1)couple2list( &((proj->diagonal[i]).first_elem) ,at1,nat1,at2,nat2);
               }
        //       freecouple2list(&((proj.diagonal[i]).first_elem) );
               scancouple((proj->diagonal[i]).first_elem);
     }
    free(skip);
    free(at1);
    free(at2);

    return; 
}
void PREFIX calc_projections(struct proj_grad_s  *proj ){
    int ncv,cv1,cv2,i,j,k,ii,jj;
    real scal,grad;
    real tmp1x,tmp1y,tmp1z;
    real tmp2x,tmp2y,tmp2z;
    struct coupling_ll *ptr;
    ncv= proj->nlist;
    for(i=0;i<colvar.nconst;i++){
        for(j=i;j<colvar.nconst;j++){
           (proj->matrix[i][j])=0.;
        }
    }
    for(i=0;i<proj->ncouples;i++){// loop over all the couples
            //printf("COUPLE %d ADDR %d\n",i,proj->couple[i].first_elem);
            cv1=proj->couple[i].cv1;
            cv2=proj->couple[i].cv2;
            scal=0.;
            // use the linked list to calculate grad(cv1) dot grad(cv2)
            ptr=proj->couple[i].first_elem;
            while(ptr!=NULL){
               tmp1x=0.;tmp1y=0.;tmp1z=0.; 
               for(j=0;j<(*ptr).nat1;j++){
                  k=(*ptr).at1[j];
                  //printf("KK %d \n",k);
                  tmp1x+=colvar.myder[cv1][k][0];
                  tmp1y+=colvar.myder[cv1][k][1];
                  tmp1z+=colvar.myder[cv1][k][2];
               }
               tmp2x=0.;tmp2y=0.;tmp2z=0.; 
               for(j=0;j<(*ptr).nat2;j++){
                  k=(*ptr).at2[j];
                  //printf("MM %d \n",k);
                  tmp2x+=colvar.myder[cv2][k][0];
                  tmp2y+=colvar.myder[cv2][k][1];
                  tmp2z+=colvar.myder[cv2][k][2];
               }
               scal+=tmp1x*tmp2x+
                     tmp1y*tmp2y+
                     tmp1z*tmp2z;
               ptr=ptr->next_elem;
            }   
            proj->matrix[cv1][cv2]=scal; 
         //   printf("II %d JJ %d VV %f \n",cv1,cv2,scal);
    }
    for(i=0;i<ncv;i++){
               cv1=proj->list[i];
               grad=0.0;
               ptr=proj->diagonal[i].first_elem;
               while(ptr!=NULL){
                   tmp1x=0.;tmp1y=0.;tmp1z=0.; 
                   for(j=0;j<(*ptr).nat1;j++){
                      k=(*ptr).at1[j];
                      //printf("NN %d \n",k);
                      tmp1x+=colvar.myder[cv1][k][0];
                      tmp1y+=colvar.myder[cv1][k][1];
                      tmp1z+=colvar.myder[cv1][k][2];
                   }
                   grad+=tmp1x*tmp1x;
                   grad+=tmp1y*tmp1y;
                   grad+=tmp1z*tmp1z;
                   ptr=ptr->next_elem;
               }
               proj->matrix[cv1][cv1]=grad; 
        //    printf("II %d JJ %d VV %f \n",cv1,cv1,grad);
     }
     for(i=0;i<ncv-1;i++){
        ii=proj->list[i];
        for(j=i+1;j<ncv;j++){
           jj=proj->list[j];
           proj->matrix[jj][ii]=proj->matrix[ii][jj];
        }
     }
 //    printf("MMM \n");
        printf("MMM ");
     for(i=0;i<ncv;i++){
        ii=proj->list[i];
        for(j=0;j<ncv;j++){
           jj=proj->list[j];
           printf(" %f ",proj->matrix[ii][jj]);
        }
     }
        printf("\n");
    return;
}
