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

void PREFIX restraint(struct mtd_data_s *mtd_data)
{
  int i_c, ncv, nth, ntrh, ntp, rpxm;             	// indexes for cycles and time

  Vhills = Vwall = 0.;                                 	// Hills and Walls potential energy initialization

  colvar.it=mtd_data->istep;                                                  // update microdynamics step
  if(colvar.it==mtd_data->istep_old){
    logical.not_same_step=0;
  } else {
    logical.not_same_step=1;
    mtd_data->istep_old=colvar.it;
  }

  ncv  = colvar.nconst;                                                  	                // number of CVs
  ntp  = (logical.not_same_step)&&(logical.print)&&(!(colvar.it%colvar.nt_print));			                // have I got to print COLVAR?
  nth  = ( (logical.not_same_step)&&(logical.do_hills)&&(!(colvar.it%hills.nt_hills))&&(!firstTime) )
         || (hills.max_height>0);             // have I got to add HILL?
  ntrh = (logical.not_same_step)&&(logical.do_hills)&&(!(colvar.it%hills.nr_hills))&&(nwalkers>1)&&(!firstTime) ;	// period in steps to read HILLS

  #ifdef GROMACS3
  set_pbc(&metapbc, mtd_data->cell);                                              // set pbc
  #endif
  #ifdef GROMACS4
  set_pbc(&metapbc, mtd_data->ePBC, mtd_data->cell);                                              // set pbc
  #endif

  if(logical.rpxm) rpxm = !(colvar.it%mtd_data->repl_ex_nst) && (logical.not_same_step);                    // have I got a replica exchange trial
  if(logical.debug){						// debug is a GROMACS global which identifies
    test_derivatives(mtd_data);		                       // the use of -debug options in mdrun
    EXIT();			        // allow for only one step of dynamics
  }

// eventually align atoms
  if(colvar.align_atoms){
     int i,iatom1,iatom2;
     real distance[3],dummy;
     for(i=0;i<colvar.align_atoms-1;i++){
       iatom1=colvar.align_list[i];
       iatom2=colvar.align_list[i+1];
       minimal_image(mtd_data->pos[iatom2],mtd_data->pos[iatom1],&dummy,distance);
       mtd_data->pos[iatom2][0]=mtd_data->pos[iatom1][0]+distance[0];
       mtd_data->pos[iatom2][1]=mtd_data->pos[iatom1][1]+distance[1];
       mtd_data->pos[iatom2][2]=mtd_data->pos[iatom1][2]+distance[2];
     }
   };

  // this cycle is intended to calculate CVs values and derivatives
  for(i_c=0;i_c<ncv;i_c++){
    colvar.ff_hills[i_c] = 0.;                 	// initialization hills forces
    cvw.fwall[i_c] = 0.;      			// initialization walls forces

    if((!colvar.on[i_c])&&(logical.rpxm)&&(!rpxm)&&(!ntp)) continue;

    switch(colvar.type_s[i_c]){
      // geometric CVs
      case 1: dist_restraint(i_c, mtd_data); break;			// DISTANCE
      case 2: mindist_restraint(i_c, mtd_data); break;               	// MINDIST
      case 3: coord_restraint(i_c, mtd_data); break;	              	// COORD
      case 4: angle_restraint(i_c, mtd_data); break;	                // ANGLE
      case 5: torsion_restraint(i_c, mtd_data); break;                  // TORSION
      case 6: alfabeta_restraint(i_c, mtd_data); break;                	// ALPHA-BETA
      // interaction CVs
      case 7: hbonds_restraint(i_c, mtd_data); break;                   // HBONDS
      case 8: dipole_restraint(i_c, mtd_data); break;       		// DIPOLE
      // conformations CVs
      case 11: radgyr_restraint(i_c, mtd_data); break;   	       	// RGYR
      case 14: rmsdtor_restraint(i_c, mtd_data); break;	               	// RMSDTOR
      case 16: dihcor_restraint(i_c, mtd_data); break;                 	// DIHEDRAL-COR
      // water CVs
      case 20: waterbridge_restraint(i_c, mtd_data); break;            	// WATERBRIDGE
      // trajectory CVs
      case 30: spath_restraint(i_c, mtd_data); break;                   // S_MAPPATH
      case 31: zpath_restraint(i_c, mtd_data); break;                   // Z_MAPPATH
      case 32: position_restraint(i_c, mtd_data); break;                // ATOM POSITION
      case 33: elstpot_restraint(i_c, mtd_data); break;                 // ELSTPOT POSITION
      case 34: puckering_restraint(i_c, mtd_data); break;               // PUCKERING
    }
  }


  mtd_data->time=colvar.it*(mtd_data->dt);

  if(logical.commitment) commit_analysis();     // commitment analysis

  // this is the really dynamics code in which we calculate hills forces and then forces on CVs.
  if(logical.do_hills){
    hills_force();						// compute hills force and energy
    if(logical.widthadapt) hills_adapt();			// to adapt gaussian width 
    if(nth)  hills_add(mtd_data);	                   	// add HILL
    if(ntrh) {
       read_hills(mtd_data,0,hills.first_read);                   	// is it time to read_hills?
       hills.first_read = 0;
    }
  }

  Vwall=soft_walls_engine(colvar.ss0,cvw.fwall);

  apply_forces(mtd_data);
//  for(i_c=0;i_c<ncv;i_c++) apply_force(i_c, mtd_data);          // add force coming from hills, restraint, walls... 

  if(logical.debug_derivatives) debug_derivatives(mtd_data,ntp);

  if(ntp) print_colvar_enercv(mtd_data->time);	        	// dump COLVAR

  if(firstTime)firstTime = 0;			                // the first PluMed step 

  if(colvar.pg.nlist!=0)calc_projections( &(colvar.pg)); 
}

//------------------------------------------------------------------------------------------

real PREFIX soft_walls_engine(real* ss0,real* force){
  // Wall potential: V_wall(s) = sigma*((s-s0+offset)/redux)**n
  // WARNING: n must be even, sigma>0, redux>0; offset>0; s0 can be positive or negative
  real uscale,lscale,uexp,lexp;
  real V;
  int i;
  V=0.;
  for(i=0;i<colvar.nconst;i++){
    if(force) force[i]=0.0;
    if(logical.upper[i]) {                                                      // if there is a soft wall on this cv
      uscale = (ss0[i]-cvw.upper[i]+cvw.uoff[i])/cvw.ueps[i];                   // calculates the position on the wall 
      if(uscale>0.) {
        uexp = (real) cvw.uexp[i];
        V+=cvw.sigma[i]*pow(uscale, cvw.uexp[i]);
        if(force) force[i]+=-(cvw.sigma[i]/cvw.ueps[i])*uexp*pow(uscale, cvw.uexp[i]-1);
      }
    }
    if(logical.lower[i]) {                                                      // if there is a soft wall on this cv
      lscale = (ss0[i]-cvw.lower[i]-cvw.loff[i])/cvw.leps[i];                   // calculates the position on the wall
      if(lscale<0.) {
        lexp = (real) cvw.lexp[i];
        V+=cvw.lsigma[i]*pow(lscale, cvw.lexp[i]);
        if(force) force[i]+=-(cvw.lsigma[i]/cvw.leps[i])*lexp*pow(lscale, cvw.lexp[i]-1);
      }
    }
  };
  return V;
};
  
//-------------------------------------------------------------------------------------------

void PREFIX steer_cv(int i_c)
{ 
  real force,tmp;
// getting starting point the first time
//  if(colvar.it==0){
  if(firstTime){
   if ( cvsteer.impose_start[i_c] == 0 ) { 
          cvsteer.pos[i_c] = cvsteer.start[i_c] = colvar.ss0[i_c]; } 
   else {    cvsteer.pos[i_c] = cvsteer.start[i_c] ; }   

   fprintf(mtd_data.fplog,"|- STEERING %d CV : STARTVALUE %f\n",i_c+1,cvsteer.start[i_c]);
   fflush(mtd_data.fplog);
 
   if(cvsteer.max[i_c] < cvsteer.start[i_c]){
    cvsteer.sign[i_c] = -1; 
   } else {
    cvsteer.sign[i_c] = +1;
   } 
   // check if you're there since the beginning
   if(cvsteer.pos[i_c]==cvsteer.max[i_c]) {
          fprintf(mtd_data.fplog,"|- STEERING %d CV ARRIVED TO TARGET POINT %f in %d STEPS\n",i_c+1,cvsteer.max[i_c],colvar.it);
          fflush(mtd_data.fplog);
   } 
#ifdef STANDALONE
   FILE *file;
   char filename[100];
   sprintf(filename, "STEER.%d.rst",i_c);
   file = fopen(filename,"w");
   fprintf(file,"%lf  %lf",cvsteer.start[i_c],cvsteer.pos[i_c]);
   fclose(file);
#endif
  } 
  else{ // increase when you're not at the first step 
#ifdef STANDALONE 
        FILE *file;
        char *str, stringa[800],filename[100];
        // open the file 
        sprintf(filename, "STEER.%d.rst",i_c); 
        file = fopen(filename,"r");
        if(file==NULL){
          char buf[1024];
          sprintf(buf,"Cannot read %s  : EXITING\n",filename);
          plumed_error(buf);
        }else{
          str = fgets(stringa, 800, file); 
          sscanf(str, "%lf %lf",&cvsteer.start[i_c],&cvsteer.pos[i_c]);  
          fprintf(mtd_data.fplog,"|- STEERING %d CV RESTARTED FROM POINT: %f  STARTED: %f\n",i_c+1,cvsteer.pos[i_c],cvsteer.start[i_c]);  
          fclose(file);
          if(cvsteer.max[i_c] < cvsteer.start[i_c]){
            cvsteer.sign[i_c] = -1; 
          } else {
            cvsteer.sign[i_c] = +1;
          } 
        }
        if((logical.not_same_step) && fabs(cvsteer.pos[i_c]-cvsteer.start[i_c])<fabs(cvsteer.max[i_c]-cvsteer.start[i_c])) {
                cvsteer.pos[i_c] += cvsteer.sign[i_c] * cvsteer.delta[i_c] / 1000.0;
        }
        sprintf(filename, "STEER.%d.rst",i_c); 
        file = fopen(filename,"w");
        fprintf(file,"%lf  %lf",cvsteer.start[i_c],cvsteer.pos[i_c]);
        fclose(file); 
#else
      if((logical.not_same_step) && fabs(cvsteer.pos[i_c]-cvsteer.start[i_c])<fabs(cvsteer.max[i_c]-cvsteer.start[i_c])) {
              cvsteer.pos[i_c] += cvsteer.sign[i_c] * cvsteer.delta[i_c] / 1000.0;
      }
#endif
  } 
  if(fabs(cvsteer.pos[i_c]-cvsteer.start[i_c])>fabs(cvsteer.max[i_c]-cvsteer.start[i_c])) {
    cvsteer.pos[i_c] = cvsteer.max[i_c];
    fprintf(mtd_data.fplog,"|- STEERING %d CV ARRIVED TO TARGET POINT %f in %d STEPS\n",i_c+1,cvsteer.max[i_c],colvar.it);  
    fflush(mtd_data.fplog); 
  }
      /* HERE PUT THE PERIODICITY YOU NEED!!!!!!!! */ 
  tmp=colvar.ss0[i_c]-cvsteer.pos[i_c];
  if(colvar.type_s[i_c]==5 ){
                   if(tmp > M_PI)
                     tmp -= 2.*M_PI;
                   if(tmp < -M_PI)
                    tmp += 2.*M_PI;
  } 
  force = -cvsteer.spring[i_c]*tmp;
  cvw.fwall[i_c] += force;
  Vwall += -0.5*force*tmp; 
}

//---------------------------------------------------------------------------------------------

void PREFIX print_colvar_enercv(real time_s)
{
    int i_c;
    static FILE *cv_file=NULL;

    if(!cv_file) cv_file = fopen((mtd_data.ionode?colfilen:"/dev/null"), "a");
    fprintf(cv_file, "%10.3f", time_s);
    for(i_c=0;i_c<colvar.nconst;i_c++) fprintf(cv_file, "   %10.5f", colvar.ss0[i_c]);
    fprintf(cv_file, "   %10.3f   %10.3f ", Vhills/mtd_data.eunit, Vwall/mtd_data.eunit);
    for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.steer[i_c]){fprintf(cv_file," RESTRAINT %d %10.5f ",i_c+1,cvsteer.pos[i_c] );} }
    fprintf(cv_file, "\n");
    fflush(cv_file);
#ifdef STANDALONE 
    fclose(cv_file);
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------

void PREFIX commit_analysis()
{
  int i, a, b;
  FILE *commit_file;

  if(firstTime){
    commit_file = fopen((mtd_data.ionode?"COMMIT":"/dev/null"), "a");
    for(i=0;i<colvar.nconst;i++) fprintf(commit_file, "   %14.7f", colvar.ss0[i]);
    fclose(commit_file);
  }

  a = 0;
  b = 0;
  for(i=0;i<colvar.nconst;i++){
    if(colvar.Amin[i]<colvar.ss0[i] && colvar.ss0[i]<colvar.Amax[i]) a++;
    if(colvar.Bmin[i]<colvar.ss0[i] && colvar.ss0[i]<colvar.Bmax[i]) b++;
  }
  if(a==colvar.commit||b==colvar.commit) {
    printf("|- SYSTEM HAS REACHED AN ENDING REGION\n");
    EXIT(); 
    logical.commitment = 0;
    commit_file = fopen((mtd_data.ionode?"COMMIT":"/dev/null"), "a");
    if(a==colvar.commit) fprintf(commit_file, " A\n");
    if(b==colvar.commit) fprintf(commit_file, " B\n");
    fclose(commit_file);
  }
}

//----------------------------------------------------------------------------------------------------------------------
void PREFIX apply_forces(struct mtd_data_s *mtd_data )
{
  real ddr, uscale, lscale, dsdt, nor, fact;
  int ix, i, iat, wall,i_c;
  #ifdef NAMD
  Vector f;   
  #endif 

// set to zero all forces
  for(i=0;i<mtd_data->natoms;i++){
    mtd_data->force[i][0] = 0.0; 
    mtd_data->force[i][1] = 0.0;
    mtd_data->force[i][2] = 0.0;
  }  

  for(i_c=0;i_c<colvar.nconst;i_c++){

    if(logical.steer[i_c]) steer_cv(i_c);			// do you want to pull this cv with a quadratic potential?

    ddr = colvar.ff_hills[i_c] + cvw.fwall[i_c];		// soft wall and hills contribution
    for(i=0;i<colvar.natoms[i_c];i++) {
      iat = colvar.cvatoms[i_c][i];
#ifdef NAMD 
       f.x=ddr*colvar.myder[i_c][i][0];
       f.y=ddr*colvar.myder[i_c][i][1];
       f.z=ddr*colvar.myder[i_c][i][2];
       addForce(iat,f);
#else
      mtd_data->force[iat][0] += ddr*colvar.myder[i_c][i][0];     // PluMeD forces
      mtd_data->force[iat][1] += ddr*colvar.myder[i_c][i][1];
      mtd_data->force[iat][2] += ddr*colvar.myder[i_c][i][2];
#endif    
    }
  }
}

