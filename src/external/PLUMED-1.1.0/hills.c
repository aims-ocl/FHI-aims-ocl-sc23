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

// This routine calculates, from the variaton of a collective variable, the width
// of the gaussian hills along the CV direction.

void PREFIX hills_adapt()
{
  real fluct, step;
  real Mss02, M2ss0;
  int i_c;

  for(i_c=0;i_c<colvar.nconst;i_c++){
   if(colvar.adapt[i_c].on){ 
// Time for recording fluctuations???
    if((logical.not_same_step) && colvar.it%colvar.adapt[i_c].stride==0 && !firstTime ){
      colvar.Mss0[i_c] += colvar.ss0[i_c];
      colvar.M2ss0[i_c] += colvar.ss0[i_c]*colvar.ss0[i_c];
    }
// Time for evaluating width???    
    if((logical.not_same_step) && colvar.it%colvar.adapt[i_c].block==0 && !firstTime ){
      step = (real) colvar.adapt[i_c].stride / colvar.adapt[i_c].block;      
      M2ss0 = colvar.M2ss0[i_c]*step;
      Mss02 = colvar.Mss0[i_c]*colvar.Mss0[i_c]*step*step;  
      if(M2ss0>Mss02) fluct = sqrt(M2ss0-Mss02);
      else fluct = 0.;
      colvar.delta_r[i_c] = colvar.adapt[i_c].widthmultiplier*fluct;
      if(colvar.delta_r[i_c] > colvar.adapt[i_c].sup) colvar.delta_r[i_c] = colvar.adapt[i_c].sup;
      if(colvar.delta_r[i_c] < colvar.adapt[i_c].inf) colvar.delta_r[i_c] = colvar.adapt[i_c].inf;
      colvar.M2ss0[i_c] = 0.;
      colvar.Mss0[i_c] = 0.;
    }
   }
  }
}

//------------------------------------------------------------------------------------------

void PREFIX hills_push(struct mtd_data_s *mtd_data,real ww, real* ss,real* delta)
{
  int icv,ncv;
  static FILE *file=NULL;
  real inv_ss0[nconst_max];
  int done,wall;

  ncv=colvar.nconst;
  if(hills.n_hills+10>hills.ntothills) hills_reallocate(mtd_data); 

  hills.ww[hills.n_hills] = ww;
  for(icv=0;icv<ncv;icv++) hills.ss0_t[hills.n_hills][icv] = ss[icv];	// new hill center
  for(icv=0;icv<ncv;icv++) colvar.delta_s[hills.n_hills][icv] = delta[icv];       // new hill width

// Multiple walkers lock file
  if(nwalkers>1){				// multiple walkers sincronization 
    done = 0;
    while(!done){				// waiting for LOCFILE to disapper 
      FILE *file2;
      file2 = fopen(locfilen, "r");		// checking LOCFILE 
      if(!file2){                                 
        file2 = fopen((mtd_data->ionode?locfilen:"/dev/null"), "w");
        fclose(file2);
	done = 1;
      }
    }
  }

// PRINT
  if(!file) file = fopen((mtd_data->ionode?hilfilen:"/dev/null"), "a");
   
  fprintf(file, "%10.3f   ", mtd_data->time);
  for(icv=0;icv<ncv;icv++) if(colvar.on[icv]) fprintf(file, "%10.5f   ", hills.ss0_t[hills.n_hills][icv]);
  for(icv=0;icv<ncv;icv++) if(colvar.on[icv]) fprintf(file, "%10.5f   ", colvar.delta_s[hills.n_hills][icv]);
  if(logical.welltemp){
   fprintf(file, "%10.5f   %4.3f \n", hills.ww[hills.n_hills]*colvar.wfactor/(colvar.wfactor-1.0)/mtd_data->eunit,colvar.wfactor);
  } else {
   fprintf(file, "%10.5f   %4.3f \n", hills.ww[hills.n_hills]/mtd_data->eunit,0.0);
  } 
  
  hills.n_hills++;                              // hills added

// flush all the time standalone: not a big deal... 
#ifdef STANDALONE 
  fclose(file);
#endif
  if(nwalkers==1) hills.read=hills.n_hills;

// WALLS
  wall = 0;

  for(icv=0;icv<ncv;icv++) {
    inv_ss0[icv] = colvar.ss0[icv];
    if(logical.ureflect[icv] && colvar.ss0[icv]>(cvw.upper[icv]-colvar.delta_r[icv]) && 
       colvar.ss0[icv]<cvw.upper[icv] && colvar.on[icv]) { inv_ss0[icv] = 2.*cvw.upper[icv]-colvar.ss0[icv]; wall=1;}
    else if(logical.lreflect[icv] && colvar.ss0[icv]<(cvw.lower[icv]+colvar.delta_r[icv]) && 
       colvar.ss0[icv]>cvw.lower[icv] && colvar.on[icv]) { inv_ss0[icv] = 2.*cvw.lower[icv]-colvar.ss0[icv]; wall=1;}
  }

  if(wall) {
    hills.ww[hills.n_hills] = hills.ww[hills.n_hills-1];
    fprintf(file, "%10.3f   ", mtd_data->time);
    for(icv=0;icv<ncv;icv++) {
      hills.ss0_t[hills.n_hills][icv] = inv_ss0[icv];      // new hill center
      if(colvar.on[icv]) fprintf(file, "%10.5f   ", hills.ss0_t[hills.n_hills][icv]);
    }
    for(icv=0;icv<ncv;icv++) {
      colvar.delta_s[hills.n_hills][icv] = colvar.delta_r[icv];       // new hill width
      if(colvar.on[icv]) fprintf(file, "%10.5f   ", colvar.delta_s[hills.n_hills][icv]);
    }
    wall=0;
    fprintf(file, "%10.5f\n", hills.ww[hills.n_hills]);
    hills.n_hills++;                              // hills added
  }

// FLUSH FILE
  fflush(file);
  if(nwalkers>1){
// we close it for multiple walkers calculation, so as to allow other processes to read it
    fclose(file);
    file=NULL;
    remove(locfilen);
  }

}

// This routine add the gaussian hills
void PREFIX hills_add(struct mtd_data_s *mtd_data)

{
  int icv,irep;
  real* all_ww;
  real** all_ss;
  real** all_delta;
  real  this_ww;
  real  this_ss[nconst_max];
  real  this_delta[nconst_max];
  int ineighbour,distance;
  int nrep,ncv;
  int myrep;
  static int last_hill_at_this_step;
  static int first=1;

  if(first) last_hill_at_this_step=colvar.it;
  first=0;

  ncv=colvar.nconst;

  if(hills.max_height>0.0) hills.wwr=hills.rate*(colvar.it-last_hill_at_this_step)*mtd_data->dt;

// ADD THE HILL
  if(logical.welltemp) {
/* since hills are added after the calculation of the bias potential, we can reuse
   the stored Vhills to decide the hills height*/
// custom external boltzman factor 
#ifdef STANDALONE 
    this_ww = hills.wwr*exp(-Vhills/(mtd_data->myboltz*(colvar.wfactor-1.0)*colvar.simtemp));
#else
    this_ww = hills.wwr*exp(-Vhills/(BOLTZ*(colvar.wfactor-1.0)*colvar.simtemp));
#endif
  } else this_ww = hills.wwr;	// new hill height 
  for(icv=0;icv<ncv;icv++) this_ss[icv]    = colvar.ss0[icv];    	// new hill center
  for(icv=0;icv<ncv;icv++) this_delta[icv] = colvar.delta_r[icv];       // new hill width

  if(hills.max_height>0.0 && this_ww<hills.max_height && (colvar.it-last_hill_at_this_step)<hills.max_stride) return;
  last_hill_at_this_step=colvar.it;
  

// add hill to the array
  if(! (logical.remd && colvar.ptmetad_neighbours>0) ) hills_push(mtd_data,this_ww,this_ss,this_delta);

// NEIGHBOURS HILLS for parallel tempering
  if(logical.remd && colvar.ptmetad_neighbours>0){
    nrep=mtd_data->nrepl;
    myrep=mtd_data->repl;
    all_ww=float_1d_array_alloc(nrep);
    all_ss=float_2d_array_alloc(nrep,ncv);
    all_delta=float_2d_array_alloc(nrep,ncv);
    for(irep=0;irep<nrep;irep++)                          all_ww[irep]=0.0;
    for(irep=0;irep<nrep;irep++) for(icv=0;icv<ncv;icv++) all_ss[irep][icv]=0.0;
    for(irep=0;irep<nrep;irep++) for(icv=0;icv<ncv;icv++) all_delta[irep][icv]=0.0;
// first broadcast on the master node of each replica
    if(mtd_data->ionode){
      all_ww[myrep]=this_ww;
      for(icv=0;icv<ncv;icv++) all_ss[myrep][icv]=this_ss[icv];
      for(icv=0;icv<ncv;icv++) all_delta[myrep][icv]=this_delta[icv];
      plumed_intersum(mtd_data,nrep, & all_ww[0]);
      plumed_intersum(mtd_data,nrep*ncv, & all_ss[0][0]);
      plumed_intersum(mtd_data,nrep*ncv, & all_delta[0][0]);
    }
// then broadcast inside each replica
    plumed_sum(mtd_data,nrep, & all_ww[0]);
    plumed_sum(mtd_data,nrep*ncv, & all_ss[0][0]);
    plumed_sum(mtd_data,nrep*ncv, & all_delta[0][0]);
    for(ineighbour=0;ineighbour<nrep;ineighbour++){
      distance=myrep-ineighbour;
      if(distance<0) distance=-distance;
      if(distance<=colvar.ptmetad_neighbours){
        this_ww=all_ww[ineighbour]*exp(-0.5*distance*distance/(colvar.ptmetad_sigma*colvar.ptmetad_sigma));;
        for(icv=0;icv<ncv;icv++)this_ss[icv]=all_ss[ineighbour][icv];
        for(icv=0;icv<ncv;icv++)this_delta[icv]=all_delta[ineighbour][icv];
        hills_push(mtd_data,this_ww,this_ss,this_delta);
      }
    }
    free_1dr_array_alloc(all_ww);
    free_2dr_array_alloc(all_ss,nrep);
    free_2dr_array_alloc(all_delta,nrep);
  }

// After hill addition, update mtd_data.Vbias
  Vhills=hills_engine(colvar.ss0,colvar.ff_hills);
}

//-----------------------------------------------------------------------------------------------

real PREFIX hills_engine_dp(int ih,real* ss0,real* dp){
  int icv,ncv;
  real dp2,diff;
  ncv=colvar.nconst;
  dp2 = 0.;
  for(icv=0;icv<ncv;icv++)if(colvar.on[icv]){
    dp[icv] = (ss0[icv]-hills.ss0_t[ih][icv])/colvar.delta_s[ih][icv];       // (cv(now)-hil(ih))/sigma
    if(colvar.type_s[icv]==5) {
      diff = ss0[icv]-hills.ss0_t[ih][icv];
      if(diff>M_PI) diff-=2.*M_PI;
      else if(diff<-M_PI) diff+=2.*M_PI;
      dp[icv] = diff/colvar.delta_s[ih][icv];
    }
    dp2 += dp[icv]*dp[icv];                  // sum dp*dp
  }
  dp2 = 0.5*dp2;
  return dp2;
};

real PREFIX hills_engine(real* ss0,real* force){
/* this routine calculate hills forces and energy in an arbitrary point.
   let us try to use it anywhere it is needed, so as to avoid errors and double coding.
   in particular, ADD HERE VARIABLES NEEDING PBC
   In GROMACS, and DL_POLY, it is parallel, and should be called with the same arguments by all the processes belongin
   to a replica
*/

  int dp2index;
  real dp2, dp[nconst_max], VhillsLast, diff;
  real Vbias;

  int nh,ncv;
  int ih,icv;

  int npe,rank;

  int nhstart; /* first hill belonging to this process */
  int nhstride; /* stride for hills belonging to this process */

  nh=hills.n_hills;
  ncv=colvar.nconst;

  if(logical.parallel_hills){
    npe=plumed_comm_size(&mtd_data);
    rank=plumed_comm_rank(&mtd_data);
  }else{
    npe=1;
    rank=0;
  };

  if(logical.do_grid){
/*
   when grid are used, it is better NOT to parallelize on hill index.
   in fact, since grid are based on a non-linear spline, splitting the hills among different processors
   leads to slighlty different forces, which means non reproducibility in parallel calculations.
   moreover, the advantage of this parallelization is only when restarting (usually one hill at a time is added)
*/
    nhstart=0;
    nhstride=1;
  }else{
    nhstart=rank;
    nhstride=npe;
  };

  Vbias=0.0;
  if(force) for(icv=0;icv<ncv;icv++) force[icv]=0.0;

// This loop is parallelized
// if logical.debug_grid is set, the actual force is calculated with the hills, but
// in a debug file we write the actual and the grid force and energy
  for(ih=nhstart;ih<nh;ih+=nhstride){
     if(logical.do_grid && ih>=grid.nhills && ih < hills.read)
        grid_addhills(&grid,hills.ww[ih],hills.ss0_t[ih],colvar.delta_s[ih],rank,npe);
     if(!logical.do_grid || logical.debug_grid || ih >= hills.read){
     dp2=hills_engine_dp(ih,ss0,dp);
     if(dp2<DP2CUTOFF){
       dp2index =  dp2*GTAB/DP2CUTOFF;
       VhillsLast = hills.ww[ih]*hills.exp[dp2index];
       Vbias += VhillsLast;                                                            // hills potential energy
       if(force) for(icv=0;icv<ncv;icv++) if(colvar.on[icv])
          force[icv] += dp[icv]/colvar.delta_s[ih][icv]*VhillsLast;	// -dU/dCV 
     }
    }
  }

  if(logical.do_grid) {
    grid.nhills = hills.read;
    if(!logical.debug_grid){
      Vbias+=grid_getstuff(&grid,ss0,force);
    } else {
      static FILE *file=NULL;
      if(!file) {
        char buf[1024];
        sprintf(buf,"DEBUG_GRID%i+%i",mtd_data.repl,plumed_comm_rank(&mtd_data));
        fprintf(stderr,"%s\n",buf);
        file=fopen(buf,"w");
        
      }
      real VbiasG;
      real forceG [nconst_max];
      int i;
      for(i=0;i<ncv;i++)  forceG[i]=0;
      VbiasG=grid_getstuff(&grid,ss0,forceG);
      for(i=0;i<ncv;i++) fprintf(file," %f",ss0[i]);
      fprintf(file,"   %f %f    ",Vbias,VbiasG);
      if(force) for(i=0;i<ncv;i++) fprintf(file," %f %f",force[i],forceG[i]);
      fprintf(file,"\n");
    }
  }

/* when using grid, each process knows the whole grid, so there is no need to reduce */
  if(logical.parallel_hills && ! logical.do_grid){
    if(force) plumed_sum(&mtd_data,ncv, & force[0]);
    plumed_sum(&mtd_data,1, & Vbias);
  }
  
  return Vbias;
};

// This routine calculate the hills contribution
void PREFIX hills_force()
{
  Vhills=hills_engine(colvar.ss0,colvar.ff_hills);
}

//-----------------------------------------------------------------------------------------------

// This routine read the HILLS file
void PREFIX read_hills(struct  mtd_data_s *mtd_data, int restart, int first_read)
{
  double dummy;//,old_factor;
  int i,j,nactive;
  long int line, ix;
  FILE *file;
  char *str, stringa[800];

  if(restart) hills.first_read=0;
  
  file = fopen(hilfilen, "r");									// open

  /****************************** UPDATING *************************/
  if(!restart && !first_read)fsetpos(file,&(hills.line_counter));

  fflush(file);											// open
  if(file==NULL) {										// not exist
          char buf[1024];
          sprintf(buf,"Cannot read HILLS file :: %s",hilfilen);
          plumed_error(buf);
  }

  line = hills.read;

  int active_to_ind[nconst_max];
  nactive=0.; //number of active cvs 
  for(i=0;i<colvar.nconst;i++){
     if(colvar.on[i]){active_to_ind[nactive]=i;nactive++;}
  }

  while(1){											// read cycle
    str = fgets(stringa, 800, file);
    if(str == NULL) break;

    // reallocate if needed during reading
    if(line+10>hills.ntothills) hills_reallocate(mtd_data);  

    j=0;
    str =(char *) strtok(stringa," \t");
    while (str != NULL)
    {
      if( j>0 && j <=nactive  ) { 
          i=active_to_ind[j-1];
          sscanf(str, "%lf", &dummy);
          hills.ss0_t[line][i] = (real) dummy; 
      //    printf("POS %d   %f ",i,  hills.ss0_t[line][i] );
      }
      else if( j>nactive && j<= 2*nactive ) { 
          i=active_to_ind[j-nactive-1];
          sscanf(str, "%lf", &dummy);
          colvar.delta_s[line][i] = (real) dummy;							// read the hills dimension
      //    printf("DELTA %d   %f ",i, colvar.delta_s[line][i]);
      }
      else if(j==2*nactive+1) { 
         sscanf(str, "%lf", &dummy);
         if(logical.welltemp){
            //str =(char *) strtok(NULL, " \t");
            //sscanf(str, "%lf", &old_factor);
            hills.ww[line] = ((real) dummy) * mtd_data->eunit * (colvar.wfactor-1.0) / colvar.wfactor;          		// read the hills height
         } else {
            hills.ww[line] = (real) dummy * mtd_data->eunit;                                                              // read the hills height
         }
      //   printf("WW   %f \n", hills.ww[line]);
      }
      str =(char *) strtok(NULL, " \t");
      j++;
    }
      
    line++;
  }
  
  fgetpos(file,&(hills.line_counter));

  fclose(file);
  hills.n_hills = line;

  if(restart){
    fprintf(mtd_data->fplog, "|- RESTARTING HILLS: # read %li \n",hills.n_hills);
  }else{
    fprintf(mtd_data->fplog, "|- UPDATING HILLS @ %lf fs from %li to %li : TOT %li HILLS read!\n", mtd_data->time,hills.read,hills.n_hills-1,hills.n_hills-hills.read);
  }
  hills.read=hills.n_hills;

//  for(i=hills.read;i<hills.n_hills;i++) fprintf(mtd_data->fplog, "UPDATING # %i HILLS %f %f \n",i,hills.ss0_t[i][0],hills.ss0_t[i][1]);
//  for(i=0;i<hills.n_hills;i++) fprintf(mtd_data->fplog, "AFTER UPDATE # %i HILLS %f %f \n",i,hills.ss0_t[i][0],hills.ss0_t[i][1]);
  fprintf(mtd_data->fplog,"\n");
  fflush(mtd_data->fplog);



}

//------------------------------------------------------------------------------------------

void PREFIX hills_reallocate( struct mtd_data_s *mtd_data) { 
      int i,j,k;
      long int oldtot;
      real  *ww_tmp;
      real  **ss0_t_tmp;
      real  **delta_s_tmp; 

      oldtot=hills.ntothills;
      hills.ntothills+=STACKDIM ;

      fprintf(mtd_data->fplog,"HILLS REALLOCATION--OLD DIMENSION: %li\n",oldtot);
      // save pointers to old arrays
      ww_tmp        = hills.ww;
      ss0_t_tmp     = hills.ss0_t;
      delta_s_tmp   = colvar.delta_s;

      // reallocate to the new stackdim
      hills.ww        = float_1d_array_alloc(hills.ntothills);
      hills.ss0_t     = float_2d_array_alloc(hills.ntothills,colvar.nconst); 
      colvar.delta_s  = float_2d_array_alloc(hills.ntothills,colvar.nconst);   

      // copy old arrays to new ones
      // WE COPY THE FULL ARRAY, NOT JUST hills.n_hills.
      // this is safer, since in reallocation during reading n_hills only set at the end
      for (i=0;i<oldtot;i++){
            hills.ww[i]=ww_tmp[i];  
            for (j=0;j<colvar.nconst;j++) hills.ss0_t[i][j]=ss0_t_tmp[i][j]; 
            for (j=0;j<colvar.nconst;j++) colvar.delta_s[i][j]=delta_s_tmp[i][j]; 
      } 

      // free old arrays
      free_1dr_array_alloc(ww_tmp);
      free_2dr_array_alloc(ss0_t_tmp,oldtot);
      free_2dr_array_alloc(delta_s_tmp,oldtot);

      fprintf(mtd_data->fplog,"HILLS REALLOCATION--NEW DIMENSION: %li\n",hills.ntothills);
}
//-------------------------------------------------------------------------------------------
// Calculate total grid dimension, allocate and initialize
void PREFIX grid_initialize(struct grid_s *grid)
{

 int i, j, ncv;
 
 grid->one2multi = NULL;
 grid->size      = 1;
 ncv             = grid->ncv;
 
 for(i=0;i<ncv;i++) {
  grid->size    *= grid->bin[i]; 
  grid->lbox[i]  = grid->max[i] - grid->min[i];  
  grid->dx[i]    = grid->lbox[i] / grid->bin[i];
 }

// allocation
 grid->pot   = float_1d_array_alloc(grid->size);  
 grid->force = float_2d_array_alloc(grid->size,ncv);

 grid->mem = (ncv+1)*grid->size*sizeof(real)/pow(1024.0,2);
 fprintf(mtd_data.fplog,"|- GRID MEMORY USAGE :: %6.2f MB \n",grid->mem);

 for(j=0;j<grid->size;j++) {
  grid->pot[j] = 0. ;
  for(i=0;i<ncv;i++) grid->force[i][j] = 0. ;
 }
}
//-------------------------------------------------------------------------------------------
// Calculate dimension of the reduced grid. This routine is called once unless
// the HILLS delta is modified during the simulation.
void PREFIX grid_resize_minigrid(struct grid_s *grid, real* delta, real cutoff)
{

 real delta_max;
 int  i, ncv;

// store cutoff for later use
 grid->cutoff   = cutoff;
 ncv            = grid->ncv;
 delta_max      = -100000.;
 grid->minisize =        1;

 for(i=0;i<ncv;i++) if(delta[grid->index[i]]>delta_max) delta_max=delta[grid->index[i]];
// this is actually half of the side of the reduced GRID 
 grid->minilbox = sqrt(2.*cutoff)*delta_max;

 for(i=0;i<ncv;i++) {
  grid->minibin[i] = floor(2.*grid->minilbox/grid->dx[i])+1;
  grid->minisize *= grid->minibin[i];
 } 

 fprintf(mtd_data.fplog,"|- UPDATING REDUCED GRID: SIZE %d pts on %d ",grid->minisize,grid->size);
 fprintf(mtd_data.fplog," DIMENSION "); for(i=0;i<ncv-1;i++) fprintf(mtd_data.fplog," %d x",grid->minibin[i]); 
 fprintf(mtd_data.fplog," %d \n",grid->minibin[ncv-1]);

// create one2multi vector
 grid_create_one2multi(grid);  

}
//-------------------------------------------------------------------------------------------
// allocate and create one2multi vector. This vector is needed
// to work efficiently on the reduced grid. Again this routine is called once unless
// the HILLS delta is modified during the simulation.
void PREFIX grid_create_one2multi(struct grid_s *grid)
{

 int i, j, k, tmpgrid, index1d, ncv;
 real mem;

 ncv   = grid->ncv;
// deallocate if previously allocated
 if(grid->one2multi) free_2di_array_alloc(grid->one2multi,grid->minisize); 

// new allocation
 grid->one2multi = int_2d_array_alloc(grid->minisize,ncv);

 mem = grid->minisize*ncv*sizeof(int)/pow(1024.0,2);
 fprintf(mtd_data.fplog,"|- GRID MEMORY USAGE :: %6.2f MB \n",grid->mem+mem);

// create one2multi
 if(ncv==1) {
  for(i=0;i<grid->minisize;i++) grid->one2multi[i][0] = i;
 } else {
  for(i=0;i<grid->minisize;i++){
   index1d=i;
   for(j=ncv-1;j>=0;j--){
    tmpgrid = 1;
    for(k=0;k<j;k++) tmpgrid*=grid->minibin[k];
    grid->one2multi[i][j] = index1d/tmpgrid;
    index1d = index1d%tmpgrid;
    if(index1d==-1) {
     grid->one2multi[i][j] -= 1;
     index1d=tmpgrid;
    }
   }
  } 
 }

}
//-------------------------------------------------------------------------------------------
// add a hills on the grid. Update potential and force
void PREFIX grid_addhills(struct grid_s *grid, real ww, real* ss, real* delta,int rank,int npe)
{

 int   i, j, ncv, flag; 
 real *xx, *dp, dp2, expo, diff_delta; 
 int  *index_nd, index_1d, dp2index;
 int *index_1d_para;
 real *pot_for_para;

 ncv  = grid->ncv;

// allocate temp array
 xx       = float_1d_array_alloc(ncv);
 dp       = float_1d_array_alloc(ncv);
 index_nd = int_1d_array_alloc(ncv); 

// preliminary checks
// 1) if the HILLS center is inside the grid
// 2) if the GRID bin size is too large
// 3) if delta is changed from previous call
 diff_delta = 0.;
 for(j=0;j<ncv;j++) {
  if((ss[grid->index[j]]<grid->min[j] || ss[grid->index[j]]>=grid->max[j]) && !grid->period[j])
   plumed_error("HILLS outside GRID. Please increase GRID size."); 
  if(grid->dx[j]>delta[grid->index[j]]/2.0) plumed_error("GRID bin size is too large compared to HILLS sigma."); 
  diff_delta += pow((grid->oldelta[j]-delta[grid->index[j]]),2); 
 }
// recalculate the dimension of the reduced grid if delta is changed
 if(sqrt(diff_delta/ncv)>0.01) {
  grid_resize_minigrid(grid,delta,DP2CUTOFF);
  for(j=0;j<ncv;j++) grid->oldelta[j] = delta[grid->index[j]];
 }

// temporary array for parallel computation
 index_1d_para=int_1d_array_alloc(grid->minisize);
 pot_for_para=float_1d_array_alloc(grid->minisize*(1+ncv));
 for(i=0;i<grid->minisize;i++) index_1d_para[i]=0;
 for(i=0;i<grid->minisize*(1+ncv);i++) pot_for_para[i]=0.0;

// add HILL to the points belonging to the reduced GRID
 for(i=rank;i<grid->minisize;i+=npe) {

  index_1d_para[i]=-1; // it means "no force on this point"

  flag=0;
  for(j=0;j<ncv;j++) {
   xx[j] = ss[grid->index[j]] - grid->minilbox + grid->dx[j] * grid->one2multi[i][j];  
   if(grid->period[j]) xx[j] -= grid->lbox[j] * rint(xx[j]/grid->lbox[j]);
   if((xx[j]<grid->min[j] || xx[j]>=grid->max[j]) && !grid->period[j])  flag=1;
   index_nd[j] = floor((xx[j]-grid->min[j])/grid->dx[j]);
  }
  if(flag==1) continue; 

// from multidimensional index to mono
  index_1d=grid_multi2one(grid,index_nd);

// add the gaussian on the GRID
  dp2 = 0.;
  for(j=0;j<ncv;j++) {
   xx[j] = grid->min[j] + grid->dx[j] * index_nd[j];
   dp[j] = xx[j] - ss[grid->index[j]];
   if(grid->period[j]) dp[j] -= grid->lbox[j] * rint(dp[j]/grid->lbox[j]); 
   dp[j] /= delta[grid->index[j]];
   dp2 += dp[j]*dp[j]; 
  }
  dp2 *= 0.5;

  if(dp2<grid->cutoff){
   dp2index =  dp2*GTAB/DP2CUTOFF;
   expo     = ww*hills.exp[dp2index];
   pot_for_para[i*(ncv+1)]=expo;
   for(j=0;j<ncv;j++) pot_for_para[i*(ncv+1)+1+j]=dp[j]/delta[grid->index[j]]*expo;
   index_1d_para[i]=index_1d;
//   grid->pot[index_1d] += expo;
//   for(j=0;j<ncv;j++) grid->force[index_1d][j] += dp[j]/delta[grid->index[j]]*expo;
  }
 
 }

 if(npe>1){ 
   plumed_sum (&mtd_data,grid->minisize*(ncv+1),pot_for_para);
   plumed_sumi(&mtd_data,grid->minisize,index_1d_para);
 }

 for(i=0;i<grid->minisize;i++) {
   if(index_1d_para[i]<0) continue;
   grid->pot[index_1d_para[i]]+=pot_for_para[i*(ncv+1)];
   for(j=0;j<ncv;j++) grid->force[index_1d_para[i]][j] += pot_for_para[i*(ncv+1)+1+j];
 }

// deallocation
 free_1dr_array_alloc(dp);
 free_1dr_array_alloc(xx);
 free_1di_array_alloc(index_nd);

 free_1dr_array_alloc(pot_for_para);
 free_1di_array_alloc(index_1d_para);

}
//-------------------------------------------------------------------------------------------
// from multidimensional index to mono dimensional
int PREFIX grid_multi2one(struct grid_s *grid, int* index_nd)
{
 int i, j, index, tmpgrid;

 index = index_nd[0];
 
 for(i=1;i<grid->ncv;i++) {
  tmpgrid = 1;
  for(j=0;j<i;j++) tmpgrid *= grid->bin[j];
  index += index_nd[i]*tmpgrid;
 }
  
 return index;

}
//-------------------------------------------------------------------------------------------
// this routine returns the bias in a certain point and optionally the forces
real PREFIX grid_getstuff(struct grid_s *grid, real* ss0, real* force)
{

 real *xx;
 int  j, *index_nd, index_1d, ncv;
 real Vbias;

 ncv  = grid->ncv;

// allocate temp array
 xx       = float_1d_array_alloc(ncv);
 index_nd = int_1d_array_alloc(ncv);

// first check if the point is inside the GRID
 for(j=0;j<ncv;j++) {
  xx[j] = ss0[grid->index[j]];
  if((xx[j]<grid->min[j] || xx[j]>=grid->max[j]) && !grid->period[j]) plumed_error("You are outside the GRID!. Please increase GRID size.");
  if(grid->period[j])  xx[j] -= grid->lbox[j] * rint(xx[j]/grid->lbox[j]);
  index_nd[j] = floor((xx[j]-grid->min[j])/grid->dx[j]);
 }

// from multidimensional index to mono
 index_1d=grid_multi2one(grid,index_nd);

 if(!logical.donot_spline){
   real f;
   real where[nconst_max];
   int  stride[nconst_max];
   real der[nconst_max];
   for(j=0;j<ncv;j++) where[j]=xx[j]-grid->min[j]-index_nd[j]*grid->dx[j];
   stride[0]=1;
   for(j=1;j<ncv;j++) stride[j]=stride[j-1]*grid->bin[j-1];
   for(j=0;j<ncv;j++) if(grid->period[j] && index_nd[j]==grid->bin[j]-1) stride[j]*=(1-grid->bin[j]);
   f=spline(ncv,grid->dx,where,& grid->pot[index_1d],&grid->force[index_1d][0],stride,der);
   Vbias=f;
   if(force) for(j=0;j<ncv;j++) force[grid->index[j]] -=der[j];
 } else {
// getting BIAS and (eventually) FORCE
   Vbias = grid->pot[index_1d];
   if(force) for(j=0;j<ncv;j++) force[grid->index[j]] += grid->force[index_1d][j];
 }

// free
 free_1dr_array_alloc(xx);
 free_1di_array_alloc(index_nd);

 return Vbias;
}



// Interpolation with a (sort of) cubic spline.
// The function is built as a sum over the nearest neighbours (i.e. 2 in 1d, 4 in 2d, 8 in 3d,...).
// Each neighbour contributes with a polynomial function which is a product of single-dimensional polynomials,
// written as functions of the distance to the neighbour in units of grid spacing
// Each polynomial is proportional to:
// (1-3x^2+2x^3)  + Q (x-2x^2+x^3)
// * its value and derivative in +1 are zero
// * its value in 0 is 1
// * its derivative in 0 is Q
// so, Q is chosen as the desired derivative at the grid point divided by the value at the grid point
// and the final function is multiplied times the value at the grid point.
//
// It works perfectly, except when the tabulated function is zero (there is a special case).
// Maybe one day I will learn the proper way to do splines...
// Giovanni

real PREFIX spline(int ndim,real *dx,real *where,real *tabf,real *tabder,int* stride,real *der){
// ndim:   dimensionality
// dx:     delta between grid points
// where:  location relative to the floor grid point (always between 0 and dx)
// tabf:   table with function, already pointed at the floor grid point
// tabder: table with minus gradients (the fastest running index is the dimension index), already pointed at the floor grid point
// stride: strides to the next point on the tabf array.
//         note that, in case of PBC, this stride should corrispond to a backward jump of (N-1) points,
//         where N is the number of points in the domain. 
//         also note that the corrisponding strides for tabder can be obtained multipling times ndim
// der:    in output, the minus gradient.

  int idim;
  int npoints,ipoint;
  real X;
  real X2;
  real X3;
  int x0[nconst_max];;
  real fd[nconst_max];
  real C[nconst_max];
  real D[nconst_max];
  int  tmp,shift;
  real f;

  npoints=1; for(idim=0;idim<ndim;idim++) npoints*=2; // npoints=2**ndim

// reset
  f=0;
  for(idim=0;idim<ndim;idim++) der[idim]=0;

// loop over neighbour points:
  for(ipoint=0;ipoint<npoints;ipoint++){

// find coordinate of neighbour point (x0) and shift
    tmp=ipoint;
    shift=0;
    for(idim=0;idim<ndim;idim++){
      x0[idim]=tmp%2; tmp/=2;
      shift+=stride[idim]*x0[idim];
    }
//fprintf(stderr,"%i\n",shift);

// reset contribution from this point:
    real ff;
    ff=1.0;

    for(idim=0;idim<ndim;idim++){
      X=fabs(where[idim]/dx[idim]-x0[idim]);
      X2=X*X;
      X3=X2*X;
      real yy;
      if(fabs(tabf[shift])<0.0000001) yy=0.0;
      else yy=tabder[shift*ndim+idim]/tabf[shift];
                                       // il - e per -derivata
      C[idim]=(1-3*X2+2*X3) - (x0[idim]?-1:1)*yy*(X-2*X2+X3)*dx[idim];
      D[idim]=( -6*X +6*X2) - (x0[idim]?-1:1)*yy*(1-4*X +3*X2)*dx[idim]; // d / dX
      D[idim]*=(x0[idim]?-1:1)/dx[idim]; // chain rule (to where)
      ff*=C[idim];
    }
    for(idim=0;idim<ndim;idim++) {
      int idim1;
      fd[idim]=D[idim];
      for(idim1=0;idim1<ndim;idim1++) if(idim1!=idim) fd[idim]*=C[idim1];
    }
    
    f+=tabf[shift]*ff;
    for(idim=0;idim<ndim;idim++) der[idim]+=tabf[shift]*fd[idim];
  }
  return f;
};




