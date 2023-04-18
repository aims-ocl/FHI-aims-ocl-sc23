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
#define EXTERNALS 1
#include "metadyn.h"
#include <assert.h>

#if defined (PLUMED_GROMACS)
void mtd_data_init (int ePBC, real *charge, real *mass, 
                    int natoms, real dt, int repl_ex_nst, int repl, 
                    int nrepl, real rte0, real rteio, const t_commrec *mcr, FILE *fplog)
{
// this interface should be consistent with that in GROMACS/mdlib/sim_utils.c
   extern void (*plumed_do_forces)(int,int,real (*pos)[3], real (*force)[3], real box[3][3]);
// set the pointer in do_forces() to the metadynamics calculator
   plumed_do_forces=meta_force_calculation;

   mtd_data.mcr = mcr;
   mtd_data.natoms = natoms;
   mtd_data.pos = float_2d_array_alloc(mtd_data.natoms,3);
   mtd_data.force = float_2d_array_alloc(mtd_data.natoms,3);
   mtd_data.vel = NULL;
   mtd_data.mass = float_1d_array_alloc(mtd_data.natoms);
   mtd_data.charge = float_1d_array_alloc(mtd_data.natoms);

   int iat,iat_dd;
   for(iat=0;iat<mtd_data.natoms;iat++){
     iat_dd=plumed_dd_index(iat);
     if(iat_dd>=0){
       mtd_data.mass[iat]=mass[iat_dd];
       mtd_data.charge[iat]=charge[iat_dd];
     } else {
       mtd_data.mass[iat]=0.0;
       mtd_data.charge[iat]=0.0;
     }
   };
#ifdef GMX_MPI
#ifdef GROMACS4
   if (DOMAINDECOMP(mtd_data.mcr)) {
// for dd, collect mass and charge
// for pd, atoms are already shared and this is not necessary
     plumed_sum(&mtd_data,mtd_data.natoms,&mtd_data.mass[0]);
     plumed_sum(&mtd_data,mtd_data.natoms,&mtd_data.charge[0]);
   }
#endif
#endif

   mtd_data.repl_ex_nst = repl_ex_nst;
   mtd_data.repl = repl;
   mtd_data.nrepl = nrepl;
   mtd_data.rte0 = rte0;
   mtd_data.rteio = rteio;
   mtd_data.dt = dt;
   mtd_data.fplog = fplog;
   mtd_data.eunit = 1.;
#if defined (GROMACS4)
   mtd_data.ePBC = ePBC;
#endif
   mtd_data.istep_old = -1;
   if(!mtd_data.ionode) mtd_data.fplog=fopen("/dev/null","w");
   logical.not_same_step=1;
   sprintf(hills.dir, ".");
}

void plumed_setstep(int istep){
  mtd_data.istep=istep;
};

// Routine which finds the local index from the global index, for gromacs domain decomposition
int plumed_dd_index(int iat){
   int cell_dd,iat_dd;
   iat_dd=iat;
#ifdef GMX_MPI
#if defined GROMACS4
   if (DOMAINDECOMP(mtd_data.mcr)) {
     cell_dd=mtd_data.mcr->dd->ga2la[iat].cell;
     if(cell_dd!=0) iat_dd=-1;
     else iat_dd=mtd_data.mcr->dd->ga2la[iat].a;
   };
#endif
#endif
   return iat_dd;
};

void meta_force_calculation(int start,int homenr,real (*pos)[3], real (*force)[3], real box[3][3])
{
   int i;
   int kk;
   int iat,iat_dd,cell_dd;
   int icv,ncv;
   rvec *pos_cvatoms;
   int nreqatoms;

   copy_mat(box, mtd_data.cell);

   ncv=colvar.nconst;

   nreqatoms=0;
   for(icv=0;icv<ncv;icv++) nreqatoms+=colvar.natoms[icv];

   nreqatoms+=colvar.align_atoms;

// compact array with only positions of atoms involved in metadynamics
   snew(pos_cvatoms,nreqatoms);

// copy relevent atoms on compact array
// in dd, each process owns part of the atoms: iat_dd provides the mapping, negative means non-local atom
   kk=0;
   for(icv=0;icv<ncv;icv++) for(i=0;i<colvar.natoms[icv];i++){
     iat_dd=plumed_dd_index(colvar.cvatoms[icv][i]);
     if(iat_dd>=start && iat_dd<start+homenr) {
       pos_cvatoms[kk][0]=pos[iat_dd][0]; pos_cvatoms[kk][1]=pos[iat_dd][1]; pos_cvatoms[kk][2]=pos[iat_dd][2];
     }else{
       pos_cvatoms[kk][0]=0.0; pos_cvatoms[kk][1]=0.0; pos_cvatoms[kk][2]=0.0;
     }
     kk++;
   }
// if an atom is aligned, we also copy it here
// this is not efficient (if you align all the atoms involved in the CVs, you will send around
// twice the data as necessary), but still should work
// (in the future, we have to implement a list of uniq atoms required by plumed, obtained as the
// union of the sets of atoms required by CV evaluation and alignment)
   for(i=0;i<colvar.align_atoms;i++){
     iat_dd=plumed_dd_index(colvar.align_list[i]);
     if(iat_dd>=start && iat_dd<start+homenr) {
       pos_cvatoms[kk][0]=pos[iat_dd][0]; pos_cvatoms[kk][1]=pos[iat_dd][1]; pos_cvatoms[kk][2]=pos[iat_dd][2];
     }else{
       pos_cvatoms[kk][0]=0.0; pos_cvatoms[kk][1]=0.0; pos_cvatoms[kk][2]=0.0;
     }
     kk++;
   }

// just to double check
   assert(kk==nreqatoms);

// collect positions
// (needed for both particle and domain decomposition)
   plumed_sum(&mtd_data,3*nreqatoms,&pos_cvatoms[0][0]);

// this is just to avoid having uninitialized elements. in practice they should not be touched
   for(i=0;i<mtd_data.natoms;i++){
     mtd_data.pos[i][0] = 0.0; mtd_data.pos[i][1] = 0.0; mtd_data.pos[i][2] = 0.0;
   }

// scatter the relevant atoms on the large mtd_data.pos array.
// ( in the future we may consider eliminating the mtd_data.pos array and working directly with pos_cvatoms,
//   which should be more efficient but will require hacking all the restraint routines.
//   similarly, a for_cvatoms could be used for forces.)
   kk=0;
   for(icv=0;icv<ncv;icv++) for(i=0;i<colvar.natoms[icv];i++){
     iat=colvar.cvatoms[icv][i];
     mtd_data.pos[iat][0]=pos_cvatoms[kk][0]; mtd_data.pos[iat][1]=pos_cvatoms[kk][1]; mtd_data.pos[iat][2]=pos_cvatoms[kk][2];
     kk++;
   }
   for(i=0;i<colvar.align_atoms;i++){
     iat=colvar.align_list[i];
     mtd_data.pos[iat][0]=pos_cvatoms[kk][0]; mtd_data.pos[iat][1]=pos_cvatoms[kk][1]; mtd_data.pos[iat][2]=pos_cvatoms[kk][2];
     kk++;
   }

// double check
   assert(kk==nreqatoms);

   restraint(&mtd_data);

// add the contribution to the gromacs forces.
// each process only updates the local atoms
   for(iat=0;iat<mtd_data.natoms;iat++){
     iat_dd=plumed_dd_index(iat);
     if(iat_dd>=start && iat_dd<start+homenr){
       force[iat_dd][0] += mtd_data.force[iat][0];
       force[iat_dd][1] += mtd_data.force[iat][1];
       force[iat_dd][2] += mtd_data.force[iat][2];
     }
   }
 
   sfree(pos_cvatoms);
}
#elif DRIVER
void mtd_data_init(int atoms, real *mass, real *charge, char *metainp, int pbc, real *box)
{
 int i;
 mtd_data.pos       = float_2d_array_alloc(atoms,3);
 mtd_data.vel       = float_2d_array_alloc(atoms,3);
 mtd_data.force     = float_2d_array_alloc(atoms,3);
 mtd_data.charge    = (real *)calloc(atoms,sizeof(real));
 mtd_data.mass      = (real *)calloc(atoms,sizeof(real));
 mtd_data.natoms    = atoms;
 mtd_data.dt        = 1;
 mtd_data.istep_old = -1;
 mtd_data.istep     = 0;  
 mtd_data.fplog     = stdout;
 mtd_data.imcon     = pbc;
 mtd_data.eunit     = 1.;
 mtd_data.ionode = 1;
 logical.not_same_step=1;
 if(mtd_data.imcon==1){
  mtd_data.cell[0]=box[0];
  mtd_data.cell[1]=box[1];
  mtd_data.cell[2]=box[2]; 
 } 

 for(i=0;i<mtd_data.natoms;i++) {
   mtd_data.mass[i]   = mass[i];
   mtd_data.charge[i] = charge[i];
 }
 strcpy(mtd_data.metaFilename, metainp);
}

void cv_calculation_(real *box, real *pos, int *ncv, real *cv)
{
 int i;
 for(i=0;i<mtd_data.natoms;i++){
 mtd_data.pos[i][0] = pos[i];
 mtd_data.pos[i][1] = pos[i + mtd_data.natoms];
 mtd_data.pos[i][2] = pos[i + 2*mtd_data.natoms];
 }

 // update cell info
 if(mtd_data.imcon==1){
  mtd_data.cell[0]=box[0];
  mtd_data.cell[1]=box[1];
  mtd_data.cell[2]=box[2];
 } 


 restraint(&mtd_data);
 // step increment
 mtd_data.istep++;

 for(i=0;i<*ncv;i++) cv[i] = colvar.ss0[i];
}
#elif STANDALONE  
void mtd_data_init(int atoms, real *mass, real *charge, int pbc, real *box, real *tstep, int *nstep, real *myboltz, real *ampli,  char *metainp )
{
 int i;
 mtd_data.pos       = float_2d_array_alloc(atoms,3);
 mtd_data.vel       = float_2d_array_alloc(atoms,3);
 mtd_data.force     = float_2d_array_alloc(atoms,3);
 mtd_data.charge    = (real *)calloc(atoms,sizeof(real));
 mtd_data.mass      = (real *)calloc(atoms,sizeof(real));
 mtd_data.natoms    = atoms;
 mtd_data.dt        = (*tstep);
 mtd_data.istep_old = -1;
 mtd_data.istep     = (*nstep);  
 mtd_data.ampli     = (*ampli);  
 mtd_data.myboltz   = (*myboltz);  
 mtd_data.imcon     = pbc;
 mtd_data.eunit     = 1.;
 mtd_data.ionode = 1;

 sprintf(hills.dir, ".");

 logical.not_same_step=1;
 if(mtd_data.imcon==1){
  mtd_data.cell[0]=box[0];
  mtd_data.cell[1]=box[1];
  mtd_data.cell[2]=box[2]; 
 } 
 if(mtd_data.istep==1){
    mtd_data.fplog = fopen("PLUMED.OUT","w");
 }else {
    //mtd_data.fplog = fopen("PLUMED.OUT","a");
    mtd_data.fplog = fopen("/dev/null","w");
 }
 for(i=0;i<mtd_data.natoms;i++) {
   mtd_data.mass[i]   = mass[i];
   mtd_data.charge[i] = charge[i];
 }
 strcpy(mtd_data.metaFilename, metainp);
}

void cv_calculation_standalone_(real *box, real *pos, real *force , real *ene)
{
 int i;
 for(i=0;i<mtd_data.natoms;i++){
 mtd_data.pos[i][0] = pos[i];
 mtd_data.pos[i][1] = pos[i + mtd_data.natoms];
 mtd_data.pos[i][2] = pos[i + 2*mtd_data.natoms];
 }

 // update cell info
 if(mtd_data.imcon==1){
  mtd_data.cell[0]=box[0];
  mtd_data.cell[1]=box[1];
  mtd_data.cell[2]=box[2];
 } 


 restraint(&mtd_data);
 // step increment
 mtd_data.istep++;

   for(i=0;i<mtd_data.natoms;i++){
       force[i                    ] = mtd_data.force[i][0];
       force[i + mtd_data.natoms  ] = mtd_data.force[i][1];
       force[i + 2*mtd_data.natoms] = mtd_data.force[i][2];
   }
   (*ene)=Vhills/mtd_data.eunit+Vwall/mtd_data.eunit;
}

#elif ACEMD
void mtd_data_init( real *charge, real *mass,
                    int natoms, real dt, int repl,
                    int nrepl, real rte0, real rteio, char *metainp, real box[3])
{
 int i;
 mtd_data.pos       = float_2d_array_alloc(natoms,3); 
 mtd_data.vel       = float_2d_array_alloc(natoms,3);
 mtd_data.force     = float_2d_array_alloc(natoms,3);
 mtd_data.charge    = charge;
 mtd_data.mass      = mass;
 mtd_data.natoms    = natoms;
 mtd_data.repl      = repl; 
 mtd_data.nrepl     = nrepl;
 mtd_data.rte0      = rte0;
 mtd_data.rteio     = rteio;
 mtd_data.dt        = dt; 
 mtd_data.istep_old = -1;
 mtd_data.istep     = 0; 
 mtd_data.imcon     = 1; 
 mtd_data.eunit     = 1.;
 mtd_data.ionode = 1;
 logical.not_same_step=1;
 sprintf(hills.dir, ".");
 mtd_data.fplog = fopen("log.file","w");
 strcpy(mtd_data.metaFilename, metainp); 
   mtd_data.cell[0]=box[0];
   mtd_data.cell[1]=box[1];
   mtd_data.cell[2]=box[2];

}

void meta_force_calculation(struct aceplug_sim_t* s )
{
   int i;


   // step increment
   mtd_data.istep++; 

   for(i=0;i<mtd_data.natoms;i++){
     mtd_data.pos[i][0] = (double) s->pos[i].x;
     mtd_data.pos[i][1] = (double) s->pos[i].y;
     mtd_data.pos[i][2] = (double) s->pos[i].z;
     //     mtd_data.vel[i][0] = (double) s->vel[i].x;
     //mtd_data.vel[i][1] = (double) s->vel[i].y;
     //mtd_data.vel[i][2] = (double) s->vel[i].z;
   }

   restraint(&mtd_data);

   for(i=0;i<mtd_data.natoms;i++){
     s->frc[i].x= s->frc[i].x+mtd_data.force[i][0];
     s->frc[i].y= s->frc[i].y+mtd_data.force[i][1];
     s->frc[i].z= s->frc[i].z+mtd_data.force[i][2];
   }
}
#elif OPEP
void mtd_data_init(int pbc, real tstep,int atoms, int repl, int nrepl, real rte0, real rteio, real *mass, char *lpath, char *logfile, char *metainp)
{
 int i;
 mtd_data.pos       = float_2d_array_alloc(atoms,3); 
 mtd_data.vel       = float_2d_array_alloc(atoms,3);
 mtd_data.force     = float_2d_array_alloc(atoms,3);
 mtd_data.charge    = (real *)calloc(atoms,sizeof(real));
 mtd_data.mass      = (real *)calloc(atoms,sizeof(real));
 mtd_data.natoms    = atoms;
 mtd_data.repl      = repl; 
 mtd_data.nrepl     = nrepl;
 mtd_data.rte0      = rte0;
 mtd_data.rteio     = rteio;
 mtd_data.dt        = tstep; 
 mtd_data.istep_old = -1;
 mtd_data.istep     = 0; 
 mtd_data.imcon     = pbc; 
 mtd_data.eunit     = 1.;
 mtd_data.ionode = 1;
#ifdef MPI
   mtd_data.comm=MPI_COMM_SELF;
   if(nrepl==1){
     mtd_data.intercomm=MPI_COMM_NULL;
   } else {
     mtd_data.intercomm=MPI_COMM_WORLD;
   };
#endif
 logical.not_same_step=1;
 for(i=0;i<mtd_data.natoms;i++) mtd_data.mass[i] = mass[i];

 strcpy(mtd_data.metaFilename, metainp); 
 if(nrepl==1){
   strcpy(hills.dir, ".");
 }else{
  strcpy(hills.dir, lpath);
 }
 strcpy(mtd_data.log, logfile); 
 mtd_data.fplog = fopen(mtd_data.log,"a"); 

}

void meta_force_calculation_(real *pos, real *force)
{
 int i;
 for(i=0;i<mtd_data.natoms;i++){
 mtd_data.pos[i][0] = pos[i];
 mtd_data.pos[i][1] = pos[i + mtd_data.natoms];
 mtd_data.pos[i][2] = pos[i + 2*mtd_data.natoms];
 }

 restraint(&mtd_data);
 // step increment
 mtd_data.istep++; 
 
 for(i=0;i<mtd_data.natoms;i++){
 force[i]                     = mtd_data.force[i][0];
 force[i + mtd_data.natoms]   = mtd_data.force[i][1];
 force[i + 2*mtd_data.natoms] = mtd_data.force[i][2];
 }
}

void share_bias_(int *rep, real *Vbias, real *Vbiasx) 
{
ptmetad_vbias(*rep,Vbias,Vbiasx);
}

void switch_fluct_(int *rep)
{
 if(logical.widthadapt) ptmetad_exchfluct(*rep);
}

void bias_exchange_(int *nrep, int *biaseed, int *ind) 
{
if(logical.rpxm) bias_exchange_traj(*nrep,biaseed,ind); 
}
#elif AMBER 
void mtd_data_init(int atoms, real dt ,real *mass, real *charge, char *metainp)
{
 int i ;
 mtd_data.pos=float_2d_array_alloc(atoms,3); 
 mtd_data.vel=float_2d_array_alloc(atoms,3);
 mtd_data.force=float_2d_array_alloc(atoms,3);
 mtd_data.charge=(real *)calloc(atoms,sizeof(real));
 for (i=0;i<atoms;i++){mtd_data.charge[i]= charge[i]/18.2223; } // amber conversion 
 mtd_data.mass=(real *)calloc(atoms,sizeof(real));
 for (i=0;i<atoms;i++){mtd_data.mass[i]=mass[i]; } 
 mtd_data.fplog = fopen("log.file","w"); 
 mtd_data.natoms = atoms;
 mtd_data.dt = dt*1000. ; // timestep from ps to fs  
 mtd_data.repl  = -1;
 mtd_data.eunit = 1.;
 mtd_data.ionode = 1;
 strcpy(mtd_data.metaFilename,metainp); 
 sprintf(hills.dir, ".");
 logical.not_same_step=1;
 mtd_data.istep_old = -1;
 mtd_data.istep = 0;
}
void meta_force_calculation__(real *box, real *rrr,real *fff, int *nn)
{
 int i;
 mtd_data.istep = (*nn);
 //cell ??
 if((box[0]==0.)&&(box[0]==box[1])&&(box[0]==box[2])) { 
    //no cell
     mtd_data.imcon=0;
 }else{
     mtd_data.imcon=1;
     mtd_data.cell[0]=box[0];
     mtd_data.cell[1]=box[1];
     mtd_data.cell[2]=box[2];
 }
 for (i=0;i<mtd_data.natoms;i++){
      mtd_data.pos[i][0]=rrr[i*3+0];
      mtd_data.pos[i][1]=rrr[i*3+1];
      mtd_data.pos[i][2]=rrr[i*3+2];
 } 
 restraint(&mtd_data);
 
 for(i=0;i<mtd_data.natoms;i++){
   // first convert in kcal
   fff[i*3+0]   += mtd_data.force[i][0]; // conversion from kcal to internal units
   fff[i*3+1]   += mtd_data.force[i][1];
   fff[i*3+2]   += mtd_data.force[i][2];
 }

}
#elif FHIAIMS
void mtd_data_init(int atoms, real dt ,real *mass, char *metainp)
{
 int i ;
 for(i=0;i<atoms;i++) {printf("%lf\t",mass[i]);} printf("\n");
 mtd_data.pos=float_2d_array_alloc(atoms,3); 
 mtd_data.vel=float_2d_array_alloc(atoms,3); // where are these from???
 mtd_data.force=float_2d_array_alloc(atoms,3);
 mtd_data.charge=(real *)calloc(atoms,sizeof(real));
 for (i=0;i<atoms;i++){mtd_data.charge[i]= 0.;} // no charges in a DFT code
 mtd_data.mass=(real *)calloc(atoms,sizeof(real));
 //printf("n atoms in meta %d\n",atoms);
 for (i=0;i<atoms;i++){mtd_data.mass[i]=mass[i];} //printf("masses in meta %lf\n",mass[i]);}
  
 mtd_data.fplog = fopen("log.file","w"); 
 mtd_data.natoms = atoms;
 mtd_data.dt = dt*1000. ; // timestep from ps to fs  
 mtd_data.repl  = -1;
 mtd_data.eunit = 1.;
 mtd_data.ionode = 1;
 strcpy(mtd_data.metaFilename,metainp); 
 sprintf(hills.dir, ".");
 logical.not_same_step=1;
 mtd_data.istep_old = -1;
 mtd_data.istep = 0;
}

void meta_force_calculation__(int *npbc, real *cell, real *rrr, real *fff, int *nn)
{
 int i,j;
 mtd_data.istep = (*nn);
 //cell ??
 printf("npbc after %d\n",*npbc);
 printf("cell vectors %lf %lf %lf\n %lf %lf %lf\n %lf %lf %lf\n", cell[0], cell[1],cell[2],cell[3],cell[4],cell[5],cell[6],cell[7],cell[8]);
 if(*npbc==0) { 
    //no cell
     mtd_data.imcon=0;
 }else if(*npbc==3) {
   for(i=0;i<=2;i++) {
     for(j=0;j<=2;j++) {
       if(i != j) {
         if(cell[3*i+j] != 0.) {
           printf("Non orthogonal cell not supported, yet\n");
           EXIT();
	 }
       }
     }	
   }
   mtd_data.imcon=1;
   mtd_data.cell[0]=cell[3*1+0];
   mtd_data.cell[1]=cell[3*2+1];
   mtd_data.cell[2]=cell[3*3+2];
 }else{
	printf("Slab (2D) and wire (1D) pbc not supported, yet\n");
	EXIT();
 }
 for (i=0;i<mtd_data.natoms;i++){
      mtd_data.pos[i][0]=rrr[i*3+0];
      mtd_data.pos[i][1]=rrr[i*3+1];
      mtd_data.pos[i][2]=rrr[i*3+2];
      printf("plumed pos %lf\t%lf\t%lf\n",rrr[i*3+0],rrr[i*3+1],rrr[i*3+2]);
 } 
 restraint(&mtd_data);
 
 for(i=0;i<mtd_data.natoms;i++){
   // first convert in kcal
   printf("plumed forces before %e\t%e\t%e\n",fff[i*3+0],fff[i*3+1],fff[i*3+2]);
   fff[i*3+0]   += mtd_data.force[i][0]; // conversion from kcal to internal units
   fff[i*3+1]   += mtd_data.force[i][1];
   fff[i*3+2]   += mtd_data.force[i][2];
   printf("plumed forces after %e\t%e\t%e\n",fff[i*3+0],fff[i*3+1],fff[i*3+2]);
 }
 	
}

/*	Below the entry points for IBM Fortran (without "_") */

void meta_force_calculation_(int *npbc, real *cell, real *rrr, real *fff, int *nn)
{
	meta_force_calculation__(npbc, cell, rrr, fff, nn);
}

void init_metadyn_(int *atoms, real *dt, real *mass, char *metainp, int ll)
{
	init_metadyn__(atoms, dt, mass, metainp, ll);
}

#elif DL_POLY
void mtd_data_init(int atoms, real dt ,real *mass, real *charge, int *imcon, real *eunit, char *metainp)
{
 int i ;
 mtd_data.pos=float_2d_array_alloc(atoms,3);
 mtd_data.vel=float_2d_array_alloc(atoms,3);
 mtd_data.force=float_2d_array_alloc(atoms,3);
 mtd_data.charge=(real *)calloc(atoms,sizeof(real));
 for (i=0;i<atoms;i++){mtd_data.charge[i]= charge[i]; }
 mtd_data.mass=(real *)calloc(atoms,sizeof(real));
 for (i=0;i<atoms;i++){mtd_data.mass[i]=mass[i]; }
 mtd_data.fplog = fopen("PLUMED.OUT","a");
 mtd_data.natoms = atoms;
 mtd_data.dt = dt ; // the timestep is kept in ps
 mtd_data.repl = -1;
 mtd_data.istep_old = -1;
 mtd_data.eunit=(*eunit);
 mtd_data.imcon=(* imcon);
#ifdef MPI
 mtd_data.intercomm=MPI_COMM_NULL;
 mtd_data.comm=MPI_COMM_WORLD;
#endif
 mtd_data.ionode = (plumed_comm_rank(&mtd_data)==0);
 strcpy(mtd_data.metaFilename,metainp);
 sprintf(hills.dir, ".");
 logical.not_same_step=1;
}
void meta_force_calculation_(real *cell, int *istep, real *xxx, real *yyy, real *zzz, real *fxx, real *fyy, real *fzz)
{
 int i;
 // cell type 
 // cell size 
 for(i=0;i<9;i++)mtd_data.cell[i]=cell[i]; 
 mtd_data.istep=(* istep); // <---paolo
 for(i=0;i<mtd_data.natoms;i++){
   mtd_data.pos[i][0] = xxx[i];
   mtd_data.pos[i][1] = yyy[i];
   mtd_data.pos[i][2] = zzz[i];
   //printf("POS %d %f %f %f\n",i, mtd_data.pos[i][0], mtd_data.pos[i][1], mtd_data.pos[i][2]);
 }

 restraint(&mtd_data);

// in DL_POLY, meta_force_calculation is called before the global sum of forces.
// thus, to avoid multiple counting, only one process is adding the force
 if(plumed_comm_rank(&mtd_data)==0) {
   for(i=0;i<mtd_data.natoms;i++){
     fxx[i]   += mtd_data.force[i][0]; 
     fyy[i]   += mtd_data.force[i][1];
     fzz[i]   += mtd_data.force[i][2];
   //printf("FORCE %d %f %f %f\n",i, mtd_data.force[i][0], mtd_data.force[i][1], mtd_data.force[i][2]);
   }
 }

}
#elif NAMD

PREFIX GlobalMasterMetaDynamics() :
    GlobalMasterEasy("metaDynamicsScript")
{
// Initialize subclass and get current config
    easy_init(config);
}

// easy_init and easy_calc are virtuals of GlobalMasterEasy: must be defined anyway
void PREFIX easy_init(const char *config) {
  int i,j,i_c;

// PluMeD initialization 
  init_metadyn();

// requesting atom to NAMD
    for(i_c=0;i_c<colvar.nconst;i_c++){
            for(i=0;i<colvar.natoms[i_c];i++){
                    j=colvar.cvatoms[i_c][i];
                    if(requestAtom(j)){
                     printf("ATOM ID %i DOESN'T EXIST\n",j);
                     EXIT(); 
                    }
            }
    }
}

void PREFIX  easy_calc() {
    Vector coord_from;
    int i,j,i_c;

// getting atom positions
    for(i_c=0;i_c<colvar.nconst;i_c++){
            for(i=0;i<colvar.natoms[i_c];i++){
                    j=colvar.cvatoms[i_c][i];
                    getPosition(j,coord_from);
                    mtd_data.pos[j][0]=coord_from.x;
                    mtd_data.pos[j][1]=coord_from.y;
                    mtd_data.pos[j][2]=coord_from.z;
            }
    }

// CV evaluation
    restraint(&mtd_data);
// step increment
    mtd_data.istep++;
}

void PREFIX mtd_data_init()
{
 int i ;
 
 // simparameters contains: cell,dt,lattice and all attributes coming from the namd input file
 SimParameters *spar=Node::Object()->simParameters; 
 // molecule contains atommass, atomcharge,isHydrogen,isOxygen,isWater,bonds_for_atoms,angle_for_atoms methods
 // and other public attributes numAtoms,numBonds,numDihedrals 
 mtd_data.natoms=Node::Object()->molecule->numAtoms; 
 // allocate the common vectors 
 mtd_data.charge=(real *)calloc(mtd_data.natoms,sizeof(real)); 
 mtd_data.mass=(real *)calloc(mtd_data.natoms,sizeof(real)); 
 mtd_data.pos=float_2d_array_alloc(mtd_data.natoms,3);
 mtd_data.vel=float_2d_array_alloc(mtd_data.natoms,3);
 mtd_data.force=float_2d_array_alloc(mtd_data.natoms,3);
 sprintf(hills.dir, ".");

 for(i=0;i<mtd_data.natoms;i++){mtd_data.mass[i]=Node::Object()->molecule->atommass(i);/* printf("ATOMMASS %d %f\n",i,mtd_data.mass[i]);*/};
 for(i=0;i<mtd_data.natoms;i++){mtd_data.charge[i]=Node::Object()->molecule->atomcharge(i);/*printf("ATOMCHARGE %d %f\n",i,mtd_data.charge[i]);*/};
 strcpy(mtd_data.metaFilename,spar->metaFilename);
 mtd_data.fplog = stdout;
 mtd_data.repl=-1;// always consider it as normal md
 mtd_data.dt=spar->dt;
 mtd_data.eunit = 1.;
 mtd_data.ionode = 1;
 mtd_data.istep_old = -1;
 mtd_data.istep = 0;
 logical.not_same_step=1;
}
void PREFIX rvec2vec(rvec rv,Vector *v)
{
v->x=rv[0];
v->y=rv[1];
v->z=rv[2];
}
#endif
// Vector operation, angle etc etc...
#if !defined(PLUMED_GROMACS)
void PREFIX oprod(const rvec a,const rvec b,rvec c)
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}
real PREFIX iprod(const rvec a,const rvec b)
{
  return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}
real PREFIX norm(const rvec a)
{
  return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}
real PREFIX norm2(const rvec a)
{
  return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
}
real PREFIX cos_angle(const rvec a,const rvec b)
{
  real   cos;
  int    m;
  real aa,bb,ip,ipa,ipb;
  ip=ipa=ipb=0.0;
  for(m=0; (m<3); m++) {
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  cos=ip/sqrt(ipa*ipb);
  if (cos > 1.0)
    return  1.0;
  if (cos <-1.0)
    return -1.0;
  return cos;
}
real PREFIX dih_angle(rvec xi, rvec xj, rvec xk, rvec xl,
               rvec r_ij,rvec r_kj,rvec r_kl,rvec m,rvec n,
               real *cos_phi,real *sign)
{
  real ipr,phi;
  real mod_rij, mod_rkj, mod_rkl;

  minimal_image(xi, xj, &mod_rij, r_ij);
  minimal_image(xk, xj, &mod_rkj, r_kj);
  minimal_image(xk, xl, &mod_rkl, r_kl); 

  oprod(r_ij,r_kj,m);               
  oprod(r_kj,r_kl,n);               
  *cos_phi=cos_angle(m,n);          
  phi=acos(*cos_phi);               
  ipr=iprod(r_ij,n);                
  (*sign)=(ipr<0.0)?-1.0:1.0;
  phi=(*sign)*phi;                 
                                   
  return phi;
}
void PREFIX clear_rvec(rvec a)
{
  a[0]=0.0;
  a[1]=0.0;
  a[2]=0.0;
}
#endif

// MPI stuff
void PREFIX plumed_sum(struct mtd_data_s *mtd_data,int nr,real r[]){
  static real *buf=NULL;
  static int nalloc=0;
  int i;
#ifdef MPI
  if (nr > nalloc) {
    nalloc = nr;
    srenew(buf,nalloc);
  }

  if(sizeof(real)==sizeof(double)){
    MPI_Allreduce(r,buf,nr,MPI_DOUBLE,MPI_SUM,mtd_data->comm);
  } else if(sizeof(real)==sizeof(float)){
    MPI_Allreduce(r,buf,nr,MPI_FLOAT,MPI_SUM,mtd_data->comm);
  } else assert(0);

  for(i=0; i<nr; i++)
    r[i] = buf[i];
#endif
};

void PREFIX plumed_sumi(struct mtd_data_s *mtd_data,int nr,int r[]){
  static int *buf=NULL;
  static int nalloc=0;
  int i;
#ifdef MPI
  if (nr > nalloc) {
    nalloc = nr;
    srenew(buf,nalloc);
  }

    MPI_Allreduce(r,buf,nr,MPI_INT,MPI_SUM,mtd_data->comm);

  for(i=0; i<nr; i++)
    r[i] = buf[i];
#endif
};

void PREFIX plumed_intersum(struct mtd_data_s *mtd_data,int nr,real r[]){
  static real *buf=NULL;
  static int nalloc=0;
  int i;
#ifdef MPI
  if (nr > nalloc) {
    nalloc = nr;
    srenew(buf,nalloc);
  }

  if(sizeof(real)==sizeof(double)){
    MPI_Allreduce(r,buf,nr,MPI_DOUBLE,MPI_SUM,mtd_data->intercomm);
  } else if(sizeof(real)==sizeof(float)){
    MPI_Allreduce(r,buf,nr,MPI_FLOAT,MPI_SUM,mtd_data->intercomm);
  } else assert(0);

  for(i=0; i<nr; i++)
    r[i] = buf[i];
#endif
};

int PREFIX plumed_comm_size(struct mtd_data_s *mtd_data){
#ifdef MPI
  int size;
  MPI_Comm_size(mtd_data->comm,&size);
  return size;
#else
  return 1;
#endif
};
int PREFIX plumed_comm_rank(struct mtd_data_s *mtd_data){
#ifdef MPI
  int rank;
  MPI_Comm_rank(mtd_data->comm,&rank);
  return rank;
#else
  return 0;
#endif
};

// random number
#if !defined(PLUMED_GROMACS)
real PREFIX rando(int *ig)
     /* generate a random number. */
{
  int  irand;

  int  m    = 100000000;
  real rm   = 100000000.0;  /* same number as m, but real format */
  int  m1   = 10000;
  int  mult = 31415821;

  real r;
  int  irandh,irandl,multh,multl;

  irand = abs(*ig) % m;

  /* multiply irand by mult, but take into account that overflow
   * must be discarded, and do not generate an error.
   */
  irandh = irand / m1;
  irandl = irand % m1;
  multh  = mult / m1;
  multl  = mult % m1;
  irand  = ((irandh*multl+irandl*multh) % m1) * m1 + irandl*multl;
  irand  = (irand + 1) % m;

  /* convert irand to a real random number between 0 and 1. */
  r = (irand / 10);
  r = r * 10 / rm;
  if ((r <= 0) || (r > 1))
    r = 0.0;
  *ig = irand;

  return r;
}
#endif
// different init_metadyn
#if defined (PLUMED_GROMACS)
void init_metadyn(int natoms, int ePBC, real *charge, real *mass, 
                  real dt, int repl_ex_nst, plumed_repl_ex repl_ex,
                  const t_commrec *mcr, FILE *fplog)
#elif defined (NAMD)
void PREFIX init_metadyn()
#elif defined (ACEMD)
void init_metadyn(int natoms, real *charge, real *mass, 
                  real dt, int repl, int nrepl, 
                  real rte0, real rteio, char *metainp, real box[3])
#elif defined (OPEP)
void init_metadyn_(int *atoms, real *ddt, int *pbc_opep, 
                   int *repl, int *nrepl,real *rte0, real *rteio, real *mass,
                   char *lpath, char *logfile, char *metainp, int ll, int mm, int jj) 
#elif defined (DL_POLY)
void init_metadyn_(int *atoms, real *ddt, real *mass, real *charge, int *imcon, real *eunit, char *metainp, int pp) 
#elif defined (DRIVER)
void init_metadyn_(int *atoms, real *mass, real *charge, int *pbc, real *box, char *metainp, int ll)
#elif defined (STANDALONE)
void init_metadyn_(int *atoms, real *mass, real *charge, int *pbc, real *box, real *tstep, int *nstep, real *myboltz, real* ampli, char *metainp, int ll)
#elif defined (AMBER)
void init_metadyn__( int *atoms, real *ddt,real *mass, real *charge ,char *metainp, int ll)
#elif defined (FHIAIMS)
void init_metadyn__(int *atoms, real *dt, real *mass, char *metainp, int ll)
#endif
{
  int i;
  long int j;
  real k;
  char stringa[800];
  FILE *fp;

#if defined (PLUMED_GROMACS)
  int repl,nrepl;
  real rte0,rteio;
  if(mcr && MASTER(mcr)){
      repl=(repl_ex_nst>0?repl_ex->repl:-1);
      nrepl=(repl_ex_nst>0?repl_ex->nrepl:1);
#ifdef GROMACS4
      rte0=(repl_ex_nst>0?repl_ex->q[0]:0);
      rteio=(repl_ex_nst>0?repl_ex->q[repl_ex->repl]:0);
#else
      rte0=(repl_ex_nst>0?repl_ex->temp[0]:0);
      rteio=(repl_ex_nst>0?repl_ex->temp[repl_ex->repl]:0);
#endif
  } else {
      repl=0;
      nrepl=0;
      rte0=0.;
      rteio=0.;
  };
   mtd_data.ionode = 1;
#ifdef MPI
// these are done here since they are needed for the plumed_sum calls
// we have to find a more consistent way to do it
   if(mcr && !MASTER(mcr)) mtd_data.ionode = 0;
   if(mcr && PAR(mcr)) {
#ifdef GROMACS4
     if(mcr->dd) mtd_data.comm=mcr->dd->mpi_comm_all;
     else mtd_data.comm=mcr->mpi_comm_mysim;
#else
     mtd_data.comm=MPI_COMM_WORLD;
#endif
   } else {
     mtd_data.comm=MPI_COMM_SELF;
   }
#ifdef GROMACS4
   if(mtd_data.ionode && mcr && mcr->ms) {
     mtd_data.intercomm=mcr->ms->mpi_comm_masters;
   } else {
     mtd_data.intercomm=MPI_COMM_NULL;
   }
#else
  if(mtd_data.ionode && nrepl>1) {
    mtd_data.intercomm=MPI_COMM_WORLD;
  } else {
    mtd_data.intercomm=MPI_COMM_NULL;
  }
#endif

  if(mcr && PAR(mcr)){
    plumed_sumi(&mtd_data,1,&repl);
    plumed_sumi(&mtd_data,1,&nrepl);
    plumed_sum(&mtd_data,1,&rte0);
    plumed_sum(&mtd_data,1,&rteio);
  }
#endif
  mtd_data_init (ePBC, charge, mass, natoms, dt, repl_ex_nst, repl, nrepl, 
                 rte0, rteio, mcr, fplog);
#elif NAMD
  mtd_data_init();
#elif ACEMD
  mtd_data_init(  charge, mass, natoms, dt, repl, nrepl, 
		  rte0, rteio, metainp,box);
#elif OPEP
  mtd_data_init(*pbc_opep,*ddt,*atoms,*repl,*nrepl,*rte0,*rteio,mass,lpath,logfile,metainp);
#elif DL_POLY
  mtd_data_init( *atoms , *ddt , mass, charge , imcon, eunit, metainp);
#elif DRIVER
  mtd_data_init( *atoms, mass, charge, metainp, *pbc, box);
#elif STANDALONE 
  mtd_data_init( *atoms, mass, charge, *pbc, box, tstep, nstep, myboltz, ampli, metainp );
#elif AMBER
   mtd_data_init(*atoms, *ddt, mass,charge, metainp);
#elif FHIAIMS
   mtd_data_init(*atoms, *dt , mass, metainp);
#endif

  read_restraint(&mtd_data);             // read META_INP

  if(colvar.nconst==0) return;                          // no CVs, no party!
  firstTime = 1;                                        // it is the first step!
  hills.ntothills = 0;
  sprintf(colfilen, "COLVAR");
  sprintf(hilfilen, "%s/HILLS", hills.dir);
  sprintf(locfilen, "%s/LOCK", hills.dir);
  if(logical.remd){                                     // in replica exchange case
    sprintf(colfilen, "%s/COLVAR%i", hills.dir, mtd_data.repl);
    sprintf(hilfilen, "%s/HILLS%i",  hills.dir, mtd_data.repl);
    sprintf(locfilen, "%s/LOCK%i",   hills.dir, mtd_data.repl);
  }
#ifdef STANDALONE 
// restore printout
    if(mtd_data.istep!=0){
            fclose(mtd_data.fplog);
            mtd_data.fplog = fopen("PLUMED.OUT","a");
            firstTime = 0;                                        // it is the first step!
            if(logical.do_hills){
                    FILE *file;
                    // if the file exists then restart
                        file = fopen(hilfilen, "r");
                        if(file!=NULL){
                                 logical.restart_hills = 1;
                        }  else {
                                 logical.restart_hills = 0;
                        }
            }
    } else {
            firstTime=1;
    }
    fflush(mtd_data.fplog);
#endif
  if(logical.do_hills){
    hills.ntothills = STACKDIM ;// dimesion of hills vector
    hills.ss0_t     = float_2d_array_alloc(hills.ntothills,colvar.nconst); 
    hills.ww        = float_1d_array_alloc(hills.ntothills);
    colvar.delta_s  = float_2d_array_alloc(hills.ntothills,colvar.nconst);
    if(logical.restart_hills) {
      read_hills(&mtd_data,1,hills.first_read);    	// if restart
    } else {                                           	// if not restart, backup old file
#if STANDALONE 
      if(nwalkers==1 &&  firstTime==1 ){
#else
      if(nwalkers==1){
#endif
        sprintf(stringa, "%s.old", hilfilen);
        rename(hilfilen, stringa);
        sprintf(stringa, "%s.old", colfilen);
        rename(colfilen, stringa);
      }

    }
    for(j=0;j<GTAB;j++) {
      k = (real) DP2CUTOFF/GTAB*j;
      hills.exp[j] = exp(-k);
    }  
  }else{
  // evtl move the old colvar file
      if( firstTime==1 ){
        sprintf(stringa, "%s.old", hilfilen);
        rename(hilfilen, stringa);
        sprintf(stringa, "%s.old", colfilen);
        rename(colfilen, stringa);
      }
  }
}
void PREFIX minimal_image(rvec pos1, rvec pos2, real *mod_rij, rvec rij)
{
#if defined (NAMD)
Vector Vect1,Vect2,rij_v;
rvec2vec(pos1,&Vect1);
rvec2vec(pos2,&Vect2);
rij_v=Node::Object()->simParameters->lattice.delta(Vect1,Vect2);
*mod_rij=rij_v.length();
rij[0]=rij_v.x;
rij[1]=rij_v.y;
rij[2]=rij_v.z;
#elif defined (PLUMED_GROMACS)
 pbc_dx(&metapbc, pos1, pos2, rij);
 *mod_rij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
#elif defined (OPEP)
real rin;
int i;
for(i=0;i<3;i++){
 rin = pos1[i] - pos2[i];
 if(mtd_data.imcon==1) {
  rij[i] = pbc_mic_(&rin);
 } else {
  rij[i] = rin;
 }
}
*mod_rij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
#elif defined (DL_POLY)
  int zero,one;
  zero=0;
  one=1;
  rij[0]=pos1[0]-pos2[0];   
  rij[1]=pos1[1]-pos2[1];   
  rij[2]=pos1[2]-pos2[2];   
  images_(&(mtd_data.imcon),&zero,&one,&one,mtd_data.cell,&rij[0],&rij[1],&rij[2]); 
  *mod_rij=sqrt(pow(rij[0],2)+pow(rij[1],2)+pow(rij[2],2));
#elif defined (ACEMD) || defined (DRIVER) || defined (AMBER) || defined (STANDALONE) || defined (FHIAIMS)
  int i;
  for(i=0;i<3;i++){
   rij[i] = pos1[i] - pos2[i];
   if(mtd_data.imcon==1){
    rij[i] -= mtd_data.cell[i]*rint(rij[i]/mtd_data.cell[i]);
   }
  }
*mod_rij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
#endif
}
void PREFIX EXIT()
{
#ifdef NAMD
CkExit();
#else
exit(1);
#endif
} 

real  **** PREFIX float_4d_array_alloc(int ii,int jj,int kk,int ll){
  real ****xx;
  real *ptr;
  int i,j,k;
  ptr=(real *)calloc(ii*jj*kk*ll,sizeof(real));
  xx=(real ****)calloc(ii,sizeof(real ***));
  for (i=0;i<ii;i++){
    xx[i]=(real ***)calloc(jj,sizeof(real **));
    for(j=0;j<jj;j++){
      xx[i][j]=(real **)calloc(kk,sizeof(real *));
      for(k=0;k<kk;k++)xx[i][j][k]= & ptr[ll*(kk*(i*jj+j)+k)];
    }
  }
  return xx;
};

real  *** PREFIX float_3d_array_alloc(int ii,int jj,int kk){
  real ***xx;
  real *ptr;
  int i,j;
  ptr=(real *)calloc(ii*jj*kk,sizeof(real));
  xx=(real ***)calloc(ii,sizeof(real **));
  for (i=0;i<ii;i++){
    xx[i]=(real **)calloc(jj,sizeof(real *));
    for(j=0;j<jj;j++) xx[i][j]= & ptr[kk*(i*jj+j)];
  };
  return xx;
};

real  ** PREFIX float_2d_array_alloc(int ii,int jj){
  real **xx;
  real *ptr;
  int i;
  ptr=(real *)calloc(ii*jj,sizeof(real));
  xx=(real **)calloc(ii,sizeof(real *));
  for (i=0;i<ii;i++)xx[i]=& ptr[jj*i];
  return xx;
};

real  * PREFIX float_1d_array_alloc(int ii){
  real *xx;
  xx=(real *)calloc(ii,sizeof(real));
  return xx;
};


int  ** PREFIX int_2d_array_alloc(int ii,int jj){
  int **xx;
  int *ptr;
  int i;
  ptr=(int *)calloc(ii*jj,sizeof(int));
  xx=(int **)calloc(ii,sizeof(int *));
  for (i=0;i<ii;i++) xx[i]= & ptr[jj*i];
  return xx;
};

int  * PREFIX int_1d_array_alloc(int ii){
  int *xx;
  xx=(int *)calloc(ii,sizeof(int));
  return xx;
};

int PREFIX free_4dr_array_alloc(real ****xx,int ii,int jj,int kk){
  int i,j;
  assert(xx[0][0][0]); free(xx[0][0][0]);
  for(i=0;i<ii;i++){
    for(j=0;j<jj;j++){
      assert(xx[i][j]);
      free(xx[i][j]);
    }
    assert(xx[i]); free(xx[i]);
  }
  free(xx);
  return 0;
};

int PREFIX free_3dr_array_alloc(real ***xx,int ii,int jj){
  int i;
  assert(xx[0][0]); free(xx[0][0]);
  for (i=0;i<ii;i++) {
    assert(xx[i]); free(xx[i]);
  }
  assert(xx); free(xx);
  return 0;
};

int PREFIX free_2dr_array_alloc(real **xx,int ii){
  assert(xx[0]); free(xx[0]);
  assert(xx); free(xx);
  return 0;
};

int PREFIX free_1dr_array_alloc(real *xx){
  free(xx);
  return 0;
};

int PREFIX free_2di_array_alloc(int **xx,int ii){
  assert(xx[0]); free(xx[0]);
  assert(xx); free(xx);
  return 0;
};

int PREFIX free_1di_array_alloc(int *xx){
  assert(xx); free(xx);
  return 0;
};

void PREFIX realquicksort ( real *v , int *ind , int left , int right ) {
        int i,last;
        if(left>=right)return;
        swap(v,ind,left,(left+right)/2);
        last=left;
        for(i=left+1;i<=right;i++){
           if ( v[i]<v[left] ) swap(v,ind,++last,i);
        }
        swap(v,ind,left,last);
        realquicksort(v,ind,left,last-1 );
        realquicksort(v,ind,last+1,right);
}
void PREFIX swap (real *v,int *ind,int i,int j){
        real temp;
        int tempi;
        temp=v[i];
        v[i]=v[j];
        v[j]=temp;
        tempi=ind[i];
        ind[i]=ind[j];
        ind[j]=tempi;
}
