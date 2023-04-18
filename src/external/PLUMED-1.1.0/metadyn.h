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
#ifndef EXTERNALS
#ifndef NAMD
   #define MYEXT extern
#else
   #define MYEXT
#endif   
#else
   #define MYEXT
#endif

// common header files
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
// NAMD header files and definition
#if defined (NAMD)
#define PREFIX GlobalMasterMetaDynamics::
#include "ComputeHomePatches.h"
#include "GlobalMaster.h"
#include "GlobalMasterEasy.h"
#include "NamdTypes.h"
#include "SimParameters.h"
#include "Molecule.h"
#include "Node.h"
#include "Vector.h" 
#include "signal.h"

// GROMACS header files
#elif defined (GROMACS3) || defined (GROMACS4)
#define PLUMED_GROMACS
#define PREFIX
#include "config.h"
#include "typedefs.h"
#include "smalloc.h"
#include "pbc.h"
#include "vec.h"
#include "physics.h"
#include "network.h"
#include "nrnb.h"
#include "bondf.h"
#include "random.h"
#include "repl_ex.h"

#ifdef GMX_MPI
#define MPI
#endif

//ACEMD START
#elif defined (ACEMD)
#define PREFIX
#define TRUE 1 
#include <aceplug.h>
//ACEMD END

//OPEP headers
#elif defined (OPEP)
#define PREFIX
#define TRUE 1 

//DRIVER HEADERS
#elif defined (DRIVER)
#define PREFIX
#define TRUE 1

//STANDALONE HEADERS
#elif defined (STANDALONE)
#define PREFIX
#define TRUE 1

//DL_POLY headers
#elif defined (DL_POLY)
#define PREFIX
#define TRUE 0 

//AMBER headers
#elif defined (AMBER)
#define PREFIX
#define TRUE 0 

//FHIAIMS headers ????
#elif defined (FHIAIMS)
#define PREFIX
#define TRUE 0 
#endif

// Eventually, we include the mpi header
#ifdef MPI
#include "mpi.h"
#endif

#if defined (GROMACS4)
// this structure is needed for GROMACS4 (it has to be consistent with the definition in repl_ex.c)
    struct gmx_repl_ex {
      int  repl;
      int  nrepl;
      real temp;
      int  type;
      real *q;
      bool bNPT;
      real *pres;
      int  *ind;
      int  nst;
      int  seed;
      int  nattempt[2];
      real *prob_sum;
      int  *nexchange;
    };
    typedef gmx_repl_ex_t plumed_repl_ex;
#endif
#if defined (GROMACS3)
// this is needed because of the inconsistency in the definition of gmx_repl_ex_t in versions 3 and 4
    typedef gmx_repl_ex_t* plumed_repl_ex;
#endif

// here we provide alternatives of gromacs macros for the other codes
#if ! defined (PLUMED_GROMACS)
typedef double real;
typedef double rvec[3];
#define snew(ptr,nelem) (ptr)= (nelem==0 ? NULL : (typeof(ptr)) calloc(nelem,sizeof(*(ptr)))) 
#define sfree(ptr) if (ptr != NULL) free(ptr)   
#define srenew(ptr,nelem) (ptr)=(typeof(ptr)) realloc(ptr,(nelem)*sizeof(*(ptr)))
#endif


// common global structures and definitions
#define nconst_max 20					// fixed maximum number of COLVARS
#define DP2CUTOFF 6.25					// sigma^2 considered for gaussian
#define GTAB 1000000					// mesh for exponential tablature
#ifndef M_PI
#define M_PI 3.14159265
#endif
#ifndef M_2PI
#define M_2PI 6.2831853 
#endif
// in case of namd kcal/mol/K 
#if defined (NAMD) || defined (AMBER) || defined (OPEP) || defined (DRIVER) || defined (ACEMD) 
#define BOLTZ 0.001987191 
// in case of dl_poly 10Joule/mol/K 
#elif defined (DL_POLY)
#define BOLTZ 0.831451115
#elif defined (FHIAIMS)
#define BOLTZ 0.00000316681535
// arbitrary choice for standalone (up to user need)
#elif defined (STANDALONE)
#define BOLTZ 1.0 
#endif
// path dimensions
#define MAXATOMS_PATH 230 
#define NMAX_PATH 8 
#define MAXFRAMES_PATH 22 
#define MAXATOMS_RMSD 230
#define MAXCHARS_PATH 40
// cmap 
#define MAXDIM_CMAP 3800 
#define MAXNUM_GROUP 10
#define MAXATOM_GROUP 30
// stackdimension for hills
#define STACKDIM  10000 

// Structure containing the parsed input file
typedef struct {
  int     nlines;     // number of lines
  int*    nwords;     // number of words (in each line)
  char*** words;      // words in each line (+ an extra empty word as a terminator(
  int     ngroups;    // number of group definitions found
  char**  groupnames; // names of groups
  int*    natoms;     // number of members (for each group)
  int**   atoms;      // members (for each group)
} t_plumed_input;


// structure for gradient projections
// linked list
struct  coupling_ll {
        int *at1,*at2,nat1,nat2;
        struct coupling_ll  *next_elem;
};
// couple of cv containing the linked list
struct  el_couple {
    int cv1,cv2;
    struct coupling_ll *first_elem;
};
struct el_diagonal {
    struct coupling_ll *first_elem;
    real *accm;
    real *acct;
};
struct proj_grad_s {
  int *list;
  int  nlist; 
  int nvar,ncouples;
  struct el_couple *couple;
  struct  el_diagonal *diagonal;
  real  **matrix;
  real  **invmatrix;
  real  volume;
  real  *averages;
  real  *prec;
};

// common NAMD/GROMACS interface
struct mtd_data_s
{
   int  natoms;
   real **pos;
   real **vel;
   real **force;
   real *charge;
   real *mass;
   int  repl;
   int  nrepl;
   int  repl_ex_nst;   
   real rte0;
   real rteio;
   real time;
   FILE *fplog;
   char metaFilename[120];
   real dt;
   int  istep;
   int  istep_old;
   real eunit;
   int    ionode;                                        // this is true only in the master node

#ifdef MPI
// communicators for parallel PLUMED
   MPI_Comm comm;       // INTRA replica
   MPI_Comm intercomm;  // INTER replica (for bias-exchange and ptmetad)
#endif

// code-specific definitions
#if defined (PLUMED_GROMACS)
   const t_commrec *mcr;  
   real cell[3][3];
#ifdef GROMACS4   
   int  ePBC;
#endif

#elif defined (OPEP)
   char log[120];
   int imcon;

#elif defined (DL_POLY)
   int imcon;
   real cell[9];

#elif defined (ACEMD) || defined (AMBER) || defined (DRIVER) || defined (STANDALONE) || defined (FHIAIMS)
   int  imcon;
   real cell[3];

#endif
#ifdef STANDALONE 
   real ampli;
   real myboltz;
#endif
};

struct logical_s
{
  int    do_hills;					// hills on/off
  int    widthadapt;
  int    restart_hills;					// restart yes/no
  int    upper[nconst_max];				// upper walls
  int    lower[nconst_max];				// lower walls
  int    steer[nconst_max];				// steering cv on/off
  int    remd;						// replica exchange parallel tempering
  int    rpxm;						// replica exchange metadynamics
  int    commitment;					// committors analysis
  int    print;						// do COLVAR 
  int    welltemp;
  int    lreflect[nconst_max];
  int    ureflect[nconst_max];
  int    debug;
  int    not_same_step;
  int    meta_inp;
  int    parallel_hills;                                // hills sum is parallelized
  int    debug_derivatives;
  int    enable_untested_features;
  int    do_grid;
  int    donot_spline;
  int    debug_grid;
};

struct adapt
{
  int    block;
  int    stride;
  int    on;
  real   inf;
  real   sup;
  real   widthmultiplier;                               // Hills width adapt
};

struct colvar_s
{
  int    nconst;					// # CVs
  int    nt_print;					// step for printing colvar and enercv
  int    it;						// PluMeD step counter
  real   delta_r [nconst_max];				// Hills width along CVs
  real   **delta_s;					// Hills width along CVs in time
  struct adapt adapt[nconst_max];                             // adaptive width structure
  real   ff_hills[nconst_max];				// force due to hills
  int    on      [nconst_max];				// hills on/off on a CV
  int    type_s  [nconst_max];				// colvar type (DIST, ...)
  int    nn      [nconst_max];				// used by for numerator exponent
  int    mm      [nconst_max];				// used by for denominator exponent
  real   r_0     [nconst_max];				// used by for binding distance
  real   d_0     [nconst_max];				// used by for binding distance
  int    intpar  [nconst_max][10];                      // array of integers (general use)
  rvec   realpar [nconst_max][10];                      // array of reals (general use)
  rvec   vecpar  [nconst_max][2];                       // array of vecors (general use)
  real   *map0   [nconst_max]; 				// inter e intra contact starting maps
  int    groups  [nconst_max];				// energy groups id, other id
  int    type    [nconst_max];				// type id (beta sheet, alpha elicas, none, ecc..)
  real   Amin    [nconst_max];				// A state for committors analysis
  real   Amax    [nconst_max];				// ""
  real   Bmin    [nconst_max];				// B ""
  real   Bmax    [nconst_max];				// ""
  int    commit;					//
  rvec   *myder  [nconst_max];				// derivatives
  real   ss0     [nconst_max];                          // CVs value 
  real   Mss0    [nconst_max];				// Hills width adapt
  real   M2ss0   [nconst_max];				// Hills width adapt
  real   wtemp;	         				// well tempered temperature
  real   simtemp;                                       // simulation temperature
  real   wfactor;                                       // welltemp factor = wtemp/simtemp
  int    list[nconst_max][4];                           // structure definition for list 
  int    natoms  [nconst_max];
  int    *cvatoms[nconst_max];
  int    logic[nconst_max];
  int    cell_pbc[nconst_max];                          // switch for applying pbc (where it applies)
  int    ptmetad_neighbours;
  real    ptmetad_sigma;
  int    align_atoms;
  int    *align_list;
  // projections 
   struct proj_grad_s  pg;
};

struct hills_s
{
  real   wwr;						// Hill height
  real   rate;						// Hill deposition rate
  real   max_height;    				// Maximum height of added hills (0 means NO maximum)
  int    max_stride;    				// Maximum stride between hills (0 means NO maximum)
  real   *ww;						// Hill height history
  long int n_hills;					// Hills added
  int    nt_hills;					// Period in step to add hills
  int    nr_hills;					// Period in step to read hills
  real   **ss0_t;					// Hills center history
  char   dir[100];					// HILLS place
  long int ntothills;					// max number of hills
  real   exp[GTAB];					// table for exponential
  long int read;
  fpos_t line_counter;
  int    first_read;
};

struct wall
{
  real   upper   [nconst_max];				// upper limit where start the wall
  real   lower   [nconst_max];				// lower limit
  real   lsigma  [nconst_max];				// lower force constant
  real   sigma   [nconst_max];				// upper force constant
  real   fwall   [nconst_max];				// force due to wall
  int    uexp    [nconst_max];				// upper softwall exponent
  int    lexp    [nconst_max];				// lower softwall exponent
  real   ueps    [nconst_max];				// redux factor for upper wall 
  real   leps    [nconst_max];				// redux factor for lower wall 
  real   uoff    [nconst_max];                          // offset for upper wall
  real   loff    [nconst_max];                          // offset for lower wall     
  int    st_inv  [nconst_max];
};

struct steer {
 real    pos     [nconst_max];				// position of umbrella
 real    delta   [nconst_max];				// increment of umbrella along cv
 real    spring  [nconst_max];				// elastic constant of umbrella
 real    max     [nconst_max];				// limit
 real    start   [nconst_max];				// start position
 int     sign    [nconst_max];
 int     impose_start [nconst_max];                      // logical for imposing a starting point
} ;

struct grid_s {
 real       min      [nconst_max];                      // GRID inferior limit
 real       max      [nconst_max];                      // GRID superior limit
 real       lbox     [nconst_max];                      // GRID bin size
 real       oldelta  [nconst_max];                      // store old HILLS delta 
 real       dx       [nconst_max];                      // GRID spacing 
 real       minilbox             ;                      // redux GRID bin size 
 real       *pot                 ;                      // array for meta bias
 real       **force              ;                      // array for forces
 real       cutoff               ;                      // store genereal DP2CUTOFF
 real       mem                  ;                      // memory info
 int        bin      [nconst_max];                      // number of bins in total GRID
 int        minibin  [nconst_max];                      // number of bins in redux GRID
 int        period   [nconst_max];                      // periodic ?
 int        index    [nconst_max];                      // to map back the id of active CV
 int        size                 ;                      // size of total GRID 
 int        minisize             ;                      // size of redux GRID
 int        **one2multi          ;                      // from 1d index to multidimensional
 int        ncv                  ;                      // number of ACTIVE CVs 
 long int   nhills               ;                      // Total number of HILLS put on GRID
} ;

struct rmsd_inpack{
       int natoms;
       real r0[3][MAXATOMS_RMSD];
       real r1[3][MAXATOMS_RMSD];
       real mass[MAXATOMS_RMSD];
};

struct rmsd_outpack{
       real r0p[3][MAXATOMS_RMSD];//centered reference  frame
       real r1p[3][MAXATOMS_RMSD];//centered running frame  
       real cmr0[3]; //center of mass of reference frame
       real cmr1[3]; //center of mass of running frame
       real err;
       real derr_dr0[3][MAXATOMS_RMSD];
       real derr_dr1[3][MAXATOMS_RMSD];
       real dderr_dr1_dr1[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
       real dderr_dr0_dr0[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
       real dderr_dr1_dr0[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
       real dderr_dr0_dr1[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
       real d[3][3];
       real dd_dr0[3][3][3][MAXATOMS_RMSD];
       real dd_dr1[3][3][3][MAXATOMS_RMSD];
};
struct rmsd_mini_outpack{
     real r0p[3][MAXATOMS_RMSD];//centered reference  frame
     real r1p[3][MAXATOMS_RMSD];//centered running frame  
     real cmr0[3]; //center of mass of reference frame
     real cmr1[3]; //center of mass of running frame
     real err;
     real derr_dr0[3][MAXATOMS_RMSD];
     real derr_dr1[3][MAXATOMS_RMSD];
     real d[3][3];
     real dd_dr0[3][3][3][MAXATOMS_RMSD];
     real dd_dr1[3][3][3][MAXATOMS_RMSD];
};
// For path in contact map space
struct group_struct{
       int number;
       int numatom[MAXNUM_GROUP];
       int index[MAXNUM_GROUP][MAXATOM_GROUP];
       int index_to_list[MAXNUM_GROUP][MAXATOM_GROUP];
       real rcm[MAXNUM_GROUP][3];
};
struct cmap_pack {
       int index1[MAXDIM_CMAP];
       int index2[MAXDIM_CMAP];       
       int index_from1[MAXDIM_CMAP];
       int index_from2[MAXDIM_CMAP];
       int atoms;
       int list[MAXATOMS_RMSD];
       int nn[MAXDIM_CMAP];
       int nd[MAXDIM_CMAP];
       int number;
       int gnumber;
       real r0[MAXDIM_CMAP];
       real weight[MAXDIM_CMAP];
       real cutoff[MAXDIM_CMAP];
       real cmap[MAXFRAMES_PATH][MAXDIM_CMAP];
       int logical_group;
       struct group_struct group;
       };
struct cmap_inpack{
       real r0[MAXATOMS_PATH][3];
       real cmap[MAXDIM_CMAP];
};
struct cmap_outpack{
       real err;
       real derr_dr0[3][MAXATOMS_PATH];
       real derr_dcm[MAXDIM_CMAP];
};
struct coordinates_frameset {
        int natoms;
        int resid[MAXATOMS_PATH],atmnum[MAXATOMS_PATH];
        char resname[MAXATOMS_PATH][MAXCHARS_PATH];
        char label[MAXATOMS_PATH][MAXCHARS_PATH];
        real pos[MAXATOMS_PATH][3];
        char residue[MAXATOMS_PATH];
        real align[MAXATOMS_PATH];    // for each atom is 1 if used for alignment, else 0
        real displace[MAXATOMS_PATH]; // for each atom is 1 if used for displacement, else 0
        int frameset_to_coord[MAXATOMS_PATH]; // for each atoms in the frameset provides the corresponding index on running coord
        int frameset_to_align[MAXATOMS_PATH]; // for each atoms in the frameset provides the corresponding index used alignment(-1 if not involved in rmsd) 
        int nalign;
        int align_to_coord[MAXATOMS_PATH];
        int align_to_frameset[MAXATOMS_PATH];
        int ndisplace;
        //real ndisplace;
        //int displaced[MAXATOMS_PATH];
};
// this is the final unti that contains all the single datas
struct hybrid_elem{
        int cvtype,cvindex;
	//
       	// distance: distance is only one float, but requires the full list of atoms so to copy the derivative 
	//
        real ref_dist;    // reference dist for this frame 
        int  nder;        // ref number of atoms in the derivative 
        int  *backtable;        // ref intdex of the atoms in the derivative (back transfer of forces) 
        // ref derivative       (copy buffer for cvs)
        rvec   *myder  ; 
};
struct hybrid_frameset{
       	// this should contain all the required data (possibly dynamically allocated)
       	// for all the structures needed in the hybrid frameset definition
        // one frame may contain many structures
        int    hbd_nelem,hbd_totalatoms;
        struct hybrid_elem **hbd_elem; 
        int    *backtable;
        rvec   *myder;
}; 

struct sz_data {
        int    number;
        real lambda;
        struct coordinates_frameset **frameset;
        struct hybrid_frameset **hbd_frameset;
        struct hybrid_frameset    hbd_running; //assumin that the number of atoms in each reference  stays constant
        char   names[MAXCHARS_PATH];
        int grad_on,umb_on,mass_on,targeted_on,sqrt_on;
        real gradmax,gradmin,gradk;
        int umblagsteps,umbblocksize,umbstride,umbcount,umbpermanency,countperm;
        real umbtolerance;
        real ****umb_block;// dimensions:  3,natoms,nframes,nblock
        real ***umb_avg;// dimensions:  3,natoms,nframes
        real ***umb_map_block;// dimensions:  tot_con,nframes,nblock
        real **umb_map_avg;// dimensions:  tot_con,nframes
        struct cmap_pack my_cmap_pack;
        char   path_type[10];
        int nneigh,*lneigh,neigh,neigh_time;
        // hybrid path
        int *lhybrid,nhybrid,*lcvhybrid; 
        real **mathybrid;
};

// NAMD CLASS DEFINITION AND PECULIAR METHODS
#ifdef NAMD
class GlobalMasterMetaDynamics : public GlobalMasterEasy 
{

 public:
 GlobalMasterMetaDynamics(); 
 virtual void easy_init(const char *);
 void easy_calc(void);
 void mtd_data_init(  );
 void init_metadyn( );
 void rvec2vec(rvec rv,Vector *v);
#elif defined(PLUMED_GROMACS)
 void mtd_data_init (int ePBC, real *charge, real *mass,
                    int natoms, real dt, int repl_ex_nst, int repl,
                    int nrepl, real rte0, real rteio, const t_commrec *mcr, FILE *fplog);

 void init_metadyn(int natoms, int ePBC, real *charge, real *mass,
                  real dt, int repl_ex_nst, plumed_repl_ex repl_ex,
                  const t_commrec *mcr, FILE *fplog);
 void ptmetad(real *Epota, real *Epotb, real *Epotba, real *Epotab, int a, int b);
 void ptmetad_vbias(int,real*,real*);
 void ptmetad_helper();
 void ptmetad_exchfluct(int);
 void bias_exchange_traj(int nrepl, int *seed, int *ind);
 void meta_force_calculation(int, int, real (*pos)[3], real (*force)[3], real box[3][3]);
 void plumed_setstep(int istep);
 int  plumed_dd_index(int iat);
#ifdef GROMACS4
#define oprod cprod
#endif
#elif OPEP
 void init_metadyn_(int *atoms, real *ddt, int *pbc_opep,
                   int *repl, int *nrepl,real *rte0, real *rteio, real *mass,
                   char *lpath, char *logfile, char *metainp, int ll, int mm, int jj);
 void meta_force_calculation_(real *pos, real *force);
 real pbc_mic_(real *rin);
 void mtd_data_init(int pbc, real tstep,int atoms, int repl, int nrepl, real rte0, real rteio, real *mass, char *lpath, char *logfile, char *metainp);
 void ptmetad_vbias(int,real*,real*);
 void share_bias_(int *rep, real *Vbias, real *Vbiasx); 
 void bias_exchange_(int *nrep, int *biaseed, int *ind);
 void bias_exchange_traj(int nrepl, int *seed, int *ind);
 void switch_fluct_(int *rep); 
 void ptmetad_exchfluct(int);
#elif ACEMD 
 void mtd_data_init( real *charge, real *mass,
                    int natoms, real dt, int repl,
		     int nrepl, real rte0, real rteio, char *metainp, real *box);
 void init_metadyn(int natoms, real *charge, real *mass, 
                   real dt, int repl, int nrepl, 
                   real rte0, real rteio, char *metainp, real box[3]);
 void meta_force_calculation(struct aceplug_sim_t* );
#elif DL_POLY
 void init_metadyn_(int *atoms, real *ddt, real *mass, real *charge, int *imcon, real *eunit, char *metain, int pp);
 void mtd_data_init(int atoms, real dt , real *mass, real *charge , int *imcon, real *eunit, char *metainp);
 void meta_force_calculation__(real *cell,real *xxx,real *yyy,real *zzz,real *fxx,real *fyy,real *fzz); 
 void images_(int *i,int *j,int *k,int *natoms,real *cell,real *xxx,real *yyy,real *zzz); 
#elif DRIVER
 void init_metadyn_(int *atoms, real *mass, real *charge, int *pbc, real *box, char *metainp, int ll);
 void mtd_data_init(int atoms, real *mass, real *charge, char *metainp, int pbc, real *box);
 void cv_calculation_(real *box, real *pos, int *ncv, real *cv);
 void ptmetad(real *Epota, real *Epotb, real *Epotba, real *Epotab, int a, int b);
 void ptmetad_vbias(int,real*,real*);
 void ptmetad_sharepot(int nrepl, int repl);
 void bias_exchange_traj(int nrepl, int *seed, int *ind);
 void ptmetad_exchflut(int repl);
 void ptmetad_exchfluct(int);
#elif STANDALONE 
 void init_metadyn_(int *atoms, real *mass, real *charge, int *pbc, real *box, real *tstep, int *nstep, real *myboltz, real *ampli, char *metainp, int ll);
 void mtd_data_init(int atoms, real *mass, real *charge, int pbc, real *box,  real *tstep, int *nstep , real* myboltz, real *ampli , char *metainp );
 void cv_calculation_standalone_(real *box, real *pos, real *force, real *ene);
#elif AMBER
 void init_metadyn__(int *atoms, real *ddt, real *mass, real *charge ,char *metainp, int ll);
 void mtd_data_init(int atoms,  real ddt , real *mass, real *charge, char *metainp);
 void meta_force_calculation__(real *box,real *rrr,real *fff, int *nn); 
#elif FHIAIMS
 void init_metadyn__();
 void init_metadyn_();
 void mtd_data_init();
 void meta_force_calculation__(); 
 void meta_force_calculation_(); 
#endif

#if ! defined (PLUMED_GROMACS)
// Vector operation (locally reimplemented for codes either than GROMACS)
 real rando(int *ig);
 void oprod(const rvec a,const rvec b,rvec c);
 real iprod(const rvec a,const rvec b);
 real norm(const rvec a);
 real norm2(const rvec a);
 real cos_angle(const rvec a,const rvec b);
 real dih_angle(rvec xi, rvec xj, rvec xk, rvec xl,
               rvec r_ij,rvec r_kj,rvec r_kl,rvec m,rvec n,
               real *cos_phi,real *sign);
 void clear_rvec(rvec a);
#endif

// common declarations
 void EXIT();
 void read_restraint(struct mtd_data_s *mtd_data);
 void read_defaults(  );
 void plumed_read_input(t_plumed_input* input,FILE* file,FILE* log);
 void plumed_clear_input(t_plumed_input*);
 int  plumed_get_group(const char *word,int **atoms,int n,t_plumed_input*,FILE *log);
 void plumed_error(const char*);
 int  plumed_atoi(const char* word);
 double plumed_atof(const char* word);
 int  plumed_get_words(char* line,char*** words);
 void plumed_sum    (struct mtd_data_s *mtd_data,int n,real* buffer);
 void plumed_sumi   (struct mtd_data_s *mtd_data,int n,int* buffer);
 void plumed_intersum(struct mtd_data_s *mtd_data,int n,real* buffer);
 int  plumed_comm_size(struct mtd_data_s *mtd_data);
 int  plumed_comm_rank(struct mtd_data_s *mtd_data);
 int  seek_word(char **word, const char *wanted);
 int  seek_word2(char **word, const char *wanted, int is);
// PBS
 void minimal_image(rvec pos1, rvec pos2, real *mod_rij, rvec rij);
// HILLS stuff
 void hills_add(struct mtd_data_s *mtd_data);
 void hills_push(struct mtd_data_s *mtd_data,real ww, real* ss,real* delta);
 real hills_engine(real*,real*);
 real hills_engine_dp(int ih,real* ss0,real* dp);
 void hills_force();
 void apply_forces(struct mtd_data_s *mtd_data );
 real soft_walls_engine(real*,real*);
 void steer_cv(int);
 void print_colvar_enercv(real time_s);
 void hills_adapt();
 void commit_analysis();
 void read_hills(struct mtd_data_s *mtd_data, int restart, int first_read);
 void hills_reallocate(struct mtd_data_s *mtd_data);
// CV routines
 void restraint(struct mtd_data_s *mtd_data);
 void test_derivatives(struct mtd_data_s *mtd_data);
 void debug_derivatives(struct mtd_data_s *mtd_data,int);
// reading...
 int  read_dist       (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_mindist    (char **word,int count,t_plumed_input *input,           FILE *fplog); 
 int  read_coord      (char **word,int count,t_plumed_input *input,           FILE *fplog); 
 int  read_angle      (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_torsion    (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_hbonds     (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_dipole     (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_rgyr       (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_waterbridge(char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_path       (char **word,int count,t_plumed_input *input,           FILE *fplog);
 int  read_position   (char **word,int count,t_plumed_input *input,           FILE *fplog); 
 int  read_elstpot    (char **word,int count,t_plumed_input *input,           FILE *fplog); 
 int  read_puckering  (char **word,int count,t_plumed_input *input,           FILE *fplog); 
// these variables read lines directly from the input, thus they need the iline pointer to increment it
 int  read_alfabeta   (char **word,int count,t_plumed_input *input,int *iline,FILE *fplog);
 int  read_rmsdtor    (char **word,int count,t_plumed_input *input,int *iline,FILE *fplog);
 int  read_dihcor     (char **word,int count,t_plumed_input *input,int *iline,FILE *fplog);
// calculating
 void dist_restraint       (int i_c, struct mtd_data_s *mtd_data);
 void mindist_restraint    (int i_c, struct mtd_data_s *mtd_data);
 void coord_restraint      (int i_c, struct mtd_data_s *mtd_data);
 void angle_restraint      (int i_c, struct mtd_data_s *mtd_data);
 void torsion_restraint    (int i_c, struct mtd_data_s *mtd_data);
 void alfabeta_restraint   (int i_c, struct mtd_data_s *mtd_data);
 void hbonds_restraint     (int i_c, struct mtd_data_s *mtd_data);
 void dipole_restraint     (int i_c, struct mtd_data_s *mtd_data);
 void radgyr_restraint     (int i_c, struct mtd_data_s *mtd_data);
 void rmsdtor_restraint    (int i_c, struct mtd_data_s *mtd_data);
 void dihcor_restraint     (int i_c, struct mtd_data_s *mtd_data);
 void waterbridge_restraint(int i_c, struct mtd_data_s *mtd_data);
 void spath_restraint      (int i_c, struct mtd_data_s *mtd_data);
 void zpath_restraint      (int i_c, struct mtd_data_s *mtd_data);
 void position_restraint   (int i_c, struct mtd_data_s *mtd_data);
 void elstpot_restraint    (int i_c, struct mtd_data_s *mtd_data);
 void puckering_restraint  (int i_c, struct mtd_data_s *mtd_data);
 real coscut(real r,real r_0, real cut);
 real dcoscut(real r,real r_0, real cut);
 real tprod (rvec a,rvec b, rvec c );
 real generate_R(int i_c,struct mtd_data_s *mtd_data,rvec* R,rvec Rp,rvec Rdp,rvec RpxRdp);
 real puckering_Zeta (rvec R,rvec RpxRdp,real mod);
 void puckering_gradZeta(rvec gradZ,int index_i,int index_j,real* z,rvec* R,rvec Rp, rvec Rdp, rvec RpxRdp,real mod);
 real puckering_Q(real *z);
 void puckering_gradQ(rvec gradQ, real* z, real Q,rvec*R, rvec Rp, rvec Rdp, rvec RpxRdp,real mod,int index);
 real puckering_phi(real *z);
 void puckering_gradphi(rvec gradphi, real* z,rvec*R, rvec Rp, rvec Rdp, rvec RpxRdp,real mod,int index);
 real puckering_theta(real * z);
 void puckering_gradtheta(rvec gradtheta, real* z,rvec*R,rvec Rp,rvec Rdp,rvec RpxRdp,real mod,int index);
// GRID stuff
real spline(int ndim,real *dx,real *where,real *tabf,real *tabder,int* stride,real *der);
void grid_initialize(struct grid_s *grid);
void grid_resize_minigrid(struct grid_s *grid, real* delta, real cutoff);
void grid_create_one2multi(struct grid_s *grid);
void grid_addhills(struct grid_s *grid, real ww, real* ss,real* delta,int rank,int npe);
real grid_getstuff(struct grid_s *grid, real* ss0,real* force);
int  grid_multi2one(struct grid_s *grid, int* index_nd);
// misc routines
 void cmap_running(struct cmap_inpack *inpack, struct cmap_pack *my_cmap_pack);
 void cmdist_eval(int frame,struct cmap_inpack *inpack,struct cmap_outpack *outpack, 
      struct cmap_pack *my_cmap_pack,int dr1_calc);
 real pow2(real x);
 void power(real x,int p,int q,real *xp,real *xq);
 void read_sz_map(struct sz_data *my_sz, char file_maps[129], char file_maps2[129],
                        char file_group[129], FILE *fplog);
 int  read_sz_rmsd(struct sz_data *my_sz, FILE *fplog);
 int  read_sz_hybrid(struct sz_data *my_sz, FILE *fplog);
 int  read_sz_coord (char *filename, struct coordinates_frameset *p, FILE *fplog);
 // hbd dist
 int  hbd_read_dist (FILE *myfile, FILE *fplog, struct hybrid_elem *  );
 int  hbd_copy_dist ( struct hybrid_elem *elem ); 
 real hbd_metrics_dist ( struct hybrid_elem *run,  struct hybrid_elem *ref );
 // hbd angle
 int  hbd_read_angle (FILE *myfile, FILE *fplog, struct hybrid_elem *  );
 int  hbd_copy_angle ( struct hybrid_elem *elem ); 
 real hbd_metrics_angle ( struct hybrid_elem *run,  struct hybrid_elem *ref );
 // hbd torsion
 int  hbd_read_torsion (FILE *myfile, FILE *fplog, struct hybrid_elem *  );
 int  hbd_copy_torsion ( struct hybrid_elem *elem ); 
 real hbd_metrics_torsion ( struct hybrid_elem *run,  struct hybrid_elem *ref );
 // hbd coord 
 int  hbd_read_coord (FILE *myfile, FILE *fplog, struct hybrid_elem *  );
 int  hbd_copy_coord ( struct hybrid_elem *elem ); 
 real hbd_metrics_coord ( struct hybrid_elem *run,  struct hybrid_elem *ref );
 // hbd target 
 int  hbd_read_target (FILE *myfile, FILE *fplog, struct hybrid_elem * , struct mtd_data_s *mtd_data );
// int  hbd_copy_target ( struct hybrid_elem *elem ); 
// real hbd_metrics_target ( struct hybrid_elem *run,  struct hybrid_elem *ref );

 int  hbd_collect_config ( struct hybrid_frameset running  );
 int  hbd_metrics ( struct hybrid_frameset *running , struct hybrid_frameset *reference , struct cmap_outpack *outpack, real **mat); 
 void msd_calculation(struct coordinates_frameset *pframeset,struct cmap_inpack *c_inpack,
                             struct cmap_outpack *c_outpack,real dmsd_dr1[3][MAXATOMS_PATH]);
 int  rmsd_pack(struct rmsd_inpack inpack,struct rmsd_outpack *outpack,int iopt,int iopt2);
 int  rmsd_mini_pack(struct rmsd_inpack inpack,struct rmsd_mini_outpack *outpack,int iopt);
 void mean_rmsd(struct sz_data *pmy_sz, real dCV_dr1[MAXFRAMES_PATH][3][MAXATOMS_PATH],
                        int i_c, FILE *fplog); 
 void mean_map(struct sz_data *pmy_sz, real dCV_dcm[MAXFRAMES_PATH][MAXDIM_CMAP],
                        int i_c, FILE *fplog);
 void dmsd_calculation(int i_c,struct coordinates_frameset *pframeset,struct cmap_inpack *c_inpack,
                             struct cmap_outpack *c_outpack,real dmsd_dr1[3][MAXATOMS_PATH]);

 // stuff for gradient projection
 void couple2list( struct coupling_ll **atlist ,int *at1,int nat1,int *at2,int nat2);
 void freecouple2list( struct coupling_ll **atlist);
 void scancouple( struct coupling_ll *atlist);
 void setup_projections(struct proj_grad_s * );   
 void calc_projections(struct proj_grad_s *);   

 void cite_please (const char* re, FILE *fplog);
 void disclaimer (FILE *fplog);
 // allocators and destructors
 real *float_1d_array_alloc(int ii);
 real **float_2d_array_alloc(int ii,int jj);
 real ***float_3d_array_alloc(int i,int j,int k);
 real ****float_4d_array_alloc(int i,int j,int k,int l);
 int  *int_1d_array_alloc(int ii);
 int  **int_2d_array_alloc(int ii,int jj);
 int free_1dr_array_alloc(real *xx); // 1d real
 int free_2dr_array_alloc(real **xx,int ii); // 2d real
 int free_3dr_array_alloc(real ***xx,int ii,int jj); // 3d read
 int free_4dr_array_alloc(real ****xx,int ii,int jj, int kk); //4d real
 int free_1di_array_alloc(int *xx);  //1d integer
 int free_2di_array_alloc(int **xx,int ii);  //2d integer
 // neighbour list tools for quicksorting
 void realquicksort ( real *v , int *ind , int left , int right ); // quicksort for neighbour list search
 void swap (real *v,int *ind,int i,int j); // used by quicksort
 void ql77_driver(real m[][4],real* lambda);
 void ql77 (int n,real *x,real *d);
/* COMMON DATA STRUCTURES  */
    MYEXT struct mtd_data_s    mtd_data;
    MYEXT struct steer         cvsteer;
    MYEXT struct wall          cvw;
    MYEXT struct logical_s     logical;
    MYEXT struct colvar_s      colvar;
    MYEXT struct hills_s       hills;
    MYEXT int   firstTime;					// first PluMed step
    MYEXT int    nwalkers;					// # walkers in parallel
    MYEXT char   colfilen[800];					// COLVAR and ENERCV files
    MYEXT char   hilfilen[800];					// HILLS file
    MYEXT char   locfilen[800];					// LOCK file
    MYEXT real   Vhills;					// Hills potential
    MYEXT real   Vwall;						// Wall potential
// path stuff
    MYEXT struct sz_data my_sz_list[NMAX_PATH];
    MYEXT int ic_to_sz[nconst_max];
    MYEXT int nsz;
    MYEXT int kill_me[nconst_max];
// GRID stuff
    MYEXT struct grid_s       grid;
#ifdef DL_POLY
    MYEXT real tstep;
#elif defined (PLUMED_GROMACS)
    MYEXT t_pbc  metapbc;
#elif NAMD
};
#endif
