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

void  PREFIX spath_restraint(int i_c, struct mtd_data_s *mtd_data) {

        int iat,i,ii,j;
        real s,ci_vec,tmp1;
        struct coordinates_frameset *pmy_coord1;
        struct hybrid_frameset *hbd_pmy;
        struct sz_data *pmy_sz;
        struct cmap_inpack inpack;
        struct cmap_outpack outpack;
	real ds_temp_dr0[MAXFRAMES_PATH][3][MAXATOMS_PATH];
	real ds_dcm[MAXFRAMES_PATH][MAXDIM_CMAP]; 
	real ds_dr0[3][MAXATOMS_PATH]; 
        int start_avg = 0; 
        real ds_dr1[MAXFRAMES_PATH][3][MAXATOMS_PATH];
        real dmsd_dr1[3][MAXATOMS_PATH];
        int  tot_con;
        real *save_err; 
        int nneigh;

        pmy_sz=&my_sz_list[ic_to_sz[i_c]];

// neigh list ?
        if(pmy_sz->neigh==1 && colvar.it%pmy_sz->neigh_time==0) {
            fprintf(mtd_data->fplog,"|- CALCULATING NEIGHBOUR LIST AT STEP %d\n",colvar.it);
            for(i=0;i< pmy_sz->number;i++)pmy_sz->lneigh[i]=i;
               nneigh=pmy_sz->number;
               save_err=(real *)malloc(pmy_sz->number*sizeof(real));
        }else {
            nneigh=pmy_sz->nneigh;
        }

// delete vectors
        for (i=0;i<3;i++){
                for (j=0;j<colvar.natoms[i_c];j++) {
                        ds_dr0[i][j]=0.;
                }
        }
        ci_vec=0.;
        s=0.;

        if(pmy_sz->umb_on && (colvar.it-pmy_sz->umblagsteps>=0)) start_avg = 1;
                                                               
        if(strcmp(pmy_sz->path_type,"HYBRID") != 0){
               for (i=0;i<colvar.natoms[i_c];i++){
                     iat = colvar.cvatoms[i_c][i];
                     inpack.r0[i][0] = mtd_data->pos[iat][0];
                     inpack.r0[i][1] = mtd_data->pos[iat][1];
                     inpack.r0[i][2] = mtd_data->pos[iat][2];
              }
        }
 
        if(strcmp(pmy_sz->path_type,"CMAP") == 0){ 
         tot_con=pmy_sz->my_cmap_pack.number+pmy_sz->my_cmap_pack.gnumber;
         cmap_running(&inpack,&pmy_sz->my_cmap_pack);
        }
        
        if(strcmp(pmy_sz->path_type,"HYBRID") == 0){ 
            //retrive values of cv and derivatives (one set)
            // assuming they are all the same for all the framesets (generally it is) 
            // in case of distance where the comp weight is on the cv evaluation and not on the 
            // distance acquire the points in CV space 
            // eg: dist between two atoms-> collect the cvval (which is the coordinate) and the  
            // eg: rmsd between two struct-> collect the actual coordinates (which is the coordinate) 
            // note: if D is the distance function from the reference  
            // d D(r_0,r)/dx= dD(r_0,r)/dr dr/dx 
            //   
            // in case of rmsd there is no need for chain rule as the derivative is directly calculated with quaternion
            // d D(r_0,r)/dx via quaternion   (it's like dr/dx=1)
            //
            //
            hbd_collect_config(pmy_sz->hbd_running);
        }

        for(ii=0;ii< nneigh;ii++){
                i=pmy_sz->lneigh[ii];
       
                if(strcmp(pmy_sz->path_type,"CMAP") == 0){
                 cmdist_eval(i,&inpack,&outpack,&pmy_sz->my_cmap_pack,start_avg); 
                } 
                if(strcmp(pmy_sz->path_type,"RMSD") == 0){
                 pmy_coord1=pmy_sz->frameset[i];
                 msd_calculation(pmy_coord1,&inpack,&outpack,dmsd_dr1);
                }
                if(strcmp(pmy_sz->path_type,"DRMS") == 0){
                 pmy_coord1=pmy_sz->frameset[i];
                 dmsd_calculation(i_c,pmy_coord1,&inpack,&outpack,dmsd_dr1);
                }
                if(strcmp(pmy_sz->path_type,"HYBRID") == 0){
                 hbd_pmy=pmy_sz->hbd_frameset[i];
                 // calculate the distance between the two frames
                 hbd_metrics(&pmy_sz->hbd_running,hbd_pmy,&outpack,pmy_sz->mathybrid);
                }
 
                //fprintf(mtd_data->fplog,"ERR %d %f \n",i,outpack.err);
                //fflush(mtd_data->fplog);
             // in case you are calculating the neigh list 
                if(pmy_sz->neigh==1 && colvar.it%pmy_sz->neigh_time==0)save_err[i]=outpack.err;

             // sqrt option
                if(pmy_sz->sqrt_on){
                   if(outpack.err<1.e-6){fprintf(mtd_data->fplog,"|-ERRROR IN ZPATH: too small error %f\n",outpack.err);EXIT();}
                   tmp1=(0.5/sqrt(outpack.err))*exp(-pmy_sz->lambda*sqrt(outpack.err));
                }else{
                   tmp1=exp(-pmy_sz->lambda*outpack.err);
                }
 
                s+=(i+1.0)*tmp1;
                ci_vec+=tmp1;
                for(j=0;j<colvar.natoms[i_c];j++){
			ds_temp_dr0[i][0][j]=(outpack.derr_dr0[0][j])*tmp1;
			ds_temp_dr0[i][1][j]=(outpack.derr_dr0[1][j])*tmp1;
			ds_temp_dr0[i][2][j]=(outpack.derr_dr0[2][j])*tmp1;
                        if((strcmp(pmy_sz->path_type,"RMSD") == 0 || strcmp(pmy_sz->path_type,"DRMS") == 0) && start_avg){
                         ds_dr1[i][0][j]=(dmsd_dr1[0][j])*tmp1*pmy_sz->lambda;
                         ds_dr1[i][1][j]=(dmsd_dr1[1][j])*tmp1*pmy_sz->lambda;
                         ds_dr1[i][2][j]=(dmsd_dr1[2][j])*tmp1*pmy_sz->lambda;
                        }
		}
                if(strcmp(pmy_sz->path_type,"CMAP") == 0 && start_avg){
                 for(j=0;j<tot_con;j++){
                         ds_dcm[i][j]=(outpack.derr_dcm[j])*tmp1*pmy_sz->lambda;
                        } 
                }
	}


        s=s/ci_vec;

	for(j=0;j<colvar.natoms[i_c];j++){
               for(ii=0;ii<nneigh;ii++){
                       i=pmy_sz->lneigh[ii];
        	       ds_dr0[0][j]+=(s-i-1.)*ds_temp_dr0[i][0][j]; 
        	       ds_dr0[1][j]+=(s-i-1.)*ds_temp_dr0[i][1][j]; 
        	       ds_dr0[2][j]+=(s-i-1.)*ds_temp_dr0[i][2][j]; 
                       if((strcmp(pmy_sz->path_type,"RMSD") == 0 || strcmp(pmy_sz->path_type,"DRMS") == 0) && start_avg){ 
                        ds_dr1[i][0][j]=ds_dr1[i][0][j]*(s-i-1.)/ci_vec;
                        ds_dr1[i][1][j]=ds_dr1[i][1][j]*(s-i-1.)/ci_vec;
                        ds_dr1[i][2][j]=ds_dr1[i][2][j]*(s-i-1.)/ci_vec;
                       }
       	    	}
                ds_dr0[0][j]=((pmy_sz->lambda)/ci_vec)*ds_dr0[0][j]; 
       	    	ds_dr0[1][j]=((pmy_sz->lambda)/ci_vec)*ds_dr0[1][j]; 
       	     	ds_dr0[2][j]=((pmy_sz->lambda)/ci_vec)*ds_dr0[2][j]; 
	}

        if(strcmp(pmy_sz->path_type,"CMAP") == 0 && start_avg){
           for(i=0;i<pmy_sz->number;i++){           
               for(j=0;j<tot_con;j++){
          	        ds_dcm[i][j]=ds_dcm[i][j]*(s-i-1.)/ci_vec;
		       }
           }
        }

        for(i=0;i<colvar.natoms[i_c];i++) {
          colvar.myder[i_c][i][0] = ds_dr0[0][i];
          colvar.myder[i_c][i][1] = ds_dr0[1][i];
          colvar.myder[i_c][i][2] = ds_dr0[2][i];
        }

        colvar.ss0[i_c]=s;


// neigh list? do quicksort and deallocate save_err
        if(pmy_sz->neigh==1 && colvar.it%pmy_sz->neigh_time==0){
             //   for(i=0;i<nneigh;i++)printf("BEFORE SORTING %d %f\n",pmy_sz->lneigh[i],save_err[i]);
                realquicksort(save_err,pmy_sz->lneigh,0,nneigh-1);
             //   for(i=0;i<nneigh;i++)printf("AFTER SORTING %d %f\n",pmy_sz->lneigh[i],save_err[i]);
                free(save_err);
        }



        if(pmy_sz->umb_on==1){
          if(strcmp(pmy_sz->path_type,"CMAP") == 0) mean_map(pmy_sz,ds_dcm,i_c,mtd_data->fplog); 
          if(strcmp(pmy_sz->path_type,"RMSD") == 0 || strcmp(pmy_sz->path_type,"DRMS") == 0) mean_rmsd(pmy_sz,ds_dr1,i_c,mtd_data->fplog); 
        }

        return;
}

// ------------------------------------------------------------------------------------------------

int PREFIX read_path(char **word, int count,t_plumed_input *input,FILE *fplog)
{

  int i,iw;
  double lambda, tol;
  double sigma = 0.0;
  char file_maps[129];
  char file_maps2[129];
  char file_group[129];
  char type[2];
  struct sz_data *my_sz;
  int path_help;
  int neigh;

  path_help=0; 
  neigh=0;

  my_sz=&(my_sz_list[nsz]);
  ic_to_sz[count]=nsz;

  my_sz->number = 0;
  my_sz->lambda = 0.;                
  my_sz->my_cmap_pack.logical_group = 0;
  my_sz->umb_on = 0;
  my_sz->grad_on = 0;
  my_sz->targeted_on = 0;
  my_sz->sqrt_on = 0;
 
// targeted md ?
  iw=seek_word(word,"TARGETED");
  if(iw>=0) {
      fprintf(fplog,"|- TARGETED MD: only one frame needed \n"); 
      my_sz->targeted_on=1;
     // my_sz->sqrt_on=1; // targeted is always square rooted
  }
// specific sqrt keyword
  iw=seek_word(word,"SQRT");
  if(iw>=0) {
      fprintf(fplog,"|- SQRT enabled: the measure will be square rooted \n"); 
      my_sz->sqrt_on=1;
  }
// type
  iw=seek_word(word,"TYPE");
  if(iw>=0) {
      sscanf(word[iw+1],"%s",my_sz->path_type);
  } else {
      fprintf(fplog,"|- NEEDED \"TYPE\" KEWORD FOR PATH\n"); 
      path_help=1; 
  }
// common stuff
  // if targeted you dont need many frames 
  if(my_sz->targeted_on==0){
     iw=seek_word(word,"NFRAMES");
     if(iw>=0){
         sscanf(word[iw+1],"%i", &my_sz->number); 
     } else {
         fprintf(fplog,"|- NEEDED \"NFRAMES\" KEYWORD FOR PATH\n"); 
         path_help=1; 
     }
  } else { 
     my_sz->number=1;
  }
  // if targeted  lambda can be simply 1.0 
  if(my_sz->targeted_on==0){
     iw=seek_word(word,"LAMBDA");
     if(iw>=0){ 
         sscanf(word[iw+1],"%lf", &lambda);
#ifdef STANDALONE
         lambda/=(mtd_data.ampli)*(mtd_data.ampli);
#endif
     } else {
         fprintf(fplog,"|- NEEDED \"LAMBDA\" KEYWORD FOR PATH\n"); 
         path_help=1; 
     }
  } else { 
     lambda=1.0;
  }
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &sigma);
             colvar.delta_r[count]  = (real) sigma; }
// rmsd parameters
  iw=seek_word(word,"FRAMESET");
  if(iw>=0) sscanf(word[iw+1],"%s", my_sz->names); 
// cmap parameters
  iw = seek_word(word,"INDEX");
  if(iw>=0) sscanf(word[iw+1],"%s", file_maps);
  iw = seek_word(word,"MAP");
  if(iw>=0) sscanf(word[iw+1],"%s", file_maps2);
  iw=seek_word(word,"GROUP");
  if(iw>=0) {
   sscanf(word[iw+1],"%s",file_group);
   my_sz->my_cmap_pack.logical_group = 1;
   }
// UMBRELLA stuff
  iw=seek_word(word,"UMB_LAG");
  if(iw>=0) {
    sscanf(word[iw+1],"%d",&(my_sz->umblagsteps));
    my_sz->umb_on=1; 
    my_sz->countperm=0;
  }
  iw=seek_word(word,"UMB_BLOCK");
  if(iw>=0) sscanf(word[iw+1],"%d",&(my_sz->umbblocksize)); 
  iw=seek_word(word,"UMB_STRIDE");
  if(iw>=0) sscanf(word[iw+1],"%d",&(my_sz->umbstride));
  iw=seek_word(word,"UMB_PERM");
  if(iw>=0) sscanf(word[iw+1],"%d",&(my_sz->umbpermanency));
  iw=seek_word(word,"UMB_TOL");
  if(iw>=0) sscanf(word[iw+1],"%lf",&tol);

// NEIGHBOUR LIST STUFF
  my_sz->neigh=0;
  my_sz->lneigh=(int *)malloc( MAXFRAMES_PATH * sizeof(int)); 
  iw=seek_word(word,"NEIGHLIST");  
  if(iw>=0) {
       sscanf(word[iw+1],"%d",&(my_sz->neigh_time)); // time for list 
       sscanf(word[iw+2],"%d",&(my_sz->nneigh));     // number of neighbours 
       if(my_sz->nneigh>=my_sz->number){
              my_sz->nneigh=my_sz->number;
              my_sz->neigh=0;
              for(i=0;i<my_sz->nneigh;i++)my_sz->lneigh[i]=i;
              neigh=1;
       }else {   
              my_sz->neigh=1;
              for(i=0;i<my_sz->nneigh;i++)my_sz->lneigh[i]=i;
              neigh=2;
       } 
 
  } else {
       my_sz->neigh=0;
       my_sz->nneigh=my_sz->number;
       for(i=0;i<my_sz->nneigh;i++)my_sz->lneigh[i]=i; 
  } 


  my_sz->lambda  = (real) lambda; 
  my_sz->umbtolerance = (real) tol;

  if(colvar.type_s[count]==30) strcpy(type,"S");
  if(colvar.type_s[count]==31) strcpy(type,"Z");

  if(path_help){
         fprintf(fplog,"|- PATH/TARGETED SYNTAX:\n");
         fprintf(fplog,"|- TYPE               : can be RMSD/CMAP/DRMSD\n");
         fprintf(fplog,"|- NFRAMES            : the number of reference structures\n");
         fprintf(fplog,"|- LAMBDA             : the common prefactor in the exponential equation for the path\n");
         fprintf(fplog,"|- SIGMA              : hills width for this cv\n");
         fprintf(fplog,"|- NEIGHLIST (opt.)   : neighlist on the closest frames    NEIGHLIST (ntimesteps) (nframes) \n");
         fprintf(fplog,"|-     \n");
         fprintf(fplog,"|- .....many other keywords are specific to RMSD/CMAP/DRMSD paths... \n");
         fprintf(fplog,"|- (RMSD)   FRAMESET :base for frameset name pippo_ will look for pippo_1.pdb pippo_2.pdb etc...  \n");
         fprintf(fplog,"|-     \n");
         fprintf(fplog,"|- e.g.\n");
         fprintf(fplog,"|- \n");
         fprintf(fplog,"|- Z_PATH   TYPE RMSD FRAMESET frame_ NFRAMES 1 LAMBDA 1.0 SIGMA 1.0 NEIGHLIST 10 5 \n");
         fprintf(fplog,"|- \n");
         fprintf(fplog,"|- S_PATH   TYPE RMSD FRAMESET frame_ NFRAMES 1 LAMBDA 1.0 SIGMA 1.0 NEIGHLIST 10 5 \n");
         fprintf(fplog,"|- \n");
         fprintf(fplog,"|- TARGETED TYPE RMSD FRAMESET frame.pdb SIGMA 1.0 \n");
         fprintf(fplog,"|- \n");
         EXIT();
  } 
  if(strcmp(my_sz->path_type,"CMAP") != 0 && strcmp(my_sz->path_type,"RMSD") != 0
     && strcmp(my_sz->path_type,"DRMS") !=0 && strcmp(my_sz->path_type,"HYBRID") !=0 ){
    fprintf(fplog,"ERROR: %s Unknown type of path!! \n",my_sz->path_type);   
    EXIT();
  }
  fprintf(fplog,"\n%1i-%s_PATH in %s space \n",count+1,type,my_sz->path_type);
  fprintf(fplog,"|--NFRAMES %i LAMBDA %f ",my_sz->number,my_sz->lambda);
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  if(strcmp(my_sz->path_type,"CMAP") == 0){
   fprintf(fplog,"|--READING CONTACT MAPS INDEX FROM FILE %s AND VALUES FROM %s\n",file_maps,file_maps2); 
   if(my_sz->my_cmap_pack.logical_group) fprintf(fplog,"|--AND GROUP FROM FILE %s \n",file_group);
  }
  if(strcmp(my_sz->path_type,"RMSD") == 0 || strcmp(my_sz->path_type,"DRMS") == 0  || strcmp(my_sz->path_type,"HYBRID") == 0 ){
   fprintf(fplog,"|--BASENAME FOR FRAMES %s \n",my_sz->names);
  }
  if(neigh==0) fprintf(fplog,"|--NEIGHBOUR LIST OFF \n");  
  if(neigh==1) fprintf(fplog,"|--NO NEED FOR NEIGHBOUR LIST: THE LIST APPEARS TOO LARGE \n");  
  if(neigh==2) {fprintf(fplog,"|--NEIGHBOUR LIST ON : TIME FOR LIST: %d TSTEPS \n", my_sz->neigh_time);
                fprintf(fplog,"|--                  : LIST LENGTH  : %d FRAMES \n", my_sz->nneigh); } 

  
  if(my_sz->number>MAXFRAMES_PATH){
                fprintf(fplog,"ERROR: Exceeded maximum number of frames in path CV\n");
                fprintf(fplog,"ERROR: increase MAXFRAMES_PATH and recompile \n");
                EXIT();
           }
  if(my_sz->number==0){
                   fprintf(fplog,"ERROR: NUMBER NOT FOUND\n ");
                   fprintf(fplog,"SYNTAX:  NAME NUMBER LAMBDA SIGMA\n ");
                   EXIT();
           }
   if(my_sz->lambda==0.){
                   fprintf(fplog,"ERROR: LAMBDA NOT FOUND\n ");
                   fprintf(fplog,"SYNTAX:  NAME NUMBER LAMBDA SIGMA\n ");
                   EXIT();
           }

   if(my_sz->umb_on){
    fprintf(fplog,"|--MEAN EVALUATION ON %s IS ACTIVE\n",type);
    fprintf(fplog,"|---LAGSTEPS %d - BLOCKSTEPS %d  - STRIDESTEPS %d  \n",my_sz->umblagsteps,my_sz->umbblocksize,my_sz->umbstride);
    fprintf(fplog,"|---PERMSTEPS %d - TOLERANCESTEPS %lf \n",my_sz->umbpermanency,my_sz->umbtolerance);
   }


// case of CMAP path

 if(strcmp(my_sz->path_type,"CMAP") == 0){

   read_sz_map(my_sz,file_maps,file_maps2,file_group,fplog);

   colvar.natoms[count]   = my_sz->my_cmap_pack.atoms; 

   snew(colvar.myder[count], colvar.natoms[count]);
   snew(colvar.cvatoms[count], colvar.natoms[count]);

   for(i=0;i<colvar.natoms[count];i++){
     colvar.cvatoms[count][i] = my_sz->my_cmap_pack.list[i];
   }
 
   if(my_sz->umb_on){
     my_sz->umb_map_block=float_3d_array_alloc(my_sz->my_cmap_pack.number+my_sz->my_cmap_pack.gnumber,my_sz->number,my_sz->umbblocksize);
     my_sz->umb_map_avg=float_2d_array_alloc(my_sz->my_cmap_pack.number+my_sz->my_cmap_pack.gnumber,my_sz->number);
     fprintf(fplog,"|---ALLOCATED for UMBRELLA STUFF %lf Kbytes \n",sizeof(real)*((real) ((my_sz->my_cmap_pack.number+
                  my_sz->my_cmap_pack.gnumber)*my_sz->number*my_sz->umbblocksize)/1000));
   }
 }

// RMSD/DRMS case 

 if(strcmp(my_sz->path_type,"RMSD") == 0 || strcmp(my_sz->path_type,"DRMS") == 0 ){

   read_sz_rmsd(my_sz,fplog);

   colvar.natoms[count]   = (my_sz->frameset[0])->natoms; 

   snew(colvar.myder[count], colvar.natoms[count]);
   snew(colvar.cvatoms[count], colvar.natoms[count]);

   for(i=0;i<colvar.natoms[count];i++){
     colvar.cvatoms[count][i] = ((my_sz->frameset[0])->atmnum[i])-1; 
   }

   if(my_sz->umb_on){
     my_sz->umb_block=float_4d_array_alloc(3,(*my_sz->frameset[0]).natoms,my_sz->number,my_sz->umbblocksize);
     my_sz->umb_avg=float_3d_array_alloc(3,(*my_sz->frameset[0]).natoms,my_sz->number);
     fprintf(fplog,"|---ALLOCATED for UMBRELLA STUFF %lf Kbytes",sizeof(real)*((real) 3*((*my_sz->frameset[0]).natoms)*(my_sz->number)*(my_sz->umbblocksize)/1000));
   }
 }
// HYBRID STRUCTURE: uses the definition of previous cvs
 if(strcmp(my_sz->path_type,"HYBRID") == 0 ){
    // read which cvs you want to hybridize  
     iw=seek_word(word,"HYBRID");
     if(iw>=0) {
        int c,ii,*t1l,t1s=0, *t2l,t2s,jj;
        ii=0;
        while (1) { 
             //      word[iw+1] is a word ???   
             ii++;
             c=*word[iw+ii];
             if (isalpha(c)!=0 ){
                 //printf("This is a word %s\n",word[iw+ii]); 
                 break; 
             } else {
       //          printf("This is a number %s\n",word[iw+ii]); 
                 if(t1s==0){ t1l=(int *)malloc(sizeof(int));
                             sscanf(word[iw+ii],"%d", &t1l[0]);
                 } else{ 
                        t2s=t1s;t2l=(int *)malloc(t2s*sizeof(int));
                        for(jj=0;jj<t1s;jj++){
                              t2l[jj]=t1l[jj];
                       }
                       free(t1l);
                       t1l=(int *)malloc((t1s+1)*sizeof(int)); 
                       for(jj=0;jj<t1s;jj++){
                              t1l[jj]=t2l[jj];
                       }
                       free(t2l);
                       sscanf(word[iw+ii],"%d", &t1l[t1s]);
                 } 
                 t1s++;
                 //for(jj=0;jj<t1s;jj++){printf("LL %d ",t1l[jj]);}; printf(" NN %d\n",t1s) ; 
             } 
        } 
        // allocate and copy a specific structure
        my_sz->nhybrid=t1s;
        my_sz->lhybrid=(int *)malloc(my_sz->nhybrid*sizeof(int));
        my_sz->lcvhybrid=(int *)malloc(my_sz->nhybrid*sizeof(int));
        for(ii=0;ii<my_sz->nhybrid;ii++){my_sz->lhybrid[ii]=t1l[ii];}  
        //for(jj=0;jj<my_sz->nhybrid;jj++){printf("LL %d ",my_sz->lhybrid[jj]);}; printf(" NN %d\n",my_sz->nhybrid) ; 
        printf("|--NUMBER OF HYBRID CVS= %d WHICH ARE ",my_sz->nhybrid);
        for(jj=0;jj<my_sz->nhybrid;jj++){
              printf("CV %d ",my_sz->lhybrid[jj]);my_sz->lhybrid[jj]--;
              if(my_sz->lhybrid[jj]>=count){printf("|--WRONG INPUT: PUT THE HYBRID CV AFTER THE ONES YOU USE TO MAKE THE HYBRID\n");EXIT();}
        };printf("\n"); 
        // set the type of the representation for each cv  
        for(jj=0;jj<my_sz->nhybrid;jj++){
            int kk=my_sz->lhybrid[jj];
            printf("CVTYPE %d\n",colvar.type_s[kk]);my_sz->lcvhybrid[jj]=colvar.type_s[kk];
            // parse the input and check if it is implemented   
            switch(my_sz->lcvhybrid[jj]){
              case 1:  printf("|- IMPLEMENTED\n") ; break;
              case 3:  printf("|- IMPLEMENTED\n") ; break;
              case 4:  printf("|- IMPLEMENTED\n") ; break;
              case 5:  printf("|- IMPLEMENTED\n") ; break;
              case 31:  printf("|- IMPLEMENTED\n") ; break;
              default:  printf("|- NOT IMPLEMENTED: now dying...\n"); EXIT();  
            }  
        }; 
        // matrix of metrics 
        my_sz->mathybrid=float_2d_array_alloc(my_sz->nhybrid,my_sz->nhybrid); 
        colvar.natoms[count]=read_sz_hybrid(my_sz,fplog); 
        snew(colvar.myder[count], colvar.natoms[count]);
        snew(colvar.cvatoms[count], colvar.natoms[count]);

        // assuming they are all the same for all the framesets (generally it is) 
        for(i=0;i<colvar.natoms[count];i++){
          colvar.cvatoms[count][i] = my_sz->hbd_frameset[0]->backtable[i]  ;
        }
         
     } else {
         printf("SOMETHING WENT WRONG \n");
     }  
//     EXIT(); 
 }
           nsz++;

           if(nsz==NMAX_PATH){
                printf("read_restraint.C: TOO MANY PATH CVS. INCREASE NMAX_PATH and recompile\n");
                EXIT();
           }


         fprintf(fplog,"\n");
         return colvar.natoms[count];
}


// ------------------------------------------------------------------------------------------------

int PREFIX read_sz_rmsd(struct sz_data *my_sz, FILE *fplog) {

        int l,i,j,k,found;
        char *str,ic[3],str2[100];
        l=0;
/*
  * allocate the pointers
 */
        my_sz->frameset=(struct coordinates_frameset **)malloc((my_sz->number)*sizeof(struct coordinates_frameset *));
        for(i=0;i< my_sz->number;i++){
                my_sz->frameset[i]=(struct coordinates_frameset *)malloc(sizeof(struct coordinates_frameset)) ;
        }

/*
 *  build names
 */
        str=&(my_sz->names[0]);
        for (i=1;i<= my_sz->number ;i++){
                strcpy(str2,my_sz->names);
                if(my_sz->targeted_on==0){
                     if(i<10){
                      ic[0]='0'+i;
                      ic[1]='\0';}
                     else if(i<100) {
                      ic[0]='0'+i/10 ;
                      ic[1]='0'+i%10 ;
                      ic[2]='\0';
                     }
                     else{
                       fprintf(fplog,"|--read_sz_input: TOO MANY FRAMES REQUIRED FOR NAME BUILDING!\n");
                       EXIT();
                     }
                     strcat(str2,ic);
                     strcat(str2,".pdb");
                }
                fprintf(fplog,"|--%s\n",str2);
                read_sz_coord(str2,my_sz->frameset[i-1],fplog);
        }

/*
 * check over that aligned atoms are the same ( this is requirement for rmsd routine
 * so you don't have to reallocate everything everytime )
 */
       for (i=0;i< my_sz->number-1;i++){
        for (j=i+1;j< my_sz->number;j++){
         if( ( *my_sz->frameset[i] ) .nalign == ( *my_sz->frameset[j] ).nalign ){
            for(k=0;k<(*my_sz->frameset[i]).nalign;k++){
               found=0;
               for(l=0;l<(*my_sz->frameset[j]).nalign;l++){
                  if( (*my_sz->frameset[i]).align_to_frameset[k]==(*my_sz->frameset[j]).align_to_frameset[l] ){found++;}
               }
               if(found==0){fprintf(fplog,"|--ERROR: ATOM %d in frameset %d not found in frameset %d\n",(*my_sz->frameset[i]).align_to_frameset[k],i,j);EXIT();}
               else if(found>1){fprintf(fplog,"|--ERROR: found multiple definition of %d in frameset %d not found in frameset %d\n",(*my_sz->frameset[i]).align_to_frameset[k],i,j);EXIT();}            }
         }
         else{fprintf(fplog,"|--ERROR : ALIGNMENT ATOMS IN THE FRAMESET %d AND %d ARE NOT THE SAME\n",i,j);EXIT();};
        }
       }
// 
// now write the backtable align_to_coord
// 
       for(l=0;l<my_sz->number;l++){// for each frameset in the set
          j=0;
          for(i=0;i<(*my_sz->frameset[l]).natoms;i++){//on all the atoms
             if((*my_sz->frameset[l]).align[i]!=0.){ // only if this atom is used in the alignment
                  (*my_sz->frameset[l]).align_to_coord[j]=i;
                  j++;
             }
             (*my_sz->frameset[l]).frameset_to_coord[i]=i;
           }
       }
       return 0;
};

// ------------------------------------------------------------------------------------------------

int PREFIX read_sz_coord (char *filename, struct coordinates_frameset *p, FILE *fplog){


    char string[400],sm[10],sm1[10];
    char *str,remark[10],end[10],atom[5],chain[3];
    real tmp0;
    FILE *fp;
    int i,l;
    double x,y,z;

    fp=fopen(filename,"r");
    if (fp == NULL)
     {fprintf(fplog,"|--UNSUCCESSFULL OPENING FOR FILE %s \n",filename);EXIT(); };
/* PDB
 
         1         2         3         4         5         6         7         8

12345678901234567890123456789012345678901234567890123456789012345678901234567890

ATOM   4150  H   ALA A 431       8.674  16.036  12.858  1.00  0.00

*/
l=0;
p->ndisplace=0;
p->nalign=0;
while(1){
  readagain:
  str=fgets(string,100,fp);
  if(str==NULL)break;
  if (feof(fp))break;
  //ATOM FIELD       
  sscanf(str,"%3s",end);if(strstr(end,"END")!=NULL){break;};
  sscanf(str,"%6s",remark);if(strstr(remark,"REMARK")!=NULL){goto readagain;};
  //ATOM     11  HA2 GLY    2      13.693  12.686  -2.859  1.00  0.00      PRO
  //ATOM     46  C          1       3.683  -2.464   1.429  1.00  1.00
  sscanf(str,"%s %d %s %s %d %lf %lf %lf %s %s",atom,&(p->atmnum[l]),(p->label[l]),(p->resname[l]),&(p->resid[l]),&x,&y,&z,sm,sm1);
  
  p->pos[l][0] = (real) x;
  p->pos[l][1] = (real) y;
  p->pos[l][2] = (real) z;
  
  #if defined (PLUMED_GROMACS)
  p->pos[l][0] /= 10.;
  p->pos[l][1] /= 10.;
  p->pos[l][2] /= 10.;
  #endif
  
  #if defined (STANDALONE)
  p->pos[l][0] *= mtd_data.ampli;
  p->pos[l][1] *= mtd_data.ampli;
  p->pos[l][2] *= mtd_data.ampli;
  #endif
  
  //fprintf(fplog,"RESID NUM %d  RESID %d  RESNAME %s LABEL %s PX %f PY %f PZ %f\n  ",p->atmnum[l],p->resid[l],p->resname[l],p->label[l],p->pos[l][0],p->pos[l][1],p->pos[l][2]);
  //alignment
  tmp0=atof(sm);
  if(tmp0==0.){
    p->align[l]=0.;
  }
  else{
    p->align[l]=tmp0;
    p->align_to_frameset[p->nalign]=l;
    p->nalign++;
  }
  //displacement
  tmp0=atof(sm1);
  if(tmp0==0.){
    p->displace[l]=0;
  }
  else{
    p->displace[l]=1;
    //p->displace[l]=tmp0;
    p->ndisplace++;
  }
  l++;
  if(l>=MAXATOMS_PATH){
          fprintf(fplog,"|--read_coord: Exceeded number of atoms per frame in path CV.\n");
          fprintf(fplog,"|--read_coord: Increase MAXATOMS_PATH and recompile.\n");
  }
}
fclose(fp);
p->natoms=l;
fprintf(fplog,"|---FOUND %d ATOMS FOR DISPLACEMENT\n",p->ndisplace);
fprintf(fplog,"|---FOUND %d ATOMS FOR ALIGNMENT\n",p->nalign);

if(p->nalign==0){
                   fprintf(fplog,"IT SEEMS YOU DONT  WANT TO ALIGN ANY ATOM\n");
                    fprintf(fplog,"Your frameset should look like:\n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"                                                         the    the    \n");
                    fprintf(fplog,"                                                         align  mea      \n");
                    fprintf(fplog,"                                                         ment   sure     \n");
                    fprintf(fplog,"                                                          |     |  \n");
                    fprintf(fplog,"                                                          V     V  \n");
                    fprintf(fplog,"ATOM   4150  H   ALA A 431       8.674  16.036  12.858  1.00  0.20\n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"\n");
                   EXIT();
                }
if(p->ndisplace==0){
                    fprintf(fplog,"IT SEEMS YOU DONT  WANT TO MEASURE THE DISPLACEMENT OF ANY ATOM\n");
                    fprintf(fplog,"Your frameset should look like:\n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"                                                         the    the    \n");
                    fprintf(fplog,"                                                         align  mea      \n");
                    fprintf(fplog,"                                                         ment   sure     \n");
                    fprintf(fplog,"                                                          |     |  \n");
                    fprintf(fplog,"                                                          V     V  \n");
                    fprintf(fplog,"ATOM   4150  H   ALA A 431       8.674  16.036  12.858  1.00  0.20\n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog," .      .    .    .  .  .          .       .       .      .     .  \n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"\n");
                    fprintf(fplog,"\n");
                    EXIT();
                   }
// find the back index ( from alignment to frameset )
l=0;
 for(i=0;i< p->natoms;i++){
   if (p->align[i]==0.){p->frameset_to_align[i]=-1 ; } // if negative there's no position in the rmsd frame
   else {p->frameset_to_align[i]=l;l++;}; // if >=0 then provides the index
 }
return 0;
};

// ------------------------------------------------------------------------------------------------

void PREFIX read_sz_map(struct sz_data *my_sz, char file_maps[129], char file_maps2[129],
                        char file_group[129], FILE *fplog){

 FILE   *fmap;
 int    i,ii,kk,jj,kkk;
 char   stringa[200],end[10],tmpstring[8],tmp2string[8],*cmstring;
 double r0,cutoff,map,weight;

// open INDEX file
   fmap=fopen(file_maps,"r");
   if(fmap == NULL){
               fprintf(fplog,"Error opening %s for reading. Program terminated. \n", file_maps);
               EXIT();
       }
   kk=0;
   ii=0;
   jj=0;
   while(1){
         cmstring=fgets(stringa,200,fmap);
         if(cmstring==NULL){break;}
         sscanf(cmstring,"%3s",end);if(strstr(end,"END")!=NULL){break;};
         sscanf(cmstring,"%7s",tmpstring);if(strstr(tmpstring,"CONTACT")!=NULL){ii=ii+1;}
         sscanf(cmstring,"%5s",tmp2string);if(strstr(tmp2string,"GROUP")!=NULL){jj=jj+1;}
         sscanf(cmstring,"%7s %d %d %d %lf %d %d %lf %lf",tmpstring,&kkk,&(my_sz->my_cmap_pack.index1[kk]),
                &(my_sz->my_cmap_pack.index2[kk]),&r0,&(my_sz->my_cmap_pack.nn[kk]),
                &(my_sz->my_cmap_pack.nd[kk]),&cutoff,&weight);


         my_sz->my_cmap_pack.index1[kk]--;
         my_sz->my_cmap_pack.index2[kk]--;
         my_sz->my_cmap_pack.r0[kk] = (real) r0;
         my_sz->my_cmap_pack.cutoff[kk] = (real) cutoff;
         my_sz->my_cmap_pack.weight[kk] = (real) weight;
         kk=kk+1;
         

         }

   my_sz->my_cmap_pack.number=ii;
   my_sz->my_cmap_pack.gnumber=jj;

   fprintf(fplog,"|--%d atomic contacts / %d group contacts \n",my_sz->my_cmap_pack.number,my_sz->my_cmap_pack.gnumber);
   fclose(fmap);

   int iii=0;
   fmap=fopen(file_maps2,"r");
   if(fmap == NULL){
     fprintf(fplog,"Error opening %s for reading. Program terminated. \n", file_maps2);
     EXIT();
       }

   while(1){
    kk=0;
    while(1){
      cmstring=fgets(stringa,200,fmap);
      if(cmstring==NULL){break;}
      sscanf(cmstring,"%3s",end);if(strstr(end,"END")!=NULL){break;};
      sscanf(cmstring,"%d %d %d %lf",&kkk,&ii,&jj,&map);
 
      my_sz->my_cmap_pack.cmap[iii][kk] = (real) map;
 
      kk=kk+1;
    }

    if(cmstring==NULL){break;}
    iii=iii+1;
    }
    if(iii!=my_sz->number){
      fprintf(fplog,"ERROR: NUMBER OF FRAMES FOUND %d IS DIFFERENT FROM EXPECTED %d \n",iii,my_sz->number);
      EXIT();
    }
    fclose(fmap);

    if(my_sz->my_cmap_pack.logical_group){
      fmap=fopen(file_group,"r");
      if(fmap == NULL){
       fprintf(fplog,"Error opening %s for reading. Program terminated. \n", file_group);
       EXIT();
      }
      while(1){
       iii=0;
       cmstring=fgets(stringa,200,fmap);
       if(cmstring==NULL){break;}
       char *result = NULL;
       result = strtok( stringa, " \t" );
       while( result != NULL ) {
        if(iii==1) jj = atoi (result);
        if(iii==2) my_sz->my_cmap_pack.group.numatom[jj-1] = atoi (result);
        if(iii>2)  my_sz->my_cmap_pack.group.index[jj-1][iii-3] = atoi (result) - 1; 
        result = strtok( NULL, " \t" );
        iii=iii+1;
       }
      }
      my_sz->my_cmap_pack.group.number=jj;
      fclose(fmap);

      printf("|--Reading GROUP file. Found %d groups \n",  my_sz->my_cmap_pack.group.number);
      for(i=0;i<my_sz->my_cmap_pack.group.number;i++){
       fprintf(fplog,"|--GROUP %d #ofATOMS %d \n",i+1,my_sz->my_cmap_pack.group.numatom[i]);
       for(ii=0;ii<my_sz->my_cmap_pack.group.numatom[i];ii++){
        fprintf(fplog,"|---AT %d \n", my_sz->my_cmap_pack.group.index[i][ii]+1);
       }
      }
     } 
     
// Calculating number of atoms involved and creating the list
// from ATOMIC contacts
            ii=0;
            if(my_sz->my_cmap_pack.number!=0){my_sz->my_cmap_pack.list[ii]=my_sz->my_cmap_pack.index1[0];}

            for(iii=1;iii<my_sz->my_cmap_pack.number;iii++){
             int flag=0;
             for(kk=0;kk<=ii;kk++){
              if(my_sz->my_cmap_pack.index1[iii]==my_sz->my_cmap_pack.list[kk]){flag=1;}
             }
             if(flag==0){
             ii=ii+1;
             my_sz->my_cmap_pack.list[ii]=my_sz->my_cmap_pack.index1[iii];
             }
            }

           for(iii=0;iii<my_sz->my_cmap_pack.number;iii++){
             int flag=0;
             for(kk=0;kk<=ii;kk++){
              if(my_sz->my_cmap_pack.index2[iii]==my_sz->my_cmap_pack.list[kk]){flag=1;}
             }
             if(flag==0){
             ii=ii+1;
             my_sz->my_cmap_pack.list[ii]=my_sz->my_cmap_pack.index2[iii];
             }
            }

            if(my_sz->my_cmap_pack.number!=0){my_sz->my_cmap_pack.atoms=ii+1;}
            else{my_sz->my_cmap_pack.atoms=0;}
// From group

 if(my_sz->my_cmap_pack.logical_group){
           for(i=0;i<my_sz->my_cmap_pack.group.number;i++){
            for(ii=0;ii<my_sz->my_cmap_pack.group.numatom[i];ii++){
             kkk=my_sz->my_cmap_pack.group.index[i][ii];
             int flag=0;
                for(kk=0;kk<my_sz->my_cmap_pack.atoms;kk++){
                 if(my_sz->my_cmap_pack.list[kk]==kkk){
                  my_sz->my_cmap_pack.group.index_to_list[i][ii]=kk;
                  flag=1;
                 }
                }
              if(flag==0){
               my_sz->my_cmap_pack.list[my_sz->my_cmap_pack.atoms]=kkk;
               my_sz->my_cmap_pack.group.index_to_list[i][ii]=my_sz->my_cmap_pack.atoms;
               my_sz->my_cmap_pack.atoms=my_sz->my_cmap_pack.atoms+1;
              }
            }
           }
  }


            fprintf(fplog,"|--Total number of atoms involved %d \n",my_sz->my_cmap_pack.atoms);

// creating connection between my_cmap_pack.list e my_cmap_pack.index_from

           for(iii=0;iii<my_sz->my_cmap_pack.number;iii++){
            for(i=0;i<my_sz->my_cmap_pack.atoms;i++){
             if(my_sz->my_cmap_pack.index1[iii]==my_sz->my_cmap_pack.list[i]){my_sz->my_cmap_pack.index_from1[iii]=i;}
             if(my_sz->my_cmap_pack.index2[iii]==my_sz->my_cmap_pack.list[i]){my_sz->my_cmap_pack.index_from2[iii]=i;}
            }
           }

    return;

}

// ------------------------------------------------------------------------------------------------

void PREFIX cmap_running(struct cmap_inpack *inpack, struct cmap_pack *my_cmap_pack){

       int mm,nn,i,j;
       real dist,xp,xq;

// atomic contact

       for(i=0;i<my_cmap_pack->number;i++){

        mm=my_cmap_pack->index_from1[i];
        nn=my_cmap_pack->index_from2[i];

        dist=sqrt(pow2(inpack->r0[mm][0]-inpack->r0[nn][0])+
                  pow2(inpack->r0[mm][1]-inpack->r0[nn][1])+
                  pow2(inpack->r0[mm][2]-inpack->r0[nn][2]));

        if (dist>my_cmap_pack->cutoff[i] || my_cmap_pack->weight[i]==0){
             inpack->cmap[i]=0.;
        }
        else{
             if (fabs(dist-my_cmap_pack->r0[i])==0.0){
              inpack->cmap[i]=(real) my_cmap_pack->nn[i]/my_cmap_pack->nd[i];
             } else { 
              power(dist/my_cmap_pack->r0[i],my_cmap_pack->nn[i],my_cmap_pack->nd[i],&xp,&xq);
              inpack->cmap[i]=(1.-xp)/(1.-xq)*my_cmap_pack->weight[i];
             }
        }
       }

       if(my_cmap_pack->logical_group){
// group contacts
// evaluating center of mass

        for(i=0;i<my_cmap_pack->group.number;i++){

         my_cmap_pack->group.rcm[i][0]=0.;
         my_cmap_pack->group.rcm[i][1]=0.;
         my_cmap_pack->group.rcm[i][2]=0.;

         for (j=0;j<my_cmap_pack->group.numatom[i];j++){
          mm=my_cmap_pack->group.index_to_list[i][j];
          my_cmap_pack->group.rcm[i][0] += inpack->r0[mm][0];
          my_cmap_pack->group.rcm[i][1] += inpack->r0[mm][1];
          my_cmap_pack->group.rcm[i][2] += inpack->r0[mm][2];
         }

         my_cmap_pack->group.rcm[i][0] /= (real) my_cmap_pack->group.numatom[i];
         my_cmap_pack->group.rcm[i][1] /= (real) my_cmap_pack->group.numatom[i];
         my_cmap_pack->group.rcm[i][2] /= (real) my_cmap_pack->group.numatom[i];

       }


// evaluating contacts

        for(i=0;i<my_cmap_pack->gnumber;i++){

         j=i+my_cmap_pack->number;
         mm=my_cmap_pack->index1[j];
         nn=my_cmap_pack->index2[j];

         dist=sqrt(pow2(my_cmap_pack->group.rcm[mm][0]-my_cmap_pack->group.rcm[nn][0])+
                   pow2(my_cmap_pack->group.rcm[mm][1]-my_cmap_pack->group.rcm[nn][1])+
                   pow2(my_cmap_pack->group.rcm[mm][2]-my_cmap_pack->group.rcm[nn][2]));

         if (dist>my_cmap_pack->cutoff[j] || my_cmap_pack->weight[j]==0){
             inpack->cmap[j]=0.;
         }
         else{
            if(fabs(dist-my_cmap_pack->r0[j])==0.0){
             inpack->cmap[j] = (real) my_cmap_pack->nn[j]/my_cmap_pack->nd[j];
            } else {
              power(dist/my_cmap_pack->r0[j],my_cmap_pack->nn[j],my_cmap_pack->nd[j],&xp,&xq);
              inpack->cmap[j]=(1.-xp)/(1.-xq)*my_cmap_pack->weight[j];
            }
         }

        }
   }
 }

// ------------------------------------------------------------------------------------------------

void PREFIX cmdist_eval(int frame,struct cmap_inpack *inpack,struct cmap_outpack *outpack,
                  struct cmap_pack *my_cmap_pack,int dr1_calc){

       int    jj,k,i,j,ii,jjj,iii;
       int    tot;

       real tmp, dist_r0;
       real tmp4_r0_0,tmp4_r0_1,tmp4_r0_2;
       real tmp1_r0,tmp2_r0,tmp3_r0;
       real tmp1,tmp2,tmp3;
       real pow_P,pow_Q;
       real R01;
       int    P1,Q1;

       tmp=0.;
       tot=my_cmap_pack->number+my_cmap_pack->gnumber;

       for(i=0;i<tot;i++){
          tmp=tmp+pow2(my_cmap_pack->cmap[frame][i]-inpack->cmap[i]);
         }
        outpack->err=tmp;


                                      /* DERIVATIVE CALCULATION:respect to running frame and frameset */

// setting derivatives to zero

        for(k=0;k<my_cmap_pack->atoms;k++){
             outpack->derr_dr0[0][k]=0.;
             outpack->derr_dr0[1][k]=0.;
             outpack->derr_dr0[2][k]=0.;
        }


           for(j=0;j<my_cmap_pack->number;j++){
            if(my_cmap_pack->weight[j]!=0){
             ii=my_cmap_pack->index_from1[j];
             jj=my_cmap_pack->index_from2[j];

              dist_r0=sqrt(pow2(inpack->r0[ii][0]-inpack->r0[jj][0])+pow2(inpack->r0[ii][1]-inpack->r0[jj][1])+pow2(inpack->r0[ii][2]-inpack->r0[jj][2]));
              tmp1_r0=inpack->cmap[j]-my_cmap_pack->cmap[frame][j];
              R01=my_cmap_pack->r0[j];
              P1=my_cmap_pack->nn[j];
              Q1=my_cmap_pack->nd[j];

              power(dist_r0/R01,P1,Q1,&pow_P,&pow_Q);

              tmp2_r0=(Q1*pow_Q*(1.-pow_P)-P1*pow_P*(1.-pow_Q))*R01/dist_r0*my_cmap_pack->weight[j];
              tmp3_r0=R01*(1.-pow_Q)*(1.-pow_Q);

              tmp4_r0_0=(inpack->r0[ii][0]-inpack->r0[jj][0])/dist_r0;
              tmp4_r0_1=(inpack->r0[ii][1]-inpack->r0[jj][1])/dist_r0;
              tmp4_r0_2=(inpack->r0[ii][2]-inpack->r0[jj][2])/dist_r0;

              tmp1=2*tmp1_r0*tmp2_r0*tmp4_r0_0/tmp3_r0;
              tmp2=2*tmp1_r0*tmp2_r0*tmp4_r0_1/tmp3_r0;
              tmp3=2*tmp1_r0*tmp2_r0*tmp4_r0_2/tmp3_r0;

              outpack->derr_dr0[0][ii]+=tmp1;
              outpack->derr_dr0[1][ii]+=tmp2;
              outpack->derr_dr0[2][ii]+=tmp3;
              outpack->derr_dr0[0][jj]-=tmp1;
              outpack->derr_dr0[1][jj]-=tmp2;
              outpack->derr_dr0[2][jj]-=tmp3;
           }
          }
// case of group contact

          if(my_cmap_pack->logical_group){
           for(j=0;j<my_cmap_pack->gnumber;j++){

            i=j+my_cmap_pack->number;
            ii=my_cmap_pack->index1[i];
            jj=my_cmap_pack->index2[i];


              dist_r0=sqrt(pow2(my_cmap_pack->group.rcm[ii][0]-my_cmap_pack->group.rcm[jj][0])+
                           pow2(my_cmap_pack->group.rcm[ii][1]-my_cmap_pack->group.rcm[jj][1])+
                           pow2(my_cmap_pack->group.rcm[ii][2]-my_cmap_pack->group.rcm[jj][2]));

              tmp1_r0=inpack->cmap[i]-my_cmap_pack->cmap[frame][i];

              R01=my_cmap_pack->r0[i];
              P1=my_cmap_pack->nn[i];
              Q1=my_cmap_pack->nd[i];

              power(dist_r0/R01,P1,Q1,&pow_P,&pow_Q);

              tmp2_r0=(Q1*pow_Q*(1.-pow_P)-P1*pow_P*(1.-pow_Q))*R01/dist_r0*my_cmap_pack->weight[i];
              tmp3_r0=R01*(1.-pow_Q)*(1.-pow_Q);

              tmp4_r0_0=(my_cmap_pack->group.rcm[ii][0]-my_cmap_pack->group.rcm[jj][0])/dist_r0;
              tmp4_r0_1=(my_cmap_pack->group.rcm[ii][1]-my_cmap_pack->group.rcm[jj][1])/dist_r0;
              tmp4_r0_2=(my_cmap_pack->group.rcm[ii][2]-my_cmap_pack->group.rcm[jj][2])/dist_r0;

              tmp1=2*tmp1_r0*tmp2_r0*tmp4_r0_0/tmp3_r0;
              tmp2=2*tmp1_r0*tmp2_r0*tmp4_r0_1/tmp3_r0;
              tmp3=2*tmp1_r0*tmp2_r0*tmp4_r0_2/tmp3_r0;

              for(jjj=0;jjj<my_cmap_pack->group.numatom[ii];jjj++){
              iii=my_cmap_pack->group.index_to_list[ii][jjj];
              outpack->derr_dr0[0][iii]+=tmp1/((real) my_cmap_pack->group.numatom[ii]);
              outpack->derr_dr0[1][iii]+=tmp2/((real) my_cmap_pack->group.numatom[ii]);
              outpack->derr_dr0[2][iii]+=tmp3/((real) my_cmap_pack->group.numatom[ii]);
              }
              for(jjj=0;jjj<my_cmap_pack->group.numatom[jj];jjj++){
              iii=my_cmap_pack->group.index_to_list[jj][jjj];
              outpack->derr_dr0[0][iii]-=tmp1/((real) my_cmap_pack->group.numatom[jj]);
              outpack->derr_dr0[1][iii]-=tmp2/((real) my_cmap_pack->group.numatom[jj]);
              outpack->derr_dr0[2][iii]-=tmp3/((real) my_cmap_pack->group.numatom[jj]);
              }

           }
          }

           if(dr1_calc){
             for(i=0;i<tot;i++){
              outpack->derr_dcm[i]=-2.*(inpack->cmap[i]-my_cmap_pack->cmap[frame][i]);
             }
           }

}

// ------------------------------------------------------------------------------------------------

real PREFIX pow2(real x){
return x*x;
}

// ------------------------------------------------------------------------------------------------

void PREFIX power(real x,int p,int q,real *xp,real *xq){
     int i;
     real tot;

     tot=1;
     if(p>=q){
       for(i=1;i<=q;i++){
        tot=tot*x;
        }
        *xq=tot;
       for(i=q+1;i<=p;i++){
        tot=tot*x;
        }
        *xp=tot;}
     else{
       for(i=1;i<=p;i++){
        tot=tot*x;
        }
        *xp=tot;
       for(i=p+1;i<=q;i++){
        tot=tot*x;
        }
        *xq=tot;}
}

// ------------------------------------------------------------------------------------------------


void  PREFIX msd_calculation(struct coordinates_frameset *pframeset,struct cmap_inpack *c_inpack,
                             struct cmap_outpack *c_outpack,real dmsd_dr1[3][MAXATOMS_PATH]){

        int j,l,k,m,n,o;

        real tmp0,tmp1,ndisplace,nalign,displace,align;
        real coeff[3][MAXATOMS_PATH];
        struct rmsd_inpack inpack;
#ifdef RMSD_FULL
	struct rmsd_outpack outpack;
#else
	struct rmsd_mini_outpack outpack;
#endif

			/* number of atoms */
	inpack.natoms=pframeset->nalign;
                           /* transfer atoms from the current set */
	for(k=0;k< pframeset->nalign;k++){
                          l=pframeset->align_to_coord[k];
                          inpack.r0[0][k]=c_inpack->r0[l][0];        
                          inpack.r0[1][k]=c_inpack->r0[l][1];        
                          inpack.r0[2][k]=c_inpack->r0[l][2];        
                         // fprintf(mtd_data.fplog,"COORD0 %d  %d %12.6f %12.6f %12.6f\n",k,l,inpack.r0[0][k],inpack.r0[1][k],inpack.r0[2][k]);
        }

                        /* transfer the atoms in the frameset */
        nalign=0.; 
        for(k=0;k<pframeset->nalign;k++){
                          l=pframeset->align_to_frameset[k];
                          inpack.r1[0][k]=pframeset->pos[l][0];
                          inpack.r1[1][k]=pframeset->pos[l][1];
                          inpack.r1[2][k]=pframeset->pos[l][2];
                          inpack.mass[k]=pframeset->align[l]; // this contains the weight in mass avg
                          nalign+=inpack.mass[k]; 
                        //  fprintf(mtd_data.fplog,"COORD1  %12.6f %12.6f %12.6f\n",inpack.r1[0][k],inpack.r1[1][k],inpack.r1[2][k]);
        }
#ifdef RMSD_FULL
        rmsd_pack(inpack,&outpack,7,1);
#else
        rmsd_mini_pack(inpack,&outpack,7);
#endif
        //cout<<"ERROR_RMSD_PACK "<<outpack.err<<endl;	
	//(*err)=outpack.err;

                            /* check rotation and translation */

        //                for(k=0;k<my_coord0.natoms;k++){
        //                     my_coord0.pos[0][k]-=outpack.cmr0[0];
        //                     my_coord0.pos[1][k]-=outpack.cmr0[1];
        //                     my_coord0.pos[2][k]-=outpack.cmr0[2];
	//		}

        //                ll=init_pdb("current.pdb"); 
        //                plot_pdb(ll,&my_coord0); 

        //                mm=init_pdb("reference.pdb"); 
        //                for(k=0;k<pframeset->natoms;k++){
        //                     pframeset->pos[0][k]-=outpack.cmr1[0];
        //                     pframeset->pos[1][k]-=outpack.cmr1[1];
        //                     pframeset->pos[2][k]-=outpack.cmr1[2];
	//		}
        //                for(k=0;k<pframeset->natoms;k++){

        //                     tmp1=pframeset->pos[0][k]*outpack.d[0][0]+
        //                          pframeset->pos[1][k]*outpack.d[0][1]+
        //                          pframeset->pos[2][k]*outpack.d[0][2];

        //                     tmp2=pframeset->pos[0][k]*outpack.d[1][0]+
        //                          pframeset->pos[1][k]*outpack.d[1][1]+
        //                          pframeset->pos[2][k]*outpack.d[1][2];

        //                     tmp3=pframeset->pos[0][k]*outpack.d[2][0]+
        //                          pframeset->pos[1][k]*outpack.d[2][1]+
        //                          pframeset->pos[2][k]*outpack.d[2][2];

	//	             pframeset->pos[0][k]=tmp1;
        //                     pframeset->pos[1][k]=tmp2;
        //                     pframeset->pos[2][k]=tmp3;
	//		}
        //                plot_pdb(mm,pframeset); 
        //                CkExit();
                            /* REAL MSD CALCULATION */

			         tmp0=0.;
                                 ndisplace=(real) pframeset->ndisplace; 
                                 // old:: nalign=double(pframeset->nalign); 

                                 for(k=0;k<pframeset->natoms;k++){
                                     for(l=0;l<3;l++){

                                        displace=(real) pframeset->displace[k]; 
                                        align=(real) pframeset->align[k]; 
                                        tmp1=0.;
                                 
                                        // contribution from rotated reference frame //
                                        for(m=0;m<3;m++){
                                           tmp1-=outpack.d[l][m]*(pframeset->pos[k][m]-outpack.cmr1[m]);
                                        }
                                 
                                        // contribution from running centered frame //
                                        j=pframeset->frameset_to_coord[k]; 
                                        tmp1+=(c_inpack->r0[j][l]-outpack.cmr0[l]); 
                                        
                                        tmp1*=displace;

                                        //printf("DISPLACED ATOM %d %f %f \n",k,tmp2,tmp1*tmp1);
                                        
                                        coeff[l][k]=tmp1;// store coefficents for derivative usage// 

                                        tmp0+=tmp1*tmp1; //squared distance added//
                                         
			             }
			         }  
                                 tmp0=tmp0/ndisplace;
                                 
			         //printf("ERRR NEW %f \n",tmp0);
                                 
                                 c_outpack->err=tmp0;

                           /* DERIVATIVE CALCULATION:respect to running frame */

                                for(k=0;k<pframeset->natoms;k++){
                                     for(l=0;l<3;l++){
                                            
                                         displace=(real) pframeset->displace[k]; 
                                         align=(real) pframeset->align[k]; 

                                         tmp1=0.;
					 
                                         tmp1 +=2.0*coeff[l][k]*displace/ndisplace ;
                                         for(o=0;o<pframeset->natoms;o++){
                                           tmp1 -=2.0*coeff[l][o]*align/(nalign*ndisplace); 
					 } 

				         j=pframeset->frameset_to_align[k]; //index of the k atom passed to the rmsd routine
                                         if(j>=0){
                                                 for(m=0;m<pframeset->natoms;m++){
                                                        displace=(real) pframeset->displace[m]; 
                                                     	for(n=0;n<3;n++){
                                                 		tmp0=0.;
                                                     		for(o=0;o<3;o++){
				 	       	                   tmp0+=outpack.dd_dr0[n][o][l][j]*(pframeset->pos[m][o]*displace-outpack.cmr1[o]);
				                        	}
                                                                tmp0*=2.0*coeff[n][m]/ndisplace;
                                                		tmp1-=tmp0;    
				 	       		}
		 		 	         }
					 }
                                         c_outpack->derr_dr0[l][k]=tmp1;
                                         //printf("DMSD %d %d %12.6f\n",l,k,dmsd_dr0[l][k]);
                                     }
                                }
				
                           /* DERIVATIVE CALCULATION:respect to frameset  */

                                for(k=0;k<pframeset->natoms;k++){
                                     for(l=0;l<3;l++){
                                            
                                         tmp1=0.;
					 
				         j=pframeset->frameset_to_align[k]; //index of the k atom passed to the rmsd routine

                                         if(j>=0){
                                                 for(m=0;m<pframeset->natoms;m++){
                                                        displace=(real) pframeset->displace[m]; 
                                                        align=(real) pframeset->align[m]; 
                                                     	for(n=0;n<3;n++){
                                                 		tmp0=0.;
                                                     		for(o=0;o<3;o++){
				 	       	                   tmp0+=outpack.dd_dr1[n][o][l][j]*
									  (pframeset->pos[m][o]-outpack.cmr1[o])*displace;
				                        	}
                                                                tmp0*=2.0*coeff[n][m]/ndisplace;
                                                		tmp1-=tmp0;    
				 	       		}
		 		 	         }
					 }

                                         displace=(real) pframeset->displace[k]; 
                                         align=(real) pframeset->align[k]; 
					 
					 tmp0=0.;
                                         for(o=0;o<3;o++){
					 	tmp0+=-(2.*coeff[o][k]/ndisplace)*displace*outpack.d[o][l];
                                         }
                                         tmp1+=tmp0;

                                         tmp0=0.;
                                         for(m=0;m<pframeset->natoms;m++){
                                         	for(o=0;o<3;o++){
							tmp0+=2.*coeff[o][m]*outpack.d[o][l]*align/(nalign*ndisplace);
						}
					 }
					 tmp1 += tmp0;

                                         dmsd_dr1[l][k]=tmp1;
                                     }
                                }
return;
};

// ------------------------------------------------------------------------------------------------


int PREFIX rmsd_pack(struct rmsd_inpack inpack,struct rmsd_outpack *outpack,int iopt,int iopt2)
{
/* declarations */
int i,j,k,l,p,ll,mm,nn,ii,ix,jx;
real rrsq,xx,yy,zz,m[4][4],rr1[4],rr0[4];
//cR double lambda[4],z[4][4],wk[20],s,q[4];
real lambda[4],s,q[4];
real dddq[3][3][4],gamma[3][3][3];
real dm_r1[4][4][3],dm_r0[4][4][3];
real dm_r1_store[4][4][3][MAXATOMS_RMSD];
real dm_r0_store[4][4][3][MAXATOMS_RMSD];
real derr_dr1_tmp[3][MAXATOMS_RMSD];
real derr_dr0_tmp[3][MAXATOMS_RMSD];
real dderr_dr1_dr1_tmp[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
real dderr_dr0_dr0_tmp[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
real dderr_dr1_dr0_tmp[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
real dderr_dr0_dr1_tmp[3][MAXATOMS_RMSD][3][MAXATOMS_RMSD];
real pi1[3][3],pi0[3][3]; 
real tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
//cR int ier,arg1,arg2,arg3;
real alpha_m1[3][3],alpha_m2[3][3],alpha_m3[3][3],alpha_m4[3][3];
xx=0.;
yy=0.;
zz=0.;
tmp1=0.;
if(iopt==5 || iopt == 7 ){
	for(i=0;i<inpack.natoms;i++){
		xx+=inpack.r0[0][i]*inpack.mass[i];
		yy+=inpack.r0[1][i]*inpack.mass[i];
		zz+=inpack.r0[2][i]*inpack.mass[i];
                tmp1+=inpack.mass[i];
	}
	xx=xx/((real) tmp1);
	yy=yy/((real) tmp1);
	zz=zz/((real) tmp1);
};
outpack->cmr0[0]=xx;
outpack->cmr0[1]=yy;
outpack->cmr0[2]=zz;
for(i=0;i<inpack.natoms;i++){
	outpack->r0p[0][i]=inpack.r0[0][i]-xx;
	outpack->r0p[1][i]=inpack.r0[1][i]-yy;
	outpack->r0p[2][i]=inpack.r0[2][i]-zz;
}
xx=0.;
yy=0.;
zz=0.;
tmp1=0.;
if(iopt==6 || iopt == 7 ){
	for(i=0;i<inpack.natoms;i++){
		xx+=inpack.r1[0][i]*inpack.mass[i];
		yy+=inpack.r1[1][i]*inpack.mass[i];
		zz+=inpack.r1[2][i]*inpack.mass[i];
                tmp1+=inpack.mass[i]; 
	};
	xx=xx/((real) tmp1);
	yy=yy/((real) tmp1);
	zz=zz/((real) tmp1);
};
outpack->cmr1[0]=xx;
outpack->cmr1[1]=yy;
outpack->cmr1[2]=zz;
for(i=0;i<inpack.natoms;i++){
	outpack->r1p[0][i]=inpack.r1[0][i]-xx;
	outpack->r1p[1][i]=inpack.r1[1][i]-yy;
	outpack->r1p[2][i]=inpack.r1[2][i]-zz;
}
// CLEAN M MATRIX
for(i=0;i<4;i++){
	for(j=0;j<4;j++){
          m[i][j]=0.;  
	}
}
// ASSIGN MATRIX ELEMENTS
for(i=0;i<inpack.natoms;i++){
	
        tmp1=sqrt(inpack.mass[i]); 
        rr1[0]=outpack->r1p[0][i]*tmp1;
        rr1[1]=outpack->r1p[1][i]*tmp1;
        rr1[2]=outpack->r1p[2][i]*tmp1;
        rr0[0]=outpack->r0p[0][i]*tmp1;
        rr0[1]=outpack->r0p[1][i]*tmp1;
        rr0[2]=outpack->r0p[2][i]*tmp1;
	
        rrsq=pow(rr0[0],2)+pow(rr0[1],2)+pow(rr0[2],2)+pow(rr1[0],2)+pow(rr1[1],2)+pow(rr1[2],2);
     
        m[0][0] +=  rrsq+2.*(-rr0[0]*rr1[0]-rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m[1][1] +=  rrsq+2.*(-rr0[0]*rr1[0]+rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m[2][2] +=  rrsq+2.*(+rr0[0]*rr1[0]-rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m[3][3] +=  rrsq+2.*(+rr0[0]*rr1[0]+rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m[0][1] += 2.*(-rr0[1]*rr1[2]+rr0[2]*rr1[1]);
        m[0][2] += 2.*( rr0[0]*rr1[2]-rr0[2]*rr1[0]);
        m[0][3] += 2.*(-rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m[1][2] -= 2.*( rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m[1][3] -= 2.*( rr0[0]*rr1[2]+rr0[2]*rr1[0]);
        m[2][3] -= 2.*( rr0[1]*rr1[2]+rr0[2]*rr1[1]);

};
m[1][0] = m[0][1];
m[2][0] = m[0][2];
m[2][1] = m[1][2];
m[3][0] = m[0][3];
m[3][1] = m[1][3];
m[3][2] = m[2][3];

// DIAGONALIZE 

ql77_driver(m,lambda);
s=1.0;
if(m[0][0]<0.)s=-1.;//correct for negative values (?)
q[0]=s*m[0][0];
q[1]=s*m[1][0];
q[2]=s*m[2][0];
q[3]=s*m[3][0];
outpack->err=sqrt(lambda[0]/((real) inpack.natoms));
if(lambda[0]==lambda[1]){
printf("DIAGONALIZATION: NON UNIQUE SOLUTION: ABORT \n");
EXIT();
}
if(iopt==0){return 0;}// JUST DIAGONALIZATION REQUIRED 

/*
 * Find the ROTATION matrix
 */
outpack->d[0][0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3]       ; 
outpack->d[1][0]=2.0*(q[1]*q[2]-q[0]*q[3]);
outpack->d[2][0]=2.0*(q[1]*q[3]+q[0]*q[2]);
outpack->d[0][1]=2.0*(q[1]*q[2]+q[0]*q[3]);
outpack->d[1][1]=q[0]*q[0]+q[2]*q[2]-q[1]*q[1]-q[3]*q[3];
outpack->d[2][1]=2.0*(q[2]*q[3]-q[0]*q[1]);
outpack->d[0][2]=2.0*(q[1]*q[3]-q[0]*q[2]);
outpack->d[1][2]=2.0*(q[2]*q[3]+q[0]*q[1]);
outpack->d[2][2]=q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
#ifdef EXTREME_DEBUG
for (i=0;i<3;i++){
printf("D_MATRIX %12.6f %12.6f %12.6f\n",outpack->d[i][0],outpack->d[i][1],outpack->d[i][2]);
}
#endif
/* 
 * first derivative in perturbation theory
 */
dddq[0][0][0]= 2.0*q[0];
dddq[1][0][0]=-2.0*q[3];
dddq[2][0][0]= 2.0*q[2];
dddq[0][1][0]= 2.0*q[3];
dddq[1][1][0]= 2.0*q[0];
dddq[2][1][0]=-2.0*q[1];
dddq[0][2][0]=-2.0*q[2];
dddq[1][2][0]= 2.0*q[1];
dddq[2][2][0]= 2.0*q[0];

dddq[0][0][1]= 2.0*q[1];
dddq[1][0][1]= 2.0*q[2];
dddq[2][0][1]= 2.0*q[3];
dddq[0][1][1]= 2.0*q[2];
dddq[1][1][1]=-2.0*q[1];
dddq[2][1][1]=-2.0*q[0];
dddq[0][2][1]= 2.0*q[3];
dddq[1][2][1]= 2.0*q[0];
dddq[2][2][1]=-2.0*q[1];

dddq[0][0][2]=-2.0*q[2];
dddq[1][0][2]= 2.0*q[1];
dddq[2][0][2]= 2.0*q[0];
dddq[0][1][2]= 2.0*q[1];
dddq[1][1][2]= 2.0*q[2];
dddq[2][1][2]= 2.0*q[3];
dddq[0][2][2]=-2.0*q[0];
dddq[1][2][2]= 2.0*q[3];
dddq[2][2][2]=-2.0*q[2];

dddq[0][0][3]=-2.0*q[3];
dddq[1][0][3]=-2.0*q[0];
dddq[2][0][3]= 2.0*q[1];
dddq[0][1][3]= 2.0*q[0];
dddq[1][1][3]=-2.0*q[3];
dddq[2][1][3]= 2.0*q[2];
dddq[0][2][3]= 2.0*q[1];
dddq[1][2][3]= 2.0*q[2];
dddq[2][2][3]= 2.0*q[3];

#ifdef EXTREME_DEBUG
printf("\n");
for(i=0;i<4;i++){
	for(j=0;j<3;j++){
		printf("MATR %12.6f %12.6f %12.6f\n",dddq[j][0][i],dddq[j][1][i],dddq[j][2][i]);
	}
        printf("\n");
}
#endif
/*
 * Build gamma 3x3x3 matrix
 */
for(i=0;i<3;i++){     //direction 
    for(j=0;j<3;j++){     //direction 
        for(k=0;k<3;k++){     //eigenvector number
            gamma[i][j][k]=0.0;
            for(l=0;l<4;l++){   //components of each eigenvector in pert. series
              if(lambda[0]==lambda[k+1]){
                 printf("FOUND DEGENERACY IN RMSD_ESS ROUTINE \n");
               /*  write(*,*)"FOUND DEGENERACY IN RMSD_ESS ROUTINE "
                 write(*,*)"I'm DYING...."
                 write(*,*)"COPYING STACK HERE "
                 write(*,*)"R0"
                 do ll=1,n
                  write(*,'(f8.3,f8.3,f8.3)')r0(1,ll),r0(2,ll),r0(3,ll)
                 enddo
                 write(*,*)"R"
                 do ll=1,n
                  write(*,'(f8.3,f8.3,f8.3)')r(1,ll),r(2,ll),r(3,ll)
                 enddo
                 stop*/
		 EXIT();} 
              else{
                gamma[i][j][k]=gamma[i][j][k]+dddq[i][j][l]*m[l][k+1]/(lambda[0]-lambda[k+1]);
	      }
	    }
	}

    }	
}
#ifdef EXTREME_DEBUG
for(i=0;i<3;i++){
	for(j=0;j<3;j++){
		printf("GAMM %12.6f %12.6f %12.6f\n",gamma[j][0][i],gamma[j][1][i],gamma[j][2][i]);
	}
        printf("\n");
}
#endif
/* 
 * Table of Derivative of the quaternion matrix respect to atom position
 */
for(i=0;i<inpack.natoms;i++){

        tmp1=(inpack.mass[i]); 
        rr1[0]=2.*outpack->r1p[0][i]*tmp1;
        rr1[1]=2.*outpack->r1p[1][i]*tmp1;
        rr1[2]=2.*outpack->r1p[2][i]*tmp1;
        rr0[0]=2.*outpack->r0p[0][i]*tmp1;
        rr0[1]=2.*outpack->r0p[1][i]*tmp1;
        rr0[2]=2.*outpack->r0p[2][i]*tmp1;
     

#ifdef EXTREME_DEBUG
        printf("ATOM %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n",rr0[0],rr0[1],rr0[2],rr1[0],rr1[1],rr1[2]);
#endif

        dm_r1 [0][0][0]=(rr1[0]-rr0[0]);
        dm_r1 [0][0][1]=(rr1[1]-rr0[1]);
        dm_r1 [0][0][2]=(rr1[2]-rr0[2]);
                      
        dm_r1 [0][1][0]=0.;
        dm_r1 [0][1][1]= rr0[2];
        dm_r1 [0][1][2]=-rr0[1];
                      
        dm_r1 [0][2][0]=-rr0[2];
        dm_r1 [0][2][1]= 0.;
        dm_r1 [0][2][2]= rr0[0];
                      
        dm_r1 [0][3][0]= rr0[1];
        dm_r1 [0][3][1]=-rr0[0];
        dm_r1 [0][3][2]= 0.;
                      
        dm_r1 [1][1][0]=(rr1[0]-rr0[0]);
        dm_r1 [1][1][1]=(rr1[1]+rr0[1]);
        dm_r1 [1][1][2]=(rr1[2]+rr0[2]);
                      
        dm_r1 [1][2][0]=-rr0[1];
        dm_r1 [1][2][1]=-rr0[0];
        dm_r1 [1][2][2]= 0.;
                      
        dm_r1 [1][3][0]=-rr0[2];
        dm_r1 [1][3][1]= 0.;
        dm_r1 [1][3][2]=-rr0[0];
                      
        dm_r1 [2][2][0]=(rr1[0]+rr0[0]);
        dm_r1 [2][2][1]=(rr1[1]-rr0[1]);
        dm_r1 [2][2][2]=(rr1[2]+rr0[2]);
                      
        dm_r1 [2][3][0]=0.;
        dm_r1 [2][3][1]=-rr0[2];
        dm_r1 [2][3][2]=-rr0[1];
                      
        dm_r1 [3][3][0]=(rr1[0]+rr0[0]);
        dm_r1 [3][3][1]=(rr1[1]+rr0[1]);
        dm_r1 [3][3][2]=(rr1[2]-rr0[2]);
/*
  derivative respec to to the other vector
 */
        dm_r0 [0][0][0]=-(rr1[0]-rr0[0]);
        dm_r0 [0][0][1]=-(rr1[1]-rr0[1]);
        dm_r0 [0][0][2]=-(rr1[2]-rr0[2]);
                      
        dm_r0 [0][1][0]=0.       ;
        dm_r0 [0][1][1]=-rr1[2];
        dm_r0 [0][1][2]=rr1[1];
                      
        dm_r0 [0][2][0]= rr1[2];      
        dm_r0 [0][2][1]= 0.;
        dm_r0 [0][2][2]=-rr1[0];
                      
        dm_r0 [0][3][0]=-rr1[1] ;     
        dm_r0 [0][3][1]= rr1[0];
        dm_r0 [0][3][2]= 0.;
                      
        dm_r0 [1][1][0]=-(rr1[0]-rr0[0]);
        dm_r0 [1][1][1]=(rr1[1]+rr0[1]);
        dm_r0 [1][1][2]=(rr1[2]+rr0[2]);
                      
        dm_r0 [1][2][0]=-rr1[1];
        dm_r0 [1][2][1]=-rr1[0];
        dm_r0 [1][2][2]= 0.;
                      
        dm_r0 [1][3][0]=-rr1[2];
        dm_r0 [1][3][1]= 0.;
        dm_r0 [1][3][2]=-rr1[0];
                      
        dm_r0 [2][2][0]=(rr1[0]+rr0[0]);
        dm_r0 [2][2][1]=-(rr1[1]-rr0[1]);
        dm_r0 [2][2][2]=(rr1[2]+rr0[2]);
                      
        dm_r0 [2][3][0]=0.;
        dm_r0 [2][3][1]=-rr1[2];
        dm_r0 [2][3][2]=-rr1[1];
                      
        dm_r0 [3][3][0]=(rr1[0]+rr0[0]);
        dm_r0 [3][3][1]=(rr1[1]+rr0[1]);
        dm_r0 [3][3][2]=-(rr1[2]-rr0[2]);
/*
 * write the diagonal
 */ 
	for(j=0;j<3;j++){

          dm_r1[1][0][j]=dm_r1[0][1][j];
          dm_r1[2][0][j]=dm_r1[0][2][j];
          dm_r1[3][0][j]=dm_r1[0][3][j];
          dm_r1[2][1][j]=dm_r1[1][2][j];
          dm_r1[3][1][j]=dm_r1[1][3][j];
          dm_r1[3][2][j]=dm_r1[2][3][j];

          dm_r0[1][0][j]=dm_r0[0][1][j];
          dm_r0[2][0][j]=dm_r0[0][2][j];
          dm_r0[3][0][j]=dm_r0[0][3][j];
          dm_r0[2][1][j]=dm_r0[1][2][j];
          dm_r0[3][1][j]=dm_r0[1][3][j];
          dm_r0[3][2][j]=dm_r0[2][3][j];
	  
          for(ll=0;ll<4;ll++){
          	for(mm=0;mm<4;mm++){
          		dm_r0_store[ll][mm][j][i]=dm_r0[ll][mm][j];
          		dm_r1_store[ll][mm][j][i]=dm_r1[ll][mm][j];
		};
	  };
 
	}
#ifdef EXTREME_DEBUG
	for(k=0;k<4;k++){
	for(l=0;l<4;l++){
	 printf("DM_R0 %12.6f %12.6f %12.6f\n",dm_r0[k][l][0],dm_r0[k][l][1],dm_r0[k][l][2]);
	}
        printf("\n"); 
        };
        for(k=0;k<4;k++){
	for(l=0;l<4;l++){
          printf("DM_R1 %12.6f %12.6f %12.6f\n",dm_r1[k][l][0],dm_r1[k][l][1],dm_r1[k][l][2]);
	}
        printf("\n"); 
        };
#endif
/*
 * pi matrix : coefficents in per theory
 */
	for(j=0;j<3;j++){
          pi1[0][j]=0.;
          pi1[1][j]=0.;
          pi1[2][j]=0.;
          pi0[0][j]=0.;
          pi0[1][j]=0.;
          pi0[2][j]=0.;
          outpack->derr_dr1 [j][i]=0.;
          outpack->derr_dr0 [j][i]=0.;

          for(k=0;k<4;k++){
            for(l=0;l<4;l++){
              outpack->derr_dr1[j][i]=outpack->derr_dr1[j][i]+q[k]*q[l]*dm_r1[l][k][j];
              outpack->derr_dr0[j][i]=outpack->derr_dr0[j][i]+q[k]*q[l]*dm_r0[l][k][j];
              for(mm=0;mm<3;mm++){
                pi0[mm][j]+=m[k][mm+1]*dm_r0[l][k][j]*q[l];
                pi1[mm][j]+=m[k][mm+1]*dm_r1[l][k][j]*q[l];  
	      };
	    };
	  };
          outpack->derr_dr1[j][i]=outpack->derr_dr1[j][i]/sqrt(4.*inpack.natoms*lambda[0]);
          outpack->derr_dr0[j][i]=outpack->derr_dr0[j][i]/sqrt(4.*inpack.natoms*lambda[0]);
	};
	for(j=0;j<3;j++){
		for (k=0;k<3;k++){
			for(l=0;l<3;l++){	    
              		outpack->dd_dr1[j][k][l][i]=0.;
              		outpack->dd_dr0[j][k][l][i]=0.;
			for(ii=0;ii<3;ii++){
                  		outpack->dd_dr1[j][k][l][i]+=gamma[j][k][ii]*pi1[ii][l]; 
                		outpack->dd_dr0[j][k][l][i]+=gamma[j][k][ii]*pi0[ii][l]; 
				}
			}
		}
	}
}
/*
 * Check arrays
 */
#ifdef EXTREME_DEBUG
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR0 %12.6f %12.6f %12.6f\n",outpack->derr_dr0[0][i],outpack->derr_dr0[1][i],outpack->derr_dr0[2][i]);
}
printf("\n");
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR1 %12.6f %12.6f %12.6f\n",outpack->derr_dr1[0][i],outpack->derr_dr1[1][i],outpack->derr_dr1[2][i]);
}
for(i=0;i<inpack.natoms;i++){
for(j=0;j<3;j++){
for(k=0;k<3;k++){
printf("DD_DR0 %12.6f %12.6f %12.6f\n",outpack->dd_dr0[j][k][0][i],outpack->dd_dr0[j][k][1][i],outpack->dd_dr0[j][k][2][i]);
}}}
for(i=0;i<inpack.natoms;i++){
for(j=0;j<3;j++){
for(k=0;k<3;k++){
printf("DD_DR1 %12.6f %12.6f %12.6f\n",outpack->dd_dr1[j][k][0][i],outpack->dd_dr1[j][k][1][i],outpack->dd_dr1[j][k][2][i]);
}}}
#endif
/*
 * Second derivative if asked
 *
 *
 */
if(iopt2==2){
/*
 *   dr0 dr
 */

for(i=0;i<3;i++){      //  r0 atom component
	for(k=0;k<inpack.natoms;k++){ // r0 atom index 
		for (j=0;j<3;j++){//  r atom component
        		for(l=0;l<inpack.natoms;l++){// r atom index
            			outpack->dderr_dr0_dr0[i][k][j][l]=0.;
            			outpack->dderr_dr0_dr1[i][k][j][l]=0.;
            			outpack->dderr_dr1_dr0[i][k][j][l]=0.;
            			outpack->dderr_dr1_dr1[i][k][j][l]=0.;
            			for(p=1;p<4;p++){//eigenvector index 
              				 tmp1=0.;
            				 tmp2=0.;
              				 tmp3=0.;
              				 tmp4=0.;
  		 	       		 for(nn=0;nn<4;nn++){//eigenvector component
  		 	       		 	for(mm=0;mm<4;mm++){//eigenvector component
                    					tmp1+=(m[nn][0]*dm_r0_store[nn][mm][i][k]*m[mm][p]);
                    					tmp2+=(m[nn][0]*dm_r1_store[nn][mm][j][l]*m[mm][p]);
                					tmp3+=(m[nn][0]*dm_r1_store[nn][mm][i][k]*m[mm][p]);
                    					tmp4+=(m[nn][0]*dm_r0_store[nn][mm][j][l]*m[mm][p]);
						}
					 }
              				 outpack->dderr_dr0_dr1[i][k][j][l]+=2.*tmp1*tmp2/(lambda[0]-lambda[p]);
              				 outpack->dderr_dr1_dr0[i][k][j][l]+=2.*tmp3*tmp4/(lambda[0]-lambda[p]);
               				 outpack->dderr_dr1_dr1[i][k][j][l]+=2.*tmp2*tmp3/(lambda[0]-lambda[p]);
               				 outpack->dderr_dr0_dr0[i][k][j][l]+=2.*tmp1*tmp4/(lambda[0]-lambda[p]);
				};

/*
 *  second order diagonal and semi-diagonal terms
 */
				if(k-l==0){ 
					if(i-j==0){ 
				
						outpack->dderr_dr1_dr1[i][k][j][l]+=2.;
              					outpack->dderr_dr0_dr0[i][k][j][l]+=2.;
						if(i==0){
               						tmp5=2.*(-pow(m[0][0],2) - pow(m[1][0],2) 
									+pow(m[2][0],2) +pow(m[3][0],2));
               						tmp6=tmp5;
						}
						if(i==1){
							tmp5=2.0*(-pow(m[0][0],2)+ pow(m[1][0],2) 
									-pow(m[2][0],2) +pow(m[3][0],2));
							tmp6=tmp5;
						}
						
              					if(i==2){
						        tmp5=2.*(-pow(m[0][0],2) + pow(m[1][0],2) 
									+ pow(m[2][0],2) -pow(m[3][0],2));
               						tmp6=tmp5;
						}

					}
					else{

             					if( i==1 && j==0 ){// dy dx 
                					tmp5=4.*(-m[0][0]*m[3][0]-m[1][0]*m[2][0]);
                					tmp6=4.*( m[0][0]*m[3][0]-m[1][0]*m[2][0]);
              					};
             					if( i==2 && j==0 ){// dz dx 
                                                        tmp5=4.*( m[0][0]*m[2][0]-m[1][0]*m[3][0]);
              						tmp6=4.*(-m[0][0]*m[2][0]-m[1][0]*m[3][0]);
						};
             					if( i==0 && j==1 ){// dx dy 
                					tmp5=4.*( m[0][0]*m[3][0]-m[1][0]*m[2][0]);
                					tmp6=4.*(-m[0][0]*m[3][0]-m[1][0]*m[2][0]);
						};
             					if( i==2 && j==1 ){// dz dx 
                					tmp5=4.*(-m[0][0]*m[1][0]-m[2][0]*m[3][0]);
                					tmp6=4.*( m[0][0]*m[1][0]-m[2][0]*m[3][0]);
						};
             					if( i==0 && j==2 ){// dx dz 
                					tmp5=4.*(-m[2][0]*m[0][0]-m[3][0]*m[1][0]);
                					tmp6=4.*( m[2][0]*m[0][0]-m[3][0]*m[1][0]);
              					};
             					if( i==1 && j==2 ){// dy dz 
                					tmp5=4.*( m[1][0]*m[0][0]-m[3][0]*m[2][0]);
                					tmp6=4.*(-m[1][0]*m[0][0]-m[3][0]*m[2][0]);
              					};

					};
            		 		outpack->dderr_dr1_dr0[i][k][j][l]+=tmp5;
             				outpack->dderr_dr0_dr1[i][k][j][l]+=tmp6;
				}; 
            			outpack->dderr_dr0_dr1[i][k][j][l]=
                                outpack->dderr_dr0_dr1[i][k][j][l]/(2.*sqrt(lambda[0]*((real) inpack.natoms)))
                                 -sqrt(((real) inpack.natoms)/lambda[0])*outpack->derr_dr0[i][k]*outpack->derr_dr1[j][l];

    			        outpack->dderr_dr1_dr0[i][k][j][l]=
          			 outpack->dderr_dr1_dr0[i][k][j][l]/(2.*sqrt(lambda[0]*((real) inpack.natoms)))
            			-sqrt(((real) inpack.natoms)/lambda[0])*outpack->derr_dr1[i][k]*outpack->derr_dr0[j][l];

        		        outpack->dderr_dr1_dr1[i][k][j][l]=
           			 outpack->dderr_dr1_dr1[i][k][j][l]/(2.*sqrt(lambda[0]*((real) inpack.natoms)))
            			-sqrt(((real) inpack.natoms)/lambda[0])*outpack->derr_dr1[i][k]*outpack->derr_dr1[j][l];

            			outpack->dderr_dr0_dr0[i][k][j][l]=
          			 outpack->dderr_dr0_dr0[i][k][j][l]/(2.*sqrt(lambda[0]*((real) inpack.natoms)))
           			 -sqrt(((real) inpack.natoms)/lambda[0])*outpack->derr_dr0[i][k]*outpack->derr_dr0[j][l];
			};
		};
	};
};
};
#ifdef EXTREME_DEBUG
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR0_DR0 %12.6f %12.6f %12.6f\n",outpack->dderr_dr0_dr0[j][i][0][l],outpack->dderr_dr0_dr0[j][i][1][l],outpack->dderr_dr0_dr0[j][i][2][l]);
		}
	}
	
}
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR1_DR1 %12.6f %12.6f %12.6f\n",outpack->dderr_dr1_dr1[j][i][0][l],outpack->dderr_dr1_dr1[j][i][1][l],outpack->dderr_dr1_dr1[j][i][2][l]);
		}
	}
	
}	
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR0_DR1 %12.6f %12.6f %12.6f\n",outpack->dderr_dr0_dr1[j][i][0][l],outpack->dderr_dr0_dr1[j][i][1][l],outpack->dderr_dr0_dr1[j][i][2][l]);
		}
	}
	
}	
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR1_DR0 %12.6f %12.6f %12.6f\n",outpack->dderr_dr1_dr0[j][i][0][l],outpack->dderr_dr1_dr0[j][i][1][l],outpack->dderr_dr1_dr0[j][i][2][l]);
		}
	}
	
}
#endif
/*
 * Now correct for cm - hard part in 2nd derivative
 *
 */
if(iopt==6 || iopt==7){
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
                derr_dr1_tmp[l][i]=outpack->derr_dr1[l][i];
			for(j=0;j<inpack.natoms;j++){
                        	derr_dr1_tmp[l][i]-=(1./((real) inpack.natoms))*outpack->derr_dr1[l][j];
			}
		}	
	}
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
			outpack->derr_dr1[l][i]=derr_dr1_tmp[l][i];
		}
	}	
}
if(iopt==5 || iopt==7){
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
                derr_dr0_tmp[l][i]=outpack->derr_dr0[l][i];
			for(j=0;j<inpack.natoms;j++){
                        	derr_dr0_tmp[l][i]-=(1./((real) inpack.natoms))*outpack->derr_dr0[l][j];
			}
		}	
	}
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
			outpack->derr_dr0[l][i]=derr_dr0_tmp[l][i];
		}
	}	
}
#ifdef EXTREME_DEBUG
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR0 %12.6f %12.6f %12.6f\n",outpack->derr_dr0[0][i],outpack->derr_dr0[1][i],outpack->derr_dr0[2][i]);
}
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR1_CM %12.6f %12.6f %12.6f\n",outpack->derr_dr1[0][i],outpack->derr_dr1[1][i],outpack->derr_dr1[2][i]);
}
#endif
if(iopt2==2){
	if(iopt==6){//dr1 correction
		for(ix=0;ix<3;ix++){	
			for(jx=0;jx<3;jx++){	
                                alpha_m1[ix][jx]=0.0;
				for(i=0;i<inpack.natoms;i++){	
					for(j=0;j<inpack.natoms;j++){	
                                        alpha_m1[ix][jx]+=outpack->dderr_dr1_dr1[ix][i][jx][j];
					}
				}
		        	alpha_m1[ix][jx]= alpha_m1[ix][jx]/(((real) inpack.natoms*inpack.natoms));
			}
		}
		for(ix=0;ix<3;ix++){	
			for(i=0;i<inpack.natoms;i++){	
				for(jx=0;jx<3;jx++){	
					for(j=0;j<inpack.natoms;j++){	
						dderr_dr1_dr1_tmp[ix][i][jx][j]=outpack->dderr_dr1_dr1[ix][i][jx][j];
						dderr_dr1_dr0_tmp[ix][i][jx][j]=outpack->dderr_dr0_dr1[ix][i][jx][j];
						dderr_dr0_dr1_tmp[ix][i][jx][j]=outpack->dderr_dr1_dr0[ix][i][jx][j];
						for(mm=0;mm<inpack.natoms;mm++){
							dderr_dr1_dr0_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*outpack->dderr_dr1_dr0[ix][i][jx][mm];
							dderr_dr0_dr1_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*outpack->dderr_dr0_dr1[ix][i][jx][mm];
							dderr_dr1_dr1_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*
						         	(outpack->dderr_dr1_dr1[ix][i][jx][mm]+outpack->dderr_dr1_dr1[ix][mm][jx][j]);
						}
						dderr_dr1_dr1_tmp[ix][i][jx][j]+=alpha_m1[ix][jx];

					}
				}
			}

		}
		for(ix=0;ix<3;ix++){	
			for(i=0;i<inpack.natoms;i++){	
				for(jx=0;jx<3;jx++){	
					for(j=0;j<inpack.natoms;j++){	
						outpack->dderr_dr1_dr1[ix][i][jx][j]=dderr_dr1_dr1_tmp[ix][i][jx][j];
						outpack->dderr_dr1_dr0[ix][i][jx][j]=dderr_dr1_dr0_tmp[ix][i][jx][j];
						outpack->dderr_dr0_dr1[ix][i][jx][j]=dderr_dr0_dr1_tmp[ix][i][jx][j];
					}
				}
			}
		}	
	}
	else if(iopt==5){//dr0 correction
		for(ix=0;ix<3;ix++){	
			for(jx=0;jx<3;jx++){	
                                alpha_m1[ix][jx]=0.0;
				for(i=0;i<inpack.natoms;i++){	
					for(j=0;j<inpack.natoms;j++){	
                                        alpha_m1[ix][jx]+=outpack->dderr_dr0_dr0[ix][i][jx][j];
					}
				}
		        	alpha_m1[ix][jx]= alpha_m1[ix][jx]/(((real) inpack.natoms*inpack.natoms));
			}
		}
		for(ix=0;ix<3;ix++){	
			for(i=0;i<inpack.natoms;i++){	
				for(jx=0;jx<3;jx++){	
					for(j=0;j<inpack.natoms;j++){	
						dderr_dr0_dr0_tmp[ix][i][jx][j]=outpack->dderr_dr0_dr0[ix][i][jx][j];
						dderr_dr1_dr0_tmp[ix][i][jx][j]=outpack->dderr_dr1_dr0[ix][i][jx][j];
						dderr_dr0_dr1_tmp[ix][i][jx][j]=outpack->dderr_dr0_dr1[ix][i][jx][j];
						for(mm=0;mm<inpack.natoms;mm++){
							dderr_dr1_dr0_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*outpack->dderr_dr1_dr0[ix][i][jx][mm];
							dderr_dr0_dr1_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*outpack->dderr_dr0_dr1[ix][i][jx][mm];
							dderr_dr0_dr0_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*
						         	(outpack->dderr_dr0_dr0[ix][i][jx][mm]+outpack->dderr_dr0_dr0[ix][mm][jx][j]);
						}
						dderr_dr0_dr0_tmp[ix][i][jx][j]+= alpha_m1[ix][jx];
					}
				}
			}
	
		}
		for(ix=0;ix<3;ix++){	
			for(i=0;i<inpack.natoms;i++){	
				for(jx=0;jx<3;jx++){	
					for(j=0;j<inpack.natoms;j++){	
						outpack->dderr_dr0_dr0[ix][i][jx][j]=dderr_dr0_dr0_tmp[ix][i][jx][j];
						outpack->dderr_dr1_dr0[ix][i][jx][j]=dderr_dr1_dr0_tmp[ix][i][jx][j];
						outpack->dderr_dr0_dr1[ix][i][jx][j]=dderr_dr0_dr1_tmp[ix][i][jx][j];
					}
				}
			}
		}	
	}
	else if(iopt==7){
		for(ix=0;ix<3;ix++){	
			for(jx=0;jx<3;jx++){	
                                alpha_m1[ix][jx]=0.0;
                                alpha_m2[ix][jx]=0.0;
                                alpha_m3[ix][jx]=0.0;
                                alpha_m4[ix][jx]=0.0;
				for(i=0;i<inpack.natoms;i++){	
					for(j=0;j<inpack.natoms;j++){	
                                     		alpha_m1[ix][jx]+=outpack->dderr_dr0_dr0[ix][i][jx][j];
                           	        	alpha_m2[ix][jx]+=outpack->dderr_dr1_dr1[ix][i][jx][j];
                               			alpha_m3[ix][jx]+=outpack->dderr_dr0_dr1[ix][i][jx][j];
                               	       	 	alpha_m4[ix][jx]+=outpack->dderr_dr1_dr0[ix][i][jx][j];
					}
				}
		        	alpha_m1[ix][jx]= alpha_m1[ix][jx]/(((real) inpack.natoms*inpack.natoms));
		        	alpha_m2[ix][jx]= alpha_m2[ix][jx]/(((real) inpack.natoms*inpack.natoms));
		        	alpha_m3[ix][jx]= alpha_m3[ix][jx]/(((real) inpack.natoms*inpack.natoms));
		        	alpha_m4[ix][jx]= alpha_m4[ix][jx]/(((real) inpack.natoms*inpack.natoms));
			}
		}
		for(ix=0;ix<3;ix++){	
			for(i=0;i<inpack.natoms;i++){	
				for(jx=0;jx<3;jx++){	
					for(j=0;j<inpack.natoms;j++){	
						dderr_dr0_dr0_tmp[ix][i][jx][j]=outpack->dderr_dr0_dr0[ix][i][jx][j];
						dderr_dr1_dr1_tmp[ix][i][jx][j]=outpack->dderr_dr1_dr1[ix][i][jx][j];
						dderr_dr1_dr0_tmp[ix][i][jx][j]=outpack->dderr_dr1_dr0[ix][i][jx][j];
						dderr_dr0_dr1_tmp[ix][i][jx][j]=outpack->dderr_dr0_dr1[ix][i][jx][j];
						for(mm=0;mm<inpack.natoms;mm++){
							dderr_dr0_dr0_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*
						         	(outpack->dderr_dr0_dr0[ix][i][jx][mm]+outpack->dderr_dr0_dr0[ix][mm][jx][j]);
							dderr_dr1_dr1_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*
						         	(outpack->dderr_dr1_dr1[ix][i][jx][mm]+outpack->dderr_dr1_dr1[ix][mm][jx][j]);
							dderr_dr0_dr1_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*
						         	(outpack->dderr_dr0_dr1[ix][i][jx][mm]+outpack->dderr_dr0_dr1[ix][mm][jx][j]);
							dderr_dr1_dr0_tmp[ix][i][jx][j]-=
								(1./((real) inpack.natoms))*
					         		(outpack->dderr_dr1_dr0[ix][i][jx][mm]+outpack->dderr_dr1_dr0[ix][mm][jx][j]);
						}
						dderr_dr0_dr0_tmp[ix][i][jx][j]+=alpha_m1[ix][jx];
						dderr_dr1_dr1_tmp[ix][i][jx][j]+=alpha_m2[ix][jx];
						dderr_dr0_dr1_tmp[ix][i][jx][j]+=alpha_m3[ix][jx];
						dderr_dr1_dr0_tmp[ix][i][jx][j]+=alpha_m4[ix][jx];
					}
				}
			}

		}
		for(ix=0;ix<3;ix++){	
			for(i=0;i<inpack.natoms;i++){	
				for(jx=0;jx<3;jx++){	
					for(j=0;j<inpack.natoms;j++){	
						outpack->dderr_dr0_dr0[ix][i][jx][j]=dderr_dr0_dr0_tmp[ix][i][jx][j];
						outpack->dderr_dr1_dr1[ix][i][jx][j]=dderr_dr1_dr1_tmp[ix][i][jx][j];
						outpack->dderr_dr1_dr0[ix][i][jx][j]=dderr_dr1_dr0_tmp[ix][i][jx][j];
						outpack->dderr_dr0_dr1[ix][i][jx][j]=dderr_dr0_dr1_tmp[ix][i][jx][j];
					}
				}
			}
		}	

	}
}
#ifdef EXTREME_DEBUG
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR0_DR0_CM %12.6f %12.6f %12.6f\n",outpack->dderr_dr0_dr0[j][i][0][l],outpack->dderr_dr0_dr0[j][i][1][l],outpack->dderr_dr0_dr0[j][i][2][l]);
		}
	}
	
}
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR1_DR1_CM %12.6f %12.6f %12.6f\n",outpack->dderr_dr1_dr1[j][i][0][l],outpack->dderr_dr1_dr1[j][i][1][l],outpack->dderr_dr1_dr1[j][i][2][l]);
		}
	}
	
}	
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR0_DR1_CM %12.6f %12.6f %12.6f\n",outpack->dderr_dr0_dr1[j][i][0][l],outpack->dderr_dr0_dr1[j][i][1][l],outpack->dderr_dr0_dr1[j][i][2][l]);
		}
	}
	
}	
for(i=0;i<inpack.natoms;i++){
	for(j=0;j<3;j++){
		for(l=0;l<inpack.natoms;l++){
			printf("DDERR_DR1_DR0_CM %12.6f %12.6f %12.6f\n",outpack->dderr_dr1_dr0[j][i][0][l],outpack->dderr_dr1_dr0[j][i][1][l],outpack->dderr_dr1_dr0[j][i][2][l]);
		}
	}
	
}
#endif

return 0;
}

// ------------------------------------------------------------------------------------------------


int PREFIX rmsd_mini_pack(struct rmsd_inpack inpack,struct rmsd_mini_outpack *outpack,int iopt)
{
/* declarations */
//cR int i,j,k,l,p,ll,mm,nn,ii,ix,jx;
int i,j,k,l,ll,mm,ii;
real rrsq,xx,yy,zz,m[4][4],rr1[4],rr0[4];
//cR double lambda[4],z[4][4],wk[20],s,q[4];
real lambda[4],s,q[4];
real dddq[3][3][4],gamma[3][3][3];
real dm_r1[4][4][3],dm_r0[4][4][3];
real dm_r1_store[4][4][3][MAXATOMS_RMSD];
real dm_r0_store[4][4][3][MAXATOMS_RMSD];
real derr_dr1_tmp[3][MAXATOMS_RMSD];
real derr_dr0_tmp[3][MAXATOMS_RMSD];
real pi1[3][3],pi0[3][3],tmp1; 
//cR double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
//cR int ier,arg1,arg2,arg3;
//cR double alpha_m1[3][3],alpha_m2[3][3],alpha_m3[3][3],alpha_m4[3][3];
xx=0.;
yy=0.;
zz=0.;
tmp1=0.;
if(iopt==5 || iopt == 7 ){
	for(i=0;i<inpack.natoms;i++){
		xx+=inpack.r0[0][i]*inpack.mass[i];
		yy+=inpack.r0[1][i]*inpack.mass[i];
		zz+=inpack.r0[2][i]*inpack.mass[i];
                tmp1+=inpack.mass[i];
	}
	xx=xx/((real) tmp1);
	yy=yy/((real) tmp1);
	zz=zz/((real) tmp1);
};
outpack->cmr0[0]=xx;
outpack->cmr0[1]=yy;
outpack->cmr0[2]=zz;
for(i=0;i<inpack.natoms;i++){
	outpack->r0p[0][i]=inpack.r0[0][i]-xx;
	outpack->r0p[1][i]=inpack.r0[1][i]-yy;
	outpack->r0p[2][i]=inpack.r0[2][i]-zz;
}
xx=0.;
yy=0.;
zz=0.;
tmp1=0.;
if(iopt==6 || iopt == 7 ){
	for(i=0;i<inpack.natoms;i++){
		xx+=inpack.r1[0][i]*inpack.mass[i];
		yy+=inpack.r1[1][i]*inpack.mass[i];
		zz+=inpack.r1[2][i]*inpack.mass[i];
                tmp1+=inpack.mass[i]; 
	};
	xx=xx/((real) tmp1);
	yy=yy/((real) tmp1);
	zz=zz/((real) tmp1);
};
outpack->cmr1[0]=xx;
outpack->cmr1[1]=yy;
outpack->cmr1[2]=zz;
for(i=0;i<inpack.natoms;i++){
	outpack->r1p[0][i]=inpack.r1[0][i]-xx;
	outpack->r1p[1][i]=inpack.r1[1][i]-yy;
	outpack->r1p[2][i]=inpack.r1[2][i]-zz;
}
// CLEAN M MATRIX
for(i=0;i<4;i++){
	for(j=0;j<4;j++){
          m[i][j]=0.;  
	}
}
// ASSIGN MATRIX ELEMENTS
for(i=0;i<inpack.natoms;i++){
	
        tmp1=sqrt(inpack.mass[i]); 
        rr1[0]=outpack->r1p[0][i]*tmp1;
        rr1[1]=outpack->r1p[1][i]*tmp1;
        rr1[2]=outpack->r1p[2][i]*tmp1;
        rr0[0]=outpack->r0p[0][i]*tmp1;
        rr0[1]=outpack->r0p[1][i]*tmp1;
        rr0[2]=outpack->r0p[2][i]*tmp1;
	
        rrsq=pow(rr0[0],2)+pow(rr0[1],2)+pow(rr0[2],2)+pow(rr1[0],2)+pow(rr1[1],2)+pow(rr1[2],2);
     
        m[0][0] +=  rrsq+2.*(-rr0[0]*rr1[0]-rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m[1][1] +=  rrsq+2.*(-rr0[0]*rr1[0]+rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m[2][2] +=  rrsq+2.*(+rr0[0]*rr1[0]-rr0[1]*rr1[1]+rr0[2]*rr1[2]);
        m[3][3] +=  rrsq+2.*(+rr0[0]*rr1[0]+rr0[1]*rr1[1]-rr0[2]*rr1[2]);
        m[0][1] += 2.*(-rr0[1]*rr1[2]+rr0[2]*rr1[1]);
        m[0][2] += 2.*( rr0[0]*rr1[2]-rr0[2]*rr1[0]);
        m[0][3] += 2.*(-rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m[1][2] -= 2.*( rr0[0]*rr1[1]+rr0[1]*rr1[0]);
        m[1][3] -= 2.*( rr0[0]*rr1[2]+rr0[2]*rr1[0]);
        m[2][3] -= 2.*( rr0[1]*rr1[2]+rr0[2]*rr1[1]);

};
m[1][0] = m[0][1];
m[2][0] = m[0][2];
m[2][1] = m[1][2];
m[3][0] = m[0][3];
m[3][1] = m[1][3];
m[3][2] = m[2][3];

// DIAGONALIZE 

ql77_driver(m,lambda);
s=1.0;
if(m[0][0]<0.)s=-1.;//correct for negative values (?)
q[0]=s*m[0][0];
q[1]=s*m[1][0];
q[2]=s*m[2][0];
q[3]=s*m[3][0];
outpack->err=sqrt(lambda[0]/((real) inpack.natoms));
if(lambda[0]==lambda[1]){
printf("DIAGONALIZATION: NON UNIQUE SOLUTION: ABORT \n");
EXIT();
}
if(iopt==0){return 0;}// JUST DIAGONALIZATION REQUIRED 

/*
 * Find the ROTATION matrix
 */
outpack->d[0][0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3]       ; 
outpack->d[1][0]=2.0*(q[1]*q[2]-q[0]*q[3]);
outpack->d[2][0]=2.0*(q[1]*q[3]+q[0]*q[2]);
outpack->d[0][1]=2.0*(q[1]*q[2]+q[0]*q[3]);
outpack->d[1][1]=q[0]*q[0]+q[2]*q[2]-q[1]*q[1]-q[3]*q[3];
outpack->d[2][1]=2.0*(q[2]*q[3]-q[0]*q[1]);
outpack->d[0][2]=2.0*(q[1]*q[3]-q[0]*q[2]);
outpack->d[1][2]=2.0*(q[2]*q[3]+q[0]*q[1]);
outpack->d[2][2]=q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
#ifdef EXTREME_DEBUG
for (i=0;i<3;i++){
printf("D_MATRIX %12.6f %12.6f %12.6f\n",outpack->d[i][0],outpack->d[i][1],outpack->d[i][2]);
}
#endif
/* 
 * first derivative in perturbation theory
 */
dddq[0][0][0]= 2.0*q[0];
dddq[1][0][0]=-2.0*q[3];
dddq[2][0][0]= 2.0*q[2];
dddq[0][1][0]= 2.0*q[3];
dddq[1][1][0]= 2.0*q[0];
dddq[2][1][0]=-2.0*q[1];
dddq[0][2][0]=-2.0*q[2];
dddq[1][2][0]= 2.0*q[1];
dddq[2][2][0]= 2.0*q[0];

dddq[0][0][1]= 2.0*q[1];
dddq[1][0][1]= 2.0*q[2];
dddq[2][0][1]= 2.0*q[3];
dddq[0][1][1]= 2.0*q[2];
dddq[1][1][1]=-2.0*q[1];
dddq[2][1][1]=-2.0*q[0];
dddq[0][2][1]= 2.0*q[3];
dddq[1][2][1]= 2.0*q[0];
dddq[2][2][1]=-2.0*q[1];

dddq[0][0][2]=-2.0*q[2];
dddq[1][0][2]= 2.0*q[1];
dddq[2][0][2]= 2.0*q[0];
dddq[0][1][2]= 2.0*q[1];
dddq[1][1][2]= 2.0*q[2];
dddq[2][1][2]= 2.0*q[3];
dddq[0][2][2]=-2.0*q[0];
dddq[1][2][2]= 2.0*q[3];
dddq[2][2][2]=-2.0*q[2];

dddq[0][0][3]=-2.0*q[3];
dddq[1][0][3]=-2.0*q[0];
dddq[2][0][3]= 2.0*q[1];
dddq[0][1][3]= 2.0*q[0];
dddq[1][1][3]=-2.0*q[3];
dddq[2][1][3]= 2.0*q[2];
dddq[0][2][3]= 2.0*q[1];
dddq[1][2][3]= 2.0*q[2];
dddq[2][2][3]= 2.0*q[3];

#ifdef EXTREME_DEBUG
printf("\n");
for(i=0;i<4;i++){
	for(j=0;j<3;j++){
		printf("MATR %12.6f %12.6f %12.6f\n",dddq[j][0][i],dddq[j][1][i],dddq[j][2][i]);
	}
        printf("\n");
}
#endif
/*
 * Build gamma 3x3x3 matrix
 */
for(i=0;i<3;i++){     //direction 
    for(j=0;j<3;j++){     //direction 
        for(k=0;k<3;k++){     //eigenvector number
            gamma[i][j][k]=0.0;
            for(l=0;l<4;l++){   //components of each eigenvector in pert. series
              if(lambda[0]==lambda[k+1]){
                 printf("FOUND DEGENERACY IN RMSD_ESS ROUTINE \n");
               /*  write(*,*)"FOUND DEGENERACY IN RMSD_ESS ROUTINE "
                 write(*,*)"I'm DYING...."
                 write(*,*)"COPYING STACK HERE "
                 write(*,*)"R0"
                 do ll=1,n
                  write(*,'(f8.3,f8.3,f8.3)')r0(1,ll),r0(2,ll),r0(3,ll)
                 enddo
                 write(*,*)"R"
                 do ll=1,n
                  write(*,'(f8.3,f8.3,f8.3)')r(1,ll),r(2,ll),r(3,ll)
                 enddo
                 stop*/
		 EXIT();} 
              else{
                gamma[i][j][k]=gamma[i][j][k]+dddq[i][j][l]*m[l][k+1]/(lambda[0]-lambda[k+1]);
	      }
	    }
	}

    }	
}
#ifdef EXTREME_DEBUG
for(i=0;i<3;i++){
	for(j=0;j<3;j++){
		printf("GAMM %12.6f %12.6f %12.6f\n",gamma[j][0][i],gamma[j][1][i],gamma[j][2][i]);
	}
        printf("\n");
}
#endif
/* 
 * Table of Derivative of the quaternion matrix respect to atom position
 */
for(i=0;i<inpack.natoms;i++){

        tmp1=(inpack.mass[i]); 
        rr1[0]=2.*outpack->r1p[0][i]*tmp1;
        rr1[1]=2.*outpack->r1p[1][i]*tmp1;
        rr1[2]=2.*outpack->r1p[2][i]*tmp1;
        rr0[0]=2.*outpack->r0p[0][i]*tmp1;
        rr0[1]=2.*outpack->r0p[1][i]*tmp1;
        rr0[2]=2.*outpack->r0p[2][i]*tmp1;
     

#ifdef EXTREME_DEBUG
        printf("ATOM %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n",rr0[0],rr0[1],rr0[2],rr1[0],rr1[1],rr1[2]);
#endif

        dm_r1 [0][0][0]=(rr1[0]-rr0[0]);
        dm_r1 [0][0][1]=(rr1[1]-rr0[1]);
        dm_r1 [0][0][2]=(rr1[2]-rr0[2]);
                      
        dm_r1 [0][1][0]=0.;
        dm_r1 [0][1][1]= rr0[2];
        dm_r1 [0][1][2]=-rr0[1];
                      
        dm_r1 [0][2][0]=-rr0[2];
        dm_r1 [0][2][1]= 0.;
        dm_r1 [0][2][2]= rr0[0];
                      
        dm_r1 [0][3][0]= rr0[1];
        dm_r1 [0][3][1]=-rr0[0];
        dm_r1 [0][3][2]= 0.;
                      
        dm_r1 [1][1][0]=(rr1[0]-rr0[0]);
        dm_r1 [1][1][1]=(rr1[1]+rr0[1]);
        dm_r1 [1][1][2]=(rr1[2]+rr0[2]);
                      
        dm_r1 [1][2][0]=-rr0[1];
        dm_r1 [1][2][1]=-rr0[0];
        dm_r1 [1][2][2]= 0.;
                      
        dm_r1 [1][3][0]=-rr0[2];
        dm_r1 [1][3][1]= 0.;
        dm_r1 [1][3][2]=-rr0[0];
                      
        dm_r1 [2][2][0]=(rr1[0]+rr0[0]);
        dm_r1 [2][2][1]=(rr1[1]-rr0[1]);
        dm_r1 [2][2][2]=(rr1[2]+rr0[2]);
                      
        dm_r1 [2][3][0]=0.;
        dm_r1 [2][3][1]=-rr0[2];
        dm_r1 [2][3][2]=-rr0[1];
                      
        dm_r1 [3][3][0]=(rr1[0]+rr0[0]);
        dm_r1 [3][3][1]=(rr1[1]+rr0[1]);
        dm_r1 [3][3][2]=(rr1[2]-rr0[2]);
/*
  derivative respec to to the other vector
 */
        dm_r0 [0][0][0]=-(rr1[0]-rr0[0]);
        dm_r0 [0][0][1]=-(rr1[1]-rr0[1]);
        dm_r0 [0][0][2]=-(rr1[2]-rr0[2]);
                      
        dm_r0 [0][1][0]=0.       ;
        dm_r0 [0][1][1]=-rr1[2];
        dm_r0 [0][1][2]=rr1[1];
                      
        dm_r0 [0][2][0]= rr1[2];      
        dm_r0 [0][2][1]= 0.;
        dm_r0 [0][2][2]=-rr1[0];
                      
        dm_r0 [0][3][0]=-rr1[1] ;     
        dm_r0 [0][3][1]= rr1[0];
        dm_r0 [0][3][2]= 0.;
                      
        dm_r0 [1][1][0]=-(rr1[0]-rr0[0]);
        dm_r0 [1][1][1]=(rr1[1]+rr0[1]);
        dm_r0 [1][1][2]=(rr1[2]+rr0[2]);
                      
        dm_r0 [1][2][0]=-rr1[1];
        dm_r0 [1][2][1]=-rr1[0];
        dm_r0 [1][2][2]= 0.;
                      
        dm_r0 [1][3][0]=-rr1[2];
        dm_r0 [1][3][1]= 0.;
        dm_r0 [1][3][2]=-rr1[0];
                      
        dm_r0 [2][2][0]=(rr1[0]+rr0[0]);
        dm_r0 [2][2][1]=-(rr1[1]-rr0[1]);
        dm_r0 [2][2][2]=(rr1[2]+rr0[2]);
                      
        dm_r0 [2][3][0]=0.;
        dm_r0 [2][3][1]=-rr1[2];
        dm_r0 [2][3][2]=-rr1[1];
                      
        dm_r0 [3][3][0]=(rr1[0]+rr0[0]);
        dm_r0 [3][3][1]=(rr1[1]+rr0[1]);
        dm_r0 [3][3][2]=-(rr1[2]-rr0[2]);
/*
 * write the diagonal
 */ 
	for(j=0;j<3;j++){

          dm_r1[1][0][j]=dm_r1[0][1][j];
          dm_r1[2][0][j]=dm_r1[0][2][j];
          dm_r1[3][0][j]=dm_r1[0][3][j];
          dm_r1[2][1][j]=dm_r1[1][2][j];
          dm_r1[3][1][j]=dm_r1[1][3][j];
          dm_r1[3][2][j]=dm_r1[2][3][j];

          dm_r0[1][0][j]=dm_r0[0][1][j];
          dm_r0[2][0][j]=dm_r0[0][2][j];
          dm_r0[3][0][j]=dm_r0[0][3][j];
          dm_r0[2][1][j]=dm_r0[1][2][j];
          dm_r0[3][1][j]=dm_r0[1][3][j];
          dm_r0[3][2][j]=dm_r0[2][3][j];
	  
          for(ll=0;ll<4;ll++){
          	for(mm=0;mm<4;mm++){
          		dm_r0_store[ll][mm][j][i]=dm_r0[ll][mm][j];
          		dm_r1_store[ll][mm][j][i]=dm_r1[ll][mm][j];
		};
	  };
 
	}
#ifdef EXTREME_DEBUG
	for(k=0;k<4;k++){
	for(l=0;l<4;l++){
	 printf("DM_R0 %12.6f %12.6f %12.6f\n",dm_r0[k][l][0],dm_r0[k][l][1],dm_r0[k][l][2]);
	}
        printf("\n"); 
        };
        for(k=0;k<4;k++){
	for(l=0;l<4;l++){
          printf("DM_R1 %12.6f %12.6f %12.6f\n",dm_r1[k][l][0],dm_r1[k][l][1],dm_r1[k][l][2]);
	}
        printf("\n"); 
        };
#endif
/*
 * pi matrix : coefficents in per theory
 */
	for(j=0;j<3;j++){
          pi1[0][j]=0.;
          pi1[1][j]=0.;
          pi1[2][j]=0.;
          pi0[0][j]=0.;
          pi0[1][j]=0.;
          pi0[2][j]=0.;
          outpack->derr_dr1 [j][i]=0.;
          outpack->derr_dr0 [j][i]=0.;

          for(k=0;k<4;k++){
            for(l=0;l<4;l++){
              outpack->derr_dr1[j][i]=outpack->derr_dr1[j][i]+q[k]*q[l]*dm_r1[l][k][j];
              outpack->derr_dr0[j][i]=outpack->derr_dr0[j][i]+q[k]*q[l]*dm_r0[l][k][j];
              for(mm=0;mm<3;mm++){
                pi0[mm][j]+=m[k][mm+1]*dm_r0[l][k][j]*q[l];
                pi1[mm][j]+=m[k][mm+1]*dm_r1[l][k][j]*q[l];  
	      };
	    };
	  };
          outpack->derr_dr1[j][i]=outpack->derr_dr1[j][i]/sqrt(4.*inpack.natoms*lambda[0]);
          outpack->derr_dr0[j][i]=outpack->derr_dr0[j][i]/sqrt(4.*inpack.natoms*lambda[0]);
	};
	for(j=0;j<3;j++){
		for (k=0;k<3;k++){
			for(l=0;l<3;l++){	    
              		outpack->dd_dr1[j][k][l][i]=0.;
              		outpack->dd_dr0[j][k][l][i]=0.;
			for(ii=0;ii<3;ii++){
                  		outpack->dd_dr1[j][k][l][i]+=gamma[j][k][ii]*pi1[ii][l]; 
                		outpack->dd_dr0[j][k][l][i]+=gamma[j][k][ii]*pi0[ii][l]; 
				}
			}
		}
	}
}
/*
 * Check arrays
 */
#ifdef EXTREME_DEBUG
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR0 %12.6f %12.6f %12.6f\n",outpack->derr_dr0[0][i],outpack->derr_dr0[1][i],outpack->derr_dr0[2][i]);
}
printf("\n");
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR1 %12.6f %12.6f %12.6f\n",outpack->derr_dr1[0][i],outpack->derr_dr1[1][i],outpack->derr_dr1[2][i]);
}
for(i=0;i<inpack.natoms;i++){
for(j=0;j<3;j++){
for(k=0;k<3;k++){
printf("DD_DR0 %12.6f %12.6f %12.6f\n",outpack->dd_dr0[j][k][0][i],outpack->dd_dr0[j][k][1][i],outpack->dd_dr0[j][k][2][i]);
}}}
for(i=0;i<inpack.natoms;i++){
for(j=0;j<3;j++){
for(k=0;k<3;k++){
printf("DD_DR1 %12.6f %12.6f %12.6f\n",outpack->dd_dr1[j][k][0][i],outpack->dd_dr1[j][k][1][i],outpack->dd_dr1[j][k][2][i]);
}}}
#endif
/*
 * Now correct for cm - hard part in 2nd derivative
 *
 */
if(iopt==6 || iopt==7){
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
                derr_dr1_tmp[l][i]=outpack->derr_dr1[l][i];
			for(j=0;j<inpack.natoms;j++){
                        	derr_dr1_tmp[l][i]-=(1./((real) inpack.natoms))*outpack->derr_dr1[l][j];
			}
		}	
	}
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
			outpack->derr_dr1[l][i]=derr_dr1_tmp[l][i];
		}
	}	
}
if(iopt==5 || iopt==7){
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
                derr_dr0_tmp[l][i]=outpack->derr_dr0[l][i];
			for(j=0;j<inpack.natoms;j++){
                        	derr_dr0_tmp[l][i]-=(1./((real) inpack.natoms))*outpack->derr_dr0[l][j];
			}
		}	
	}
	for(l=0;l<3;l++){
		for(i=0;i<inpack.natoms;i++){
			outpack->derr_dr0[l][i]=derr_dr0_tmp[l][i];
		}
	}	
}
#ifdef EXTREME_DEBUG
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR0 %12.6f %12.6f %12.6f\n",outpack->derr_dr0[0][i],outpack->derr_dr0[1][i],outpack->derr_dr0[2][i]);
}
for(i=0;i<inpack.natoms;i++){
printf("DERR_DR1_CM %12.6f %12.6f %12.6f\n",outpack->derr_dr1[0][i],outpack->derr_dr1[1][i],outpack->derr_dr1[2][i]);
}
#endif
return 0;
}

// ------------------------------------------------------------------------------------------------

void PREFIX ql77_driver(real m[][4],real* lambda){
int i,j;
real ll[16];
for(i=0;i<4;i++){
 for(j=0;j<4;j++){
	 ll[4*i+j]=m[i][j]; }};
#ifdef EXTREME_DEBUG
for(i=0;i<4;i++){
         printf("M_MATRIX %12.6f %12.6f %12.6f %12.6f\n",m[i][0],m[i][1],m[i][2],m[i][3]) ;  
}
#endif
ql77(4,ll,lambda);
#ifdef EXTREME_DEBUG
printf("EIGV %f %f %f %f\n",lambda[0],lambda[1],lambda[2],lambda[3]);
#endif
//back to square representation: columns have eigenvectors
for(j=0;j<4;j++){
	for(i=0;i<4;i++){
		 m[i][j]=ll[4*j+i]; }
	};
};
void PREFIX ql77 (int n,real *x,real *d)
{
  int i,j,k,l,ni;
  real *e,h,g,f,b,s,p,r,c,absp;
  real totwork,work;

  const real eps=7.e-14,tol=1.e-30;

  e=(real *)malloc(n*sizeof(real));

  totwork = 0;
  for(ni=1; ni<n; ni++)
    totwork += pow((real) n-ni,3);

  work=0;
  for(ni=1; (ni < n); ni++) {
    i=n-ni;
    l=i-1;
    h=0.0;                                                          
    g=x[i+n*(i-1)];
    if (l > 0) {
      for(k=0; (k<l); k++) 
	h=h+x[i+n*k]*x[i+n*k];
      s=h+g*g;
      if (s < tol)
	h=0.0;
      else if (h > 0.0) { 
	l=l+1;
	f=g;
	g=sqrt(s);
	g=-g;
	h=s-f*g;                                                           
	x[i+n*(i-1)]=f-g;
	f=0.0;
	for(j=0; (j<l); j++) {
	  x[j+n*i]=x[i+n*j]/h;
	  s=0.0;
	  for(k=0; (k<=j); k++)
	    s=s+x[j+n*k]*x[i+n*k];
	  for(k=j+1; (k<l); k++)
	    s=s+x[k+n*j]*x[i+n*k];
	  e[j]=s/h;
	  f=f+s*x[j+n*i];
	}
	f=f/(h+h);
	for(j=0; (j<l); j++) 
	  e[j]=e[j]-f*x[i+n*j];
	for(j=0; (j<l); j++) {
	  f=x[i+n*j];
	  s=e[j];
	  for(k=0; (k<=j); k++)
	    x[j+n*k]=x[j+n*k]-f*e[k]-x[i+n*k]*s;
	}
      }
    }
    d[i]=h;
    e[i-1]=g;

    work += pow((real) n-ni,3);
  }

  /*
   *  accumulation of transformation matrix and intermediate d vector
   */
  
  d[0]=x[0];
  x[0]=1.0;

  work=0;
  for(i=1; (i<n); i++) {
    if (d[i] > 0.0) {
      for(j=0; (j<i); j++) {
	s=0.0;
	for(k=0; (k<i); k++) 
	  s=s+x[i+n*k]*x[k+n*j];
	for(k=0; (k<i); k++)
	  x[k+n*j]=x[k+n*j]-s*x[k+n*i];
      }
    }
    d[i]=x[i+n*i];
    x[i+n*i]=1.0;
    for(j=0; (j<i); j++) {
      x[i+n*j]=0.0;
      x[j+n*i]=0.0;
    }
    work += pow((real) i,3);
  }

  /*
   *  ql iterates
   */

  b=0.0;
  f=0.0;
  e[n-1]=0.0;
  totwork += pow((real) n,3);
  work=0;
  for(l=0; (l<n); l++) {
    h=eps*(fabs(d[l])+fabs(e[l]));
    if (h > b) 
      b=h;                                                   
    for(j=l; (j<n); j++) {
      if(fabs(e[j]) <= b) 
	break;
    }
    if (j != l) { 
      do {
	g=d[l];
	p=(d[l+1]-g)*0.5/e[l];
	r=sqrt(p*p+1.0);
	if(p < 0.0)
	  p=p-r;
	else
	  p=p+r;
	d[l]=e[l]/p;
	h=g-d[l];                                                     
	for(i=l+1; (i<n); i++)
	  d[i]=d[i]-h;                                                       
	f=f+h;                                                             
	p=d[j];
	c=1.0;
	s=0.0;
	for(ni=l; (ni<j); ni++) {
	  i=l+j-1-ni;
	  g=c*e[i];
	  h=c*p;
	  if(fabs(p) >= fabs(e[i])) {
	    c=e[i]/p;
	    r=sqrt(c*c+1.0);
	    e[i+1]=s*p*r;
	    s=c/r;
	    c=1.0/r;
	  } else {
	    c=p/e[i];
	    r=sqrt(c*c+1.0);
	    e[i+1]=s*e[i]*r;
	    s=1.0/r;
	    c=c/r;
	  }
	  p=c*d[i]-s*g;
	  d[i+1]=h+s*(c*g+s*d[i]);
	  for(k=0; (k<n); k++) {
	    h=x[k+n*(i+1)];
	    x[k+n*(i+1)]=x[k+n*i]*s+h*c;
	    x[k+n*i]=x[k+n*i]*c-h*s;
	  }
	}
	e[l]=s*p;
	d[l]=c*p;
      } while (fabs(e[l]) > b); 
    }
    d[l]=d[l]+f;

    work += pow((real) n-l,3);
  }

  /*
   *  put eigenvalues and eigenvectors in 
   *  desired ascending order
   */
  
  for(i=0; (i<n-1); i++) {
    k    = i;
    p    = d[i];
    absp = fabs(d[i]);
    for(j=i+1; (j<n); j++) {
      if(fabs(d[j]) < absp) {
	k    = j;
	p    = d[j];
	absp = fabs(d[j]);
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for(j=0; (j<n); j++) {
	p        = x[j+n*i];
	x[j+n*i] = x[j+n*k];
	x[j+n*k] = p;
      }
    }
  }

	  free(e);
	}

//-----------------------------------------------------------------------------------------------------

void PREFIX dmsd_calculation(int i_c,struct coordinates_frameset *pframeset,struct cmap_inpack *c_inpack,
                             struct cmap_outpack *c_outpack,real dmsd_dr1[3][MAXATOMS_PATH]){

    int i,j,ix;
    rvec rij,rij0;
    real mod_rij,mod_rij0, drmsd, fact;

    drmsd = 0.;
    for(i=0;i<colvar.natoms[i_c];i++) for(ix=0;ix<3;ix++) {c_outpack->derr_dr0[ix][i] = dmsd_dr1[ix][i] = 0.;} 

  
    for(i=0;i<colvar.natoms[i_c]-1;i++) {
     for(j=i+1;j<colvar.natoms[i_c];j++) {
      minimal_image(c_inpack->r0[i], c_inpack->r0[j], &mod_rij, rij);
      minimal_image(pframeset->pos[i], pframeset->pos[j], &mod_rij0, rij0);
      drmsd += (mod_rij-mod_rij0)*(mod_rij-mod_rij0);
      for(ix=0;ix<3;ix++) {
        c_outpack->derr_dr0[ix][i] +=  rij[ix]*(mod_rij-mod_rij0)/mod_rij;
        c_outpack->derr_dr0[ix][j] += -rij[ix]*(mod_rij-mod_rij0)/mod_rij;
        dmsd_dr1[ix][i] += -rij0[ix]*(mod_rij-mod_rij0)/mod_rij0;
        dmsd_dr1[ix][j] +=  rij0[ix]*(mod_rij-mod_rij0)/mod_rij0;
      }
     }
    } 

   fact = 2./((real)colvar.natoms[i_c]*((real)colvar.natoms[i_c]-1.));

   c_outpack->err = drmsd * fact; 

   for(i=0;i<colvar.natoms[i_c];i++) {
    c_outpack->derr_dr0[0][i] *= 2.*fact;
    c_outpack->derr_dr0[1][i] *= 2.*fact;
    c_outpack->derr_dr0[2][i] *= 2.*fact; 
    dmsd_dr1[0][i] *= 2.*fact;
    dmsd_dr1[1][i] *= 2.*fact;
    dmsd_dr1[2][i] *= 2.*fact;
   }

}

int PREFIX read_sz_hybrid(struct sz_data *my_sz, FILE *fplog) {

        int l,i,j,k,jj,kk,found;
        char *str,ic[3],str2[100];
        l=0;
/*
  * allocate the pointers: fully dynamical
 */
        my_sz->hbd_frameset=(struct hybrid_frameset **)malloc((my_sz->number)*sizeof(struct hybrid_frameset *));
        for(i=0;i< my_sz->number;i++){
                my_sz->hbd_frameset[i]=(struct hybrid_frameset *)malloc(sizeof(struct hybrid_frameset)) ;
                my_sz->hbd_frameset[i]->hbd_nelem=my_sz->nhybrid;
                my_sz->hbd_frameset[i]->hbd_elem=(struct hybrid_elem **)malloc((my_sz->hbd_frameset[i]->hbd_nelem)*sizeof(struct hybrid_elem * ));
                my_sz->hbd_frameset[i]->hbd_totalatoms=0;   
                // nelem=number of element in the hybrid representation
                for (j=0;j<my_sz->hbd_frameset[i]->hbd_nelem;j++){
                        my_sz->hbd_frameset[i]->hbd_elem[j]=(struct hybrid_elem *)malloc(sizeof(struct hybrid_elem));
                        // type of the cv in the list
                        my_sz->hbd_frameset[i]->hbd_elem[j]->cvtype=my_sz->lcvhybrid[j]; 
                        // progressive index in the list of the cvs
                        my_sz->hbd_frameset[i]->hbd_elem[j]->cvindex=my_sz->lhybrid[j]; 
                } 
        }

/*
 *  build names
 */
        str=&(my_sz->names[0]);
        for (i=1;i<= my_sz->number ;i++){
                strcpy(str2,my_sz->names);
                if(i<10){
                 ic[0]='0'+i;
                 ic[1]='\0';}
                else if(i<100) {
                 ic[0]='0'+i/10 ;
                 ic[1]='0'+i%10 ;
                 ic[2]='\0';
                }
                else{
                  fprintf(fplog,"|--read_sz_input: TOO MANY FRAMES REQUIRED FOR NAME BUILDING!\n");
                  EXIT();
                }
                strcat(str2,ic);
                strcat(str2,".dat");
                fprintf(fplog,"|--%s\n",str2);
                //read_sz_coord(str2,my_sz->frameset[i-1],fplog);
                FILE *fp;
                // test if the file exists 
                fp=fopen(str2,"r");
                if (fp == NULL)
                {fprintf(fplog,"|--UNSUCCESSFULL OPENING FOR FILE %s \n",str2);EXIT(); }; 
                // loop over the needed cv: open the file with the right reader and allocate accordingly 
                for (j=0;j<my_sz->nhybrid;j++){
                      switch(my_sz->lcvhybrid[j]) {
                         // read the distance, get back with the file pointer
                         case 1:  printf("CALLING DIST READER\n");
                                 my_sz->hbd_frameset[i-1]->hbd_totalatoms+=hbd_read_dist(fp,fplog,my_sz->hbd_frameset[i-1]->hbd_elem[j]); break;
                         case 3:  printf("CALLING COORD READER\n");
                                 my_sz->hbd_frameset[i-1]->hbd_totalatoms+=hbd_read_coord(fp,fplog,my_sz->hbd_frameset[i-1]->hbd_elem[j]); break;
                         case 4:  printf("CALLING ANGLE READER\n");
                                 my_sz->hbd_frameset[i-1]->hbd_totalatoms+=hbd_read_angle(fp,fplog,my_sz->hbd_frameset[i-1]->hbd_elem[j]); break;
                         case 5:  printf("CALLING TORSION READER\n");
                                 my_sz->hbd_frameset[i-1]->hbd_totalatoms+=hbd_read_torsion(fp,fplog,my_sz->hbd_frameset[i-1]->hbd_elem[j]); break;
                         case 31:  printf("CALLING TARGET READER\n");
                                 my_sz->hbd_frameset[i-1]->hbd_totalatoms+=hbd_read_target(fp,fplog,my_sz->hbd_frameset[i-1]->hbd_elem[j], &mtd_data); break;

                         default: printf("READER NOT IMPLEMENTED\n"); EXIT(); break;
                      } 
                }
                printf("TOTAL NUMBER OF ATOMS IN FRAMESET %d IS %d \n",i-1,my_sz->hbd_frameset[i-1]->hbd_totalatoms);
                // final lines should be the MATRIX
                char string[400],matrix[10],end[10],*str;
                while(1){
                     readagain:
                     str=fgets(string,100,fp);  
                     sscanf(str,"%3s",end);if(strstr(end,"END")!=NULL){break;};
                     sscanf(str,"%6s",matrix);  
                     if(strstr(matrix,"MATRIX")!=NULL){
                            for(jj=0;jj<my_sz->nhybrid ;jj++){ 
                     		str=fgets(string,100,fp);  
                                //printf("MATRIX:   %s",string);
                                char *pp=strtok(string," \t\n");
                                my_sz->mathybrid[jj][0]=atof(pp);//printf("XXX %d %d %s\n",jj,0,pp);
                                for(kk=1;kk<my_sz->nhybrid;kk++){
                                    char *pp=strtok(NULL," \t\n");//printf("XXX %d %d  %s\n",jj,kk,pp);
                                    my_sz->mathybrid[jj][kk]=atof(pp);
                                }
                            }
     	             }
		}
                // print the matrix
                for(jj=0;jj<my_sz->nhybrid;jj++){
                   printf("MATRIX "); 
                   for(kk=0;kk<my_sz->nhybrid;kk++){
                     	printf(" %12.6f ",my_sz->mathybrid[jj][kk]); 
		   }  	
                   printf("\n");
                }
        }

///// 
///// now write the backtable from elem 
///// 
        for (i=0;i< my_sz->number ;i++){
                l=0;
                // a global backtable to keep track of all the index of all atoms
                snew(my_sz->hbd_frameset[i]->backtable,my_sz->hbd_frameset[i]->hbd_totalatoms);
                printf("CHECK ON TOTAL NUMBER OF ATOMS\n");
                for (j=0;j<my_sz->hbd_frameset[i]->hbd_nelem;j++){
                   struct  hybrid_elem *ee;
                   ee=(my_sz->hbd_frameset[i]->hbd_elem[j]) ; 
                   for (k=0;k<ee->nder;k++){
                     my_sz->hbd_frameset[i]->backtable[l]= ee->backtable[k]; 
                     l++;
                   }
	 	}	
                printf("TOTAL NUMBER OF ATOMS IN FRAMESET %d IS %d \n",i,l);
                snew(my_sz->hbd_frameset[i]->myder,my_sz->hbd_frameset[i]->hbd_totalatoms); 
	}		
        // allocate the running structure: needed to calculate the reference before comparing to the stored value 
        my_sz->hbd_running.hbd_totalatoms=my_sz->hbd_frameset[0]->hbd_totalatoms;
        my_sz->hbd_running.hbd_nelem=my_sz->hbd_frameset[0]->hbd_nelem;
        snew(my_sz->hbd_running.myder,my_sz->hbd_frameset[0]->hbd_totalatoms);
        my_sz->hbd_running.hbd_elem=(struct hybrid_elem **)malloc((my_sz->hbd_running.hbd_nelem)*sizeof(struct hybrid_elem * ));
        for (j=0;j<my_sz->hbd_running.hbd_nelem;j++){
                        my_sz->hbd_running.hbd_elem[j]=(struct hybrid_elem *)malloc(sizeof(struct hybrid_elem));
                        my_sz->hbd_running.hbd_elem[j]->cvtype=my_sz->lcvhybrid[j];
                        my_sz->hbd_running.hbd_elem[j]->cvindex=my_sz->lhybrid[j];
                        my_sz->hbd_running.hbd_elem[j]->nder=my_sz->hbd_frameset[0]->hbd_elem[j]->nder;
                        printf("ELEM:  CVTYPE %d CVINDEX %d\n",my_sz->hbd_running.hbd_elem[j]->cvtype+1,my_sz->hbd_running.hbd_elem[j]->cvindex);
                        snew(my_sz->hbd_running.hbd_elem[j]->myder,colvar.natoms[my_sz->hbd_running.hbd_elem[j]->cvindex]);
        }
        printf("HBD_RUNNING: NELEM %d NATOM %d\n",my_sz->hbd_running.hbd_nelem, my_sz->hbd_running.hbd_totalatoms);
       return my_sz->hbd_frameset[0]->hbd_totalatoms;
};

// ------------------------------------------------------------------------------------------------


int PREFIX  hbd_read_dist (  FILE *myfile, FILE *fplog, struct hybrid_elem *elem ){
    FILE *initpos;
    char string[400],remark[10],end[10],*str,type[10];
    initpos=myfile;
    int done=0;
    int rdone=0;
    double uno;
    // scan a remark line and check if the type corresponds (zero level of idiotproofing)
    while(1){
       readagain:
       str=fgets(string,100,myfile);
       //printf("LINE %s ",string);
       sscanf(str,"%3s",end);if(strstr(end,"END")!=NULL){break;}; 
       sscanf(str,"%6s",remark);
       if(strstr(remark,"REMARK")!=NULL){
            sscanf(str,"%s %s",remark,type); 
            if(strstr(type,"DIST")!=NULL){printf("TYPE CONFIRMED: DIST \n");rdone=1;}
            goto readagain;
       };
       if(!done){
          sscanf(str,"%lf",&uno);done=1; elem->ref_dist = (real) uno; printf("FOUND REFERENCE DISTANCE %f\n",elem->ref_dist); 
       }
    } 
    if(rdone==0){printf("MISSING REMARK FIELD IN PATH FRAMES FOR CV CONFIRM\n");EXIT();}; 
    if(done==0){printf("MISSING END FIELD IN PATH FRAMES\n");EXIT();}; 
    printf("READ CORRECTLY ENDED\n");
    // access colvar for making a correct allocation of the derivative 
    elem->nder=colvar.natoms[elem->cvindex]; 
    printf("NATOMS INVOLVED IN CV %d INDEX :",elem->nder); 
    snew(elem->myder, elem->nder);   
    snew(elem->backtable, elem->nder);   
    int i;  
    for(i=0;i<elem->nder;i++){
        // backtable containt the index of the atoms on which you calculate the property
        elem->backtable[i]=colvar.cvatoms[elem->cvindex][i]; 
        printf("%d ",elem->backtable[i]+1);
    }
    printf("\n");
    return elem->nder ;
};
int PREFIX  hbd_read_angle (  FILE *myfile, FILE *fplog, struct hybrid_elem *elem ){
    FILE *initpos;
    char string[400],remark[10],end[10],*str,type[10];
    initpos=myfile;
    int done=0;
    int rdone=0;
    double uno;
    // scan a remark line and check if the type corresponds (zero level of idiotproofing)
    while(1){
       readagain:
       str=fgets(string,100,myfile);
       //printf("LINE %s ",string);
       sscanf(str,"%3s",end);if(strstr(end,"END")!=NULL){break;}; 
       sscanf(str,"%6s",remark);
       if(strstr(remark,"REMARK")!=NULL){
            sscanf(str,"%s %s",remark,type); 
            if(strstr(type,"ANGLE")!=NULL){printf("TYPE CONFIRMED: ANGLE \n");rdone=1;}
            goto readagain;
       };
       if(!done){
          sscanf(str,"%lf",&uno);done=1; elem->ref_dist = (real) uno; printf("FOUND REFERENCE ANGLE %f\n",elem->ref_dist); 
       }
    } 
    if(rdone==0){printf("MISSING REMARK FIELD IN PATH FRAMES FOR CV CONFIRM\n");EXIT();}; 
    if(done==0){printf("MISSING END FIELD IN PATH FRAMES\n");EXIT();}; 
    printf("READ CORRECTLY ENDED\n");
    // access colvar for making a correct allocation of the derivative 
    elem->nder=colvar.natoms[elem->cvindex]; 
    printf("NATOMS INVOLVED IN CV %d INDEX :",elem->nder); 
    snew(elem->myder, elem->nder);   
    snew(elem->backtable, elem->nder);   
    int i;  
    for(i=0;i<elem->nder;i++){
        // backtable containt the index of the atoms on which you calculate the property
        elem->backtable[i]=colvar.cvatoms[elem->cvindex][i]; 
        printf("%d ",elem->backtable[i]+1);
    }
    printf("\n");
    return elem->nder ;
};
int PREFIX  hbd_read_torsion (  FILE *myfile, FILE *fplog, struct hybrid_elem *elem ){
    FILE *initpos;
    char string[400],remark[10],end[10],*str,type[10];
    initpos=myfile;
    int done=0;
    int rdone=0;
    double uno;
    // scan a remark line and check if the type corresponds (zero level of idiotproofing)
    while(1){
       readagain:
       str=fgets(string,100,myfile);
       //printf("LINE %s ",string);
       sscanf(str,"%3s",end);if(strstr(end,"END")!=NULL){break;}; 
       sscanf(str,"%6s",remark);
       if(strstr(remark,"REMARK")!=NULL){
            sscanf(str,"%s %s",remark,type); 
            if(strstr(type,"TORSION")!=NULL){printf("TYPE CONFIRMED: TORSION \n");rdone=1;}
            goto readagain;
       };
       if(!done){
          sscanf(str,"%lf",&uno);done=1; elem->ref_dist = (real) uno; printf("FOUND REFERENCE TORSION %f\n",elem->ref_dist); 
       }
    } 
    if(rdone==0){printf("MISSING REMARK FIELD IN PATH FRAMES FOR CV CONFIRM\n");EXIT();}; 
    if(done==0){printf("MISSING END FIELD IN PATH FRAMES\n");EXIT();}; 
    printf("READ CORRECTLY ENDED\n");
    // access colvar for making a correct allocation of the derivative 
    elem->nder=colvar.natoms[elem->cvindex]; 
    printf("NATOMS INVOLVED IN CV %d INDEX :",elem->nder); 
    snew(elem->myder, elem->nder);   
    snew(elem->backtable, elem->nder);   
    int i;  
    for(i=0;i<elem->nder;i++){
        // backtable containt the index of the atoms on which you calculate the property
        elem->backtable[i]=colvar.cvatoms[elem->cvindex][i]; 
        printf("%d ",elem->backtable[i]+1);
    }
    printf("\n");
    return elem->nder ;
};
int PREFIX  hbd_read_coord (  FILE *myfile, FILE *fplog, struct hybrid_elem *elem ){
    FILE *initpos;
    char string[400],remark[10],end[10],*str,type[10];
    initpos=myfile;
    int done=0;
    int rdone=0;
    double uno;
    // scan a remark line and check if the type corresponds (zero level of idiotproofing)
    while(1){
       readagain:
       str=fgets(string,100,myfile);
       //printf("LINE %s ",string);
       sscanf(str,"%3s",end);if(strstr(end,"END")!=NULL){break;}; 
       sscanf(str,"%6s",remark);
       if(strstr(remark,"REMARK")!=NULL){
            sscanf(str,"%s %s",remark,type); 
            if(strstr(type,"COORD")!=NULL){printf("TYPE CONFIRMED: COORD \n");rdone=1;}
            goto readagain;
       };
       if(!done){
          sscanf(str,"%lf",&uno); done=1; elem->ref_dist = (real) uno; printf("FOUND REFERENCE COORD %f\n",elem->ref_dist); 
       }
    } 
    if(rdone==0){printf("MISSING REMARK FIELD IN PATH FRAMES FOR CV CONFIRM\n");EXIT();}; 
    if(done==0){printf("MISSING END FIELD IN PATH FRAMES\n");EXIT();}; 
    printf("READ CORRECTLY ENDED\n");
    // access colvar for making a correct allocation of the derivative 
    elem->nder=colvar.natoms[elem->cvindex]; 
    printf("NATOMS INVOLVED IN CV %d INDEX :",elem->nder); 
    snew(elem->myder, elem->nder);   
    snew(elem->backtable, elem->nder);   
    int i;  
    for(i=0;i<elem->nder;i++){
        // backtable containt the index of the atoms on which you calculate the property
        elem->backtable[i]=colvar.cvatoms[elem->cvindex][i]; 
        printf("%d ",elem->backtable[i]+1);
    }
    printf("\n");
    return elem->nder ;
};
int PREFIX  hbd_read_target (  FILE *myfile, FILE *fplog, struct hybrid_elem *elem , struct mtd_data_s *mtd_data){
    FILE *initpos;
    char string[400],remark[10],end[10],*str,type[10];
    initpos=myfile;
    int done=0;
    int rdone=0;
    int type_of_target,i_c;
    struct sz_data *pmy_sz;
    i_c=elem->cvindex; 
    pmy_sz=&my_sz_list[ic_to_sz[i_c]];
    printf("|- TARGET TYPE IS %s\n",pmy_sz->path_type); 
     
    // scan a remark line and check if the type corresponds (zero level of idiotproofing)
    while(1){
       readagain:
       str=fgets(string,100,myfile);
       //printf("LINE %s ",string);
       sscanf(str,"%3s",end);if(strstr(end,"END")!=NULL){break;}; 
       sscanf(str,"%6s",remark);
       if(strstr(remark,"REMARK")!=NULL){
            sscanf(str,"%s %s",remark,type); 
            if(strstr(type,"TARGET")!=NULL){printf("TYPE CONFIRMED: TARGET \n");rdone=1;}
            goto readagain;
       };
       if(!done){
            // map or pdb ? 
            if(strcmp(pmy_sz->path_type,"RMSD") == 0){
               printf("|- READER NOT YET IMPLEMENTED FOR TYPE %s \n",pmy_sz->path_type);
                 EXIT();
            } else {
               printf("|- READER NOT YET IMPLEMENTED FOR TYPE %s \n",pmy_sz->path_type);
                 EXIT();
            }
            //sscanf(str,"%lf",&elem->ref_dist);done=1; printf("FOUND REFERENCE TARGET VALUE %f\n",elem->ref_dist); 
       }
    } 
    if(rdone==0){printf("MISSING REMARK FIELD IN PATH FRAMES FOR CV CONFIRM\n");EXIT();}; 
    if(done==0){printf("MISSING END FIELD IN PATH FRAMES\n");EXIT();}; 
    printf("READ CORRECTLY ENDED\n");
    // access colvar for making a correct allocation of the derivative 
    elem->nder=colvar.natoms[elem->cvindex]; 
    printf("NATOMS INVOLVED IN CV %d INDEX :",elem->nder); 
    snew(elem->myder, elem->nder);   
    snew(elem->backtable, elem->nder);   
    int i;  
    for(i=0;i<elem->nder;i++){
        // backtable containt the index of the atoms on which you calculate the property
        elem->backtable[i]=colvar.cvatoms[elem->cvindex][i]; 
        printf("%d ",elem->backtable[i]+1);
    }
    printf("\n");
    return elem->nder ;
};

int  PREFIX  hbd_collect_config ( struct hybrid_frameset running  ) {
int i;
//   printf("NOW COLLECTING REFERENCE VALUES AND DERIVATIVE FOR EACH STRUCT \n");
   for (i=0;i< running.hbd_nelem;i++){
//           printf("ELEM CV %d TYPE %d\n",running.hbd_elem[i]->cvindex,running.hbd_elem[i]->cvtype);
           // retrieve the value and the derivative:
           switch (running.hbd_elem[i]->cvtype){  
                // in case you have a distance you need the value and the derivatives : use chain rules
               case 1:  hbd_copy_dist(running.hbd_elem[i]); break ;
                // in case you have a rmsd you just have to retrieve the configurations : dont use chain rules. The metrics does it all 
               case 3:  hbd_copy_coord(running.hbd_elem[i]); break ;
                // in case you have a rmsd you just have to retrieve the configurations : dont use chain rules. The metrics does it all 
               case 4:  hbd_copy_angle(running.hbd_elem[i]); break ;
                // in case you have a rmsd you just have to retrieve the configurations : dont use chain rules. The metrics does it all 
               case 5:  hbd_copy_torsion(running.hbd_elem[i]); break ;
                // in case you have a rmsd you just have to retrieve the configurations : dont use chain rules. The metrics does it all 
                default:  printf("|- NOT IMPLEMENTED: now dying...\n"); EXIT();
           }
           //printf("RUNNING VAL %f\n",running.hbd_elem[i]->ref_dist);
   }
//   printf("DONE!\n");
   return 1; 
}

int  PREFIX  hbd_copy_dist ( struct hybrid_elem *elem ) {

          int i,cvindex=elem->cvindex;

          dist_restraint(cvindex, &mtd_data);  
         
          elem->ref_dist=colvar.ss0[cvindex];
          for(i=0;i<colvar.natoms[cvindex];i++){
               elem->myder[i][0]=colvar.myder[cvindex][i][0]; 
               elem->myder[i][1]=colvar.myder[cvindex][i][1]; 
               elem->myder[i][2]=colvar.myder[cvindex][i][2]; 
          } 
          return 1; 
}
int  PREFIX  hbd_copy_coord ( struct hybrid_elem *elem ) {

          int i,cvindex=elem->cvindex;

          coord_restraint(cvindex, &mtd_data);  
         
          elem->ref_dist=colvar.ss0[cvindex];
          for(i=0;i<colvar.natoms[cvindex];i++){
               elem->myder[i][0]=colvar.myder[cvindex][i][0]; 
               elem->myder[i][1]=colvar.myder[cvindex][i][1]; 
               elem->myder[i][2]=colvar.myder[cvindex][i][2]; 
          } 
          return 1; 
}
int  PREFIX  hbd_copy_angle ( struct hybrid_elem *elem ) {

          int i,cvindex=elem->cvindex;

          angle_restraint(cvindex, &mtd_data);  
         
          elem->ref_dist=colvar.ss0[cvindex];
          for(i=0;i<colvar.natoms[cvindex];i++){
               elem->myder[i][0]=colvar.myder[cvindex][i][0]; 
               elem->myder[i][1]=colvar.myder[cvindex][i][1]; 
               elem->myder[i][2]=colvar.myder[cvindex][i][2]; 
          } 
          return 1; 
}
int  PREFIX  hbd_copy_torsion ( struct hybrid_elem *elem ) {

          int i,cvindex=elem->cvindex;

          torsion_restraint(cvindex, &mtd_data);  
         
          elem->ref_dist=colvar.ss0[cvindex];
          for(i=0;i<colvar.natoms[cvindex];i++){
               elem->myder[i][0]=colvar.myder[cvindex][i][0]; 
               elem->myder[i][1]=colvar.myder[cvindex][i][1]; 
               elem->myder[i][2]=colvar.myder[cvindex][i][2]; 
          } 
          return 1; 
}
int  PREFIX  hbd_metrics (struct hybrid_frameset *running, struct hybrid_frameset *reference , struct cmap_outpack *outpack, real **mat) {
          int n=running->hbd_nelem;
          real dist=0.;
          // do only diagonal terms for the time being...
          int i,jj; 
         real *store=(real *)malloc(n*sizeof(real)); 
         int  *start=(int *)malloc(n*sizeof(int)); 
//          real store [10]; 
//          int start [10]; 
          for(i=0;i<n;i++){
                   switch(running->hbd_elem[i]->cvtype){
                       case 1: 		store[i]=hbd_metrics_dist ( running->hbd_elem[i] , reference->hbd_elem[i]);
                                        break;
                       case 3: 		store[i]=hbd_metrics_coord ( running->hbd_elem[i] , reference->hbd_elem[i]);
                                        break;
                       case 4: 		store[i]=hbd_metrics_angle ( running->hbd_elem[i] , reference->hbd_elem[i]);
                                        break;
                       case 5: 		store[i]=hbd_metrics_torsion ( running->hbd_elem[i] , reference->hbd_elem[i]);
                                        break;
                       default: 	printf("|- NOT IMPLEMENTED: now dying...\n"); EXIT();
                   }

          }
          // calculate the matrix with nondiagonal representation
          // first the diagonal contributions 
          int j=0;
          for(i=0;i<n;i++){
                   start[i]=j;
                   for(jj=0;jj<running->hbd_elem[i]->nder;jj++){
                       outpack->derr_dr0[0][j+jj]=2.*mat[i][i]*store[i]*reference->hbd_elem[i]->myder[jj][0];
                       outpack->derr_dr0[1][j+jj]=2.*mat[i][i]*store[i]*reference->hbd_elem[i]->myder[jj][1];
                       outpack->derr_dr0[2][j+jj]=2.*mat[i][i]*store[i]*reference->hbd_elem[i]->myder[jj][2];
                       //printf("OUTPACK %d VAL %f %f %f\n",j+jj,outpack->derr_dr0[0][j+jj],outpack->derr_dr0[1][j+jj],outpack->derr_dr0[2][j+jj]);
                       // MISSING THE DERIVATIVE OF THE MATRIX!!!!!!!!!
	           } 
                   dist+=store[i]*store[i]*mat[i][i];
                   j+=jj;
          }
          // non diagonal elements
          for(i=0;i<n-1;i++){
                   for(j=i+1;j<n;j++){
                       dist+=2.*mat[i][j]*store[i]*store[j];
                       int l=start[i];
                       for(jj=0;jj<running->hbd_elem[i]->nder;jj++){
                           outpack->derr_dr0[0][l+jj]+=2.*mat[i][j]*store[j]*reference->hbd_elem[i]->myder[jj][0];
                           outpack->derr_dr0[1][l+jj]+=2.*mat[i][j]*store[j]*reference->hbd_elem[i]->myder[jj][1];
                           outpack->derr_dr0[2][l+jj]+=2.*mat[i][j]*store[j]*reference->hbd_elem[i]->myder[jj][2];
                         //  printf("NOND %d VAL %f %f %f\n",l+jj,2.*mat[i][j]*store[j]*reference->hbd_elem[i]->myder[jj][0],2.*mat[i][j]*store[j]*reference->hbd_elem[i]->myder[jj][1],2.*mat[i][j]*store[j]*reference->hbd_elem[i]->myder[jj][2]);
                       // MISSING THE DERIVATIVE OF THE MATRIX!!!!!!!!!
	               } 
                       l=start[j];
                       for(jj=0;jj<running->hbd_elem[j]->nder;jj++){
                           outpack->derr_dr0[0][l+jj]+=2.*mat[i][j]*store[i]*reference->hbd_elem[j]->myder[jj][0];
                           outpack->derr_dr0[1][l+jj]+=2.*mat[i][j]*store[i]*reference->hbd_elem[j]->myder[jj][1];
                           outpack->derr_dr0[2][l+jj]+=2.*mat[i][j]*store[i]*reference->hbd_elem[j]->myder[jj][2];
                         //  printf("NOND %d VAL %f %f %f\n",l+jj,2.*mat[i][j]*store[i]*reference->hbd_elem[j]->myder[jj][0],2.*mat[i][j]*store[i]*reference->hbd_elem[j]->myder[jj][1],2.*mat[i][j]*store[i]*reference->hbd_elem[j]->myder[jj][2]);
                       // MISSING THE DERIVATIVE OF THE MATRIX!!!!!!!!!
	               } 
                   }
          }

          free(store);
          free(start);
          outpack->err=dist; 
          return 1; 
};

real PREFIX hbd_metrics_dist ( struct hybrid_elem *run,  struct hybrid_elem *ref ) {
          real dist=0.;
          int i;
          dist=sqrt(pow2(run->ref_dist-ref->ref_dist));
          //transfer the derivative in the reference element so you have 
          //one block of derivative each element 
          real one=run->ref_dist-ref->ref_dist>0.?1.:-1.;
          for (i=0;i<ref->nder;i++ ){
                   ref->myder[i][0]=run->myder[i][0]*one;
                   ref->myder[i][1]=run->myder[i][1]*one;
                   ref->myder[i][2]=run->myder[i][2]*one;
          }
          return dist;       
} 
real PREFIX hbd_metrics_coord ( struct hybrid_elem *run,  struct hybrid_elem *ref ) {
          int i;
          real dist=0.;
          dist=sqrt(pow2(run->ref_dist-ref->ref_dist));
          //transfer the derivative in the reference element so you have 
          //one block of derivative each element 
          real one=run->ref_dist-ref->ref_dist>0.?1.:-1.;
          for (i=0;i<ref->nder;i++ ){
                   ref->myder[i][0]=run->myder[i][0]*one;
                   ref->myder[i][1]=run->myder[i][1]*one;
                   ref->myder[i][2]=run->myder[i][2]*one;
          }
          return dist;       
}
real PREFIX hbd_metrics_angle ( struct hybrid_elem *run,  struct hybrid_elem *ref ) {
          int i;
          real dist=0.;
          // need periodic representation
           // 0.5(1-cos(angle-a_0)) 
          dist=M_PI*0.5*(1.0-cos(run->ref_dist-ref->ref_dist));
          //transfer the derivative in the reference element so you have 
          //one block of derivative each element 
          real one=M_PI*0.5*sin(run->ref_dist-ref->ref_dist);
          for (i=0;i<ref->nder;i++ ){
                   ref->myder[i][0]=run->myder[i][0]*one;
                   ref->myder[i][1]=run->myder[i][1]*one;
                   ref->myder[i][2]=run->myder[i][2]*one;
          }
          return dist;       
}
real PREFIX hbd_metrics_torsion ( struct hybrid_elem *run,  struct hybrid_elem *ref ) {
          real dist=0.;
          int i;
          // need periodic representation
           // 0.5(1-cos(angle-a_0)) 
          dist=M_PI*0.5*(1.0-cos(run->ref_dist-ref->ref_dist));
          //transfer the derivative in the reference element so you have 
          //one block of derivative each element 
          real one=M_PI*0.5*sin(run->ref_dist-ref->ref_dist);
          for (i=0;i<ref->nder;i++ ){
                   ref->myder[i][0]=run->myder[i][0]*one;
                   ref->myder[i][1]=run->myder[i][1]*one;
                   ref->myder[i][2]=run->myder[i][2]*one;
          }
          return dist;       
}
