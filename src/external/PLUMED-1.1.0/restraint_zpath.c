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

void  PREFIX zpath_restraint(int i_c, struct mtd_data_s *mtd_data) {

        int iat,i,j;
        real z,ci_vec,tmp1;
        struct coordinates_frameset *pmy_coord1;
        struct sz_data *pmy_sz;
        struct hybrid_frameset *hbd_pmy;
        struct cmap_inpack inpack;
        struct cmap_outpack outpack;
        real dz_dr0[3][MAXATOMS_PATH];
        real dz_dcm[MAXFRAMES_PATH][MAXDIM_CMAP];
        real dz_dr1[MAXFRAMES_PATH][3][MAXATOMS_PATH];
        real dmsd_dr1[3][MAXATOMS_PATH];
        int   start_avg = 0;
        int    tot_con;    
        real *save_err; 
        int nneigh,ii;

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



        for (i=0;i<3;i++){
                for (j=0;j<colvar.natoms[i_c];j++) {
                        dz_dr0[i][j]=0.;
                }
        }

        ci_vec=0.;

   	for (i=0;i<colvar.natoms[i_c];i++){
            iat = colvar.cvatoms[i_c][i];
       	    inpack.r0[i][0] = mtd_data->pos[iat][0];
            inpack.r0[i][1] = mtd_data->pos[iat][1];
            inpack.r0[i][2] = mtd_data->pos[iat][2];
	}



        if(pmy_sz->umb_on && (colvar.it-pmy_sz->umblagsteps>=0)) start_avg = 1;  

        if(strcmp(pmy_sz->path_type,"CMAP") == 0){
         tot_con=pmy_sz->my_cmap_pack.number+pmy_sz->my_cmap_pack.gnumber;
         cmap_running(&inpack,&pmy_sz->my_cmap_pack);
        }

        if(strcmp(pmy_sz->path_type,"HYBRID") == 0){ 
            hbd_collect_config(pmy_sz->hbd_running);
        }

        for(ii=0;ii< nneigh;ii++){
                i=pmy_sz->lneigh[ii];

                //fprintf(mtd_data->fplog,"PRE RMSD %d %d\n",i,pmy_sz->lneigh[ii]);
                //fflush(mtd_data->fplog);
 
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
                 hbd_metrics(&pmy_sz->hbd_running,hbd_pmy,&outpack,pmy_sz->mathybrid);
                }
 
             //   fprintf(mtd_data->fplog,"ERR %d %f \n",i,outpack.err);
                fflush(mtd_data->fplog);
 
             // in case you are calculating the neigh list 
                if(pmy_sz->neigh==1 && colvar.it%pmy_sz->neigh_time==0)save_err[i]=outpack.err;

             // sqrt option
                if(pmy_sz->sqrt_on){
                   if(outpack.err<1.e-6){fprintf(mtd_data->fplog,"|-ERRROR IN ZPATH: too small error %f\n",outpack.err);EXIT();}
                   tmp1=(0.5/sqrt(outpack.err))*exp(-pmy_sz->lambda*sqrt(outpack.err));
                }else{
                   tmp1=exp(-pmy_sz->lambda*outpack.err);
                }
                ci_vec+=tmp1;
                for(j=0;j<colvar.natoms[i_c];j++){
                        dz_dr0[0][j]+=tmp1*outpack.derr_dr0[0][j];
                        dz_dr0[1][j]+=tmp1*outpack.derr_dr0[1][j];
                        dz_dr0[2][j]+=tmp1*outpack.derr_dr0[2][j];
                        if((strcmp(pmy_sz->path_type,"RMSD") == 0 || strcmp(pmy_sz->path_type,"DRMS") == 0) && start_avg){ 
                         dz_dr1[i][0][j]=(dmsd_dr1[0][j])*tmp1;
                         dz_dr1[i][1][j]=(dmsd_dr1[1][j])*tmp1;
                         dz_dr1[i][2][j]=(dmsd_dr1[2][j])*tmp1;
                        }
                }
                if(strcmp(pmy_sz->path_type,"CMAP") == 0 && start_avg){ 
                 for(j=0;j<tot_con;j++){
                          dz_dcm[i][j]=(outpack.derr_dcm[j])*tmp1;
                        }
                 }
        }


        z=-(1./pmy_sz->lambda)*log(ci_vec);

        for(j=0;j<colvar.natoms[i_c];j++){
                dz_dr0[0][j]=(1./ci_vec) *dz_dr0[0][j];
                dz_dr0[1][j]=(1./ci_vec) *dz_dr0[1][j];
                dz_dr0[2][j]=(1./ci_vec) *dz_dr0[2][j];
                 if((strcmp(pmy_sz->path_type,"RMSD") == 0 || strcmp(pmy_sz->path_type,"DRMS") == 0) && start_avg){ 
                  for(ii=0;ii<nneigh;ii++){
                       i=pmy_sz->lneigh[ii];
                       dz_dr1[i][0][j]=dz_dr1[i][0][j]/ci_vec;
                       dz_dr1[i][1][j]=dz_dr1[i][1][j]/ci_vec;
                       dz_dr1[i][2][j]=dz_dr1[i][2][j]/ci_vec;
                  }
                }
        }

         if(strcmp(pmy_sz->path_type,"CMAP") == 0 && start_avg){ 
           for(ii=0;ii<nneigh;ii++){
               i=pmy_sz->lneigh[ii];
               for(j=0;j<tot_con;j++){
                  dz_dcm[i][j]=dz_dcm[i][j]/ci_vec;
               }
           }
        }


        for(i=0;i<colvar.natoms[i_c];i++) {
          colvar.myder[i_c][i][0] = dz_dr0[0][i];
          colvar.myder[i_c][i][1] = dz_dr0[1][i];
          colvar.myder[i_c][i][2] = dz_dr0[2][i];
        }

        colvar.ss0[i_c]=z;

        // neigh list? do quicksort and deallocate save_err
        if(pmy_sz->neigh==1 && colvar.it%pmy_sz->neigh_time==0){
             //   for(i=0;i<nneigh;i++)printf("BEFORE SORTING %d %f\n",pmy_sz->lneigh[i],save_err[i]);
                realquicksort(save_err,pmy_sz->lneigh,0,nneigh-1);
             //   for(i=0;i<nneigh;i++)printf("AFTER SORTING %d %f\n",pmy_sz->lneigh[i],save_err[i]);
                free(save_err);
        }


        if(pmy_sz->umb_on==1){
          if(strcmp(pmy_sz->path_type,"CMAP") == 0) mean_map(pmy_sz,dz_dcm,i_c,mtd_data->fplog);
          if(strcmp(pmy_sz->path_type,"RMSD") == 0 || strcmp(pmy_sz->path_type,"DRMS") == 0) mean_rmsd(pmy_sz,dz_dr1,i_c,mtd_data->fplog);
        }

 
        return;
}


// ------------------------------------------------------------------------------------------------

 void PREFIX mean_rmsd(struct sz_data *pmy_sz, real dCV_dr1[MAXFRAMES_PATH][3][MAXATOMS_PATH],
                        int i_c, FILE *fplog){

   real tmp1,tmp2,tmp3,tmp4;
   int i,j,k,l;  
   FILE   *fp;


          if(logical.steer[i_c] != 1){
            fprintf(fplog,"ERROR:  UMBRELLA MUST BE USED TOGHETER WITH STEER KEYWORD \n ");
            EXIT();
          } 

          if(!colvar.it){
                   kill_me[i_c]=0;
                   for(i=0;i<3;i++){
                     for(j=0;j<(*pmy_sz->frameset[0]).natoms;j++){
                       for(k=0;k< pmy_sz->number;k++){
                         for(l=0;l<pmy_sz->umbblocksize;l++){
                          pmy_sz->umb_block[i][j][k][l]=0.;
                         }
                       }
                     }
                   }  
                   for(i=0;i<3;i++){
                     for(j=0;j<(*pmy_sz->frameset[0]).natoms;j++){
                       for(k=0;k< pmy_sz->number;k++){
                          pmy_sz->umb_avg[i][j][k]=0.;
                       }
                     }
                   }  
                   pmy_sz->umbcount=0;
          }
          // END OF INITIZIALIZATION
          if(  (colvar.it-pmy_sz->umblagsteps>=0)  &&  ((colvar.it-pmy_sz->umblagsteps)%pmy_sz->umbstride==0)  ){
             pmy_sz->umbcount++;
             l=pmy_sz->umbcount%pmy_sz->umbblocksize;
             //cout<<"UMBCOUNT "<<pmy_sz->umbcount<<endl;
             for(i=0;i<3;i++){
                for(j=0;j<(*pmy_sz->frameset[0]).natoms;j++){
                  for(k=0;k< pmy_sz->number;k++){
                     tmp1=2.*cvsteer.spring[i_c]*(colvar.ss0[i_c]-cvsteer.pos[i_c])*dCV_dr1[k][i][j];    
                     pmy_sz->umb_avg[i][j][k]      +=tmp1 ; 
                     pmy_sz->umb_block[i][j][k][l]  =pmy_sz->umb_avg[i][j][k]/((real) pmy_sz->umbcount) ; 
                  }
                }
             }  
             /*evaluate convergence only when the block is full*/ 
             if( (pmy_sz->umbcount>=pmy_sz->umbblocksize) && (pmy_sz->umbcount%pmy_sz->umbblocksize == 0) ){
                     
                     tmp3=0.;
                     tmp4=0.;
                     for(i=0;i<3;i++){
                         for(j=0;j<(*pmy_sz->frameset[0]).natoms;j++){
                             for(k=0;k< pmy_sz->number;k++){
                               tmp3+=pow(pmy_sz->umb_avg[i][j][k]/((real) pmy_sz->umbcount),2);    
                             }
                         }
                     } 
                     tmp3=sqrt(tmp3);
                     tmp2=0.;
                     for(l=0;l<pmy_sz->umbblocksize;l++){
                       tmp1=0.;
                       tmp4=0.; 
                       for(i=0;i<3;i++){
                         for(j=0;j<(*pmy_sz->frameset[0]).natoms;j++){
                             for(k=0;k< pmy_sz->number;k++){
                               tmp1+=pmy_sz->umb_block[i][j][k][l]*(pmy_sz->umb_avg[i][j][k]/((real) pmy_sz->umbcount));    
                               tmp4+=pow(pmy_sz->umb_block[i][j][k][l],2);    
                             }
                         }
                       } 
                       tmp4=sqrt(tmp4);
                       tmp2+=pow(tmp1/(tmp3*tmp4),2);
                     } 
                     tmp2=tmp2/((real) pmy_sz->umbblocksize); // mean of the projections
                     /* FILE PRINTOUT */
                     if(pmy_sz->umbcount==pmy_sz->umbblocksize){
                      if(colvar.type_s[i_c] == 30) fp=fopen((mtd_data.ionode?"umbrella_evol_s.dat":"/dev/null"),"w");
                      if(colvar.type_s[i_c] == 31) fp=fopen((mtd_data.ionode?"umbrella_evol_a.dat":"/dev/null"),"w");
                     } else {
                      if(colvar.type_s[i_c] == 30) fp=fopen((mtd_data.ionode?"umbrella_evol_s.dat":"/dev/null"),"a");
                      if(colvar.type_s[i_c] == 31) fp=fopen((mtd_data.ionode?"umbrella_evol_a.dat":"/dev/null"),"a");
                     }
 
                     if(fp!=NULL){
                          for(k=0;k<pmy_sz->number ;k++){
                                fprintf(fp,"%% TRAJ %d\n",k+1);
                                for(j=0;j<(*pmy_sz->frameset[0]).natoms;j++){
                                    fprintf(fp,"%15.8e %15.8e %15.8e\n",pmy_sz->umb_avg[0][j][k]/((real) pmy_sz->umbcount),
                                                                        pmy_sz->umb_avg[1][j][k]/((real) pmy_sz->umbcount),
                                                                        pmy_sz->umb_avg[2][j][k]/((real) pmy_sz->umbcount));
				}
                          }
                          fprintf(fp,"\n\n");
                          fclose(fp);
                     }
                     else { fprintf(fplog,"ERROR IN OPENING THE UMBRELLA FILE \n"); }

                     // criteria for stop
                     if(colvar.type_s[i_c] == 30)  fprintf(fplog,"PROJECTION_SS VALUE %lf PERM %d \n",sqrt((1.-tmp2)*(1.-tmp2)),pmy_sz->countperm);
                     if(colvar.type_s[i_c] == 31)  fprintf(fplog,"PROJECTION_ZZ VALUE %lf PERM %d \n",sqrt((1.-tmp2)*(1.-tmp2)),pmy_sz->countperm); 
                     if( sqrt((1.-tmp2)*(1.-tmp2))<pmy_sz->umbtolerance){
                       pmy_sz->countperm++;
                       if(pmy_sz->countperm >= pmy_sz->umbpermanency){
                              kill_me[i_c]=1;
                              /* Decide wether kill the calculation or wait for others  check kill_me vector*/
                              i=0;
                              j=0;
                              for(k=0;k<nsz;k++){
                                    if(kill_me[k]>=0)i++;
                                    if(kill_me[k]>0)j++;
                              }
                              // IF KILL_ME == NUMBER OF VARIABLES ON WHICH CALCULATE THE VALUE THEN DUMP AND KILL
                              if(i==j) {
                               fprintf(fplog,"** VECTOR CONVERGED !! Exiting... \n");
                               EXIT(); 
                              } 
                       }  
                     }
                     else{
                       pmy_sz->countperm=0;
                       kill_me[i_c]=0;
                     }
           
             } 

          }

      return;          
}

// ------------------------------------------------------------------------------------------------

 void PREFIX mean_map(struct sz_data *pmy_sz, real dCV_dcm[MAXFRAMES_PATH][MAXDIM_CMAP],
                        int i_c, FILE *fplog){

   real tmp1,tmp2,tmp3,tmp4;
   int    i,j,k,l;  
   FILE   *fp;
   int    tot_con;



          if(logical.steer[i_c] != 1){
            fprintf(fplog,"ERROR:  UMBRELLA MUST BE USED TOGHETER WITH STEER KEYWORD \n ");
            EXIT();
          } 

          tot_con=pmy_sz->my_cmap_pack.number+pmy_sz->my_cmap_pack.gnumber;
          
          if(!colvar.it){
                   kill_me[i_c]=0;
                     for(j=0;j<tot_con;j++){
                       for(k=0;k<pmy_sz->number;k++){
                         for(l=0;l<pmy_sz->umbblocksize;l++){
                          pmy_sz->umb_map_block[j][k][l]=0.;
                         }
                       }
                     }
                     for(j=0;j<tot_con;j++){
                       for(k=0;k<pmy_sz->number;k++){
                          pmy_sz->umb_map_avg[j][k]=0.;
                       }
                     }
                   pmy_sz->umbcount=0;
          }
          // END OF INITIALIZATION
          if(  (colvar.it-pmy_sz->umblagsteps>=0)  &&  ((colvar.it-pmy_sz->umblagsteps)%pmy_sz->umbstride==0)  ){
                                  
             pmy_sz->umbcount++;
             l=pmy_sz->umbcount%pmy_sz->umbblocksize;

                for(j=0;j<tot_con;j++){
                  for(k=0;k<pmy_sz->number;k++){
                     tmp1=cvsteer.spring[i_c]*(colvar.ss0[i_c]-cvsteer.pos[i_c])*dCV_dcm[k][j];    
                     pmy_sz->umb_map_avg[j][k]      +=tmp1 ; 
                     pmy_sz->umb_map_block[j][k][l]  =pmy_sz->umb_map_avg[j][k]/((real) pmy_sz->umbcount) ; 
                  }
                }
             /*evaluate convergence only when the block is full*/ 
             if( (pmy_sz->umbcount>=pmy_sz->umbblocksize) && (pmy_sz->umbcount%pmy_sz->umbblocksize == 0) ){
                     tmp3=0.;
                     tmp4=0.;
                         for(j=0;j<tot_con;j++){
                             for(k=0;k<pmy_sz->number;k++){
                               tmp3+=pow(pmy_sz->umb_map_avg[j][k]/((real) pmy_sz->umbcount),2);    
                             }
                         }
                     tmp3=sqrt(tmp3);
                     tmp2=0.;
                     for(l=0;l<pmy_sz->umbblocksize;l++){
                       tmp1=0.;
                       tmp4=0.; 
                         for(j=0;j<tot_con;j++){
                             for(k=0;k<pmy_sz->number;k++){
                               tmp1+=pmy_sz->umb_map_block[j][k][l]*(pmy_sz->umb_map_avg[j][k]/((real) pmy_sz->umbcount));    
                               tmp4+=pow(pmy_sz->umb_map_block[j][k][l],2);    
                             }
                         }
                       tmp4=sqrt(tmp4);
                       tmp2+=pow(tmp1/(tmp3*tmp4),2);
                     } 
                     tmp2=tmp2/((real) pmy_sz->umbblocksize); // mean of the projections
                     /* FILE PRINTOUT */
                     if(pmy_sz->umbcount==pmy_sz->umbblocksize){
                      if(colvar.type_s[i_c] == 30) fp=fopen("umbrella_evol_s.dat","w");
                      if(colvar.type_s[i_c] == 31) fp=fopen("umbrella_evol_a.dat","w");
                     } else {
                      if(colvar.type_s[i_c] == 30) fp=fopen("umbrella_evol_s.dat","a");
                      if(colvar.type_s[i_c] == 31) fp=fopen("umbrella_evol_a.dat","a");
                     }
 
                     if(fp!=NULL){
                          for(k=0;k<pmy_sz->number ;k++){
                                fprintf(fp,"%% TRAJ %d\n",k+1);
                                for(j=0;j<tot_con;j++){
                                    fprintf(fp,"%15.8e \n",pmy_sz->umb_map_avg[j][k]/((real) pmy_sz->umbcount));
				}
                          }
                          fprintf(fp,"\n\n");
                          fclose(fp);
                     }
                     else { fprintf(fplog,"ERROR IN OPENING THE UMBRELLA FILE \n"); }

                     // criteria for stop
                     if(colvar.type_s[i_c] == 30)  fprintf(fplog,"PROJECTION_SS VALUE %lf PERM %d \n",sqrt((1.-tmp2)*(1.-tmp2)),pmy_sz->countperm);
                     if(colvar.type_s[i_c] == 31)  fprintf(fplog,"PROJECTION_ZZ VALUE %lf PERM %d \n",sqrt((1.-tmp2)*(1.-tmp2)),pmy_sz->countperm);

                     if( sqrt((1.-tmp2)*(1.-tmp2))<pmy_sz->umbtolerance){
                       pmy_sz->countperm++;
                       if(pmy_sz->countperm >= pmy_sz->umbpermanency){
                              kill_me[i_c]=1;
                              /* Decide wether kill the calculation or wait for others  check kill_me vector*/
                              i=0;
                              j=0;
                              for(k=0;k<nsz;k++){
                                    if(kill_me[k]>=0)i++;
                                    if(kill_me[k]>0)j++;
                              }
                              // IF KILL_ME == NUMBER OF VARIABLES ON WHICH CALCULATE THE VALUE THEN DUMP AND KILL
                              if(i==j) {
                               fprintf(fplog,"** VECTOR CONVERGED !! Exiting... \n");
                               EXIT(); 
                              }
                       }  
                     }
                     else{
                       pmy_sz->countperm=0;
                       kill_me[i_c]=0;
                     }
           
             } 

          }

      return;          
}
