#include "opencl_util_mpi.h" // 注意一定要在 pass_mod_var 之前，以免被宏污染
#include "pass_mod_var.h"
#include "save_load_var.h"
#include "sum_up_whole_potential_shanghui.h"
#include "cmake_help.h"
#include <alloca.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>


#define CL_TARGET_OPENCL_VERSION 200
#ifdef APPLE
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#include <CL/cl_ext.h>
#endif


#undef b0
#undef b2
#undef b4
#undef b6
#undef a_save

extern int MV(mpi_tasks, myid);
extern int MV(mpi_tasks, mpi_platform_relative_id);
#define myid MV(mpi_tasks, myid)
#define mpi_platform_relative_id MV(mpi_tasks, mpi_platform_relative_id)
#define mpi_task_per_gpu MV(mpi_tasks, mpi_task_per_gpu)
#define mpi_per_node MV(mpi_tasks, mpi_per_node)
#define DEVICE_ID 0

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

size_t platformId = 0; // INFO choose a platform

cl_platform_id *platforms;
cl_program program;
cl_context context;
cl_command_queue cQ;
cl_uint numOfPlatforms;
const char *file_names[] = {CURRENT_DIRECTORY "/sum_up.cl", 
                            CURRENT_DIRECTORY "/integrate_first_order_rho.cl"};
                            // "/public/home/autopar/FHI-aims-test/fhi-aims_MPE_O3_local_index_final/src/DFPT_OpenCL/integration_points.cl"};
const char *kernel_names[] = {"sum_up_whole_potential_shanghui_sub_t_"
                              , "integrate_first_order_rho_sub_tmp2_"
                              , "integrate_first_order_h_sub_tmp2_"
                              , "sum_up_whole_potential_shanghui_pre_proc_"};
                              // "c_integration_points2"};

#define number_of_files (sizeof(file_names) / sizeof(file_names[0]))
#define number_of_kernels (sizeof(kernel_names) / sizeof(kernel_names[0]))

cl_kernel kernels[number_of_kernels];
char *buffer[number_of_files];
size_t sizes[number_of_files];

// const char *buffer[] = {ocl_source_sumup, ocl_source_rho_h};
// size_t sizes[] = {sizeof(ocl_source_sumup), sizeof(ocl_source_rho_h)};

size_t localSize[] = {128}; // 可能被重新设置 !!!
size_t globalSize[] = {128 * 128};  // 可能被重新设置 !!!

size_t globalSize_sum_up_pre_proc[1];
size_t localSize_sum_up_pre_proc[1] = {64};
const int i_center_tile_size_default = 256;

extern double *Fp_function_spline_slice;
extern double *Fpc_function_spline_slice;
extern double *Fp;

RHO_PARAM rho_param;
H_PARAM H_param;

static int opencl_init_finished = 0;
static int opencl_common_buffer_init_finished = 0;
static int remember_arg_index_1 = -10;
static int remember_arg_Fp_function_spline = -10;
static int remember_arg_n_my_batches_work_sumup = -10;
static int remember_arg_Fp_max_grid = -10;
static int sum_up_first_begin_finished = 0;
static int sum_up_begin_0_finished = 0;
static int rho_first_begin_finished = 0;
static int H_first_begin_finished = 0;
static int h_begin_0_finished = 0;

struct CL_BUF_COM_T {
  cl_mem species;
  cl_mem empty;
  cl_mem centers_hartree_potential;
  cl_mem center_to_atom;
  cl_mem species_center;
  cl_mem coords_center;
  cl_mem l_hartree;
  cl_mem n_grid;
  cl_mem n_radial;
  cl_mem r_grid_min;
  cl_mem r_grid_inc;
  cl_mem log_r_grid_inc;
  cl_mem scale_radial;
  cl_mem r_radial;
  cl_mem r_grid;
  cl_mem n_cc_lm_ijk;
  cl_mem index_cc;
  cl_mem index_ijk_max_cc;
  cl_mem b0;
  cl_mem b2;
  cl_mem b4;
  cl_mem b6;
  cl_mem a_save;
  cl_mem Fp_function_spline_slice;
  cl_mem Fpc_function_spline_slice;
  // ---
  cl_mem rho_multipole_index;
  cl_mem compensation_norm;
  cl_mem compensation_radius;
  cl_mem rho_multipole_h_p_s;
  cl_mem multipole_radius_free;
  // rho global
  cl_mem perm_basis_fns_spl;
  cl_mem outer_radius_sq;
  cl_mem basis_fn;
  cl_mem basis_l;
  cl_mem atom_radius_sq;
  cl_mem basis_fn_start_spl;
  cl_mem basis_fn_atom;
  cl_mem basis_wave_ordered;
  cl_mem basis_kinetic_ordered;
  cl_mem Cbasis_to_basis;
  cl_mem Cbasis_to_center;
  cl_mem centers_basis_integrals; // 可能因为宏展开出问题
  cl_mem index_hamiltonian;
  cl_mem position_in_hamiltonian;
  cl_mem column_index_hamiltonian;
  // pbc_lists_coords_center 即 coords_center，只是为了规避一点点重名
  cl_mem center_to_cell;
  // loop helper
  cl_mem point_to_i_batch;
  cl_mem point_to_i_index;
  cl_mem valid_point_to_i_full_point;
  cl_mem index_cc_aos;
  cl_mem i_center_to_centers_index;
  // sum up
  // cl_mem partition_tab_std;
  // cl_mem delta_v_hartree;               // (n_full_points_work)
  // cl_mem rho_multipole;                 // (n_full_points_work)
  cl_mem centers_rho_multipole_spl;     // (l_pot_max+1)**2, n_max_spline, n_max_radial+2, n_atoms)
  cl_mem centers_delta_v_hart_part_spl; // (l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid, n_atoms)
  // cl_mem adap_outer_radius_sq;          // (n_atoms)
  // cl_mem multipole_radius_sq;           // (n_atoms)
  // cl_mem l_hartree_max_far_distance;    // (n_atoms)
  // cl_mem outer_potential_radius;        // (0:l_pot_max, n_atoms)
  // cl_mem multipole_c;                   // (n_cc_lm_ijk(l_pot_max), n_atoms)
  // sum up tmp
  cl_mem angular_integral_log; // per block (l_pot_max + 1) * (l_pot_max + 1) * n_max_grid
  cl_mem Fp; // global_size * (l_pot_max + 2) * n_centers_hartree_potential)
  cl_mem coord_c;
  cl_mem coord_mat;
  cl_mem rest_mat;
  cl_mem vector;
  cl_mem delta_v_hartree_multipole_component;
  cl_mem rho_multipole_component;
  cl_mem ylm_tab;
  // sum_up batches
  // cl_mem batches_size_sumup;
  // cl_mem batches_points_coords_sumup;
  // rho tmp
  cl_mem dist_tab_sq__;
  cl_mem dist_tab__;
  cl_mem dir_tab__;
  cl_mem atom_index__;
  cl_mem atom_index_inv__;
  cl_mem i_basis_fns__;
  cl_mem i_basis_fns_inv__;
  cl_mem i_atom_fns__;
  cl_mem spline_array_start__;
  cl_mem spline_array_end__;
  cl_mem one_over_dist_tab__;
  cl_mem rad_index__;
  cl_mem wave_index__;
  cl_mem l_index__;
  cl_mem l_count__;
  cl_mem fn_atom__;
  cl_mem zero_index_point__;
  cl_mem wave__;
  cl_mem first_order_density_matrix_con__;
  cl_mem i_r__;
  cl_mem trigonom_tab__;
  cl_mem radial_wave__;
  cl_mem spline_array_aux__;
  cl_mem aux_radial__;
  cl_mem ylm_tab__;
  cl_mem dylm_dtheta_tab__;
  cl_mem scaled_dylm_dphi_tab__;
  cl_mem kinetic_wave__;
  cl_mem grid_coord__;
  cl_mem H_times_psi__;
  cl_mem T_plus_V__;
  cl_mem contract__;
  cl_mem wave_t__;
  cl_mem first_order_H_dense__;
  // rho batches
  cl_mem batches_size_rho;
  cl_mem batches_batch_n_compute_rho;
  cl_mem batches_batch_i_basis_rho;
  cl_mem batches_points_coords_rho;
  // H batches
  cl_mem batches_size_H;
  cl_mem batches_batch_n_compute_H;
  cl_mem batches_batch_i_basis_H;
  cl_mem batches_points_coords_H;
  // rho
  cl_mem basis_l_max__;
  cl_mem n_points_all_batches__;
  cl_mem n_batch_centers_all_batches__;
  cl_mem batch_center_all_batches__;
  cl_mem batch_point_to_i_full_point__;
  cl_mem ins_idx_all_batches__;
  cl_mem first_order_rho__;
  cl_mem first_order_density_matrix__;
  cl_mem partition_tab__;
  // cl_mem tmp_rho__; // only for swcl
  // H
  cl_mem batches_batch_i_basis_h__;
  cl_mem partition_all_batches__;
  cl_mem first_order_H__;
  cl_mem local_potential_parts_all_points__;
  cl_mem local_first_order_rho_all_batches__;
  cl_mem local_first_order_potential_all_batches__;
  cl_mem local_dVxc_drho_all_batches__;
  cl_mem local_rho_gradient__;
  cl_mem first_order_gradient_rho__;
} cl_buf_com;

typedef struct CL_BUF_SUMUP_T {
  // sum up
  cl_mem partition_tab_std;
  cl_mem delta_v_hartree;               // (n_full_points_work)
  cl_mem rho_multipole;                 // (n_full_points_work)
  // cl_mem centers_rho_multipole_spl;     // (l_pot_max+1)**2, n_max_spline, n_max_radial+2, n_atoms)
  // cl_mem centers_delta_v_hart_part_spl; // (l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid, n_atoms)
  cl_mem adap_outer_radius_sq;          // (n_atoms)
  cl_mem multipole_radius_sq;           // (n_atoms)
  cl_mem l_hartree_max_far_distance;    // (n_atoms)
  cl_mem outer_potential_radius;        // (0:l_pot_max, n_atoms)
  cl_mem multipole_c;                   // (n_cc_lm_ijk(l_pot_max), n_atoms)
  // sum_up batches
  cl_mem batches_size_sumup;
  cl_mem batches_points_coords_sumup;

  cl_mem point_to_i_batch;
  cl_mem point_to_i_index;
  cl_mem valid_point_to_i_full_point;
} CL_BUF_SUMUP;

CL_BUF_SUMUP cl_buf_sumup[16];

#define IF_ERROR_EXIT(cond, err_code, str)                                                                             \
  if (cond) {                                                                                                          \
    printf("Error! rank%d, %s:%d, %s\nError_code=%d\n%s\n", myid, __FILE__, __LINE__, __FUNCTION__, err_code, str);   \
    fflush(stdout);                                                                                                    \
    exit(-1);                                                                                                          \
  }

#define _CHK_(size1, size2)                                                                                            \
  if ((size1) != (size2)) {                                                                                            \
    printf("Error! rank%d, %s:%d, %s\n", myid, __FILE__, __LINE__, __FUNCTION__);                                     \
    fflush(stdout);                                                                                                    \
    exit(-1);                                                                                                          \
  }
#define _NL_ status = func(&endl, sizeof(char), 1, file_p);
#define _FWD_(type, var)                                                                                               \
  status = func(&var, sizeof(type), 1, file_p);                                                                        \
  if (print) {                                                                                                         \
    printf("rank%d, %s = %d\n", myid, #var, var);                                                                      \
  }
#define _FW_(type, var, size, var_out, cl_mem_flag)                                                                    \
  cl_buf_com.var_out = clCreateBuffer(context, cl_mem_flag, sizeof(type) * (size), var, &error);                       \
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clCreateBuffer failed");

#define _FWV_(type, var, size, var_out, cl_mem_flag)                                                                    \
  var_out = clCreateBuffer(context, cl_mem_flag, sizeof(type) * (size), var, &error);                       \
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clCreateBuffer failed");

#define clSetKernelArgWithCheck(kernel, arg_index, arg_size, arg_value) {\
  cl_int error = clSetKernelArg(kernel, arg_index, arg_size, arg_value); \
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed");  }

void loadProgramSource(const char **files, size_t length, char **buffer, size_t *sizes) {
  /* Read each source file (*.cl) and store the contents into a temporary datastore */
  for (size_t i = 0; i < length; i++) {
    FILE *file = fopen(files[i], "r");
    if (file == NULL) {
      printf("Couldn't read the program file\n");
      fflush(stdout);
      exit(1);
    }
    fseek(file, 0, SEEK_END);
    sizes[i] = ftell(file);
    rewind(file); // reset the file pointer so that 'fread' reads from the front
    buffer[i] = (char *)malloc(sizes[i] + 1);
    buffer[i][sizes[i]] = '\0';
    size_t status = fread(buffer[i], sizeof(char), sizes[i], file);
    IF_ERROR_EXIT(status != sizes[i], (int)status, "fread filed, read num != file size");
    fclose(file);
  }
}

void opencl_c_general_init(){

}

void displayDeviceDetails(cl_device_id id, cl_device_info param_name, const char* paramNameAsStr) ; 

void displayPlatformInfo(cl_platform_id id,
                         cl_platform_info param_name,
                         const char* paramNameAsStr) {
    cl_int error = 0;
    size_t paramSize = 0;
    error = clGetPlatformInfo( id, param_name, 0, NULL, &paramSize );
    char* moreInfo = (char*)alloca( sizeof(char) * paramSize);
    error = clGetPlatformInfo( id, param_name, paramSize, moreInfo, NULL );
    if (error != CL_SUCCESS ) {
        perror("Unable to find any OpenCL platform information");
        return;
    }
    printf("%s: %s\n", paramNameAsStr, moreInfo);
}

void displayDeviceInfo(cl_platform_id id, 
                       cl_device_type dev_type) {
    /* OpenCL 1.1 device types */
    cl_int error = 0;
    cl_uint numOfDevices = 0;

    /* Determine how many devices are connected to your platform */
    error = clGetDeviceIDs(id, dev_type, 0, NULL, &numOfDevices);
    if (error != CL_SUCCESS ) { 
        perror("Unable to obtain any OpenCL compliant device info");
        exit(1);
    }
    cl_device_id* devices = (cl_device_id*) alloca(sizeof(cl_device_id) * numOfDevices);

    /* Load the information about your devices into the variable 'devices' */
    error = clGetDeviceIDs(id, dev_type, numOfDevices, devices, NULL);
    if (error != CL_SUCCESS ) { 
        perror("Unable to obtain any OpenCL compliant device info");
        exit(1);
    }
    printf("Number of detected OpenCL devices: %d\n", numOfDevices);
    /* We attempt to retrieve some information about the devices. */
    for(int i = 0; i < numOfDevices; ++ i ) {
        displayDeviceDetails( devices[i], CL_DEVICE_TYPE, "CL_DEVICE_TYPE" );
        displayDeviceDetails( devices[i], CL_DEVICE_NAME, "CL_DEVICE_NAME" );
        displayDeviceDetails( devices[i], CL_DEVICE_VENDOR, "CL_DEVICE_VENDOR" );
        displayDeviceDetails( devices[i], CL_DEVICE_VENDOR_ID, "CL_DEVICE_VENDOR_ID" );
        displayDeviceDetails( devices[i], CL_DEVICE_MAX_MEM_ALLOC_SIZE, "CL_DEVICE_MAX_MEM_ALLOC_SIZE" );
        displayDeviceDetails( devices[i], CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, "CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE" );
        displayDeviceDetails( devices[i], CL_DEVICE_GLOBAL_MEM_SIZE, "CL_DEVICE_GLOBAL_MEM_SIZE" );
        displayDeviceDetails( devices[i], CL_DEVICE_MAX_COMPUTE_UNITS, "CL_DEVICE_MAX_COMPUTE_UNITS" );
        displayDeviceDetails( devices[i], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS" );
        displayDeviceDetails( devices[i], CL_DEVICE_MAX_WORK_ITEM_SIZES, "CL_DEVICE_MAX_WORK_ITEM_SIZES" );
        displayDeviceDetails( devices[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, "CL_DEVICE_MAX_WORK_GROUP_SIZE" );
        displayDeviceDetails( devices[i], CL_DEVICE_TOPOLOGY_AMD, "CL_DEVICE_TOPOLOGY_AMD" );
    }
}

void displayDeviceDetails(cl_device_id id,
                          cl_device_info param_name, 
                          const char* paramNameAsStr) {
  cl_int error = 0;
  size_t paramSize = 0;

  error = clGetDeviceInfo( id, param_name, 0, NULL, &paramSize );
  if (error != CL_SUCCESS ) {
    perror("Unable to obtain device info for param\n");
    return;
  }

  /* the cl_device_info are preprocessor directives defined in cl.h */
  switch (param_name) {
    case CL_DEVICE_TYPE: {
            cl_device_type* devType = (cl_device_type*) alloca(sizeof(cl_device_type) * paramSize);
            error = clGetDeviceInfo( id, param_name, paramSize, devType, NULL );
            if (error != CL_SUCCESS ) {
                perror("Unable to obtain device info for param\n");
                return;
            }
            switch (*devType) {
              case CL_DEVICE_TYPE_CPU : printf("CPU detected\n");break;
              case CL_DEVICE_TYPE_GPU : printf("GPU detected\n");break;
              case CL_DEVICE_TYPE_ACCELERATOR : printf("Accelerator detected\n");break;
              case CL_DEVICE_TYPE_DEFAULT : printf("default detected\n");break;
            }
            }break;
    case CL_DEVICE_VENDOR_ID : 
    case CL_DEVICE_MAX_COMPUTE_UNITS : 
    case CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS : {
            cl_uint* ret = (cl_uint*) alloca(sizeof(cl_uint) * paramSize);
            error = clGetDeviceInfo( id, param_name, paramSize, ret, NULL );
            if (error != CL_SUCCESS ) {
                perror("Unable to obtain device info for param\n");
                return;
            }
            switch (param_name) {
                case CL_DEVICE_VENDOR_ID: printf("\tVENDOR ID: 0x%x\n", *ret); break;
                case CL_DEVICE_MAX_COMPUTE_UNITS: printf("\tMaximum number of parallel compute units: %d\n", *ret); break;
                case CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS: printf("\tMaximum dimensions for global/local work-item IDs: %d\n", *ret); break;
            }
         }break;
    case CL_DEVICE_MAX_WORK_ITEM_SIZES : {
            cl_uint maxWIDimensions;
            size_t* ret = (size_t*) alloca(sizeof(size_t) * paramSize);
            error = clGetDeviceInfo( id, param_name, paramSize, ret, NULL );

            error = clGetDeviceInfo( id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &maxWIDimensions, NULL );
            if (error != CL_SUCCESS ) {
                perror("Unable to obtain device info for param\n");
                return;
            }
            printf("\tMaximum number of work-items in each dimension: ( ");
            for(cl_int i =0; i < maxWIDimensions; ++i ) {
                printf("%lu ", ret[i]);
            }
            printf(" )\n");
            }break;
    case CL_DEVICE_MAX_WORK_GROUP_SIZE : {
            size_t* ret = (size_t*) alloca(sizeof(size_t) * paramSize);
            error = clGetDeviceInfo( id, param_name, paramSize, ret, NULL );
            if (error != CL_SUCCESS ) {
                perror("Unable to obtain device info for param\n");
                return;
            }
            printf("\tMaximum number of work-items in a work-group: %lu\n", *ret);
            }break;
    case CL_DEVICE_NAME :
    case CL_DEVICE_VENDOR : {
            char data[48];
            error = clGetDeviceInfo( id, param_name, paramSize, data, NULL );
            if (error != CL_SUCCESS ) {
                perror("Unable to obtain device name/vendor info for param\n");
                return;
            }
            switch (param_name) {
                case CL_DEVICE_NAME : printf("\tDevice name is %s\n", data);break;
                case CL_DEVICE_VENDOR : printf("\tDevice vendor is %s\n", data);break;
            }
    } break;
    case CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE: {
            cl_uint* size = (cl_uint*) alloca(sizeof(cl_uint) * paramSize);
            error = clGetDeviceInfo( id, param_name, paramSize, size, NULL );
            if (error != CL_SUCCESS ) {
                perror("Unable to obtain device name/vendor info for param\n");
                return;
            }
            printf("\tDevice global cacheline size: %d bytes\n", (*size)); break;
    } break;
    case CL_DEVICE_GLOBAL_MEM_SIZE:
    case CL_DEVICE_MAX_MEM_ALLOC_SIZE: {
            cl_ulong* size = (cl_ulong*) alloca(sizeof(cl_ulong) * paramSize);
            error = clGetDeviceInfo( id, param_name, paramSize, size, NULL );
            if (error != CL_SUCCESS ) {
                perror("Unable to obtain device name/vendor info for param\n");
                return;
            }
            switch (param_name) {
                case CL_DEVICE_GLOBAL_MEM_SIZE: printf("\tDevice global mem: %ld mega-bytes\n", (*size)>>20); break;
                case CL_DEVICE_MAX_MEM_ALLOC_SIZE: printf("\tDevice max memory allocation: %ld mega-bytes\n", (*size)>>20); break; 
            }
    } break;
    case CL_DEVICE_TOPOLOGY_AMD: {
            cl_device_topology_amd topology;
            error = clGetDeviceInfo (id, param_name, sizeof(cl_device_topology_amd), &topology, NULL);
            if (error != CL_SUCCESS ) {
                perror("Unable to obtain device cl_device_topology_amd info for param\n");
                return;
            }
            if (topology.raw.type == CL_DEVICE_TOPOLOGY_TYPE_PCIE_AMD) {
                printf("INFO: Topology: PCI[ B#%d, D#%d, F#%d ]\n", (int)topology.pcie.bus
                    , (int)topology.pcie.device, (int)topology.pcie.function);
            }
    } break;

  } //end of switch
         
}

void opencl_init_() {
  if(opencl_init_finished)
    return;
  opencl_init_finished = 1;
  cl_int error;

  if(myid == 0){
    printf("opencl_util_debug_io = %d\n", opencl_util_debug_io);
  }

  // 得到平台数量
  error = clGetPlatformIDs(0, NULL, &numOfPlatforms);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "Unable to find any OpenCL platforms");
  platforms = (cl_platform_id *)alloca(sizeof(cl_platform_id) * numOfPlatforms);
  // printf("Number of OpenCL platforms found: %d\n", numOfPlatforms);
  IF_ERROR_EXIT(numOfPlatforms <= platformId, numOfPlatforms, "The selected platformId is out of range");
  error = clGetPlatformIDs(numOfPlatforms, platforms, NULL);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "Unable to find any OpenCL platforms");
  if(myid == 0){
    printf("===================================================\n");
    for(cl_uint i = 0; i < numOfPlatforms; ++i) {
      displayPlatformInfo( platforms[i], CL_PLATFORM_PROFILE, "CL_PLATFORM_PROFILE" );
      displayPlatformInfo( platforms[i], CL_PLATFORM_VERSION, "CL_PLATFORM_VERSION" );
      displayPlatformInfo( platforms[i], CL_PLATFORM_NAME,    "CL_PLATFORM_NAME" );
      displayPlatformInfo( platforms[i], CL_PLATFORM_VENDOR,  "CL_PLATFORM_VENDOR" );
      displayPlatformInfo( platforms[i], CL_PLATFORM_EXTENSIONS, "CL_PLATFORM_EXTENSIONS" );
      // Assume that we don't know how many devices are OpenCL compliant, we locate everything !
      displayDeviceInfo( platforms[i], CL_DEVICE_TYPE_ALL );
      printf("--------------------------------\n");     
    }
    printf("===================================================\n");
  }

  // int device_id = DEVICE_ID;
  // 选择指定平台上的指定设备
  // int device_id = 0;
  int device_id = mpi_platform_relative_id / mpi_task_per_gpu;

  if(myid <= mpi_per_node){
    printf("rank%d, will choose platform %d, device %d\n", myid, (int)platformId, (int)device_id);
  }

  // 得到选定平台设备数量
  cl_uint numOfDevices = 0;
  error = clGetDeviceIDs(platforms[platformId], CL_DEVICE_TYPE_ALL, 0, NULL, &numOfDevices);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "Unable to obtain any OpenCL compliant device info");
  cl_device_id *devices = (cl_device_id *)alloca(sizeof(cl_device_id) * numOfDevices);
  IF_ERROR_EXIT(numOfDevices <= device_id, numOfDevices, "numOfDevices <= device_id, could not choose");
  // 加载选定平台设备信息
  error = clGetDeviceIDs(platforms[platformId], CL_DEVICE_TYPE_ALL, numOfDevices, devices, NULL);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "Unable to obtain any OpenCL compliant device info");
  // printf("myid=%d, Number of detected OpenCL devices: %d, choose %d\n", myid, numOfDevices, device_id);
  // fflush(stdout);
  // cl_device_topology_amd topology;
  // error = clGetDeviceInfo (devices[device_id], CL_DEVICE_TOPOLOGY_AMD, sizeof(cl_device_topology_amd), &topology, NULL);
  // IF_ERROR_EXIT(error != CL_SUCCESS, error, "Unable to obtain device cl_device_topology_amd info for param");
  // printf("myid=%d, Number of detected OpenCL devices: %d, choose %d, Topology: PCI[ B#%d, D#%d, F#%d ]\n", 
  //   myid, numOfDevices, device_id, (int)topology.pcie.bus , (int)topology.pcie.device, (int)topology.pcie.function);
  // 创建 OpenCL 上下文
  cl_context_properties ctx[3] = {CL_CONTEXT_PLATFORM, (cl_context_properties)platforms[platformId], 0};
  context = clCreateContext(ctx, numOfDevices, devices, NULL, NULL, &error);
  // context = clCreateContext(ctx, 1, &devices[device_id], NULL, NULL, &error);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "Can't create a valid OpenCL context");
  // 读入一个或多个文件，仅文件操作，OpenCL 无关
  loadProgramSource(file_names, number_of_files, buffer, sizes);

  // 创建程序对象
  program = clCreateProgramWithSource(context, number_of_files, (const char **)buffer, sizes, &error);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "Can't create the OpenCL program object");
  // 编译程序对象
  char *program_log;
  size_t log_size;
  char clBuildOption[512];
  int N_PERIODIC_OR_USE_HARTREE_NON_PERIODIC_EWALD;
  if(n_periodic > 0 || use_hartree_non_periodic_ewald){
    N_PERIODIC_OR_USE_HARTREE_NON_PERIODIC_EWALD = 1;
  }else{
    N_PERIODIC_OR_USE_HARTREE_NON_PERIODIC_EWALD = 0;
  }
	sprintf(clBuildOption, "-DL_POT_MAX=%d -DHARTREE_FP_FUNCTION_SPLINES%d -DN_PERIODIC_OR_USE_HARTREE_NON_PERIODIC_EWALD%d"
    " -DMACRO_n_centers_integrals=%d -DMACRO_n_max_compute_atoms=%d -DMACRO_n_centers=%d"
    " -DLOCALSIZE_SUM_UP_PRE_PROC=%d", 
    l_pot_max, hartree_fp_function_splines, N_PERIODIC_OR_USE_HARTREE_NON_PERIODIC_EWALD, n_centers_integrals,
    n_max_compute_atoms, n_centers, (int)localSize_sum_up_pre_proc[0]);
  if(myid == 0){
	  printf("clBuildOption=\"%s\"\n", clBuildOption);
    printf("n_periodic = %d, use_hartree_non_periodic_ewald = %d\n", n_periodic, use_hartree_non_periodic_ewald);
  }
  error = clBuildProgram(program, 1, &devices[device_id], clBuildOption, NULL, NULL);
  if (error != CL_SUCCESS) {
    // 编译失败时打印错误
    clGetProgramBuildInfo(program, devices[device_id], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
    program_log = (char *)malloc(log_size + 1);
    program_log[log_size] = '\0';
    clGetProgramBuildInfo(program, devices[device_id], CL_PROGRAM_BUILD_LOG, log_size + 1, program_log, NULL);
    printf("\n=== clBuildProgram ERROR ===\n\n%s\n=============\n", program_log);
    fflush(stdout);
    free(program_log);
    exit(1);
  }
  // 构建 kernel
  for (int i = 0; i < number_of_kernels; i++) {
    kernels[i] = clCreateKernel(program, kernel_names[i], &error);
    char tmp_str[256];
    sprintf(tmp_str, "clCreateKernel failed, kernel name: %s", kernel_names[i]);
    IF_ERROR_EXIT(kernels[i] == NULL || error != CL_SUCCESS, error, tmp_str);
  }
  // 构建命令队列
  cQ = clCreateCommandQueueWithProperties(context, devices[device_id], 0, &error);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "Unable to create command-queue");
  if(myid == 0){
    printf("rank%d, %s finished\n", myid, __FUNCTION__);
    fflush(stdout);
  }
}

void opencl_finish() {
  if(!opencl_init_finished)
    return;
  opencl_init_finished = 0;
  clReleaseCommandQueue(cQ);
  for (int i = 0; i < number_of_kernels; i++)
    clReleaseKernel(kernels[i]);
  for (int i = 0; i < number_of_files; i++)
    free(buffer[i]);
  clReleaseProgram(program);
  clReleaseContext(context);
}

void opencl_common_buffer_init_() {
  if(opencl_common_buffer_init_finished)
    return;
  opencl_common_buffer_init_finished = 1;
  cl_int error;
  unsigned int arg_index = 0;
  _FW_(int, MV(geometry, species), n_atoms, species, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(pbc_lists, centers_hartree_potential), n_centers_hartree_potential, centers_hartree_potential,
       CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(pbc_lists, center_to_atom), n_centers, center_to_atom, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(pbc_lists, species_center), n_centers, species_center, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, MV(pbc_lists, coords_center), 3 * n_centers, coords_center, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(species_data, l_hartree), n_species, l_hartree, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(grids, n_grid), n_species, n_grid, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(grids, n_radial), n_species, n_radial, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);

  _FW_(double, MV(grids, r_grid_min), n_species, r_grid_min, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, MV(grids, r_grid_inc), n_species, r_grid_inc, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, MV(grids, log_r_grid_inc), n_species, log_r_grid_inc, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, MV(grids, scale_radial), n_species, scale_radial, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(analytic_multipole_coefficients, n_cc_lm_ijk), (l_max_analytic_multipole + 1), n_cc_lm_ijk,
       CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(analytic_multipole_coefficients, index_cc), n_cc_lm_ijk(l_max_analytic_multipole) * 6, index_cc,
       CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(analytic_multipole_coefficients, index_ijk_max_cc), 3 * (l_max_analytic_multipole + 1), index_ijk_max_cc,
       CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, MV(hartree_potential_real_p0, b0), pmaxab + 1, b0, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, MV(hartree_potential_real_p0, b2), pmaxab + 1, b2, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, MV(hartree_potential_real_p0, b4), pmaxab + 1, b4, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, MV(hartree_potential_real_p0, b6), pmaxab + 1, b6, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, MV(hartree_potential_real_p0, a_save), pmaxab + 1, a_save, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  // _FW_(double, Fp_function_spline_slice, (lmax_Fp + 1) * 4 * (Fp_max_grid+1), Fp_function_spline_slice,
  //      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  // _FW_(double, Fpc_function_spline_slice, (lmax_Fp + 1) * 4 * (Fp_max_grid+1), Fpc_function_spline_slice,
  //      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  // rho global
  _FW_(int, MV(basis, perm_basis_fns_spl), n_basis_fns, perm_basis_fns_spl, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, MV(basis, outer_radius_sq), n_basis_fns, outer_radius_sq, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(basis, basis_fn), n_basis, basis_fn, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(basis, basis_l), n_basis, basis_l, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, MV(basis, atom_radius_sq), n_species, atom_radius_sq, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(basis, basis_fn_start_spl), n_species, basis_fn_start_spl, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(basis, basis_fn_atom), n_basis_fns * n_atoms, basis_fn_atom, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, MV(basis, basis_wave_ordered), n_basis_fns * n_max_spline * n_max_grid, basis_wave_ordered, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);  // 进程间可能不同
  _FW_(double, MV(basis, basis_kinetic_ordered), n_basis_fns * n_max_spline * n_max_grid, basis_kinetic_ordered, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);  // 进程间可能不同

  _FW_(int, MV(pbc_lists, cbasis_to_basis), n_centers_basis_T, Cbasis_to_basis, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(pbc_lists, cbasis_to_center), n_centers_basis_T, Cbasis_to_center, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(pbc_lists, centers_basis_integrals), n_centers_basis_integrals, centers_basis_integrals, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, MV(pbc_lists, index_hamiltonian), 2 * index_hamiltonian_dim2 * n_basis, index_hamiltonian, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR); // 进程间可能不同
  _FW_(int, MV(pbc_lists, position_in_hamiltonian), position_in_hamiltonian_dim1 * position_in_hamiltonian_dim2, position_in_hamiltonian, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);  // 进程间可能不同
  _FW_(int, MV(pbc_lists, column_index_hamiltonian), column_index_hamiltonian_size, column_index_hamiltonian, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);  // 进程间可能不同

  _FW_(int, MV(pbc_lists, center_to_cell), n_centers, center_to_cell, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);

  // if(ctrl_use_sumup_pre_c_cl_version){
  _FW_(double, MV(grids, r_radial), n_max_radial * n_species, r_radial, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, MV(grids, r_grid), n_max_grid * n_species, r_grid, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, NULL, n_atoms, rho_multipole_index, CL_MEM_READ_ONLY);
  if(compensate_multipole_errors){
    _FW_(int, MV(hartree_potential_storage, compensation_norm), n_atoms, compensation_norm, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
    _FW_(int, MV(hartree_potential_storage, compensation_radius), n_atoms, compensation_radius, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  } else {
    _FW_(int, NULL, 1, compensation_norm, CL_MEM_READ_ONLY);
    _FW_(int, NULL, 1, compensation_radius, CL_MEM_READ_ONLY);
  }

  _FW_(double, MV(species_data, multipole_radius_free), n_species, multipole_radius_free, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  // }
  // ------
}

void opencl_common_buffer_free_() {
  if(!opencl_common_buffer_init_finished)
    return;
  opencl_common_buffer_init_finished = 0;
  cl_int error;
  unsigned int arg_index = 0;

  clReleaseMemObject(cl_buf_com.species);
  clReleaseMemObject(cl_buf_com.centers_hartree_potential);
  clReleaseMemObject(cl_buf_com.center_to_atom);
  clReleaseMemObject(cl_buf_com.species_center);
  clReleaseMemObject(cl_buf_com.coords_center);
  clReleaseMemObject(cl_buf_com.l_hartree);
  clReleaseMemObject(cl_buf_com.n_grid);
  clReleaseMemObject(cl_buf_com.n_radial);

  clReleaseMemObject(cl_buf_com.r_grid_min);
  clReleaseMemObject(cl_buf_com.r_grid_inc);
  clReleaseMemObject(cl_buf_com.log_r_grid_inc);
  clReleaseMemObject(cl_buf_com.scale_radial);
  clReleaseMemObject(cl_buf_com.n_cc_lm_ijk);
  clReleaseMemObject(cl_buf_com.index_cc);
  clReleaseMemObject(cl_buf_com.index_ijk_max_cc);
  clReleaseMemObject(cl_buf_com.b0);
  clReleaseMemObject(cl_buf_com.b2);
  clReleaseMemObject(cl_buf_com.b4);
  clReleaseMemObject(cl_buf_com.b6);
  clReleaseMemObject(cl_buf_com.a_save);

  // rho global
  clReleaseMemObject(cl_buf_com.perm_basis_fns_spl);
  clReleaseMemObject(cl_buf_com.outer_radius_sq);
  clReleaseMemObject(cl_buf_com.basis_fn);
  clReleaseMemObject(cl_buf_com.basis_l);
  clReleaseMemObject(cl_buf_com.atom_radius_sq);
  clReleaseMemObject(cl_buf_com.basis_fn_start_spl);
  clReleaseMemObject(cl_buf_com.basis_fn_atom);
  clReleaseMemObject(cl_buf_com.basis_wave_ordered);
  clReleaseMemObject(cl_buf_com.basis_kinetic_ordered);

  clReleaseMemObject(cl_buf_com.Cbasis_to_basis);
  clReleaseMemObject(cl_buf_com.Cbasis_to_center);
  clReleaseMemObject(cl_buf_com.centers_basis_integrals);
  clReleaseMemObject(cl_buf_com.index_hamiltonian);
  clReleaseMemObject(cl_buf_com.position_in_hamiltonian);
  clReleaseMemObject(cl_buf_com.column_index_hamiltonian);

  clReleaseMemObject(cl_buf_com.center_to_cell);

  // if(ctrl_use_sumup_pre_c_cl_version){
  clReleaseMemObject(cl_buf_com.r_radial);
  clReleaseMemObject(cl_buf_com.r_grid);
  clReleaseMemObject(cl_buf_com.rho_multipole_index);
  clReleaseMemObject(cl_buf_com.compensation_norm);
  clReleaseMemObject(cl_buf_com.compensation_radius);
  // clReleaseMemObject(cl_buf_com.rho_multipole_h_p_s);
  clReleaseMemObject(cl_buf_com.multipole_radius_free);
  // }
  // ------
}

void sum_up_pre_processing_init_(){
  opencl_init_();
  opencl_common_buffer_init_();

  globalSize_sum_up_pre_proc[0] = 256 * localSize_sum_up_pre_proc[0];  // static !
  // globalSize_sum_up_pre_proc[0] = n_atoms * localSize_sum_up_pre_proc[0];  // dynamic !

  int arg_index = 0;
  cl_int error;
  
  _FW_(double, NULL, (l_pot_max + 1) * (l_pot_max + 1) * n_max_grid * (globalSize_sum_up_pre_proc[0] / localSize_sum_up_pre_proc[0]), 
      angular_integral_log, CL_MEM_READ_WRITE);

  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(int), &n_max_radial);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(int), &l_pot_max);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(int), &n_max_spline);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(int), &n_hartree_grid);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(int), &n_atoms);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(int), &n_max_grid);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(int), &Adams_Moulton_integrator);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(int), &compensate_multipole_errors);

  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.species);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.l_hartree);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.multipole_radius_free);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.n_grid);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.n_radial);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.r_grid_inc);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.scale_radial);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.r_grid);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.r_radial);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.rho_multipole_h_p_s);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.rho_multipole_index);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.compensation_norm);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.compensation_radius);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.centers_hartree_potential);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.center_to_atom);

  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), &cl_buf_com.angular_integral_log);
}

void sum_up_pre_processing_part_(int* n_coeff_hartree_, int* i_center_begin_, int* i_center_end_, 
  cl_mem* centers_rho_multipole_spl, cl_mem* centers_delta_v_hart_part_spl, cl_mem* i_center_to_centers_index, 
  int debug){
  struct timeval start, end;
  gettimeofday(&start, NULL);
  // cl_buf_com.centers_rho_multipole_spl;
  // _FW_(double, NULL, (l_pot_max + 1) * (l_pot_max + 1) * n_max_spline * (n_max_radial + 2) * centers_tile_size, centers_rho_multipole_spl, CL_MEM_READ_WRITE);
  // _FW_(double, NULL, (l_pot_max + 1) * (l_pot_max + 1) * n_coeff_hartree * n_hartree_grid * centers_tile_size, centers_delta_v_hart_part_spl, CL_MEM_READ_WRITE);

  int arg_index = 24;
  cl_int error;
  
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), centers_rho_multipole_spl);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), centers_delta_v_hart_part_spl);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(int), n_coeff_hartree_);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(int), i_center_begin_);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(int), i_center_end_);
  clSetKernelArgWithCheck(kernels[3], arg_index++, sizeof(cl_mem), i_center_to_centers_index);

  error = clEnqueueNDRangeKernel(cQ, kernels[3], 1, NULL, globalSize_sum_up_pre_proc, localSize_sum_up_pre_proc, 0, NULL, NULL);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clEnqueueNDRangeKernel failed")

  // clEnqueueReadBuffer(cQ, centers_rho_multipole_spl, CL_TRUE, 0, 
  //   sizeof(double) * (l_pot_max + 1) * (l_pot_max + 1) * n_max_spline * (n_max_radial + 2) * (*i_center_end_ - *i_center_begin_),
  //   param_centers_rho_multipole_spl, 1, &event, NULL);
  // clEnqueueReadBuffer(cQ, centers_delta_v_hart_part_spl, CL_TRUE, 0, 
  //   sizeof(double) * (l_pot_max + 1) * (l_pot_max + 1) * n_coeff_hartree * n_hartree_grid * (*i_center_end_ - *i_center_begin_),
  //   param_centers_delta_v_hart_part_spl, 1, &event, NULL);
  clFinish(cQ);

  gettimeofday(&end, NULL);
  long time_use = 1000000 * ( end.tv_sec - start.tv_sec ) + end.tv_usec - start.tv_usec;
  // if(myid == 0 && debug)
  //   printf("    rank%d, %s: %lf seconds\n", myid, __func__, time_use/1000000.0);
}

void sum_up_pre_processing_finish(){
  clReleaseMemObject(cl_buf_com.angular_integral_log);
}

void sum_up_first_begin(){
  if(sum_up_first_begin_finished)
    return;
  sum_up_first_begin_finished = 1;

  opencl_init_();
  opencl_common_buffer_init_();

  int arg_index;
  cl_int error;
  int fast_ylm = 1;
  int new_ylm = 0;

  // TODO 这两没有释放
  // TODO WARNING 进程间协作可能没有处理这个
  if(Fp_max_grid == 0){
    _FW_(double, NULL, 1, Fp_function_spline_slice, CL_MEM_READ_ONLY);
    _FW_(double, NULL, 1, Fpc_function_spline_slice, CL_MEM_READ_ONLY);
  } else {
    _FW_(double, MV(hartree_f_p_functions, fp_function_spline), (lmax_Fp + 1) * n_max_spline * (Fp_max_grid), Fp_function_spline_slice,
         CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
    _FW_(double, MV(hartree_f_p_functions, fpc_function_spline), (lmax_Fp + 1) * n_max_spline * (Fp_max_grid), Fpc_function_spline_slice,
         CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  }

  arg_index = 11;
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &n_centers_hartree_potential);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &n_periodic);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &n_max_radial);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &l_pot_max);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &n_max_spline);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &n_hartree_grid);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &n_species);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &n_atoms);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &n_centers);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &n_max_batch_size);
  remember_arg_n_my_batches_work_sumup = arg_index;
  arg_index += 2;
  // =~error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &n_my_batches_work_sumup);
  // =~error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &n_full_points_work_sumup);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &use_hartree_non_periodic_ewald);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &hartree_fp_function_splines);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &fast_ylm);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &new_ylm);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &l_max_analytic_multipole);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &n_hartree_atoms);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &hartree_force_l_add);
  remember_arg_Fp_max_grid = arg_index;
  arg_index += 5;
  // =~error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &Fp_max_grid);
  // =~error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &lmax_Fp);
  // =~error = clSetKernelArg(kernels[0], arg_index++, sizeof(double), &Fp_grid_min);
  // =~error = clSetKernelArg(kernels[0], arg_index++, sizeof(double), &Fp_grid_inc);
  // =~error = clSetKernelArg(kernels[0], arg_index++, sizeof(double), &Fp_grid_max);
  arg_index = 35;
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.species);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.centers_hartree_potential);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.center_to_atom);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.species_center);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.coords_center);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.l_hartree);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.n_grid);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.n_radial);
  remember_arg_index_1 = arg_index;
  // error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.batches_size_s);
  // error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.batches_points_coords_s);
  arg_index += 2;
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.r_grid_min);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.log_r_grid_inc);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.scale_radial);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.n_cc_lm_ijk);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.index_cc);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.index_ijk_max_cc);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.b0);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.b2);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.b4);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.b6);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.a_save);
  remember_arg_Fp_function_spline = arg_index;

}

void sum_up_begin_0_() {
  if(sum_up_begin_0_finished)
    return;
  sum_up_begin_0_finished = 1;

  sum_up_first_begin();
  cl_int error;

  // _FW_(double, sum_up_param.centers_rho_multipole_spl, (l_pot_max + 1) * (l_pot_max + 1) * n_max_spline * (n_max_radial + 2) * n_atoms,
  //      centers_rho_multipole_spl, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  // _FW_(double, sum_up_param.centers_delta_v_hart_part_spl, (l_pot_max + 1) * (l_pot_max + 1) * n_coeff_hartree * n_hartree_grid * n_atoms,
  //      centers_delta_v_hart_part_spl, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);  // CL_MEM_USE_HOST_PTR

}

void sum_up_begin() {
  sum_up_begin_0_();
  sum_up_begin_0_finished = 0;

  long time_uses[32];
  char *time_infos[32];
  size_t time_index = 0;
  struct timeval start, end;
  gettimeofday(&start, NULL);

  cl_int error;
  unsigned int arg_index = 0;
  // unsigned int arg_index_centers_rho_multipole_spl = 0;
  unsigned int arg_index_i_center_begin = 0;

  // int i_center_tile_size = n_centers_hartree_potential;
  int i_center_tile_size = MIN(i_center_tile_size_default, n_centers_hartree_potential);

  size_t localSize[] = {256};
  size_t globalSize[] = {256 * 128};

  // printf("index_cc size %d\n", n_cc_lm_ijk(l_max_analytic_multipole));
  int *index_cc_aos = (int *)malloc(sizeof(int) * n_cc_lm_ijk(l_max_analytic_multipole) * 6);
  for(int i=1; i<=n_cc_lm_ijk(l_max_analytic_multipole); i++){
    index_cc_aos[(i-1)*4] = index_cc(i, 3, n_cc_lm_ijk(l_max_analytic_multipole));
    index_cc_aos[(i-1)*4+1] = index_cc(i, 4, n_cc_lm_ijk(l_max_analytic_multipole));
    index_cc_aos[(i-1)*4+2] = index_cc(i, 5, n_cc_lm_ijk(l_max_analytic_multipole));
    index_cc_aos[(i-1)*4+3] = index_cc(i, 6, n_cc_lm_ijk(l_max_analytic_multipole));
  }
  _FW_(int, index_cc_aos, n_cc_lm_ijk(l_max_analytic_multipole) * 6, index_cc_aos, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);

  // TODO 这一块应该也可以塞到 sum_up_begin_0_ 中，256*128*11*11x8 也才几十兆，
  // sumup tmp
  _FW_(double, NULL, globalSize[0] * (l_pot_max + 2), Fp, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * 3 * (l_pot_max + 1), coord_c, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * 1, coord_mat, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * 1, rest_mat, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * 1, vector, CL_MEM_READ_WRITE);  // 留下这个，不过这个理应已经废弃
  _FW_(double, NULL, globalSize[0] * 1, delta_v_hartree_multipole_component,
       CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * 1, rho_multipole_component, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * 1, ylm_tab, CL_MEM_READ_WRITE);
  // sumup tmp
  arg_index = 63;
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.Fp);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.coord_c);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.coord_mat);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.rest_mat);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.vector);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.delta_v_hartree_multipole_component);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.rho_multipole_component);
  error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.ylm_tab);
  arg_index_i_center_begin = arg_index;

  // TODO 这一块也可以放前面去
  int max_spl_atom = 1;
  for(int i_center = 1; i_center <= n_centers_hartree_potential; i_center++){
      int current_center = MV(pbc_lists, centers_hartree_potential)[i_center - 1];
      int current_spl_atom = MV(pbc_lists, center_to_atom)[current_center - 1];
      max_spl_atom = max_spl_atom > current_spl_atom ? max_spl_atom : current_spl_atom;
  }

  int *spl_atom_to_i_center = (int*)malloc(sizeof(int) * (max_spl_atom+1));
  int *i_center_to_centers_index = (int*)malloc(sizeof(int) * n_centers_hartree_potential);

  for(int i_center_tile = 0; i_center_tile < n_centers_hartree_potential; i_center_tile += i_center_tile_size){

    for(int i=0; i<(max_spl_atom+1); i++)
      spl_atom_to_i_center[i] = -1;

    for(int i_center_ = i_center_tile; i_center_ < MIN(i_center_tile + i_center_tile_size, n_centers_hartree_potential); i_center_++){
      int i_center = i_center_ + 1;
      int current_center = MV(pbc_lists, centers_hartree_potential)[i_center_];
      int current_spl_atom = MV(pbc_lists, center_to_atom)[current_center - 1];
      
      
      if(spl_atom_to_i_center[current_spl_atom] != -1){
        i_center_to_centers_index[i_center_] = spl_atom_to_i_center[current_spl_atom];
      } else {
        i_center_to_centers_index[i_center_] = i_center_ - i_center_tile;
        spl_atom_to_i_center[current_spl_atom] = i_center_ - i_center_tile;
      }
    }
  }

  _FW_(int, i_center_to_centers_index, n_centers_hartree_potential, i_center_to_centers_index, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  error = clSetKernelArg(kernels[0], arg_index_i_center_begin+2, sizeof(cl_mem), &cl_buf_com.i_center_to_centers_index);

  error = clEnqueueWriteBuffer(cQ, cl_buf_com.rho_multipole_index, CL_TRUE, 0, sizeof(int) * n_atoms, 
    MV(hartree_potential_storage, rho_multipole_index), 0, NULL, NULL);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clEnqueueWriteBuffer failed");

  error = clSetKernelArg(kernels[0], remember_arg_Fp_function_spline, sizeof(cl_mem), &cl_buf_com.Fp_function_spline_slice);  // TODO 赋值位置可能得改
  error = clSetKernelArg(kernels[0], remember_arg_Fp_function_spline+1, sizeof(cl_mem), &cl_buf_com.Fpc_function_spline_slice); // TODO 赋值位置可能得改

  _FW_(double, NULL, (l_pot_max + 1) * (l_pot_max + 1) * n_max_spline * (n_max_radial + 2) * i_center_tile_size, centers_rho_multipole_spl, CL_MEM_READ_WRITE);
  _FW_(double, NULL, (l_pot_max + 1) * (l_pot_max + 1) * n_coeff_hartree * n_hartree_grid * i_center_tile_size, centers_delta_v_hart_part_spl, CL_MEM_READ_WRITE);
  // _FW_(double, sum_up_param.centers_rho_multipole_spl, (l_pot_max + 1) * (l_pot_max + 1) * n_max_spline * (n_max_radial + 2) * n_atoms,
  //      centers_rho_multipole_spl, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  // _FW_(double, sum_up_param.centers_delta_v_hart_part_spl, (l_pot_max + 1) * (l_pot_max + 1) * n_coeff_hartree * n_hartree_grid * n_atoms,
  //      centers_delta_v_hart_part_spl, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);  // CL_MEM_USE_HOST_PTR

  if(MV(hartree_potential_storage, use_rho_multipole_shmem)){
    _FW_(double, MV(hartree_potential_storage, rho_multipole_shmem_ptr), 
      ((l_pot_max+1)*(l_pot_max+1) * (n_max_radial+2) * n_rho_multipole_atoms), rho_multipole_h_p_s, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  } else {
    _FW_(double, MV(hartree_potential_storage, rho_multipole), 
      ((l_pot_max+1)*(l_pot_max+1) * (n_max_radial+2) * n_rho_multipole_atoms), rho_multipole_h_p_s, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  }
  sum_up_pre_processing_init_();

  // ----------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------

  // int vid = 0;
#ifdef S_DEBUG
  init_opencl_util_mpi_();
  int vnum = 1;
#else
  int vnum = MV(mpi_tasks, mpi_task_per_gpu);
#endif

  #undef n_my_batches_work_sumup
  // #undef n_max_batch_size
  #undef n_full_points_work_sumup
  #undef Fp_max_grid
  #undef lmax_Fp
  #undef Fp_grid_min
  #undef Fp_grid_inc
  #undef Fp_grid_max
  #undef batches_size_sumup
  
  int i_full_points[16] = {0};
  int i_valid_points[16] = {0};

  for(int vid = 0; vid < vnum; vid++){

    // printf("-------------vid=%d--------------\n", vid);

    _FWV_(int, ocl_util_vars_all[vid].batches_size_sumup, ocl_util_vars_all[vid].n_my_batches_work_sumup, cl_buf_sumup[vid].batches_size_sumup
        , CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR); // FLAG 不同
    _FWV_(double, ocl_util_vars_all[vid].batches_points_coords_sumup, 3 * n_max_batch_size * ocl_util_vars_all[vid].n_my_batches_work_sumup
        , cl_buf_sumup[vid].batches_points_coords_sumup, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR); // FLAG 不同

    // sumup
    _FWV_(double, ocl_util_vars_all[vid].partition_tab, ocl_util_vars_all[vid].n_full_points_work_sumup, cl_buf_sumup[vid].partition_tab_std, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
    _FWV_(double, ocl_util_vars_all[vid].delta_v_hartree, ocl_util_vars_all[vid].n_full_points_work_sumup, cl_buf_sumup[vid].delta_v_hartree, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR); // 输出
    _FWV_(double, ocl_util_vars_all[vid].rho_multipole, ocl_util_vars_all[vid].n_full_points_work_sumup, cl_buf_sumup[vid].rho_multipole, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);   // 输出

    _FWV_(double, ocl_util_vars_all[vid].adap_outer_radius_sq, n_atoms, cl_buf_sumup[vid].adap_outer_radius_sq, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
    _FWV_(double, ocl_util_vars_all[vid].multipole_radius_sq, n_atoms, cl_buf_sumup[vid].multipole_radius_sq, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
    _FWV_(int, ocl_util_vars_all[vid].l_hartree_max_far_distance, n_atoms, cl_buf_sumup[vid].l_hartree_max_far_distance, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
    _FWV_(double, ocl_util_vars_all[vid].outer_potential_radius, (l_pot_max + 1) * n_atoms, cl_buf_sumup[vid].outer_potential_radius, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
    _FWV_(double, ocl_util_vars_all[vid].multipole_c, n_cc_lm_ijk(l_pot_max) * n_atoms, cl_buf_sumup[vid].multipole_c, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);

    // TODO 这一块也想办法提前面去
    int *point_to_i_batch = (int *)malloc(sizeof(int) * ocl_util_vars_all[vid].n_full_points_work_sumup);
    int *point_to_i_index = (int *)malloc(sizeof(int) * ocl_util_vars_all[vid].n_full_points_work_sumup);
    int *valid_point_to_i_full_point = (int *)malloc(sizeof(int) * ocl_util_vars_all[vid].n_full_points_work_sumup);

    for (int i_batch = 1; i_batch <= ocl_util_vars_all[vid].n_my_batches_work_sumup; i_batch++) {
      for (int i_index = 1; i_index <= ocl_util_vars_all[vid].batches_size_sumup[i_batch-1]; i_index++) {
        point_to_i_batch[i_full_points[vid]] = i_batch;
        point_to_i_index[i_full_points[vid]] = i_index;
        i_full_points[vid]++;
        if(ocl_util_vars_all[vid].partition_tab[i_full_points[vid]-1] > 0.0){
          valid_point_to_i_full_point[i_valid_points[vid]] = i_full_points[vid]-1;
          i_valid_points[vid]++;
        }
      }
    }
    if(i_full_points[vid] > ocl_util_vars_all[vid].n_full_points_work_sumup){
      printf("rank%d, i_full_points=%d, n_full_points_work_sumup=%d\n", myid+vid, i_full_points[vid], ocl_util_vars_all[vid].n_full_points_work_sumup);
      fflush(stdout);
      char save_file_name[64];
      sprintf(save_file_name, "mdata_outer_rank%d_%d.bin", myid+vid, 101010);
      m_save_load(save_file_name, 0, 1);
      m_save_load_sumup(save_file_name, 0, 1);
      fflush(stdout);
      exit(-1);
    }

    // TODO 验证一下大小，看看这一块也能不能放前面去，4 * n_full_points_work_sumup / (1024 * 1024) < 20 就行
    // loop helper

    _FWV_(int, point_to_i_batch, ocl_util_vars_all[vid].n_full_points_work_sumup, cl_buf_sumup[vid].point_to_i_batch, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
    _FWV_(int, point_to_i_index, ocl_util_vars_all[vid].n_full_points_work_sumup, cl_buf_sumup[vid].point_to_i_index, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
    _FWV_(int, valid_point_to_i_full_point, ocl_util_vars_all[vid].n_full_points_work_sumup, cl_buf_sumup[vid].valid_point_to_i_full_point, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);

    free(point_to_i_batch);
    free(point_to_i_index);
    free(valid_point_to_i_full_point);
  }

  gettimeofday(&end, NULL);
  time_uses[time_index] = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  time_infos[time_index++] = "sumup writebuf and setarg";
  if(myid < 8 && opencl_util_debug_io == 1){
    printf("rank%d, %s: %lf seconds\n", myid, time_infos[time_index-1], time_uses[time_index-1]/1000000.0);
    fflush(stdout);
  }
  gettimeofday(&start, NULL);

  cl_event event;
  long time_preproc = 0;
  long time_main_kernel = 0;

  for(int i_center_tile = 0; i_center_tile < n_centers_hartree_potential; i_center_tile += i_center_tile_size){
    // TODO 改用 kernel 版的 preprocessing，注意设备端数组指针的使用

    int i_center_begin = i_center_tile;
    int i_center_end = MIN(i_center_tile + i_center_tile_size, n_centers_hartree_potential);

    struct timeval start, end;
    gettimeofday(&start, NULL);

    // OpenCL version
    sum_up_pre_processing_part_(&n_coeff_hartree, &i_center_begin, &i_center_end, 
      &cl_buf_com.centers_rho_multipole_spl, &cl_buf_com.centers_delta_v_hart_part_spl, 
      &cl_buf_com.i_center_to_centers_index, 1);

    gettimeofday(&end, NULL);
    time_preproc += 1000000 * ( end.tv_sec - start.tv_sec ) + end.tv_usec - start.tv_usec;

    error = clSetKernelArg(kernels[0], arg_index_i_center_begin, sizeof(int), &i_center_begin);
    error = clSetKernelArg(kernels[0], arg_index_i_center_begin+1, sizeof(int), &i_center_end);

    // -------------------------------------

    for(int vid = 0; vid < vnum; vid++){

      error = clSetKernelArg(kernels[0], remember_arg_n_my_batches_work_sumup, sizeof(int), &ocl_util_vars_all[vid].n_my_batches_work_sumup);
      error = clSetKernelArg(kernels[0], remember_arg_n_my_batches_work_sumup+1, sizeof(int), &ocl_util_vars_all[vid].n_full_points_work_sumup);

      error = clSetKernelArg(kernels[0], remember_arg_Fp_max_grid, sizeof(int), &ocl_util_vars_all[vid].Fp_max_grid);
      error = clSetKernelArg(kernels[0], remember_arg_Fp_max_grid+1, sizeof(int), &ocl_util_vars_all[vid].lmax_Fp);
      error = clSetKernelArg(kernels[0], remember_arg_Fp_max_grid+2, sizeof(double), &ocl_util_vars_all[vid].Fp_grid_min);
      error = clSetKernelArg(kernels[0], remember_arg_Fp_max_grid+3, sizeof(double), &ocl_util_vars_all[vid].Fp_grid_inc);
      error = clSetKernelArg(kernels[0], remember_arg_Fp_max_grid+4, sizeof(double), &ocl_util_vars_all[vid].Fp_grid_max);

      error = clSetKernelArg(kernels[0], remember_arg_index_1, sizeof(cl_mem), &cl_buf_sumup[vid].batches_size_sumup);
      error = clSetKernelArg(kernels[0], remember_arg_index_1+1, sizeof(cl_mem), &cl_buf_sumup[vid].batches_points_coords_sumup);

      // sumup specific
      arg_index = 0;
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &sum_up_param.forces_on);
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_sumup[vid].partition_tab_std);
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_sumup[vid].delta_v_hartree);
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_sumup[vid].rho_multipole);
      // arg_index_centers_rho_multipole_spl = arg_index+4;
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.centers_rho_multipole_spl);
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.centers_delta_v_hart_part_spl);
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_sumup[vid].adap_outer_radius_sq);
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_sumup[vid].multipole_radius_sq);
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_sumup[vid].l_hartree_max_far_distance);
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_sumup[vid].outer_potential_radius);
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_sumup[vid].multipole_c);

      arg_index = 58;
      // 有 partition_tab 预判： i_valid_points ，没有： i_full_points
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(int), &i_valid_points[vid]); // TODO 注意这个 ！！！, 有没有 partition_tab 的预判不同 !
      IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed")
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_sumup[vid].point_to_i_batch);
      IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed")
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_sumup[vid].point_to_i_index);
      IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed")
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_sumup[vid].valid_point_to_i_full_point);
      IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed")
      error = clSetKernelArg(kernels[0], arg_index++, sizeof(cl_mem), &cl_buf_com.index_cc_aos);
      IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed")

      // -------------------------------------

      gettimeofday(&start, NULL);

      error = clEnqueueNDRangeKernel(cQ, kernels[0], 1, NULL, globalSize, localSize, 0, NULL, &event);
      IF_ERROR_EXIT(error != CL_SUCCESS, error, "clEnqueueNDRangeKernel failed")
      clFinish(cQ);

      gettimeofday(&end, NULL);
      long time_use = 1000000 * ( end.tv_sec - start.tv_sec ) + end.tv_usec - start.tv_usec;
      time_main_kernel += time_use;
    }
    // if(myid == 0)
    //   printf("    rank%d, %s: %lf seconds\n", myid, "sumup-kernel-tile", time_use/1000000.0);
  }

  // // WARNING: event
  // clEnqueueReadBuffer(cQ, cl_buf_com.delta_v_hartree, CL_TRUE, 0, sizeof(double) * n_full_points_work_sumup, sum_up_param.delta_v_hartree, 0, NULL, NULL);
  // clEnqueueReadBuffer(cQ, cl_buf_com.rho_multipole, CL_TRUE, 0, sizeof(double) * n_full_points_work_sumup, sum_up_param.rho_multipole, 0, NULL, NULL);

for(int vid = 0; vid < vnum; vid++){
  clEnqueueReadBuffer(cQ, cl_buf_sumup[vid].delta_v_hartree, CL_TRUE, 0, sizeof(double) * ocl_util_vars_all[vid].n_full_points_work_sumup, ocl_util_vars_all[vid].delta_v_hartree, 1, &event, NULL);
  clEnqueueReadBuffer(cQ, cl_buf_sumup[vid].rho_multipole, CL_TRUE, 0, sizeof(double) * ocl_util_vars_all[vid].n_full_points_work_sumup, ocl_util_vars_all[vid].rho_multipole, 1, &event, NULL);
}
  clFinish(cQ);
  // delta_v_hartree 和 rho_multipole 靠 use_host_ptr 直接同步

  gettimeofday(&end, NULL);
  time_uses[time_index] = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  time_infos[time_index++] = "sumup kernel and readbuf";
  if(myid < 8 && opencl_util_debug_io == 1){
    printf("rank%d, %s: %lf seconds, preproc %lf s, main kernel %lf s\n", myid, time_infos[time_index-1], time_uses[time_index-1]/1000000.0,
      time_preproc/1000000.0, time_main_kernel/1000000.0);
    fflush(stdout);
  }
  gettimeofday(&start, NULL);

  free(index_cc_aos);

  // free(local_centers_rho_multipole_spl);
  // free(local_centers_delta_v_hart_part_spl);
  free(spl_atom_to_i_center);
  free(i_center_to_centers_index);

  sum_up_pre_processing_finish();

for(int vid = 0; vid < vnum; vid++){
  clReleaseMemObject(cl_buf_sumup[vid].batches_size_sumup);
  clReleaseMemObject(cl_buf_sumup[vid].batches_points_coords_sumup);

  clReleaseMemObject(cl_buf_sumup[vid].partition_tab_std);
  clReleaseMemObject(cl_buf_sumup[vid].delta_v_hartree);
  clReleaseMemObject(cl_buf_sumup[vid].rho_multipole);
  clReleaseMemObject(cl_buf_sumup[vid].adap_outer_radius_sq);
  clReleaseMemObject(cl_buf_sumup[vid].multipole_radius_sq);
  clReleaseMemObject(cl_buf_sumup[vid].l_hartree_max_far_distance);
  clReleaseMemObject(cl_buf_sumup[vid].outer_potential_radius);
  clReleaseMemObject(cl_buf_sumup[vid].multipole_c);

  clReleaseMemObject(cl_buf_sumup[vid].point_to_i_batch);
  clReleaseMemObject(cl_buf_sumup[vid].point_to_i_index);
  clReleaseMemObject(cl_buf_sumup[vid].valid_point_to_i_full_point);
}

  clReleaseMemObject(cl_buf_com.index_cc_aos);

  clReleaseMemObject(cl_buf_com.centers_rho_multipole_spl);
  clReleaseMemObject(cl_buf_com.centers_delta_v_hart_part_spl);

  clReleaseMemObject(cl_buf_com.Fp);
  clReleaseMemObject(cl_buf_com.coord_c);
  clReleaseMemObject(cl_buf_com.coord_mat);
  clReleaseMemObject(cl_buf_com.rest_mat);
  clReleaseMemObject(cl_buf_com.vector);
  clReleaseMemObject(cl_buf_com.delta_v_hartree_multipole_component);
  clReleaseMemObject(cl_buf_com.rho_multipole_component);
  clReleaseMemObject(cl_buf_com.ylm_tab);

  clReleaseMemObject(cl_buf_com.rho_multipole_h_p_s);
  // char save_file_name[64];
  // sprintf(save_file_name, "sumup_check_rank%d_%d.bin", myid, sumup_c_count);
  // FILE *file_p = fopen(save_file_name, "w");
  // for (int i = 0; i < n_full_points_work_stable; i++) {
  //   fprintf(file_p, "%6d, %.13lf, %.13lf\n", i, sum_up_param.delta_v_hartree[i], sum_up_param.rho_multipole[i]);
  // }
  // fclose(file_p);

  if (myid == 0)
    m_save_check_sumup_(sum_up_param.delta_v_hartree, sum_up_param.rho_multipole);

  // printf("End\n");

  gettimeofday(&end, NULL);
  time_uses[time_index] = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  time_infos[time_index++] = "sumup write check file(for debug)";
  if(myid < 8 && opencl_util_debug_io == 1){
    printf("rank%d, %s: %lf seconds\n", myid, time_infos[time_index-1], time_uses[time_index-1]/1000000.0);
    fflush(stdout);
  }
  // for(size_t i=0; i<time_index; i++){
  //   printf("rank%d, %s: %lf seconds\n", myid, time_infos[i], time_uses[i]/1000000.0);
  // }
}

void sum_up_final_end(){
  clReleaseMemObject(cl_buf_com.Fp_function_spline_slice);
  clReleaseMemObject(cl_buf_com.Fpc_function_spline_slice);
}

static int rho_pass_vars_count = -1;
void rho_pass_vars_(
  int* l_ylm_max,
  int* n_local_matrix_size,
  int* n_basis_local,
  int* perm_n_full_points,
  int* first_order_density_matrix_size,
  int* basis_l_max,
  int* n_points_all_batches,
  int* n_batch_centers_all_batches,
  int* batch_center_all_batches,
  int* batch_point_to_i_full_point,
  int* ins_idx_all_batches,
  double* first_order_rho,
  double* first_order_density_matrix,
  double* partition_tab
){
  rho_pass_vars_count++;

  rho_param.l_ylm_max = *l_ylm_max;
  rho_param.n_local_matrix_size = *n_local_matrix_size;
  rho_param.n_basis_local = *n_basis_local;
  rho_param.perm_n_full_points = *perm_n_full_points;
  rho_param.first_order_density_matrix_size = *first_order_density_matrix_size;
  rho_param.basis_l_max = basis_l_max;
  rho_param.n_points_all_batches = n_points_all_batches;
  rho_param.n_batch_centers_all_batches = n_batch_centers_all_batches;
  rho_param.batch_center_all_batches = batch_center_all_batches;
  rho_param.batch_point_to_i_full_point = batch_point_to_i_full_point;
  rho_param.ins_idx_all_batches = ins_idx_all_batches;
  rho_param.first_order_rho = first_order_rho;
  rho_param.first_order_density_matrix = first_order_density_matrix;
  rho_param.partition_tab = partition_tab;
  char save_file_name[64];
  if((myid == 0 || myid == 3) && rho_pass_vars_count <= 1){
    sprintf(save_file_name, "mdata_outer_rank%d_%d.bin", myid, rho_pass_vars_count);
    m_save_load_rho(save_file_name, 0, 1);
  }
}

static int H_pass_vars_count = -1;
void h_pass_vars_(
  int* j_coord_,
  int* n_spin_,
  int* l_ylm_max_,
  int* n_basis_local_,
  int* n_matrix_size_,
  int* basis_l_max,
  int* n_points_all_batches,
  int* n_batch_centers_all_batches,
  int* batch_center_all_batches,
  int* ins_idx_all_batches,
  int* batches_batch_i_basis_h,
  double* partition_all_batches,
  double* first_order_H,
  double* local_potential_parts_all_points,
  double* local_first_order_rho_all_batches,
  double* local_first_order_potential_all_batches,
  double* local_dVxc_drho_all_batches,
  double* local_rho_gradient,
  double* first_order_gradient_rho
){
  H_pass_vars_count++;

  H_param.j_coord = *j_coord_;
  H_param.n_spin = *n_spin_;
  H_param.l_ylm_max = *l_ylm_max_;
  H_param.n_basis_local = *n_basis_local_;
  H_param.n_matrix_size = *n_matrix_size_;
  H_param.basis_l_max = basis_l_max;
  H_param.n_points_all_batches = n_points_all_batches;
  H_param.n_batch_centers_all_batches = n_batch_centers_all_batches;
  H_param.batch_center_all_batches = batch_center_all_batches;
  H_param.ins_idx_all_batches = ins_idx_all_batches;
  H_param.batches_batch_i_basis_h = batches_batch_i_basis_h;
  H_param.partition_all_batches = partition_all_batches;
  H_param.first_order_H = first_order_H;
  H_param.local_potential_parts_all_points = local_potential_parts_all_points;
  H_param.local_first_order_rho_all_batches = local_first_order_rho_all_batches;
  H_param.local_first_order_potential_all_batches = local_first_order_potential_all_batches;
  H_param.local_dVxc_drho_all_batches = local_dVxc_drho_all_batches;
  H_param.local_rho_gradient = local_rho_gradient;
  H_param.first_order_gradient_rho = first_order_gradient_rho;
  char save_file_name[64];
  if((myid == 0 || myid == 3) && H_pass_vars_count <= 4){
    sprintf(save_file_name, "mdata_outer_rank%d_%d.bin", myid, H_pass_vars_count);
    m_save_load_H(save_file_name, 0, 1);
  }
}

void rho_first_begin(){
  if(rho_first_begin_finished)
    return;
  rho_first_begin_finished = 1;
  opencl_init_();
  opencl_common_buffer_init_();

  int arg_index;
  cl_int error = 0;
  int fast_ylm = 1;
  int new_ylm = 0;

  arg_index = 13;
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &n_centers);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &n_centers_integrals);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &n_max_compute_fns_ham);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &n_basis_fns);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &n_centers_basis_I);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &n_max_grid);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &n_max_compute_atoms);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &n_max_compute_ham);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &n_max_compute_dens);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &n_max_batch_size);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &index_hamiltonian_dim2);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &position_in_hamiltonian_dim1);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &position_in_hamiltonian_dim2);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &column_index_hamiltonian_size);
  _CHK_(error, CL_SUCCESS);

  arg_index = 29;
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.center_to_atom);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.species_center);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.center_to_cell);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.Cbasis_to_basis);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.Cbasis_to_center);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.centers_basis_integrals);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.index_hamiltonian);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.position_in_hamiltonian);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.column_index_hamiltonian);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.coords_center);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.n_grid);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.r_grid_min);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.log_r_grid_inc);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.perm_basis_fns_spl);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.outer_radius_sq);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.basis_fn);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.basis_l);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.atom_radius_sq);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.basis_fn_start_spl);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.basis_fn_atom);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.basis_wave_ordered);
  _CHK_(error, CL_SUCCESS);
}

void rho_begin(){
  rho_first_begin();
  rho_first_begin_finished = 0;

  long time_uses[32];
  char *time_infos[32];
  size_t time_index = 0;
  struct timeval start, end;
  gettimeofday(&start, NULL);

  cl_int error;
  int arg_index;

  size_t localSize[] = {256}; // 覆盖前面的设置
  size_t globalSize[] = {256 * 128};  // 覆盖前面的设置
  if(myid == 0)
    printf("n_max_batch_size=%d\n", n_max_batch_size);
  // rho tmp
  _FW_(double, NULL, globalSize[0] * (n_max_compute_atoms), dist_tab_sq__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (n_max_compute_atoms), dist_tab__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (3 * n_max_compute_atoms), dir_tab__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_atoms), atom_index__, CL_MEM_READ_WRITE);
  // _FW_(int, NULL, globalSize[0] * (n_centers_integrals), atom_index_inv__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, 1, atom_index_inv__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, 1, i_basis_fns__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_basis_fns * (n_max_compute_atoms+1)), i_basis_fns_inv__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, 1, i_atom_fns__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_atoms), spline_array_start__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_atoms), spline_array_end__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, 1, one_over_dist_tab__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, 1, rad_index__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_fns_ham), wave_index__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, 1, l_index__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_fns_ham), l_count__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_fns_ham), fn_atom__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_ham), zero_index_point__, CL_MEM_READ_WRITE);
  // _FW_(double, NULL, (globalSize[0] * (n_max_compute_ham) + 128), wave__, CL_MEM_READ_WRITE); // 多加 128 为了避免 TILE 后越界
  _FW_(double, NULL, (globalSize[0]/localSize[0]) * ((n_max_batch_size+127)/128*128) * ((n_max_compute_ham+127)/128*128) + 256 + 16 * n_max_compute_ham, 
                wave__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, (globalSize[0]/localSize[0]) * (n_max_compute_dens * n_max_compute_dens) + 256 + 16 * n_max_compute_ham, 
                first_order_density_matrix_con__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (n_max_compute_atoms), i_r__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (4 * n_max_compute_atoms), trigonom_tab__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (n_max_compute_fns_ham), radial_wave__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (n_basis_fns), spline_array_aux__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (n_max_compute_atoms * n_basis_fns), aux_radial__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * ((rho_param.l_ylm_max + 1) * (rho_param.l_ylm_max + 1) * n_max_compute_atoms), ylm_tab__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * ((rho_param.l_ylm_max + 1) * (rho_param.l_ylm_max + 1) * n_max_compute_atoms), dylm_dtheta_tab__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, ((rho_param.l_ylm_max + 1) * (rho_param.l_ylm_max + 1) * n_max_compute_atoms), scaled_dylm_dphi_tab__, CL_MEM_READ_WRITE); // 反正这个和上面那个暂时无用
  // _FW_(double, NULL, (globalSize[0]/localSize[0]) * ((n_max_batch_size+127)/128*128) * ((n_max_compute_ham+127)/128*128) + 128, 
  //               tmp_rho__, CL_MEM_READ_WRITE);  // only for swcl
  // rho tmp
  arg_index = 54;
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.dist_tab_sq__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.dist_tab__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.dir_tab__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.atom_index__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.atom_index_inv__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.i_basis_fns__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.i_basis_fns_inv__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.i_atom_fns__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.spline_array_start__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.spline_array_end__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.one_over_dist_tab__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.rad_index__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.wave_index__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.l_index__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.l_count__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.fn_atom__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.zero_index_point__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.wave__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.first_order_density_matrix_con__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.i_r__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.trigonom_tab__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.radial_wave__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.spline_array_aux__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.aux_radial__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.ylm_tab__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.dylm_dtheta_tab__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.scaled_dylm_dphi_tab__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &max_n_batch_centers);
  // clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.tmp_rho__);  // only for swcl

  // clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(double) * 1024, NULL); // local mem

  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed");

  // rho param
  _FW_(int, rho_param.basis_l_max, n_species, basis_l_max__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, rho_param.n_points_all_batches, n_my_batches_work_rho, n_points_all_batches__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, rho_param.n_batch_centers_all_batches, n_my_batches_work_rho, n_batch_centers_all_batches__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, rho_param.batch_center_all_batches, max_n_batch_centers * n_my_batches_work_rho, batch_center_all_batches__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, rho_param.batch_point_to_i_full_point, n_max_batch_size * n_my_batches_work_rho, batch_point_to_i_full_point__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  if(rho_param.n_basis_local > 0){
    _FW_(int, rho_param.ins_idx_all_batches, rho_param.n_basis_local * n_my_batches_work_rho, ins_idx_all_batches__, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);
    _FW_(double, rho_param.first_order_rho, rho_param.perm_n_full_points, first_order_rho__, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);  // 输出
  } else {
    _FW_(int, NULL, 1, ins_idx_all_batches__, CL_MEM_READ_WRITE);  // 废弃
    _FW_(double, rho_param.first_order_rho, n_full_points_work_rho, first_order_rho__, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);  // 输出
  }
  _FW_(double, rho_param.first_order_density_matrix, rho_param.first_order_density_matrix_size, first_order_density_matrix__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, rho_param.partition_tab, n_full_points_work_rho, partition_tab__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);

  // rho param
  arg_index = 0;
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &rho_param.l_ylm_max);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &rho_param.n_local_matrix_size);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &rho_param.n_basis_local);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &rho_param.first_order_density_matrix_size);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.basis_l_max__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.n_points_all_batches__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.n_batch_centers_all_batches__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.batch_center_all_batches__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.batch_point_to_i_full_point__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.ins_idx_all_batches__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.first_order_rho__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.first_order_density_matrix__);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.partition_tab__);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed");

  // rho batches
  _FW_(int, MV(opencl_util, batches_size_rho), n_my_batches_work_rho, batches_size_rho, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  _FW_(int, MV(opencl_util, batches_batch_n_compute_rho), n_my_batches_work_rho, batches_batch_n_compute_rho, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  // _FW_(int, MV(opencl_util, batches_batch_i_basis_rho), n_centers_basis_I * n_my_batches_work_rho, batches_batch_i_basis_rho, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  _FW_(int, MV(opencl_util, batches_batch_i_basis_rho), n_max_compute_dens * n_my_batches_work_rho, batches_batch_i_basis_rho, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  _FW_(double, MV(opencl_util, batches_points_coords_rho), 3 * n_max_batch_size * n_my_batches_work_rho, batches_points_coords_rho, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  // rho batches
  arg_index = 50;
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.batches_size_rho);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.batches_batch_n_compute_rho);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.batches_batch_i_basis_rho);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(cl_mem), &cl_buf_com.batches_points_coords_rho);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed");

  arg_index = 27;
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &n_my_batches_work_rho);
  clSetKernelArgWithCheck(kernels[1], arg_index++, sizeof(int), &n_full_points_work_rho);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed");

  gettimeofday(&end, NULL);
  time_uses[time_index] = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  time_infos[time_index++] = "rho writebuf and setarg";
  if(myid < 8 && opencl_util_debug_io == 1){
    printf("rank%d, %s: %lf seconds\n", myid, time_infos[time_index-1], time_uses[time_index-1]/1000000.0);
    fflush(stdout);
  }
  gettimeofday(&start, NULL);

  // printf("start rho kernel\n");
  // fflush(stdout);

  cl_event event;
  error = clEnqueueNDRangeKernel(cQ, kernels[1], 1, NULL, globalSize, localSize, 0, NULL, &event);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clEnqueueNDRangeKernel failed")
  if(rho_param.n_basis_local > 0){
    clEnqueueReadBuffer(cQ, cl_buf_com.first_order_rho__, CL_TRUE, 0, sizeof(double) * rho_param.perm_n_full_points, rho_param.first_order_rho, 1, &event, NULL);
  } else {
    clEnqueueReadBuffer(cQ, cl_buf_com.first_order_rho__, CL_TRUE, 0, sizeof(double) * n_full_points_work_rho, rho_param.first_order_rho, 1, &event, NULL);
  }
  clFinish(cQ);

  gettimeofday(&end, NULL);
  time_uses[time_index] = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  time_infos[time_index++] = "rho kernel and readbuf";
  if(myid < 8 && opencl_util_debug_io == 1){
    printf("rank%d, %s: %lf seconds\n", myid, time_infos[time_index-1], time_uses[time_index-1]/1000000.0);
    fflush(stdout);
  }
  gettimeofday(&start, NULL);


  clReleaseMemObject(cl_buf_com.batches_size_rho);
  clReleaseMemObject(cl_buf_com.batches_batch_n_compute_rho);
  clReleaseMemObject(cl_buf_com.batches_batch_i_basis_rho);
  clReleaseMemObject(cl_buf_com.batches_points_coords_rho);

  clReleaseMemObject(cl_buf_com.basis_l_max__);
  clReleaseMemObject(cl_buf_com.n_points_all_batches__);
  clReleaseMemObject(cl_buf_com.n_batch_centers_all_batches__);
  clReleaseMemObject(cl_buf_com.batch_center_all_batches__);
  clReleaseMemObject(cl_buf_com.batch_point_to_i_full_point__);
  clReleaseMemObject(cl_buf_com.ins_idx_all_batches__);
  clReleaseMemObject(cl_buf_com.first_order_rho__);
  clReleaseMemObject(cl_buf_com.first_order_density_matrix__);
  clReleaseMemObject(cl_buf_com.partition_tab__);

  clReleaseMemObject(cl_buf_com.dist_tab_sq__);
  clReleaseMemObject(cl_buf_com.dist_tab__);
  clReleaseMemObject(cl_buf_com.dir_tab__);
  clReleaseMemObject(cl_buf_com.atom_index__);
  clReleaseMemObject(cl_buf_com.atom_index_inv__);
  clReleaseMemObject(cl_buf_com.i_basis_fns__);
  clReleaseMemObject(cl_buf_com.i_basis_fns_inv__);
  clReleaseMemObject(cl_buf_com.i_atom_fns__);
  clReleaseMemObject(cl_buf_com.spline_array_start__);
  clReleaseMemObject(cl_buf_com.spline_array_end__);
  clReleaseMemObject(cl_buf_com.one_over_dist_tab__);
  clReleaseMemObject(cl_buf_com.rad_index__);
  clReleaseMemObject(cl_buf_com.wave_index__);
  clReleaseMemObject(cl_buf_com.l_index__);
  clReleaseMemObject(cl_buf_com.l_count__);
  clReleaseMemObject(cl_buf_com.fn_atom__);
  clReleaseMemObject(cl_buf_com.zero_index_point__);
  clReleaseMemObject(cl_buf_com.wave__);
  clReleaseMemObject(cl_buf_com.first_order_density_matrix_con__);
  clReleaseMemObject(cl_buf_com.i_r__);
  clReleaseMemObject(cl_buf_com.trigonom_tab__);
  clReleaseMemObject(cl_buf_com.radial_wave__);
  clReleaseMemObject(cl_buf_com.spline_array_aux__);
  clReleaseMemObject(cl_buf_com.aux_radial__);
  clReleaseMemObject(cl_buf_com.ylm_tab__);
  clReleaseMemObject(cl_buf_com.dylm_dtheta_tab__);
  clReleaseMemObject(cl_buf_com.scaled_dylm_dphi_tab__);
  // clReleaseMemObject(cl_buf_com.tmp_rho__);  // only for swcl

  if (myid == 0)
    m_save_check_rho_(rho_param.first_order_rho);

  // printf("End\n");

  gettimeofday(&end, NULL);
  time_uses[time_index] = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  time_infos[time_index++] = "rho write check file(for debug)";
  if(myid < 8 && opencl_util_debug_io == 1){
    printf("rank%d, %s: %lf seconds\n", myid, time_infos[time_index-1], time_uses[time_index-1]/1000000.0);
    fflush(stdout);
  }
  // for(size_t i=0; i<time_index; i++){
  //   if(myid < 8)
  //     printf("rank%d, %s: %lf seconds\n", myid, time_infos[i], time_uses[i]/1000000.0);
  // }
}

void H_first_begin(){
  if(H_first_begin_finished)
    return;
  H_first_begin_finished = 1;
  opencl_init_();
  opencl_common_buffer_init_();

  int arg_index;
  cl_int error = 0;
  int fast_ylm = 1;
  int new_ylm = 0;

  arg_index = 19;
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &n_centers);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &n_centers_integrals);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &n_max_compute_fns_ham);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &n_basis_fns);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &n_centers_basis_I);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &n_max_grid);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &n_max_compute_atoms);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &n_max_compute_ham);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &n_max_compute_dens);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &n_max_batch_size);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &index_hamiltonian_dim2);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &position_in_hamiltonian_dim1);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &position_in_hamiltonian_dim2);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &column_index_hamiltonian_size);
  _CHK_(error, CL_SUCCESS);

  arg_index = 35;
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.center_to_atom);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.species_center);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.center_to_cell);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.Cbasis_to_basis);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.Cbasis_to_center);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.centers_basis_integrals);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.index_hamiltonian);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.position_in_hamiltonian);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.column_index_hamiltonian);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.coords_center);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.n_grid);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.r_grid_min);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.log_r_grid_inc);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.perm_basis_fns_spl);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.outer_radius_sq);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.basis_fn);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.basis_l);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.atom_radius_sq);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.basis_fn_start_spl);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.basis_fn_atom);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.basis_wave_ordered);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.basis_kinetic_ordered);  // new
  _CHK_(error, CL_SUCCESS);
}


void h_begin_0_(){
  if(h_begin_0_finished)
    return;
  h_begin_0_finished = 1;

  H_first_begin();
  H_first_begin_finished = 0;

  cl_int error;

  // H param
  _FW_(int, H_param.basis_l_max, n_species, basis_l_max__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, H_param.n_points_all_batches, n_my_batches_work_h, n_points_all_batches__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, H_param.n_batch_centers_all_batches, n_my_batches_work_h, n_batch_centers_all_batches__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(int, H_param.batch_center_all_batches, max_n_batch_centers * n_my_batches_work_h, batch_center_all_batches__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  // _FW_(int, H_param.batch_point_to_i_full_point, n_max_batch_size * n_my_batches_work_h, batch_point_to_i_full_point__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  if(H_param.n_basis_local > 0){
    _FW_(int, H_param.ins_idx_all_batches, H_param.n_basis_local * n_my_batches_work_h, ins_idx_all_batches__, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);
  } else {
    _FW_(int, NULL, 1, ins_idx_all_batches__, CL_MEM_READ_WRITE);  // 废弃
  }
  // _FW_(int, H_param.batches_batch_i_basis_h, (n_centers_basis_I * n_my_batches_work_h), batches_batch_i_basis_h__, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);  // 输出
  // _FW_(int, H_param.batches_batch_i_basis_h, 1, batches_batch_i_basis_h__, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);  // 输出
  _FW_(int, NULL, 1, batches_batch_i_basis_h__, CL_MEM_READ_WRITE);
  _FW_(double, H_param.partition_all_batches, (n_max_batch_size * n_my_batches_work_h), partition_all_batches__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, H_param.first_order_H, (H_param.n_matrix_size * H_param.n_spin), first_order_H__, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);
  // _FW_(double, H_param.local_potential_parts_all_points, (H_param.n_spin * n_full_points_work_h), local_potential_parts_all_points__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, H_param.local_potential_parts_all_points, 1, local_potential_parts_all_points__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, H_param.local_first_order_rho_all_batches, (H_param.n_spin * n_max_batch_size * n_my_batches_work_h), local_first_order_rho_all_batches__, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  _FW_(double, H_param.local_first_order_potential_all_batches, (n_max_batch_size * n_my_batches_work_h), local_first_order_potential_all_batches__, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  _FW_(double, H_param.local_dVxc_drho_all_batches, (3 * n_max_batch_size * n_my_batches_work_h), local_dVxc_drho_all_batches__, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  _FW_(double, H_param.local_rho_gradient, (3 * H_param.n_spin * n_max_batch_size), local_rho_gradient__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
  _FW_(double, H_param.first_order_gradient_rho, (3 * H_param.n_spin * n_max_batch_size), first_order_gradient_rho__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);

  // H batches
  _FW_(int, MV(opencl_util, batches_size_h), n_my_batches_work_h, batches_size_H, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  _FW_(int, MV(opencl_util, batches_batch_n_compute_h), n_my_batches_work_h, batches_batch_n_compute_H, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  _FW_(int, MV(opencl_util, batches_batch_i_basis_h), n_max_compute_dens * n_my_batches_work_h, batches_batch_i_basis_H, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
  _FW_(double, MV(opencl_util, batches_points_coords_h), 3 * n_max_batch_size * n_my_batches_work_h, batches_points_coords_H, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR);
}

void H_begin(){
  h_begin_0_();
  h_begin_0_finished = 0;

  long time_uses[32];
  char *time_infos[32];
  size_t time_index = 0;
  struct timeval start, end;
  gettimeofday(&start, NULL);

  cl_int error;
  int arg_index;

  size_t localSize[] = {256}; // 覆盖前面的设置
  size_t globalSize[] = {256 * 128};  // 覆盖前面的设置
  // printf("n_max_batch_size=%d\n", n_max_batch_size);
  // H tmp
  _FW_(double, NULL, globalSize[0] * (n_max_compute_atoms), dist_tab_sq__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (n_max_compute_atoms), dist_tab__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (3 * n_max_compute_atoms), dir_tab__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_atoms), atom_index__, CL_MEM_READ_WRITE);
  // _FW_(int, NULL, globalSize[0] * (n_centers_integrals), atom_index_inv__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, 1, atom_index_inv__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, 1, i_basis_fns__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_basis_fns * (n_max_compute_atoms+1)), i_basis_fns_inv__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, 1, i_atom_fns__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_atoms), spline_array_start__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_atoms), spline_array_end__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, 1, one_over_dist_tab__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, 1, rad_index__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_fns_ham), wave_index__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, 1, l_index__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_fns_ham), l_count__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_fns_ham), fn_atom__, CL_MEM_READ_WRITE);
  _FW_(int, NULL, globalSize[0] * (n_max_compute_ham), zero_index_point__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, (globalSize[0]/localSize[0]) * ((n_max_batch_size+127)/128*128) * ((n_max_compute_ham+127)/128*128) + 256 + 16 * n_max_compute_ham, 
                wave__, CL_MEM_READ_WRITE); // 多加 128 为了避免 TILE 后越界, 长宽按 128 对齐
  _FW_(double, NULL, 1, first_order_density_matrix_con__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (n_max_compute_atoms), i_r__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (4 * n_max_compute_atoms), trigonom_tab__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (n_max_compute_fns_ham), radial_wave__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (n_basis_fns), spline_array_aux__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * (n_max_compute_atoms * n_basis_fns), aux_radial__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * ((H_param.l_ylm_max + 1) * (H_param.l_ylm_max + 1) * n_max_compute_atoms), ylm_tab__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, globalSize[0] * ((H_param.l_ylm_max + 1) * (H_param.l_ylm_max + 1) * n_max_compute_atoms), dylm_dtheta_tab__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, 100, scaled_dylm_dphi_tab__, CL_MEM_READ_WRITE); // 反正这个和上面那个暂时无用
  _FW_(double, NULL, 1, kinetic_wave__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, (globalSize[0]/localSize[0]) * ((n_max_batch_size+127)/128*128)+ 256, grid_coord__, CL_MEM_READ_WRITE); // 长宽按 128 对齐
  // _FW_(double, NULL, globalSize[0] * n_max_compute_ham * H_param.n_spin, H_times_psi__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, 1, H_times_psi__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, 1, T_plus_V__, CL_MEM_READ_WRITE);
  // _FW_(double, NULL, globalSize[0] * (n_max_compute_atoms * n_basis_fns), T_plus_V__, CL_MEM_READ_WRITE);
  // _FW_(double, NULL, (globalSize[0]/localSize[0]) * ((n_max_batch_size+127)/128*128) * n_max_compute_ham, contract__, CL_MEM_READ_WRITE);
  // _FW_(double, NULL, (globalSize[0]/localSize[0]) * ((n_max_batch_size+127)/128*128) * n_max_compute_ham, wave_t__, CL_MEM_READ_WRITE);
  // _FW_(double, NULL, (globalSize[0]/localSize[0]) * n_max_compute_ham * n_max_compute_ham * H_param.n_spin, first_order_H_dense__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, 1, contract__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, 1, wave_t__, CL_MEM_READ_WRITE);
  _FW_(double, NULL, (globalSize[0]/localSize[0] + 1) * n_max_compute_ham * n_max_compute_ham + 128 * H_param.n_spin, first_order_H_dense__, CL_MEM_READ_WRITE);
  // H tmp
  arg_index = 60;
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.dist_tab_sq__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.dist_tab__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.dir_tab__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.atom_index__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.atom_index_inv__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.i_basis_fns__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.i_basis_fns_inv__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.i_atom_fns__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.spline_array_start__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.spline_array_end__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.one_over_dist_tab__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.rad_index__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.wave_index__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.l_index__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.l_count__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.fn_atom__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.zero_index_point__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.wave__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.first_order_density_matrix_con__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.i_r__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.trigonom_tab__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.radial_wave__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.spline_array_aux__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.aux_radial__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.ylm_tab__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.dylm_dtheta_tab__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.scaled_dylm_dphi_tab__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.kinetic_wave__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.grid_coord__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.H_times_psi__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.T_plus_V__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.contract__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.wave_t__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.first_order_H_dense__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &max_n_batch_centers);

  // clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(double) * 1024, NULL); // local mem

  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed");

  // H param
  arg_index = 0;
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &H_param.j_coord);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &H_param.n_spin);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &H_param.l_ylm_max);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &H_param.n_basis_local);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &H_param.n_matrix_size);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.basis_l_max__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.n_points_all_batches__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.n_batch_centers_all_batches__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.batch_center_all_batches__);
  // clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.batch_point_to_i_full_point__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.ins_idx_all_batches__);

  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.batches_batch_i_basis_h__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.partition_all_batches__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.first_order_H__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.local_potential_parts_all_points__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.local_first_order_rho_all_batches__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.local_first_order_potential_all_batches__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.local_dVxc_drho_all_batches__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.local_rho_gradient__);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.first_order_gradient_rho__);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed");

  // H batches
  arg_index = 57;
  // clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.batches_size_H);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.batches_batch_n_compute_H);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.batches_batch_i_basis_H);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(cl_mem), &cl_buf_com.batches_points_coords_H);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed");

  arg_index = 33;
  // int test_batch = 1;
  // clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &test_batch);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &n_my_batches_work_h);
  clSetKernelArgWithCheck(kernels[2], arg_index++, sizeof(int), &n_full_points_work_h);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clSetKernelArg failed");

  gettimeofday(&end, NULL);
  time_uses[time_index] = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  time_infos[time_index++] = "H writebuf and setarg";
  if(myid < 8 && opencl_util_debug_io == 1){
    printf("rank%d, %s: %lf seconds\n", myid, time_infos[time_index-1], time_uses[time_index-1]/1000000.0);
    fflush(stdout);
  }
  gettimeofday(&start, NULL);

  // printf("start H kernel\n");
  // fflush(stdout);

  cl_event event;
  error = clEnqueueNDRangeKernel(cQ, kernels[2], 1, NULL, globalSize, localSize, 0, NULL, &event);
  IF_ERROR_EXIT(error != CL_SUCCESS, error, "clEnqueueNDRangeKernel failed")
  clEnqueueReadBuffer(cQ, cl_buf_com.first_order_H__, CL_TRUE, 0, sizeof(double) * (H_param.n_matrix_size * H_param.n_spin), H_param.first_order_H, 1, &event, NULL);
  clFinish(cQ);


  gettimeofday(&end, NULL);
  time_uses[time_index] = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  time_infos[time_index++] = "H kernel and readbuf";
  if(myid < 8 && opencl_util_debug_io == 1){
    printf("rank%d, %s: %lf seconds\n", myid, time_infos[time_index-1], time_uses[time_index-1]/1000000.0);
    fflush(stdout);
  }
  gettimeofday(&start, NULL);

  clReleaseMemObject(cl_buf_com.batches_size_H);
  clReleaseMemObject(cl_buf_com.batches_batch_n_compute_H);
  clReleaseMemObject(cl_buf_com.batches_batch_i_basis_H);
  clReleaseMemObject(cl_buf_com.batches_points_coords_H);

  clReleaseMemObject(cl_buf_com.basis_l_max__);
  clReleaseMemObject(cl_buf_com.n_points_all_batches__);
  clReleaseMemObject(cl_buf_com.n_batch_centers_all_batches__);
  clReleaseMemObject(cl_buf_com.batch_center_all_batches__);
  // clReleaseMemObject(cl_buf_com.batch_point_to_i_full_point__);
  clReleaseMemObject(cl_buf_com.ins_idx_all_batches__);

  clReleaseMemObject(cl_buf_com.batches_batch_i_basis_h__);
  clReleaseMemObject(cl_buf_com.partition_all_batches__);
  clReleaseMemObject(cl_buf_com.first_order_H__);
  clReleaseMemObject(cl_buf_com.local_potential_parts_all_points__);
  clReleaseMemObject(cl_buf_com.local_first_order_rho_all_batches__);
  clReleaseMemObject(cl_buf_com.local_first_order_potential_all_batches__);
  clReleaseMemObject(cl_buf_com.local_dVxc_drho_all_batches__);
  clReleaseMemObject(cl_buf_com.local_rho_gradient__);
  clReleaseMemObject(cl_buf_com.first_order_gradient_rho__);

  clReleaseMemObject(cl_buf_com.dist_tab_sq__);
  clReleaseMemObject(cl_buf_com.dist_tab__);
  clReleaseMemObject(cl_buf_com.dir_tab__);
  clReleaseMemObject(cl_buf_com.atom_index__);
  clReleaseMemObject(cl_buf_com.atom_index_inv__);
  clReleaseMemObject(cl_buf_com.i_basis_fns__);
  clReleaseMemObject(cl_buf_com.i_basis_fns_inv__);
  clReleaseMemObject(cl_buf_com.i_atom_fns__);
  clReleaseMemObject(cl_buf_com.spline_array_start__);
  clReleaseMemObject(cl_buf_com.spline_array_end__);
  clReleaseMemObject(cl_buf_com.one_over_dist_tab__);
  clReleaseMemObject(cl_buf_com.rad_index__);
  clReleaseMemObject(cl_buf_com.wave_index__);
  clReleaseMemObject(cl_buf_com.l_index__);
  clReleaseMemObject(cl_buf_com.l_count__);
  clReleaseMemObject(cl_buf_com.fn_atom__);
  clReleaseMemObject(cl_buf_com.zero_index_point__);
  clReleaseMemObject(cl_buf_com.wave__);
  clReleaseMemObject(cl_buf_com.first_order_density_matrix_con__);
  clReleaseMemObject(cl_buf_com.i_r__);
  clReleaseMemObject(cl_buf_com.trigonom_tab__);
  clReleaseMemObject(cl_buf_com.radial_wave__);
  clReleaseMemObject(cl_buf_com.spline_array_aux__);
  clReleaseMemObject(cl_buf_com.aux_radial__);
  clReleaseMemObject(cl_buf_com.ylm_tab__);
  clReleaseMemObject(cl_buf_com.dylm_dtheta_tab__);
  clReleaseMemObject(cl_buf_com.scaled_dylm_dphi_tab__);
  clReleaseMemObject(cl_buf_com.kinetic_wave__);
  clReleaseMemObject(cl_buf_com.grid_coord__);
  clReleaseMemObject(cl_buf_com.H_times_psi__);
  clReleaseMemObject(cl_buf_com.T_plus_V__);
  clReleaseMemObject(cl_buf_com.contract__);
  clReleaseMemObject(cl_buf_com.wave_t__);
  clReleaseMemObject(cl_buf_com.first_order_H_dense__);

  if (myid == 0)
    m_save_check_h_(H_param.first_order_H, &(H_param.n_spin), &(H_param.n_matrix_size));

  // printf("End\n");

  // opencl_common_buffer_free_();

  gettimeofday(&end, NULL);
  time_uses[time_index] = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  time_infos[time_index++] = "H write check file(for debug)";
  if(myid < 8 && opencl_util_debug_io == 1){
    printf("rank%d, %s: %lf seconds\n", myid, time_infos[time_index-1], time_uses[time_index-1]/1000000.0);
    fflush(stdout);
  }
  // for(size_t i=0; i<time_index; i++){
  //   if(myid < 8)
  //     printf("rank%d, %s: %lf seconds\n", myid, time_infos[i], time_uses[i]/1000000.0);
  // }
}