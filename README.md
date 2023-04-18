## 1. Artifact Identification

### 1.1. Contributions of the paper and how they will be justified with our artifacts.

This work re-designs the quantum perturbation calculations within FHI-aims (a all-electron full-potential framework), with a portable and scalable OpenCL implementation for the heterogeneous many-core supercomputer, with following major innovations:

1. A portable quantum perturbation implementation has been developed using OpenCL, i.e., a cross-platform unified programming framework, to enable efficient all-electron all-potential simulations across various supercomputers.
2. Major scaling limitations in quantum perturbation has been removed, by a locality-enhancing task mapping strategy. It enables neighbouring atoms to be simulated in the same MPI process, leading to significantly reduced per-process memory consumption and more efficient memory accesses.
3. Efficient collective communication has been enabled, by introducing a packed hierarchical communication scheme, which reduces the total number of collective communications by packing several of them together, and accelerates each communication by applying a hierarchical approach for inter-node and intra-node data synthesizing. 
4. A set of performance-portable OpenCL optimizations have been implemented to improve the efficiencies of compute units and memory sub-system for various heterogeneous processors, by reordering OpenCL kernel invocations, memory accesses and arithmetic operations. 
5. We perform a full *ab initio* simulation of a real material containing 3006 atoms and compare it with experimental results. We have evaluated our OpenCL-accelerated quantum perturbation simulation on two advanced supercomputers, and results show that our implementation exhibits remarkably performance on various material systems, enabling the system scale to 200,000 atoms with all-electron precision. 

Our artifacts include an optimized version, and several sub-optimized versions with certain optimization excluded. By running the optimized version, contributions 1 and 5 could be justified directly. By comparing results of the optimized version with certain sub-optimized version, contributions 2-4 could be testified with corresponding speedup, respectively. 

### 1.2. Software architecture of our artifacts.

All versions of artifacts provided leverage the same software architecture (also the same as the open source FHI-aims-ocl-sc23 (this project)) as follows.
1. DFT phase. It serves to provide data for the DFPT phase, with scf_solver as its entry point.
2. DFPT phase. It completes quantum perturbation calculation by iteratively calling following functions, until self-consistency has been reached. In our research, we only focus on the polar_reduce_memory version of DFPT, with cpscf_solver_polar_reduce_memory as its entry point.
    1. DM. It calculates the response of the density matrix, with evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory or evaluate_first_order_DM_polar_reduce_memory as its entry point, depending on the runtime settings.
    2. Sumup (OpenCL-accelerated). It calculates the response of the electronic density, with sum_up_whole_potential_shanghui_dielectric as its entry point, launching 2 OpenCL kernels sum_up_whole_potential_shanghui_pre_proc_ and sum_up_whole_potential_shanghui_sub_t_.
    3. Rho (OpenCL-accelerated). It calculates the response of the total electrostatic potential, with integrate_first_order_rho_polar_reduce_memory as its entry point, launching 1 OpenCL kernel integrate_first_order_rho_sub_t_.
    4. H (OpenCL-accelerated). It calculates the response of the Kohn-Sham Hamiltonian matrix, with integrate_first_order_H_polar_reduce_memory as its entry point, launching 1 OpenCL kernel integrate_first_order_h_sub_t_.

### 1.3. To what extent our artifacts contribute to reproductivity of experiments.

All results in our evaluation, could be reproduced by our provided artifacts, given enough computing nodes on the two HPC systems. However, speedup results could be reproduced even with a few computing nodes.

## 2. Reproducibility of Experiments (outline)

All calculations and scalability tests were run with FHI-aims on two supercomputers, the new generation Sunway suptercomputer and an AMD-GPU-accelerated supercomputer. The experiment workflow includes the following three steps. Note that different versions of the code and different test cases (including input files) should be chosen for different evaluations.

1. Build: Install and compile the libraies and the code. Aftering compiling, their will be an binary call aims.191127.scalapack.mpi.x in the bin directory of the project.
  
2. Execution: Run the binary on a cluster with specific input file. First, go to the the testcase directory which contains control file (control.in) and the geometry file (geometry.in). Then submit the calculation task to the job management system of the cluster.
  
3. Data Processing: Use a specific python script to extract relevant information including the execution time of some subprocedures and so on from the output file, which is generated during execution.
  
The installation of the software is expected to take 30 minutes. The execution time of the binary ranges from 10 minutes to 4 hours according to different testcases. The extraction takes about 2 minutes.

The key result of this workflow is the output file generated during execution, which includes execution time of different parts, the value of some runtime variable corresponding to memory consumption and the result of the DFPT calculation.

The numerical results found in the article include execution time of some subprocedures, speedups between and scalability efficiency, which can be directly obtained, calculated from the above results. The figures of the evaluation part in the article is visual version of the numerical results.

## How to build

Prepare:

1. Install or load mathlibs to provide LAPACK and SCALAPACK support. We use MKL on the AMD-GPU-accelerated supercomputer and swMathlib on the new generation Sunway suptercomputer.
2. Install or load a OpenCL driver to provide OpenCL support. We use ROCM on the AMD-GPU-accelerated supercomputer and SWCL on the new generation Sunway suptercomputer.
3. Install or load a compiler supporting MPI and C+Frotran. We use intelmpi and gcc on the AMD-GPU-accelerated supercomputer and swmpi and gcc on the new generation Sunway suptercomputer.
3. Set `PATH` and `LD_LIBRARY_PATH` for the software and libraries above.
4. Modify the `FC`, `CC`, `LAPACK` and `SCALAPACK` option in file `src/Makefile` based on the above settings.
5. Modify the `mpi_task_per_gpu` and `mpi_task_per_gpu` option in file `src/DFPT_OpenCL/opencl_util.f90` to adapt to your operating environment.

Compiler:
```sh
cd src
make allclean
make scalapack.mpi -jN # N is the thread you want to use for compiling
```

Result:
You can get an executable file `bin/aims.191127.scalapack.mpi.x`.

## How to run

Testcases are stored in the directory `testcases-SC`. `example-sbatch.sh` is an example slurm script.

An example:

Obtain the speedup of gpu relative to cpu with the testcase $H(C_2H_4)_{1667}H$, on a computer cluster using slurm as the job management system

```
cd testcases-SC/H(C2H4)nH/1w_atoms_mini_local_index
sbatch sbatch.sh
```

It takes aboue 40 minutes to finish the task. `slurm-xxx.out` in the same directoy is the output file, which is just an file redirected from stdout.

## Data Processing

`get_info.py` is an example script to help user get information from the output file, including the execution time of some subprocedures and so on. You need to modify it according to the path of the output file and the information you need. Then just use command `python3 get_info.py` to get the information you need.