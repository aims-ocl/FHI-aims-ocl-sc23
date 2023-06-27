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
  
3. Data Processing: Manually or use a specific script to extract relevant information including the execution time of some subprocedures and so on from the output file, which is generated during execution.
  

The software installation process will take 5-30 minutes, depending on whether the dependency library is ready and the number of cores available for compilation. The execution time of the binary ranges from 10 minutes to 4 hours according to different testcases. The extraction takes about 2 minutes.

The key result of this workflow is the output file generated during execution, which includes execution time of different parts, the value of some runtime variable corresponding to memory consumption and the result of the DFPT calculation.

The numerical results found in the article include execution time of some subprocedures, speedups between and scalability efficiency, which can be directly obtained, calculated from the above results. The figures of the evaluation part in the article is visual version of the numerical results.

## The Code

Our main work on OpenCL is in directory `src/DFPT_OpenCL`.

For sunway support, pleash checkout to the opencl-develop-sw branch.

For more information about FHI-aims, please refer [Home - FHI-aims](https://fhi-aims.org/).

## How to run

### Hardware Requirement

The core part of this work is written based on OpenCL, so your hardware devices need to support OpenCL.

All calculations and scalability tests were run with FHI-aims on two supercomputers, the new generation Sunway suptercomputer and an AMD-GPU-accelerated supercomputer with a special version of AMD GPU (between MI50 and MI60, with 64 computing units like mi60, but only 16GB of memory like MI50.).

The opencl-related code of this repository is written and optimized for these two specific architectures, and has not been tested on other platforms at present. In order to achieve better reproduction effect, please try to test it on the same platform.

On other platforms, in order to optimize, you may need to try to adjust parameters including workgroup size and blocking. NVIDIA GPU is temporarily unsupported due to some OpenCL semantics, and we will add support in the future to enhance the portability of our work.

### Software Requirement

1. Install or load mathlibs to provide LAPACK and SCALAPACK support. We use MKL on the AMD-GPU-accelerated supercomputer and xMath on the new generation Sunway suptercomputer.
2. Install or load a OpenCL driver to provide OpenCL support. We use ROCM on the AMD-GPU-accelerated supercomputer and SWCL on the new generation Sunway suptercomputer.
3. Install or load a compiler supporting MPI and C+Frotran. We use intelmpi and gcc on the AMD-GPU-accelerated supercomputer and swmpi and gcc on the new generation Sunway suptercomputer.

### Compile

Check out to （没写）branch for specific 。。。。

1. Set `PATH` and `LD_LIBRARY_PATH` for the software and libraries above. We recommend that you package the environment variable design as an `env.sh` for convenience，so that you just need to `source env.sh` later. We provide the examples under `testcases-SC/envs-example/<example-name>/env.sh`.
  
2. Copy initial_cache.example.cmake as initial_cache.cmake. Modify the `CMAKE_Fortran_COMPILER`, `CMAKE_C_COMPILER`, `CMAKE_Fortran_FLAGS`, `CMAKE_C_FLAGS` options and so on in file based on the your environment. We provide the examples under `testcases-SC/envs-example/<example-name>/initial_cache.cmake`.
  

Then we can compile the code as following. We also privide examples of build scripts under `testcases-SC/envs-example/<example-name>/build.sh`.

```bash
source env.sh
cmake -C ./initial_cache.cmake -S . -B build -DCMAKE_PREFIX_PATH="" -DOpenCL_FOUND=True -DOpenCL_LIBRARY=<YOUR_PATH>/opencl/lib/ -DOpenCL_INCLUDE_DIR=<YOUR>PATH>/opencl/include/
cmake --build build -j4
```

It takes about 5-30 minutes to compile, depending on the peformance of the CPU and the number of CPU cores you use. After successful compilation, you will be able to see `build/aims.191127.scalapack.mpi.x`.

### DataSets

We mainly provide four kinds of test cases under testcases-SC, among which Polypeptide is used for basic test because it is simple, time-consuming and representative.

1. H(C2H4)nH, 6n+2 atoms, used as large cases.
  
2. HIV-1, 49 atoms, recommend to run with 1 cpu core, takes about 5 minutes.
  
3. Polypeptide, 312 atoms, recomment to run with 8 cpu cores, takes about 4 minutes.
  
4. RBD, 3006 atoms, recommend to run with 32/64 cpu cores, takes 30 minutes.
  

In each cases, we provide `control.in` and `geometry.in` which is needed by fhi-aims. We also provide `sbatch.sh` and `test-sw.sh` as examples of test scripts.

For more $H(C_2 H_4)_nH$ cases, we provide `testcases-SC/H(C2H4)nH/c2h4_generate.F` to generate. Compile it and then run it with n as the first argument to generate the case you need. For small cases (6n + 2 <= 5000), you should use the same control.in as `testcases-SC/H(C2H4)nH/c2h4-3k/control.in`. For large cases, you should use the same control.in as `testcases-SC/H(C2H4)nH/1w_atoms_mini_local_index/control.in`. Here's an example to generate $H(C_2 H_4)_{5000}H$.

We only support DFPT polar_reduce_memory now and plan to support other similar DFPT systems later.

```bash
cd testcases-SC/H(C2H4)nH
mkdir newcases
gfortran ../c2h4_generate.F
./a.out 5000 > geometry.in
cp ../1w_atoms_mini_local_index/control.in .
```

### Run (basically)

Use Polypeptide as the testcase.

We need to set some environment variables to enable some functions. To enable opencl acceleration, set FHI_AIMS_OCL to ON. To compare it with CPU version, turn it to OFF. The MPI_PER_NODE should be set the same as tasks-per-node and GPU_PER_NODE should also be set to suit the actual situation (for sunway version, they should always be set to 1 since each main cpu core is equipped with an acceleration unit). Here's an example setting environment variables.

```bash
export FHI_AIMS_OCL_RHO_MULTIPOLE_TILE_SIZE=512
export FHI_AIMS_OCL_RHO_MULTIPOLE_SHMEM=OFF
export FHI_AIMS_OCL=ON
export FHI_AIMS_OCL_DFPT_POLAR_RHO=ON
export FHI_AIMS_OCL_DFPT_DIEL_SUMUP=ON
export FHI_AIMS_OCL_DFPT_POLAR_H=ON
export MPI_PER_NODE=1
export GPU_PER_NODE=1
export FHI_AIMS_OCL_DEBUG_IO=ON
```

Modify some settings including file path in the example `sbatch.sh` or `test-sw.sh` to run fhi-aims on your computer or cluster. Here's a complete example to run.

```bash
# remember to set the same environment variables during compilation including PATH and LD_LIBRARY_PATH !
cd ./testcases-SC/Polypeptide
sbatch.sh
# the output file will be like slurm-xxxxx.out
```

All output will be output to stdout. We are mainly concerned with time information like the following. For simplicity, we set CPSCF run 2 iteration (with j_coord=1...3, totally 6 DFPT iterations). In testcases with less CPSCF iterations, the main time-consuming of the program is other parts, so we just concern the time of an iteration. The first iterations will always run on CPU only. We recommend you to see the time under `CPSCF working for j_coord = 2` and `Begin CP-self-consistency iteration # 2` . (Use a text editor to find the second occurrence of `End CPSCF iteration # 2`.) With FHI_AIMS_OCL=ON and FHI_AIMS_OCL=OFF，you can two output file and compare the time.

```
  End CPSCF iteration #     2                  :  max(cpu_time)    wall_clock(cpu1)
  | Time for this iteration                    :        5.320 s           5.320 s
  | first_order_density                        :        1.010 s           1.003 s
  | first_order_potential                      :        2.020 s           2.013 s
  | first_order_H                              :        1.890 s           1.897 s
  | first_order_DM                             :        0.400 s           0.402 s
  | Solution of Sternheimer eqns.              :        0.000 s           0.000 s
```

### Specific Experiments

#### Locality-Enhancing Task Mapping

For convenience, this experiments use a very small molecular system HIV-1 which only needs 1 cpu core and 1 acceleration unit (like a gpu). We only focus on the time of two kernels, which is small in this case but will be larger in more complex cases. It will take about 6 minutes. Then you can use grep to get each kernel time.

With locality-enhancing

```bash
cd ./testcases-SC/HIV-1
sbatch.sh
grep 'rank0, rho kernel and readbuf' <OUTPUT_FILE>
grep 'rank0, H kernel and readbuf' <OUTPUT_FILE>
```

Without locality-enhancing

```bash
# remember to set the same environment variables during compilation including PATH and LD_LIBRARY_PATH !
cd testcases-SC/HIV-1-without-local-enhance
sbatch.sh
grep 'rank0, rho kernel and readbuf' <OUTPUT_FILE>
grep 'rank0, H kernel and readbuf' <OUTPUT_FILE>
```

Then you can compare the time as Figure9 (b).

#### Packed Hierarchical Collective Communication

This optimization can only be reflected when the calculation scale is large (e.g. 8-256 nodes).

We can use these environment variables to control the optimization. Modify these variables in run script (e.g. sbatch_template.sh).

Baseline

```bash
export FHI_AIMS_OCL_RHO_MULTIPOLE_TILE_SIZE=1
export FHI_AIMS_OCL_RHO_MULTIPOLE_SHMEM=OFF
```

Packed

```bash
export FHI_AIMS_OCL_RHO_MULTIPOLE_TILE_SIZE=512
export FHI_AIMS_OCL_RHO_MULTIPOLE_SHMEM=OFF
```

Packed Hierarchical

```bash
export FHI_AIMS_OCL_RHO_MULTIPOLE_TILE_SIZE=512
export FHI_AIMS_OCL_RHO_MULTIPOLE_SHMEM=ON
```

Then you can test it.

```bash
# remember to set the same environment variables during compilation including PATH and LD_LIBRARY_PATH !
cd testcases-SC/H(C2H4)nH/3w_atoms_mini_local_index
# modify node_num and for count in {1..4} in the script to adjust the number of nodes used.
# the script will help you to test a testcases using different number of cpu and gpu
../test-sbatch.sh
 grep 'Hartree multipole update' <OUTPUT_FILE> # use the first one
```

To cut the test time, you can add `stop` after the following code in `src/scf_solver.f90` and them recompile.

```fortran
call output_times(deffmt, "Hartree multipole update", &
        &                 time_hartree_multi, clock_time_hartree_multi,OL_norm)
```

#### Scalability

```bash
# remember to set the same environment variables during compilation including PATH and LD_LIBRARY_PATH !
cd testcases-SC/H(C2H4)nH/3w_atoms_mini_local_index # or other cases 
# modify node_num and for count in {1..4} in the script to adjust the number of nodes used.
# the script will help you to test a testcases using different number of cpu and gpu
../test-sbatch.sh
```

Then you can get the time following the previous part "Run (basically)".

Some other experiments still needs further effort since it's a little difficult to change the unoptimized part from the new version that has been adjusted.