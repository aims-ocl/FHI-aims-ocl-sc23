# GPU Code Notes

The following is a set of notes that I (WPH) wrote during a GPU code refactoring from 22-23 February, 2018 to align the various GPU-accelerated algorithms into one consistent style.  Hence, the disproportionate emphasis given to a style guide.

## Code Structure

The general flow of the GPU acceleration of an algorithm is:
```fortran
call ALGORITHM_create_gpu()
call set_VARIABLE_gpu()
call OPERATION_gpu()
call get_VARIABLE_gpu()
call ALGORITHM_destroy_gpu()
```

This naming convention has been chosen to roughly mimic the cuBLAS naming convention (cublasCreate(), cublasSetVector(), operations, cublasGetVector(), cublasDestroy()) while still being intuitive I hope to people that don't know GPU programming.

## Style Used

- 3 spaces for indentation, 6 spaces for continuation.
- 80 characters per line, but I'm not fanatical about this.
- To be consistent, all argument lists for functions are broken onto separate lines, even in cases where they could fit into one line.
- Subroutine calls across multiple lines are considered as line continuations (i.e. each additional line is prefaced by six spaces), but there is no restriction to where line breaks are placed.
  - Generally, I start the parameter list on the next line.
- All alignment should be "fixed" and independent of the contents of the code.   It is my personal experience that alignments based on anything within the code (data types, variable names, etc.) may look pretty for small bits of code, but quickly become unmanagable as the size of the code grows.
- All subroutines and variables accessible on the Fortran side should use lowercase snake_case to match aims' style.
- All subroutines and variables not visible to the Fortran side should use lowerCamelCase to match CUDA's style.  The exception is device pointers, which have "dev_" prepended to their name, and macros, which use SCREAMING_SNAKE_CASE.
- Whenever possible, names in the CUDA code should match the corresponding names in the main body of the FHI-aims source code.
- GPU-accelerated subroutines visible from FHI-aims should have a _gpu() suffix to clearly distinguish them in the main FHI-aims code.
  - In the GPU source code, there's no reason to abstract away CUDA or cuBLAS functionality, as this makes the code less readable.  It is assumed that people looking at the GPU source code should know what those two libraries are.
  - Within the main body of the FHI-aims source code, people will have no idea what a CUDA or a cuBLAS is.  Everyone knows what a GPU is.
- Fortran logicals are passed as integers, and the zero-if-false convention is used.
- Data movement from device to host should always be explicit.

# TODO
- Standardize whether array dimensions are stored in the namespace of an algorithm or whether they're passed into subroutines as arguments
- Standardize host-to-device data movement throughout algorithms
  - Sometimes it's done explicitly in the code via set_VARIABLE_gpu() calls, sometimes it's done as part of other subroutines.
  - Sometimes the number of elements are supplied, sometimes they're not.
- Implement CUDA stream in density and forces.
  - The density and density gradient calculations are completely independent of one another, as are the forces and analytical stress tensor calculations.  They could be put into different streams of their respective algorithm.
- Rename aimsSetVector, aimsGetVector, aimsSetMatrix, and aimsGetMatrix to something a little less vague
- Make generalized versions of batch integration subroutines
- Add doxygen or Sphinx documentation.
