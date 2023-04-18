#ifndef CUDA_MACRO_H
#define CUDA_MACRO_H

/* The following macro is used to control, via preprocessor directives,  whether 
 * an underscore should be added to the end of names of functions which will be
 * called by the main Fortran code.
 * A nice side benefit is that it clearly indicates within the code which
 * functions will be called by the main aims code base (their names will be
 * wrapped by this macro), and which are completely internal to the CUDA
 * implementation (their names will not be wrapped by this macro).
 * Generally, an underscore should be added for all compiler suites except IBM
 * XL. */
#ifdef  Add_
#define FORTRAN(name) name##_
#else
#define FORTRAN(name) name
#endif

#endif
