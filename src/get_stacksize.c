#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

void get_stacksize_c(int *maxstack, int *currentstack, int *err) {


  struct rlimit rlim;
  *err = getrlimit(RLIMIT_STACK, &rlim);

  
  *maxstack = rlim.rlim_max;
  *currentstack = rlim.rlim_cur;

}

