#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>

void unlimit_stack_(void) {
 struct rlimit rlim = { RLIM_INFINITY, RLIM_INFINITY };
 if ( setrlimit(RLIMIT_STACK, &rlim) == -1 ) {
  perror("setrlimit error");
  exit(1);
 }
}
