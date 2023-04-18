#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void change_directory_c(const char* path, int* err) 
{
  *err = chdir(path);
}
