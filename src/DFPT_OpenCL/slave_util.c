#include <slave.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define UP_ALIGN_PTR(ptr) ((void *)(((long)ptr + 63) & ~63))

#define MAX_ALLOCATE_SIZE 1024 * 1024 * 10
__thread char dynamic_allocate_buf[MAX_ALLOCATE_SIZE];
__thread_local void *dynamic_allocate_pointer;
__thread_local int dynamic_allocate_pointer_inited = 0;

void init_allocate() { dynamic_allocate_pointer = dynamic_allocate_buf; }
void init_allocate_first() { 
  if(dynamic_allocate_pointer_inited == 0){
    dynamic_allocate_pointer = dynamic_allocate_buf;
    dynamic_allocate_pointer_inited = 1;
  }
}

void *thread_allocate(size_t size) {
  void *ret = dynamic_allocate_pointer;
  dynamic_allocate_pointer += size;
  dynamic_allocate_pointer = UP_ALIGN_PTR(dynamic_allocate_pointer);
  return ret;
}

void free_allocate_pointer(void *free_addr) { dynamic_allocate_pointer = free_addr; }
