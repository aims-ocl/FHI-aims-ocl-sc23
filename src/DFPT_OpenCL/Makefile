# source /public/software/compiler/intel/oneapi/mkl/2021.3.0/env/vars.sh
# source /public/software/compiler/intel/oneapi/mpi/2021.3.0/env/vars.sh
# 方便起见一般都会先定义编译器链接器
# CC = mpicc
CC = gcc
# LD = gcc

# 正则表达式表示目录下所有.c文件，相当于：SRCS = main.c a.c b.c
# SRCS = $(wildcard *.c)
SRCS = opencl_util_c.c save_load_var.c integrate_first_order_rho.c sum_up_whole_potential_shanghui.c main.c opencl_util_mpi_m.c

# OBJS表示SRCS中把列表中的.c全部替换为.o，相当于：OBJS = main.o a.o b.o
OBJS = $(patsubst %c, %o, $(SRCS))

# 可执行文件的名字
TARGET = a.out

# .PHONE伪目标，具体含义百度一下一大堆介绍
.PHONY:all clean

# 要生成的目标文件
all: $(TARGET)

$(TARGET): $(SRCS)
	# $(CC) -o $@ $^ -lOpenCL -lm -g -O2
	# $(CC) -o $@ $^ -I/opt/rocm/opencl/include -lm -L/opt/rocm/opencl/lib/x86_64 -lOpenCL -g -O2
	$(CC) -o $@ $^ -I/public/software/apps/AI/compiler/rocm-4.0.1/opencl/include -lm -L/public/software/apps/AI/compiler/rocm-4.0.1/opencl/lib -lOpenCL -g -O3 -DDEBUG_NO_MPI=1
# 第一行依赖关系：冒号后面为依赖的文件，相当于Hello: main.o a.o b.o
# 第二行规则：$@表示目标文件，$^表示所有依赖文件，$<表示第一个依赖文件
# $(TARGET): $(OBJS)
# 	$(LD) -o $@ $^ -lm -lOpenCL

# # 上一句目标文件依赖一大堆.o文件，这句表示所有.o都由相应名字的.c文件自动生成
# %.o:%.c
# 	$(CC) -c $^ -O2 -g -Wall

# make clean删除所有.o和目标文件
clean:
	rm -f $(OBJS) $(TARGET)

run: $(TARGET)
	# srun -n 1 -p dongsheng --gres=dcu:1 ./a.out > log.txt
	srun -n 1 -p normal --gres=dcu:1 ./a.out | tee log.txt
	# srun -n 1 -p normal --gres=dcu:1 valgrind --log-file=valReport2 --leak-check=full --show-reachable=yes --leak-resolution=low ./a.out | tee log.txt

# gcc main.c opencl_util_c.c sum_up_whole_potential_shanghui.c save_load_var.c  -I/opt/rocm/opencl/include -lOpenCL -L/opt/rocm/opencl/lib/x86_64 -lm
# gcc sum_up_1.c -I/opt/rocm/opencl/include -lOpenCL -L/opt/rocm/opencl/lib/x86_64
# srun -n 1 -p normal --gres=dcu:1 ./a.out