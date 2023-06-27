# 方便起见一般都会先定义编译器链接器
CC=swgcc
LD=swgcc
OCL=/home/export/base/shisuan/swyjs/cq/swcl/install/bin/swcl
OCL_INCLUDE=/home/export/base/shisuan/swyjs/cq/swcl/install/lib/swcl/host/include
OCL_LIB=/home/export/base/shisuan/swyjs/cq/swcl/install/lib/swcl
IFOR=tmp.conb.cl

sdir := swcl

# 正则表达式表示目录下所有.c文件，相当于：SRCS = main.c a.c b.c
HOSTSRCS = opencl_util_c.c save_load_var.c integrate_first_order_rho.c sum_up_whole_potential_shanghui.c main.c opencl_util_mpi_m.c $(sdir)/slave.info.c
SLAVESRCS = $(sdir)/slave.c  slave_util.c

# OBJS表示SRCS中把列表中的.c全部替换为.o，相当于：OBJS = main.o a.o b.o
# HOSTOBJS = $(patsubst %c, %host.o, $(HOSTSRCS))
HOSTOBJS = $(patsubst %c, %o, $(HOSTSRCS))
SLAVEOBJS = $(patsubst %c, %slave.o, $(SLAVESRCS))
# SLAVEOBJS = $(patsubst %c, %o, $(SLAVESRCS))
OBJS = $(patsubst %c, %o, $(SRCS))

# 可执行文件的名字
TARGET = a.out

# .PHONE伪目标，具体含义百度一下一大堆介绍
.PHONY:all clean

# 要生成的目标文件
all: $(TARGET)

tmp.conb.cl: sum_up.cl integrate_first_order_rho.cl
	cat $^ > tmp.conb.cl

swcl-slave: $(sdir)/slave.c

$(sdir):
	mkdir -p $@

$(sdir)/slave.c: $(IFOR) | $(sdir)
	-cd $(sdir) && $(OCL) ../$< -only-cl2slave

$(sdir)/slave.info.c: $(IFOR) | $(sdir)
	cd $(sdir) && $(OCL) ../$< -only-cl2slave

$(TARGET): $(HOSTOBJS) $(SLAVEOBJS)
	$(CC) -o $@ $^ /home/export/online1/mdt00/shisuan/swyjs/wzk/share/lib/sw_slave_dgemm.o  \
	 -lstdc++ -lm -mhybrid -O2 -lm_slave -L $(OCL_LIB) -lswclapi


%.o:%.c
	$(CC) -o $@ -c $^ -O2 -g -Wall -mhost -I $(OCL_INCLUDE) -mieee -Wno-comment -Wno-unused-variable -DS_DEBUG

%.slave.o:%.c
	$(CC) -o $@ -c $^ -O2 -c -msimd -lm -g -mslave -mieee

# make clean删除所有.o和目标文件
clean:
	rm *.out -f
	rm -f $(OBJS) $(HOSTOBJS) $(SLAVEOBJS) $(TARGET) *.o
	rm $(sdir) -rf

run: $(TARGET)
	bsub -b -o bg.out -m 1 -J aims-swcl -n 1 -I -cgsp 64 -share_size 13000 -host_stack 1024 -priv_size 16 ./a.out 
	# head -n 10000 sumup_check_rank0_-1.bin > sumup_check_rank0_10000.bin

# gcc main.c opencl_util_c.c sum_up_whole_potential_shanghui.c save_load_var.c  -I/opt/rocm/opencl/include -lOpenCL -L/opt/rocm/opencl/lib/x86_64 -lm
# gcc sum_up_1.c -I/opt/rocm/opencl/include -lOpenCL -L/opt/rocm/opencl/lib/x86_64
# srun -n 1 -p normal --gres=dcu:1 ./a.out
# code --diff ./rho_check_rank0_-1.bin /home/export/online1/mdt00/shisuan/swyjs/wzk/test/DFPT_OpenCL_20230223/test_data/rho_check_rank0_1_cl.bin