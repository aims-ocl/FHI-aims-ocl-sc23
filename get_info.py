# python .\get_info.py | Tee-Object -FilePath log.get_info.1.txt

import sys
import os
import numpy as np

search_group_1 = [
    [1, "Total time                                 :", " s", " s"],
    [1, "Total time for Hartree multipole sum", " s", " s"],
    [5, "Time for this iteration", " s", " s"],
    [5, "| first_order_density", " s", " s"],
    [5, "| first_order_potential", " s", " s"],
    [5, "| first_order_H", " s", " s"],
    [5, "| first_order_DM", " s", " s"],
    [2, "rank0, rho kernel and readbuf", ":", "s"],
    [7, " myid=           0  time_work=", " ", ":"],
    [7, " myid=           7  time_work=", " ", ":"],
    [2, "rank0, sumup kernel and readbuf", ":", "s"],
    [2, "rank0, sumup kernel and readbuf", "preproc", "s"],
    [2, "rank0, sumup kernel and readbuf", "main kernel", "s"],
    [3, " myid=           0  sumup-main-time_work=", " ", ":"],
    [3, " myid=           7  sumup-main-time_work=", " ", ":"],
    [3, "rank0, H kernel and readbuf", ":", "s"],
    [8, " myid=           0  time_work=", " ", ":"],
    [8, " myid=           7  time_work=", " ", ":"],
    [3, "| Time for this iteration", " s", " s"],
    [1, "| first_order_density", " s", " s"],
    [1, "| first_order_potential", " s", " s"],
    [1, "| first_order_H", " s", " s"],
    [1, "| first_order_DM", " s", " s"],
    [2, " myid=           0  time_work=", " ", ":"],
    [2, " myid=           7  time_work=", " ", ":"],
    [1, " myid=           0  sumup-main-time_work=", " ", ":"],
    [3, " myid=           0  time_work=", " ", ":"],
    [3, " myid=           7  time_work=", " ", ":"],
    [1, "Time summed over all CPUs for potential: real work", " ", "s"],
    [2, "Time summed over all CPUs for potential: real work", " ", "s"],
    [1, " myid=           0  partition_grid:part1", "=", ":"]

]

search_group_2 = [
    [" myid=           0  sumup-time_work-before-mpi_barrier", "=", ":"],
    [" myid=           0  time of update_hartree_potential_p2_shanghui", "=", ":"]
]

file_list = [
    "./strong_weak_test_all2all/1.5w/sbatch_node4_ntask128_gpu16/slurm-6654765.out",
    "./strong_weak_test_all2all/1.5w/sbatch_node8_ntask256_gpu32/slurm-6654766.out",
    "./strong_weak_test_all2all/1.5w/sbatch_node16_ntask512_gpu64/slurm-6654767.out",
    "./strong_weak_test_all2all/1.5w/sbatch_node32_ntask1024_gpu128/slurm-6654768.out",
    "./strong_weak_test_all2all/1.5w/sbatch_node64_ntask2048_gpu256/slurm-6654769.out",
]

# index must be >= 1, counted from 1


def find_str(fstr, index, substr):
    i = 1
    ret_index = fstr.find(substr)
    tmpstr = fstr[max(0, ret_index-10): min(len(fstr)-1, ret_index+20)]
    while i < index and ret_index != -1:
        ret_index = fstr.find(substr, ret_index+1)
        tmpstr = fstr[max(0, ret_index-10): min(len(fstr)-1, ret_index+20)]
        i = i+1
    return ret_index


def index2linestr(fstr, index):
    end_index = fstr.find("\n", index)
    if end_index == -1:
        return fstr[index:]
    else:
        return fstr[index:end_index]


if __name__ == "__main__":
    all_save_str_c = []
    all_save_str_c2 = []
    for i in range(0, len(search_group_2)):
        all_save_str_c2.append([])
    # filename_in = sys.argv[1]
    # if True:
    for filename_in in file_list:
        (filepath, tempfilename) = os.path.split(filename_in)
        (filename, extension) = os.path.splitext(tempfilename)
        filename_out = os.path.splitext(
            filename_in)[0] + "_info" + os.path.splitext(filename_in)[1]
        print("\n\n=====================================================")
        print("filename_in='{}'\nfilename_out='{}'".format(
            filename_in, filename_out), flush=True)
        with open(filename_in, "r") as f_in:
            str_in = f_in.read()
            save_str = []
            save_str_c = [filename.split('-')[-1]]

            for index, substr, s1, s2 in search_group_1:
                s_index = find_str(str_in, index, substr)
                s_count = str_in.count(substr)
                if s_index == -1:
                    print("s_index==-1, index={}, substr={}".format(index, substr), flush=True)
                    exit(-1)
                    # continued
                s_line = index2linestr(str_in, s_index)
                pstr = "{} : {}".format(s_line, s_count)
                print(pstr, flush=True)
                save_str.append(pstr)
                tmp_str = pstr[len(substr):]
                # print(tmp_str)
                c1 = tmp_str.find(s1)
                c2 = tmp_str.find(s2, c1+1)
                tmp_str2 = tmp_str[c1+len(s1):c2].strip()
                # print(tmp_str2)
                save_str_c.append(tmp_str2)

            sg2_index = -1
            for substr, s1, s2 in search_group_2:
                sg2_index = sg2_index + 1
                count = 1
                s_index = find_str(str_in, 1, substr)

                num_list = []

                while s_index != -1:
                    s_line = index2linestr(str_in, s_index)
                    pstr = "{} : {}".format(s_line, s_count)
                    print(pstr, flush=True)
                    save_str.append(pstr)

                    tmp_str = pstr[len(substr):]
                    c1 = tmp_str.find(s1)
                    c2 = tmp_str.find(s2, c1+1)
                    tmp_str2 = tmp_str[c1+len(s1):c2].strip()
                    tmp_num = float(tmp_str2)
                    num_list.append(tmp_num)

                    count = count + 1
                    s_index = find_str(str_in, count, substr)
                
                statistic_data = "{},{},{}".format(max(num_list), min(num_list), np.mean(num_list))
                all_save_str_c2[sg2_index].append(statistic_data)
                print(sg2_index)
                print(statistic_data)

            all_save_str_c.append(save_str_c)
    for sl in all_save_str_c:
        print(", ".join(sl), flush=True)
    tmp_index = -1
    for tc in all_save_str_c2:
        tmp_index = tmp_index + 1
        print("-------{}------".format(search_group_2[tmp_index][0]))
        print("\n".join(tc), flush=True)
    
