# python .\get_info.py | Tee-Object -FilePath log.get_info.1.txt

import sys
import os
import numpy as np

search_group_1 = [
]

# search_group_2 = [
#     ["update_hartree_potential_p1_without_old:part1", "=", ":"],
#     ["update_hartree_potential_p1_without_old:part2", "=", ":"],
#     ["update_hartree_potential_p1_without_old:part3", "=", ":"],
#     ["update_hartree_potential_p1_without_old:part4", "=", ":"],
#     ["update_hartree_potential_p1_without_old:part5", "=", ":"],
#     ["Hartree multipole update", ":", "s"],
# ]
search_group_2 = [
    ["update_hartree_potential_p1_shanghui:part1", "=", ":"],
    ["update_hartree_potential_p1_shanghui:part2", "=", ":"],
    ["update_hartree_potential_p1_shanghui:part3", "=", ":"],
    ["update_hartree_potential_p1_shanghui:part4", "=", ":"],
    ["update_hartree_potential_p1_shanghui:part5", "=", ":"],
    ["Hartree multipole update", ":", "s"],
]

file_list = [

    # "./strong_weak_test_all2all/3w_atoms_mini_local_index/sbatch_node128_ntask4096_gpu512/slurm-6719064.out",

    # "./strong_weak_test_all2all/6w_atoms_mini_local_index/sbatch_node16_ntask512_gpu64/slurm-6719076.out",

    # basic
    "./3w_atoms_mini_local_index/log.athread.sumup.1680598813425.txt",  # 256
    "./3w_atoms_mini_local_index/log.athread.sumup.1680598813710.txt",  # 512
    "./3w_atoms_mini_local_index/log.athread.sumup.1680598813803.txt",  # 1024
    "./3w_atoms_mini_local_index/log.athread.sumup.1680598813875.txt",  # 2048
    "./3w_atoms_mini_local_index/log.athread.sumup.1680612444660.txt",  # 4096

    "./6w_atoms_mini_local_index/log.athread.sumup.1680605558555.txt",  # 512
    "./6w_atoms_mini_local_index/log.athread.sumup.1680605572714.txt",  # 1024
    "./6w_atoms_mini_local_index/log.athread.sumup.1680605572810.txt",  # 2048
    "./6w_atoms_mini_local_index/log.athread.sumup.1680612601797.txt",  # 4096
    "./6w_atoms_mini_local_index/log.athread.sumup.1680612601935.txt",  # 8192

    # update512
    "./3w_atoms_mini_local_index/log.athread.sumup.1680599641997.txt",  # 256
    "./3w_atoms_mini_local_index/log.athread.sumup.1680599642413.txt",  # 512
    "./3w_atoms_mini_local_index/log.athread.sumup.1680599642506.txt",  # 1024
    "./3w_atoms_mini_local_index/log.athread.sumup.1680599642577.txt",  # 2048
    "./3w_atoms_mini_local_index/log.athread.sumup.1680612487244.txt",  # 4096

    "./6w_atoms_mini_local_index/log.athread.sumup.1680605581450.txt",  # 512
    "./6w_atoms_mini_local_index/log.athread.sumup.1680653250561.txt",  # 1024
    "./6w_atoms_mini_local_index/log.athread.sumup.1680653250683.txt",  # 2048
    "./6w_atoms_mini_local_index/log.athread.sumup.1680612579700.txt",  # 4096
    "./6w_atoms_mini_local_index/log.athread.sumup.1680612580969.txt",  # 8192

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

            sg2_index = -1
            for substr, s1, s2 in search_group_2:
                sg2_index = sg2_index + 1
                count = 1
                s_index = find_str(str_in, 1, substr)

                num_list = []

                while s_index != -1:
                    s_line = index2linestr(str_in, s_index)
                    pstr = "{} : {}".format(s_line, s_index)
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
                
                if len(num_list) != 0:
                    statistic_data = "{},{},{}".format(max(num_list), min(num_list), np.mean(num_list))
                else:
                    statistic_data = "{},{},{}".format(-1, -1, -1)

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
    
