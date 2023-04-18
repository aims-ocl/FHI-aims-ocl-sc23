# -*- coding: utf-8 -*-
import h5py 
from numpy import *

# Python script for converting the hdf5 file containing the Momentummatrix into 
# an ascii file. mommat.h5-->mommat.dat
# Be aware of filesizes!

# numpy needed (available in ubuntu repositories)
# h5py needed (available in ubuntu repositories)

f=h5py.File('mommat.h5', 'r')
list_of_names = []
f.visit(list_of_names.append)
ind_mom=list_of_names.index('Momentummatrix')
dim=int32(asarray(str(f.listitems()[ind_mom][1]).strip().split('(')[1].split(')')[0].split(',')))
dip=f.require_dataset(f.listnames()[ind_mom],tuple(dim),float)
dip=dip[...]
ind_k=list_of_names.index('k_points')
dim_k=int32(asarray(str(f.listitems()[ind_k][1]).strip().split('(')[1].split(')')[0].split(',')))
k_points=f.require_dataset(f.listnames()[ind_k],tuple(dim_k),float)
k_points=k_points[...]
ind_win=list_of_names.index('Energy_window')
dim_win=int32(asarray(str(f.listitems()[ind_win][1]).strip().split('(')[1].split(')')[0].split(',')))
win=f.require_dataset(f.listnames()[ind_win],tuple(dim_win),float)
win=win[...]
f.close()
write_file=open('mommat.dat','w')
write_file.writelines('# column 1: i state\n')
write_file.writelines('# column 2: j state\n')
write_file.writelines('# column 3: Re(<\\varphi_i|\\nabla_x|\\varphi_j>)\n')
write_file.writelines('# column 4: Im(<\\varphi_i|\\nabla_x|\\varphi_j>)\n')
write_file.writelines('# column 5: Re(<\\varphi_i|\\nabla_y|\\varphi_j>)\n')
write_file.writelines('# column 6: Im(<\\varphi_i|\\nabla_y|\\varphi_j>)\n')
write_file.writelines('# column 7: Re(<\\varphi_i|\\nabla_z|\\varphi_j>)\n')
write_file.writelines('# column 8: Im(<\\varphi_i|\\nabla_z|\\varphi_j>)\n')
write_file.writelines('#\n')
for k in range(int(k_points[:,:,:,0].max())):
  num=0   
  ind_k=k_points[:,:,:,0]==k
  write_file.writelines('# k_point '+str(int(k+1))+' at '+str(float(k_points[ind_k,1]))+
                       ' '+str(float(k_points[ind_k,2]))+' '+str(float(k_points[ind_k,3]))+'\n')
  for l in arange(win[0],win[0]+win[1]-win[0]+1):
    for m in arange(l, win[0]+win[1]-win[0]+1):      
      write_file.writelines(str(int(l))+' '+str(int(m))+' '+str(dip[ind_k,num,:][0][0])+
                        ' '+str(dip[ind_k,num,:][0][1])+' '+str(dip[ind_k,num,:][0][2])+
                        ' '+str(dip[ind_k,num,:][0][3])+' '+str(dip[ind_k,num,:][0][4])+
                        ' '+str(dip[ind_k,num,:][0][5])+'\n')
      num = num+1
write_file.close()
