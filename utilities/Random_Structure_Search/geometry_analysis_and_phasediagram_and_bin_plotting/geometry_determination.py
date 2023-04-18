# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 10:28:21 2014

@author: jank
"""

#from subprocess import call, Popen, PIPE
import os
import re
import string
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpat
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator, NullLocator, AutoLocator
import matplotlib.image as mpimg
from matplotlib._png import read_png
from pylab import *

class GeoDeterm:
  def get_line_nums(self,filename):
    with open(filename) as f:
      for i, l in enumerate(f):
        pass
    return i + 1

  def readData(self):
    f=open("end_geo_files","r")
    self.endfiles = f.readlines()
    f.close()
    f=open("all-lowest_tot_en.dat")
    self.lowest_tot_en = f.readlines()
    f.close()
    lowest_tot_file = self.lowest_tot_en[0].split(" ")[4]
    self.lowest_tot_en = np.float(self.lowest_tot_en[0].split(" ")[2])
    self.lowest_tot_file = re.sub(r'run.light', 'end_geo.in', lowest_tot_file).strip()

  def makeDotPlot(self):
    plt.rc('text', usetex=True)
    plt.rc('font', family="serif", serif="Palatino")
    template_color = 'black'
    f=open("all-converged.dat")
    content = f.readlines()
    f.close()

    print("\n  Generating DotPlot first...\n")

    x = [line.split('   ')[0] for line in content]
    y = [line.split('   ')[1] for line in content]
    structs = [line.split('    '[2]) for line in content]

    strnum = ""
    gmin = 0
    gmax = -10000000
    k=0

    for i in y:
      k=k+1
      if (float(i)<gmin):
        gmin = float(i)
        strnum = structs[k]
      if (float(i)>gmax):
        gmax = float(i)

    gmax = gmax - gmin

    for i in range(0,k):
      x[i] = i+1
      y[i] = float(y[i]) - gmin
      if (float(y[i]==0)):
       xann = i+1

    for i in range(10):
      if(k<=pow(10,i)):
        xlrng = pow(10,i-1)
        break

    for i in range(10):
      if(gmax<=pow(10,i)):
        ylrng = pow(10,i-1)
        break

    locymin = 100
    for i in range(k):
      if(y[i]>0.5):
        if(y[i]<locymin):
          locxmin = i+1
          locymin = y[i]

    locymax = 0
    for i in range(k):
      if(y[i]<0.5):
        if(y[i]>locymax):
          locxmax = i+1
          locymax = y[i]

    figa = plt.figure()

    ax1a = figa.add_subplot(111)

    majxloc = MultipleLocator(xlrng)
    majxfmt = FormatStrFormatter('%d')
    minxloc = AutoMinorLocator(10)

    majyloc = MultipleLocator(ylrng)
    majyfmt = FormatStrFormatter('%d')
    minyloc = AutoMinorLocator(10) #MultipleLocator(10)

    ax1a.spines['bottom'].set_color(template_color)
    ax1a.spines['top'].set_color(template_color) 
    ax1a.spines['right'].set_color(template_color)
    ax1a.spines['left'].set_color(template_color)

    ax1a.tick_params(axis='x', colors=template_color)
    ax1a.tick_params(axis='y', colors=template_color)

    ax1a.yaxis.label.set_color(template_color)
    ax1a.xaxis.label.set_color(template_color)

    ax1a.tick_params(axis='both', which='major', color=template_color, labelsize=14)

    plt.title(r"{\Large Random Structure Search Results} \newline {\hspace*{6mm}\large SiC $000\bar{1}$ surface Si-substructure search}", fontsize=18, color=template_color)

    ax1a.set_xlabel(r"Structure Number", color=template_color, fontsize=16)
    ax1a.set_ylabel(r"$\Delta$E$_{tot}$ [eV]", color=template_color, fontsize=16)

    plt.xlim((0,k+1))
    plt.ylim((-0.45,gmax+0.5))

    ax1a.plot(locxmin, locymin, "o", color="k", markerfacecolor="w", markersize=10)
    ax1a.plot(locxmax, locymax, "o", color="k", markerfacecolor="w", markersize=10)

    ax1a.plot(xann, 0,"8", color="k", markerfacecolor="w", markersize=13)
    ax1a.plot(x, y, "bo", markersize=2, label=r"$\Delta$E$_{tot}$")

    ax1a.xaxis.set_major_locator(majxloc)
    ax1a.xaxis.set_major_formatter(majxfmt)
    ax1a.xaxis.set_minor_locator(minxloc)

    ax1a.yaxis.set_major_locator(majyloc)
    ax1a.yaxis.set_major_formatter(majyfmt)
    ax1a.yaxis.set_minor_locator(minyloc)

    plt.annotate("", xy=(locxmin, locymin), xycoords='data',
                 xytext=(locxmax, locymax), textcoords='data',
                 arrowprops=dict(arrowstyle="<->", linestyle="dashed",
                 color="0.7", #grey level or color name
                 shrinkA=9, shrinkB=9,
                 patchA=None, patchB=None,
                 connectionstyle="arc3, rad=0.1", #kind of arrow
                 ))

    locediff = abs(locymax-locymin)

    outstr = r"$\Delta\epsilon$="+repr(int(locediff*1000)/1000)+"eV"

    plt.text( locxmin+(locxmax-locxmin)/2,
              locymin+(locymax-locymin)/2,
              outstr, fontsize=10)

    annstr = r"{\large E$_{min}^{"+strnum[0]+"}$ = \small "+repr(int(gmin*100)/100)+"eV}"

    plt.annotate(annstr, xy=(xann,0), xytext=(xann-3,-0.4), fontsize=10)

    plt.savefig('dots.png', dpi=170)

  def processData(self,geo_thresh):
    total_end_geos = 0
    lnum = self.get_line_nums(self.endfiles[1].split("    ")[0])
    for line in self.endfiles:
      total_end_geos += 1

    anum = lnum - 3
    self.num_atoms = anum

    counted_si = 0
    counted_c = 0
    counted_h = 0
    lj=0

    name = self.endfiles[0].split("    ")[0]
    f = open(name,"r")
    content = f.readlines()
    f.close()
    for cont in content:
      cont = cont.split("%")[0]
      entrs = cont.split()
      if(len(entrs)==0):
        continue
      if(len(entrs)==4):
        lj+=1
      if(len(entrs)==5):
        if(entrs[4]=='Si'):
          counted_si += 1
        if(entrs[4]=='C'):
          counted_c += 1
        if(entrs[4]=='H'):
          counted_h += 1

    if(lj!=0 and lj!=3):
      print("       *****      ERROR      *****")
      print("        Geometry File Format Error")
      print("  Lattice vectors must be 3 or not present")
      quit()

    geo_code = ["Null" for x in range(0,total_end_geos)]
    for i in range(0,total_end_geos):
      geo_code[i] = self.endfiles[i].split("    ")[0]

    print("\n   I found %1s geometry files.\n" % total_end_geos) # with %1s lines each. (%1s lattice vectors)\n" % (total_end_geos,lnum,lj))
    print("      In the geometry files I found %3.1i Si, %3.1i C and %2.1i H atoms." % (counted_si, counted_c, counted_h))
    print("      Calculation done by +- 1 cell periodic expansion within x-y plane.")
    print("      Resulting distances are %1d in each file yielding in %1d total." % ( 5*anum*(anum-1)/2, 5*anum*(anum-1)*total_end_geos/2 ) )
    print("      Threshold for geometry difference deviations is %3s Angstrøm [1Å = 0.1nm]." % geo_thresh )
    print("      Overall total lowest energy found is %1s eV" % self.lowest_tot_en)

    lat_data = np.empty((3,3))
    lat_data.fill(np.nan)
    geo_data = np.empty((anum,3))
    geo_data.fill(np.nan)

    num_dists = int(5 * anum * (anum -1)/2)
    dist_tab = np.empty((total_end_geos,num_dists))
    dist_tab.fill(np.nan)

    en_tot = np.empty((total_end_geos))
    en_tot.fill(np.nan)

    print("\n\n   Reading files and computing geometry differences now...")

    k=0
    for line in self.endfiles:
      line = line.split("%")[0]
      vals = line.split()
      en_tot[k] = vals[1]
      k+=1

    k=0
    for line in self.endfiles:
      name = line.split("    ")[0]
      f = open(name,"r")
      content = f.readlines()
      f.close()
      m=0
      lj = 0
      gj = 0
      for cont in content:
        cont = cont.split("%")[0]
        entrs = cont.split()
        if(len(entrs)==0):
          continue
        if(len(entrs)==4):
          lat_data[lj][0] = entrs[1]
          lat_data[lj][1] = entrs[2]
          lat_data[lj][2] = entrs[3]
          lj += 1
        if(len(entrs)==5):
          geo_data[gj][0] = entrs[1]
          geo_data[gj][1] = entrs[2]
          geo_data[gj][2] = entrs[3]
          gj += 1

      for i in range(0,anum):
        for j in range(i+1,anum):
          dist_tab[k][m] = np.sqrt( \
               np.power(geo_data[i][0] - geo_data[j][0],2) +
               np.power(geo_data[i][1] - geo_data[j][1],2) +
               np.power(geo_data[i][2] - geo_data[j][2],2) )
          m += 1
          dist_tab[k][m] = np.sqrt( \
               np.power(geo_data[i][0] - geo_data[j][0]+lat_data[0][0],2) +
               np.power(geo_data[i][1] - geo_data[j][1]+lat_data[0][1],2) +
               np.power(geo_data[i][2] - geo_data[j][2],2) )
          m += 1
          dist_tab[k][m] = np.sqrt( \
               np.power(geo_data[i][0] - geo_data[j][0]-lat_data[0][0],2) +
               np.power(geo_data[i][1] - geo_data[j][1]-lat_data[0][1],2) +
               np.power(geo_data[i][2] - geo_data[j][2],2) )
          m += 1
          dist_tab[k][m] = np.sqrt( \
               np.power(geo_data[i][0] - geo_data[j][0]+lat_data[1][0],2) +
               np.power(geo_data[i][1] - geo_data[j][1]+lat_data[1][1],2) +
               np.power(geo_data[i][2] - geo_data[j][2],2) )
          m += 1
          dist_tab[k][m] = np.sqrt( \
               np.power(geo_data[i][0] - geo_data[j][0]-lat_data[1][0],2) +
               np.power(geo_data[i][1] - geo_data[j][1]-lat_data[1][1],2) +
               np.power(geo_data[i][2] - geo_data[j][2],2) )
          m += 1
      k+=1

    area = self.getSurfaceArea(lat_data)

    del geo_data

    hita = np.empty((20*total_end_geos))
    hita.fill(np.nan)
    hitb = np.empty((20*total_end_geos))
    hitb.fill(np.nan)

    for k in range(0,total_end_geos):
      dist_tab[k].sort()

    print("\n   Comparing all geometry differences to each other...")

    m = 0
    for k in range(0,total_end_geos):
      if(np.mod(k,10)==0):
        print("    - Working... (k = %3s)" % k)
      for l in range(k+1,total_end_geos):
        identical = 1
        while (identical==1):
          for i in range(0,num_dists):
            if( (dist_tab[k][i]!="nan") and (dist_tab[l][i]!="nan") ):
              if(np.abs(dist_tab[k][i]-dist_tab[l][i])>=geo_thresh):
                identical = 0
                break
          if(identical==1):
            hita[m] = k
            hitb[m] = l
            m+=1
            print("     ☺ Hit at %1d %1d" % (k, l))
            identical = 0

    f=open('matchtes','w')
    print("\n   After that I would state that:")
    for i in range(0,m):
      if (hita[i] != "nan"):
        if(hita[i]!=hitb[i]):
          print("   - Structures %1s and %1s seem to be matching with dE= %1s meV" %
            ( geo_code[int(hita[i])], geo_code[int(hitb[i])], np.abs(en_tot[int(hita[i])]-en_tot[int(hitb[i])])*float(100.0) ) )
          outstr = geo_code[int(hita[i])]+"   "+geo_code[int(hitb[i])]+"\n"
          f.write(outstr)
    f.close()

    self.all_hits = [["Null" for i in range(3)] for j in range(2*m)]
    k=0
    for i in range(m):
      known_eq = 0
      known_a = 0
      known_b = 0
      if(hita[i]==hitb[i]):
        for j in range(2*m):
          if(self.all_hits[j][0]==geo_code[int(hita[i])]):
            known_eq = 1
        if(known_eq==0):
          self.all_hits[k][0] = geo_code[int(hita[i])]
          self.all_hits[k][1] = en_tot[int(hita[i])]
          self.all_hits[k][2] = en_tot[int(hita[i])]-self.lowest_tot_en
          k+=1
      else:
        for j in range(2*m):
          if(self.all_hits[j][0]==geo_code[int(hita[i])]):
            known_a = 1
          if(self.all_hits[j][0]==geo_code[int(hitb[i])]):
            known_b = 1
        if(known_a==0):
          self.all_hits[k][0] = geo_code[int(hita[i])]
          self.all_hits[k][1] = en_tot[int(hita[i])]
          self.all_hits[k][2] = en_tot[int(hita[i])]-self.lowest_tot_en
          k+=1
        if(known_b==0):
          self.all_hits[k][0] = geo_code[int(hitb[i])]
          self.all_hits[k][1] = en_tot[int(hitb[i])]
          self.all_hits[k][2] = en_tot[int(hitb[i])]-self.lowest_tot_en
          k+=1

    self.num_hits = k

    return area, self.all_hits, self.lowest_tot_en, counted_si, counted_c, counted_h

  def getSurfaceArea(self,lat_data):
    area =   lat_data[0][1]*lat_data[1][2]-lat_data[0][2]*lat_data[1][1] \
           + lat_data[0][2]*lat_data[1][0]-lat_data[0][0]*lat_data[1][2] \
           + lat_data[0][0]*lat_data[1][1]-lat_data[1][0]*lat_data[0][1]
    return area

  def makePlot(self,num_bins):
    plt.rc('text', usetex=True)
    plt.rc('font', family="serif", serif="Palatino")

    min_de = 1000000
    max_de = 0

    for i in range(self.num_hits):
      if(min_de>self.all_hits[i][2]):
        min_de = self.all_hits[i][2]
      if(max_de<self.all_hits[i][2]):
        max_de = self.all_hits[i][2]

    plot_range = max_de-min_de
    bin_width = np.float(plot_range)/np.float(num_bins)

    x_index = np.empty((num_bins))
    x_index.fill(np.nan)

    bins = np.empty((num_bins))
    bins.fill(0)

    structure_file = ['' for i in range(num_bins)]

    for i in range(num_bins):
      x_index[i] = i+1
      for k in range(self.num_hits):
        if(self.all_hits[k][0]!="Null"):
          if(np.float(self.all_hits[k][2])>=i*bin_width+min_de and np.float(self.all_hits[k][2])<(i+1)*bin_width+min_de):
            bins[i] += 1
            structure_file[i] += self.all_hits[k][0]+'   '

    for i in range(10):
     if(num_bins<=pow(10,i)):
      x_scale_over = 0.3 * pow(10,i-1)
      x_ticks_div_ma = pow(10,i-1)
      x_ticks_div_mi = pow(10,i-2)
      break

    barwdth = 0.15

    fig = plt.figure()
    it = fig.add_subplot(111)
    plt.title(r"{\huge Canned Results}\\ {\centering\Large Random Structure Search}")
    plt.xlabel(r"Structure Bin {\large ["+("{0:.3g}".format(bin_width))+" meV per Bin]}", fontsize=16)
    plt.ylabel(r"$\sum$ Structures {\large [$\sharp$ per Bin]}", fontsize=16)  
    plt.xlim((1-x_scale_over,np.amax(x_index)+x_scale_over))
    plt.ylim((0.0,np.amax(bins)+0.5))

    annstr = r"{\Large E}$_{min}^{tot}$= "+repr((self.lowest_tot_en*1000)/1000)+" eV"

    it.text(0.6, 0.9, annstr, verticalalignment='bottom', horizontalalignment='left', transform=it.transAxes, color='black', fontsize=14)

    it.xaxis.set_major_locator(MultipleLocator(x_ticks_div_ma))
    it.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    it.xaxis.set_minor_locator(MultipleLocator(x_ticks_div_mi))
    it.yaxis.set_major_locator(MultipleLocator(1))
    it.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    it.yaxis.set_minor_locator(NullLocator())

    self.RenderSlabPicture()

    call = 'convert it.png -transparent white slab_geo.png'
    os.system(call)
    call = 'rm -f it.png pic.* go2d* input'
    os.system(call)

    slab_geo = read_png('slab_geo.png')
    xmin = 0 - 0.15 * num_bins
    xmax = 0 + 0.35 * num_bins
    ymin = np.amax(bins) - 0.55 * np.amax(bins)
    ymax = np.amax(bins) - 0.05 * np.amax(bins)
    it.imshow(slab_geo, extent=[xmin,xmax,ymin,ymax], origin='upper', alpha=1, aspect='auto', interpolation='none')

    plt.bar(x_index-barwdth/2, bins, barwdth, color='blue', edgecolor='blue')
    plt.tight_layout()
    plt.savefig('bins.png',dpi=170)

  def RenderSlabPicture(self):
    mged_output_resolution = "-w 1280 -n 1024"

    f = open(self.lowest_tot_file,'r')
    content = f.readlines()
    f.close()

    print ("\n   We have ",self.num_atoms," atoms in the input file.\n")

    vals = [[0 for x in range(7)] for x in range(self.num_atoms)]

    h_sph = [0 for x in range(self.num_atoms)]
    c_sph = [0 for x in range(self.num_atoms)]
    n_sph = [0 for x in range(self.num_atoms)]
    o_sph = [0 for x in range(self.num_atoms)]
    au_sph = [0 for x in range(self.num_atoms)]
    si_sph = [0 for x in range(self.num_atoms)]

    i=0
    for cont in content:
      cont = cont.split("%")[0]
      entrs = cont.split()
      if(entrs[0]=='lattice_vector'):
        continue
      if(entrs[0]=='constrain_relaxation'):
        continue
      vals[i][1] = np.float(entrs[1])
      vals[i][2] = np.float(entrs[2])
      vals[i][3] = np.float(entrs[3])

      if entrs[4] == "H":
        vals[i][0] = "Hydrogen"
        vals[i][4] = 0.25 #radius of shpere
      if entrs[4] == "C":
        vals[i][0] = "Carbon"
        vals[i][4] = 0.9 #radius of shpere
      if entrs[4] == "O":
        vals[i][0] = "Oxygen"
        vals[i][4] = 0.6 #radius of shpere
      if entrs[4] == "N":
        vals[i][0] = "Nitrogen"
        vals[i][4] = 0.65 #radius of shpere
      if entrs[4] == "Au":
        vals[i][0] = "Aurum"
        vals[i][4] = 1.5 #radius of shpere
      if entrs[4] == "Si":
        vals[i][0] = "Silicon"
        vals[i][4] = 1.1 #radius of shpere

      i += 1

    h_present = 0
    c_present = 0
    n_present = 0
    o_present = 0
    au_present = 0
    si_present = 0

    sphere_cnt=0
    y_cnt=0

    call = 'rm -f it.png pic.g go2d* input'
    os.system(call)
    f = open('input','w')
    f.write("units mm\n")
    for i in range(0,self.num_atoms):
      if vals[i][0] == "Hydrogen":
        h_present = 1
        f.write('%-5s%-1i%-8s %14s %14s %14s %14s\n' % ('in sp',y_cnt,'.s sph',vals[i][1],vals[i][2],vals[i][3],vals[i][4]) )
        h_sph[i] = y_cnt+1
        y_cnt += 1
      if vals[i][0] == "Carbon":
        c_present = 1
        f.write('%-5s%-1i%-8s %14s %14s %14s %14s\n' % ('in sp',y_cnt,'.s sph',vals[i][1],vals[i][2],vals[i][3],vals[i][4]) )
        c_sph[i] = y_cnt+1
        y_cnt += 1
      if vals[i][0] == "Nitrogen":
        n_present = 1
        f.write('%-5s%-1i%-8s %14s %14s %14s %14s\n' % ('in sp',y_cnt,'.s sph',vals[i][1],vals[i][2],vals[i][3],vals[i][4]) )
        n_sph[i] = y_cnt+1
        y_cnt += 1
      if vals[i][0] == "Oxygen":
        o_present = 1
        f.write('%-5s%-1i%-8s %14s %14s %14s %14s\n' % ('in sp',y_cnt,'.s sph',vals[i][1],vals[i][2],vals[i][3],vals[i][4]) )
        o_sph[i] = y_cnt+1
        y_cnt += 1
      if vals[i][0] == "Aurum":
        au_present = 1
        f.write('%-5s%-1i%-8s %14s %14s %14s %14s\n' % ('in sp',y_cnt,'.s sph',vals[i][1],vals[i][2],vals[i][3],vals[i][4]) )
        au_sph[i] = y_cnt+1
        y_cnt += 1
      if vals[i][0] == "Silicon":
        si_present = 1
        f.write('%-5s%-1i%-8s %14s %14s %14s %14s\n' % ('in sp',y_cnt,'.s sph',vals[i][1],vals[i][2],vals[i][3],vals[i][4]) )
        si_sph[i] = y_cnt+1
        y_cnt += 1

    if h_present == 1: 
      h_region = "r region1.r"
      for i in range(0,self.num_atoms):
        if h_sph[i] != 0:
          h_region += " u sp%s.s" % (h_sph[i]-1)
      h_region += "\n"
      f.write(h_region)
      f.write("draw region1.r\n")
      f.write("comb h_net u region1.r\n")
      f.write("mater h_net \"plastic tr=0.0 re=0.1 sp=0.8 di=0.3 sh=7 ex=0.0 ri=1.0\" 155 155 255 1\n")

    if c_present == 1:
      c_region = "r region2.r"
      for i in range(0,self.num_atoms):
        if c_sph[i] != 0:
          c_region += " u sp%s.s" % (c_sph[i]-1)
      c_region += "\n"
      f.write(c_region)
      f.write("draw region2.r\n")
      f.write("comb c_net u region2.r\n")
      f.write("mater c_net \"plastic tr=0.0 re=0.25 sp=0.45 di=0.25 sh=6 ex=0.0 ri=1.65 em=0.0\" 125 125 125 1\n")

    if n_present == 1:
      n_region = "r region3.r"
      for i in range(0,self.num_atoms):
        if n_sph[i] != 0:
          n_region += " u sp%s.s" % (n_sph[i]-1)
      n_region += "\n"
      f.write(n_region)
      f.write("draw region3.r\n")
      f.write("comb n_net u region3.r\n")
      f.write("mater n_net \"plastic tr=0.0 re=0.1 sp=0.8 di=0.3 sh=7 ex=0.0 ri=1.0\" 010 015 175 1\n")

    if o_present == 1:
      o_region = "r region4.r"
      for i in range(0,self.num_atoms):
        if o_sph[i] != 0:
          o_region += " u sp%s.s" % (o_sph[i]-1)
      o_region += "\n"
      f.write(o_region)
      f.write("draw region4.r\n")
      f.write("comb o_net u region4.r\n")
      f.write("mater o_net \"plastic tr=0.0 re=0.1 sp=0.8 di=0.3 sh=7 ex=0.0 ri=1.0\" 220 010 010 1\n")

    if au_present == 1:
      au_region = "r region5.r"
      for i in range(0,self.num_atoms):
        if au_sph[i] != 0:
          au_region += " u sp%s.s" % (au_sph[i]-1)
      au_region += "\n"
      f.write(au_region)
      f.write("draw region5.r\n")
      f.write("comb au_net u region5.r\n")
      f.write("mater au_net \"plastic tr=0.0 re=0.8 sp=0.8 di=0.3 sh=10 ex=0.0 ri=0.8 em=0.8\" 074 187 050 1\n")

    if si_present == 1:
      si_region = "r region6.r"
      for i in range(0,self.num_atoms):
        if si_sph[i] != 0:
          si_region += " u sp%s.s" % (si_sph[i]-1)
      si_region += "\n"
      f.write(si_region)
      f.write("draw region6.r\n")
      f.write("comb si_net u region6.r\n")
      f.write("mater si_net \"plastic tr=0.0 re=0.6 sp=0.4 di=0.7 sh=05 ex=0.0 ri=0.8 em=0.8\" 142 035 035 1\n")

    f.write("in sp1000000.s sph 25.0 18.0 -23.0 1.0\n")
    f.write("r lights1.r u sp1000000.s\n")
    f.write("draw lights1.r\n")
    f.write("comb light u lights1.r\n")
    f.write("mater lights1.r \"light f=3.0\" 225 225 255 0\n")

    f.write("Z\n")

    if h_present == 1:
      f.write("draw h_net\n")
    if c_present == 1:
      f.write("draw c_net\n")
    if n_present == 1:
      f.write("draw n_net\n")
    if o_present == 1:
      f.write("draw o_net\n")
    if au_present == 1:
      f.write("draw au_net\n")
    if si_present == 1:
      f.write("draw si_net\n")

    f.write("draw light\n")

    xang = -100
    yang = 0
    zang = 0

    f.write("rot %5.3f %5.3f %5.3f\n" % (xang, yang, zang) )
    f.write("autoview\n")
    f.write("zoom 2.5\n")
    f.write("eye_pt 0.0 0.0 15.0\n")
    f.write("saveview go2d %1s\n" % (mged_output_resolution) )
    f.close()

    call = 'rm -f it.png pic.* go2d*'
    os.system(call)
    call = 'cat input | mged -c pic.g'
    os.system(call)
    call = 'sed -i s/"rt -M"/"rt -M -P 2 -C255\/255\/255"/ go2d'
    os.system(call)
    call = 'chmod 700 go2d'
    os.system(call)
    call = './go2d'
    os.system(call)
    call = 'pix-png '+mged_output_resolution+' -o it.png go2d.pix '
    os.system(call)

class PhaseDiagram:
  def readData(self):
    f=open("../1x1_clean_geometry")
    self.clean_1x1_geo = f.readlines()
    f.close()
    f=open("../1x1_clean_energy")
    self.eslab_clean = np.float(f.readline())
    f.close
    f=open("all-lowest_tot_en.dat")
    self.lowest_tot_en = f.readlines()
    self.lowest = np.float(self.lowest_tot_en[0].split(" ")[2])
    f.close()
    f=open("1/1/geometry.in")
    self.firstgeo = f.readlines()
    f.close()

  def getSurfaceArea(self,lat_data):
    area =   lat_data[0][1]*lat_data[1][2]-lat_data[0][2]*lat_data[1][1] \
           + lat_data[0][2]*lat_data[1][0]-lat_data[0][0]*lat_data[1][2] \
           + lat_data[0][0]*lat_data[1][1]-lat_data[1][0]*lat_data[0][1]
    return area

  def setRefVals(self):
    plt.rc('text', usetex=True)
    plt.rc('font', family="serif", serif="Palatino")

    self.esic = -17824.738061031 * 0.5
    self.uc_grit = -1036.55209
    self.uc_dmnd = -1036.496202
    self.uc_sic = -1037.0521
    self.usi = -7875.3233
    self.uh = -15.874026543

    self.grit_line_x = [(self.uc_grit-self.uc_grit) for i in range(-10,11)]
    self.grit_line_y = [i for i in range(-10,11)]

    self.dmnd_line_x = [(self.uc_dmnd-self.uc_grit) for i in range(-10,11)]
    self.dmnd_line_y = [i for i in range(-10,11)]

    self.si_line_x = [(self.uc_sic-self.uc_grit) for i in range(-10,11)]
    self.si_line_y = [i for i in range(-10,11)]

    self.zero_line_x = [i for i in range(-1,2)]
    self.zero_line_y = [0 for i in range(-1,2)]

  def setPlotParams(self):
    fig = plt.figure()

    self.ax1 = fig.add_subplot(111)
    plt.title(r"{\Large Phase Diagram for the SiC-$000\bar{1}$ Surface Energy} \newline {\large $~~~~~~~~~~~~~~~~~~~~~~\gamma_C=\frac{1}{\textrm{\small A}}\left(\right.$E$_{\textrm{\small slab}}-$N$_{\textrm{\footnotesize Si}}\mu_{\textrm{\footnotesize Si}}-$N$_{\textrm{\footnotesize C}}\mu_{\textrm{\footnotesize C}}-$N$\left._{\textrm{\footnotesize H}}\mu_{\textrm{\footnotesize H}}\right)$}",
             fontsize=6, horizontalalignment='center')

    self.ax1.set_xlabel(r"$\mu_\textrm{\small C}-$E$_\textrm{\small C}^{\textrm{\small 1$\times$1}}$ [eV]", fontsize=16)
    self.ax1.set_ylabel(r"$\gamma_C$ [eV/SiC$_{1\times1}$]", fontsize=16)
    self.ax1.tick_params(axis='both', which='major', labelsize=14)

    majxloc = MultipleLocator(0.1)
    majxfmt = FormatStrFormatter('%2.1f')
    minxloc = AutoMinorLocator(10)

    majyloc = MultipleLocator(0.5)
    majyfmt = FormatStrFormatter('%1.1f')
    minyloc = AutoMinorLocator(10)

    self.ax1.xaxis.set_major_locator(majxloc)
    self.ax1.xaxis.set_major_formatter(majxfmt)
    self.ax1.xaxis.set_minor_locator(minxloc)

    self.ax1.yaxis.set_major_locator(majyloc)
    self.ax1.yaxis.set_major_formatter(majyfmt)
    self.ax1.yaxis.set_minor_locator(minyloc)

    plt.xlim((-0.57,0.11))
    plt.ylim((-1.23,1.23))

    self.myaspect = (0.11+0.57)/(1.23+1.23)

  def PlotRefVals(self):
    self.ax1.plot(self.grit_line_x, self.grit_line_y, ":", color="k", linewidth=0.8)
    plt.text( 0.005, 0.85, r"Graphite", rotation=90, fontsize=12)

    self.ax1.plot(self.dmnd_line_x, self.dmnd_line_y, ":", color="k", linewidth=0.8)
    plt.text( self.uc_dmnd-self.uc_grit+0.005, 0.85, r"Diamond", rotation=90, fontsize=12)

    self.ax1.plot(self.si_line_x, self.dmnd_line_y, ":", color="k", linewidth=0.8)
    plt.text( self.uc_sic-self.uc_grit-0.015, 0.85, r"SiC Bulk", rotation=90, fontsize=12)

    x = np.arange(-1.0, 1.0, 0.001)
    self.ax1.fill_between(x, -2, 2, where = x >= self.uc_dmnd-self.uc_grit, facecolor='black', alpha=0.10, interpolate=False)
    self.ax1.fill_between(x, -2, 2, where = x >= self.uc_grit-self.uc_grit, facecolor='black', alpha=0.05, interpolate=False)
    self.ax1.fill_between(x, -2, 2, where = x <= self.uc_sic -self.uc_grit, facecolor='black', alpha=0.15, interpolate=False)

    self.ax1.plot(self.zero_line_x, self.zero_line_y, "-", color="k", linewidth=0.75)
    plt.text( -0.49, -0.08, r"$1\times 1_\textrm{\footnotesize 0Si$\otimes$0C}$", rotation=0, fontsize=10)

  def gamma_c(self,i,n,a1x1,esl1,esl2,nsi1,nsi2,nc1,nc2,nh1,nh2,esic,uc_grit):
    a = 1/a1x1 * ( esl1 - nsi1 * (esic-self.uc_grit) - nc1 * (self.uc_grit) - nh1 * self.uh )
    b = 1/(n*a1x1) * ( esl2 - nsi2 * (esic-(self.uc_grit+i)) - nc2 * (self.uc_grit + i) - nh2 * self.uh )
    result = (b-a)*a1x1
    return result

  def PlotCalcVals(self, area, emin, counted_si, counted_c, counted_h):
    self.readData()

    data = np.empty((3,3))
    data.fill(np.nan)

    k=0
    for cont in self.clean_1x1_geo:
      cont = cont.split("%")[0]
      entrs = cont.split()
      if(len(entrs)==0):
        continue
      if(len(entrs)==4):
        data[k][0] = entrs[1]
        data[k][1] = entrs[2]
        data[k][2] = entrs[3]
        k+=1

    data_geo = np.empty((3,3))
    data_geo.fill(np.nan)

    k=0
    for cont in self.firstgeo:
      cont = cont.split("%")[0]
      entrs = cont.split()
      if(len(entrs)==0):
        continue
      if(len(entrs)==4):
        data_geo[k][0] = entrs[1]
        data_geo[k][1] = entrs[2]
        data_geo[k][2] = entrs[3]
        k+=1

    xax = np.empty((3))
    xax.fill(np.nan)
    xax = [ i for i in range(-1,2)]

    #############################
    # Define the reference slab #
    #############################
    a1x1 = self.getSurfaceArea(data)
    esl1 = self.eslab_clean
    nsi1 = 4
    nc1 = 4
    nh1 = 1

    #######################
    # SiC 2x2 1AdAtoms Si #
    #######################
    esl2 = -150528.780565918
    n = 4
    nsi2 = n*4 + 1
    nc2 = n*4
    nh2 = n*1

    g_line1 = np.empty((3))
    g_line1.fill(np.nan)
    g_line1 = [self.gamma_c(i,n,a1x1,esl1,esl2,nsi1,nsi2,nc1,nc2,nh1,nh2,self.esic,self.uc_grit) for i in range(-1,2)]
    self.ax1.plot(xax,g_line1, "-", color="k", linewidth=1)
    x_txt = -0.4
    y_txt = self.gamma_c(x_txt,n,a1x1,esl1,esl2,nsi1,nsi2,nc1,nc2,nh1,nh2,self.esic,self.uc_grit)
    ang = np.arctan(g_line1[1]-g_line1[0])*180/np.pi*self.myaspect
    plt.text( x_txt, y_txt-0.06, r"$2\times 2_{\textrm{\footnotesize 1Si$\otimes$0C}}$", rotation=ang, fontsize=10)

    ########################
    # SiC 3x3 13AdAtoms Si #
    ########################
    esl2 = -583777.283232368
    n = 9
    nsi2 = n*6 + 13
    nc2 = n*6
    nh2 = n*1

    g_line2 = np.empty((3))
    g_line2.fill(np.nan)
    g_line2 = [self.gamma_c(i,n,a1x1,esl1,esl2,nsi1,nsi2,nc1,nc2,nh1,nh2,self.esic,self.uc_grit) for i in range(-1,2)]
    self.ax1.plot(xax,g_line2, "-", color="k", linewidth=1)
    x_txt = -0.28
    y_txt = self.gamma_c(x_txt,n,a1x1,esl1,esl2,nsi1,nsi2,nc1,nc2,nh1,nh2,self.esic,self.uc_grit)
    ang = np.arctan(g_line2[1]-g_line2[0])*180/np.pi*self.myaspect
    plt.text( x_txt, y_txt+0.01, r"$3\times 3_{\textrm{\footnotesize 13Si$\otimes$0C}}$", rotation=ang, fontsize=10)

    ########################
    # SiC Calculation Data #
    ########################
    esl2 = self.lowest
    area = self.getSurfaceArea(data_geo)
    n = area/a1x1
    nsi2 = counted_si
    nc2 = counted_c
    nh2 = counted_h

    g_line3 = np.empty((3))
    g_line3.fill(np.nan)
    g_line3 = [self.gamma_c(i,n,a1x1,esl1,esl2,nsi1,nsi2,nc1,nc2,nh1,nh2,self.esic,self.uc_grit) for i in range(-1,2)]
    self.ax1.plot(xax,g_line3, "-", color="k", linewidth=1)
    x_txt = -0.48
    y_txt = self.gamma_c(x_txt,n,a1x1,esl1,esl2,nsi1,nsi2,nc1,nc2,nh1,nh2,self.esic,self.uc_grit)
    ang = np.arctan(g_line3[1]-g_line3[0])*180/np.pi*self.myaspect
    plt.text( x_txt, y_txt-0.00, r"$2\times 2_{\textrm{\footnotesize 5Si$\otimes$0C}}$", rotation=ang, fontsize=10)

#    inset = axes([0.45, 0.15, 0.3, 0.25], axisbg='w')
#    inset.xaxis.set_minor_locator(MultipleLocator(0.01))
#    inset.yaxis.set_minor_locator(MultipleLocator(0.01))

#    x = np.arange(-1.0, 1.0, 0.0001)
#    inset.fill_between(x, -2, 2, where=x>=self.uc_grit-self.uc_grit, facecolor='black', alpha=0.05, interpolate=False)
#    setp(inset,xlim=(-0.05,0.025), xticks=[-0.04,-0.02,0,0.02], ylim=(-0.15,-0.05), yticks=[-0.15,-0.10,-0.05])
#    plot(self.grit_line_x,self.grit_line_y, ":", color='k', linewidth=1)
#    plot(xax,g_line1, color='k', linewidth=1)
#    plot(xax,g_line2, color='k', linewidth=1)
#    plot(xax,g_line3, color='k', linewidth=1)
#    self.ax1.annotate(r"", xy=(0, -0.095), xytext=(-0.03, -0.3),
#                    arrowprops=dict(arrowstyle='fancy', facecolor='w', connectionstyle='arc3,rad=-0.2'), color='k')
#    inset.annotate(r"$3\times 3_{\textrm{\footnotesize 13Si$\otimes$0C}}$", xy=(-0.04,self.gamma_3x3_13si(-0.045)), xytext=(-0.03,-0.13),
#                 arrowprops=dict(arrowstyle='simple', facecolor='k', connectionstyle='arc3,rad=-0.4', linewidth=0.1), color='k', fontsize=10)
#    inset.annotate(r"$2\times 2_{\textrm{\footnotesize 5Si$\otimes$0C}}$", xy=(-0.005,self.gamma_2x2_5si(-0.005)), xytext=(-0.0135,-0.07),
#                 arrowprops=dict(arrowstyle='simple', facecolor='k', connectionstyle='arc3,rad=0.4', linewidth=0.1), color='k', fontsize=10)
#    inset.annotate(r"$2\times 2_{\textrm{\footnotesize 1Si$\otimes$0C}}$", xy=(0.0175,self.gamma_2x2_1si(0.00175)), xytext=(0.001,-0.11),
#                 arrowprops=dict(arrowstyle='simple', facecolor='k', connectionstyle='arc3,rad=0.4', linewidth=0.1), color='k', fontsize=10)

  def AddGeoPics(self):
    sic_geo = read_png('../sic_geo.png')
    ite_geo = read_png('../ite_geo.png')
    dia_geo = read_png('../dia_geo.png')

    self.ax1.imshow(sic_geo, extent=[-0.53,-0.47,0.95,1.15], origin='upper', alpha=1, aspect='auto', interpolation='none')
    self.ax1.imshow(ite_geo, extent=[-0.03,0.03,0.95,1.15], origin='upper', alpha=1, aspect='auto', interpolation='none')
    self.ax1.imshow(dia_geo, extent=[0.04,0.10,0.95,1.15], origin='upper', alpha=1, aspect='auto', interpolation='none')

if(__name__=="__main__"):
  this = GeoDeterm()

  this.makeDotPlot()

  this.readData()

  area = np.nan
  emin = np.nan
  nsi = np.nan
  nc = np.nan
  nh = np.nan
  geo_thresh = 1E-2
  area, geo_hits, emin, nsi, nc, nh = this.processData(geo_thresh)

  num_bins = 1000
  this.makePlot(num_bins)

  pd = PhaseDiagram()

  pd.setPlotParams()
  pd.setRefVals()
  pd.PlotRefVals()

  pd.PlotCalcVals(area, emin, nsi, nc, nh)

  pd.AddGeoPics()

  plt.savefig("phase.png",dpi=170)
