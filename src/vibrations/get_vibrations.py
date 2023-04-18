"""
  modified by Simon Erker
  
  Python version of numerical_vibrations.pl based on numpy and scipy
 
  The script creates  folders containing the aims input files for single 
  calculations (mode = 0). These are the 6M  calculations for the M 
  unconstrained atoms in the geometry.in in the execution path. The atoms are 
  displaced by +/-delta (default 0.0025\AA) in x, y, and z. The control.in 
  is copied to all folders. After perfoming the DFT calculations the 
  vibrational modes can be calculated (mode = 1). mode = 2 combines the two 
  other modes (only usefule if aims-binary is provided). Only the atoms without 
  constrain_relaxation .true. in geometry.in will be displaced. 
  (-d <delta> is optional input see help)

  > get_vibrations.py [options] <name> <mode>

  Output:
 {name}.xyz - molecular file that contains the structure, the eigenmodes, 
              their frequencies, and their IR intensity.
 ir.{name}.dat - IR intensities

 This script is considered beta stage. Please use with causion. The idea 
 behind this script is the application to vibrations of molecular adsorbates 
 on (large) surfaces. Typically, only a  few atoms have to be displaced and 
 turn over times are too large to submit only one calculation (loads of CPUs 
 could be required). The pearl script produces a lot of temporary output 
 files and it does a lot of loop evalutions for large systems (Can take 24h 
 for 500 atoms). No additional fortran routine required. The post-processing 
 is based on numpy/scipy. Plotting the IR-spectrum requires matplotlib.

"""

from scipy import constants, signal
from numpy import array, zeros, ones, float64, append, sqrt, arange, newaxis, sum
import copy, os, shutil
 
USAGE = """%prog [options] <name> <mode>

<name> will be used as prefix for all output files.
<mode> selects preprocessing (and calculation if a path to FHI-aims binary is 
       is provides) and postprocessing - calcultion of vibrational modes from 
       forces obtained from DFT calculations.
       
"""

def element(elem):
    """Element symbols"""
    data = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", 
            "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", 
            "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", 
            "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", 
            "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", 
            "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", 
            "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", 
            "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", 
            "Pa", "U"]
    try: 
        return data[elem-1] 
    except TypeError: 
        return data.index(elem.capitalize())+1
        
def mass(elem):
    """Atomic masses"""
    data = [1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.0107, 14.0067, 
            15.9994, 18.9984032, 20.1797, 22.9897, 24.305, 26.9815, 28.0855, 
            30.973762, 32.065, 35.453, 39.948, 39.0983, 40.078, 44.9559, 
            47.867, 50.9415, 51.9961, 54.938, 55.845, 58.9332, 58.6934, 
            63.546, 65.409, 69.723, 72.64, 74.9216, 78.96, 79.904, 83.8, 
            85.4678, 87.62, 88.90585, 91.224, 92.9064, 95.94, 98, 101.07, 
            102.9055, 106.42, 107.8682, 112.411, 114.818, 118.71, 121.76,
            127.6, 126.90447, 131.293, 132.9055, 137.327, 138.9055, 140.116, 
            140.9077, 144.24, 145, 150.36, 151.964, 157.25, 158.9253, 162.5, 
            164.9303, 167.259, 168.9342, 173.04, 174.967, 178.49, 180.9479, 
            183.84, 186.207, 190.23, 192.217, 195.078, 196.9665, 200.59,
            204.3833, 207.2, 208.9804, 209, 210, 222, 223, 226,227, 232.0381, 
            231.0359, 238.02891]
    try:
        return data[elem-1]
    except TypeError:
        return data[element(elem)-1]
 
class Atom:
    """Atom object"""
    def __init__(self, kind, coord, constraint=None, force=None):
        self.kind = kind
        self.coord = array(coord,float64)
        self.constraint = constraint if constraint else False
        self.force = array(force,float64) if force else zeros(3)
    def Z(self):
        return element(self.kind)
    def mass(self):
        return mass(self.Z())
    def copy(self):
        return copy.deepcopy(self)
    def to_str(self, elem_symb=True):
        elem = (element if elem_symb else str)(self.Z())
        if self.constraint:
          constraint = ' constrain_relaxation .true.\n'
        else:
          constraint = ''
        line='atom{0:14.8f}'.format(self.coord[0])
        line=line+'{0:14.8f}'.format(self.coord[1])
        line=line+'{0:14.8f}'.format(self.coord[2])
        line=line+'{0:>4}\n'.format(self.kind)
        line=line+constraint
        return line
class structure:
    """Structure of atoms object"""
    def __init__(self, atoms=None, vacuum_level=None, desc=None,\
                 periodic=None):
        self.atoms = atoms if atoms else array([])
        self.desc = desc if desc else ""
        self.vacuum_level = vacuum_level
        self.constrained=array([],bool)
        self.lattice_vector=zeros([0,3])
        self.periodic = periodic if periodic else False
        
    @classmethod
    def copy(self):
        return copy.deepcopy(self)
    def n(self):
        return len(self.atoms)
    def mass(self):
        return sum([mass(a.Z) for a in self.atoms])
    def join(self, add):
        self.atoms=append(self.atoms,add.copy())
        self.constrained=append(self.constrained,add.constraint)
    def to_str(self):
      line=''
      if self.periodic:
        line=line+'set_vacuum_level{0:14.8f}\n'.format(self.vacuum_level)
        for i in range(3):
          line=line+'lattice_vector{0:14.8f}'.format(self.lattice_vector[i,0])
          line=line+'{0:14.8f}'.format(self.lattice_vector[i,1])
          line=line+'{0:14.8f}\n'.format(self.lattice_vector[i,2])
      for j in range(self.n()):
        line=line+self.atoms[j].to_str()
      return line

def split_line(lines):
  """Split input line"""
  line_array=array(lines.strip().split(' '))
  line_vals=line_array[line_array!='']
  return line_vals

def replace_submission(template_job, name, counter, filename):
    """
     Prepare submission script
     Only jobname and output file will be replaced
     Feel free to add more
    """
    template_job=template_job.replace('<jobname>',name+"_"+str(counter))
    template_job=template_job.replace('<outfile>',filename)
    job_out=open('job.sh','w')
    job_out.write(template_job)
    job_out.close()

def lorentz(pi,g,x):
    """
      Lorentzian function
    """
    x0 = (x.max()-x.min())-x.min()
    lr = ((g**2)/((x-0.5*x0)**2+g**2))
    return lr - lr.min()
    
def main():
  import optparse
  from numpy import sum

  # Parse command line
  parser = optparse.OptionParser(usage=USAGE)
  parser.add_option("-p", "--plot", action="store_true",
                    help="Generate pdf with IR-spectrum, broadened with Lorentzian")
  parser.add_option("-i", "--info", action="store_true",
                      help="Set up/ Calculate vibrations & quit")
  parser.add_option("-s", "--suffix", action="store",
                    help="Call suffix for binary e.g. 'mpirun -n 4 '",
                    default='')
  parser.add_option("-r", "--run", action="store",
                    help="path to FHI-aims binary",default='')
  parser.add_option("-x", "--relax", action="store_true",
                    help="Relax initial geometry")
  parser.add_option("-m", "--molden", action="store_true",
                    help="Output in molden format")
  parser.add_option("-w", "--distort", action="store_true",
                    help="Output geometry distorted along imaginary modes")
  parser.add_option("-t", "--submit", action="store",
                    help="""\
Path to submission script, string <jobname>
will be replaced by name + counter, string 
                            <outfile> will be replaced by filename""")
  parser.add_option("-d", "--delta", action="store", type="float",
                    help="Displacement", default=0.0025)
  parser.add_option("-b", "--broadening", action="store", type="float",
                    help="Broadening for IR-spectrum in cm^{-1}", default=5)

  options, args = parser.parse_args()
  if options.info:
      print(__doc__)
      sys.exit(0)
  if len(args) != 2:
      parser.error("Need exactly two arguments")
  
  AIMS_CALL=options.suffix+' '+options.run
  hessian_thresh = -1
  name=args[0]
  mode=args[1] 
  delta=options.delta
  broadening=options.broadening

  run_aims=False
  if options.run!='': run_aims=True
  
  submit_script = options.submit is not None
  
  if options.plot:
    import matplotlib as mpl
    mpl.use('Agg') 
    from pylab import figure

  if options.plot or mode=='1' or mode=='2':
    from pylab import savetxt, transpose, eig, argsort, sort, sign, pi, dot, sum, linspace, argmin, r_, convolve
    
  # Constant from scipy.constants
  bohr=constants.value('Bohr radius')*1.e10
  hartree=constants.value('Hartree energy in eV')
  at_u=constants.value('atomic mass unit-kilogram relationship')
  eV=constants.value('electron volt-joule relationship')
  c=constants.value('speed of light in vacuum')
  Ang=1.0e-10
  hbar=constants.value('Planck constant over 2 pi')
  Avo=constants.value('Avogadro constant')
  kb=constants.value('Boltzmann constant in eV/K')
  pi=constants.pi
  hessian_factor   = eV/(at_u*Ang*Ang) 
  grad_dipole_factor=(eV/(1./(10*c)))/Ang  #(eV/Ang -> D/Ang)
  ir_factor = 1
  
  # Asign all filenames
  inputgeomerty = 'geometry.in.'+name
  inputcontrol  = 'control.in.'+name
  atomicmasses  = 'masses.'+name+'.dat'; 
  xyzfile       = name+'.xyz';
  moldenname    =name+'.molden';
  hessianname   = 'hessian.'+name+'.dat'; 
  graddipolename   = 'grad_dipole.'+name+'.dat'; 
  irname   = name+'.ir'; 
  deltas=array([-delta,delta])
  coeff=array([-1,1])
  c_zero = - 1. / (2. * delta)


  f=open('control.in','r')                   # read control.in template
  template_control=f.read()
  f.close

  if submit_script:
    f=open(options.submit,'r')               # read submission script template
    template_job=f.read()
    f.close

  folder=''                                  # Dummy
  ########### Central Point ##################################################
  if options.relax and (mode=='0' or mode=='2'):
    # First relax input geometry
    filename=name+'.out'
    folder=name+'_relaxation' 
    if not os.path.exists(folder): os.mkdir(folder)            # Create folder
    shutil.copy('geometry.in', folder+'/geometry.in')          # Copy geometry
    new_control=open(folder+'/control.in','w')
    new_control.write(template_control+'relax_geometry trm 1E-3\n') # Relax!
    new_control.close()
    os.chdir(folder)                             # Change directoy
    print('Central Point')
    if run_aims:
      os.system(AIMS_CALL+' > '+filename)       # Run aims and pipe the output into a file named 'filename'
    if submit_script: replace_submission(template_job, name, 0, filename)
    os.chdir('..') 

  ############################################################################
  # Check for relaxed geometry
  if os.path.exists(folder+'/geometry.in.next_step'):  
    geometry=open(folder+'/geometry.in.next_step','r')
  else:
    geometry=open('geometry.in','r')
    
  # Read input geometry  
  n_line=0
  struc=structure()
  lines=geometry.readlines()

  for line in lines:
    n_line= n_line+1
    if line.rfind('set_vacuum_level')!=-1:   # Vacuum Level
      struc.vacuum_level=float(split_line(line)[-1])
    if line.rfind('lattice_vector')!=-1:    # Lattice vectors and periodic
      lat=split_line(line)[1:]
      struc.lattice_vector=append(struc.lattice_vector,float64(array(lat))[newaxis,:],axis=0)
      struc.periodic=True
    if line.rfind('atom')!=-1:              # Set atoms
      line_vals=split_line(line)
      at=Atom(line_vals[-1],line_vals[1:-1])
      if n_line<len(lines):
        nextline=lines[n_line]
        if nextline.rfind('constrain_relaxation')!=-1: # constrained?
            at=Atom(line_vals[-1],line_vals[1:-1],True)
        else:
            at=Atom(line_vals[-1],line_vals[1:-1])
      struc.join(at)
  geometry.close()  
  n_atoms= struc.n()
  n_constrained=n_atoms-sum(struc.constrained)
  print('Total number of atoms found: '+str(n_atoms))
  print('Number of atoms considered in fequency calcualtion: '+str(n_constrained))

  # Atomic mass file
  mass_file=open(atomicmasses,'w')
  mass_vector=zeros([0])
  for at_unconstrained in struc.atoms[struc.constrained==False]:
      mass_vector=append(mass_vector,ones(3)*1./sqrt(at_unconstrained.mass()))
      line='{0:10.5f}'.format(at_unconstrained.mass())
      for i in range(3):
          line=line+'{0:11.4f}'.format(at_unconstrained.coord[i])
      line=line+'{0:}\n'.format(at_unconstrained.kind)
      mass_file.writelines(line)
  mass_file.close()

  # Init
  dip = zeros([n_constrained*3,3])
  hessian = zeros([n_constrained*3,n_constrained*3])
  index=0
  counter=1
  
  # Set up / Read folders for displaced atoms
  for atom in arange(n_atoms)[struc.constrained==False]:
    for coord in arange(3):
      for delta in deltas:
        #filename=name+'.i_atom_'+str(atom)+'.i_coord_'+str(coord)+'.displ_'+str(delta)+'.out'
        filename='aims.out'
        folder=name+'.i_atom_'+str(atom)+'.i_coord_'+str(coord)+'.displ_'+str(delta)

        if mode=='0' or mode=='2':   # Put new geometry and control.in into folder
          struc_new=copy.deepcopy(struc)
          struc_new.atoms[atom].coord[coord]=struc_new.atoms[atom].coord[coord]+delta
          geoname='geometry.i_atom_'+str(atom)+'.i_coord_'+str(coord)+\
                  '.displ_'+str(delta)+'.in'
          if not os.path.exists(folder): os.mkdir(folder)
          new_geo=open(folder+'/geometry.in','w')
          newline='#\n# temporary structure-file for finite-difference calculation of forces\n'
          newline=newline+'# displacement {0:8.4f} of \# atom '.format(delta)+'{0:5} direction {1:5}\n#\n'.format(atom,coord)
          new_geo.writelines(newline+struc_new.to_str())
          new_geo.close()
          new_control=open(folder+'/control.in','w')
          template_control=template_control.replace('relax_geometry','#relax_geometry')
          if struc.periodic:
            new_control.write(template_control+'compute_forces .true. \n'+'final_forces_cleaned '+'.true. \n'+'evaluate_work_function \n')
          else:
            new_control.write(template_control+'compute_forces .true. \n'+'final_forces_cleaned '+'.true. \n'+'output dipole \n')
          new_control.close()
          os.chdir(folder)                                   # Change directoy
          print('Processing atom: '+str(atom+1)+'/'+str(n_atoms)+', coord.: '+str(coord+1)+'/'+str(3)+', delta: '+str(delta))
          if run_aims:                           
            os.system(AIMS_CALL+' > '+filename) # Run aims and pipe the output into a file named 'filename'
          if submit_script: replace_submission(template_job, name, counter, filename)
          #os.system('sbatch -N 4 job.sh') # Mind the environment variables <======================================== Modify according to you system
          counter=counter+1
          os.chdir('..') 

        if mode=='1' or mode=='2':   # Read output 
          forces_reached=False
          atom_count=0
          data=open(folder+'/'+filename)
          for line in data.readlines():
            if struc.periodic:
              if line.rfind('| Dipole correction potential jump')!=-1:
                dip_jump = float(split_line(line)[-2]) #Periodic
            else:
              if line.rfind('| Total dipole moment [eAng]')!=-1:
                dip_jump = float64(split_line(line)[-3:]) #Cluster
            if forces_reached and atom_count<n_atoms: #Read Forces
              struc.atoms[atom_count].force=float64(split_line(line)[2:])
              atom_count=atom_count+1
              if atom_count==n_atoms:
                forces_reached=False
            if line.rfind('Total atomic forces')!=-1:
              forces_reached=True
          data.close()
          if struc.periodic:
            dip[index,2]=dip[index,2]+dip_jump*coeff[deltas==delta]*c_zero
          else:
            dip[index,:]=dip[index,:]+dip_jump*coeff[deltas==delta]*c_zero

          forces=array([])
          for at_unconstrained in struc.atoms[struc.constrained==False]:
              forces=append(forces,coeff[deltas==delta]*at_unconstrained.force)
          hessian[index,:]=hessian[index,:]+forces*c_zero
      
      index=index+1

  if mode=='1' or mode=='2': # Calculate vibrations
    print('Entering hessian diagonalization')
    print('Number of atoms                = '+str(n_atoms))
    print('Name of Hessian input file     = '+hessianname)
    print('Name of grad dipole input file = '+graddipolename)
    print('Name of Masses  input file     = '+atomicmasses)
    print('Name of XYZ output file        = '+xyzfile)
    print('Threshold for Matrix elements  = '+str(hessian_thresh))
    if (hessian_thresh < 0.0): print('     All matrix elements are taken into account by default\n')
    savetxt(hessianname,hessian)
    savetxt(graddipolename,dip)

    mass_mat=mass_vector[:,newaxis]*mass_vector[newaxis,:]
    hessian[abs(hessian)<hessian_thresh]=0.0
    hessian=hessian*mass_mat*hessian_factor
    hessian=(hessian+transpose(hessian))/2.
    # Diagonalize hessian (scipy)
    print('Solving eigenvalue system for Hessian Matrix')
    freq, eig_vec = eig(hessian)
    print('Done ... ')
    eig_vec=eig_vec[:,argsort(freq)]
    freq=sort(sign(freq)*sqrt(abs(freq)))
    ZPE=hbar*(freq)/(2.0*eV)
    freq = (freq)/(200.*pi*c)
    
    grad_dipole = dip * grad_dipole_factor
    eig_vec = eig_vec*mass_vector[:,newaxis]*ones(len(mass_vector))[newaxis,:]
    infrared_intensity = sum(dot(transpose(grad_dipole),eig_vec)**2,axis=0)*ir_factor
    reduced_mass=sum(eig_vec**2,axis=0)
    norm = sqrt(reduced_mass)
    eig_vec = eig_vec/norm
    
    # The rest is output, xyz, IR,...
    print('Results\n')
    print('List of all frequencies found:')
    print('Mode number      Frequency [cm^(-1)]   Zero point energy [eV]   IR-intensity [D^2/Ang^2]')
    for i in range(len(freq)):
      print('{0:11}{1:25.8f}{2:25.8f}{3:25.8f}'.format(i+1,freq[i],ZPE[i],infrared_intensity[i]))
    print('\n')
    print('Summary of zero point energy for entire system:')
    print('| Cumulative ZPE               = {0:15.8f} eV'.format(sum(ZPE)))
    print('| without first six eigenmodes = {0:15.8f} eV\n'.format(sum(ZPE)-sum(ZPE[:6])))
    print('Stability checking - eigenvalues should all be positive for a stable structure. ')
    print('The six smallest frequencies should be (almost) zero:')
    string=''
    for zz in freq[:6]: string=string+'{0:25.8f}'.format(zz)
    print(string)
    print('Compare this with the largest eigenvalue, ')
    print('{0:25.8f}'.format(freq[-1]))
    
    nums=arange(n_atoms)[struc.constrained==False]
    nums2=arange(n_atoms)[struc.constrained]
    newline=''
    newline_ir='[INT]\n'
    if options.molden:
      newline_molden='[Molden Format]\n[GEOMETRIES] XYZ\n'
      newline_molden=newline_molden+'{0:6}\n'.format(n_atoms)+'\n'
      for i_atoms in range(n_constrained):
        newline_molden=newline_molden+'{0:6}'.format(struc.atoms[nums[i_atoms]].kind)
        for i_coord in range(3):
          newline_molden=newline_molden+'{0:10.4f}'.format(struc.atoms[nums[i_atoms]].coord[i_coord])
        newline_molden=newline_molden+'\n'
      newline_molden=newline_molden+'[FREQ]\n'   
      for i in range(len(freq)):
        newline_molden=newline_molden+'{0:10.3f}\n'.format(freq[i])
      newline_molden=newline_molden+'[INT]\n' 
      for i in range(len(freq)):
        newline_molden=newline_molden+'{0:17.6e}\n'.format(infrared_intensity[i])
      newline_molden=newline_molden+'[FR-COORD]\n'
      newline_molden=newline_molden+'{0:6}\n'.format(n_atoms)+'\n'
      for i_atoms in range(n_constrained):
        newline_molden=newline_molden+'{0:6}'.format(struc.atoms[nums[i_atoms]].kind)
        for i_coord in range(3):
          newline_molden=newline_molden+'{0:10.4f}'.format(struc.atoms[nums[i_atoms]].coord[i_coord]/bohr)
        newline_molden=newline_molden+'\n'
      newline_molden=newline_molden+'[FR-NORM-COORD]\n'
    
    for i in range(len(freq)):
      newline=newline+'{0:6}\n'.format(n_atoms)
      if freq[i]>0:
        newline=newline+'stable frequency at '
      elif freq[i]<0:
        newline=newline+'unstable frequency at '
        if options.distort and freq[i]<-50:
          struc_new=copy.deepcopy(struc)
          for i_atoms in range(n_constrained):
            for i_coord in range(3):
              struc_new.atoms[i_atoms].coord[i_coord]=\
              struc_new.atoms[i_atoms].coord[i_coord]+eig_vec[(i_atoms)*3+i_coord,i]
          geoname=name+'.distorted.vibration_'+str(i+1)+'.geometry.in'
          new_geo=open(geoname,'w')
          newline_geo='#\n# distorted structure-file for based on eigenmodes\n'
          newline_geo=newline_geo+'# vibration {0:5} :{1:10.3f} 1/cm\n#\n'.format(i+1,freq[i])
          new_geo.writelines(newline_geo+struc_new.to_str())
          new_geo.close()
      elif freq[i]==0:
        newline=newline+'translation or rotation '
      newline=newline+'{0:10.3f} 1/cm IR int. is '.format(freq[i])
      newline=newline+'{0:10.4e} D^2/Ang^2; red. mass is '.format(infrared_intensity[i])
      newline=newline+'{0:5.3f} a.m.u.; force const. is '.format(1.0/reduced_mass[i])
      newline=newline+'{0:5.3f} mDyne/Ang.\n'.format(((freq[i]*(200*pi*c))**2)*(1.0/reduced_mass[i])*at_u*1.e-2)
      if options.molden: newline_molden=newline_molden+'vibration {0:6}\n'.format(i+1)
      for i_atoms in range(n_constrained):
        newline=newline+'{0:6}'.format(struc.atoms[nums[i_atoms]].kind)
        for i_coord in range(3):
          newline=newline+'{0:10.4f}'.format(struc.atoms[nums[i_atoms]].coord[i_coord])
        for i_coord in range(3):
          newline=newline+'{0:10.4f}'.format(eig_vec[(i_atoms)*3+i_coord,i])
          if options.molden: newline_molden=newline_molden+'{0:10.4f}'.format(eig_vec[(i_atoms)*3+i_coord,i]/bohr)
        newline=newline+'\n'
        if options.molden: newline_molden=newline_molden+'\n'
      for i_atoms in range(n_atoms-n_constrained):
        newline=newline+'{0:6}'.format(struc.atoms[nums2[i_atoms]].kind)
        for i_coord in range(3):
          newline=newline+'{0:10.4f}'.format(struc.atoms[nums2[i_atoms]].coord[i_coord])
        for i_coord in range(3):
          newline=newline+'{0:10.4f}'.format(0.0)
        newline=newline+'\n'
      newline_ir=newline_ir+'{0:11}{1:25.8f}{2:25.8f}{3:25.8f}\n'.format(i+1,freq[i],ZPE[i],infrared_intensity[i])
    xyz=open(xyzfile,'w')
    xyz.writelines(newline)
    xyz.close()
    ir=open(irname,'w')
    ir.writelines('#Mode number      Frequency [cm^(-1)]   Zero point energy [eV]   IR-intensity [D^2/Ang^2]')
    ir.writelines(newline_ir)
    ir.close()
    
    if options.molden:
      molden=open(moldenname,'w')
      molden.writelines(newline_molden)
      molden.close()
    
    if (mode=='1' or mode=='2') and options.plot:
      x=linspace(freq.min()-500,freq.max()+500,1000)
      z=zeros(len(x))
      for i in range(len(freq)):
          z[argmin(abs(x-freq[i]))]=infrared_intensity[i]
      window_len=150
      lorentzian=lorentz(pi,broadening,arange(250))#signal.gaussian(window_len,broadening)
      s=r_[z[window_len-1:0:-1],z,z[-1:-window_len:-1]]
      z_convolve=convolve(lorentzian/lorentzian.sum(),s,mode='same')[window_len-1:-window_len+1]
      fig=figure(0)
      ax=fig.add_subplot(111)
      ax.plot(x,z_convolve,'r',lw=2)
      ax.set_xlim([freq.min()-500,freq.max()+500])
      ax.set_ylim([-0.01,ax.get_ylim()[1]])
      ax.set_yticks([])
      ax.set_xlabel('Frequency [1/cm]',size=20)
      ax.set_ylabel('Intensity [a.u.]',size=20)
      fig.savefig(name+'_IR_spectrum.pdf')
      
    print('\n Done. ')
    
if __name__ == "__main__":
    main()
