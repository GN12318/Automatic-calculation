from ase import Atoms, Atom
from ase.io import read
import pymatgen, os, math, shutil
from pymatgen.io.vasp.inputs import Kpoints, Poscar
def automatic_gamma_density(file, desired_directory, dens, output_directory):

       
    flie = file
    filepath = desired_directory + file 
    output_KPOINTS = output_directory + f"{file}.kp"
    dens = dens  
        
    # Read file
    p = Poscar.from_file(filepath)
        
    #length of the reciprocal lattice vector
    latt = p.structure.lattice 
    R = latt.reciprocal_lattice
    V = [R.a, R.b, R.c]  #recp length
    num_div = [int(round(max(1, n/(2*dens*math.pi)))) for n in V]
        
    """for n in B:
        N = round(max(1, n/(2*dens*math.pi)))   V Wang:Vaspkit DOI:10.1016/j.cpc.2021.108033"""
               
    style = Kpoints.supported_modes.Gamma
    comment = f"generate with grid density = {dens:.02f} "
    num_kpts = 0
    #write kpoints file    
    with open(output_KPOINTS, 'w') as file:
        file.write(f"{comment}\n")
        file.write(f"{num_kpts}\n")
        file.write(f"{style}\n")
        file.write(f"{num_div[0]}" + ' ' + f"{num_div[1]}" + ' ' + f"{num_div[2]}\n")
        file.write('0.0' + ' ' + '0.0' + ' ' + '0.0')
        
    return Kpoints(comment, num_kpts, style, [num_div], (0, 0, 0))
	
mydir = '221001'
file_path = "/public3/home/scg5072/workdir/ALSI/TEST/POSCAR/"
k_outpath = "/public3/home/scg5072/workdir/ALSI/KOUT/"
files = os.listdir(file_path)
for file in files:
    atoms = read(file_path + file, format = 'vasp')
    print(atoms)
    MK = automatic_gamma_density(file, file_path, 0.03, k_outpath)
    print(MK)
    
    kps = MK.kpts[0]
#     print(kps)
         
    # build calc
    calc = Vasp(
            istart = 0,
            icharg = 2,
            ibrion = -1,
            isif = 2,
            ismear = 0,
            sigma = 0.1,
            nsw = 0,
            ediff = 1E-5,
            ediffg = -1E-3,
            encut = 450,
            lwave = False, 
            lcharg = False,
            gamma = True,
	    xc = 'PBE',
            kpts = kps,
            directory = mydir) 
    
    atoms.calc = calc
    atoms.get_potential_energy()
    
    print(atoms.get_potential_energy())
    os.rename('/public3/home/scg5072/workdir/ALSI/OOUT/' + 'OUTCAR', '/public3/home/scg5072/workdir/ALSI/OOUT/' + f"{file}.out")
