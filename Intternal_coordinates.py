import numpy as np
import sys

def distance(x1, x2): #function for calculating bond length, the parameteres x1,x2 are the coordiantes presnet in an array
	return np.linalg.norm(x2-x1) #np.linalg.norm this function under numpy calculate the magnitude of a vector

# Function to calculate the bond angle between three atoms	
def angle (x1, x2, x3):
	v1= x1-x2 # Vector from atom 2 to atom 1
	v2= x3-x2 # Vector from atom 3 to atom 2
	cos_theta= np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)) # Cosine of the angle
	
	return np.degrees((np.arccos(cos_theta))) # angle to degrees

# Function to calculate the dihedral between three atoms
def dihedral(x1, x2, x3, x4):
	v1 = x2 - x1  # Vector from atom 2 to atom 1
	v2 = x3 - x2  # Vector from atom 3 to atom 2
	v3 = x4 - x3  # Vector from atom 4 to atom 3

	normal1 = np.cross(v1, v2) #another vector that we get by doing cross product
	normal2 = np.cross(v2, v3)

	cos_phi = np.dot(normal1, normal2) / (np.linalg.norm(normal1) * np.linalg.norm(normal2))
    
	return np.degrees(np.arccos(cos_phi))

def read_pdb(filename):
	coords={}
	connectivity={}
	with open(filename,'r') as file:
		for line in file:
			if line.startswith(('ATOM','HETATM')):
				atomID = int(line[6:11].strip())
				atom=line[13:16].strip()
				x=float(line[30:38].strip())
				y=float(line[38:46].strip())
				z=float(line[46:54].strip())
				coords[atomID]=(atom,np.array([x,y,z]))
			if line.startswith('CONECT'):
				parts=line.split()
				atomID=int(parts[1])
				bonded_atom=[int(bond) for bond in parts[2:]]
				connectivity[atomID]=bonded_atom
			
	return coords,connectivity

def analysis(pdb_file, input_file, output):
	coords,connectivity = read_pdb(pdb_file)

	with open(input_file, 'r') as infile, open(output, 'w') as outfile:
		for line in infile:
			atom_ids = [int(i) for i in line.split()]

			if len(atom_ids) == 2:
				if atom_ids[1] in connectivity.get(atom_ids[0],[]):
					x1 = coords[atom_ids[0]][1]
					x2 = coords[atom_ids[1]][1]
					bond_length = distance(x1, x2)
					outfile.write(f"Bond length between atom {atom_ids[0]} and {atom_ids[1]}: {bond_length:.3f}\n\n")
				else:
					outfile.write(f"Atoms {atom_ids[0]} and {atom_ids[1]} are not directly connected.\n\n")
                   
                
			elif len(atom_ids) == 3:
				if atom_ids[0] in connectivity.get(atom_ids[1], []) and atom_ids[2] in connectivity.get(atom_ids[1], []) :
					x1 = coords[atom_ids[0]][1]
					x2 = coords[atom_ids[1]][1]
					x3 = coords[atom_ids[2]][1]
					bond_angle = angle(x1, x2,x3)
					outfile.write(f"Bond angle between atom {atom_ids[0]} , {atom_ids[1]} and {atom_ids[2]}: {bond_angle:.3f}\n\n")
				else:
					outfile.write(f"Atoms {atom_ids[0]} , {atom_ids[1]} and {atom_ids[2]} are not connected in continous fation.\n\n")
               
			elif len(atom_ids) == 4:
				x1 = coords[atom_ids[0]][1]
				x2 = coords[atom_ids[1]][1]
				x3 = coords[atom_ids[2]][1]
				x4 = coords[atom_ids[3]][1]
				dihedral_angle = dihedral(x1, x2,x3,x4)
				
				if atom_ids[0] in connectivity.get(atom_ids[1], []) and atom_ids[2] in connectivity.get(atom_ids[1], []) and atom_ids[3] in connectivity.get(atom_ids[2]) :
					outfile.write(f"Dihedral angle between atom {atom_ids[0]} , {atom_ids[1]} ,{atom_ids[2]} and {atom_ids[3]}: {dihedral_angle:.3f}  (The entered atom IDs or their sequence form a connected structure.\n\n")
				
				
					
					
				elif atom_ids[0] in connectivity.get(atom_ids[1], []) and atom_ids[2] in connectivity.get(atom_ids[1], []) :
					outfile.write(f"Dihedral angle between atom {atom_ids[0]} , {atom_ids[1]} ,{atom_ids[2]} and {atom_ids[3]}: {dihedral_angle:.3f}  (atom {atom_ids[0]} , {atom_ids[1]} and {atom_ids[2]} are connected).\n\n")
					
				
				elif atom_ids[0] in connectivity.get(atom_ids[1], []) and atom_ids[2] in connectivity.get(atom_ids[3], []) :
					
					outfile.write(f"Dihedral angle between atom {atom_ids[0]} , {atom_ids[1]} ,{atom_ids[2]} and {atom_ids[3]}: {dihedral_angle:.3f}  (atom {atom_ids[0]} and {atom_ids[1]} are connected , {atom_ids[2]} and {atom_ids[3]} are connected).\n\n")
					
					
				else:
			
					outfile.write(f"Dihedral angle between atom {atom_ids[0]} , {atom_ids[1]} ,{atom_ids[2]} and {atom_ids[3]}: {dihedral_angle:.3f}  (The entered atom IDs or their sequence do not form a connected structure)\n\n")
					
					


pdb_filename = sys.argv[1]
input_filename = sys.argv[2]
output_filename = sys.argv[3]
analysis(pdb_filename, input_filename, output_filename)

	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
	
	
		
			 
        
