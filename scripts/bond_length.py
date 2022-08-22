import sys , math
import numpy as np

#common used atoms radius
A_radi = {  'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71,
  'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 'Br': 1.14, 'I' : 1.33, 'He': 0.30,
  'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 'Na': 1.02,
  'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75,
  'Ti': 0.86, 'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
  'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
  'Se': 1.17, 'Kr': 1.03, 'X' : 0.00}
#threshold factor
k = 1.2
# The three following functions are just for printing DATA
def print_bond_graph():
	print("Bond Graph")
	for i in range(n):
		print('{num:>4}{name:>2}'.format(name = atomic_graph[i+1][0], num = i+1), ':' , end = ' ')
		for j in atomic_graph[i+1][1:]:
			print(j,end=' ')
		print()
	print()

def print_original_info():
	print("Original Info")
	for i in atoms:
		print('{J:>4}'.format(J=i[0])	 ,end = '')
		for j in i[1:]:
			print('{I:>10}'.format(I=j) ,end='')
		print()
	print()

def print_bonds():
	print('Bonds')
	fkey = bonds.keys()
	for i in fkey:
		print ('{:<8}{:>8}{:>12.6f}'.format('('+str(i)+')' , bonds[i][1] , bonds[i][0]))

# Reads file from input 
def read(a : str):
	f = open(a , 'r')
	n = int(f.readline().strip('\n') )
	molecule_name = f.readline().strip('\n')
	for i in range(n):
		a = f.readline()
		if a:
			atoms.append(a.strip('\n').split())
		else: break
	return n , molecule_name
# Calculats length of bond from given coordinations
def length(a : list ,b : list):
	R = 0.0
	for i in range(3):
		x = float(b[i])	 - float(a[i])
		R += x**2
	return np.sqrt(R)
# Checks if atoms can form bond in such distance
def if_bond(a,b):
	RR = length(a[1:] , b[1:])
	if  RR <= k*(A_radi[a[0]]+A_radi[b[0]]) 	:
		return RR , True
	else:
		return RR , False
# Checks all DATA and applies the necessary function 
def surf(atoms):
	for i in range(n):
		for j in range(n-1,i,-1):
			key = bonds.keys()
			L , flag = if_bond(atoms[i] , atoms[j])
			if flag:
				atomic_graph[i+1].append(j+1)
				atomic_graph[j+1].append(i+1)
				if not((str(i+1)+'-'+str(j+1)) or (str(j+1) +'-'+str(i+1))) in key:
					bonds[str(i+1)+'-'+str(j+1)] = [L ,'({a}-{b})'.format(a = atoms[i][0] , b = atoms[j][0])]

atoms = list()
n , molecule_name = read(sys.argv[1])
atomic_graph = dict([(i+1,[atoms[i][0]]) for i in range(n)])
bonds = dict()


surf(atoms)
print_original_info()
print_bond_graph()
print_bonds()