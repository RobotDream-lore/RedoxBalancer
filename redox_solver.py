import numpy as np
import re
from collections import OrderedDict 
import math
regex2 = r"(?P<molecule_c>^\d)"
regex = r"((?P<atom>[A-Z]{1}[a-z]*)(?P<atom_c>\d*))"
r = "H2 + O2 -> H2O"
rr = "Na2CO3"

class Molecule():
    def __init__(self,s) -> None:
        self.input_molecule_string = s
        print(self.input_molecule_string)
        self.atoms = []
        self.__pre_filter()
        self.out_coeff = 1
    
    def __pre_filter(self):
        self.molecule_coeff = re.findall(regex2,self.input_molecule_string)
        
        try:
            self.minimum_molecule_formula = self.input_molecule_string[len(str(self.molecule_coeff[0])):]
        except:
            self.minimum_molecule_formula = self.input_molecule_string
        
        print(self.minimum_molecule_formula)
        
        for _,atom,n in re.findall(regex,self.input_molecule_string):
            if n:
                pass
            else:
                n = 1
            self.atoms.append((atom,n))
            
        self.atoms = dict(self.atoms)

    def return_output(self):
        if self.out_coeff == 1:
            return self.minimum_molecule_formula
        else:
            o = str(int(self.out_coeff))
            o += self.minimum_molecule_formula
            return o

class Factor():
    
    def __init__(self,s) -> None:
        self.factor_string = s
        self.total_atoms = {}
        self.__pre_filter()
        
    def __pre_filter(self):
        self.molecules = list(map(lambda x: Molecule(x),list(map(lambda x: x.strip(),self.factor_string.split('+')))))
        self.__compute_total_atoms()

    def __compute_total_atoms(self):
        for m in self.molecules:
            for k,v in m.atoms.items():
                if k in self.total_atoms:
                    self.total_atoms[k] += int(v)
                else:
                    self.total_atoms[k] = int(v)
        self.total_atoms = dict(OrderedDict(sorted(self.total_atoms.items())))
        return self.total_atoms
    
    def return_output(self):
        o = ""
        for i,m in enumerate(self.molecules):
            o += m.return_output()
            if i < len(self.molecules)-1:
                o += ' + '
        return o


class Redox():
    def __init__(self,s) -> None:
        self.input_redox_string = s
        self.factors_string = []
        self.factors = []
        self.bilanced_costants = []
        self.__pre_filter()
    
    def __pre_filter(self):
        #Crea una lista con 2 elementi ogni elemento è una lista n con n numero di molecole nel fattore di sinistra o destra
        for i in self.input_redox_string.split('->'):
            self.factors.append(Factor(i))
        self.__create_matrix()
    
    def __solver_equations(self,U):
        # find the eigenvalues and eigenvector of U(transpose).U
        e_vals, e_vecs = np.linalg.eig(np.dot(U.T, U))  
        # extract the eigenvector (column) associated with the minimum eigenvalue
        return e_vecs[:, np.argmin(e_vals)] 

    def solve(self):
        solution = self.__solver_equations(self.matrix)
        s = solution / np.amin(solution)
        o = s
        for i in range(2,10):
            
            f = np.array(list(map(lambda x: math.modf(float(x)), o)))
            print(o)
            v = np.array(list(map(lambda x: True if (abs(x[0]) < 0.1 or abs(x[0]) > 0.9) else False ,f)))
            print(v)
            # if (o<0.98).any() or not((1 - np.absolute(f) < 0.01).all()):
            #     o = np.array(s*i)
            #     # print(o)
            #     # print(o<0.98)
            if not v.all():
                o = np.array(s*i)
        
        o = list(map(round,o))
        #print(o)
        self.bilanced_costants = o
        self.update_costants()

    
    def update_costants(self):
        for i,f in enumerate(self.factors):
            for j,m in enumerate(f.molecules):
                self.factors[i].molecules[j].out_coeff = self.bilanced_costants[i*len(self.factors[0].molecules) + j]
    
    def __create_matrix(self):
        matrix = []
        
        for j,f in enumerate(self.factors):
            print(f.total_atoms)
            for i,m in enumerate(f.molecules):
                tmp = []
                for k,v in f.total_atoms.items():
                    #print(k,v)
                    if k in m.atoms:
                        if j == 0:
                            tmp.append(int(m.atoms[k]))
                        else:
                            tmp.append(-int(m.atoms[k]))
                    else:
                        tmp.append(0)
                
                matrix.append(tmp)
        self.matrix = np.array(matrix).T
    
    def return_ouput(self):
        o = ""
        for i,f in enumerate(self.factors):
            o += f.return_output()
            if i < len(self.factors) - 1:
                o += ' -> '
        return o
       
            
redox = Redox(r)
redox.solve()
print(redox.return_ouput())