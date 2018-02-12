import numpy as np
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol
from siesta_utils.conversions import kcaltoeV 
from simtk.openmm import LocalEnergyMinimizer


def coords_to_pdb(coords, n_mol):
    
    n_systems = int(len(coords)/3/n_mol)

    cnt = 0

    labels = ['O   ','H1  ','H2  ']
    labeloh = ['O', 'H', 'H']

    with open('tmp.pdb','w') as pdbfile:
        
        for system in range(n_systems):
            mol_cnt = 1
            atm_cnt = 1
            
            pdbfile.write('MODEL{:9d}\n'.format(system + 1))
            
            for molecule in range(n_mol):
                for l, loh in zip(labels,labeloh):
                    pdbfile.write('HETATM{:5d}  '.format(atm_cnt) + l + 'HOH A{:4d}    '.format(mol_cnt))
                    pdbfile.write('{:7.4f} {:7.4f} {:7.4f}  1.00  0.00'.format(coords[cnt,0],
                                                                              coords[cnt,1],
                                                                              coords[cnt,2]) + ' '*11 + loh + '  \n')
                    cnt += 1
                    atm_cnt += 1
                    
                pdbfile.write('HETATM{:5d}  '.format(atm_cnt) + 'M   ' + 'HOH A{:4d}    '.format(mol_cnt))
                pdbfile.write('{:7.4f} {:7.4f} {:7.4f}  1.00  0.00 \t \t \n'.format(0,0,0))
                atm_cnt += 1
                mol_cnt += 1

            pdbfile.write('TER   {:5d}  '.format(atm_cnt) + 'M   ' + 'HOH A{:4d}    \n'.format(mol_cnt))
            pdbfile.write('ENDMDL \n')

    return 'tmp.pdb'
