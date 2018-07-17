#!/usr/bin/env python
'''
Driver that allows for the calculation of potential energies of water systems 
provided in an xyz-file
'''
from ase import Atoms
from mbpol_calculator import *
import sys
import argparse
import ase.io as io
import ipyparallel as ipp
import time
import pandas as pd

if __name__ == '__main__':

    # Parse Command line
    parser = argparse.ArgumentParser(description='Use Mbpol to calculate energies')
    parser.add_argument('xyz_file', action="store")
    parser.add_argument('out_file', action="store")
    parser.add_argument('-a', metavar='a', type=float, nargs = '?',
     default=0.0, help='Sidelength of cubic box')
    parser.add_argument('-pbc', action=('store_true'), help='Use PBC?')
    parser.add_argument('-ipprofile', metavar= 'ipprofile', type=str, nargs='?',
     default='None',help='Which ipcluster profile to use (default: None)')
    args = parser.parse_args()

    start_time = time.time()
    print('Reading from file {}'.format(args.xyz_file))

    with open(args.xyz_file, 'r') as xyz:
        n_atoms = int(xyz.readline())
        n_mol = int(n_atoms / 3)
        if args.a == 0.0:
            try:
                aline = xyz.readline().split()
                if len(aline) == 3:
                    a, b, c = aline[0], aline[1], aline[2]
                elif len(aline) == 1:
                    a = b = c = aline[0]
                else:
                    raise TypeError
            except TypeError:
                print('type Error')
                a = b = c =  0.0
        else:
            a = b = c =  args.a

    h2o = Atoms('{}OHH'.format(n_mol),
                positions = io.read(args.xyz_file, index = 0).get_positions(),
                cell = [a, b, c],
                pbc = args.pbc)

    reconnect_monomers(h2o)
    reconnect_monomers(h2o)
    h2o.calc = MbpolCalculator(h2o)

    # Find number of systems
    n_systems = len(io.read(args.xyz_file, index = ':'))

    # Find mpi settings
    indices = list(range(n_systems))
    pot_energies = []
    for i in indices:
        print(i)
        h2o.set_positions(io.read(args.xyz_file, index = i).get_positions())
        reconnect_monomers(h2o)
        reconnect_monomers(h2o)
        pot_energies.append(h2o.get_potential_energy())           

    pd.DataFrame(pot_energies).to_csv(args.out_file, index = None, header = None)
    print('Calculation took {} s'.format(time.time()-start_time))
    print(pot_energies)
