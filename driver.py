from ase import Atoms
from mbpol_calculator import *

if __name__ == '__main__':

    a = 15.646
    h2o = Atoms('128OHH',
                positions = np.genfromtxt('./128.csv', delimiter = ','),
                cell = [a, a, a],
                pbc = True)

    # As provided, 128.csv contains monomers that are split up across unit cells
    # boundaries. MBPol/OpenMM cannot deal with that which is why reconnet_monomers()
    # has to be used. 
    reconnect_monomers(h2o)
    h2o.calc = MbpolCalculator(h2o)
    print('Epot[eV] = {}'.format(h2o.get_potential_energy()))
