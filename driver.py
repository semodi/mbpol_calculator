from ase import Atoms
from mbpol_calculator import *

if __name__ == '__main__':

    a = 15.646
    h2o = Atoms('128OHH',
                positions = np.genfromtxt('./128.csv', delimiter = ','),
                cell = [a, a, a],
                pbc = True)

    reconnect_monomers(h2o)
    h2o.calc = MbpolCalculator(h2o)
    print(h2o.get_potential_energy())
