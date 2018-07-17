from ase import Atoms
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol
import numpy as np
from ase import units as ase_units
from .utils import coords_to_pdb
import os
eVtokcal = 23.06035
kcaltoeV = 1/eVtokcal

kilocalorie_per_mole_per_angstrom = unit.kilocalorie_per_mole/unit.angstrom

class MbpolCalculator:

    def __init__(self, atoms):

        # Create a PDB file that can be read by OpenMM
        coords_to_pdb(atoms.get_positions(),int(atoms.get_number_of_atoms()/3))
        pdb = app.PDBFile("./tmp.pdb")
        os.remove('./tmp.pdb')

        cell = atoms.get_cell()

        # Make sure bond lengths are within MBPol range
        distances = atoms.get_all_distances()
        for i in range(0,len(distances),3):
            for j in [1,2]:
                if distances[i,i+j] > 2.5:
                    print('WARNING: rOH > 2.5 A!!!!!! ({})'.format(distances[i,i+j]))
                    if atoms.get_pbc()[0]: print('Try: reconnect_monomers()')

        self.atoms = atoms
        
        if np.count_nonzero(cell - np.diag(np.diag(cell))) != 0:
            raise Exception('Only orthorhombic unit cells supported')

        if np.alltrue(atoms.get_pbc()):
            self.nonbondedMethod = app.PME
            if np.allclose(np.diag(cell),np.array([0,0,0])):
                raise Exception('Specify unit cell for PBC')
            self.boxSize = tuple(np.diag(cell)) * unit.angstrom
            pdb.topology.setUnitCellDimensions(self.boxSize)
            self.nonbondedCutoff = cell[0,0] / 2 * unit.angstrom
        elif np.alltrue([not a for a in atoms.get_pbc()]):
            self.nonbondedMethod = app.CutoffNonPeriodic
            self.nonbondedCutoff = 1e4 * unit.angstrom
        else:
            raise Exception('Only non-PBC or PBC in all euclidean directions supported')

        self.forcefield = app.ForceField(mbpol.__file__.replace("mbpol.py", "mbpol.xml"))
        self.system = self.forcefield.createSystem(pdb.topology,
                                         nonbondedMethod=self.nonbondedMethod,
                                         nonbondedCutoff=self.nonbondedCutoff)

        self.platform = mm.Platform.getPlatformByName('Reference')
        integrator = mm.VerletIntegrator(0.5) # Never used but needed for context
        self.simulation = app.Simulation(pdb.topology, self.system, integrator, self.platform)
        self.set_positions(atoms) #ASE: Angstrom , OMM : nm
        self.state = self.simulation.context.getState(getForces=True, getEnergy=True)
        self.last_coordinates = np.array(pdb.positions.value_in_unit(unit.angstrom))
        self.last_coordinates = np.delete(self.last_coordinates,
                                          np.arange(3,len(self.last_coordinates),4), axis = 0).reshape(-1,3)

    def calculation_required(self, atoms, quantities = None):
        return (not np.all(self.last_coordinates == atoms.positions))

    def set_positions(self, atoms):
        # Convert positions as stored in ase Atoms instance to OpenMM/MBPol format
        # which includes a 4th position per molecule for the virtual site "M"
        pos = np.zeros([int(len(atoms.positions)/3*4),3]).reshape(-1,4,3)
        pos[:,:3,:] = atoms.positions.reshape(-1,3,3)
        self.simulation.context.setPositions(pos.reshape(-1,3)/10) #ASE: Angstrom , OMM : nm
        self.simulation.context.computeVirtualSites()

    def get_forces(self, atoms):
        if self.calculation_required(atoms):
            self.set_positions(atoms)
            self.state = self.simulation.context.getState(getForces=True, getEnergy=True)
            self.last_coordinates = np.array(atoms.positions)
        forces = np.array(self.state.getForces().value_in_unit(kilocalorie_per_mole_per_angstrom))*kcaltoeV
        # Delete forces on virtual sites
        return np.delete(forces, np.arange(3,len(forces),4), axis = 0)

    def get_potential_energy(self, atoms, force_consistent = False):
        if self.calculation_required(atoms):
            self.set_positions(atoms)
            self.state = self.simulation.context.getState(getEnergy=True, getForces=True)
            self.last_coordinates = np.array(atoms.positions)
        return self.state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)*kcaltoeV

    def get_stress(self, atoms):
        return np.zeros([3,3])
#        raise Exception('get_stress() not implemented')

def reconnect_monomers(atoms, rc = 2.5):
    """ Reconnect Hydrogen and Oxygen that are split across unit cell
        boundaries and are therefore situated on opposite ends of the unit
        cell
    """
    boxsize = np.diag(atoms.get_cell())
    pos0 = np.array(atoms.positions)

    for i,_ in enumerate(atoms.get_positions()[::3]):

        if atoms.get_distance(i*3,i*3+1) > rc:
            d = atoms.positions[i*3] - atoms.positions[i*3+1]
            which = np.where(np.abs(d) > 5)[0]
            for w in which:
                pos0[i*3+1, w] += d[w]/np.abs(d[w]) * boxsize[w]
        if atoms.get_distance(i*3,i*3+2) > rc:
            d = atoms.positions[i*3] - atoms.positions[i*3+2]
            which = np.where(np.abs(d) > 5)[0]
            for w in which:
                pos0[i*3+2, w] += d[w]/np.abs(d[w]) * boxsize[w]

    atoms.set_positions(pos0)
    return atoms
