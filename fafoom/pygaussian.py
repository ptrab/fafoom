#   Copyright 2015 Adriana Supady
#
#   This file is part of fafoom.
#
#   Fafoom is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   Fafoom is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with fafoom.  If not, see <http://www.gnu.org/licenses/>.

'''Wrapper for Gaussian'''
from __future__ import division
import glob
import os
import subprocess

from utilities import sdf2xyz
from utilities import pse

hartree2eV = 27.21138602

class GaussianObject():
    '''Create and handle Gaussian objects.'''
    def __init__(self, commandline, memory="256mb", chargemult="0 1", nprocs=1):
        """Initialize the GaussianObject.

        Args(required):
            commandline
        Args(optional):
            memory (defaul="256mb")
            chargemult (default="0 1")
            nprocs (default=1)
        Raises:
            KeyError: if the commandline is not defined
        """
        self.commandline = commandline
        self.memory = memory
        self.chargemult = chargemult
        self.nprocs = nprocs

    def generate_input(self, sdf_string):
        """Create input files for Gaussian.
        Args:
            sdf_string (str)
        """
        xyz_string = sdf2xyz(sdf_string)
        coord = xyz_string.split('\n')
        string1 = '%nprocshared='+str(self.nprocs)+"\n"
        string2 = '%mem='+str(self.memory)+"\n"
        string3 = '#p '+str(self.commandline)+"\n"
        #string4 = 'scf=(maxcycles=N,conver=M)'+"\n"
        string5 = "\n fafoom\n"
        string6 = "\n"+str(self.chargemult)+"\n"
        with open('gaussian_job.gjf', 'w') as f:
            f.write(string1)
            f.write(string2)
            f.write(string3)
            #f.write(string4)
            f.write(string5)
            f.write(string6)
            f.write('\n'.join(coord[2:]))
            f.write('\n\n')
        f.close()

    def run_gaussian(self):
        """Run gaussian and write output to 'gaussian_job.log'.
        The optimized geometry is written to 'orca_molecule.xyz'.

        Warning: this function uses subprocessing to invoke the run.
        The subprocess's shell is set to TRUE.
        Args:
            none
        Raises:
            OSError: if gaussian_job.gjf not present in the working directory
        """
        success = False
        if os.path.exists('gaussian_job.gjf') is False:
            raise OSError("Required input file not present.")
        gaussian = subprocess.Popen(str("g09 gaussian_job.gjf"),
            stdout=subprocess.PIPE, shell=True)
        out = subprocess.Popen(
            ['cat'], stdin=gaussian.stdout,
            stdout=open('gaussian_job.log', 'w'), shell=True)
        out.wait()

        searchfile = open("gaussian_job.log", "r")

        s0 = "Stationary point found"
        s = "SCF Done"
        not_conv = True
        for line in searchfile:
            if s0 in line:
                not_conv = False
        searchfile.close()
        if not_conv:
            killfile = open("kill.dat", "w")
            killfile.close()
        else:
            searchfile = open("gaussian_job.log", "r")
            for line in searchfile:
                if s in line:
                    print line
                    energy_tmp = float(line.split()[4].split('\n')[0])
                    print energy_tmp
            searchfile.close()
            self.energy = energy_tmp

            #parsing ugly gaussian output for the geometry
            with open("gaussian_job.log") as f:
                lines = f.readlines()
                string = 0
                for num, line in enumerate(lines, 1):
                    line = line.split('\n')[0]
                    #numbers of atoms
                    if 'NAtoms' in line:
                        natoms = line.split()[1]
                        xyzstring = natoms+"\n\n"
                    #find the position of the geometry
                    if 'Standard orientation' in line:
                        start=num+4
                        end=start+int(natoms)
                        geom = lines[start:end]
                #produce the xyz-"input"
                for line in geom:
                    line = line.split('\n')[0]
                    centernr, atomicnr, atype, xcoord, ycoord, zcoord = line.split()
                    xyzstring += pse(int(atomicnr))+' '+xcoord+' '+ycoord+' '+zcoord+'\n'

            self.xyz_string_opt = xyzstring
            success = True
        return success

    def get_energy(self):
        """Get the energy of the molecule.

        Returns:
            energy (float) in eV
        Raises:
            AttributeError: if energy hasn't been calculated yet
        """
        if not hasattr(self, 'energy'):
            raise AttributeError("The calculation wasn't performed yet.")
        else:
            return hartree2eV*self.energy

    def get_xyz_string_opt(self):
        """Get the optimized xyz string.

        Returns:
            optimized xyz string (str)
        Raises:
            AttributeError: if the optimization hasn't been performed yet
        """
        if not hasattr(self, 'xyz_string_opt'):
            raise AttributeError("The calculation wasn't performed yet.")
        else:
            return self.xyz_string_opt

    def clean(self):
        """Clean the working direction after the orca calculation has been
        completed.
        """
        for f in glob.glob("orca_molecule.*"):
            os.remove(f)
