'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import os
import sys
import random
import locale
import numpy
import sasmol.sascalc as sascalc
import sassie.simulate.energy.dihedral_energy as dihedral_energy
import sassie.simulate.monomer_monte_carlo.dihedral_rotate as dihedral_rotate

#	STEP
#
#	01/18/2011	--	initial coding			:	jc
#
# LC	 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	STEP the module that contains the base classes that 
	handle setup and analysis routines for each Monte Carlo, ahem, step.
	
'''


class Setup():

    def __init__(self):
        pass

    def chooser(self, coor, m1, pairdat, dtheta, numranges, reslow, numcont, dihedral_parameters, beta, residue_rotation_indices, residue_rotation_mask, nonbondflag, first_last_resid, molecule_type, seed_object):
        '''
        This method is responsible select a range, a residue, get
        appropriate masks & indices, theta parameters, and choose
        an angle to vary.
        '''

        # pick a range
        if(seed_object != -1):
            ran_num = seed_object.rand()
            idum = int((numranges) * ran_num)
        else:
            idum = int((numranges) * random.random())
        # set the low range and number of contiguous for this range
        low = reslow[idum]
        num = numcont[idum]
        # pick a residue from this range
        q0a = numpy.arange(numcont[idum]) + low
        if(seed_object != -1):
            ran_num = seed_object.rand()
            ndum = int(len(q0a) * ran_num)
        else:
            ndum = int(len(q0a) * random.random())
        q0 = q0a[ndum]

        # get indices for this residue
        indices = residue_rotation_indices[q0]
        # get rotation mask for this residue
        this_mask = numpy.array(residue_rotation_mask[q0])

        # this max dtheta for this region
        this_dtheta = dtheta[idum]

        # pick a random angle for the list of dihedral angles for this molecule
        # type
        number_of_dihedrals = len(dihedral_parameters)

        search = True

        if(seed_object != -1):
            ran_num = seed_object.rand()
            if(molecule_type == 'rna'):
                # chose an angle type that does not include the delta dihedral
                # angle
                while(search):
                    this_dihedral = int(number_of_dihedrals * ran_num)
                    if(this_dihedral != 3):
                        if(q0 == first_last_resid[0]):
                            if(this_dihedral != 0):
                                search = False
                        elif(q0 == first_last_resid[1]):
                            if(this_dihedral != 4 and this_dihedral != 5):
                                search = False
                        else:
                            search = False
            else:
                this_dihedral = int(number_of_dihedrals * ran_num)
        else:
            if(molecule_type == 'rna'):
                # chose an angle type that does not include the delta dihedral
                # angle
                while(search):
                    this_dihedral = int(number_of_dihedrals * random.random())
                    if(this_dihedral != 3):
                        if(q0 == first_last_resid[0]):
                            if(this_dihedral != 0):
                                search = False
                        elif(q0 == first_last_resid[1]):
                            if(this_dihedral != 4 and this_dihedral != 5):
                                search = False
                        else:
                            search = False
            else:
                this_dihedral = int(number_of_dihedrals * random.random())

        # set specific dihedral parameters in this range (idum) for this
        # residue (ndum)
        if(molecule_type == 'protein'):
            resphi = dihedral_parameters[0]
            respsi = dihedral_parameters[1]
            parm = []
            parm = [resphi[idum][1][ndum], respsi[idum][1][ndum]]

            if(q0 == first_last_resid[0]):
                an = 'psi'
            elif(q0 == first_last_resid[1]):
                an = 'phi'
            elif(this_dihedral == 0):
                an = 'phi'
            else:
                an = 'psi'

        elif(molecule_type == 'rna'):
            resalpha = dihedral_parameters[0]
            resbeta = dihedral_parameters[1]
            resgamma = dihedral_parameters[2]
            resdelta = dihedral_parameters[3]
            resepsilon = dihedral_parameters[4]
            reseta = dihedral_parameters[5]
            parm = []
            parm = [resalpha[idum][1][ndum], resbeta[idum][1][ndum], resgamma[idum][1][
                ndum], resdelta[idum][1][ndum], resepsilon[idum][1][ndum], reseta[idum][1][ndum]]

            rna_dihedrals = ['alpha', 'beta',
                             'gamma', 'delta', 'epsilon', 'eta']
            an = rna_dihedrals[this_dihedral]

        vdi, vdf, trial_theta = self.find_angle(coor, an, indices, this_dihedral, q0, this_dtheta,
                                                beta, parm, this_mask, nonbondflag, first_last_resid, molecule_type, seed_object)

        pairdat[0] = an
        pairdat[1] = q0
        pairdat[2] = trial_theta

        return vdi, vdf, indices, this_mask

    def find_angle(self, coor, an, indices, this_dihedral, q0, this_dtheta, beta, parm, this_mask, nonbondflag, first_last_resid, molecule_type, seed_object):
        '''
        Method to search for a given dihedral angle which is chosen such that the
        angle has an valid dihedral energy.

        '''

       # 	angle_value=sascalc.Measure.calculate_dihedral(coor,this_sub_mask)
        angle_value = dihedral_rotate.measure(
            coor, indices, an, this_mask, q0, first_last_resid, molecule_type)
        search = 1
        while(search == 1):
            if(seed_object != -1):
                ran_num = seed_object.rand()
                trial_theta = this_dtheta * (2.0 * ran_num - 1)

            else:
                trial_theta = this_dtheta * (2.0 * random.random() - 1)
            try:
                theta = self.check_angle(angle_value, trial_theta)
            except:
                print('BUG IN STEP ANGLE')
                print('this_theta = ', this_theta)
                print('trial_theta = ', trial_theta)
                print('>>> setting theta to 0.0')
                theta = 0.0

            search, vdi, vdf = dihedral_energy.calc(
                this_dihedral, angle_value, theta, parm, beta, nonbondflag, seed_object)

        return vdi, vdf, trial_theta

    def check_angle(self, angle_value, trial_theta):
        '''
        Method checks the trial_theta such that it is between -180 and 180 degrees

        '''
        sum_theta = angle_value + trial_theta

        while sum_theta >= 180.0:
            # print 'sum_theta > 180.0: sum_theta = ',sum_theta
            sum_theta -= 360.0
            # print 'sum_theta > 180.0: sum_theta = ',sum_theta
        while sum_theta < -180.0:
            # print 'sum_theta < 0.0: sum_theta = ',sum_theta
            sum_theta += 360.0
            # print 'sum_theta < 0.0: sum_theta = ',sum_theta

#        if (sum_theta < 0.0):
#            #print 'sum_theta < 0.0: sum_theta = ',sum_theta
#            sum_theta = abs((sum_theta % 180.0))
#            #print 'sum_theta < 0.0: sum_theta = ',sum_theta
#

#        outfile = open("check_angle.txt","a")
#        outfile.write("%f\n" % (sum_theta))
#        outfile.close()

        return sum_theta

