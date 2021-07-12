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
import os,string,sys

#       WRITE_APBS_INPUT
#
#       10/1/2007       --      initial coding                  :       jc
#       09/16/2012      --      adapated from namd input 	:       jc
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
        WRITE_APBS_INPUT is the function to write a generic input
	file for the program APBS (version 1.8).

	This module is called by APBS.PY

	Writes input file for vacuum and solvent polar electrostatics 

	REFERENCE:

      Baker NA, Sept D, Joseph S, Holst MJ, McCammon JA. Electrostatics of 
	nanosystems: application to microtubules and the ribosome. 
	Proc. Natl. Acad. Sci. USA 98, 10037-10041 2001. 

	M. Holst and F. Saied, Multigrid solution of the Poisson-Boltzmann equation. 
	J. Comput. Chem. 14, 105-113, 1993.

	M. Holst and F. Saied, Numerical solution of the nonlinear Poisson-Boltzmann 
	equation: Developing more robust and efficient methods. 
	J. Comput. Chem. 16, 337-364, 1995.

	M. Holst, Adaptive numerical treatment of elliptic systems on manifolds. 
	Advances in Computational Mathematics 15, 139-191, 2001. 

	R. Bank and M. Holst, A New Paradigm for Parallel Adaptive Meshing Algorithms. 
	SIAM Review 45, 291-323, 2003.

'''

def write_apbs_input(maximum_dimensions,temperature,inputfilename,ion_charge,ion_conc,ion_radius):

        '''
        WRITE_APBS_INPUT is the function to write a generic input
	file for the program APBS (version 1.8).
       

        INPUT:  variable descriptions:
      
		min_max		: list of minimum and maximum dimensions
		temperature	:	temperature

        OUTPUT:
                inputfilename:          name of input file that function will write

        '''

	max_x = maximum_dimensions[0]
	max_y = maximum_dimensions[1]
	max_z = maximum_dimensions[2]

	outfile=open(inputfilename,'w')

	
#	outfile.write('%s\n' % ('exclude                 scaled1-4'))

	outfile.write('%s\n' % ('read'))
	outfile.write('\t%s\n' % ('mol pqr junk.pqr'))
	outfile.write('%s\n' % ('end'))

	outfile.write('%s\n' % ('elec'))
	outfile.write('\t%s\n' % ('mg-auto'))
	outfile.write('\t%s\n' % ('dime 353 193 289'))  ### this needs to be variable

	outfile.write('\t%s\t%f\t%f\t%f\n' % ('cglen', max_x*1.70,max_y*1.70,max_z*1.70)) ### why 170%
	outfile.write('\t%s\t%f\t%f\t%f\n' % ('fglen', max_x*1.30,max_y*1.30,max_z*1.30)) ### why 130%

	outfile.write('\t%s\n' % ('cgcent mol 1'))
	outfile.write('\t%s\n' % ('fgcent mol 1'))
	outfile.write('\t%s\n' % ('mol 1'))
	outfile.write('\t%s\n' % ('lpbe'))
	outfile.write('\t%s\n' % ('bcfl sdh'))

	outfile.write('\t%s\t%f\t%f\t%f\n' % ('ion',ion_charge,ion_conc,ion_radius))
	outfile.write('\t%s\t%f\t%f\t%f\n' % ('ion',-ion_charge,ion_conc,ion_radius))

	outfile.write('\t%s\n' % ('pdie 2.0000'))
	outfile.write('\t%s\n' % ('sdie 78.5400'))
	outfile.write('\t%s\n' % ('srfm smol'))
	outfile.write('\t%s\n' % ('chgm spl2'))
	outfile.write('\t%s\n' % ('sdens 10.00'))
	outfile.write('\t%s\n' % ('srad 1.40'))
	outfile.write('\t%s\n' % ('swin 0.30'))
	outfile.write('\t%s\t%f\n' % ('temp',temperature))
	outfile.write('\t%s\n' % ('calcenergy total'))
	outfile.write('\t%s\n' % ('calcforce no'))
	#outfile.write('\t%s\n' % ('write pot dx pot'))
	outfile.write('%s\n' % ('end'))
		
	outfile.write('%s\n' % ('elec'))
	outfile.write('\t%s\n' % ('mg-auto'))
	outfile.write('\t%s\n' % ('dime 353 193 289'))  ### this needs to be variable

	outfile.write('\t%s\t%f\t%f\t%f\n' % ('cglen', max_x*1.70,max_y*1.70,max_z*1.70))  ### why 170%
	outfile.write('\t%s\t%f\t%f\t%f\n' % ('fglen', max_x*1.30,max_y*1.30,max_z*1.30))  ### why 130%

	outfile.write('\t%s\n' % ('cgcent mol 1'))
	outfile.write('\t%s\n' % ('fgcent mol 1'))
	outfile.write('\t%s\n' % ('mol 1'))
	outfile.write('\t%s\n' % ('lpbe'))
	outfile.write('\t%s\n' % ('bcfl sdh'))

	outfile.write('\t%s\t%f\t%f\t%f\n' % ('ion',ion_charge,ion_conc,ion_radius))
	outfile.write('\t%s\t%f\t%f\t%f\n' % ('ion',-ion_charge,ion_conc,ion_radius))

	outfile.write('\t%s\n' % ('pdie 2.0000'))
	outfile.write('\t%s\n' % ('sdie 2.0000'))
	outfile.write('\t%s\n' % ('srfm smol'))
	outfile.write('\t%s\n' % ('chgm spl2'))
	outfile.write('\t%s\n' % ('sdens 10.00'))
	outfile.write('\t%s\n' % ('srad 1.40'))
	outfile.write('\t%s\n' % ('swin 0.30'))
	outfile.write('\t%s\t%f\n' % ('temp',temperature))
	outfile.write('\t%s\n' % ('calcenergy total'))
	outfile.write('\t%s\n' % ('calcforce no'))
	#outfile.write('\t%s\n' % ('write pot dx pot'))
	outfile.write('%s\n' % ('end'))
		
	outfile.write('%s\n' % ('print elecEnergy 1 end'))    
	outfile.write('%s\n' % ('print elecEnergy 2 end'))    
	outfile.write('%s\n' % ('print elecEnergy 2 - 1 end'))    
	outfile.write('%s\n' % ('quit'))
    
	outfile.close()

	return 



if __name__=='__main__': 
        print 'running as a main process'

	import numpy
	min_max =  [numpy.array([-72.92417908, -44.44118118, -46.919384  ]), numpy.array([ 82.26298523,  30.79856491,  73.71658325])]
	maximum_dimensions = [min_max[1][0] - min_max[0][0],min_max[1][1] - min_max[0][1], min_max[1][2] - min_max[0][2]]
	print 'min_max = ',min_max
	print 'maximum_dimensions = ',maximum_dimensions
	temperature = 300.0
	inputfilename = 'dum.in'

	write_apbs_input(maximum_dimensions,temperature,inputfilename)
       
else:
        print 'running as a spawned process'

