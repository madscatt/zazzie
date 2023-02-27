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
import string,locale,numpy,math
import sasmol.sasmath as sasmath

def angle(coor,angleparam):
	angleenergy=0.0 ; ubenergy=0.0 
	deg2rad=numpy.pi/180.0

	for i in xrange(len(angleparam)):
		vub=0.0	
		ttrio=angleparam[i][0]	
		atm1=locale.atoi(ttrio[0])	
		atm2=locale.atoi(ttrio[1])	
		atm3=locale.atoi(ttrio[2])	
		coor1=coor[atm1]; coor2=coor[atm2] ; coor3=coor[atm3]
		theta=sasmath.calc_angle(coor1,coor2,coor3)
		ktheta=locale.atof(angleparam[i][1][0])	
		theta0=deg2rad*locale.atof(angleparam[i][1][1])	

		if(len(angleparam[i][1])>2):
			kub=locale.atof(angleparam[i][1][2])
			so=locale.atof(angleparam[i][1][3])
			v = coor3-coor1
			s = math.sqrt(sum(v*v))
			#s=sascalc.calc_dist(coor1,coor3)
			vub=kub*(s-so)**2.0
	
		vangle=ktheta*(theta-theta0)**2.0
		angleenergy=angleenergy+vangle
		ubenergy=ubenergy+vub

	return angleenergy,ubenergy

