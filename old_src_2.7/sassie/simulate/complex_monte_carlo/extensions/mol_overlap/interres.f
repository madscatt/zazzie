	subroutine mover(c1,na1,c2,na2,cut,np2,atom_list)
c	subroutine mover(c1,na1,c2,na2,cut,r,ni)
	double precision c1(na1,3)
	double precision c2(na2,3)
	double precision cut
	integer na1,na2,ni,npairs,i,j,k,count,qex,np2
c	integer r(ni,2)
	integer atom_list(np2)
	double precision x1,y1,z1,x2,y2,z2,dist
	double precision dx2,dy2,dz2

cf2py	intent(in) :: c1,c2,cut,atom_list
cf2py	intent(out):: atom_list
cf2py	intent(hide):: na1,na2,ni,npairs,count,i,j,k
cf2py	intent(hide):: x1,y1,z1,x2,y2,z2,dist,qex,np2

c	npairs = (na1+na2)*((na1+na2)-1)/2

c	do 30 i=1,npairs
c		atom_list(i) = 0
c  30  continue

	count = 1
	do 200,i=1,na1
		x1 = c1(i,0)	
		y1 = c1(i,1)	
		z1 = c1(i,2)	
		do 100,j=1,na2
			x2 = c2(j,0)	
			y2 = c2(j,1)	
			z2 = c2(j,2)
			qex=0	
c			do 50,k=1,ni
c        1         2         3         4         5           
c2345678901234567890123456789012345678901234567890123456789 
c			  if((i.eq.r(k,0)).and.(j.eq.r(k,1))) then
c              	    qex=1	
c                    endif
c  50			continue
			if(qex.eq.0) then
				dx2 = (x2-x1)*(x2-x1)
				dy2 = (y2-y1)*(y2-y1)
				dz2 = (z2-z1)*(z2-z1)
				dist = sqrt(dx2+dy2+dz2)
				if(dist.lt.cut) then
					atom_list(count) = 1	
					count = count + 1 
				endif
			endif
  100 	continue
  200	continue

	write(6,*) 'FORTAN: COUNT = ',count

	end
