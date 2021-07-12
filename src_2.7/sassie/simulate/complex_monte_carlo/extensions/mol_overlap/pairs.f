	subroutine pairs(avdw,natoms,npairs,cutoffarray)
	double precision avdw(natoms,2)
	double precision cutoffarray(npairs)
	integer natoms,npairs,count,i,j
	double precision r1,r2

cf2py	intent(in) :: avdw,cutoffarray
cf2py	intent(out):: cutoffarray
cf2py	intent(hide)::natoms,npairs
cf2py	intent(hide)::r1,r2

	count = 1
	do 200,i=1,natoms-1
		r1 = avdw(i,2)   	   
		do 100,j=i+1,natoms
			r2 = avdw(j,2)
			if (count .eq. 1) then
				write(6,*) (r1+r2)
			end if
			cutoffarray(count) = (r1+r2)
			count = count + 1 
  100 	continue
  200	continue

	end
