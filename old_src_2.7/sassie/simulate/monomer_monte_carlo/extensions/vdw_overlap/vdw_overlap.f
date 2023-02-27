        subroutine overlap(coor1,natoms1,cutoffarray,npairs,fact,check)
        double precision coor1(natoms1,3)
        double precision cutoffarray(npairs),fact
        integer natoms1,check,count
        double precision x1,y1,z1,x2,y2,z2,diff2,dist

cf2py	intent(in) :: coor1,cutoffarray,fact
cf2py	intent(out):: check
cf2py	intent(hide)::natoms1,npairs
cf2py	intent(hide)::x1,y1,z1,x2,y2,z2,diff2,dist

        count = 1
	  check = 0
        do 200,i=1,natoms1-1
           x1=coor1(i,1)
           y1=coor1(i,2)
           z1=coor1(i,3)
           do 100,j=i+1,natoms1
             x2=coor1(j,1)
             y2=coor1(j,2)
             z2=coor1(j,3)
             diff2=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
             dist=sqrt(diff2)
             if(dist.lt.(fact*cutoffarray(count))) then
              check=1
              exit
             endif
             count = count + 1

  100        continue

             if(check==1) then
              exit
             endif
  200   continue

        end


