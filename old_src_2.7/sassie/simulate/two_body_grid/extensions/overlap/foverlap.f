        subroutine foverlap(coor1,coor2,natoms1,natoms2,cut,check)
        double precision coor1(natoms1,3)
        double precision coor2(natoms2,3)
        double precision cut
        integer natoms1,natoms2,check
        double precision x1,y1,z1,x2,y2,z2,diff2,dist

cf2py	intent(in) :: coor1,coor2,cut
cf2py	intent(out):: check
cf2py	intent(hide)::natoms1,natoms2
cf2py	intent(hide)::x1,y1,z1,x2,y2,z2,diff2,dist

c        check = 0

        do 200,i=1,natoms1
           x1=coor1(i,1)
           y1=coor1(i,2)
           z1=coor1(i,3)
           do 100,j=1,natoms2
             x2=coor2(j,1)
             y2=coor2(j,2)
             z2=coor2(j,3)
             diff2=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
             dist=sqrt(diff2)
             if(dist.lt.cut) then
              check=1
              exit
             endif
  100        continue

             if(check==1) then
              exit
             endif
  200   continue

        end
