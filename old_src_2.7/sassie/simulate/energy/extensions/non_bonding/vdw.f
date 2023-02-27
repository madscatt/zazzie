
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine fvdw(coor,vdwp,nbcutoff,nbmargin,natoms,vdwenergy)
        integer natoms,i,j
        double precision coor(natoms,3)
        double precision vdwp(natoms,2)
C
        double precision rij,epsi,epsj,rmi,rmj,eps,rminij
        double precision x1,y1,z1,x2,y2,z2,dx2,dy2,dz2,switchscale,vij
        double precision vdwenergy,nbcutoff,nbmargin
        double precision arg1,arg2,arg3,rij6,rij12

cf2py intent(in) :: coor,vdwp,natoms,nbcutoff,nbmargin
cf2py intent(out):: vdwenergy
cf2py intent(hide):: i,j,rij,epsi,epsj,rmi,rmj,eps,rminij
cf2py intent(hide):: arg1,arg2,arg3,rij6,rij12,switchscale,vij

        vdwenergy = 0d0
        switchscale = 1d0
        do 200,i=1,natoms
            epsi = vdwp(i,1)
            rmi = vdwp(i,2)
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,j=i+1,natoms
               epsj = vdwp(j,1)
               rmj = vdwp(j,2)
               x2=coor(j,1)
               y2=coor(j,2)
               z2=coor(j,3)

               dx2=(x1-x2)*(x1-x2)
               dy2=(y1-y2)*(y1-y2)
               dz2=(z1-z2)*(z1-z2)
               rij=sqrt(dx2+dy2+dz2)

               eps = sqrt(epsi*epsj)
               rminij = rmi+rmj

C call getswitch(rij,nbcutoff,nbmargin,switchscale)

C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

               if(rij . LE . nbcutoff) then
                  switchscale=1d0
               else if(rij .GT. nbcutoff .AND. rij .LE. nbmargin) then
                  arg1 = (nbmargin-rij)*(nbmargin-rij)
                  arg2 = (nbmargin+2.0*rij-3.0*nbcutoff)
                  arg3 = ((nbmargin-nbcutoff)**3.0)
                  switchscale = (arg1*arg2)/arg3
               else if (rij .GT. nbmargin) then
                  switchscale=0.0
               end if
               rij6 = (rminij/rij)**6.0
               rij12 = rij6 * rij6
               vij=eps*(rij12-2.0*rij6)*switchscale

               vdwenergy=vdwenergy+vij

  100   continue
  200   continue

        end


