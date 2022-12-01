
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012
C
        subroutine inter(coor1,coor2,charge1,charge2,vdwp1,vdwp2,
     cswitchd,nbcutoff,slength,er,boxl,sflag,vdwflag,natoms1,natoms2,
     cnonbondenergy)
      
        integer natoms1,natoms2,i,j,np,sflag,vdwflag
        double precision coor1(natoms1,3) 
        double precision coor2(natoms2,3) 
        double precision charge1(natoms1)
        double precision charge2(natoms2)
        double precision vdwp1(natoms1,2)
        double precision vdwp2(natoms2,2)

        double precision rij,pi,qconv,qconv2,eps,qi,qj,er
        double precision epsi,epsj,rmi,rmj,epsv,rminij,slength
        double precision ang2m,jtokjpmol,kjtokcal,conv
        double precision x2,y2,z2,dx2,dy2,dz2,switchscale,vij
        double precision elenergy,vdwenergy,switchd,nbcutoff
        double precision arg1,arg2,arg3,rij6,rij12,boxl,invboxl
        double precision nonbondenergy

C        write(*,*) 'fortran: switchd = ',switchd
C        write(*,*) 'fortran: nbcutoff = ',nbcutoff

cf2py intent(in) :: coor1,charge1,natoms1,switchd,nbcutoff,er,slength
cf2py intent(in) :: coor2,charge2,natoms2,vdwp1,vdwp2,boxl,sflag
cf2py intent(in) :: vdwflag
cf2py intent(out):: nonbondenergy
cf2py intent(hide):: i,j,pi,qconv,qconv2,eps,ang2m,jtokjpmol,kjtokcal
cf2py intent(hide):: conv,np,x2,y2,z2,dx2,dy2,dz2,rij,switchscale,vij
cf2py intent(hide):: arg1,arg2,arg3,rij6,rij12,q1,q2,invboxl
cf2py intent(hide):: epsi,epsj,rmi,rmj,epsv,rminij

C eV --> C
        qconv=1.620177E-19
C q*q --> C^2
        qconv2=qconv*qconv
C epsilon --> C^2/(J*m)
        eps=8.854187816E-12
C angstrom --> m
        ang2m=1E-10
C J --> KJ/mol
        jtokjpmol=6.022137E+20
C KJ --> kcal
        kjtokcal=1d0/4.184

C qconv2=1.0

        pi=dacos(-1d0)
        conv=(qconv2/(4.0*pi*eps*ang2m))*jtokjpmol*kjtokcal
C
        conv=332.4683

        np=0

C  qi*qj/(er*rij)

        invboxl = 1d0/boxl

        elenergy = 0d0
        vdwenergy = 0d0
        switchscale = 1d0
        nonbondingenergy = 0d0

        do 200,i=1,natoms1
            qi=charge1(i)
            epsi = vdwp1(i,1) 
            rmi = vdwp1(i,2) 
            x1=coor1(i,1)
            y1=coor1(i,2)
            z1=coor1(i,3)
            do 100,j=1,natoms2
               qj=charge2(j)
               epsj = vdwp2(j,1) 
               rmj = vdwp2(j,2) 
               x2=coor2(j,1)
               y2=coor2(j,2)
               z2=coor2(j,3)

               dx = x2-x1
               dy = y2-y1
               dz = z2-z1

               dx = dx - boxl*anint(dx*invboxl)
               dy = dy - boxl*anint(dy*invboxl)
               dz = dz - boxl*anint(dz*invboxl)

               dx2=dx*dx
               dy2=dy*dy
               dz2=dz*dz

               rij=sqrt(dx2+dy2+dz2)

               if(rij . LT. 1E-2) then
C                  write(*,*) 'rij is screwed = ',rij
C                  call flush()
                  exit
               endif

               epsv = sqrt(epsi*epsj)
               rminij = rmi+rmj

               if(rij . LE . switchd) then
                  switchscale=1d0
               else if(rij .GT. switchd .AND. rij .LE. nbcutoff) then 
                  arg1 = (nbcutoff-rij)*(nbcutoff-rij)
                  arg2 = (nbcutoff+2.0*rij-3.0*switchd)
                  arg3 = ((nbcutoff-switchd)**3.0)
                  switchscale = (arg1*arg2)/arg3
               else if (rij .GT. nbcutoff) then
                  switchscale=0d0
               end if

               if(sflag .EQ. 0) then
                  vij=((qi*qj)/(er*rij))*switchscale
               else if(slag .EQ. 1) then
                  switchscale = exp(-rij/slength)
                  vij=((qi*qj)/(er*rij))*switchscale
               end if

               elenergy=elenergy+vij
             
C               write(*,*) 'el = ',elenergy
 
               rij6 = (rminij/rij)**6.0
               rij12 = rij6 * rij6
               
               if(vdwflag .EQ. 1) then
                  vij=epsv*(rij12-2.0*rij6)*switchscale
               else
                  vij=epsv*(rij12)*switchscale
               endif

               vdwenergy=vdwenergy+vij
C               write(*,*) 'vdw = ',vdwenergy

               np=np+1

  100   continue
  200   continue

        elenergy = elenergy*conv

        nonbondenergy = elenergy + vdwenergy

        end


