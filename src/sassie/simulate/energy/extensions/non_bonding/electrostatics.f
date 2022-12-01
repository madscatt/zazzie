
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine fel(coor,charge,nbcutoff,nbmargin,natoms,elenergy)
        integer natoms,i,j,np
        double precision coor(natoms,3) 
        double precision charge(natoms)
        double precision rij,pi,qconv,qconv2,eps,e1,qi,qj
        double precision ang2m,jtokjpmol,kjtokcal,conv
        double precision x2,y2,z2,dx2,dy2,dz2,switchscale,vij
        double precision elenergy,nbcutoff,nbmargin
        double precision arg1,arg2,arg3
        
cf2py intent(in) :: coor,charge,natoms,nbcutoff,nbmargin
cf2py intent(out):: elenergy
cf2py intent(hide):: i,j,pi,qconv,qconv2,eps,ang2m,jtokjpmol,kjtokcal
cf2py intent(hide):: conv,np,x2,y2,z2,dx2,dy2,dz2,rij,switchscale,vij
cf2py intent(hide):: arg1,arg2,arg3


C  to call this from python
C       
C    import sys ; sys.path.append('./')
C    import electrostatics
C
C    elenergy = elecrostatics.fel(coor,charge,nbcutoff,nbmargin,natoms)
C
C    sudo python setup_electrostatics.py build
C    cp build/lib*/electrostatics.o (or .so) to directory that you are 
C    running in
C

C  qi*qj/(e1*rij)
        elenergy = 0.0
        e1 = 1d0

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

        write(*,*) 'in fortran: conv = ',conv

        do 200,i=1,natoms
            qi=charge(i)
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,j=i+1,natoms
               qj=charge(j)
               x2=coor(j,1)
               y2=coor(j,2)
               z2=coor(j,3)

               dx2=(x1-x2)*(x1-x2)
               dy2=(y1-y2)*(y1-y2)
               dz2=(z1-z2)*(z1-z2)
               rij=sqrt(dx2+dy2+dz2)
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
               vij=((qi*qj)/(e1*rij))*switchscale

               elenergy=elenergy+vij
               np=np+1
  100   continue
  200   continue

        elenergy = elenergy*conv

        end


