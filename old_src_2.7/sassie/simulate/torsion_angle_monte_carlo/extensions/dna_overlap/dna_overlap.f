C
C Author:  Steven C. Howell
C Purpose: Calculate interaction energy between every pair of atoms
C Created: January 2014
C
C $Id: dna_overlap.f 3249 2016-06-27 03:44:03Z schowell $
C
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine wca_nbyn(coor,r0,rN,c0,cN,w,wca,n)
        integer n,i,j,r0,rN,c0,cN
        double precision coor(n,3)
        double precision wca(n,n)
        double precision cutoff,w,rij
        double precision x2,y2,z2,dx2,dy2,dz2
        double precision x1,y1,z1

cf2py intent(in) :: coor,r0,rN,c0,cN,w,wca,n
cf2py intent(in,out):: wca
cf2py intent(hide):: i,j,cutoff,rij
cf2py intent(hide):: x2,y2,z2,dx2,dy2,dz2,rij
cf2py intent(hide):: x1,y1,z1

C    to call this from python:
C
C    import sys ; sys.path.append('./')
C    import electrostatics
C
C    elenergy = collision.wca_nbyn(coor,charge,nbcutoff,nbmargin,n)
C
C
C    to setup and incorporate into python:
C
C    sudo python setup_collision.py build
C    cp build/lib*/collision.o (or .so) to directory that you are
C    running in

        cutoff = 2.**(1./6.)*w

        do 200,i=r0,rN
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,j=c0,cN
               x2=coor(j,1)
               y2=coor(j,2)
               z2=coor(j,3)

               dx2=(x1-x2)*(x1-x2)
               dy2=(y1-y2)*(y1-y2)
               dz2=(z1-z2)*(z1-z2)
               rij=sqrt(dx2+dy2+dz2)

               if(rij . LT . cutoff) then
                  wca(i,j)=(w/rij)**12.-(w/rij)**6.+0.25
               else
                  wca(i,j)=0
               end if
  100   continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine wca_nbym(coor_r,coor_c,r0,rN,c0,cN,w,wca,n,m)
        integer n,m,c0,cN,r0,rN,i,j
        double precision coor_r(n,3),coor_c(m,3)
        double precision wca(n,m)
        double precision cutoff,w,rij
        double precision x2,y2,z2,dx2,dy2,dz2
        double precision x1,y1,z1

cf2py intent(in) :: coor_r,coor_c,r0,rN,c0,cN,w,wca,n,m
cf2py intent(in,out):: wca
cf2py intent(hide):: i,j,cutoff,rij
cf2py intent(hide):: x2,y2,z2,dx2,dy2,dz2,rij
cf2py intent(hide):: x1,y1,z1

        cutoff = 2.**(1./6.)*w

        do 200,i=r0,rN
            x1=coor_r(i,1)
            y1=coor_r(i,2)
            z1=coor_r(i,3)
            do 100,j=c0,cN
               x2=coor_c(j,1)
               y2=coor_c(j,2)
               z2=coor_c(j,3)

               dx2=(x1-x2)*(x1-x2)
               dy2=(y1-y2)*(y1-y2)
               dz2=(z1-z2)*(z1-z2)
               rij=sqrt(dx2+dy2+dz2)

               if(rij . LT . cutoff) then
                  wca(i,j)=(w/rij)**12.-(w/rij)**6.+0.25
               else
                  wca(i,j)=0
               end if
  100   continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine wca_d(coor,tbead,w,wca,n)
        integer n,tbead,i,j
        double precision coor(n,3)
        double precision wca(n,n)
        double precision cutoff,w,rij
        double precision dx2,dy2,dz2,dx,dy,dz
        double precision x1,y1,z1,x2,y2,z2

cf2py intent(in) :: coor,tbead,w,wca,n
cf2py intent(in,out):: wca
cf2py intent(hide):: i,j,cutoff,rij
cf2py intent(hide):: x2,y2,z2,dx2,dy2,dz2,rij
cf2py intent(hide):: x1,y1,z1

        cutoff = 2.**(1./6.)*w

        do 200,i=tbead,n
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,j=1,i-1
               x2=coor(j,1)
               y2=coor(j,2)
               z2=coor(j,3)

               dx = x2-x1
               dy = y2-y1
               dz = z2-z1

               dx2=dx*dx
               dy2=dy*dy
               dz2=dz*dz
               rij=sqrt(dx2+dy2+dz2)

               if(rij . LT . cutoff) then
                  wca(i,j)=(w/rij)**12.-(w/rij)**6.+0.25
               else
                  wca(i,j)=0
               end if
  100   continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine overlap1(coor,natoms,cutoff,check)
        double precision coor(natoms,3)
        double precision cutoff
        integer natoms,check,count,i,j
        double precision x1,y1,z1,x2,y2,z2,diff2,dist

cf2py intent(in) :: coor,cutoff
cf2py intent(out):: check
cf2py intent(hide)::natoms,i,j
cf2py intent(hide)::x1,y1,z1,x2,y2,z2,diff2,dist

        count = 1
        check = 0
        do 200,i=1,natoms-1
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,j=i+1,natoms
                x2=coor(j,1)
                y2=coor(j,2)
                z2=coor(j,3)
                diff2=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
                dist=sqrt(diff2)
                if(dist . LT . cutoff) then
                    check=1
                    exit
                endif
            count = count + 1

  100       continue

            if(check==1) then
                exit
            endif
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine overlap_skip_n(coor,natoms,nskip,cutoff,check,i,j)
        double precision coor(natoms,3)
        double precision cutoff
        integer natoms,check,count,i,j,nskip
        double precision x1,y1,z1,x2,y2,z2,diff2,dist

cf2py intent(in) :: coor,cutoff,nskip
cf2py intent(out):: check,i,j
cf2py intent(hide)::natoms
cf2py intent(hide)::x1,y1,z1,x2,y2,z2,diff2,dist

        count = 1
        check = 0
        do 200,i=1,natoms-nskip-1
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,j=i+nskip+1,natoms
                x2=coor(j,1)
                y2=coor(j,2)
                z2=coor(j,3)
                diff2=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
                dist=sqrt(diff2)
                if(dist . LT . cutoff) then
C                    write (*,*) "collision between beads"
C                    write (*,*) dist, i, j
C                    write (*,*) "b:", j
C                    write (*,*) "dist",dist
                    check=1
                    exit
                endif
            count = count + 1

  100       continue

            if(check==1) then
                exit
            endif
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012
        subroutine overlap2(coor1a,coor1b,n1a,n1b,cutoff,check,i,j)
        double precision coor1a(n1a,3),coor1b(n1b,3)
        double precision cutoff
        integer n1a,n1b,check,count,i,j
        double precision x1,y1,z1,x2,y2,z2,diff2,dist

cf2py intent(in) :: coor1a,coor1b,cutoff
cf2py intent(out):: check,i,j
cf2py intent(hide)::n1a,n1b
cf2py intent(hide)::x1,y1,z1,x2,y2,z2,diff2,dist

        count = 1
        check = 0
        do 200,i=1,n1a
            x1=coor1a(i,1)
            y1=coor1a(i,2)
            z1=coor1a(i,3)
            do 100,j=1,n1b
                x2=coor1b(j,1)
                y2=coor1b(j,2)
                z2=coor1b(j,3)
                diff2=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
                dist=sqrt(diff2)
                if(dist . LT . cutoff) then
C                    write (*,*) "collision between beads"
C                    write (*,*) dist, i, j
C                    write (*,*) "b:", j
C                    write (*,*) "dist",dist
                    check=1
                    exit
                endif
            count = count + 1

  100       continue

            if(check==1) then
                exit
            endif
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine overlap_dist(coor1a,coor1b,dist,n1a,n1b)
        double precision coor1a(n1a,3),coor1b(n1b,3),dist(n1a*n1b)
        integer n1a,n1b,i,j
        double precision x1,y1,z1,x2,y2,z2,diff2

cf2py intent(in) :: coor1a,coor1b,dist
cf2py intent(in,out):: dist
cf2py intent(hide)::n1a,n1b,i,j
cf2py intent(hide)::x1,y1,z1,x2,y2,z2,diff2

        do 200,i=1,n1a
            x1=coor1a(i,1)
            y1=coor1a(i,2)
            z1=coor1a(i,3)
            do 100,j=1,n1b
                x2=coor1b(j,1)
                y2=coor1b(j,2)
                z2=coor1b(j,3)
                diff2=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
                dist((i-1)*n1a+j)=sqrt(diff2)
  100       continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012
C
        subroutine coulomb(coor,charge,e,t,switchd,nbcutoff,n,el,el2)

        integer n,i,j,np
        double precision coor(n,3)
        double precision charge(n)
        double precision e,t,switchd,nbcutoff
        double precision el,el2

        double precision pi,qconv,qconv2,eps,ang2m,jtokjpmol
        double precision kjtokcal,conv,conv2
        double precision xi,yi,zi,xj,yj,zj
        double precision dx,dy,dz,dx2,dy2,dz2
        double precision rij,switchscale,vij
        double precision arg1,arg2,arg3,qi,qj,kb

cf2py intent(in) :: coor,charge,e,t,n,switchd,nbcutoff
cf2py intent(out):: el,el2
cf2py intent(hide):: i,j,np,pi,qconv,qconv2,eps,ang2m,jtokjpmol
cf2py intent(hide):: kjtokcal,conv,conv2
cf2py intent(hide):: xi,yi,zi,xj,yj,zj
cf2py intent(hide):: dx,dy,dz,dx2,dy2,dz2
cf2py intent(hide):: rij,switchscale,vij
cf2py intent(hide):: arg1,arg2,arg3,qi,qj,kb

C eV --> C                        (yes)
        qconv=1.602177E-19
C q*q --> C^2                     (yes)
        qconv2=qconv*qconv
C epsilon --> C^2/(J*m)           (yes)
        eps=8.854187817E-12
C angstrom --> m                  (yes)
        ang2m=1E-10
C J --> KJ/mol                    (yes)
        jtokjpmol=6.02214129E+20
C KJ --> kcal                     (yes)
        kjtokcal=1d0/4.184

C kB --> (J/K):
        kb=1.3806488E-23

        pi=dacos(-1d0)
C        conv=(qconv2/(4.0*pi*eps*ang2m*e))*jtokjpmol*kjtokcal
C        conv=332.4683
         conv=332.06/e

C       conversion to unitless energy (t: temperature in kelvin)
        conv2=qconv2/(4.0*pi*eps*ang2m*e*kb*t)

        np=0

C  bppbi*bppbj/rij

        el = 0d0
        el2 = 0d0
        switchscale = 1d0

        do 200,i=1,n-3
            qi=charge(i)
            xi=coor(i,1)
            yi=coor(i,2)
            zi=coor(i,3)
            do 100,j=i+4,n
               qj=charge(j)
               xj=coor(j,1)
               yj=coor(j,2)
               zj=coor(j,3)

               dx = xj-xi
               dy = yj-yi
               dz = zj-zi

               dx2=dx*dx
               dy2=dy*dy
               dz2=dz*dz
               rij=sqrt(dx2+dy2+dz2)

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

               vij=((qi*qj)/(rij))*switchscale

               el=el+vij

C               write(*,*) 'el = ',el

               np=np+1

  100   continue
  200   continue

        el2 = el*conv2
        el = el*conv

        end

C123456789012345678901234567890123456789012345678901234567890123456789012
C
        subroutine s_debye(coor,charge,e,t,ld,switchd,nbcutoff,n,energy)

        integer n,i,j
        double precision coor(n,3)
        double precision charge(n)
        double precision e,t,ld,switchd,nbcutoff,energy

C        double precision pi,q,q2,e0,ang2m,jtokjpmol,kb
C        double precision kjtokcal,conv,conv2
        double precision xi,yi,zi,xj,yj,zj
        double precision dx,dy,dz,dx2,dy2,dz2
        double precision rij,switchscale,uij,unitless
        double precision arg1,arg2,arg3,qi,qj

C cf2py intent(in) :: coor,charge,ld,e,t,n,switchd,nbcutoff
C cf2py intent(out):: el,el2
cf2py intent(in) :: coor,charge,e,t,ld,switchd,nbcutoff,n
cf2py intent(out):: energy
cf2py intent(hide):: i,j
cf2py intent(hide):: xi,yi,zi,xj,yj,zj
cf2py intent(hide):: dx,dy,dz,dx2,dy2,dz2
cf2py intent(hide):: rij,switchscale,uij,unitless
cf2py intent(hide):: arg1,arg2,arg3,qi,qj
C cf2py intent(hide):: pi,q,q2,e0,ang2m,jtokjpmol,kb
C cf2py intent(hide):: kjtokcal,conv,conv2

C coor: bead/atom coordinates (Angstroms)
C charge: bead/atom charge (Coulombs)
C e: dielectric constant
C t: temperature (Kelvin)
C ld: Debye length (Angstroms)
C switchd: distance to start transitioning to zero
C nbcutoff: distance to cutoff energy to zero
C n: number of coordinates

C constant --> units
C elementary charge --> C         (yes)
C        q=1.6021766208E-19
C q*q --> C^2                     (yes)
C        q2=q*q
C electric constant --> F/m       (yes)
C        e0=8.854187817E-12
C angstrom --> m                  (yes)
C        ang2m=1E-10
C J --> KJ/mol                    (yes)
C        jtokjpmol=6.02214129E+20
C KJ --> kcal                     (yes)
C        kjtokcal=1d0/4.184
C
C kB --> J/K:
C        kb=1.38064852E-23
C
C        pi=dacos(-1d0)
C
C        conv=(q2/(4.0*pi*e0*ang2m*e))*jtokjpmol*kjtokcal
C        conv=332.4683
C        kcalpmol=332.06/e
C
C       conversion to unitless energy (t: temperature in kelvin)
C       unitless=q2/(4.0*pi*e0*ang2m*e*kb*t) C = 167101.00216036927

        unitless=167101.002160369273042306303977966309/(e*t)
C        write(*,*) 'unitless= ',unitless

C        ekcalpmol = 0d0
        energy = 0d0
        switchscale = 1d0

        do 200,i=1,n+3
            qi=charge(i)
            xi=coor(i,1)
            yi=coor(i,2)
            zi=coor(i,3)
            do 100,j=i+4,n
               qj=charge(j)
               xj=coor(j,1)
               yj=coor(j,2)
               zj=coor(j,3)

               dx = xj-xi
               dy = yj-yi
               dz = zj-zi

               dx2=dx*dx
               dy2=dy*dy
               dz2=dz*dz
               rij=sqrt(dx2+dy2+dz2)

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

               uij=((qi*qj)/(rij))*exp(-rij/ld)*switchscale

               energy=energy+uij

C               write(*,*) 'energy = ',energy

  100   continue
  200   continue

        energy = energy*unitless

C        ekcalpmol = el*kcalpmol

        end
