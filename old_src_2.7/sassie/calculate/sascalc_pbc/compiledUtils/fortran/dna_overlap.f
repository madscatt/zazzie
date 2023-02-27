C
C Author:  Steven C. Howell
C Purpose: Calculate interaction energy between every pair of atoms
C Created: January 2014
C
C $Id: dna_overlap.f 2545 2015-06-20 04:51:24Z schowell $
C
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine debye(coor,Q,B,I_tmp,qI,box,PBC,n,numQ)
        integer n,i,j,numQ,qI
        double precision coor(n,3)
        double precision Q,Rjk,qr,b_k,b_j,q_dot_r
        double precision B(numQ)
        double precision x2,y2,z2,dx2,dy2,dz2
        double precision x1,y1,z1,box,x12,y12,z12
        double precision I_tmp
        logical PBC

cf2py intent(in) :: coor,Q,B,qI,box,PBC,n,numQ
cf2py intent(in,out):: I_tmp
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

        I_tmp=0
        do 200,j=1,n
            x1=coor(j,1)
            y1=coor(j,2)
            z1=coor(j,3)
            b_j = B(n*qI+j)
            do 100,k=1,n
               x2=coor(k,1)
               y2=coor(k,2)
               z2=coor(k,3)
               


              if(PBC) THEN
                   x12=(x1-x2)
                   y12=(y1-y2)
                   z12=(z1-z2)

                   x12=x12-box*dnint(x12/box) ! periodic boundary
                   y12=y12-box*dnint(y12/box) ! conditions
                   z12=z12-box*dnint(z12/box)
                   Rjk=dsqrt(x12**2+y12**2+z12**2)
               else

                   dx2=(x1-x2)*(x1-x2)
                   dy2=(y1-y2)*(y1-y2)
                   dz2=(z1-z2)*(z1-z2)
                   Rjk=dsqrt(dx2+dy2+dz2)
               endif
               

               qr = Q*Rjk
               if(qr.NE.0) THEN
                   b_k = B(n*qI+k)
                   I_tmp = I_tmp + b_j*b_k*dsin(qr)/qr
               endif
               if(qr.eq.0) THEN
                   b_k=B(n*qI+k)
                   I_tmp = I_tmp+b_j*b_k
                endif
            
  100   continue
  200   continue

        end
c         1         2         3         4         5         6         7
c123456789012345678901234567890123456789012345678901234567890123456789012

c The Radial Distribution Function
c D. Frenkel & B. Smit: Understanding Molecular Simulation 2nd Ed
c Algorithm 7, pg 86

        subroutine gr(switch,coor,g,box,delg,ngr,nhis,npart)
        integer ngr,nhis,npart,switch,i,j
        double precision coor(npart,3),box,g(nhis)
        double precision nid,vb,r,rho,pi
        double precision x1,y1,z1,x2,y2,z2,x12,y12,z12

cf2py intent(in) :: coor,ngr,box,nhis,npart,delg
cf2py intent(in,out):: g
cf2py intent(hide):: i,j,cutoff,rij
cf2py intent(hide):: x2,y2,z2,dx2,dy2,dz2,rij
cf2py intent(hide):: x1,y1,z1

        ngr=ngr+1
        do 200 i=1,npart-1 ! loop over all atom pairs
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100 j=i+1,npart
                x2=coor(j,1)
                y2=coor(j,2)
                z2=coor(j,3)

                x12=x1-x2
                y12=y1-y2
                z12=z1-z2

                x12=x12-box*dnint(x12/box) ! periodic boundary
                y12=y12-box*dnint(y12/box) ! conditions
                z12=z12-box*dnint(z12/box)

                rij=dsqrt(x12**2+y12**2+z12**2)

                if (rij.lt.box/2) then ! only within half the box
                    ig=int(rij/delg)
                    g(ig)=g(ig)+2
                endif
100           continue
200       continue
        return
        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine calc_full_2(dist,gv,Q,B,qI,I_tmp,n,bsize,gv_size)
        integer n,i,j,qI,bsize,gv_size
        double precision dist(n,n,3)
        double precision gv(gv_size,3)
        double precision B(bsize) 
        double precision I_tmp, b_j,b_k,q_dot_r
        double precision cutoff,w,rij
        double precision x2,y2,z2,dx2,dy2,dz2
        double precision x1,y1,z1,Q,q_(3)

cf2py intent(in) :: distances,gold_vects,Q,B,qI,n,bsize,gv_size 
cf2py intent(in,out):: I_tmp
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

        do 300,z=1,gv_size
            q_(1) = gv(z,1)*Q
            q_(2) = gv(z,2)*Q
            q_(3) = gv(z,3)*Q
            do 200,j=1,n
                b_j = B(n*qI+j)
                b_j = 1
                do 100,k=1,n
                    !if(j.EQ.k) THEN
                    !    I_tmp = I_tmp + b_j*b_j
                    !endif
                    q_dot_r=q_(1)*dist(j,k,1)+q_(2)*dist(j,k,2)!&
                    !    + q(3)*distances(j,k,3)
                    q_dot_r = q_dot_r + q_(3)*dist(j,k,3)
                    b_k = B(n*qI+k)
                    b_k = 1
                    I_tmp = I_tmp +b_j*b_k*(dsin(q_dot_r))

  100   continue
  200   continue
  300   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine calc_full(distances,q,B,qI,I_tmp,n,bsize)
        integer n,i,j,qI,bsize
        double precision distances(n,n,3)
        double precision q(3), B(bsize) 
        double precision I_tmp, b_j,b_k,q_dot_r
        double precision cutoff,w,rij
        double precision x2,y2,z2,dx2,dy2,dz2
        double precision x1,y1,z1

cf2py intent(in) :: coor,r0,rN,c0,cN,w,wca,n
cf2py intent(in,out):: I_tmp
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


        do 200,j=1,n
            b_j = B(n*qI+j)
            do 100,k=1,n
                if(j.EQ.k) THEN
                    I_tmp = I_tmp + b_j*b_j
                endif
                q_dot_r = q(1)*distances(j,k,1)+q(2)*distances(j,k,2) !&
!                    + q(3)*distances(j,k,3)
                q_dot_r = q_dot_r + q(3)*distances(j,k,3)
                b_k = B(n*qI+k)
                I_tmp = I_tmp + b_j*b_k*dcos(q_dot_r)

  100   continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine scatteringExp(coor,q,wca,n)
        integer n,i,j,r0,rN,c0,cN
        double precision coor(n,3)
        double precision Q 
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


        do 200,j=1,n
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,k=1,n
               x2=coor(j,1)
               y2=coor(j,2)
               z2=coor(j,3)

               dx2=(x1-x2)
               dy2=(y1-y2)
               dz2=(z1-z2)
              ! q_dot_r = qx *dx2 + qy*dy2+qz*dz2
              ! dexp(q_dot_r)
              ! Rjk=sqrt(dx2+dy2+dz2)

  100   continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine distance_vectors(coor,wca,n)
        integer n,i,j,r0,rN,c0,cN
        double precision coor(n,3)
        double precision wca(n,n,3)
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


        do 200,i=1,n
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,j=1,n
               x2=coor(j,1)
               y2=coor(j,2)
               z2=coor(j,3)
               wca(i,j,1)=(x1-x2) 
               wca(i,j,2)=(y1-y2)
               wca(i,j,3)=(z1-z2)

  100   continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine distance_vectors_pbc(coor,box,wca,n)
        integer n,i,j,r0,rN,c0,cN
        double precision coor(n,3)
        double precision wca(n,n,3)
        double precision cutoff,w,rij
        double precision x2,y2,z2,dx2,dy2,dz2
        double precision x1,y1,z1,x12,y12,z12

cf2py intent(in) :: coor,box,n
cf2py intent(in,out):: wca
cf2py intent(hide):: i,j,cutoff,rij
cf2py intent(hide):: x2,y2,z2,dx2,dy2,dz2,rij
cf2py intent(hide):: x1,y1,z1,x12,y12,z12

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


        do 200,i=1,n
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,j=1,n
               x2=coor(j,1)
               y2=coor(j,2)
               z2=coor(j,3)
               x12 = x1-x2
               y12 = y1-y2
               z12 = z1-z2
               wca(i,j,1)= x12 - x12*dnint(x12/box)
               wca(i,j,2)= y12 - y12*dnint(y12/box)
               wca(i,j,3)= z12 - z12*dnint(z12/box)

  100   continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine g_hist(coor,wca,n)
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

C    For PBC
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


        do 200,i=1,n
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,j=1,n
               x2=coor(j,1)
               y2=coor(j,2)
               z2=coor(j,3)

               dx2=(x1-x2)*(x1-x2)
               dy2=(y1-y2)*(y1-y2)
               dz2=(z1-z2)*(z1-z2)
               wca(i,j)=sqrt(dx2+dy2+dz2)

  100   continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine distances(coor,wca,n)
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


        do 200,i=1,n
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,j=1,n
               x2=coor(j,1)
               y2=coor(j,2)
               z2=coor(j,3)

               dx2=(x1-x2)*(x1-x2)
               dy2=(y1-y2)*(y1-y2)
               dz2=(z1-z2)*(z1-z2)
               wca(i,j)=sqrt(dx2+dy2+dz2)

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
