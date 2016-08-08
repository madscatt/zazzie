
        subroutine update_gr(x,y,z,boxl,nbins,deltag,natoms,tgr)
        double precision boxl
        double precision x(natoms),y(natoms),z(natoms)
        integer natoms,i,j,nbins,ig
        double precision tgr(nbins)
        double precision rxij,ryij,rzij,rijsq,deltag

cf2py	intent(in)  x,y,z,boxl,nbins,deltag
cf2py	intent(out) tgr
cf2py	intent(hide) natoms,i,j,ig
cf2py   depend(nbins) tgr
cf2py	intent(hide) rxij,ryij,rzij,rijsq,rcutsq

        rcutsq=(boxl/2.0)**2.0

        do 200,i=1,natoms-1
                do 300,j=i+1,natoms
                        rxij=x(j)-x(i)
                        ryij=y(j)-y(i)
                        rzij=z(j)-z(i)
         
                        rxij=rxij-boxl*(ANINT(rxij/boxl))
                        ryij=ryij-boxl*(ANINT(ryij/boxl))
                        rzij=rzij-boxl*(ANINT(rzij/boxl))
                         
                        rijsq=rxij*rxij+ryij*ryij+rzij*rzij

                        if (rijsq .lt. rcutsq) then
                                ig = int(sqrt(rijsq)/deltag)
                                tgr(ig) = tgr(ig) + 2.0
                        endif

  300   continue
  200   continue

        end
  
