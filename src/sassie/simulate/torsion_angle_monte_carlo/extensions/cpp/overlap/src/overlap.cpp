#include <overlap.h>

#include <iostream>

int
overlap(const double * const coor, const double * const r, const int natoms, const double scale, const int * const bond_table)
{
        int flag_overlap = 0;
        double xi,yi,zi,ri,xj,yj,zj,rj,dist ;
        int i,j;
        for (i=0; i<natoms-1 ; ++i)
        {
            xi=coor[3*i];
            yi=coor[3*i+1];
            zi=coor[3*i+2];
            ri=r[i];
            for (j=i+1; j<natoms ; ++j)
            {
                if (bond_table && bond_table[(j*(j-1))/2+i]) continue;
                //if (bond_table && bond_table[(j*(j-1))/2+i]) {printf("test bonds: %d %d\n",i,j); continue;}
                xj=coor[3*j];
                yj=coor[3*j+1];
                zj=coor[3*j+2];
                rj=r[j];
                dist=sqrt((xj-xi)*(xj-xi)+(yj-yi)*(yj-yi)+(zj-zi)*(zj-zi)) ;
                //if (i==2 && j==3) std::cout<<"ZHL ri: "<<ri<<" rj: "<<rj<<" dist: "<<dist<<" scale: "<<scale<<std::endl;
                if (dist < ((ri+rj)*scale))
                {
                    //printf("Overlap found between atomic index: %d %d\n",i,j);
                    flag_overlap=1 ;
                    break ;
                }
            }
            if(flag_overlap==1)	break ;
        }

        return flag_overlap;
}
