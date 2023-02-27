from cffi import FFI
ffibuilder = FFI()

ffibuilder.cdef("""
    double sQ(double**, double *, int, int);
""")
ffibuilder.set_source("_cLoops",
"""
    double sQ(double** coor,double* q,int x, int y){
        int j, k; 
        double dotX,dotY,dotZ;
        double sum = 0.0;
    
        for (j=0; j<x; j++){
            for (k=0; k<y; k++){
                dotX = q[0]*(coor[j][0]-coor[k][0]);
                dotY = q[1]*(coor[j][1]-coor[k][1]);
                dotZ = q[2]*(coor[j][2]-coor[k][2]);
                sum = sum+cos(dotX+dotY+dotZ);
            }
        }
        return(sum);
}
""")
if __name__=="__main__":
    ffibuilder.compile(verbose=True)
