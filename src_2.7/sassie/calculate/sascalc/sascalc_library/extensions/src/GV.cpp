//
// Hailiang Zhang
// NIST & UTK
//

#include "GV.h"
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <string.h>



/*
//////////////////////////////////////////////////////////
/// constructor 
/// currently not implemented due to overhead of B array preparation
//////////////////////////////////////////////////////////
sascalc::GV::
GV(const SasMol & mol, const ScVars & scvar)
{
}
*/

//////////////////////////////////////////////////////////
/// constructor 
/// patched to take the python-preprocess B array as the input
/// assume that the B-array memory layout follows what has been prepared in sascalc_util.py
//////////////////////////////////////////////////////////
sascalc::GV::
GV(const int Natoms, const ScVars & scvar, const double *const B_neutron_array, const double *const B_xray_array):
SasCalc(Natoms, scvar, B_neutron_array, B_xray_array),
_gv_method(scvar._gv_method),
_gv_parameter(scvar._gv_parameter)
{
}

//////////////////////////////////////////////////////////
/// Set the golden vectors
//////////////////////////////////////////////////////////
double *
sascalc::GV::
_getGV(int Ngv) const
{
    // assign Ngv
    Ngv = Ngv;

    // setup GV
    if (Ngv%2==0)
    {
        std::cout<<"The number of golden vectors should be an odd integer, and so it will be reset to be: "<<Ngv+1<<std::endl;
        ++Ngv;
    }
    double *gv = new double[3*Ngv];
    
    const double phi_inv = 2.0/(1+sqrt(5.)); // golden ratio
    double cos_theta, sin_theta, phi;
    double qx,qy,qz;
    int igv;
    const int rank = Ngv/2;
    for (int i=-rank; i<=rank; i++)
    {   
        sin_theta = cos(asin(2.0*i/Ngv));
        cos_theta = 2.0*i/Ngv;
        phi = 2*M_PI*i*phi_inv;
        igv = i + rank;
        gv[igv] = sin_theta*cos(phi);
        gv[Ngv+igv] = sin_theta*sin(phi);
        gv[2*Ngv+igv] = cos_theta;
    }

    // return
    return gv;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// get Ngv_start based on the tolerance value
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int
sascalc::GV::
_tolerance_to_ngv(const double tolerance) const
{
    if (tolerance>=0.05) return 11;
    if (tolerance<0.05 && tolerance>=0.04) return 13;
    if (tolerance<0.04 && tolerance>=0.03) return 15;
    if (tolerance<0.03 && tolerance>=0.02) return 19;
    if (tolerance<0.02 && tolerance>=0.01) return 31;
    if (tolerance<0.01) return 35;
}

//////////////////////////////////////////////////////////
/// Golden vector internal calculator
//  for single item and single frame using fixed method
//////////////////////////////////////////////////////////
void
sascalc::GV::
_calculate_singleFrame_fixed(const int *const flag_skip, const int Nitems, const double *const coor, double *const Iq, const int Ngv)
{
    // locals
    int iatom,iq,igv;
    double x, y, z, b;
    double qx,qy,qz;
    double q_dot_r;
    double sine, cose;
    double *Ireal = new double[Nitems];
    double *Iimag = new double[Nitems];
    double qmag;
    int item;

    double *gv = _getGV(Ngv);

    // summation
    for (iq=0; iq<_Nq; ++iq)
    {
        for (item=0; item<Nitems; ++item) if (!flag_skip[item]) Iq[item*_Nq + iq] = 0.0;
        qmag = _q[iq];
        //std::cout<<"Q: "<<qmag<<std::endl;
        for (igv=0; igv<Ngv; ++igv)
        {
            qx = qmag * gv[igv];
            qy = qmag * gv[Ngv+igv];
            qz = qmag * gv[2*Ngv+igv];
            for (item=0; item<Nitems; ++item)
            {
                Ireal[item] = 0.0;
                Iimag[item] = 0.0;
            }
            for (iatom=0; iatom<_Natoms; ++iatom)
            {
                // get coordinate
                x = coor[iatom];
                y = coor[_Natoms + iatom];
                z = coor[_Natoms*2 + iatom];

                item = 0;
                // neutron
                for (auto &B_neutron : _B_neutron_vector)
                {
                    if (!flag_skip[item]) 
                    {
                        auto B_descriptor = B_neutron.first;
                        auto B = B_neutron.second;
                        b = B[iatom];
                        q_dot_r = qx*x+qy*y+qz*z;
                        sine = sin(q_dot_r);
                        cose = cos(q_dot_r);
                        Ireal[item] += b*cose;
                        Iimag[item] += b*sine;
                    }
                    //if (item==0 && iq==1 && igv==1 && iatom==1) std::cout<<"("<<qx<<" "<<qy<<" "<<qz<<"), ("<<x<<" "<<y<<" "<<z<<") "<<b<<std::endl;
                    ++item;
                }
                // xray
                for (auto &B_xray : _B_xray_vector)
                {
                    if (!flag_skip[item]) 
                    {
                        auto B_descriptor = B_xray.first;
                        auto B = B_xray.second;
                        b = B[_Natoms*iq+iatom];
                        q_dot_r = qx*x+qy*y+qz*z;
                        sine = sin(q_dot_r);
                        cose = cos(q_dot_r);
                        Ireal[item] += b*cose;
                        Iimag[item] += b*sine;
                    }
                    ++item;
                }
            }
            for (item=0; item<Nitems; ++item)  if (!flag_skip[item]) Iq[item*_Nq + iq] += Ireal[item]*Ireal[item] + Iimag[item]*Iimag[item];
        }
        for (item=0; item<Nitems; ++item)  if (!flag_skip[item]) Iq[item*_Nq + iq] /= Ngv;
        //for (item=0; item<Nitems; ++item)  if (item==1) printf("%8.3f %16.8f, ",qmag,Iq[item*_Nq + iq]); std::cout<<std::endl;
    }

    // clean up
    delete [] gv;
    delete [] Ireal;
    delete [] Iimag;
}

//////////////////////////////////////////////////////////
/// Golden vector internal calculator
//  for single item and single frame using converged method
//////////////////////////////////////////////////////////
void
sascalc::GV::
_calculate_singleFrame_converge(const double *const coor, double *const Iq, const double tolerance)
{
    const int Nitems = _B_neutron_vector.size()+_B_xray_vector.size();

    const int Ngv_start = _tolerance_to_ngv(tolerance);

    double * Iq_run_ave = new double[Nitems*_Nq];
    double * Iq_run_ave_previous = new double[Nitems*_Nq];

    int i, item;
    int *flag_skip = new int[Nitems];
    int *Ngvs = new int[Nitems];
    double *err = new double[Nitems];
    for (item=0; item<Nitems; ++item)
    {
        Ngvs[item] = Ngv_start;
        flag_skip[item] = false;
    }
    int offset;
    for (int count=1, Ngv=Ngv_start; count<=100; ++count, Ngv+=2)
    {
        _calculate_singleFrame_fixed(flag_skip, Nitems, coor, Iq, Ngv);
        for (item=0; item<Nitems; ++item)
        {
            if (flag_skip[item]) continue;
            offset = item*_Nq;
            if (count==1) for (i=0; i<_Nq; ++i) Iq_run_ave[offset + i] = Iq[offset + i];
            else
            {
                err[item] = 0.0;
                for (i=0; i<_Nq; ++i)
                {
                    Iq_run_ave_previous[offset+i]=Iq_run_ave[offset+i];
                    Iq_run_ave[offset+i] = (Iq_run_ave[offset+i]*(count-1)+Iq[offset+i])/count;
                    err[item] += fabs(Iq_run_ave[offset+i]-Iq_run_ave_previous[offset+i])/Iq_run_ave_previous[offset+i];
                }
                err[item] /= _Nq;
                if (err[item]<=tolerance) flag_skip[item] = true;
            }
            Ngvs[item] += 2;
        }
    }
    for (item=0; item<Nitems; ++item) Ngvs[item] -= 2;

    // clean up
    delete [] Iq_run_ave;
    delete [] Iq_run_ave_previous;
    delete [] flag_skip;
    delete [] Ngvs;
    delete [] err;
}

int replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

//////////////////////////////////////////////////////////
/// Golden vector calculator
//////////////////////////////////////////////////////////
sascalc::ScResults *
sascalc::GV::
calculate(const double *const coor, const int Nframes)
{
    // initialize a SasResults object
    sascalc::ScResults * pResults = new sascalc::ScResults(Nframes, _Nq, _Qmax);

    const int Nitems = _B_neutron_vector.size()+_B_xray_vector.size();
    double * Iq = new double[Nframes * Nitems * _Nq];
    int item, count;
    int *flag_skip = new int[Nitems];
    for (item=0; item<Nitems; ++item) flag_skip[item] = false;

    // loop over frames
    int offset = 0;
    static int frame_tot = 0;
    for (int frame=0; frame<Nframes; ++frame)
    {
        offset = frame*Nitems*_Nq;
        if (_gv_method==std::string("fixed")) _calculate_singleFrame_fixed(flag_skip, Nitems, coor + frame*_Natoms*3, Iq+offset, int(_gv_parameter));
        else if (_gv_method==std::string("converge")) _calculate_singleFrame_converge(coor + frame*_Natoms*3, Iq+offset, _gv_parameter);

        // save results
        item = 0;
        std::stringstream ss;
        ss<<"frame-"<<frame_tot++<<" ";
        count = 0;
        for (auto &B_neutron : _B_neutron_vector)
        {
            std::string B_descriptor = std::string("neutron ") + ss.str() + B_neutron.first;
            pResults->push_result(std::make_pair(std::string(B_descriptor), Iq + offset + item*_Nq));

            int flag_complete = false;
            char *Bs = strdup(B_descriptor.c_str());
            char * pch;
            pch = strtok (Bs, " ");
            int count1 = 0;
            while (pch != NULL)
            {
                if (count1==3 && std::string(pch)==std::string("complete")) flag_complete=true;
                ++count1;
                pch = strtok (NULL, " ");
            }
            if (flag_complete)
            {
                const double I0 = _neutron_I0[count];
                const double * Iq_data = Iq + offset + item*_Nq;
                const double I0data = *(Iq_data);
                double * Inorm  = new double[_Nq];
                double * Ierr  = new double[_Nq];
                for (int k=0; k<_Nq; ++k){
                    Inorm[k] = Iq_data[k]*(I0/I0data);
                    Ierr[k] = Inorm[k]*0.1;
                }
                std::string tmp = B_neutron.first;
                replace(tmp, "complete", "normalized");
                B_descriptor = std::string("neutron ") + ss.str() + tmp;
                pResults->push_result(std::make_pair(std::string(B_descriptor), Inorm));
                replace(tmp, "normalized", "error");
                B_descriptor = std::string("neutron ") + ss.str() + tmp;
                pResults->push_result(std::make_pair(std::string(B_descriptor), Ierr));
            }

            ++item;
            ++count;
        }
        count = 0;
        for (auto &B_xray : _B_xray_vector)
        {
            std::string B_descriptor = std::string("x-ray ") + ss.str() + B_xray.first;
            pResults->push_result(std::make_pair(std::string(B_descriptor), Iq + offset + item*_Nq));

            int flag_complete = false;
            char *Bs = strdup(B_descriptor.c_str());
            char * pch;
            pch = strtok (Bs, " ");
            int count1 = 0;
            while (pch != NULL)
            {
                if (count1==3 && std::string(pch)==std::string("complete")) flag_complete=true;
                ++count1;
                pch = strtok (NULL, " ");
            }
            if (flag_complete)
            {
                const double I0 = _xray_I0[count];
                const double * Iq_data = Iq + offset + item*_Nq;
                const double I0data = *(Iq_data);
                double * Inorm  = new double[_Nq];
                double * Ierr  = new double[_Nq];
                for (int k=0; k<_Nq; ++k){
                    Inorm[k] = Iq_data[k]*(I0/I0data);
                    Ierr[k] = Inorm[k]*0.1;
                }
                std::string tmp = B_xray.first;
                replace(tmp, "complete", "normalized");
                B_descriptor = std::string("x-ray ") + ss.str() + tmp;
                pResults->push_result(std::make_pair(std::string(B_descriptor), Inorm));
                replace(tmp, "normalized", "error");
                B_descriptor = std::string("x-ray ") + ss.str() + tmp;
                pResults->push_result(std::make_pair(std::string(B_descriptor), Ierr));
            }

            ++item;
            ++count;
        }
    }

    delete [] flag_skip;

    // return
    return pResults;
}

//////////////////////////////////////////////////////////
/// GV deallocation
//////////////////////////////////////////////////////////
sascalc::GV::
~GV()
{
}
