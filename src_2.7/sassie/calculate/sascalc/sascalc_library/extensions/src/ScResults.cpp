//
// Hailiang Zhang
// NIST & UTK
//

#include "ScResults.h"
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <H5Cpp.h>
#include <set>

//////////////////////////////////////////////////////////
/// constructor 
//////////////////////////////////////////////////////////
sascalc::ScResults::
ScResults(const int Nframes, const int Nq, const double Qmax):
_Nframes(Nframes),
_Nq(Nq),
_Qmax(Qmax)
{
    _Q = new double[_Nq];
    for (int iq=0; iq<_Nq; ++iq) _Q[iq] = iq*(_Qmax/(_Nq-1));
}

//////////////////////////////////////////////////////////
/// save the results
//////////////////////////////////////////////////////////
const sascalc::ScResults &
sascalc::ScResults::
save(const std::string &group, const std::string &runname) const
{
    // setup directory
    std::string run_dir = runname+"/"; // it's important that directory name is postfixed with a slash
    static int count = 0;
    if (count==0)
    {
        if (!opendir(run_dir.c_str()))
        {
            mkdir(run_dir.c_str(), S_IRWXU);
            mkdir((run_dir+"sascalc/").c_str(), S_IRWXU);
        }
        else 
        {
            system((std::string("exec rm -rf ")+run_dir+std::string("/")).c_str());
            mkdir(run_dir.c_str(), S_IRWXU);
            mkdir((run_dir+"sascalc/").c_str(), S_IRWXU);
        }
    }

    H5::DataSet dataset;
    H5::DataSpace fspace;

    // create the hdf5 file
    const std::string filename(run_dir+"sascalc/results_raw.h5");
    const std::string h5file_name(filename);
    H5::H5File * h5file;
    if (std::ifstream(filename))
        h5file = new H5::H5File(h5file_name, H5F_ACC_RDWR);
    else
        h5file = new H5::H5File(h5file_name, H5F_ACC_TRUNC);
    const std::string group_path  = std::string("/")+group;
    h5file->createGroup(group_path);
    H5::Attribute attribute;
    H5::DataSpace attribute_space = H5::DataSpace(H5S_SCALAR);
    H5::StrType attribute_strtype = H5::StrType(H5::PredType::C_S1, 256);

    // Q
    hsize_t * fdimq = new hsize_t(1);
    fdimq[0] = _Nq;
    fspace = H5::DataSpace(1, fdimq);
    H5::FloatType datatype( H5::PredType::NATIVE_DOUBLE );
    if (count==0)
    {
        dataset = H5::DataSet(h5file->createDataSet("Q", datatype, fspace));
        dataset.write(_Q, datatype, fspace);
        attribute = dataset.createAttribute(H5std_string("Help"), attribute_strtype, attribute_space);
        attribute.write(attribute_strtype, H5std_string("Q values"));
    }

    // Q
    hsize_t * fdim = new hsize_t(2);
    fdim[0] = _Nq;

    std::set<std::string> sset;
    std::set<std::basic_string<char> >::iterator it;
    for (auto &Iq_data : _Iq_vector)
    {
        auto B_descriptor = Iq_data.first;
        char *Bs = strdup(B_descriptor.c_str());
        char *pch = strtok (Bs," ");
        std::string path;
        std::string full_path;
        do {
            path += std::string("/")+pch;
            pch = strtok (NULL, " ");
            full_path = group_path+path;
            it = sset.find(full_path);
            if (it == sset.end())
            {
                sset.insert(full_path);
                h5file->createGroup(full_path);
            }
        } while (pch != NULL);

        auto Iq = Iq_data.second;
        fspace = H5::DataSpace(1, fdim);
        H5::FloatType datatype( H5::PredType::NATIVE_DOUBLE );
        dataset = H5::DataSet(h5file->createDataSet(group_path+path+"/data", datatype, fspace));
        dataset.write(Iq, datatype, fspace);
        attribute = dataset.createAttribute(H5std_string("Help"), attribute_strtype, attribute_space);
        attribute.write(attribute_strtype, H5std_string("I(q) values"));
    }

    delete h5file;
    delete fdim;
    delete fdimq;

    ++count;

    return *this;
}

/*
//////////////////////////////////////////////////////////
/// epilogure
/// re-organize and save the results for the end-user
//////////////////////////////////////////////////////////
const sascalc::ScResults &
sascalc::ScResults::
epilogue(const std::string &group) const
{
    // setup directory
    std::string run_dir = "run_0/"; // it's important that directory name is postfixed with a slash

    H5::DataSet dataset;
    H5::DataSpace fspace;

    // create the hdf5 file
    const std::string filename(run_dir+"sascalc/results.h5");
    const std::string h5file_name(filename);
    H5::H5File * h5file;
    if (std::ifstream(filename))
        h5file = new H5::H5File(h5file_name, H5F_ACC_RDWR);
    else
        h5file = new H5::H5File(h5file_name, H5F_ACC_TRUNC);
    const std::string group_path  = std::string("/")+group;
    h5file->createGroup(group_path);
    H5::Attribute attribute;
    H5::DataSpace attribute_space = H5::DataSpace(H5S_SCALAR);
    H5::StrType attribute_strtype = H5::StrType(H5::PredType::C_S1, 256);

    // Q
    hsize_t * fdim = new hsize_t(2);
    fdim[0] = _Nq;

    std::map<std::string, std::vector<double>> smap;
    for (auto &Iq_data : _Iq_vector)
    {
        auto B_descriptor = Iq_data.first;
        char *Bs = strdup(B_descriptor.c_str());
        char *pch = strtok (Bs," ");
        std::string path;
        std::string full_path;
        do {
            path += std::string("/")+pch;
            pch = strtok (NULL, " ");
            full_path = group_path+path;
            if (smap.find(full_path) == smap.end())
            {
                sset.insert(full_path);
                h5file->createGroup(full_path);
            }
        } while (pch != NULL);

        auto Iq = Iq_data.second;
        fspace = H5::DataSpace(1, fdim);
        H5::FloatType datatype( H5::PredType::NATIVE_DOUBLE );
        dataset = H5::DataSet(h5file->createDataSet(group_path+path+"/data", datatype, fspace));
        dataset.write(Iq, datatype, fspace);
        attribute = dataset.createAttribute(H5std_string("Help"), attribute_strtype, attribute_space);
        attribute.write(attribute_strtype, H5std_string("I(q) values"));
    }

    delete h5file;
    delete fdim;
    delete fdimq;

    return *this;
}
*/



//////////////////////////////////////////////////////////
/// ScResults deallocation
//////////////////////////////////////////////////////////
sascalc::ScResults::
~ScResults()
{
    delete [] _Q;
}

//////////////////////////////////////////////////////////
/// output
//////////////////////////////////////////////////////////
std::ostream &
sascalc::
operator<<(std::ostream &os, const sascalc::ScResults & results)
{
    for (auto &Iq_data : results._Iq_vector)
    {
        auto B_descriptor = Iq_data.first;
        auto Iq = Iq_data.second;
        os<<B_descriptor<<std::endl;
        for (int i=0; i<3; ++i) os<<Iq[i]<<" "; os<<std::endl;
    }

    return os;
}
