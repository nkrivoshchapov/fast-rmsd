#include "conformation.h"
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

Conformation::Conformation() :
    is_duplicate(false) {
    //BOOST_LOG_TRIVIAL(info) << "Running empty constructor";
}

Conformation::Conformation(const std::string& filename) :
    is_duplicate(false) {
    //BOOST_LOG_TRIVIAL(info) << "Initializing conformation from " << filename;
    std::ifstream infile(filename);
    myname = filename;
    std::string line;
    unsigned int linecount = 0;
    double x,y,z;
    std::string sym;
    while (std::getline(infile, line))
    {
        //BOOST_LOG_TRIVIAL(info) << line;
        std::istringstream iss(line);
        if(linecount == 0) {
            if (!(iss >> n_atoms)) {
                //BOOST_LOG_TRIVIAL(error) << "Unable to read input file. " << line;
                return;
            } else {
                xyz_m = gsl_matrix_alloc(n_atoms, 3);
                temp_xyz = gsl_matrix_alloc(n_atoms, 3);
                //BOOST_LOG_TRIVIAL(info) << "Natoms = " << n_atoms;
            }
        }

        if(linecount == 1) {
            if (!(iss >> energy)) {
                //BOOST_LOG_TRIVIAL(error) << "Unable to read input file. " << line;
                return;
            } else {
                //BOOST_LOG_TRIVIAL(info) << "Energy = " << energy;
            }
        }

        if((linecount >= 2) && (linecount <= 2+n_atoms)){
            if (!(iss >> sym >> x >> y >> z)) {
                //BOOST_LOG_TRIVIAL(error) << "Unable to read input file. " << line;
                return;
            } else {
                gsl_matrix_set(xyz_m,linecount-2,0,x);
                gsl_matrix_set(xyz_m,linecount-2,1,y);
                gsl_matrix_set(xyz_m,linecount-2,2,z);
                //BOOST_LOG_TRIVIAL(info) << "x = " << x
//                                        << " y = " << y
//                                        << " z = " << z
//                                        << " sym = " << sym;
            }
        }
        linecount++;
    }
    c_m = gsl_matrix_alloc(3, 3);
    v_m = gsl_matrix_alloc(3, 3);
    rot_m = gsl_matrix_alloc(3, 3);
    work_v = gsl_vector_alloc(3);
    s_v = gsl_vector_alloc(3);
    //BOOST_LOG_TRIVIAL(info) << "Coordinate matrix = \n" << print_matrix(xyz_m);
}

//TODO use move semantics
Conformation::Conformation(const Conformation& other) {
    //BOOST_LOG_TRIVIAL(info) << "Creating a copy of  " << other.myname;
    is_duplicate = other.is_duplicate;
    myname = std::string("Copy of " + other.myname);
    energy = other.energy;
    n_atoms = other.n_atoms;
    xyz_m = gsl_matrix_alloc(n_atoms, 3);
    gsl_matrix_memcpy(xyz_m,other.xyz_m);
    temp_xyz = gsl_matrix_alloc(n_atoms, 3);
    c_m = gsl_matrix_alloc(3, 3);
    v_m = gsl_matrix_alloc(3, 3);
    rot_m = gsl_matrix_alloc(3, 3);
    work_v = gsl_vector_alloc(3);
    s_v = gsl_vector_alloc(3);
}

Conformation::~Conformation() {
    //BOOST_LOG_TRIVIAL(info) << "Running destructor of  " << myname;
    gsl_matrix_free(xyz_m);
    gsl_matrix_free(c_m);
    gsl_matrix_free(v_m);
    gsl_matrix_free(rot_m);
    gsl_matrix_free(temp_xyz);
    gsl_vector_free(s_v);
    gsl_vector_free(work_v);
}

std::string Conformation::print_matrix(const gsl_matrix *m)
{
    std::string res;
    for (size_t i = 0; i < m->size1; i++) {
        for (size_t j = 0; j < m->size2; j++) {
            res +=  boost::lexical_cast<std::string>(gsl_matrix_get(m, i, j));
            res += " ";
        }
        res += "\n";
    }
    return res;
}

std::string Conformation::print_vector(const gsl_vector *v)
{
    std::string res;
    for (size_t i = 0; i < v->size; i++) {
        res +=  boost::lexical_cast<std::string>(gsl_vector_get(v, i));
        res += " ";
    }
    return res;
}

void Conformation::prepare() {
    double aver_x = 0, aver_y = 0, aver_z = 0;
    for (size_t i = 0; i < xyz_m->size1; i++) {
        aver_x += gsl_matrix_get(xyz_m, i, 0);
        aver_y += gsl_matrix_get(xyz_m, i, 1);
        aver_z += gsl_matrix_get(xyz_m, i, 2);
    }
    aver_x /= n_atoms;
    aver_y /= n_atoms;
    aver_z /= n_atoms;
    //BOOST_LOG_TRIVIAL(info) << "Calculated centroid for " << myname <<
//            " = " << aver_x << " " << aver_y << " " << aver_z;

    for (size_t i = 0; i < xyz_m->size1; i++) {
        gsl_matrix_set(xyz_m, i, 0,gsl_matrix_get(xyz_m, i, 0) - aver_x);
        gsl_matrix_set(xyz_m, i, 1,gsl_matrix_get(xyz_m, i, 1) - aver_y);
        gsl_matrix_set(xyz_m, i, 2,gsl_matrix_get(xyz_m, i, 2) - aver_z);
    }
    //BOOST_LOG_TRIVIAL(info) << "Coordinate matrix after subtraction = \n" << print_matrix(xyz_m);
}

inline double Conformation::det3x3(const gsl_matrix* m){
    return  gsl_matrix_get(m, 0, 0)*gsl_matrix_get(m, 1, 1)*gsl_matrix_get(m, 2, 2)+
            gsl_matrix_get(m, 1, 0)*gsl_matrix_get(m, 2, 1)*gsl_matrix_get(m, 0, 2)+
            gsl_matrix_get(m, 2, 0)*gsl_matrix_get(m, 1, 2)*gsl_matrix_get(m, 0, 1)-
            gsl_matrix_get(m, 2, 0)*gsl_matrix_get(m, 1, 1)*gsl_matrix_get(m, 0, 2)-
            gsl_matrix_get(m, 1, 0)*gsl_matrix_get(m, 0, 1)*gsl_matrix_get(m, 2, 2)-
            gsl_matrix_get(m, 2, 1)*gsl_matrix_get(m, 1, 2)*gsl_matrix_get(m, 0, 0);
}

double Conformation::rmsd(const Conformation& other) {
    gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,xyz_m,other.xyz_m,0,c_m);
    //BOOST_LOG_TRIVIAL(info) << "C = A.T * B = \n" << print_matrix(c_m);

    gsl_linalg_SV_decomp(c_m, v_m, s_v, work_v); // c_m -> U
    //BOOST_LOG_TRIVIAL(info) << "SVD done! U = \n" << print_matrix(c_m);
    //BOOST_LOG_TRIVIAL(info) << "V = \n" << print_matrix(v_m);
    //BOOST_LOG_TRIVIAL(info) << "S = \n" << print_vector(s_v);

    double crit = det3x3(c_m)*det3x3(v_m);
    //BOOST_LOG_TRIVIAL(info) << "det(U) = \n" << det3x3(c_m);
    //BOOST_LOG_TRIVIAL(info) << "det(V) = \n" << det3x3(v_m);
    //BOOST_LOG_TRIVIAL(info) << "Criteria = " << crit;
    if(crit < 0){
//        gsl_vector_set(s_v,2,-gsl_vector_get(s_v,2));
//        //BOOST_LOG_TRIVIAL(info) << "Corrected S = \n" << print_vector(s_v);
        for (size_t i = 0; i < 3; i++) {
            gsl_matrix_set(c_m,i,2,-gsl_matrix_get(c_m,i,2));
        }
        //BOOST_LOG_TRIVIAL(info) << "Corrected U = \n" << print_matrix(c_m);
    }
    gsl_blas_dgemm(CblasNoTrans, CblasTrans,1, c_m, v_m,0, rot_m);
    //BOOST_LOG_TRIVIAL(info) << "Rotation matrix = U * V = \n" << print_matrix(rot_m);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1,xyz_m,rot_m,0,temp_xyz);
    //BOOST_LOG_TRIVIAL(info) << "Rotated A = \n" << print_matrix(temp_xyz);

    gsl_matrix_sub(temp_xyz, other.xyz_m);
    //BOOST_LOG_TRIVIAL(info) << "Difference = \n" << print_matrix(temp_xyz);

    double res = 0;
    for (size_t i = 0; i < temp_xyz->size1; i++) {
        for (size_t j = 0; j < temp_xyz->size2; j++) {
            res += pow(gsl_matrix_get(temp_xyz, i, j), 2);
        }
    }
    res /= n_atoms;
    res = sqrt(res);
    //BOOST_LOG_TRIVIAL(info) << "Calced RMSD = \n" << res;
    return res;
}
