#ifndef MOLECULE_H
#define MOLECULE_H

#define BOOST_LOG_DYN_LINK 1
#include <boost/log/trivial.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

class Conformation {
    public:
        Conformation();
        Conformation(const Conformation&);
        explicit Conformation(const std::string&);
        ~Conformation();
        void prepare();
        double rmsd(const Conformation&);
        gsl_matrix *xyz_m;
        unsigned int n_atoms;

    private:
        gsl_matrix *c_m, *v_m, *rot_m, *temp_xyz;
        double energy;
        gsl_vector *work_v, *s_v;
        std::string myname;
        std::string print_matrix(const gsl_matrix*);
        std::string print_vector(const gsl_vector*);
        static inline double det3x3(const gsl_matrix* m);
};


#endif //MOLECULE_H
