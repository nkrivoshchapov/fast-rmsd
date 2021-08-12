#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include "conformation.h"
#include <list>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <chrono>
#include <fstream>
#include <sstream>
#include <omp.h>

#define RMSD_THRESHHOLD 0.03
#define ENERGY_THRESHHOLD 0.01 // kcal/mol
#define HtoKCAL 627.509474

int main() {
    BOOST_LOG_TRIVIAL(info) << "Starting the script";
    const std::string target_path( "../../filter_set" );
    const boost::regex my_filter( ".*\.xyz" );
    typedef std::list<Conformation> conflist ;
    conflist conformations;
    std::list<boost::filesystem::path> testfiles;
    boost::filesystem::directory_iterator end_itr;
    boost::smatch what;
    for( boost::filesystem::directory_iterator i(target_path); i != end_itr; ++i ) {
        BOOST_LOG_TRIVIAL(info) << "Checking out = " << i->path().string();
        if( !boost::regex_match( i->path().filename().string(), what, my_filter ) ) continue;
        testfiles.push_back(i->path());
    }

    testfiles.sort();

    
    for (auto & testfile : testfiles) {
        conformations.emplace_back(Conformation(testfile.string()));
        BOOST_LOG_TRIVIAL(info) << "Added file (sorted) = " << testfile.string();
    }

    for (auto & conformation : conformations) {
        conformation.prepare();
    }

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    unsigned int nconf = conformations.size(), idx = 0;
    #pragma omp parallel
    #pragma omp single
    {
        for(auto confA = conformations.begin(); confA != conformations.end(); ++confA) {
            if(confA->is_duplicate)
                continue;

            #pragma omp task firstprivate(confA)
            {
                BOOST_LOG_TRIVIAL(info) << "Thread " << omp_get_thread_num() << " Done " << idx << "/" << nconf;
                for (auto confB = std::next(confA, 1); confB != conformations.end(); ++confB) {
                    if((!confB->is_duplicate) && (abs(confA->energy - confB->energy) * HtoKCAL < ENERGY_THRESHHOLD) && (confA->rmsd(*confB) < RMSD_THRESHHOLD)) {
                        BOOST_LOG_TRIVIAL(info) << "Thread " << omp_get_thread_num() << " Erasing " << confA->myname;
                        confA->is_duplicate = true;
                        break;
                    }
                }
                idx++;
            }
        }
        #pragma omp taskwait
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    BOOST_LOG_TRIVIAL(info) << "Time Elapsed: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " ms";
    return 0;
}

#pragma clang diagnostic pop