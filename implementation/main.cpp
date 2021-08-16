#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <iostream>
#include "conformation.h"
#include <list>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <algorithm>
#include <execution>
#include <chrono>
#include <fstream>
#include <sstream>
#include <omp.h>

#define RMSD_THRESHHOLD 0.3
#define ENERGYDIFF_THRESHHOLD 0.0001
#define ENERGY_THRESHHOLD 30 // kcal/mol
#define HtoKCAL 627.509474

int main() {
    BOOST_LOG_TRIVIAL(info) << "Starting the script";
    const std::string target_path( "../../filter_set" );
    const boost::regex my_filter( ".*noh\.xyz" );
    boost::filesystem::directory_iterator end_itr;
    boost::smatch what;
    unsigned int num_files = 0;
    for( boost::filesystem::directory_iterator i(target_path); i != end_itr; ++i ) {
        BOOST_LOG_TRIVIAL(info) << "Checking out = " << i->path().string();
        if( !boost::regex_match( i->path().filename().string(), what, my_filter ) ) continue;
        num_files++;
    }

    auto * testfiles = new boost::filesystem::path [num_files];
    auto * conformations = new Conformation[num_files];
    int count = 0;
    for(boost::filesystem::directory_iterator i(target_path); i != end_itr; ++i) {
        BOOST_LOG_TRIVIAL(info) << "Loading = " << i->path().string();
        if( !boost::regex_match( i->path().filename().string(), what, my_filter ) ) continue;
        testfiles[count] = i->path();
        count++;
    }
    std::sort(std::execution::par_unseq, testfiles, testfiles+num_files);

    #pragma omp parallel for
    for(unsigned int i = 0; i < num_files; ++i){
        BOOST_LOG_TRIVIAL(info) << "Thread " <<  omp_get_thread_num() <<" Added file (sorted) = " << testfiles[i].string();
        conformations[i] = Conformation(testfiles[i].string());
        conformations[i].prepare();
        BOOST_LOG_TRIVIAL(info) << "Thread " <<  omp_get_thread_num() <<" Confname = " << conformations[i].myname;
    }
    delete[] testfiles;

    double minener = conformations[0].energy;
    for(unsigned int i = 0; i < num_files; ++i){
        if(conformations[i].energy < minener) {
            minener = conformations[i].energy;
        }
    }

    for(unsigned int i = 0; i < num_files; ++i){
        if(abs(conformations[i].energy - minener) * HtoKCAL > ENERGY_THRESHHOLD){
            conformations[i].is_duplicate = true;
        }
    }

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    #pragma omp parallel for collapse(2)
    for(unsigned int i = 0; i < num_files; ++i) {
        for (unsigned int j = 0; j < num_files; ++j) {
            if((conformations[i].is_duplicate) || (i == j))
                continue;
            if(/*(!conformations[j].is_duplicate) && */(abs(conformations[i].energy - conformations[j].energy) < ENERGYDIFF_THRESHHOLD) && (conformations[i].rmsd(conformations[j]) < RMSD_THRESHHOLD))
            {
                #pragma omp critical
                    conformations[i].is_duplicate = true;
                BOOST_LOG_TRIVIAL(info) << "Thread " << omp_get_thread_num() << " Erasing " << conformations[i].myname;
            }
        }
    }

    std::ofstream myfile;
    myfile.open ("unique_conformers.txt");
    for (unsigned int i = 0; i < num_files; ++i) {
        if(!conformations[i].is_duplicate)
            myfile << conformations[i].myname << "\n";
    }
    myfile.close();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    BOOST_LOG_TRIVIAL(info) << "Time Elapsed: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " ms";
    delete[] conformations;
    return 0;
}

#pragma clang diagnostic pop