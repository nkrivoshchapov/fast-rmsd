#include <iostream>
#include "conformation.h"
#include <list>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <chrono>

#define RMSD_THRESHHOLD 0.3
#define ENERGYDIFF_THRESHHOLD 0.0001
#define ENERGY_THRESHHOLD 30 // kcal/mol
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

    double minener = conformations.begin()->energy;
    for (auto & conformation : conformations) {
        if(conformation.energy < minener) {
            minener = conformation.energy;
        }
        conformation.prepare();
    }

    unsigned int nconf = conformations.size(), idx = 0;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (auto confA = conformations.begin(); confA != conformations.end(); ++confA) {
        BOOST_LOG_TRIVIAL(info) << "Done " << idx << "/" << nconf;
        if(abs(confA->energy - minener) * HtoKCAL > ENERGY_THRESHHOLD) {
            confA = conformations.erase(confA);
            confA--;
            continue;
        }

        for (auto confB = std::next(confA, 1); confB != conformations.end(); ++confB) {
            if((abs(confA->energy - confB->energy) < ENERGYDIFF_THRESHHOLD) && (confA->rmsd(*confB) < RMSD_THRESHHOLD)) {
                BOOST_LOG_TRIVIAL(info) << "Erasing " << confA->myname;
                confA = conformations.erase(confA);
                confA--;
                break;
            }
        }
        idx++;
    }

    std::ofstream myfile;
    myfile.open ("unique_conformers.txt");
    for (auto & conformation : conformations) {
        myfile << conformation.myname << "\n";
    }
    myfile.close();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    BOOST_LOG_TRIVIAL(info) << "Time Elapsed: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " ms";
    return 0;
}
