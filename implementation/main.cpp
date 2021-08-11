#include <iostream>
#include "conformation.h"
#include <list>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <chrono>
#include <fstream>
#include <sstream>

int main() {
    BOOST_LOG_TRIVIAL(info) << "Starting the script";
    const std::string target_path( "../../testset" );
    const boost::regex my_filter( ".*\.xyz" );
    std::list<Conformation> conformations;
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

    typedef std::list<double> OneDList;
    typedef std::list<std::list<double>> TwoDList;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    TwoDList result;
    for (auto & confA : conformations) {
        OneDList curline;
        for (auto &confB : conformations) {
            curline.push_back(confA.rmsd(confB));
        }
        result.push_back(curline);
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    BOOST_LOG_TRIVIAL(info) << "Time Elapsed on RMSDs: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " ms";

    std::ofstream myfile;
    myfile.open ("../../results_cpp.dat");
    for (auto & curline : result) {
        for (auto & curitem : curline) {
            myfile << curitem << " ";
        }
        myfile << "\n";
    }
    myfile.close();
    return 0;
}
