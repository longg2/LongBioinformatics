// This is the boiler plate
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <numeric>
//#include <unistd.h> // Only for debugging
//using namespace std;

void usage(){
	std::cerr << "C++ implementation of the DepthStatistics.awk script"<< std::endl ;
	std::cerr << "It should be noted that the AWK script is faster than this..."<< std::endl ;
	std::cerr <<	"DepthStatistics <FILENAME>" << std::endl;
}

void summaryStats(std::vector<int>& depths, std::string& seq, std::string& fileName){ 	// This is what will do the LCA analysis
	// First we need to figure out the threshold
	double stats[5]; // These are where all of the stats will go. We 
	
	// Need to get the sum real quick
	double sum = std::accumulate(depths.begin(), depths.end(), 0.0);
	double mean = sum/depths.size();
	stats[0] = mean;

	// Now for the SD
	double sumofsquares = std::inner_product(depths.begin(), depths.end(), depths.begin(),0.0);
	double sd = std::sqrt((sumofsquares / depths.size()) - std::pow((sum / depths.size()),2));
	stats[1] = sd;

	// Getting the 95% Error
	stats[2] = 1.96 * sd/std::sqrt(depths.size());

	// Getting the CV
	stats[3] = sd/mean;

	// Getting the PCov
	double tmp = depths.size() - std::count(depths.begin(), depths.end(), 0);
	stats[4] = tmp/depths.size();

	// Here I'm cutting the fileName so that it's only the part of interest
	size_t found = fileName.find_last_of("/") ; // In case it's in a folder
	if(found != std::string::npos){
		fileName = fileName.substr(found, std::string::npos);
	}
	found = fileName.find_last_of(".") ;
	fileName = fileName.substr(0, found);
	
	// Printing out the results
	std::cout << fileName << "\t" << seq;
	for(int i = 0; i < 5; i++){
		std::cout << "\t" << stats[i];
	}
	std::cout << "\n";
};

void readDepth(std::ifstream& file, std::string& sampleName){
	std::vector<int> chromDepths; 
	//ifstream inFile(file);
	std::string entry,tmp,seqID,prevSeqID; // Entry is the line, seqID is the first ID and prevSeqID is the last one
	int dep;
	bool first = true;
	while (std::getline(file, entry)) // Iterate through the file
	{
		std::istringstream iss(entry); // read the line
		std::getline(iss, seqID, '\t'); // Getting the sequence

		if(first){
			prevSeqID = seqID;
			first = false;
		}

		if(seqID != prevSeqID){
			summaryStats(chromDepths, prevSeqID, sampleName);
			chromDepths.clear(); // Emptying the vector after we've used it
		}

		// If they're the same ID as before
		std::getline(iss, tmp, '\t'); // Need to skip the position!
		std::getline(iss, tmp, '\t'); // Getting the count
		dep = std::stoi(tmp);
		chromDepths.push_back(dep);
		prevSeqID = seqID;
	}

		// This is the final bit here, accounting for anything we might have missed
		summaryStats(chromDepths, prevSeqID, sampleName);
}

int main(int argc, char* argv[]){ // The main program runs here
	std::string fileName; // The blast file name
	if(argc != 2){ // Test if we have all the variables
	        usage();
	        return 1;
	}
	fileName = argv[1];

	std::ifstream depthFile; // Creating an ifstream object
	depthFile.open(fileName); // Opening the file

	// Now we need to check if the file exists
	if(depthFile){
		readDepth(depthFile, fileName); // Reading in the blast file
	}else{
		std::cerr << "File: " << fileName << "doesn't exist!" << std::endl;
	        usage();
	        return 1;
	}

	return 0;
}
