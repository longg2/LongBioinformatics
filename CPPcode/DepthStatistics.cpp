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
using namespace std;

void usage(){
	cerr << "C++ implementation of the DepthStatistics.awk script"<< endl ;
	cerr << "It should be noted that the AWK script is faster than this..."<< endl ;
	cerr <<	"DepthStatistics <FILENAME>" << endl;
}

void summaryStats(vector<int>& depths, string& seq, string& fileName){ 	// This is what will do the LCA analysis
	// First we need to figure out the threshold
	double stats[5]; // These are where all of the stats will go. We 
	
	// Need to get the sum real quick
	double sum = accumulate(depths.begin(), depths.end(), 0.0);
	double mean = sum/depths.size();
	stats[0] = mean;

	// Now for the SD
	double sumofsquares = inner_product(depths.begin(), depths.end(), depths.begin(),0.0);
	double sd = sqrt((sumofsquares / depths.size()) - pow((sum / depths.size()),2));
	stats[1] = sd;

	// Getting the 95% Error
	stats[2] = 1.96 * sd/sqrt(depths.size());

	// Getting the CV
	stats[3] = sd/mean;

	// Getting the PCov
	double tmp = depths.size() - count(depths.begin(), depths.end(), 0);
	stats[4] = tmp/depths.size();

	// Here I'm cutting the fileName so that it's only the part of interest
	size_t found = fileName.find_last_of("/") ; // In case it's in a folder
	if(found != string::npos){
		fileName = fileName.substr(found, string::npos);
	}
	found = fileName.find_last_of(".") ;
	fileName = fileName.substr(0, found);
	
	// Printing out the results
	cout << fileName << "\t" << seq;
	for(int i = 0; i < 5; i++){
		cout << "\t" << stats[i];
	}
	cout << "\n";
};

void readDepth(ifstream& file, string& sampleName){
	vector<int> chromDepths; 
	//ifstream inFile(file);
	string entry,tmp,seqID,prevSeqID; // Entry is the line, seqID is the first ID and prevSeqID is the last one
	int dep;
	bool first = true;
	while (getline(file, entry)) // Iterate through the file
	{
		istringstream iss(entry); // read the line
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
		dep = stoi(tmp);
		chromDepths.push_back(dep);
		prevSeqID = seqID;
	}

		// This is the final bit here, accounting for anything we might have missed
		summaryStats(chromDepths, prevSeqID, sampleName);
}

int main(int argc, char* argv[]){ // The main program runs here
	string fileName; // The blast file name
	if(argc != 2){ // Test if we have all the variables
	        usage();
	        return 1;
	}
	fileName = argv[1];

	ifstream depthFile; // Creating an ifstream object
	depthFile.open(fileName); // Opening the file

	// Now we need to check if the file exists
	if(depthFile){
		readDepth(depthFile, fileName); // Reading in the blast file
	}else{
		cout << "File: " << fileName << "doesn't exist!" << endl;
	        usage();
	        return 1;
	}

	return 0;
}
