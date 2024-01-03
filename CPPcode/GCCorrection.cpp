// Trying to simplify GC Corrections
//

// This is the boiler plate
#include <iostream>
#include <vector>
#include <map>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <algorithm>
#include <cmath>
#include <numeric>
//#include <unistd.h> // Only for debugging
using namespace std;

void usage(){
	cout << "C++ implementation of the GC correction that I often end up doing"<< endl ;
	cout << "PLAN: If I'm giving it genes it'll just do it one group at a time, otherwise, every xKb"<< endl ;
	cout <<	"GCCorrection <DepthFile> <FastaFile> <Float or Integer>" << endl;
}

//class DNASeq{
//		string sequence, name; // Want to leave the sequence and name unchanged throughout the whole thing
//	public:
//		void insert (string, string); // Inserting the sequence
//		int len () {return (sequence.length());} // Length of the sequence
//
//
//};
//
//float gcContentCalc(string& DNA){
//	// Need to get the sum of the GCs
//	int count = 0;
//
//	// The Gs
//	for(int i = 0; (i = DNA.find('G',i)) != string::npos; i++){
//		count++;
//	}
//	
//	// The Cs
//	for(int i = 0; (i = DNA.find('C',i)) != string::npos; i++){
//		count++;
//	}
//}

class DNA{ // My DNA Structure. Might want to convert to a class to save space?
		unsigned char byteSeq ;  // Converting the sequence to bytes
	public:
		void init (string);
		string name; // Name of the sequence
		string seq; // The actual sequence
		unsigned int len = seq.length(); // Length of the sequence
	//float gcContent = gcContentCalc(seq); // The GC Content
};

void DNA::init(string seq){
	std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper); // Converting string to uppercase
	for(auto &&a : seq){
		switch (a){ // Need to assign the bits
			case 'A':
				0b00 ;
			case 'T':
				0b11 ;
			case 'C':
				0b01 ;
			case 'G':
				0b10 ;

		}
	}
}
//struct DNA_Seq{ // My DNA Structure. Might want to convert to a class to save space?
//	string name; // Name of the sequence
//	string seq; // The actual sequence
//	unsigned int len = seq.length(); // Length of the sequence
//	//float gcContent = gcContentCalc(seq); // The GC Content
//};

struct SlideWindow{ // My DNA Structure. Might want to convert to a class to save space?
	int pos;
	float mean;
	//float gcContent = gcContentCalc(seq); // The GC Content
};

int catchNA(string test){
	const string testVal = "NA";
	
	if(test == testVal){
		return 0;
	}else{
		return stoi(test);
	}
}

void summaryStats(std::vector<int>& depths, std::string& seq){ 	// This is what will do the LCA analysis
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
//	size_t found = fileName.find_last_of("/") ; // In case it's in a folder
//	if(found != std::string::npos){
//		fileName = fileName.substr(found, std::string::npos);
//	}
//	found = fileName.find_last_of(".") ;
//	fileName = fileName.substr(0, found);
//	
//	// Printing out the results
//	std::cout << fileName << "\t" << seq;
//	for(int i = 0; i < 5; i++){
//		std::cout << "\t" << stats[i];
//	}
//	std::cout << "\n";
};

slidingWindow(vector<int>& listNums, float& k){
	int slideSize; 
	int pos;
	float sumDepths = 0;
	float mean;

	if(k < 1){ // First, we need to convert k into an integer if less than zero
		slideSize = static_cast<int>(round(listNums.size() * k));
	}
//TODO: I'M RIGHT HERE
	for(auto i : listNums){
		if(pos % slideSize == 0){ // If we've reached the size I want
			mean = sumDepths / pos; 
			sumDepths = 0;
		}
	}

}

int readFasta(ifstream& file, map<string,DNA_Seq>& fastaSeqs){ // Need to figure out proper understanding of void and getting variables
	//vector<DNA_Seq>  fastaSeqs;
	bool first=true; // Want to test if this is really a fasta file! Will need to read the first line
	string entry, name, sequence;
	DNA_Seq fasta;
	
	while (getline(file, entry)) // Iterating through the file
	{
		istringstream iss(entry); // reading the line
		if(first){ // If on the first line
			first = false;
			if(entry[0] == '>'){ // If it's indeed a fasta file
				name = entry.substr(1, string::npos); // Don't want to save the '>'
			}else{ // if it isn't
				return 1;
				break;
			}
		}else if(entry[0] == '>'){ // If we're pulling out the name, we need to first save what we currently have 
			fasta.name = name;
			fasta.seq = sequence;
			fasta.len = sequence.length();

			fastaSeqs.insert({name, fasta});
			name = entry.substr(1, string::npos);
			sequence = ""; 
		}else{
			sequence = sequence.append(entry); // Appending the sequence that we have

		}
	}

	// Once we're out of the loop, we need to save what we currently have in our loop
	fastaSeqs.insert({name, fasta});
	return 0;
}

void readDepth(std::ifstream& file, map<string, Stats>& summStats, float& windowSize){
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
			slidingWindow(chromDepths, windowSize);
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
		double sum = std::accumulate(depths.begin(), depths.end(), 0.0);
		double mean = sum/depths.size();
		stats[0] = mean;
}

int main(int argc, char* argv[]){ // The main program runs here
	string depthFileName, fastaFileName; // The blast file name
	float windowSize;
	 if(argc != 4){ // Test if we have all the variables
		 usage();
		 return 1;
	 }
	 depthFileName = argv[1];
	 fastaFileName = argv[2];
	 windowSize = stof(argv[3]);

	 // Opening the files
	 ifstream depthFile;
	 depthFile.open(depthFileName);

	 ifstream fastaFile;
	 fastaFile.open(fastaFileName);

	 if(depthFile && fastaFile){ // Maybe?
		cerr << "It read both files!" << endl;		// 
		map<string, DNA_Seq> fastaSeqs;
		map<string, Stats> depthStats;
		readFasta(fastaFile, fastaSeqs);
		readDepth(depthFile, depthStats, windowSize);

		for(auto const& x : fastaSeqs){
			cout << ">" << x.first << endl;
			cout << x.second.seq << endl;
		}
		return 0;
	 }else{ // If both files aren't opened, need to figure out which is the problem child
		 if(depthFile){
		 	cerr << "File: " << fastaFileName << "doesn't exist!" << endl;
		 	usage();
		 	return 1;
		 }else{
		 	cerr << "File: " << depthFileName << "doesn't exist!" << endl;
		 	usage();
		 	return 1;
		 }
	 }

	 return 0;
}
