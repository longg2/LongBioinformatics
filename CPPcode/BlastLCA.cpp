// This is going to take a while to figure out, but, I want to do this in C++
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
//#include <unistd.h> // Only for debugging
using namespace std;

void usage(){
	cout << "C++ implementation of the BlastLCA.R script"<< endl ;
	cout << "This version assumes that the input file is sorted by sequence name"<< endl ;
	cout << "If that is not the case, please use BlastLCAUnsorted instead"<< endl ;
	cout <<	"BlastLCA <FILENAME> <Detection threshold > 0.5 and <= 1>" << endl;
}

int catchNA(string test){
	const string testVal = "NA";
	
	if(test == testVal){
		return 0;
	}else{
		return stoi(test);
	}
}

void printRecord(string seqID, array<int,7> record){ 	// printing out the record
	string tmp;
	tmp = seqID;
	for(int i = 0; i < 7; i++){
		tmp = tmp + "\t" + to_string(record[i]);
	}
	cout << tmp << endl;
};

int threshSearch(vector<int> counts, int threshLCA){
	int i;
	for(i = 0; i < counts.size(); i++){
		//cerr << threshLCA << " " << counts[i] << endl;
		if(counts[i] >= threshLCA){
			return i;
		}
	}
	return -1	;
};

array<int,7> array8to7(array<int,8> begArray){
	array<int,7> endArray;
		for(int i = 1;i < 8;i++){
			endArray[i - 1] = begArray[i];
		}
		return endArray;
}

array<int,7> LCA(vector<array<int,8>>& record, float& th){ 	// This is what will do the LCA analysis
	// First we need to figure out the threshold
	int counts,i,threshLCA;
	int totalHits = 0;
	vector<int>::iterator it;
	vector<int> allCounts;
	
	// First, if there's only one element in the vector, might as well print it
	if(record.size() == 1){	return array8to7(record.at(0));	}

//	cerr << "Multiple records found" << endl;
	for(auto& a : record){ // Quickly getting the total count size
		totalHits += a[0];
//		cerr << totalHits << "\t" << a[0] << endl;
		allCounts.push_back(a[0]);
	}

	threshLCA = ceil(totalHits * th);
	//cerr << "ThreshLCA = " << threshLCA << endl;

	// Looking to see if we can get away at the species level
	i = threshSearch(allCounts, threshLCA)	;
	//cerr << "After the search: i = " << i << ", allCounts Size = "<< allCounts.size() << endl;

	if(i > -1){ // if there's a species level hit
		//cerr << "Found a one off!" << endl;
		if(record.at(i)[7] > 0){ // If the taxon is NA, we're not interested....
			return array8to7(record.at(i));
		}
	};
	
	//cerr << "Now entering the loop" << endl;
	
	for(int t = 6; t != 0; t--){ // traversing the taxonomy backwards
		//cerr << t << endl;
			// Need to figure out the taxa
		map<int, int> taxCounts; // This will have the counts
		map<int, int> taxIndex; // This will have the first index of our hit

		for(int i = 0; i < record.size(); i++){ // Getting the taxa counts at our current taxonomic level. 6 = Genus -> 1 = Kingdom
			if(taxCounts.empty()){ // If it's empty, then we need to insert one
				taxCounts.insert({record.at(i)[t],record.at(i)[0]});
				taxIndex.insert({record.at(i)[t],i});
			}else{
				if(taxCounts.count(record.at(i)[t])){ // If key in the map
					int tmpCount = taxCounts.at(record.at(i)[t]) + record.at(i)[0];
					taxCounts[record.at(i)[t]] = tmpCount;

				}else{ // We need to insert the new taxa
					taxCounts.insert({record.at(i)[t],record.at(i)[0]});
					taxIndex.insert({record.at(i)[t],i});
				}
			}
		}

		// Finding if an element exists that's larger than our threshold
		auto it = find_if(taxCounts.begin(), taxCounts.end(), 
				[threshLCA](const pair<int,int> & x) -> bool  // This here is a lambda going through my map
				{return x.second >= threshLCA;}
				);

		//cerr << taxCounts.size() << "\t" << it->first << "\t" << it->second << endl;
		//sleep(2);
		// Testing if we found our hit
		if(taxCounts.count(it->first)){
			int key = it->first; // This is where the hit was found
			auto it = taxIndex.find(key);
			auto a = record.at(it->second); // This pull out the first 

			if(a[t] > 0){ // If the taxon is NA, we're not interested....
				for(int s = t + 1; s < 8; s++){ // Filling in the zeroes
					a[s] = 0;
				}
				return array8to7(a);
			}
		}
	}

	// Now if we've found nothing at all!
	array<int,7> struckOut;
	struckOut.fill(0); // Filling it with blanks
	return struckOut;
};

void readBlast(ifstream& file, float& thresh){
	vector<array<int,8>> blastFile; // The two dimensional array
	array<int,7> taxonResults; // The array output from LCA

	string entry, seqID, prevSeqID, tmp; // Entry is the line, seqID is the first ID and prevSeqID is the last one
	bool first = true;

	while (getline(file, entry)) // Iterate through the file
	{
		istringstream iss(entry); // read the line

		array<int,8> rec; // Create the record
		std::getline(iss, tmp, '\t'); // Getting the count
		rec[0] = stoi(tmp);
		std::getline(iss, seqID, '\t'); // Getting the sequence

		if(first){
			prevSeqID = seqID;
			first = false;
		}

		// If we've arrived at a new sequence, we need to process what we had previously
		if(seqID != prevSeqID){
	//		array<int,7> taxonResults; // The array output from LCA

			taxonResults = LCA(blastFile, thresh);
			printRecord(prevSeqID, taxonResults); // Printing the results
			blastFile.clear(); // Emptying the vector after we've used it
		}

		for(int i = 1; i < 8; i++){
			std::getline(iss, tmp, '\t'); // Getting the taxa
			rec[i] = catchNA(tmp);
		}
		
		blastFile.push_back(rec);
		prevSeqID = seqID;
	}
	
	taxonResults = LCA(blastFile, thresh);
	printRecord(prevSeqID, taxonResults); // Printing the results
}

int main(int argc, char* argv[]){ // The main program runs here
	string fileName; // The blast file name
	 float thresh; // The blast file name
	 if(argc != 3){ // Test if we have all the variables
		 usage();
		 return 1;
	 }
	 fileName = argv[1];
	 thresh = atof(argv[2]);
	 if(thresh > 1 | thresh <= 0.5f){ // Test if the threshold works
		 cout << "The threshold must be between 0.5 and 1" << endl;
		 usage();
		 return 1;

	 }

	 // Opening the file
	 ifstream blastFile;
	 blastFile.open(fileName);

	 if(blastFile){
		cout << "Sequence\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies" << endl;
	 	readBlast(blastFile,thresh); // Reading in the blast file
	 }else{
		 cout << "File: " << fileName << "doesn't exist!" << endl;
		 usage();
		 return 1;
	 }
	 


	 return 0;
}
