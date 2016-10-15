//============================================================================
// Name        : MEOIMap.cpp
//
//arg 1: Type of matching to use
//arg 2: Element(s) of interest file
//arg 3: Input file of reads
//arg 4: Output directory
//============================================================================


#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>

#include <seqan/index.h>
#include <seqan/seq_io.h>

#include "CGReadAndMapFilePair.hpp"
#include "CGRead.h"
#include "ReadFinders.h"

////////////////////////////////         ARGUMENTS
const int TYPE_ARG=1;
const int EOI_ARG=2;
const int INPUT_READS_ARG=3;
const int OUTPUT_DIR_ARG=4;
const int MIN_GAP1_ARG=5;
const int MIN_GAP2_ARG=7;
const int MIN_GAP3_ARG=9;
const int MAX_GAP1_ARG=6;
const int MAX_GAP2_ARG=8;
const int MAX_GAP3_ARG=10;

//TODO find a better place for this
/**std::ostream& operator<<(std::ostream& out, ReadFinders& readFinders){
	readFinders.writeFoundLocations(out);
	return out;
}**/

int main(int argc, const char * argv[]) {


	//TODO - verbosity and logging
	//These control the outputs - may create verbosity arg at later point

	const bool DEBUG = false;
	const bool VERBOSE = false;

	//verbose, debug and error output targets
	//  for relatively easy replacement later
	std::ostream& vout = std::cout;
	std::ostream& dout = std::cout;
	std::ostream& eout = std::cerr;

	//TODO end of really hack log stuff - change me

	//stat bins for summary
	long mapCounts[11]={0,0,0,0,0,0,0,0,0,0,0};

	std::string usage = "viraMap [Type of matching: Exact] [EOI file] [Input file of read files] [Output directory] ([min gap 1] [max gap 1] [min gap 2] [max gap 2] [min gap3] [max gap 3])";

	//Get system clock info for start/end duration
	//      from http://en.cppreference.com/w/cpp/chrono
	//TODO eclipse doesn't like these lines - blue borders to tell them from real errors

	std::chrono::time_point<std::chrono::system_clock> timeStart;
	std::chrono::time_point<std::chrono::system_clock> timeEnd;
	timeStart = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(timeStart);
	//TODO eclipse doesn't like these lines - blue borders to tell them from real errors

	if (VERBOSE) vout << "viraMap started at: "<< std::ctime(&start_time) << std::endl;


	// ****** CHECK ARG NUMBER

	if (argc != 5 && argc!=11){
		eout << "Incorrect number of arguments"<< std::endl;
		eout << usage << std::endl;
		exit(1);
	}

	// ****** ARG ******  Matching Type to use
	std::string matchType(argv[TYPE_ARG]);

	//TODO Figure out how to switch between modes.

	/*	if( matchType == "Exact"){
		typedef seqan::IndexEsa TIndexSpec;
	}
	else{
		typedef seqan::IndexEsa TIndexSpec;
	}*/


	// ****** ARG ******  EOI file
	//read pattern from file into index for searching
	seqan::SeqFileIn eoiFile;

	if (!open(eoiFile, argv[EOI_ARG]))
	{
		eout << "Could not open Element(s) of Interest file" <<argv[EOI_ARG]<<std::endl;
		exit(1);
	}

	if(DEBUG) dout << "EOI file looks good. Arg: "<< EOI_ARG<<std::endl;

	// ****** ARG ***** Input file of reads


	//fs::path p (argv[INPUT_READS_ARG]);
	//fs::path p ("/Users/BethLocke/Dropbox/Thesis/Software/viraMap/data/testReads.tsv");
	std::string p(argv[INPUT_READS_ARG]);
	std::string laneID;

	//should throw exception if does not
	//fs::exists(p);

	if( p.find("reads_GS") != std::string::npos ){
		//get the lane ID and .tsv, is after the "reads" and is 20 characters long (includes preceeding _)
		if(DEBUG) dout<<"--try to substr"<<std::endl;
		laneID = p.substr(  p.find("reads_GS")+5 ,20);
		if(DEBUG) dout<<"--LaneID found: "<<laneID<<std::endl;
	}

	else{
		eout<<"Input file does not follow reads_GS convention: ";
		eout<< p <<std::endl;
		exit(1);
	}

	CGReadAndMapFileTrio * rmFilePair;
	rmFilePair = new CGReadAndMapFileTrio(p);

	if(!rmFilePair->isOpen()){
		eout<<"File or mapping file could not be opened, or not a file: ";
		eout<< p <<std::endl;
		exit(1);
	}

	if(DEBUG) dout<<"Read and map files look good. Arg: "<< INPUT_READS_ARG<<std::endl;

	// ****** ARG **** Output directory

	std::string outprefix (argv[OUTPUT_DIR_ARG]);
	std::string outsuffix ("_viraMap.tsv");

	//add suffix different from input suffix - don't want to mistakenly overwrite originals!
	std::string readOutName = outprefix + "reads" + laneID + outsuffix;
	std::string mapOutName = outprefix + "mapping" + laneID + outsuffix;
	std::string locOutName = outprefix + "viraMapping" + laneID + outsuffix;
	std::string sumOutName = outprefix + "summary"+laneID + outsuffix;

	std::ofstream readOut(readOutName.c_str());
	std::ofstream mapOut(mapOutName.c_str());
	std::ofstream locOut(locOutName.c_str());
	std::ofstream sumOut(sumOutName.c_str());

	if(!readOut.is_open() || !mapOut.is_open() || !locOut.is_open() || !sumOut.is_open()){
		eout << laneID << ": Could not open output files: "<< readOutName <<" "<< mapOutName <<" "<<locOutName <<" "<<sumOutName <<" " <<std::endl;
		exit(1);
	}

	if(DEBUG) dout<<"Output directory looks good. Arg: "<< OUTPUT_DIR_ARG<<std::endl;


	//defaults
	int ng1 = -3;
	int ng2 = 0;
	int ng3 = 5;
	int xg1= -1;
	int xg2 = 2;
	int xg3 = 7;

	//*** Gap arguments
	if (argc == 11){
		std::istringstream ss(argv[MIN_GAP1_ARG]);
		if (!(ss >> ng1)){
			eout << "Invalid number for min gap 1 specified as argument: "<< MIN_GAP1_ARG<< " " << argv[MIN_GAP1_ARG] << '\n';
			exit(1);
		}
		ss.clear();
		ss.str(argv[MAX_GAP1_ARG]);
		if (!(ss >> xg1)){
			eout << "Invalid number for min gap 1 specified as argument: "<< MAX_GAP1_ARG<< " " << argv[MAX_GAP1_ARG] << '\n';
			exit(1);
		}
		ss.clear();
		ss.str(argv[MIN_GAP2_ARG]);
		if (!(ss >> ng2)){
			eout << "Invalid number for min gap 1 specified as argument: "<< MIN_GAP2_ARG<< " " << argv[MIN_GAP2_ARG] << '\n';
			exit(1);
		}
		ss.clear();
		ss.str(argv[MAX_GAP2_ARG]);
		if (!(ss >> xg2)){
			eout << "Invalid number for min gap 1 specified as argument: "<< MAX_GAP2_ARG<< " " << argv[MAX_GAP2_ARG] << '\n';
			exit(1);
		}
		ss.clear();
		ss.str(argv[MIN_GAP3_ARG]);
		if (!(ss >> ng3)){
			eout << "Invalid number for min gap 1 specified as argument: "<< MIN_GAP3_ARG<< " " << argv[MIN_GAP3_ARG] << '\n';
			exit(1);
		}
		ss.clear();
		ss.str(argv[MAX_GAP3_ARG]);
		if (!(ss >> xg3)){
			eout << "Invalid number for min gap 1 specified as argument: "<< MAX_GAP3_ARG<< " " << argv[MAX_GAP3_ARG] << '\n';
			exit(1);
		}
	}


	// ************************  END OF PARSE AND VERIFY ARGUMENTS

	sumOut << "viraMap started at: "<< std::ctime(&start_time) << std::endl;

	//create read finder
	ReadFinders readFinders(eoiFile, ng1,xg1,ng2,xg2,ng3,xg3);

	//int readFlag;
	int readID=0;
	const CGRead * read;

	sumOut << laneID<<": File: "<<p<< std::endl;
	sumOut << laneID<<"-Read in"<<std::endl;

	int tempFound4 =0;
	if(rmFilePair->isOpen()){

		if(VERBOSE) vout<<"read header: \n" << rmFilePair->readFileHeader();
		if(VERBOSE) vout<<"mapping header: \n" << rmFilePair->mapFileHeader();

		while ( rmFilePair->hasNext() ) { // keep reading until end-of-file

			//go to the next (or the first on the first time through)
			read = &(rmFilePair->next());

			//get ref to the next read
			//read = rmFilePair->read();

			//parse out half dnbs
			//readFlag = rmFilePair->readFlag(); //readFile >> flag;

			if(DEBUG) dout << "Read: "<< read->readSegment(2) <<std::endl;

			//TODO - add different modes if(EXACT ||  ){
			if (read->avgQualityLeft() >= 18 && read->avgQualityRight() >= 18) {

				//read = ;
				readFinders.searchForRead(read);

				//*******************************************************   Current

				if (readFinders.found()){

					//get read flag and start creating all output files
					//locFile << readFinders.foundLocations();

					//output read
					readOut << readID << "\t"<< rmFilePair->readEntry() <<std::endl;

					//output mapping
					const std::vector<std::string> mappingEntry = rmFilePair->mapEntry();
					for (std::vector<std::string>::const_iterator cit = mappingEntry.begin(); cit != mappingEntry.end(); ++cit){
						mapOut << readID <<"\t"<< cit->data()<<std::endl;
					}

					//output viral locations
					readFinders.writeFoundLocations(locOut, readID);

					//records read flag
					mapCounts[rmFilePair->readFlag()]++;

					tempFound4++;

				}

			}
			else{
				if(DEBUG) dout << "\tQuality too low" <<std::endl;
				//used to generate quality banks for testing
				//if( rmFilePair->readEntry().find("\'")==std::string::npos)
				//	std::cout << rmFilePair->readEntry() <<std::endl;
			}

			readID++; //always increment - this way reads are trackable between samples.
		}//end while file is good

	}//end if map and read are open
	else{
		sumOut<<laneID<<": Could not open one or both of map or read file" <<std::endl;
	}

	delete rmFilePair;


	//***********  Write to Summary file

	sumOut <<"Reads Found: " << readID << std::endl;
	sumOut <<"Map Flag Counts:"<<std::endl;
	sumOut << "0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10"<<std::endl;
	sumOut << mapCounts[0]<<"\t"<<
			mapCounts[1]<<"\t"<<
			mapCounts[2]<<"\t"<<
			mapCounts[3]<<"\t"<<
			mapCounts[4]<<"\t"<<
			mapCounts[5]<<"\t"<<
			mapCounts[6]<<"\t"<<
			mapCounts[7]<<"\t"<<
			mapCounts[8]<<"\t"<<
			mapCounts[9]<<"\t"<<
			mapCounts[10]<<"\t" <<std::endl;

	readFinders.writeProfile(sumOut);

	sumOut <<"Found the anchor ever: "<< readFinders.tempFoundL2Ever <<std::endl;
	sumOut <<"Found the anchor total: "<< readFinders.tempFoundL2Total <<std::endl;
	sumOut <<"Found all 4 left: "<< tempFound4 <<std::endl;

	//Find and report end time and total duration.
	//      from: http://en.cppreference.com/w/cpp/chrono
	//TODO eclipse doesn't like these lines - blue borders to tell them from real errors
	//
	timeEnd = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsedSeconds = timeEnd-timeStart;
	std::time_t end_time = std::chrono::system_clock::to_time_t(timeEnd);

	if (VERBOSE) std::cout << "finished computation at " << std::ctime(&end_time) << "elapsed time: " << elapsedSeconds.count() << "s\n";
	sumOut << "finished computation at " << std::ctime(&end_time) << "elapsed time: " << elapsedSeconds.count() << "s\n";
	//
	//TODO eclipse doesn't like these lines - blue borders to tell them from real errors

	sumOut <<"Analysis Complete for this read/map file pair"<<std::endl;
	sumOut.close();
	locOut.close();
	readOut.close();
	mapOut.close();
	return 0;
}


