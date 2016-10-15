/*
 * ReadFinder.cpp
 *
 *  Created on: Aug 9, 2016
 *      Author: BethLocke
 */

#include "ReadFinders.h"

/**
 * Reads the file and creates an index of from the sequences from it.
 * It will use exact matching.
 */
ReadFinders::ReadFinders(seqan::SeqFileIn& haystacksFile){
	//:dout("/Users/BethLocke/Dropbox/Thesis/Software/viraMap/data/testOutput/test_debugFinder.txt") {


	//defaults from white papers CG
	_min_gap1 = -3;
	_min_gap2 = 0;
	_min_gap3 = 5;

	_max_gap1 = -1;
	_max_gap2 = 2;
	_max_gap3 = 7;

	init(haystacksFile);

}


ReadFinders::ReadFinders(seqan::SeqFileIn& haystacksFile, int minGap1,
		int maxGap1, int minGap2, int maxGap2, int minGap3, int maxGap3){
	//	//:dout("/Users/BethLocke/Dropbox/Thesis/Software/viraMap/data/testOutput/test_debugFinder.txt") {

	_min_gap1 = minGap1;
	_min_gap2 = minGap2;
	_min_gap3 = minGap3;

	_max_gap1 = maxGap1;
	_max_gap2 = maxGap2;
	_max_gap3 = maxGap3;

	init(haystacksFile);
}

void ReadFinders::init(seqan::SeqFileIn& haystacksFile) {


	//if using output file for debug stream
	/*if(!dout.is_open()){
		eout << "Could not open test log file "<<std::endl;
		exit(1);
	}*/

	//*****  Read the input sequence file into the index
	try
	{
		//reads all records from file into the two string sets
		seqan::readRecords(_ids, _haystacks, haystacksFile);
	}

	catch (seqan::Exception const & e)
	{
		eout << "ERROR: " << e.what() << std::endl;
		exit(1);
	}

	//print all ids
	/*for( unsigned int i =0; i< seqan::length(_ids); ++i){
		std::cout << i << " "<<_ids[i] <<std::endl;
	}*/

	_rcHaystacks = _haystacks;
	seqan::reverseComplement(_rcHaystacks);

	//build the index from the haystacks
	TIndex findex(_haystacks);

	for (int i=0;i< NUMFINDERS; ++i )
		// - not sure why the Finder doesn't take the index but the haystacks but following docs: http://seqan.readthedocs.io/en/master/Tutorial/Algorithms/PatternMatching/IndexedPatternMatching.html
		_finders[i] = new TFinder(_haystacks);

	//build the index from the reverse complimented haystacks
	TIndex rcindex(_rcHaystacks);

	for (int i=0;i< NUMFINDERS; ++i )
		// - not sure why the Finder doesn't take the index but the haystacks but following docs: http://seqan.readthedocs.io/en/master/Tutorial/Algorithms/PatternMatching/IndexedPatternMatching.html
		_rcFinders[i] = new TFinder(_rcHaystacks);


	if(PROFILE){
		for (TIDIterator it = seqan::begin(_ids); it != seqan::end(_ids); ++it)
		{
			_leftFoundIn.insert(std::pair<TID,long>((*it) ,0));
			_rightFoundIn.insert(std::pair<TID, long>((*it),0));
		}
	}
	if(DEBUG) dout<<"Input sequence file read into index"<<std::endl;

	/* NOTE: Indices in SeqAn are built on demand.
	 * That means that the index tables are not build when the constructor is called,
	 * but when we search for a pattern for the first time.
	 */
}


/**
 * Will look for the given read in the haystacks.
 *
 */
//TODO Throw exception
void ReadFinders::searchForRead(const CGRead * read) {
	this->clear();
	doSearch(read, true);
	doSearch(read, false);
}


//finders must match
void ReadFinders::doSearch(const CGRead * read, bool isForward) {

	//teehee
	TFinder** findersToUse;
	if (isForward)
		findersToUse = _finders;
	else
		findersToUse = _rcFinders;


	std::map<int, std::vector<int>> foundSegmentLocations[8];
	//for(int i=0;i< length(_ids);++i){

	//}
	//	std::multimap<int,int> foundSegmentLocations[8];

	// Used for finding locations

	//will iterate over all l2 entries
	std::map<int, std::vector<int>>::iterator l2_key_it;
	//will only iterate over matching vector at a certain key
	std::vector<int>::iterator l3_it;
	std::vector<int>::iterator l2_it;
	std::vector<int>::iterator l1_it;
	std::vector<int>::iterator l0_it;

	//std::pair <std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> l3_it_bounds;
	//std::pair <std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> l1_it_bounds;
	//std::pair <std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> l0_it_bounds;

	//keep track of found gaps
	int foundGap1=0;
	int foundGap2=0;
	int foundGap3=0;

	bool popL0=false;
	bool popR0=false;

	// LEFT

	/*
	 *     L0       L1       L2      L3
	 *         gap1     gap2    gap3
	 */
	if( seqan::find(*(findersToUse[L2]), read->readSegment(L2)) ){
		tempFoundL2Ever++;
		tempFoundL2Total++;
		if(DEBUG) dout << "\tFOUND [" << beginPosition(*(findersToUse[L2])) << ',' << endPosition(*(findersToUse[L2])) << ")\t" << infix(*(findersToUse[L2])) << std::endl;

		//read out of the Position pair into the multimap for searching later

		foundSegmentLocations[L2][getValueI1(beginPosition(*(findersToUse[L2])))].push_back( getValueI2(beginPosition(*(findersToUse[L2]))) );

		//Load every location where this read is found
		while( seqan::find(*(findersToUse[L2]), read->readSegment(L2)) ){

			foundSegmentLocations[L2][getValueI1(beginPosition(*(findersToUse[L2]))) ].push_back( getValueI2(beginPosition(*(findersToUse[L2]))) );

			tempFoundL2Total++;
			if(DEBUG) dout << ".";
		}
		if(DEBUG) dout << std::endl;

		//check other finders
		if ( seqan::find( *(findersToUse[L3]), read->readSegment(L3)) &&
				seqan::find( *(findersToUse[L1]), read->readSegment(L1)) &&
				seqan::find(*(findersToUse[L0]), read->readSegment(L0)) )
		{
			//std::cout << "Found l3 l1 l0 ";
			//add the first found locations
			foundSegmentLocations[L3][getValueI1(beginPosition(*(findersToUse[L3]))) ].push_back( getValueI2(beginPosition(*(findersToUse[L3]))) );
			foundSegmentLocations[L1][getValueI1(beginPosition(*(findersToUse[L1]))) ].push_back( getValueI2(beginPosition(*(findersToUse[L1]))) );
			foundSegmentLocations[L0][getValueI1(beginPosition(*(findersToUse[L0]))) ].push_back( getValueI2(beginPosition(*(findersToUse[L0]))) );

			//Load every other location where this read is found
			while( seqan::find(*(findersToUse[L3]), read->readSegment(L3)) ){
				foundSegmentLocations[L3][getValueI1(beginPosition(*(findersToUse[L3]))) ].push_back( getValueI2(beginPosition(*(findersToUse[L3]))) );
				if(DEBUG) dout << ".";
			}
			//Load every other location where this read is found
			while( seqan::find(*(findersToUse[L1]), read->readSegment(L1)) ){
				foundSegmentLocations[L1][getValueI1(beginPosition(*(findersToUse[L1]))) ].push_back( getValueI2(beginPosition(*(findersToUse[L1]))) );
				if(DEBUG) dout << ".";
			}


			//LEAVE OUT L0 - as it is tiny and will be found many times - only do as verification if 3 are found in conjuction

			if(DEBUG) dout << std::endl;

			//Output all found locations:
			if(DEBUG){
				std::map<int,std::vector<int>>::iterator it;
				dout << "Locations found:\n\tL0 -forward?"<<isForward<<"\n\t\t" <<foundSegmentLocations[L0].size()<<std::endl;

				//for (it=foundLocations[L0].begin(); it!=foundLocations[L0].end(); ++it)
				//	dout << "\t\t"<<(*it).first << " : " << (*it).second << '\n';
				dout << "Locations found:\n\tL1:\n";
				for (it=foundSegmentLocations[L1].begin() ; it!=foundSegmentLocations[L1].end(); ++it)
					dout << "\t\t"<< it->first<< std::endl;// << " : " << it->second.size() << '\n';

				dout << "Locations found:\n\tL2:\n";
				for (it=foundSegmentLocations[L2].begin(); it!=foundSegmentLocations[L2].end(); ++it)
					dout << "\t\t"<<it->first << std::endl;//<< " : " << it->second.size() << '\n';

				dout << "Locations found:\n\tL3:\n";
				for (it=foundSegmentLocations[L3].begin(); it!=foundSegmentLocations[L3].end(); ++it)
					dout << "\t\t"<<it->first << std::endl;//<< " : " << it->second.size() << '\n';
			}

			/*for(l2_it = _foundSegmentLocations[L2].begin(); l2_it!=_foundSegmentLocations[L2].end(); ++l2_it){
				std::cout << l2_it->first << " " << l2_it->second <<std::endl;
			}*/
			for(l2_key_it = foundSegmentLocations[L2].begin(); l2_key_it != foundSegmentLocations[L2].end(); ++l2_key_it){
				//check if the other two have an entry for that key
				int id_key = l2_key_it->first;

				//Quick sanity check
				//the first element of the pair is the start
				l3_it = foundSegmentLocations[L3][id_key].begin();
				l1_it = foundSegmentLocations[L1][id_key].begin();

				//if our iterators have anything in them, then there is a set of results for that key in the multimap
				if( l3_it != foundSegmentLocations[L3][id_key].end() && l1_it != foundSegmentLocations[L1][id_key].end()){


					//check all locations for this key in l2
					for (l2_it = l2_key_it->second.begin(); l2_it != l2_key_it->second.end(); ++l2_it){

						//restart looking at this id
						l3_it = foundSegmentLocations[L3][id_key].begin();
						l1_it = foundSegmentLocations[L1][id_key].begin();

						int l2_start = (*l2_it);

						//we have a location in the same org
						if(DEBUG) dout <<"Found putative L!"<<std::endl;

						int l1_found=-1;

						//until it is found or we run out - the second element of the boundary pair
						for (;l1_found < 0 && l1_it != foundSegmentLocations[L1][id_key].end(); ++l1_it){
							int l1_end = (*l1_it) + 10;
							foundGap2 = l2_start-l1_end;
							if (foundGap2 <= _max_gap2 && foundGap2 >=_min_gap2)
								l1_found=(*l1_it);
						}
						if(DEBUG && l1_found>0) dout << "**FOUND MATCHING L1"<<std::endl;

						int l3_found=-1;

						//until it is found or we run out
						for (;l3_found < 0 && l3_it != foundSegmentLocations[L3][id_key].end(); ++l3_it){
							int l2_end = l2_start + 10;
							int l3_start = (*l3_it);
							foundGap3 = l3_start-l2_end;
							if (foundGap3 <= _max_gap3 && foundGap3 >= _min_gap3)
								l3_found=l3_start;
						}

						if(DEBUG && l3_found>0) dout << "**FOUND MATCHING L3"<<std::endl;

						int l0_found=-1;
						if(l1_found>0 && l3_found>0) {

							//if we haven't already, finish loading every location where this 5mer read is found
							if(!popL0){
								//foundSegmentLocations[L0].size() == 1)  Does not work
								popL0=true;
								while( seqan::find(*(findersToUse[L0]), read->readSegment(L0)) ){
									foundSegmentLocations[L0][ getValueI1(beginPosition(*(findersToUse[L0]))) ].push_back( getValueI2(beginPosition(*(findersToUse[L0]))) );
								}
							}

							//verify there is an l0

							//the first element of the pair is the start
							l0_it = foundSegmentLocations[L0][id_key].begin();

							for (;l0_found < 0 && l0_it != foundSegmentLocations[L0][id_key].end();
									++l0_it){
								int l0_end = (*l0_it) + 5;
								foundGap1 = l1_found-l0_end;
								if (foundGap1 <= _max_gap1 && foundGap1 >=_min_gap1){
									l0_found=(*l0_it);
									//std::cout << id_key<< " - " << l0_it->first <<": "<< l0_it->second <<std::endl;
								}
							}

							//FOUND EVERYTHING
							if(l0_found>0){
								_foundLeft=true;

								//add to found locations
								std::string foundID = std::string(toCString(_ids[id_key]));
								//std::cout << isForward<< "- " <<id_key << " -- "<<foundID << " " << l0_found <<std::endl;
								if (!isForward)
									//count in from the last index length -1 (not the full length)
									l0_found = length(_haystacks[id_key]) -1- l0_found;
								_foundReadLocations.push_back( ReadLocation( foundID, true,isForward,l0_found,foundGap1, foundGap2, foundGap3));

								if(PROFILE) this->_leftFoundIn[_ids[id_key]]++;
								if(DEBUG) dout << "******* FOUND FULL MATCH"<<std::endl;
							}
						}
					}
				}
			}
		}
	}
	else{
		if(DEBUG) dout<< "\tLeft arm not found. Is forward? "<<isForward<<std::endl;
	}
	//std::cout<<std::endl;

	// RIGHT
	/*
	 *     R3       R2       R1      R0
	 *         gap3     gap2    gap1
	 */


	//will iterate over all l2 entries
	std::map<int, std::vector<int>>::iterator r2_key_it;
	//will only iterate over matching vector at a certain key
	std::vector<int>::iterator r3_it;
	std::vector<int>::iterator r2_it;
	std::vector<int>::iterator r1_it;
	std::vector<int>::iterator r0_it;

	//keep track of found gaps
	foundGap1=0;
	foundGap2=0;
	foundGap3=0;


	if( seqan::find(*(findersToUse[R2]), read->readSegment(R2)) ){

		if(DEBUG) dout << "\tFOUND [" << beginPosition(*(findersToUse[R2])) << ',' << endPosition(*(findersToUse[R2])) << ")\t" << infix(*(findersToUse[R2])) << std::endl;
		//std::cout << "Found r2 ";
		//read out of the Position pair into the multimap for searching later
		foundSegmentLocations[R2][ getValueI1(beginPosition(*(findersToUse[R2]))) ].push_back( getValueI2(beginPosition(*(findersToUse[R2]))) );

		//Load every location where this read is found
		while( seqan::find(*(findersToUse[R2]), read->readSegment(R2)) ){
			foundSegmentLocations[R2][getValueI1(beginPosition(*(findersToUse[R2]))) ].push_back( getValueI2(beginPosition(*(findersToUse[R2]))) );
			if(DEBUG) dout << ".";
		}
		if(DEBUG) dout << std::endl;
		//check other finders
		if ( seqan::find( *(findersToUse[R3]), read->readSegment(R3)) &&
				seqan::find( *(findersToUse[R1]), read->readSegment(R1)) &&
				seqan::find(*(findersToUse[R0]), read->readSegment(R0)) )
		{
			//std::cout << "Found r3 r1 r0 ";
			//add the first found locations
			foundSegmentLocations[R3][getValueI1(beginPosition(*(findersToUse[R3]))) ].push_back( getValueI2(beginPosition(*(findersToUse[R3]))) );
			foundSegmentLocations[R1][getValueI1(beginPosition(*(findersToUse[R1]))) ].push_back( getValueI2(beginPosition(*(findersToUse[R1]))) );
			foundSegmentLocations[R0][getValueI1(beginPosition(*(findersToUse[R0]))) ].push_back( getValueI2(beginPosition(*(findersToUse[R0]))) );

			//Load every other location where this read is found
			while( seqan::find(*(findersToUse[R3]), read->readSegment(R3)) ){
				foundSegmentLocations[R3][ getValueI1(beginPosition(*(findersToUse[R3]))) ].push_back( getValueI2(beginPosition(*(findersToUse[R3]))) );
				if(DEBUG) dout << ".";
			}
			//Load every other location where this read is found
			while( seqan::find(*(findersToUse[R1]), read->readSegment(R1)) ){
				foundSegmentLocations[R1][ getValueI1(beginPosition(*(findersToUse[R1]))) ].push_back( getValueI2(beginPosition(*(findersToUse[R1]))) );
				if(DEBUG) dout << ".";
			}
			//LEAVE OUT L0 - as it is tiny and will be found many times - only do as verification if 3 are found in conjuction

			if(DEBUG) dout << std::endl;

			//Output all found locations:
			if(DEBUG){
				std::map<int,std::vector<int>>::iterator it;
				dout << "Locations found:\n\tR0 -forward?"<<isForward<<":\n\t\t" <<foundSegmentLocations[R0].size()<<std::endl;

				//for (it=foundLocations[L0].begin(); it!=foundLocations[L0].end(); ++it)
				//	dout << "\t\t"<<(*it).first << " : " << (*it).second << '\n';
				dout << "Locations found:\n\tR1:\n";
				for (it=foundSegmentLocations[R1].begin(); it!=foundSegmentLocations[R1].end(); ++it)
					dout << "\t\t"<<it->first << std::endl;//" : " << it->second << '\n';
				dout << "Locations found:\n\tR2:\n";
				for (it=foundSegmentLocations[R2].begin(); it!=foundSegmentLocations[R2].end(); ++it)
					dout << "\t\t"<<it->first << std::endl;// << it->second << '\n';
				dout << "Locations found:\n\tR3:\n";
				for (it=foundSegmentLocations[R3].begin(); it!=foundSegmentLocations[R3].end(); ++it)
					dout << "\t\t"<<it->first  << std::endl;//<< " : " << it->second << '\n';
			}


			for(r2_key_it = foundSegmentLocations[R2].begin(); r2_key_it!=foundSegmentLocations[R2].end(); ++r2_key_it){
				//check if the other two have an entry for that key
				int id_key = r2_key_it->first;


				//the first element of the pair is the start
				r3_it = foundSegmentLocations[R3][id_key].begin();
				r1_it = foundSegmentLocations[R3][id_key].begin();


				//if our iterators have anything in them, then there is a set of results for that key in the multimap
				if( r3_it != foundSegmentLocations[R3][id_key].end() && r1_it != foundSegmentLocations[R1][id_key].end()){



					//check all locations for this key in l2
					for (r2_it = r2_key_it->second.begin(); r2_it != r2_key_it->second.end(); ++r2_it){

						//restart looking at this id
						r3_it = foundSegmentLocations[R3][id_key].begin();
						r1_it = foundSegmentLocations[R1][id_key].begin();

						int r2_start = (*r2_it);

						int r2_end = r2_start + 10;

						//we have a location in the same org
						if(DEBUG) dout <<"Found putative!"<<std::endl;


						int r1_found=-1;

						//until it is found or we run out
						for (;r1_found < 0 && r1_it != foundSegmentLocations[R1][id_key].end(); ++r1_it){

							foundGap2 = (*r1_it) -r2_end;
							if (foundGap2 <= _max_gap2 && foundGap2 >=_min_gap2)
								r1_found=(*r1_it);
						}

						if(DEBUG && r1_found>0) dout << "**FOUND MATCHING R1"<<std::endl;

						int r3_found=-1;

						//until it is found or we run out
						for (;r3_found < 0 && r3_it != foundSegmentLocations[R3][id_key].end(); ++r3_it){
							int r3_end = (*r3_it) + 10;
							foundGap3 = r2_start-r3_end;
							if (foundGap3 <= _max_gap3 && foundGap3 >= _min_gap3)
								r3_found=(*r3_it);
						}

						if(DEBUG && r3_found>0) dout << "**FOUND MATCHING R3"<<std::endl;

						int r0_found=-1;
						if(r1_found>0 && r3_found>0) {

							//At this point, if we haven't already, finish loading every location where this 5mer read is found
							if(!popR0){
								//foundSegmentLocations[R0].size() == 1){  does not work
								popR0=true;
								while( seqan::find(*(findersToUse[R0]), read->readSegment(R0)) ){
									foundSegmentLocations[R0][getValueI1(beginPosition(*(findersToUse[R0]))) ].push_back( getValueI2(beginPosition(*(findersToUse[R0]))) );
								}
							}
							//verify there is an r0


							//the first element of the pair is the start
							r0_it = foundSegmentLocations[R0][id_key].begin();
							for (;r0_found < 0 && r0_it !=foundSegmentLocations[R0][id_key].end();
									++r0_it){
								int r1_end =  r1_found + 10;
								foundGap1 = (*r0_it) - r1_end;
								if (foundGap1 <= _max_gap1 && foundGap1 >=_min_gap1)
									r0_found=(*r0_it);
							}

							//FOUND EVERYTHING
							if(r0_found>0){
								_foundRight=true;

								//add to found locations
								std::string foundID = std::string(toCString(_ids[id_key]));
								//std::cout << isForward<< "- " <<id_key << " -- "<<foundID << " " << r3_found <<std::endl;
								if (!isForward)
									//count in from the last index length -1 (not the full length)
									r3_found = length(_haystacks[id_key]) -1- r3_found;
								_foundReadLocations.push_back( ReadLocation( foundID, false, isForward,r3_found,foundGap1, foundGap2, foundGap3));
								if(PROFILE) this->_rightFoundIn[_ids[id_key]]++;
								if(DEBUG) dout << "******* FOUND FULL MATCH"<<std::endl;


							}
						}
					}

				}
			}

		}
	}
	else{
		if(DEBUG) dout<< "\tRight arm not found.  Is forward? "<<isForward<<std::endl;
	}
	//std::cout<<std::endl;
}

/**
 * Returns true if either arm of the read was found
 */
bool ReadFinders::found() {
	return _foundLeft||_foundRight;
}


void ReadFinders::clear() {
	////clear the finders and location multimaps
	for (int i=0; i<NUMFINDERS; ++i){
		seqan::clear(*(_finders[i]));
		seqan::clear(*(_rcFinders[i]));
	}
	_foundReadLocations.clear();
	_foundLeft=false;
	_foundRight=false;
}

void ReadFinders::writeFoundLocations(std::ostream& out, int readID){

	for (std::vector<ReadLocation>::iterator it= _foundReadLocations.begin(); it< _foundReadLocations.end(); ++it){
		out << readID<<"\t"
				<< it->isLeft() << "\t"
				<< it->isForward() << "\t"
				<< it->getViraId() << "\t"
				<< it->getStart() << "\t"
				<< it->getGap1() << "\t"
				<< it->getGap2() << "\t"
				<< it->getGap3() << "\t"
				<< std::endl;
	}

}


void ReadFinders::writeProfile(std::ostream& out) {
	if(PROFILE){
		out<<"ID"<<std::endl;
		for (TIDIterator it = seqan::begin(_ids); it != seqan::end(_ids); ++it)
			out << "\t"<< *it ;
		out<<std::endl<<"Left Found Count";
		for (TIDIterator it = seqan::begin(_ids); it != seqan::end(_ids); ++it)
			out << "\t"<< _leftFoundIn[*it] ;
		out<<std::endl<<"Right Found Count";
		for (TIDIterator it = seqan::begin(_ids); it != seqan::end(_ids); ++it)
			out << "\t"<< _rightFoundIn[*it] ;
		out<<std::endl;
	}
	else{
		out<<"Profile flag in ReadFinders.h not set to true - no profile generated"<<std::endl;
	}
}



const std::vector<ReadLocation>& ReadFinders::foundLocations() {
	return _foundReadLocations;
}


/**
 * Deconstructor
 *  Finders allocated in constructor and cleared each time a new search initiated
 */
ReadFinders::~ReadFinders() {
	for (int i=0; i< NUMFINDERS;++i){
		delete _finders[i];
	}
}
