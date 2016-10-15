/*
 * ReadFinder.h
 *
 *  Created on: Aug 9, 2016
 *      Author: BethLocke
 */

#ifndef READFINDERS_H_
#define READFINDERS_H_

#include <seqan/index.h>

#include "CGReadAndMapFilePair.hpp"
#include "ReadLocation.h"

typedef seqan::CharString TID;
typedef seqan::StringSet<TID> TIDSet;
typedef seqan::Iterator<TIDSet>::Type TIDIterator;

typedef seqan::Dna5String TInputGenome;
typedef seqan::StringSet<TInputGenome> THaystacks;
typedef seqan::Index< THaystacks, seqan::IndexEsa<> > TIndex;
typedef seqan::Finder< TIndex > TFinder;

class ReadFinders {
public:

	const static int NUMFINDERS=8;

	//Temp used for testing - generates stats.
	int tempFoundL2Ever =0;
	int tempFoundL2Total =0;

	ReadFinders(seqan::SeqFileIn& haystacksFile);
	ReadFinders(seqan::SeqFileIn& haystacksFile, int minGap1,int maxGap1,int minGap2,int maxGap2,int minGap3,int maxGap3);

	void searchForRead(const CGRead * read);
	bool found();

	const std::vector<ReadLocation>& foundLocations();
	void writeFoundLocations(std::ostream& out, int readID);
	void writeProfile(std::ostream& out);

	virtual ~ReadFinders();

private:
	void clear();
	void init(seqan::SeqFileIn& haystacksFile);
	void doSearch(const CGRead * read, bool isForward);

	//These control the outputs - may create verbosity arg at later point
	//TODO - verbosity
	const bool DEBUG = false;
	const bool VERBOSE = false;
	const bool PROFILE = true;

	//verbose, debug and error output targets
	//  for relatively easy replacement later
	std::ostream& vout = std::cout;
	std::ostream& dout = std::cout;
	//std::ofstream dout;// = std::ofstream();
	std::ostream& eout = std::cout;

	const int L0=0;
	const int L1=1;
	const int L2=2;
	const int L3=3;

	const int R0=7;
	const int R1=6;
	const int R2=5;
	const int R3=4;

	int _min_gap1;
	int _min_gap2;
	int _min_gap3;

	int _max_gap1;
	int _max_gap2;
	int _max_gap3;

	TIDSet _ids;
	std::map<TID, long> _leftFoundIn;
	std::map<TID, long> _rightFoundIn;
	THaystacks _haystacks;
	TIndex _index;
	THaystacks _rcHaystacks;

	TFinder* _finders[NUMFINDERS];
	TFinder* _rcFinders[NUMFINDERS];

	//std::string _read;

	std::vector< ReadLocation > _foundReadLocations;


	bool _foundLeft=false;
	bool _foundRight=false;


};


#endif /* READFINDERS_H_ */
