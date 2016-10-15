/*
 * CGRead.h
 *
 *  Created on: Aug 10, 2016
 *      Author: BethLocke
 */

#ifndef CGREAD_H_
#define CGREAD_H_
#include <iostream>
#include <seqan/seq_io.h>

class CGRead {
public:
	CGRead();
	//CGRead(const CGRead&);
	CGRead(std::string read, std::string qual);
	//CGRead(char lineRead[]);

	int avgQuality() const;
	int avgQualityLeft() const;
	int avgQualityRight() const;

	//getters from read file entry
	seqan::String<char> readSegment(int i) const;
	seqan::String<char> maskedReadSegment( int i, int q) const;

	virtual ~CGRead();

private:
	std::string _read;
	std::string _qual;
};

#endif /* CGREAD_H_ */
