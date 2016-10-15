/*
 * CGReadFile.hpp
 *
 *  Created on: Aug 21, 2012
 *      Author: BethLocke
 */

#ifndef CGREADANDMAPFILEPAIR_HPP_
#define CGREADANDMAPFILEPAIR_HPP_


#include <iostream>
#include <fstream>

#include <seqan/sequence.h>

#include "CGRead.h"


class CGReadAndMapFileTrio {

public:
	CGReadAndMapFileTrio(std::string readFName);
	~CGReadAndMapFileTrio();

	bool isGood();
	bool isOpen();
	bool hasNext();

	std::string readFileHeader();
	std::string mapFileHeader();

	const CGRead& next();

	const CGRead& read(){ return _read; }

	int readFlag();

	//string versions of the whole entry in both read and map files.
	//map file entry may be multiple lines.
	const std::string& readEntry();
	const std::vector<std::string>& mapEntry();

private:
	std::string _readHeader;
	std::string _mapHeader;

	int _readFlag;
	int _mapFlag;

	CGRead _read;
	std::string _readEntry;
	std::vector<std::string> _mapEntry;

	std::ifstream _readFile;
	std::ifstream _mapFile;

};


#endif /* CGREADFILE_HPP_ */
