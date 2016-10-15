/*
 * CGReadAndMapFilePair.cpp
 *
 *  Created on: Aug 21, 2012
 *      Author: BethLocke
 */

#include <string>
#include "CGReadAndMapFilePair.hpp"

CGReadAndMapFileTrio::CGReadAndMapFileTrio(std::string readFilePath){
	_readFile.open(readFilePath.c_str());
	//std::cout<<"reading: " <<readFilePath<<std::endl;
	//get corresponding map file name
	int readPos = readFilePath.find("reads_GS");
	std::string mapFilePath(readFilePath.substr(0,readPos));
	mapFilePath+="mapping";
	mapFilePath+=readFilePath.substr(readPos+5);

	//read out headers on both input files.
	_mapFile.open(mapFilePath.c_str());
	//std::cout<<mapFilePath<<std::endl;
	char lineRead[256]; //used to read in header lines and rest of line after flag read out.

	//read out the headers
	if(_mapFile.is_open() && _readFile.is_open()){

		//read first 15lines (header)
		for(int i=0; i<15;i++){
			_readFile.getline(lineRead, 256);
			_readHeader.append(lineRead);
			_readHeader.append("\n");
		}
		for(int i=0;i<13;++i){
			_mapFile.getline(lineRead,256);
			_mapHeader.append(lineRead);
			_mapHeader.append("\n");
		}
	}

}

CGReadAndMapFileTrio::~CGReadAndMapFileTrio(){
	_readFile.close();
	_mapFile.close();
}

bool CGReadAndMapFileTrio::isGood(){
	//check status of both files
	return ( _readFile.good() && _mapFile.good() );
}

bool CGReadAndMapFileTrio::isOpen(){
	//check status of both files
	return ( _readFile.is_open() && _mapFile.is_open() );
}

bool CGReadAndMapFileTrio::hasNext(){
	//peek for an int - check if next char is a valid digit
	int nextFlag =  _readFile.peek()-48;
	return isGood() && nextFlag>=0 && nextFlag<=9;
}


std::string CGReadAndMapFileTrio::readFileHeader(){
	return _readHeader;
}


std::string CGReadAndMapFileTrio::mapFileHeader(){
	return _mapHeader;
}

//advance to the next entry in both read and map files
const CGRead& CGReadAndMapFileTrio::next(){
	//the file pointers are already at the position of the next read
	char lineRead[256];

	_readEntry.clear();
	_mapEntry.clear();

	//read the current read entry from the read file
	//_readFlag=_readFile.peek()-48;  // misses the 10 flag
	_readFile >> _readFlag; //reads whatever int it finds
	_readFile.getline(lineRead,256);

	//should be predictably the tab followed by read followed by tab
	std::string toParse(lineRead);

	//reconstruct the original line for posterity(output to file again later)

	//c++ 11
	_readEntry += std::to_string(_readFlag);
	_readEntry += lineRead;

	//the read comes after the flag character and a tab, and is 70 bp long
	_read = CGRead(toParse.substr(1,70),toParse.substr(72,70));

	//reset to default values
	_mapFlag = -1;

	//if the read flag is 5,6,9 or 10, there will be no mapping entries, and we should not read anything from the mapping file.
	//as per CG documentation
	if (_readFlag != 5 && _readFlag != 6 && _readFlag != 9 && _readFlag != 10 ){
		//there is something to read in, until an odd flag is reached.
		while((_mapFile.peek()-48) %2 ==0){
			_mapFile.getline(lineRead,256);
			_mapEntry.push_back(lineRead);
		}
		//read in the last entry
		//mapflag should always be odd
		_mapFlag = _mapFile.peek()-48;
		_mapFile.getline(lineRead,256);
		_mapEntry.push_back(lineRead);
	}
	return _read;
}



int CGReadAndMapFileTrio::readFlag(){
	return _readFlag;
}


const std::string& CGReadAndMapFileTrio::readEntry(){
	return _readEntry;
}

const std::vector<std::string>& CGReadAndMapFileTrio::mapEntry(){
	return _mapEntry;
}

//end of file
