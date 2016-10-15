/*
 * CGRead.cpp
 *
 *  Created on: Aug 10, 2016
 *      Author: BethLocke
 */

#include "CGRead.h"



CGRead::CGRead(){
	_read=std::string();
	_qual=std::string();
}

/*
 * This is what it should be doing
 * CGRead::CGRead(const CGRead& toCopy) {
	//just copy the strings (what would happen automatically anyway
	_read = toCopy._read;
	_qual = toCopy._qual;
}*/

CGRead::CGRead(std::string read, std::string qual) {
	_read = read;
	_qual = qual;
}

// more efficent
//follows read Entry format from the CG file
//Doesn't allow 10s
/*CGRead::CGRead(char readEntry[]) {
	_read = std::string(readEntry[2],70);
	_qual = std::string(readEntry[73],70);
}*/

int CGRead::avgQuality() const{
	int avg=0;

	for (int i=0;i<(int)_qual.length();i++){
		char c = _qual.at(i);
		avg += ((int)c - 33);
	}
	avg = avg/70;

	return avg;
}

int CGRead::avgQualityLeft() const{
	int avg=0;
	int l = (int)_qual.length();
	for (int i=0;i<l/2;i++){
		avg += ((int)_qual.at(i) - 33);
	}

	avg = avg/(l/2);

	return avg;
}
int CGRead::avgQualityRight() const{
	int avg=0;
	int l = (int)_qual.length();
	for (int i=l/2;i<l;i++){
		avg += ((int)_qual.at(i) - 33);
	}
	avg = avg/(l/2);

	return avg;
}

seqan::String<char> CGRead::readSegment(int i) const{
	if (_read.length() !=70)
		return "";


	switch(i){
	case 0:
		return _read.substr(0,5);
	case 1:
		return _read.substr(5,10);
	case 2:
		return _read.substr(15,10);
	case 3:
		return _read.substr(25,10);
	case 4:
		return _read.substr(35,10);
	case 5:
		return _read.substr(45,10);
	case 6:
		return _read.substr(55,10);
	case 7:
		return _read.substr(65,10);
	}
	//will default to the empty string
	return "";
}

seqan::String<char> CGRead::maskedReadSegment( int i, int q) const{
	std::string read = "";

	switch(i){
	case 0:
		read = _read.substr(0,5);
		for (int i=0;i < 5;i++){
			if((int)_qual.at(i) < q + 33){
				read.at(i)='.';
			}

		}
		return read;
	case 1:
		read = _read.substr(5,10);
		for (int i=5;i < 15;i++){
			if((int)_qual.at(i) < q + 33){
				read.at(i-5)='.';
			}

		}
		return read;
	case 2:


		read = _read.substr(15,10);
		for (int i=15;i < 25;i++){
			if((int)_qual.at(i) < q + 33){
				read.at(i-15)='.';
			}

		}
		return read;
	case 3:

		read = _read.substr(25,10);
		for (int i=25;i < 35;i++){
			if((int)_qual.at(i) < q + 33){
				read.at(i-25)='.';
			}

		}
		return read;
	case 4:
		read = _read.substr(35,10);
		for (int i=35;i < 45;i++){
			if((int)_qual.at(i) < q + 33){
				read.at(i-35)='.';
			}

		}
		return read;
	case 5:
		read = _read.substr(45,10);
		for (int i=45;i < 55;i++){
			if((int)_qual.at(i) < q + 33){
				read.at(i-45)='.';
			}

		}
		return read;
	case 6:
		read = _read.substr(55,10);
		for (int i=55;i < 65;i++){
			if((int)_qual.at(i) < q + 33){
				read.at(i-55)='.';
			}

		}
		return read;
	case 7:

		read = _read.substr(65,5);
		for (int i=65;i < 70;i++){
			if((int)_qual.at(i) < q + 33){
				read.at(i-65)='.';
			}

		}
		return read;
	}
	//will default to the empty string
	return read;
}



CGRead::~CGRead() {
	// TODO Auto-generated destructor stub
}

