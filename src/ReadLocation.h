/*
 * ReadLocation.h
 *
 *  Created on: Aug 30, 2016
 *      Author: BethLocke
 */

#ifndef READLOCATION_H_
#define READLOCATION_H_

#include <iostream>

class ReadLocation {
public:
	ReadLocation();
	ReadLocation(std::string& id, bool left,bool forward,int start, int gap1,int gap2, int gap3);
	virtual ~ReadLocation();

	inline bool isForward() const { return _forward; }
	inline int getGap1() const {return _gap1; }
	inline int getGap2() const {return _gap2;}
	inline int getGap3() const {return _gap3;}
	inline int getStart() const {return _start;}
	inline bool isLeft() const {return _left;}
	inline const std::string& getViraId() const {return _viraID;}

private:
	//false is right
	bool _left;
	//false is reverse
	bool _forward;
	std::string _viraID;
	int _start;
	int _gap1;
	int _gap2;
	int _gap3;
};

#endif /* READLOCATION_H_ */
