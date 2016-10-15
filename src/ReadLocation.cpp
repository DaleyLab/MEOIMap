/*
 * ReadLocation.cpp
 *
 *  Created on: Aug 30, 2016
 *      Author: BethLocke
 */

#include "ReadLocation.h"

ReadLocation::ReadLocation() {
	// TODO Auto-generated constructor stub

}

ReadLocation::ReadLocation(std::string& id, bool left, bool forward, int start, int gap1,
		int gap2, int gap3) {
	this->forward = forward;
	this->left = left;
	this->viraID = id;
	this->start = start;
	this->gap1 = gap1;
	this->gap2 = gap2;
	this->gap3 = gap3;
}

ReadLocation::~ReadLocation() {
	// TODO Auto-generated destructor stub
}

