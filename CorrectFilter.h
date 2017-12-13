/*
 * CorrectFilter.h
 *
 *  Created on: Nov 18, 2014
 *      Author: e4k2
 */

#ifndef CORRECTFILTER_H_
#define CORRECTFILTER_H_

#include "TestCondition.h"
#include "Contact.h"

class ContactMap;

/*
 * Used for testing the filters. This is the reference filter which uses the ContactMap from the native structure,
 * so it should give perfect filtering
 */
class CorrectFilter : public TestCondition
{
private:
	ContactMap& refContactMap; // from Main.cpp
	double minDistCutoff; // set to DISTCUTOFF in main.cpp
	double fractCutoff;

public:
	CorrectFilter(ContactMap& refMap, double minDistCutoff, double fractCutoff);
	CorrectFilter(const CorrectFilter& cf);
	virtual ~CorrectFilter();
	friend void swap(CorrectFilter& cf1, CorrectFilter& cf2);
	CorrectFilter& operator=(CorrectFilter cf);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
	virtual void reset();
};

#endif /* CORRECTFILTER_H_ */
