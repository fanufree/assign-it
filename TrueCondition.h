/*
 * TrueCondition.h
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#ifndef TRUECONDITION_H_
#define TRUECONDITION_H_

#include "TestCondition.h"

class TrueCondition: public TestCondition
{
public:
	TrueCondition();
	TrueCondition(const TrueCondition& ac);
	virtual ~TrueCondition();
	friend void swap(TrueCondition& ac1, TrueCondition& ac2);
	TrueCondition& operator=(TrueCondition ac);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
};

#endif /* TRUECONDITION_H_ */
