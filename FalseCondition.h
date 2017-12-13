/*
 * FalseCondition.h
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#ifndef FALSECONDITION_H_
#define FALSECONDITION_H_

#include "TestCondition.h"

class FalseCondition: public TestCondition
{
public:
	FalseCondition();
	FalseCondition(const FalseCondition& ac);
	virtual ~FalseCondition();
	friend void swap(FalseCondition& ac1, FalseCondition& ac2);
	FalseCondition& operator=(FalseCondition ac);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
};

#endif /* FALSECONDITION_H_ */
