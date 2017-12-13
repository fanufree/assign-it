/*
 * IfElseCondition.h
 *
 *  Created on: Jun 18, 2014
 *      Author: e4k2
 */

#ifndef IfElseCondition_H_
#define IfElseCondition_H_

#include "TestCondition.h"

class IfElseCondition: public TestCondition
{
private:
	TestCondition* testCond; // cannot be NULL
	TestCondition* ifBody;  // can be NULL
	TestCondition* elseBody; // can be NULL


public:
	IfElseCondition(TestCondition* testCond, TestCondition* ifBody, TestCondition* elseBody);
	IfElseCondition(const IfElseCondition& ifc);
	virtual ~IfElseCondition();
	friend void swap(IfElseCondition& ifc1, IfElseCondition& ifc2);
	IfElseCondition& operator=(IfElseCondition ifc);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
};

#endif /* IfElseCondition_H_ */
