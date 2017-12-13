/*
 * NumAssignmentsCondition.h
 *
 *  Created on: Oct 16, 2015
 *      Author: e4k2
 */

#ifndef NUMASSIGNMENTSCONDITION_H_
#define NUMASSIGNMENTSCONDITION_H_

#include "TestCondition.h"

class NumAssignmentsCondition: public TestCondition
{
	bool isLR; // true if count LR assignments, else count SR assignments
	int numAssigned; // num assigned LR or SR depending on isLR
	double fract;  // scale factor on numAssignable
	int numAssignable; // test returns true if numAssigned >= fract*numAssignable; else returns false
	bool retValue; // the fixed return value of test

public:
	NumAssignmentsCondition(bool isLR, int numAssigned, double fract, int numAssignable);
	NumAssignmentsCondition(const NumAssignmentsCondition& nac);
	virtual ~NumAssignmentsCondition();
	friend void swap(NumAssignmentsCondition& nac1, NumAssignmentsCondition& nac2);
	NumAssignmentsCondition& operator=(NumAssignmentsCondition nac);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
};

#endif /* NUMASSIGNMENTSCONDITION_H_ */
