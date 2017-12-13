/*
 * AssignedCondition.h
 *
 *  Created on: Jun 20, 2014
 *      Author: e4k2
 */

#ifndef ASSIGNEDCONDITION_H_
#define ASSIGNEDCONDITION_H_

#include "TestCondition.h"
#include "CSProtein.h"

/**
 * For checking for isolated assignment contacts
 */
class AssignedCondition: public TestCondition
{
private:
	int assignMatrix[MAXPROSIZE][MAXPROSIZE]; // indexed by res1-1, res2-1; elements that are out of bounds should have value 0
	                                          // assumes matrix is symmetric
	int countCutoff;
	int windowSize; // inclusive; 0=only consider assignments between i and j; 1=also consider i+/-1 and j+/-1, ...


public:
	AssignedCondition(int assignments[MAXPROSIZE][MAXPROSIZE], int countCut, int window);
	AssignedCondition(const AssignedCondition& ac);
	virtual ~AssignedCondition();
	friend void swap(AssignedCondition& ac1, AssignedCondition& ac2);
	AssignedCondition& operator=(AssignedCondition ac);
	virtual bool test(PCAssignment* pca); // returns false if count of assignments between res i, j +/- window < countCutoff
	virtual void print(int indent);
};

#endif /* ASSIGNEDCONDITION_H_ */
