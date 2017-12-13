/*
 * AssignedCondition2.h
 *
 *  Created on: Nov 4, 2015
 *      Author: e4k2
 */

#ifndef ASSIGNEDCONDITION2_H_
#define ASSIGNEDCONDITION2_H_

#include "TestCondition.h"
#include "CSProtein.h"

class AssignedCondition2: public TestCondition
{
private:
	int assignMatrix[MAXPROSIZE][MAXPROSIZE]; // indexed by res1-1, res2-1; elements that are out of bounds should have value 0
	                                          // assumes matrix is symmetric
	int countCutoffSR;
	int windowSizeSR; // inclusive; 0=only consider assignments between i and j; 1=also consider i+/-1 and j+/-1, ...
	int countCutoffLR;
	int windowSizeLR;

public:
	AssignedCondition2(int assignments[MAXPROSIZE][MAXPROSIZE], int countCutSR, int windowSR, int countCutLR, int windowLR);
	AssignedCondition2(const AssignedCondition2& ac);
	virtual ~AssignedCondition2();
	friend void swap(AssignedCondition2& ac1, AssignedCondition2& ac2);
	AssignedCondition2& operator=(AssignedCondition2 ac);
	virtual bool test(PCAssignment* pca); // returns false if count of assignments between res i, j +/- window < countCutoff
	virtual void print(int indent);
};

#endif /* ASSIGNEDCONDITION2_H_ */
