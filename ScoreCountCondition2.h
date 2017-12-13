/*
 * ScoreCountCondition2.h
 *
 *  Created on: Oct 13, 2015
 *      Author: e4k2
 */

#ifndef SCORECOUNTCONDITION2_H_
#define SCORECOUNTCONDITION2_H_

#include "TestCondition.h"

class ScoreCountCondition2: public TestCondition
{
private:
	int* termsSR;
	int* termsLR;
	int lengthSR;
	int lengthLR;
	double* boundsSR;
	double* boundsLR;
	int countSR;
	int countLR;

public:
	ScoreCountCondition2(int termsArrSR[], int lengthSR, double boundsArrSR[], int countSR,
			int termsArrLR[], int lengthLR, double boundsArrLR[], int countLR);
	ScoreCountCondition2(const ScoreCountCondition2& st);
	virtual ~ScoreCountCondition2();
	friend void swap(ScoreCountCondition2& st1, ScoreCountCondition2& st2);
	ScoreCountCondition2& operator=(ScoreCountCondition2 st);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
};

#endif /* SCORECOUNTCONDITION2_H_ */
