/*
 * ScoreCountCondition1.h
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#ifndef SCORECOUNTCONDITION1_H_
#define SCORECOUNTCONDITION1_H_

#include "TestCondition.h"

/**
 * test(pca) returns true if at least COUNT (inclusive) terms satisfy their bounds (>= bounds)
 */
class ScoreCountCondition1: public TestCondition
{
private:
	int* terms; // see ScoreTermCondition
	int length;
	double* bounds;
	int count;

public:
	ScoreCountCondition1(int termsArr[], int length, double boundsArr[], int count); // the array is copied
	ScoreCountCondition1(const ScoreCountCondition1& st);
	virtual ~ScoreCountCondition1();
	friend void swap(ScoreCountCondition1& st1, ScoreCountCondition1& st2);
	ScoreCountCondition1& operator=(ScoreCountCondition1 st);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
};

#endif /* SCORECOUNTCONDITION1_H_ */
