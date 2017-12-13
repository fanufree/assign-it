/*
 * ScoreTermCondition2.h
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#ifndef SCORETERMCONDITION2_H_
#define SCORETERMCONDITION2_H_

#include "TestCondition.h"

class PCAssignment;

class ScoreTermCondition2: public TestCondition
{
private:
	int* termsSR; // terms[i]: 0=cs, 1=str,2=intensity, 3=sym, 4=interres, 5=net, 6=netStr, 7=ambig, 8=db, 9=total, other_values=flag
	int* termsLR;
	int lengthSR;
	int lengthLR;
	double* boundsSR; // lower bound/minimum value on score; test(pca) returns true if terms >= bounds; false otherwise
	                 // use small positive number like 0.00001 for bounds if want terms to be > 0
	double* boundsLR;

public:
	ScoreTermCondition2(int termsArrSR[], int lengthSR, double boundsArrSR[], int termsArrLR[], int lengthLR, double boundsArrLR[]); // the array is copied
	ScoreTermCondition2(const ScoreTermCondition2& st);
	virtual ~ScoreTermCondition2();
	friend void swap(ScoreTermCondition2& st1, ScoreTermCondition2& st2);
	ScoreTermCondition2& operator=(ScoreTermCondition2 st);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
};

#endif /* SCORETERMCONDITION2_H_ */
