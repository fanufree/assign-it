/*
 * ScoreTermCondition1.h
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#ifndef SCORETERMCONDITION1_H_
#define SCORETERMCONDITION1_H_

#include "TestCondition.h"

class PCAssignment;

class ScoreTermCondition1: public TestCondition
{
private:
	int* terms; // terms[i]: 0=cs, 1=str,2=intensity, 3=sym, 4=interres, 5=net, 6=netStr, 7=ambig, 8=db, 9=total, other_values=flag
	int length;
	double* bounds; // lower bound/minimum value on score; test(pca) returns true if terms >= bounds; false otherwise
	                 // use small positive number like 0.00001 for bounds if want terms to be > 0

public:
	ScoreTermCondition1(int termsArr[], int length, double boundsArr[]); // the array is copied
	ScoreTermCondition1(const ScoreTermCondition1& st);
	virtual ~ScoreTermCondition1();
	friend void swap(ScoreTermCondition1& st1, ScoreTermCondition1& st2);
	ScoreTermCondition1& operator=(ScoreTermCondition1 st);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
};

#endif /* SCORETERMCONDITION1_H_ */
