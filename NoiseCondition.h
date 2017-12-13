/*
 * NoiseCondition.h
 *
 *  Created on: Jul 24, 2015
 *      Author: e4k2
 */

#ifndef NOISECONDITION_H_
#define NOISECONDITION_H_

#include "TestCondition.h"

class NoiseCondition : public TestCondition
{
public:
	// test(pca) returns true if numNOEPeaks/numResidues <= ratioTolerance
	NoiseCondition(int numNOEPeaks, int numResidues, double ratioTolerance);
	NoiseCondition(const NoiseCondition& nc);
	virtual ~NoiseCondition();
	friend void swap(NoiseCondition& nc1, NoiseCondition& nc2);
	NoiseCondition& operator=(NoiseCondition nc);
	virtual bool test(PCAssignment* pca); // returns true if numNOEPeaks/numResidues <= ratioTolerance
	virtual void print(int indent);
	bool testValue; // returned by test
};

#endif /* NOISECONDITION_H_ */
