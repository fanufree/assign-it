/*
 * InStructureConditionFract.h
 *
 *  Created on: Jun 27, 2014
 *      Author: e4k2
 */

#ifndef INSTRUCTURECONDITIONFRACT_H_
#define INSTRUCTURECONDITIONFRACT_H_

#include "TestCondition.h"
#include "ContactMap.h"

/**
 * Similar to InStructureCondition but with different type of cutoff
 * Only Long Range contacts (seqsep >= 6) are considered (test always returns true for assignments with seqsep < 6)
 * test(PCA) assumes that the pca's arguments are passed in decreasing score order starting from the highest scoring pca
 * test(PCA) return false when
 * fraction(restraints/assignment_possibilities considered so far where min dist in templates > distCutoff) > fractCutoff
 * assuming currentCount > minCount. If currentCount <= minCount, test(PCA) always returns true
 * where currentCount is the current number of assignments tested so far
 * fraction = currentNumViols/currentCount
 * As soon as fraction > fractCutoff happens, test(PCA) will then always return false unless reset is called
 */
class InStructureConditionFract: public TestCondition
{
private:
	int minCount;
	int currentNumViols;
	int currentCount;
	double distCutoff;
	double fractCutoff;
	ContactMap& contactMap;
	int maxCount;
	int seqSep;
	int numTimes; // num times test(pca) is called since the last call to reset()

public:
	InStructureConditionFract(int minCount, double distCutoff, double fractCutoff, ContactMap& contactMap,
			int maxCount, int seqSep);
	InStructureConditionFract(const InStructureConditionFract& iscf);
	virtual ~InStructureConditionFract();
	friend void swap(InStructureConditionFract& iscf1, InStructureConditionFract& iscf2);
	InStructureConditionFract& operator=(InStructureConditionFract iscf);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
	virtual void reset(); // resets currentNumViols, currentCount
};

#endif /* INSTRUCTURECONDITIONFRACT_H_ */
