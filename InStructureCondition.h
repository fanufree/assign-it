/*
 * InStructureCondition.h
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#ifndef INSTRUCTURECONDITION_H_
#define INSTRUCTURECONDITION_H_

#include "TestCondition.h"
#include "Contact.h"
#include <tr1/unordered_map>
#include <tr1/array>

class ContactMap;

/**
 * Keeps track of the current state of looping through all the assignments/assignment possibilities
 * test(pca) starts returning false when (based on constructor)
 * count(restraints/assignment_possibilities where min dist in templates > distCutoff) > countCutoff;
 * Only Long Range contacts (seqsep >= 6) are considered
 * Assignments should be sorted in decreasing order by score before calling test repeatedly
 */
class InStructureCondition: public TestCondition
{
private:
	double distCutoff;
	int countCutoff;
	int currentCount; // keeps track of current number of distance violations
	                  // count(restraints where min dist in templates > distCutoff)
	ContactMap& contactMap; // from Main.cpp

public:
	/**
	 * test long range contacts only (seqsep>=6) and keep rest of contacts (seqsep < 6)
	 */
	InStructureCondition(double distCutoff, int countCutoff, ContactMap& contactMap);
	InStructureCondition(const InStructureCondition& isc);
	virtual ~InStructureCondition();
	friend void swap(InStructureCondition& isc1, InStructureCondition& isc2);
	InStructureCondition& operator=(InStructureCondition isc);
	virtual bool test(PCAssignment* pca); // returns false if condition (a) is reached
	virtual void print(int indent);
	virtual void reset(); // resets the currentCount
};

#endif /* INSTRUCTURECONDITION_H_ */
