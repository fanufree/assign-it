/*
 * InStructureCutCondition.h
 *
 *  Created on: Sep 12, 2014
 *      Author: e4k2
 */

#ifndef INSTRUCTURECUTCONDITION_H_
#define INSTRUCTURECUTCONDITION_H_

#include "TestCondition.h"

class ContactMap;

/**
 * Once numViols > fractNumassignedCut*NumAssignments (or num assignment possibilities - which one to use is detected automatically),
 * then test() will always return false
 * Otherwise test() will return true
 * Each violation is the negation of (dist <= distCutoff in >= fractCutoff of templates)
 */
class InStructureCutCondition: public TestCondition
{
private:
	double fractCutoff;
	double fractNumAssignedCut; // cutoff for number assignments or assignment possibilities
	int numAssignments; // equal to number of assignments or assignment possibilities (input from user)
	ContactMap& contactMap; // from Main.cpp
	int numViols; // keeps track of current number of distance violations where each violation is
	               // negation of (dist <= distCutoff in >= fractCutoff of templates)
	              // distCutoff is assumed to have been used to generate contactMap

public:
	InStructureCutCondition(double fractCutoff, double fractNumAssignedCut,
			int numAssignments, ContactMap& contactMap);
	InStructureCutCondition(const InStructureCutCondition& isc);
	virtual ~InStructureCutCondition();
	friend void swap(InStructureCutCondition& isc1, InStructureCutCondition& isc2);
	InStructureCutCondition& operator=(InStructureCutCondition isc);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
	virtual void reset(); // resets numViols
};

#endif /* INSTRUCTURECUTCONDITION_H_ */
