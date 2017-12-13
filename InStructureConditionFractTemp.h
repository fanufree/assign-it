/*
 * InStructureConditionFractTemp.h
 *
 *  Created on: Jun 29, 2014
 *      Author: e4k2
 */

#ifndef INSTRUCTURECONDITIONFRACTTEMP_H_
#define INSTRUCTURECONDITIONFRACTTEMP_H_

#include "TestCondition.h"
#include <list>

class CSProtein;
class ContactMap;

using namespace std;

/**
 * Does not require assigments to be sorted by score
 * test(pca) returns false when
 * fraction of templates satisifying restraint with dist < distCutoff is < fractTempCutoff (true if >= fractTempCutoff)
 * Only contacts with seqsep > 1 are considered; test(pca) always returns true for contacts with seqsep <= 1
 */
class InStructureConditionFractTemp: public TestCondition
{
private:
	double distCutoff;
	double fractTempCutoff;
	ContactMap* contactMap;

public:
	// structures must contain the 3D coordinates of the models
	InStructureConditionFractTemp(double distCutoff, double fractTempCutoff, list<CSProtein*>& structures);
	InStructureConditionFractTemp(const InStructureConditionFractTemp& is);
	virtual ~InStructureConditionFractTemp();
	friend void swap(InStructureConditionFractTemp& is1, InStructureConditionFractTemp& is2);
	InStructureConditionFractTemp& operator=(InStructureConditionFractTemp is);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
};

#endif /* INSTRUCTURECONDITIONFRACTTEMP_H_ */
