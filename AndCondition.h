/*
 * AndCondition.h
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#ifndef ANDCONDITION_H_
#define ANDCONDITION_H_

#include "TestCondition.h"

class PCAssignment;

/**
 * Note that condition b is not evaluated in test(pca) if condition a evaluates to false
 * Therefore, for InStructureCondition and related conditions, if you want to always evaluate the condition
 * (so that you can record some "state" of the condition like currentCount in InStructureConditoin),
 * then put that condition as condition a instead of b
 */
class AndCondition : public TestCondition
{
private:
	TestCondition* a; // must be non-null
	TestCondition* b; // must be non-null
public:
	AndCondition(TestCondition* aa, TestCondition* bb);
	AndCondition(const AndCondition& ac);
	virtual ~AndCondition();
	friend void swap(AndCondition& ac1, AndCondition& ac2);
	AndCondition& operator=(AndCondition ac);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
};

#endif /* ANDCONDITION_H_ */
