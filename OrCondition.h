/*
 * OrCondition.h
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#ifndef ORCONDITION_H_
#define ORCONDITION_H_

#include "TestCondition.h"

/**
 * Note in test(pca), condition b is not evaluated if condition a already evaluates to true
 * Therefore, for InStructureCondition and related conditions, if you want to always evaluate the condition
 * (so that you can record some "state" of the condition like currentCount in InStructureConditoin),
 * then put that condition as condition a instead of b
 */
class OrCondition: public TestCondition
{
private:
	TestCondition* a; // must be non-null
	TestCondition* b; // must be non-null
public:
	OrCondition(TestCondition* aa, TestCondition* bb);
	OrCondition(const OrCondition& ac);
	virtual ~OrCondition();
	friend void swap(OrCondition& ac1, OrCondition& ac2);
	OrCondition& operator=(OrCondition ac);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
};

#endif /* ORCONDITION_H_ */
