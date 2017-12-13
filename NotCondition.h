/*
 * NotCondition.h
 *
 *  Created on: Sep 12, 2014
 *      Author: e4k2
 */

#ifndef NOTCONDITION_H_
#define NOTCONDITION_H_

#include "TestCondition.h"

class NotCondition: public TestCondition
{
private:
	TestCondition* cond; // cannot be null

public:
	NotCondition(TestCondition* c);
	NotCondition(const NotCondition& nc);
	virtual ~NotCondition();
	friend void swap(NotCondition& nc1, NotCondition& nc2);
	NotCondition& operator=(NotCondition nc);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
};

#endif /* NOTCONDITION_H_ */
