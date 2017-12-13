/*
 * IfElseCondition.cpp
 *
 *  Created on: Jun 18, 2014
 *      Author: e4k2
 */

#include "IfElseCondition.h"
#include <stddef.h>
#include <algorithm>
#include <cstdio>



IfElseCondition::IfElseCondition(TestCondition* testCond, TestCondition* ifBody, TestCondition* elseBody) : testCond(testCond), ifBody(ifBody), elseBody(elseBody)
{
}

IfElseCondition::IfElseCondition(const IfElseCondition& ifc) : testCond(ifc.testCond), ifBody(ifc.ifBody), elseBody(ifc.elseBody)
{
}

IfElseCondition::~IfElseCondition()
{
	delete testCond;
	if (ifBody != NULL)
		delete ifBody;
	if (elseBody != NULL)
		delete elseBody;
}

void swap(IfElseCondition& ifc1, IfElseCondition& ifc2)
{
	using namespace std;
	swap(ifc1.testCond, ifc2.testCond);
	swap(ifc1.ifBody, ifc2.ifBody);
	swap(ifc1.elseBody, ifc2.elseBody);
}

IfElseCondition& IfElseCondition::operator=(IfElseCondition ifc)
{
	swap(*this, ifc);
	return *this;
}

bool IfElseCondition::test(PCAssignment* pca)
{
	if (testCond->test(pca))
	{
		if (ifBody != NULL)
			return ifBody->test(pca);
		else
			return true;
	}
	else
	{
		if (elseBody != NULL)
			return elseBody->test(pca);
		else
			return false;
	}
	return false;
}

void IfElseCondition::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("IF\n");
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("TEST\n");
	testCond->print(indent+1);
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("BODY\n");
	if (ifBody != NULL)
		ifBody->print(indent+1);
	else
	{
		for (int i = 0; i < indent; i++)
			printf("\t");
		printf("\tEmpty\n");
	}
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("ELSE\n");
	if (elseBody != NULL)
		elseBody->print(indent+1);
	else
	{
		for (int i = 0; i < indent; i++)
			printf("\t");
		printf("\tEmpty\n");
	}
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("ENDIF\n");
}
