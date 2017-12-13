/*
 * NotCondition.cpp
 *
 *  Created on: Sep 12, 2014
 *      Author: e4k2
 */

#include "NotCondition.h"
#include <stddef.h>
#include <algorithm>
#include <cstdio>

NotCondition::NotCondition(TestCondition* c) : cond(c)
{
}

NotCondition::NotCondition(const NotCondition& nc) : cond(nc.cond)
{
}

NotCondition::~NotCondition()
{
}

void swap(NotCondition& nc1, NotCondition& nc2)
{
	using namespace std;
	swap(nc1.cond, nc2.cond);
}

NotCondition& NotCondition::operator =(NotCondition nc)
{
	swap(*this, nc);
	return *this;
}

bool NotCondition::test(PCAssignment* pca)
{
	return !(cond->test(pca));
}

void NotCondition::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("Not\n");
	cond->print(indent+1);
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("ENDNOT\n");
}
