/*
 * OrCondition.cpp
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#include "OrCondition.h"
#include <cstdio>
#include <algorithm>


OrCondition::OrCondition(TestCondition* aa, TestCondition* bb) : a(aa), b(bb)
{
}

OrCondition::OrCondition(const OrCondition& ac) : a(ac.a), b(ac.b)
{
}

OrCondition::~OrCondition()
{
	delete a;
	delete b;
}

void swap(OrCondition& ac1, OrCondition& ac2)
{
	using namespace std;
	swap(ac1.a, ac2.a);
	swap(ac1.b, ac2.b);
}

OrCondition& OrCondition::operator=(OrCondition ac)
{
	swap(*this,ac);
	return *this;
}

bool OrCondition::test(PCAssignment* pca)
{
	return a->test(pca) || b->test(pca);
}

void OrCondition::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("OR\n");
	a->print(indent+1);
	b->print(indent+1);
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("End OR\n");
}
