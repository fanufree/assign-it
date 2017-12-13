/*
 * AndCondition.cpp
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#include "AndCondition.h"
#include <cstdio>
#include <algorithm>

AndCondition::AndCondition(TestCondition* aa, TestCondition* bb) : a(aa), b(bb)
{
}

AndCondition::AndCondition(const AndCondition& ac) : a(ac.a), b(ac.b)
{
}

AndCondition::~AndCondition()
{
	delete a;
	delete b;
}

void swap(AndCondition& ac1, AndCondition& ac2)
{
	using namespace std;
	std::swap(ac1.a, ac2.a);
	std::swap(ac1.b, ac2.b);
}

AndCondition& AndCondition::operator=(AndCondition ac)
{
	swap(*this,ac);
	return *this;
}

bool AndCondition::test(PCAssignment* pca)
{
	return a->test(pca) && b->test(pca);
}

void AndCondition::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("AND\n");
	a->print(indent+1);
	b->print(indent+1);
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("End And\n");
}
