/*
 * TrueCondition.cpp
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#include "TrueCondition.h"
#include <cstdio>
#include <algorithm>


TrueCondition::TrueCondition()
{
}

TrueCondition::TrueCondition(const TrueCondition& ac)
{
}

TrueCondition::~TrueCondition()
{
}

void swap(TrueCondition& ac1, TrueCondition& ac2)
{
}

TrueCondition& TrueCondition::operator=(TrueCondition ac)
{
	swap(*this,ac);
	return *this;
}

bool TrueCondition::test(PCAssignment* pca)
{
	return true;
}

void TrueCondition::print(int indent)
{
	for (int i = 0; i < indent; i++) {
		printf("\t");
	}
	printf("True\n");
}

