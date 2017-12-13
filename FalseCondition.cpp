/*
 * FalseCondition.cpp
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#include "FalseCondition.h"
#include <cstdio>
#include <algorithm>

FalseCondition::FalseCondition()
{
}

FalseCondition::FalseCondition(const FalseCondition& ac)
{
}

FalseCondition::~FalseCondition()
{
}

void swap(FalseCondition& ac1, FalseCondition& ac2)
{
}

FalseCondition& FalseCondition::operator=(FalseCondition ac)
{
	swap(*this,ac);
	return *this;
}

bool FalseCondition::test(PCAssignment* pca)
{
	return false;
}

void FalseCondition::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("False\n");
}

