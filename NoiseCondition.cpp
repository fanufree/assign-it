/*
 * NoiseCondition.cpp
 *
 *  Created on: Jul 24, 2015
 *      Author: e4k2
 */

#include "NoiseCondition.h"
#include <cstdio>
#include <algorithm>

NoiseCondition::NoiseCondition(int numNOEPeaks, int numResidues, double ratioTolerance)
{
	double r = double(numNOEPeaks)/double(numResidues);
	if (r <= ratioTolerance)
		testValue = true;
	else
		testValue = false;
}

NoiseCondition::NoiseCondition(const NoiseCondition& nc) : testValue(nc.testValue)
{
}

NoiseCondition::~NoiseCondition()
{
}

void swap(NoiseCondition& nc1, NoiseCondition& nc2)
{
	using std::swap;
	swap(nc1.testValue, nc2.testValue);
}

NoiseCondition& NoiseCondition::operator=(NoiseCondition nc)
{
	swap(*this, nc);
	return *this;
}

bool NoiseCondition::test(PCAssignment* pca)
{
	return testValue;
}

void NoiseCondition::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	if (testValue)
		printf("NOISE CONDITION true\n");
	else
		printf("NOISE CONDITION false\n");
}
