/*
 * NumAssignmentsCondition.cpp
 *
 *  Created on: Oct 16, 2015
 *      Author: e4k2
 */

#include "NumAssignmentsCondition.h"
#include "PCAssignment.h"
#include <cmath>
#include <algorithm>
#include <cstdio>

NumAssignmentsCondition::NumAssignmentsCondition(bool isLR, int numAssigned, double fract, int numAssignable) :
              isLR(isLR), numAssigned(numAssigned), fract(fract), numAssignable(numAssignable)
{
	int upperbound = (int)round(fract*numAssignable);
	if (numAssigned >= upperbound)
		retValue = true;
	else
		retValue = false;
}

NumAssignmentsCondition::NumAssignmentsCondition(const NumAssignmentsCondition& nac) : isLR(nac.isLR),
		numAssigned(nac.numAssigned), fract(nac.fract), numAssignable(nac.numAssignable), retValue(nac.retValue)
{
}

NumAssignmentsCondition::~NumAssignmentsCondition()
{
}

void swap(NumAssignmentsCondition& nac1, NumAssignmentsCondition& nac2)
{
	using namespace std;
	std::swap(nac1.isLR, nac2.isLR);
	std::swap(nac1.numAssigned, nac2.numAssigned);
	std::swap(nac1.numAssignable, nac2.numAssignable);
	std::swap(nac1.retValue, nac2.retValue);
}

NumAssignmentsCondition& NumAssignmentsCondition::operator=(NumAssignmentsCondition nac)
{
	swap(*this, nac);
	return *this;
}

bool NumAssignmentsCondition::test(PCAssignment* pca)
{
	return retValue;
}

void NumAssignmentsCondition::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	int upperbound = (int)round(fract*numAssignable);
	if (isLR)
	{
		if (retValue)
			printf("NumAssignments: true isLR %d >= %d (%f * %d)\n",
					numAssigned,upperbound,fract,numAssignable);
		else
			printf("NumAssignments: false isLR %d >= %d  (%f * %d)\n",
					numAssigned,upperbound,fract,numAssignable);
	}
	else
	{
		if (retValue)
			printf("NumAssignments: true isSR %d >= %d (%f * %d)\n",
					numAssigned,upperbound,fract,numAssignable);
		else
			printf("NumAssignments: false isSR %d >= %d  (%f * %d)\n",
					numAssigned,upperbound,fract,numAssignable);
	}
}
