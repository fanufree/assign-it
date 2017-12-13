/*
 * AssignedCondition.cpp
 *
 *  Created on: Jun 20, 2014
 *      Author: e4k2
 */

#include "AssignedCondition.h"
#include "PCAssignment.h"
#include <cstdio>

AssignedCondition::AssignedCondition(int assignments[MAXPROSIZE][MAXPROSIZE], int countCut, int window)
{
	for (int i = 0; i < MAXPROSIZE; i++)
	{
		for (int j = 0; j < MAXPROSIZE; j++)
		{
			assignMatrix[i][j] = assignments[i][j];
		}
	}
	countCutoff = countCut;
	windowSize = window;
}

AssignedCondition::AssignedCondition(const AssignedCondition& ac)
{
	for (int i = 0; i < MAXPROSIZE; i++)
	{
		for (int j = 0; j < MAXPROSIZE; j++)
		{
			assignMatrix[i][j] = ac.assignMatrix[i][j];
		}
	}
	countCutoff = ac.countCutoff;
	windowSize = ac.windowSize;
}

AssignedCondition::~AssignedCondition()
{
}

void swap(AssignedCondition& ac1, AssignedCondition& ac2)
{
	using namespace std;
	swap(ac1.assignMatrix, ac2.assignMatrix);
	swap(ac1.countCutoff, ac2.countCutoff);
	swap(ac1.windowSize, ac2.windowSize);
}

AssignedCondition& AssignedCondition::operator =(AssignedCondition ac)
{
	swap(*this, ac);
	return *this;
}

bool AssignedCondition::test(PCAssignment* pca)
{
	int res1;
	int res2;
	pca->pc.getResPair(res1, res2);
	int count = 0;
	int i1 = res1-1;
	int i2 = res2-1;
	for (int w1 = -windowSize; w1 <= windowSize; w1++)
	{
		if (i1+w1 < 0)
			continue;
		if (i1+w1 >= MAXPROSIZE)
			continue;

		for (int w2 = -windowSize; w2 <= windowSize; w2++)
		{
			if (i2+w2 < 0)
				continue;
			if (i2+w2 >= MAXPROSIZE)
				continue;
			count += assignMatrix[i1+w1][i2+w2];
		}
	}
	if (count < countCutoff)
		return false;
	else
		return true;
}

void AssignedCondition::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("AssignedCondition: CountCutoff: %d  WindowSize: %d\n",countCutoff, windowSize);
}
