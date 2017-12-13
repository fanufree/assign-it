/*
 * AssignedCondition2.cpp
 *
 *  Created on: Nov 4, 2015
 *      Author: e4k2
 */

#include "AssignedCondition2.h"
#include "PCAssignment.h"
#include "Main.h"
#include <cstdio>
#include <cmath>

AssignedCondition2::AssignedCondition2(int assignments[MAXPROSIZE][MAXPROSIZE], int countCutSR, int windowSR, int countCutLR, int windowLR)
{
	for (int i = 0; i < MAXPROSIZE; i++)
	{
		for (int j = 0; j < MAXPROSIZE; j++)
		{
			assignMatrix[i][j] = assignments[i][j];
		}
	}
	countCutoffSR = countCutSR;
	windowSizeSR = windowSR;

	countCutoffLR = countCutLR;
	windowSizeLR = windowLR;
}

AssignedCondition2::AssignedCondition2(const AssignedCondition2& ac)
{
	for (int i = 0; i < MAXPROSIZE; i++)
	{
		for (int j = 0; j < MAXPROSIZE; j++)
		{
			assignMatrix[i][j] = ac.assignMatrix[i][j];
		}
	}
	countCutoffSR = ac.countCutoffSR;
	windowSizeSR = ac.windowSizeSR;
	countCutoffLR = ac.countCutoffLR;
	windowSizeLR = ac.windowSizeLR;
}

AssignedCondition2::~AssignedCondition2()
{
}

void swap(AssignedCondition2& ac1, AssignedCondition2& ac2)
{
	using namespace std;
	swap(ac1.assignMatrix, ac2.assignMatrix);
	swap(ac1.countCutoffSR, ac2.countCutoffSR);
	swap(ac1.windowSizeSR, ac2.windowSizeSR);

	swap(ac1.countCutoffLR, ac2.countCutoffLR);
	swap(ac1.windowSizeLR, ac2.windowSizeLR);
}

AssignedCondition2& AssignedCondition2::operator =(AssignedCondition2 ac)
{
	swap(*this, ac);
	return *this;
}

bool AssignedCondition2::test(PCAssignment* pca)
{
	int res1;
	int res2;
	pca->pc.getResPair(res1, res2);
	int count = 0;
	int i1 = res1-1;
	int i2 = res2-1;
	if (abs(res2-res1) < LRSEQSEP)
	{
		for (int w1 = -windowSizeSR; w1 <= windowSizeSR; w1++)
		{
			if (i1+w1 < 0)
				continue;
			if (i1+w1 >= MAXPROSIZE)
				continue;

			for (int w2 = -windowSizeSR; w2 <= windowSizeSR; w2++)
			{
				if (i2+w2 < 0)
					continue;
				if (i2+w2 >= MAXPROSIZE)
					continue;
				count += assignMatrix[i1+w1][i2+w2];
			}
		}
		if (count < countCutoffSR)
			return false;
		else
			return true;
	}
	else
	{
		for (int w1 = -windowSizeLR; w1 <= windowSizeLR; w1++)
		{
			if (i1+w1 < 0)
				continue;
			if (i1+w1 >= MAXPROSIZE)
				continue;

			for (int w2 = -windowSizeLR; w2 <= windowSizeLR; w2++)
			{
				if (i2+w2 < 0)
					continue;
				if (i2+w2 >= MAXPROSIZE)
					continue;
				count += assignMatrix[i1+w1][i2+w2];
			}
		}
		if (count < countCutoffLR)
			return false;
		else
			return true;

	}
}

void AssignedCondition2::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("AssignedCondition2: CountCutoffSR: %d  WindowSizeSR: %d  CountCutoffLR: %d  WindowSizeLR: %d\n",
			countCutoffSR, windowSizeSR, countCutoffLR, windowSizeLR);
}

