/*
 * ScoreCountCondition2.cpp
 *
 *  Created on: Oct 13, 2015
 *      Author: e4k2
 */

#include "ScoreCountCondition2.h"
#include "Score.h"
#include "PCAssignment.h"
#include "Main.h"
#include <algorithm>
#include <cstdio>

ScoreCountCondition2::ScoreCountCondition2(int termsArrSR[], int lengthSR, double boundsArrSR[], int countSR,
		int termsArrLR[], int lengthLR, double boundsArrLR[], int countLR)
{
	this->lengthSR = lengthSR;
	this->termsSR = new int[lengthSR];
	this->boundsSR = new double[lengthSR];
	std::fill(this->termsSR, this->termsSR+lengthSR, -1);
	std::copy(termsArrSR, termsArrSR+lengthSR, termsSR);
	std::copy(boundsArrSR, boundsArrSR+lengthSR, boundsSR);
	this->countSR = countSR;

	this->lengthLR = lengthLR;
	this->termsLR = new int[lengthLR];
	this->boundsLR = new double[lengthLR];
	std::fill(this->termsLR, this->termsLR+lengthLR, -1);
	std::copy(termsArrLR, termsArrLR+lengthLR, termsLR);
	std::copy(boundsArrLR, boundsArrLR+lengthLR, boundsLR);
	this->countLR = countLR;

}

ScoreCountCondition2::ScoreCountCondition2(const ScoreCountCondition2& st)
{
	this->lengthSR = st.lengthSR;
	this->termsSR = new int[lengthSR];
	this->boundsSR = new double[lengthSR];
	std::fill(this->termsSR, this->termsSR+lengthSR, -1);
	std::copy(st.termsSR, st.termsSR+lengthSR, termsSR);
	std::copy(st.boundsSR, st.boundsSR+lengthSR, boundsSR);
	this->countSR = st.countSR;

	this->lengthLR = st.lengthLR;
	this->termsLR = new int[lengthLR];
	this->boundsLR = new double[lengthLR];
	std::fill(this->termsLR, this->termsLR+lengthLR, -1);
	std::copy(st.termsLR, st.termsLR+lengthLR, termsLR);
	std::copy(st.boundsLR, st.boundsLR+lengthLR, boundsLR);
	this->countLR = st.countLR;
}

ScoreCountCondition2::~ScoreCountCondition2()
{
	delete[] termsSR;
	delete[] boundsSR;
	delete[] termsLR;
	delete[] boundsLR;
}

void swap(ScoreCountCondition2& st1, ScoreCountCondition2& st2)
{
	using namespace std;

	int* tempTerms1SR = st1.termsSR;
	double* tempBounds1SR = st1.boundsSR;

	std::swap(st1.lengthSR, st2.lengthSR);
	st1.termsSR = st2.termsSR;
	st2.termsSR = tempTerms1SR;
	st1.boundsSR = st2.boundsSR;
	st2.boundsSR = tempBounds1SR;
	std::swap(st1.countSR, st2.countSR);

	int* tempTerms1LR = st1.termsLR;
	double* tempBounds1LR = st1.boundsLR;

	std::swap(st1.lengthLR, st2.lengthLR);
	st1.termsLR = st2.termsLR;
	st2.termsLR = tempTerms1LR;
	st1.boundsLR = st2.boundsLR;
	st2.boundsLR = tempBounds1LR;
	std::swap(st1.countLR, st2.countLR);
}

ScoreCountCondition2& ScoreCountCondition2::operator=(ScoreCountCondition2 st)
{
	swap(*this, st);
	return *this;
}

bool ScoreCountCondition2::test(PCAssignment* pca)
{
	Score& s = pca->score;
	int c = 0;
	int seqsep = pca->pc.getSeqSep();
	if (seqsep < LRSEQSEP)
	{
		for (int i = 0; i < lengthSR; i++)
		{
			if (termsSR[i] > 10000)
			{
				unsigned int flag = termsSR[i]-10000;
				if (s.getScore(flag) >= boundsSR[i])
					c++;
				continue;
			}
			switch (termsSR[i])
			{
			case 0:
				if (s.cs >= boundsSR[i])
					c++;
				break;
			case 1:
				if (s.str >= boundsSR[i])
					c++;
				break;
			case 2:
				if (s.intensity >= boundsSR[i])
					c++;
				break;
			case 3:
				if (s.sym >= boundsSR[i])
					c++;
				break;
			case 4:
				if (s.interres >= boundsSR[i])
					c++;
				break;
			case 5:
				if (s.net >= boundsSR[i])
					c++;
				break;
			case 6:
				if (s.netStr >= boundsSR[i])
					c++;
				break;
			case 7:
				if (s.ambig >= boundsSR[i])
					c++;
				break;
			case 8:
				if (s.db >= boundsSR[i])
					c++;
				break;
			case 9:
				if (s.total >= boundsSR[i])
					c++;
				break;
			default:
				printf("Invalid terms\n");
				exit(-1);
			}
		}
		if (c >= countSR)
			return true;
		else
			return false;
	}
	else
	{
		for (int i = 0; i < lengthLR; i++)
		{
			if (termsLR[i] > 10000)
			{
				unsigned int flag = termsLR[i]-10000;
				if (s.getScore(flag) >= boundsLR[i])
					c++;
				continue;
			}
			switch (termsLR[i])
			{
			case 0:
				if (s.cs >= boundsLR[i])
					c++;
				break;
			case 1:
				if (s.str >= boundsLR[i])
					c++;
				break;
			case 2:
				if (s.intensity >= boundsLR[i])
					c++;
				break;
			case 3:
				if (s.sym >= boundsLR[i])
					c++;
				break;
			case 4:
				if (s.interres >= boundsLR[i])
					c++;
				break;
			case 5:
				if (s.net >= boundsLR[i])
					c++;
				break;
			case 6:
				if (s.netStr >= boundsLR[i])
					c++;
				break;
			case 7:
				if (s.ambig >= boundsLR[i])
					c++;
				break;
			case 8:
				if (s.db >= boundsLR[i])
					c++;
				break;
			case 9:
				if (s.total >= boundsLR[i])
					c++;
				break;
			default:
				printf("Invalid terms\n");
				exit(-1);
			}
		}
		if (c >= countLR)
			return true;
		else
			return false;
	}
	return true;
}

void ScoreCountCondition2::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("ScoreCountBounds SR, CountThreshold=%d, NumConditions=%d\n",countSR,lengthSR);
	for (int i = 0; i < lengthSR; i++)
	{
		if (termsSR[i] > 10000)
		{
			unsigned int flag = termsSR[i]-10000;
			for (int j = 0; j < indent; j++)
				printf("\t");
			printf("\tflag %u, %8.6f\n",flag, boundsSR[i]);
			continue;
		}
		else
		{
			switch (termsSR[i])
			{
			case 0:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tcs,%8.6f\n",boundsSR[i]);
				break;
			case 1:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tstr,%8.6f\n",boundsSR[i]);
				break;
			case 2:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tintensity,%8.6f\n",boundsSR[i]);
				break;
			case 3:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tsym,%8.6f\n",boundsSR[i]);
				break;
			case 4:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tinterres,%8.6f\n",boundsSR[i]);
				break;
			case 5:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tnet,%8.6f\n",boundsSR[i]);
				break;
			case 6:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tnetStr,%8.6f\n",boundsSR[i]);
				break;
			case 7:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tambig,%8.6f\n",boundsSR[i]);
				break;
			case 8:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tbias,%8.6f\n",boundsSR[i]);
				break;
			case 9:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\ttotal,%8.6f\n",boundsSR[i]);
				break;
			default:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tInvalid terms\n");
				exit(-1);
			}
		}
	}

	// LR
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("ScoreCountBounds LR, CountThreshold=%d, NumConditions=%d\n",countLR,lengthLR);
	for (int i = 0; i < lengthLR; i++)
	{
		if (termsLR[i] > 10000)
		{
			unsigned int flag = termsLR[i]-10000;
			for (int j = 0; j < indent; j++)
				printf("\t");
			printf("\tflag %u, %8.6f\n",flag, boundsLR[i]);
			continue;
		}
		else
		{
			switch (termsLR[i])
			{
			case 0:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tcs,%8.6f\n",boundsLR[i]);
				break;
			case 1:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tstr,%8.6f\n",boundsLR[i]);
				break;
			case 2:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tintensity,%8.6f\n",boundsLR[i]);
				break;
			case 3:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tsym,%8.6f\n",boundsLR[i]);
				break;
			case 4:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tinterres,%8.6f\n",boundsLR[i]);
				break;
			case 5:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tnet,%8.6f\n",boundsLR[i]);
				break;
			case 6:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tnetStr,%8.6f\n",boundsLR[i]);
				break;
			case 7:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tambig,%8.6f\n",boundsLR[i]);
				break;
			case 8:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tbias,%8.6f\n",boundsLR[i]);
				break;
			case 9:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\ttotal,%8.6f\n",boundsLR[i]);
				break;
			default:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tInvalid terms\n");
				exit(-1);
			}
		}
	}
}
