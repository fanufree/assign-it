/*
 * ScoreTermCondition2.cpp
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#include "ScoreTermCondition2.h"
#include "Score.h"
#include "PCAssignment.h"
#include "Main.h"
#include <algorithm>
#include <cstdio>


ScoreTermCondition2::ScoreTermCondition2(int termsArrSR[], int lengthSR, double boundsArrSR[], int termsArrLR[], int lengthLR, double boundsArrLR[])
{
	this->lengthSR = lengthSR;
	this->termsSR = new int[lengthSR];
	this->boundsSR = new double[lengthSR];
	std::fill(this->termsSR, this->termsSR+lengthSR, -1);
	std::copy(termsArrSR, termsArrSR+lengthSR, termsSR);
	std::copy(boundsArrSR, boundsArrSR+lengthSR, boundsSR);

	this->lengthLR = lengthLR;
	this->termsLR = new int[lengthLR];
	this->boundsLR = new double[lengthLR];
	std::fill(this->termsLR, this->termsLR+lengthLR, -1);
	std::copy(termsArrLR, termsArrLR+lengthLR, termsLR);
	std::copy(boundsArrLR, boundsArrLR+lengthLR, boundsLR);
}

ScoreTermCondition2::ScoreTermCondition2(const ScoreTermCondition2& st)
{
	this->lengthSR = st.lengthSR;
	this->termsSR = new int[lengthSR];
	this->boundsSR = new double[lengthSR];
	std::fill(this->termsSR, this->termsSR+lengthSR, -1);
	std::copy(st.termsSR, st.termsSR+lengthSR, termsSR);
	std::copy(st.boundsSR, st.boundsSR+lengthSR, boundsSR);

	this->lengthLR = st.lengthLR;
	this->termsLR = new int[lengthLR];
	this->boundsLR = new double[lengthLR];
	std::fill(this->termsLR, this->termsLR+lengthLR, -1);
	std::copy(st.termsLR, st.termsLR+lengthLR, termsLR);
	std::copy(st.boundsLR, st.boundsLR+lengthLR, boundsLR);
}

ScoreTermCondition2::~ScoreTermCondition2()
{
	delete[] termsSR;
	delete[] boundsSR;
	delete[] termsLR;
	delete[] boundsLR;
}

void swap(ScoreTermCondition2& st1, ScoreTermCondition2& st2)
{
	using namespace std;

	int* tempTerms1SR = st1.termsSR;
	double* tempBounds1SR = st1.boundsSR;

	std::swap(st1.lengthSR, st2.lengthSR);
	st1.termsSR = st2.termsSR;
	st2.termsSR = tempTerms1SR;
	st1.boundsSR = st2.boundsSR;
	st2.boundsSR = tempBounds1SR;

	int* tempTerms1LR = st1.termsLR;
	double* tempBounds1LR = st1.boundsLR;

	std::swap(st1.lengthLR, st2.lengthLR);
	st1.termsLR = st2.termsLR;
	st2.termsLR = tempTerms1LR;
	st1.boundsLR = st2.boundsLR;
	st2.boundsLR = tempBounds1LR;
}

ScoreTermCondition2& ScoreTermCondition2::operator=(ScoreTermCondition2 st)
{
	swap(*this, st);
	return *this;
}

bool ScoreTermCondition2::test(PCAssignment* pca)
{
	Score& s = pca->score;
	int seqsep = pca->pc.getSeqSep();
	if (seqsep < LRSEQSEP)
	{
		for (int i = 0; i < lengthSR; i++)
		{
			if (termsSR[i] > 10000)
			{
				unsigned int flag = termsSR[i]-10000;
				if (s.getScore(flag) < boundsSR[i])
					return false;
				continue;
			}
			else
			{
				switch (termsSR[i])
				{
				case 0:
					if (s.cs < boundsSR[i])
						return false;
					break;
				case 1:
					if (s.str < boundsSR[i])
						return false;
					break;
				case 2:
					if (s.intensity < boundsSR[i])
						return false;
					break;
				case 3:
					if (s.sym < boundsSR[i])
						return false;
					break;
				case 4:
					if (s.interres < boundsSR[i])
						return false;
					break;
				case 5:
					if (s.net < boundsSR[i])
						return false;
					break;
				case 6:
					if (s.netStr < boundsSR[i])
						return false;
					break;
				case 7:
					if (s.ambig < boundsSR[i])
						return false;
					break;
				case 8:
					if (s.db < boundsSR[i])
						return false;
					break;
				case 9:
					if (s.total < boundsSR[i])
						return false;
					break;
				default:
					printf("Invalid terms\n");
					exit(-1);
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < lengthLR; i++)
		{
			if (termsLR[i] > 10000)
			{
				unsigned int flag = termsLR[i]-10000;
				if (s.getScore(flag) < boundsLR[i])
					return false;
				continue;
			}
			else
			{
				switch (termsLR[i])
				{
				case 0:
					if (s.cs < boundsLR[i])
						return false;
					break;
				case 1:
					if (s.str < boundsLR[i])
						return false;
					break;
				case 2:
					if (s.intensity < boundsLR[i])
						return false;
					break;
				case 3:
					if (s.sym < boundsLR[i])
						return false;
					break;
				case 4:
					if (s.interres < boundsLR[i])
						return false;
					break;
				case 5:
					if (s.net < boundsLR[i])
						return false;
					break;
				case 6:
					if (s.netStr < boundsLR[i])
						return false;
					break;
				case 7:
					if (s.ambig < boundsLR[i])
						return false;
					break;
				case 8:
					if (s.db < boundsLR[i])
						return false;
					break;
				case 9:
					if (s.total < boundsLR[i])
						return false;
					break;
				default:
					printf("Invalid terms\n");
					exit(-1);
				}
			}
		}
	}
	return true;
}


void ScoreTermCondition2::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("ScoreTermBounds SR\n");
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

	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("ScoreTermBounds LR\n");
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
