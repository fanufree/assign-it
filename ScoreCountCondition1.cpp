/*
 * ScoreCountCondition1.cpp
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#include "ScoreCountCondition1.h"

#include "Score.h"
#include "PCAssignment.h"
#include <algorithm>
#include <cstdio>


ScoreCountCondition1::ScoreCountCondition1(int termsArr[], int length, double boundsArr[], int count)
{
	this->length = length;
	this->terms = new int[length];
	this->bounds = new double[length];
	std::fill(this->terms, this->terms+length, -1);
	std::copy(termsArr, termsArr+length, terms);
	std::copy(boundsArr, boundsArr+length, bounds);
	this->count = count;
}

ScoreCountCondition1::ScoreCountCondition1(const ScoreCountCondition1& st)
{
	this->length = st.length;
	this->terms = new int[length];
	this->bounds = new double[length];
	std::fill(this->terms, this->terms+length, -1);
	std::copy(st.terms, st.terms+length, terms);
	std::copy(st.bounds, st.bounds+length, bounds);
	this->count = st.count;
}

ScoreCountCondition1::~ScoreCountCondition1()
{
	delete[] terms;
	delete[] bounds;
}

void swap(ScoreCountCondition1& st1, ScoreCountCondition1& st2)
{
	using namespace std;

	int* tempTerms1 = st1.terms;
	double* tempBounds1 = st1.bounds;

	std::swap(st1.length, st2.length);
	st1.terms = st2.terms;
	st2.terms = tempTerms1;
	st1.bounds = st2.bounds;
	st2.bounds = tempBounds1;
	std::swap(st1.count, st2.count);
}

ScoreCountCondition1& ScoreCountCondition1::operator=(ScoreCountCondition1 st)
{
	swap(*this, st);
	return *this;
}

bool ScoreCountCondition1::test(PCAssignment* pca)
{
	Score& s = pca->score;
	int c = 0;
	for (int i = 0; i < length; i++)
	{
		if (terms[i] > 10000)
		{
			unsigned int flag = terms[i]-10000;
			if (s.getScore(flag) >= bounds[i])
				c++;
			continue;
		}
		switch (terms[i])
		{
		case 0:
			if (s.cs >= bounds[i])
				c++;
			break;
		case 1:
			if (s.str >= bounds[i])
				c++;
			break;
		case 2:
			if (s.intensity >= bounds[i])
				c++;
			break;
		case 3:
			if (s.sym >= bounds[i])
				c++;
			break;
		case 4:
			if (s.interres >= bounds[i])
				c++;
			break;
		case 5:
			if (s.net >= bounds[i])
				c++;
			break;
		case 6:
			if (s.netStr >= bounds[i])
				c++;
			break;
		case 7:
			if (s.ambig >= bounds[i])
				c++;
			break;
		case 8:
			if (s.db >= bounds[i])
				c++;
			break;
		case 9:
			if (s.total >= bounds[i])
				c++;
			break;
		default:
			printf("Invalid terms\n");
			exit(-1);
		}
	}
	if (c >= count)
		return true;
	else
		return false;
}


void ScoreCountCondition1::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("ScoreCountBounds, CountThreshold=%d, NumConditions=%d\n",count,length);
	for (int i = 0; i < length; i++)
	{
		if (terms[i] > 10000)
		{
			unsigned int flag = terms[i]-10000;
			for (int j = 0; j < indent; j++)
				printf("\t");
			printf("\tflag %u, %8.6f\n",flag, bounds[i]);
			continue;
		}
		else
		{
			switch (terms[i])
			{
			case 0:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tcs,%8.6f\n",bounds[i]);
				break;
			case 1:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tstr,%8.6f\n",bounds[i]);
				break;
			case 2:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tintensity,%8.6f\n",bounds[i]);
				break;
			case 3:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tsym,%8.6f\n",bounds[i]);
				break;
			case 4:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tinterres,%8.6f\n",bounds[i]);
				break;
			case 5:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tnet,%8.6f\n",bounds[i]);
				break;
			case 6:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tnetStr,%8.6f\n",bounds[i]);
				break;
			case 7:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tambig,%8.6f\n",bounds[i]);
				break;
			case 8:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\tbias,%8.6f\n",bounds[i]);
				break;
			case 9:
				for (int j = 0; j < indent; j++)
					printf("\t");
				printf("\ttotal,%8.6f\n",bounds[i]);
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
