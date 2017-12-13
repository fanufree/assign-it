/*
 * ScoreTermCondition1.cpp
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#include "ScoreTermCondition1.h"
#include "Score.h"
#include "PCAssignment.h"
#include <algorithm>
#include <cstdio>


ScoreTermCondition1::ScoreTermCondition1(int termsArr[], int length, double boundsArr[])
{
	this->length = length;
	this->terms = new int[length];
	this->bounds = new double[length];
	std::fill(this->terms, this->terms+length, -1);
	std::copy(termsArr, termsArr+length, terms);
	std::copy(boundsArr, boundsArr+length, bounds);
}

ScoreTermCondition1::ScoreTermCondition1(const ScoreTermCondition1& st)
{
	this->length = st.length;
	this->terms = new int[length];
	this->bounds = new double[length];
	std::fill(this->terms, this->terms+length, -1);
	std::copy(st.terms, st.terms+length, terms);
	std::copy(st.bounds, st.bounds+length, bounds);
}

ScoreTermCondition1::~ScoreTermCondition1()
{
	delete[] terms;
	delete[] bounds;
}

void swap(ScoreTermCondition1& st1, ScoreTermCondition1& st2)
{
	using namespace std;

	int* tempTerms1 = st1.terms;
	double* tempBounds1 = st1.bounds;

	std::swap(st1.length, st2.length);
	st1.terms = st2.terms;
	st2.terms = tempTerms1;
	st1.bounds = st2.bounds;
	st2.bounds = tempBounds1;
}

ScoreTermCondition1& ScoreTermCondition1::operator=(ScoreTermCondition1 st)
{
	swap(*this, st);
	return *this;
}

bool ScoreTermCondition1::test(PCAssignment* pca)
{
	Score& s = pca->score;
	for (int i = 0; i < length; i++)
	{
		if (terms[i] > 10000)
		{
			unsigned int flag = terms[i]-10000;
			if (s.getScore(flag) < bounds[i])
				return false;
			continue;
		}
		else
		{
			switch (terms[i])
			{
			case 0:
				if (s.cs < bounds[i])
					return false;
				break;
			case 1:
				if (s.str < bounds[i])
					return false;
				break;
			case 2:
				if (s.intensity < bounds[i])
					return false;
				break;
			case 3:
				if (s.sym < bounds[i])
					return false;
				break;
			case 4:
				if (s.interres < bounds[i])
					return false;
				break;
			case 5:
				if (s.net < bounds[i])
					return false;
				break;
			case 6:
				if (s.netStr < bounds[i])
					return false;
				break;
			case 7:
				if (s.ambig < bounds[i])
					return false;
				break;
			case 8:
				if (s.db < bounds[i])
					return false;
				break;
			case 9:
				if (s.total < bounds[i])
					return false;
				break;
			default:
				printf("Invalid terms\n");
				exit(-1);
			}
		}
	}
	return true;
}


void ScoreTermCondition1::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("ScoreTermBounds\n");
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
