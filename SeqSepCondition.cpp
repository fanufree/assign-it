/*
 * SeqSepCondition.cpp
 *
 *  Created on: Jun 29, 2014
 *      Author: e4k2
 */

#include "SeqSepCondition.h"
#include "PseudoContact.h"
#include "PCAssignment.h"
#include <cstdio>
#include <algorithm>

SeqSepCondition::SeqSepCondition(int seqSep) : seqSepCutoff(seqSep)
{
}

SeqSepCondition::SeqSepCondition(const SeqSepCondition& ssc) : seqSepCutoff(ssc.seqSepCutoff)
{
}

SeqSepCondition::~SeqSepCondition()
{
}

void swap(SeqSepCondition& ssc1, SeqSepCondition& ssc2)
{
	using namespace std;
	swap(ssc1.seqSepCutoff, ssc2.seqSepCutoff);
}

SeqSepCondition& SeqSepCondition::operator=(SeqSepCondition ssc)
{
	swap(*this, ssc);
	return *this;
}

bool SeqSepCondition::test(PCAssignment* pca)
{
	PseudoContact& pc = pca->pc;
	int res1 = 0;
	int res2 = 0;
	pc.getResPair(res1,res2);
	int seqSep = abs(res1-res2);
	if (seqSep < seqSepCutoff)
		return true;
	else
		return false;
}

void SeqSepCondition::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("SeqSep: < %d\n",seqSepCutoff);
}
