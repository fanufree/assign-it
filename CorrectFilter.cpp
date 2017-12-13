/*
 * CorrectFilter.cpp
 *
 *  Created on: Nov 18, 2014
 *      Author: e4k2
 */

#include "CorrectFilter.h"
#include "ContactMap.h"
#include "PseudoContact.h"
#include "PCAssignment.h"
#include <cstdio>
#include <algorithm>
#include <tr1/array>
#include <tr1/unordered_set>

CorrectFilter::CorrectFilter(ContactMap& refMap, double minDistCutoff, double fractCutoff) :
                             refContactMap(refMap), minDistCutoff(minDistCutoff), fractCutoff(fractCutoff)
{
}

CorrectFilter::CorrectFilter(const CorrectFilter& cf) : refContactMap(cf.refContactMap), minDistCutoff(cf.minDistCutoff), fractCutoff(cf.fractCutoff)
{
}

CorrectFilter::~CorrectFilter()
{
}

void swap(CorrectFilter& cf1, CorrectFilter& cf2)
{
	using namespace std;
	swap(cf1.refContactMap, cf2.refContactMap);
	swap(cf1.minDistCutoff, cf2.minDistCutoff);
	swap(cf1.fractCutoff, cf2.fractCutoff);
}

CorrectFilter& CorrectFilter::operator=(CorrectFilter cf)
{
	swap(*this, cf);
	return *this;
}

bool CorrectFilter::test(PCAssignment* pca)
{
	tr1::unordered_set<Contact>& contacts = pca->pc.getContacts();
	for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
	{
		const Contact& c = *itC;
		tr1::array<double,ContactMap::CMENTRYSIZE>& entry = refContactMap[c];
		double minDist = entry[ContactMap::MINDISTINDEX];
		double fractStr = entry[ContactMap::FRACSTRUCINDEX];
		if (minDist <= minDistCutoff && fractStr >= fractCutoff)
			return true;
	}
	return false;
}

void CorrectFilter::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("CorrectFilter: minDistCutoff=%6.3f fractCutoff=%4.2f\n",minDistCutoff, fractCutoff);
}

void CorrectFilter::reset()
{
}
