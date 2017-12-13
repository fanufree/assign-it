/*
 * InStructureCondition.cpp
 *
 *  Created on: Jun 19, 2014
 *      Author: e4k2
 */

#include "InStructureCondition.h"
#include "ContactMap.h"
#include "PseudoContact.h"
#include "PCAssignment.h"
#include "Main.h"
#include <cstdio>
#include <algorithm>
#include <tr1/unordered_set>

InStructureCondition::InStructureCondition(double distCutoff, int countCutoff, ContactMap& cm) :
                          distCutoff(distCutoff), countCutoff(countCutoff), currentCount(0), contactMap(cm)
{
}

InStructureCondition::InStructureCondition(const InStructureCondition& isc) :
		distCutoff(isc.distCutoff), countCutoff(isc.countCutoff), currentCount(isc.currentCount), contactMap(isc.contactMap)
{
}

InStructureCondition::~InStructureCondition()
{
}

void swap(InStructureCondition& isc1, InStructureCondition& isc2)
{
	using namespace std;
	swap(isc1.distCutoff, isc2.distCutoff);
	swap(isc1.countCutoff, isc2.countCutoff);
	swap(isc1.currentCount, isc2.currentCount);
	swap(isc1.contactMap, isc2.contactMap);
}

InStructureCondition& InStructureCondition::operator=(InStructureCondition isc)
{
	swap(*this, isc);
	return *this;
}

bool InStructureCondition::test(PCAssignment* pca)
{
	if (currentCount > countCutoff) // all subsequent assignments are removed once cutoff is reached
		return false;

	PseudoContact& pc = pca->pc;
	int res1 = 0;
	int res2 = 0;
	pc.getResPair(res1,res2);
	int diff = abs(res1-res2);
	if (diff < 6)
		return true; // short range always returns true unless currentCount > countCutoff

	// get min dist in templates
	double bestMinDist = INVALIDDISTANCE;
	tr1::unordered_set<Contact>& contacts = pca->pc.getContacts();
	for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
	{
		const Contact& c = *itC;
		tr1::array<double,ContactMap::CMENTRYSIZE>& entry = contactMap[c];
		double minDist = entry[ContactMap::MINDISTINDEX];
		if (minDist < bestMinDist)
			bestMinDist = minDist;
	}
	if (bestMinDist > distCutoff) // distance violation
	{
		currentCount++;
		if (currentCount > countCutoff)
			return false;
		else
			return true;
	}
	return true;
}

void InStructureCondition::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("InStructureCondition: distCutoff=%6.3f countCutoff=%d currentCount=%d\n",
			distCutoff, countCutoff, currentCount);
}

void InStructureCondition::reset()
{
	currentCount = 0;
}
