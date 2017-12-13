/*
 * InStructureConditionFract.cpp
 *
 *  Created on: Jun 27, 2014
 *      Author: e4k2
 */

#include "InStructureConditionFract.h"
#include "ContactMap.h"
#include "PseudoContact.h"
#include "PCAssignment.h"
#include "Main.h"
#include <cstdio>
#include <algorithm>
#include <tr1/unordered_set>

InStructureConditionFract::InStructureConditionFract(int minCount, double distCutoff,
		double fractCutoff, ContactMap& contactMap, int maxCount, int seqSep) : minCount(minCount),
		currentNumViols(0), currentCount(0), distCutoff(distCutoff), fractCutoff(fractCutoff),
		contactMap(contactMap), maxCount(maxCount), seqSep(seqSep), numTimes(0)
{
}

InStructureConditionFract::InStructureConditionFract(const InStructureConditionFract& iscf) :
		minCount(iscf.minCount), currentNumViols(iscf.currentNumViols), currentCount(iscf.currentCount),
		distCutoff(iscf.distCutoff), fractCutoff(iscf.fractCutoff), contactMap(iscf.contactMap),
		maxCount(iscf.maxCount), seqSep(iscf.seqSep), numTimes(iscf.numTimes)
{
}

InStructureConditionFract::~InStructureConditionFract()
{
}

void swap(InStructureConditionFract& iscf1, InStructureConditionFract& iscf2)
{
	using namespace std;
	swap(iscf1.minCount, iscf2.minCount);
	swap(iscf1.currentNumViols, iscf2.currentNumViols);
	swap(iscf1.currentCount, iscf2.currentCount);
	swap(iscf1.distCutoff, iscf2.distCutoff);
	swap(iscf1.fractCutoff, iscf2.fractCutoff);
	swap(iscf1.contactMap, iscf2.contactMap);
	swap(iscf1.maxCount, iscf2.maxCount);
	swap(iscf1.seqSep, iscf2.seqSep);
	swap(iscf1.numTimes, iscf2.numTimes);
}

InStructureConditionFract& InStructureConditionFract::operator=(InStructureConditionFract iscf)
{
	swap(*this, iscf);
	return *this;
}

bool InStructureConditionFract::test(PCAssignment* pca)
{
	if (numTimes > maxCount)
		return false;
	if (currentCount > 0)
	{
		double curFract = double(currentNumViols)/double(currentCount);
		if (currentCount > minCount && curFract > fractCutoff)
			return false; // once cutoff has been reached all subsequent assignments should be removed
	}
	else if (fractCutoff == 0)
	{
		return false;
	}

	PseudoContact& pc = pca->pc;
	int res1 = 0;
	int res2 = 0;
	pc.getResPair(res1,res2);
	int diff = abs(res1-res2);
	if (diff < seqSep)
	{
		numTimes++;
		if (numTimes > maxCount)
			return false;
		return true; // all short range contacts are allowed if cutoff or max has not been reached
	}

	// can assert here pca is long range contact

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
	currentCount++;
	if (bestMinDist > distCutoff) // distance violation
	{
		currentNumViols++;
		double curFract = double(currentNumViols)/double(currentCount);
		if (currentCount > minCount && curFract > fractCutoff)
		{
			return false;
		}
	}
	numTimes++;
	if (numTimes > maxCount)
		return false;
	return true;
}


void InStructureConditionFract::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("InStructureConditionFract: minCount=%d currentNumViols=%d currentCount=%d distCutoff=%6.3f fractCutoff=%6.4f maxCount=%d seqSep=%d\n",
			minCount, currentNumViols, currentCount, distCutoff, fractCutoff, maxCount, seqSep);
}

void InStructureConditionFract::reset()
{
	currentNumViols = 0;
	currentCount = 0;
	numTimes = 0;
}
