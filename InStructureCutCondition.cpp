/*
 * InStructureCutCondition.cpp
 *
 *  Created on: Sep 12, 2014
 *      Author: e4k2
 */

#include "InStructureCutCondition.h"
#include "ContactMap.h"
#include "PCAssignment.h"
#include <cstdio>
#include <algorithm>
#include <tr1/unordered_set>

InStructureCutCondition::InStructureCutCondition(double fractCutoff, double fractNumAssignedCut,
		int numAssignments, ContactMap& contactMap) : fractCutoff(fractCutoff),
		  fractNumAssignedCut(fractNumAssignedCut), numAssignments(numAssignments), contactMap(contactMap), numViols(0)
{
}

InStructureCutCondition::InStructureCutCondition(const InStructureCutCondition& isc) : fractCutoff(isc.fractCutoff),
		fractNumAssignedCut(isc.fractNumAssignedCut), numAssignments(isc.numAssignments),
		contactMap(isc.contactMap), numViols(isc.numViols)
{
}

InStructureCutCondition::~InStructureCutCondition()
{
}

void swap(InStructureCutCondition& isc1, InStructureCutCondition& isc2)
{
	using namespace std;
	swap(isc1.fractCutoff, isc2.fractCutoff);
	swap(isc1.fractNumAssignedCut, isc2.fractNumAssignedCut);
	swap(isc1.numAssignments, isc2.numAssignments);
	swap(isc1.contactMap, isc2.contactMap);
	swap(isc1.numViols, isc2.numViols);
}

InStructureCutCondition& InStructureCutCondition::operator=(InStructureCutCondition isc)
{
	swap(*this, isc);
	return *this;
}

bool InStructureCutCondition::test(PCAssignment* pca)
{
	if (numViols > fractNumAssignedCut*numAssignments)
		return false;

	// get fract structures that meet distance cutoff
	// return true if any contact in pc meets this
	tr1::unordered_set<Contact>& contacts = pca->pc.getContacts();
	bool ok = false;
	for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
	{
		const Contact& c = *itC;
		tr1::array<double,ContactMap::CMENTRYSIZE>& entry = contactMap[c];
		double fract = entry[ContactMap::FRACSTRUCINDEX];
		if (fract >= fractCutoff)
		{
			ok = true;
			break;
		}
	}
	if (!ok)
	{
		numViols++;
		if (numViols > fractNumAssignedCut*numAssignments)
			return false;
		else
			return true;
	}
	else
	{
		return true;
	}
}

void InStructureCutCondition::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("InStructureCutCondition: fractCutoff=%6.4f fractNumAssignedCut=%6.4f numAllowViols=%d numViols=%d\n",
			 fractCutoff, fractNumAssignedCut, numAssignments, numViols);
}

void InStructureCutCondition::reset()
{
	numViols = 0;
}



