/*
 * InStructureWindow.cpp
 *
 *  Created on: Nov 19, 2014
 *      Author: e4k2
 */

#include "InStructureWindow.h"
#include "ContactMap.h"
#include "PseudoContact.h"
#include "PCAssignment.h"
#include "Main.h"
#include <cstdio>
#include <algorithm>
#include <tr1/unordered_set>

InStructureWindow::InStructureWindow(int winLen, double fractViol, ContactMap& contactMap, double distCutoff) :
     winLen(winLen), fractViol(fractViol), contactMap(contactMap), distCutoff(distCutoff), window(winLen,false), numViol(0), index(0)
{
}

InStructureWindow::InStructureWindow(const InStructureWindow& isw) :
     winLen(isw.winLen), fractViol(isw.fractViol), contactMap(isw.contactMap), distCutoff(isw.distCutoff), window(isw.window), numViol(isw.numViol), index(isw.index)
{
}

InStructureWindow::~InStructureWindow()
{
}

void swap(InStructureWindow& isw1, InStructureWindow& isw2)
{
	using namespace std;
	swap(isw1.winLen, isw2.winLen);
	swap(isw1.fractViol, isw2.fractViol);
	swap(isw1.contactMap, isw2.contactMap);
	swap(isw1.distCutoff, isw2.distCutoff);
	swap(isw1.window, isw2.window);
	swap(isw1.numViol, isw2.numViol);
	swap(isw1.index, isw2.index);
}

InStructureWindow& InStructureWindow::operator=(InStructureWindow isw)
{
	swap(*this, isw);
	return *this;
}

bool InStructureWindow::test(PCAssignment* pca)
{
	double fract = double(numViol)/double(winLen);
	if (fract > fractViol)
		return false;

	PseudoContact& pc = pca->pc;
	int res1 = 0;
	int res2 = 0;
	pc.getResPair(res1,res2);
	int diff = abs(res1-res2);
	if (diff < 2)
		return true; // all short range contacts are allowed if cutoff has not been reached

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
	bool replace = window[index]; // the value in window to be replaced
	if (bestMinDist > distCutoff) // distance violation
	{
		window[index] = true;
		if (!replace)
		{
			numViol++;
			fract = double(numViol)/double(winLen);
			if (fract > fractViol)
				return false;
		}
	}
	else
	{
		window[index] = false;
		if (replace)
			numViol--;
	}
	index++;
	if (index % winLen == 0)
		index = 0;
	return true;
}

void InStructureWindow::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("InStructureWindow: winLen=%d fractViol=%6.3f numViol=%d index=%d\n",winLen,fractViol,numViol,index);
}

void InStructureWindow::reset()
{
	window.assign(winLen,false);
	numViol = 0;
	index = 0;
}
