/*
 * InStructureConditionFractTemp.cpp
 *
 *  Created on: Jun 29, 2014
 *      Author: e4k2
 */

#include "InStructureConditionFractTemp.h"
#include <tr1/unordered_set>
#include <tr1/array>
#include <cstdio>
#include <cstdlib>
#include <stddef.h>
#include "ContactMap.h"
#include "CSProtein.h"
#include "Residue.h"
#include "Atom.h"
#include "PseudoContact.h"
#include "PCAssignment.h"

InStructureConditionFractTemp::InStructureConditionFractTemp(double distCutoff, double fractTempCutoff, list<CSProtein*>& structures)
	: distCutoff(distCutoff), fractTempCutoff(fractTempCutoff)
{
	// make contact map based on parameters
	contactMap = new ContactMap();
	CSProtein* csa = structures.front();
	double numStructuresF = (double)(structures.size());
	for (int r1 = 1; r1 <= csa->size; r1++)
	{
		Residue* res1 = (*csa)[r1];
		for (AtomIterator itX1 = res1->begin('X'); itX1 != res1->end('X'); itX1++)
		{
			Atom* x1Atom = *itX1;
			if (x1Atom->getNumProtons() < 1)
				continue; // no proton
			bool isMethylX1 = x1Atom->isMethyl();
			Atom* hx1Atom = NULL; // used to ensure we count methyl groups only once
			for (HIterator itHX1 = x1Atom->beginH(); itHX1 != x1Atom->endH(); itHX1++)
			{
				if (hx1Atom != NULL && isMethylX1)
					break; // count the protons of a methyl group as one pseudoproton
				hx1Atom = *itHX1;
				for (int r2 = r1; r2 <= csa->size; r2++)
				{
					Residue* res2 = (*csa)[r2];
					// NOEs can be in same residue
					// iterator over all heavy atoms, then get their H's
					for (AtomIterator itX2 = res2->begin('X'); itX2 != res2->end('X'); itX2++)
					{
						Atom* x2Atom = *itX2;
						if (x2Atom->getNumProtons() < 1)
							continue;
						bool isMethylX2 = x2Atom->isMethyl();
						Atom* hx2Atom = NULL;
						for (HIterator itHX2 = x2Atom->beginH(); itHX2 != x2Atom->endH(); itHX2++)
						{
							if (hx2Atom != NULL && isMethylX2)
								break;
							hx2Atom = *itHX2;
							if (r1 == r2 && hx1Atom >= hx2Atom)
								continue; // avoid double counting and self contacts

							Contact contact(res1,x1Atom,hx1Atom,res2,x2Atom,hx2Atom);

							int contactCount = 0;

							for (list<CSProtein*>::const_iterator itS = structures.begin(); itS != structures.end(); itS++)
							{
								CSProtein* structure = *itS;
								double dist = 0;
								if (structure->inContact(contact,distCutoff,dist)) // inContact considers methyl's by taking min dist;
									contactCount++;
							}
							double fractionContact = double(contactCount)/numStructuresF;
							tr1::array<double,ContactMap::CMENTRYSIZE> entry = ContactMap::makeEntry(fractionContact, 0, 0, 0, 0);// only fractionContact used
							contactMap->add(contact, entry);
							Contact cr = contact.reverse();
							contactMap->add(cr, entry);
						} // end for each h attached to N or C of res2
					} // end for each N or C atom of res1
				} // end for eah res2
			} // end for each h attached to N or C of res1
		} // end for each N or C atom of res1
	} // end for each res1
}

InStructureConditionFractTemp::InStructureConditionFractTemp(const InStructureConditionFractTemp& is)
	: distCutoff(is.distCutoff), fractTempCutoff(is.fractTempCutoff), contactMap(is.contactMap)
{
}

InStructureConditionFractTemp::~InStructureConditionFractTemp()
{
	delete contactMap;
}

void swap(InStructureConditionFractTemp& is1, InStructureConditionFractTemp& is2)
{
	using namespace std;
	swap(is1.distCutoff, is2.distCutoff);
	swap(is1.fractTempCutoff, is2.fractTempCutoff);
	swap(is1.contactMap, is2.contactMap);
}

InStructureConditionFractTemp& InStructureConditionFractTemp::operator=(InStructureConditionFractTemp is)
{
	swap(*this, is);
	return *this;
}

bool InStructureConditionFractTemp::test(PCAssignment* pca)
{
	PseudoContact& pc = pca->pc;
	int res1 = 0;
	int res2 = 0;
	pc.getResPair(res1,res2);
	int diff = abs(res1-res2);
	if (diff < 2)
		return true;

	tr1::unordered_set<Contact>& contacts = pca->pc.getContacts();
	for (tr1::unordered_set<Contact>::iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
	{
		const Contact& c = *itC;
		tr1::array<double,ContactMap::CMENTRYSIZE>& entry = (*contactMap)[c];
		if (entry[ContactMap::FRACSTRUCINDEX] >= fractTempCutoff)
			return true;
	}
	return false;
}

void InStructureConditionFractTemp::print(int indent)
{
	for (int i = 0; i < indent; i++)
		printf("\t");
	printf("InStructureConditionFractTemp: distCutoff=%6.3f fractTempCutoff=%6.4f\n", distCutoff, fractTempCutoff);
}
