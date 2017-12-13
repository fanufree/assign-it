/*
 * ContactMap.cpp
 *
 *  Created on: Jun 20, 2014
 *      Author: e4k2
 */

#include "ContactMap.h"
#include "ContactMapIterator.h"
#include <algorithm>

ContactMap::ContactMap()
{
}

ContactMap::ContactMap(const ContactMap& cm) : contactMap(cm.contactMap)
{
}

ContactMap::~ContactMap()
{
}

void swap(ContactMap& cm1, ContactMap& cm2)
{
	using namespace std;
	swap(cm1.contactMap, cm2.contactMap);
}


ContactMap& ContactMap::operator=(ContactMap cm)
{
	swap(*this, cm);
	return *this;
}

void ContactMap::add(Contact& c, tr1::array<double,CMENTRYSIZE>& entry)
{
	contactMap[c] = entry;
}

tr1::array<double,ContactMap::CMENTRYSIZE>& ContactMap::operator[](const Contact& c)
{
	return contactMap[c];
}

tr1::array<double, ContactMap::CMENTRYSIZE> ContactMap::makeEntry(double fractionContact, double minDist, double fractionContactSphere, double avgDist, double stdevDist)
{
	tr1::array<double,CMENTRYSIZE> entry;
	entry[FRACSTRUCINDEX] = fractionContact;
	entry[MINDISTINDEX] = minDist;
	entry[FRACSTRUCSPHINDEX] = fractionContactSphere;
	entry[AVGINDEX] = avgDist;
	entry[STDEVINDEX] = stdevDist;
	return entry;
}

double ContactMap::get(Contact& c, int fieldIndex)
{
	return contactMap[c][fieldIndex];
}

ContactMapIterator ContactMap::begin() const
{
	return ContactMapIterator(contactMap.begin(),this);
}

ContactMapIterator ContactMap::end() const
{
	return ContactMapIterator(contactMap.end(),this);
}

ContactMapIterator ContactMap::find(const Contact& c)
{
	tr1::unordered_map< Contact, tr1::array<double,CMENTRYSIZE> >::iterator itFind = contactMap.find(c);
	return ContactMapIterator(itFind,this);
}

int ContactMap::getNumContacts() const
{
	return contactMap.size();
}

int ContactMap::getNumContacts(double distCutoff) const
{
	using namespace std;
	int count = 0;
	for (tr1::unordered_map< Contact, tr1::array< double, CMENTRYSIZE > >::const_iterator it = contactMap.begin(); it != contactMap.end(); ++it)
	{
		const tr1::array<double,CMENTRYSIZE>& arr = it->second;
		if (arr[MINDISTINDEX] <= distCutoff)
			count++;
	}
	return count;
}
