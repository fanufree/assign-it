/*
 * ContactMap.h
 *
 *  Created on: Jun 20, 2014
 *      Author: e4k2
 */

#ifndef CONTACTMAP_H_
#define CONTACTMAP_H_

#include "Contact.h"
#include <tr1/array>

using namespace std;

class ContactMapIterator;

/**
 * Contains all contacts (distance determined by caller) added by the user using makeEntry and add
 */

class ContactMap
{
public:
	static const int CMENTRYSIZE = 5; // contact map entry size
	static const int FRACSTRUCINDEX = 0; // index for fraction of input structures with contact (within DISTCUTOFF)
	static const int MINDISTINDEX = 1; // index for minimum distance among the structures for this contact; can be INVALIDDISTANCE if coordinates do not exist
	static const int FRACSTRUCSPHINDEX = 2; //  index for fraction of input structures with contact, allowing for intersectSphere
	static const int AVGINDEX = 3; // avg dist; if a structure is missing coordinates, the avg dist is set to INVALIDDISTANCE
	static const int STDEVINDEX = 4;
	tr1::unordered_map< Contact, tr1::array<double,CMENTRYSIZE> > contactMap;

	ContactMap();
	ContactMap(const ContactMap& cm);
	virtual ~ContactMap();
	friend void swap(ContactMap& cm1, ContactMap& cm2);
	ContactMap& operator=(ContactMap cm);
	void add(Contact& c, tr1::array<double,CMENTRYSIZE>& entry); // definition of contact depends on caller
	tr1::array<double,CMENTRYSIZE>& operator[](const Contact& c);
	static tr1::array<double, CMENTRYSIZE> makeEntry(double fractionContact, double minDist, double fractionContactSphere, double avgDist, double stdevDist);
	double get(Contact& c, int fieldIndex); // fieldIndex = must be one of the above static const ints (no checking is done)
	ContactMapIterator find(const Contact& c);
	int getNumContacts() const;
	int getNumContacts(double distCutoff) const; // num contacts where mindist <= distCutoff

	friend class ContactMapIterator;
	ContactMapIterator begin() const;
	ContactMapIterator end() const;
};

#endif /* CONTACTMAP_H_ */
