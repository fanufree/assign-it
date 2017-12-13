/*
 * ContactMapIterator.h
 *
 *  Created on: Jun 20, 2014
 *      Author: e4k2
 */

#ifndef CONTACTMAPITERATOR_H_
#define CONTACTMAPITERATOR_H_

#include "Contact.h"
#include "ContactMap.h"
#include <tr1/unordered_map>
#include <tr1/array>


/**
 * Read-only iterator. Cannot be used to edit the contact map
 */
class ContactMapIterator
{
public:
	ContactMapIterator(tr1::unordered_map< Contact, tr1::array<double,ContactMap::CMENTRYSIZE> >::const_iterator i, const ContactMap* c);
	~ContactMapIterator();
	ContactMapIterator(const ContactMapIterator& i);
	friend void swap(ContactMapIterator& it1, ContactMapIterator& it2);
	ContactMapIterator& operator=(ContactMapIterator i);
	ContactMapIterator& operator++();
	ContactMapIterator operator++(int); // postfix
	bool operator==(const ContactMapIterator& i) const;
	bool operator!=(const ContactMapIterator& i) const;
	const Contact& first();
	const tr1::array< double,ContactMap::CMENTRYSIZE >& second();
	// pair< Contact,tr1::array< double,ContactMap::CMENTRYSIZE > >& operator->();
private:
	tr1::unordered_map< Contact, tr1::array<double,ContactMap::CMENTRYSIZE> >::const_iterator it;
	const ContactMap* cm;
};

#endif /* CONTACTMAPITERATOR_H_ */
