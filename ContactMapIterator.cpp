/*
 * ContactMapIterator.cpp
 *
 *  Created on: Jun 20, 2014
 *      Author: e4k2
 */

#include "ContactMapIterator.h"

ContactMapIterator::ContactMapIterator(tr1::unordered_map< Contact, tr1::array<double,ContactMap::CMENTRYSIZE> >::const_iterator i,
		 const ContactMap* c) : it(i), cm(c)
{
}

ContactMapIterator::ContactMapIterator(const ContactMapIterator& i) : it(i.it), cm(i.cm)
{
}

ContactMapIterator::~ContactMapIterator()
{
}

void swap(ContactMapIterator& it1, ContactMapIterator& it2)
{
	using namespace std;
	swap(it1.it, it2.it);
	swap(it1.cm, it2.cm);
}

ContactMapIterator& ContactMapIterator::operator=(ContactMapIterator i)
{
	swap(*this, i);
	return *this;
}

ContactMapIterator& ContactMapIterator::operator++()
{
	it++;
	return *this;
}

ContactMapIterator ContactMapIterator::operator++(int) // postfix
{
	ContactMapIterator cmIt = *this;
	++*this;
	return cmIt;
}

bool ContactMapIterator::operator==(const ContactMapIterator& i) const
{
	return (it==i.it) && (cm==i.cm);
}

bool ContactMapIterator::operator!=(const ContactMapIterator& i) const
{
	return !(*this==i);
}

const Contact& ContactMapIterator::first()
{
	return it->first;
}

const tr1::array< double,ContactMap::CMENTRYSIZE >& ContactMapIterator::second()
{
	return it->second;
}

//pair< Contact,tr1::array< double,ContactMap::CMENTRYSIZE > >& ContactMapIterator::operator->()
//{
//	return it->;
//}
