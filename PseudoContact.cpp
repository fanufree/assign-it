/*
 * PseudoContact.cpp
 *
 *  Created on: 2012-05-03
 *      Author: e4k2
 */

#include <algorithm>
#include <cassert>
#include "PseudoContact.h"
#include "Utilities.h"


PseudoContact::PseudoContact() : contacts(), res1(-1), res2(-1)
{
}

PseudoContact::PseudoContact(const Contact& c) : contacts(), res1(c.r1->num), res2(c.r2->num)
{
	contacts.insert(c);
}

PseudoContact::~PseudoContact()
{
}

PseudoContact::PseudoContact(const PseudoContact& pc) : contacts(pc.contacts), res1(pc.res1), res2(pc.res2)
{
}

void swap(PseudoContact& first, PseudoContact& second)
{
	using std::swap;
	swap(first.contacts,second.contacts);
	swap(first.res1,second.res1);
	swap(first.res2,second.res2);
}

//PseudoContact::PseudoContact(PseudoContact&& pc) : contacts()
//{
//	swap(*this,pc);
//}

PseudoContact& PseudoContact::operator=(PseudoContact pc)
{
	swap(*this,pc);
	return *this;
}

void PseudoContact::add(const Contact& c)
{
	contacts.insert(c);
	if (res1 < 0)
	{
		res1 = c.r1->num;
		res2 = c.r2->num;
	}
	else
	{
		assert(c.r1->num == res1 && c.r2->num == res2);
	}
}

int PseudoContact::numContacts() const
{
	return contacts.size();
}

void PseudoContact::print() const
{
	for (tr1::unordered_set<Contact>::const_iterator it = contacts.begin(); it != contacts.end(); ++it)
	{
		it->print();
	}
}

bool PseudoContact::operator==(const PseudoContact& pc) const
{
	if (contacts.size() != pc.contacts.size())
		return false;
	for (tr1::unordered_set<Contact>::const_iterator it = pc.contacts.begin(); it != pc.contacts.end(); ++it)
	{
		if (contacts.find(*it) == contacts.end())
			return false;
	}
	return true;
}

bool PseudoContact::areSymmetric(const PseudoContact& pc) const
{
	if (contacts.size() != pc.contacts.size())
		return false;
	for (tr1::unordered_set<Contact>::const_iterator it = pc.contacts.begin(); it != pc.contacts.end(); ++it)
	{
		Contact c = it->reverse();
		if (contacts.find(c) == contacts.end())
			return false;
	}
	return true;
}

bool PseudoContact::operator!=(const PseudoContact& pc) const
{
	return !(*this == pc);
}

tr1::unordered_set<Contact>& PseudoContact::getContacts()
{
	return contacts;
}

Contact PseudoContact::getOneContact()
{
	return *(contacts.begin());
}

// returns res1 <= res2
void PseudoContact::getResPairOrdered(int& res1, int& res2) const
{
	if (this->res1 <= this->res2)
	{
		res1 = this->res1;
		res2 = this->res2;
	}
	else
	{
		res1 = this->res2;
		res2 = this->res1;
	}
}

void PseudoContact::getResPair(int& res1, int& res2) const
{
	res1 = this->res1;
	res2 = this->res2;
}

void PseudoContact::setReverse()
{
	std::swap(res1,res2);
	vector<Contact> temp(contacts.begin(),contacts.end());
	contacts.clear();
	for (vector<Contact>::iterator it = temp.begin(); it != temp.end(); ++it)
	{
		contacts.insert(it->reverse());
	}
}

bool PseudoContact::overlaps(const PseudoContact& pc) const
{
	if (res1 == pc.res1 && res2 == pc.res2)
	{
		for (tr1::unordered_set<Contact>::const_iterator it1 = contacts.begin(); it1 != contacts.end(); ++it1)
		{
			const Contact& c1 = *it1;
			for (tr1::unordered_set<Contact>::const_iterator it2 = pc.contacts.begin(); it2 != pc.contacts.end(); ++it2)
			{
				const Contact& c2 = *it2;
				if (c1 == c2)
					return true;
			}
		}
	}
	else if (res1 == pc.res2 && res2 == pc.res1)
	{
		for (tr1::unordered_set<Contact>::const_iterator it1 = contacts.begin(); it1 != contacts.end(); ++it1)
		{
			const Contact& c1 = *it1;
			for (tr1::unordered_set<Contact>::const_iterator it2 = pc.contacts.begin(); it2 != pc.contacts.end(); ++it2)
			{
				const Contact& c2 = it2->reverse();
				if (c1 == c2)
					return true;
			}
		}
	}
	return false;
}

void PseudoContact::addAll(const PseudoContact& pc)
{
	if (res1 == pc.res1 && res2 == pc.res2)
	{
		contacts.insert(pc.contacts.begin(),pc.contacts.end());
	}
	else if (res1 == pc.res2 && res2 == pc.res1)
	{
		for (tr1::unordered_set<Contact>::const_iterator it2 = pc.contacts.begin(); it2 != pc.contacts.end(); ++it2)
		{
			const Contact& c2 = it2->reverse();
			contacts.insert(c2);
		}
	}
}

int PseudoContact::getSeqSep() const
{
	return abs(res1-res2);
}
