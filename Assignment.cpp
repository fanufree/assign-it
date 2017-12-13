/*
 * Assignment.cpp
 *
 *  Created on: 2012-05-06
 *      Author: e4k2
 */

#include <cstdio>
#include "Assignment.h"
#include "NOE.h"
#include "Contact.h"

Assignment::Assignment() : noe(0), contact(), score()
{
}

Assignment::~Assignment()
{
}

Assignment::Assignment(NOE* n, const Contact& c, const Score& s) : noe(n), contact(c), score(s)
{
}

void swap(Assignment& first, Assignment& second)
{
	std::swap(first.noe,second.noe);
	swap(first.contact,second.contact);
	swap(first.score,second.score);
}

Assignment::Assignment(const Assignment& a) : noe(a.noe), contact(a.contact), score(a.score)
{
}

//Assignment::Assignment(Assignment&& a) : noe(0), contact(), score()
//{
//	swap(*this,a);
//}

Assignment& Assignment::operator=(Assignment a)
{
	swap(*this,a);
	return *this;
}

void Assignment::print() const
{
	printf("Assignment: ");
	noe->print();
	printf("\t");
	contact.print();
	printf("\t");
	score.print();
}

bool Assignment::operator==(const Assignment& a) const
{
	if (noe == a.noe && contact == a.contact)
		return true;
	else
		return false;
}

bool Assignment::operator!=(const Assignment& a) const
{
	return !(*this == a);
}

