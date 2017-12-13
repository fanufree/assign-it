/*
 * Assignment.h
 *
 *  Created on: 2012-05-06
 *      Author: e4k2
 */

#ifndef ASSIGNMENT_H_
#define ASSIGNMENT_H_

#include "Contact.h"
#include "NOE.h"
#include "Score.h"

using namespace std;


class Assignment
{
public:
	NOE* noe;
	Contact contact;
	Score score;
	Assignment();
	Assignment(NOE* noe, const Contact& c, const Score& s);
	Assignment(const Assignment& a);
	// Assignment(Assignment&& a);
	friend void swap(Assignment& first, Assignment& second);
	Assignment& operator=(Assignment a);
	void print() const;
	virtual ~Assignment();
	bool operator==(const Assignment& a) const;
	bool operator!=(const Assignment& a) const;
};

namespace std
{
	namespace tr1
	{
		template<>
		class hash<Assignment>
		{
			public:
			size_t operator()(const Assignment& a) const
			{
				// hash(i)=i*2654435761 mod 2^32
				return hash<long>()(long(a.noe)) ^ hash<Contact>()(a.contact);
			}
		};
	}
}

#endif /* ASSIGNMENT_H_ */
