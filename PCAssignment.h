/*
 * PCAssignment.h
 *
 *  Created on: 2012-05-09
 *      Author: e4k2
 */

#ifndef PCASSIGNMENT_H_
#define PCASSIGNMENT_H_

#include "PseudoContact.h"
#include "Score.h"

class NOE;

/**
 * PseudoContact - NOE assignment
 */
class PCAssignment
{
public:
	NOE* noe;
	PseudoContact pc;
	Score score;

	PCAssignment();
	PCAssignment(NOE* n, const PseudoContact& p, const Score& s);
	PCAssignment(const PCAssignment& a);
	friend void swap(PCAssignment& first, PCAssignment& second);
	// PCAssignment(PCAssignment&& a);
	PCAssignment& operator=(PCAssignment a);
	void print() const;
	virtual ~PCAssignment();
	bool operator==(const PCAssignment& a) const;
	bool operator!=(const PCAssignment& a) const;
	void setReversePC(); // reverses the direction of the contacts in pc. If this is done, the direction of the chemical shifts in noe
	            // will not match the order of the chemical shifts of the atoms in pc, so care must be taken to avoid functions
	           // that compare the chemical shifts
};

namespace std
{
	namespace tr1
	{
		template<>
		class hash<PCAssignment>
		{
			public:
			size_t operator()(const PCAssignment& a) const
			{
				// hash(i)=i*2654435761 mod 2^32
				return hash<long>()(long(a.noe)) ^ hash<PseudoContact>()(a.pc);
			}
		};
	}
}

struct PCAssignmentSortByScore
{
	bool operator()(PCAssignment* a, PCAssignment *b)
	{
		return a->score.total < b->score.total;
	}
};

#endif /* PCASSIGNMENT_H_ */
