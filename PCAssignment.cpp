/*
 * PCAssignment.cpp
 *
 *  Created on: 2012-05-09
 *      Author: e4k2
 */

#include <cstdio>
#include "PCAssignment.h"
#include "PseudoContact.h"
#include "Score.h"
#include "NOE.h"

PCAssignment::PCAssignment() : noe(0), pc(), score()
{
}

PCAssignment::PCAssignment(NOE* n, const PseudoContact& p, const Score& s) : noe(n), pc(p), score(s)
{
}

PCAssignment::PCAssignment(const PCAssignment& a) : noe(a.noe), pc(a.pc), score(a.score)
{
}

void swap(PCAssignment& first, PCAssignment& second)
{
	std::swap(first.noe,second.noe);
	swap(first.pc,second.pc);
	swap(first.score,second.score);
}

//PCAssignment::PCAssignment(PCAssignment&& a) : noe(0), pc(), score()
//{
//	swap(*this,a);
//}

PCAssignment& PCAssignment::operator=(PCAssignment a)
{
	swap(*this,a);
	return *this;
}

void PCAssignment::print() const
{
	printf("BEGIN ASSIGNMENT\n");
	noe->print();
	pc.print();
	score.print();
	printf("END ASSIGNMENT\n");
}

PCAssignment::~PCAssignment()
{
}

bool PCAssignment::operator==(const PCAssignment& a) const
{
	return noe == a.noe && pc == a.pc;
}

bool PCAssignment::operator!=(const PCAssignment& a) const
{
	return !(*this == a);
}

void PCAssignment::setReversePC()
{
	pc.setReverse();
}
