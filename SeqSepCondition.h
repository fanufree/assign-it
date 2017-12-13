/*
 * SeqSepCondition.h
 *
 *  Created on: Jun 29, 2014
 *      Author: e4k2
 */

#ifndef SEQSEPCONDITION_H_
#define SEQSEPCONDITION_H_

#include "TestCondition.h"

class SeqSepCondition: public TestCondition
{
	int seqSepCutoff; // returns true if sequence separation is < seqSepCutoff
	                  // if 0, always returns false
public:
	SeqSepCondition(int seqSep);
	SeqSepCondition(const SeqSepCondition& ssc);
	virtual ~SeqSepCondition();
	friend void swap(SeqSepCondition& ssc1, SeqSepCondition& ssc2);
	SeqSepCondition& operator=(SeqSepCondition ssc);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
};

#endif /* SEQSEPCONDITION_H_ */
