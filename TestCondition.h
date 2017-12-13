/*
 * TestCondition.h
 *
 *  Created on: Jun 18, 2014
 *      Author: e4k2
 */

#ifndef TESTCONDITION_H_
#define TESTCONDITION_H_

class PCAssignment;

class TestCondition
{
public:
	TestCondition();
	virtual ~TestCondition();

	virtual bool test(PCAssignment* pca) = 0; // returns true if should keep pca; false = eliminate pca
	virtual void print(int indent) = 0;
	virtual void reset(); // resets the state of this test; i.e. after testing all the pca's, we do the tests again, so we need to reset the state
};

#endif /* TESTCONDITION_H_ */
