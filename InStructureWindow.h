/*
 * InStructureWindow.h
 *
 *  Created on: Nov 19, 2014
 *      Author: e4k2
 */

#ifndef INSTRUCTUREWINDOW_H_
#define INSTRUCTUREWINDOW_H_

#include "TestCondition.h"
#include "ContactMap.h"
#include <vector>

/*
 * Discard all further assignments if num violations in current window of assignments is > fraction of violations, where each violation is
 * min dist among the structures > DISTCUTOFF
 * Initially, the window consists of winLen non-violations
 * Ignores assignments with seq sep < 2
 */
class InStructureWindow : public TestCondition
{
private:
	int winLen;
	double fractViol;
	ContactMap& contactMap;
	double distCutoff;
	vector<bool> window; // [i] = true if it is a violation; wraps around
	int numViol; // num violations in window
	int index; // index of next entry to add; 0-based index into window

public:
	InStructureWindow(int winLen, double fractViol, ContactMap& contactMap, double distCutoff);
	InStructureWindow(const InStructureWindow& isw);
	virtual ~InStructureWindow();
	friend void swap(InStructureWindow& isw1, InStructureWindow& isw2);
	InStructureWindow& operator=(InStructureWindow isw);
	virtual bool test(PCAssignment* pca);
	virtual void print(int indent);
	virtual void reset();
};

#endif /* INSTRUCTUREWINDOW_H_ */
