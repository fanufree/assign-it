/*
 * Interval.h
 *
 *  Created on: 2012-10-12
 *      Author: e4k2
 */

#ifndef INTERVAL_H_
#define INTERVAL_H_

#include <vector>

using namespace std;

class Interval
{
public:
	Interval();
	Interval(double left, double right);
	Interval(const Interval& i);
	virtual ~Interval();
	friend void swap(Interval& i1, Interval& i2);
	// Interval(Interval&& i);
	Interval& operator=(Interval i);
	void print();
	bool intersect(const Interval& in); // true if this and in intersect
	bool intersect(const Interval& in, Interval& ret); // ret is the intersection if
	                         // this and in intersect; otherwise it is unmodified
	                         // Returns true if this and in intersect
	                         // ret == this is a possible input
	static void getMaximumClique(vector<Interval*>& input, Interval& output,
			int& maxCliqueSize, vector<Interval*>& intsMaxClique); // intsMaxClique = intervals
	                             // from input in max clique
	                             // Side-effects: input is sorted by left end point
	// Assumes ref and input are sorted in increasing order
	static void calibrate(const vector<double>& ref, const vector<double>& input,
			double matchH, double searchH, double& offsetH, bool debug); // return the calibration offset to add to input
	//static void test1();
	//static void test2();
	double l; // left end point
	double r; // right end point
	double d; // to be annotated by calibrate
};


struct IntervalSortByLeftEnd
{
	bool operator()(Interval* a, Interval* b)
	{
		if (a->l != b->l)
			return a->l < b->l;
		else
			return a->r < b->r; // break ties by right end point
	}
};

struct IntervalSortByD
{
	bool operator()(Interval* a, Interval* b)
	{
		return a->d < b->d;
	}
};

#endif /* INTERVAL_H_ */
