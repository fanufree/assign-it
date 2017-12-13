/*
 * Interval.cpp
 *
 *  Created on: 2012-10-12
 *      Author: e4k2
 */

#include <algorithm>
#include <cstdio>
#include <list>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <tr1/unordered_map>
#include "Interval.h"

Interval::Interval() : l(0), r(0), d(0)
{
}

Interval::Interval(const Interval& i) : l(i.l), r(i.r), d(i.d)
{
}

Interval::Interval(double left, double right) : l(left), r(right), d(0)
{
}

void swap(Interval& i1, Interval& i2)
{
	using std::swap;
	swap(i1.l,i2.l);
	swap(i1.r,i2.r);
	swap(i1.d,i2.d);
}

//Interval::Interval(Interval&& i) :  l(0), r(0), d(0)
//{
//	swap(*this,i);
//}

Interval& Interval::operator=(Interval i)
{
	swap(*this,i);
	return *this;
}

Interval::~Interval()
{
}


void Interval::print()
{
	printf("[%6.4f, %6.4f]\n",l,r);
}

bool Interval::intersect(const Interval& in)
{
	if (in.l <= r && l <= in.r)
	{
		return true;
	}
	return false;
}

bool Interval::intersect(const Interval& in, Interval& ret)
{
	if (!intersect(in))
		return false;
	/*
	 * Default case
	 * in         L------R
	 * this           L------R
	 */
	double left = l;
	double right = in.r;
	if (l < in.l)
	{
		if (r < in.r)
		{
			/*
			 * in           L--------R
			 * this     L------R
			 */
			left = in.l;
			right = r;
		}
		else
		{
			/*
			 * in           L----R
			 * this     L-----------R
			 */
			left = in.l;
			right = in.r;
		}
	}
	else if (r < in.r)
	{
		/*
		 * in      L-------------------R
	     * this         L-----R
		 */
		left = l;
		right = r;
	}
	ret.l = left;
	ret.r = right;
	return true;
}

void Interval::getMaximumClique(vector<Interval*>& input, Interval& output,
			int& maxCliqueSize, vector<Interval*>& intsMaxClique)
{
	sort(input.begin(),input.end(),IntervalSortByLeftEnd());
	maxCliqueSize = 1;
	vector<Interval*>::const_iterator it1Max = input.begin()+1; // used to retrieve the max clique
	vector<Interval*>::const_iterator it2MaxStart = input.begin(); // used to retrieve the max clique
	vector<Interval*>::const_iterator it2Start = input.begin(); // used to speed up the intersection calculations
	Interval inter; // not used

	for (vector<Interval*>::const_iterator it1 = input.begin()+1; it1 != input.end(); ++it1)
	{
		Interval* i1 = *it1;
		bool firstTimeSet = false;
		int cliqueSize = 1; // includes i1
		// check intersection with previous intervals (by left end point);
		for (vector<Interval*>::const_iterator it2 = it2Start; it2 != it1; ++it2)
		{
			Interval* i2 = *it2;
			if (i1->intersect(*i2,inter))
			{
				cliqueSize++;
				if (!firstTimeSet)
				{
					it2Start = it2;
					firstTimeSet = true;
				}
			}
		}
		if (cliqueSize > maxCliqueSize)
		{
			maxCliqueSize = cliqueSize;
			it1Max = it1;
			it2MaxStart = it2Start;
		}
	}
	if (maxCliqueSize > 1)
	{
		intsMaxClique.clear(); // make sure empty
		Interval* i1 = *it1Max;
		intsMaxClique.push_back(i1);
		output = *i1;
		for (vector<Interval*>::const_iterator it2 = it2MaxStart; it2 != it1Max; ++it2)
		{
			Interval* i2 = *it2;
			if (i1->intersect(*i2,inter))
			{
				intsMaxClique.push_back(i2);
				output.intersect(inter,output);
			}
		}
	}
	else if (maxCliqueSize == 1)
	{
		output = *(intsMaxClique.front());
	}
}

/*
 * ref should be the assigned chemical shifts
 * input should be the NOEs
 * Assumes ref and input are sorted in increasing order
 */
void Interval::calibrate(const vector<double>& ref, const vector<double>& input, double matchH, double searchH, double& offsetH, bool debug)
{
	vector<Interval*> intervalMatches;
	offsetH = 0;
	// for efficiency reasons, we will walk through nIn<-noes until nIn > hRef, then we update hRef one step
	// nIn is set to the first point of intersection of the previous hRef since it can also intersect with
	// this hRef in addition to the previous hRef
	vector<double>::const_iterator itInputStart = input.begin();
	for (vector<double>::const_iterator it = ref.begin(); it != ref.end(); ++it) // chem shift assign
	{
		Interval refIn(-searchH,searchH);
		double hRef = *it;
		bool firstTimeSet = false; // update itInputStart to the first hIn of intersection
		// next hRef will use this hIn as the starting point
		for (vector<double>::const_iterator itIn = itInputStart; itIn != input.end(); ++itIn) // noes
		{
			double hIn = *itIn;
			double dH = fabs(hRef-hIn);
			double tol = matchH+searchH;
			if (dH <= tol)
			{
				Interval inIn(hIn-hRef-matchH, hIn-hRef+matchH);
				Interval* intersect = new Interval;
				if (refIn.intersect(inIn,*intersect))
				{
					intersect->d = hRef-hIn;
					intervalMatches.push_back(intersect);
					if (!firstTimeSet)
					{
						itInputStart = itIn;
						firstTimeSet = true;
					}
				}
				else
				{
					delete intersect;
					break; // no need to traverse input further since hIn will keep increasing
				}
			}
		}
	}
	if (intervalMatches.size() > 0)
	{
		Interval output;
		int maxCliqueSize = 0;
		vector<Interval*> intsMaxClique;
		getMaximumClique(intervalMatches, output, maxCliqueSize, intsMaxClique);
		if (maxCliqueSize > 0)
		{
			tr1::unordered_map<int, int> counts;
			int bestBin = -1; // for finding the most common bin
			int bestBinCount = 0;
			for (vector<Interval*>::iterator it = intsMaxClique.begin(); it != intsMaxClique.end(); ++it)
			{
				Interval* iv = *it;
				// printf("%f\t%f\t\%f\n", iv->d, iv->l, iv->r);
				int v = (int)round(iv->d*100); // group by 0.0X value
				// printf("%d %f\n",v,iv->d);
				tr1::unordered_map<int,int>::iterator itC = counts.find(v);
				int count = 1;
				if (itC != counts.end())
					count = counts[v]+1;
				counts[v] = count;
				if (count > bestBinCount)
				{
					bestBinCount = count;
					bestBin = v;
				}
			}
			// find the most common bin
			vector<double> ds;
			ds.reserve(intsMaxClique.size()/2);
			for (vector<Interval*>::iterator it = intsMaxClique.begin(); it != intsMaxClique.end(); ++it)
			{
				Interval* iv = *it;
				int v = (int)round(iv->d*100);
				if (v == bestBin)
					ds.push_back(iv->d);
			}
			sort(ds.begin(),ds.end());
			// printf("Most common bin is %d bestCount=%d\n",bestBin,bestBinCount);
			offsetH = ds[ds.size()/2];

			// get first quartile displacement vector; filter outliers
//			sort(intsMaxClique.begin(),intsMaxClique.end(),IntervalSortByD());
//			double midpt = intsMaxClique[maxCliqueSize/2]->d;
//			if (midpt > 0)
//				offsetH = intsMaxClique[maxCliqueSize/4]->d;
//			else
//				offsetH = intsMaxClique[3*maxCliqueSize/4]->d;
			if (debug)
				printf("calibration with %6.3f\n",offsetH);
		}
		else
		{
			if (debug)
				printf("Unable to calibrate. No clique\n");
		}
	}
	else
	{
		if (debug)
			printf("Unable to calibrate. No chemical shift matches possible\n");
	}
	for (vector<Interval*>::iterator it = intervalMatches.begin(); it != intervalMatches.end(); ++it)
	{
		delete *it;
	}
}

/**
MaxClique interval
[6.0000, 6.0000]
Intervals in max clique of size 4
[6.0000, 6.0000]
[2.0000, 6.0000]
[3.0000, 8.0000]
[5.0000, 10.0000]
 */
//void Interval::test1()
//{
//	vector<Interval*> input;
//	input.push_back(new Interval(1,4));
//	input.push_back(new Interval(2,6));
//	input.push_back(new Interval(3,8));
//	input.push_back(new Interval(5,10));
//	input.push_back(new Interval(7,9));
//	input.push_back(new Interval(9,10));
//	input.push_back(new Interval(6,6));
//	Interval output;
//	int maxCliqueSize = 0;
//	vector<Interval*> intsMaxClique;
//	Interval::getMaximumClique(input,output,maxCliqueSize,intsMaxClique);
//	printf("MaxClique interval\n");
//	output.print();
//	printf("Intervals in max clique of size %d\n",maxCliqueSize);
//	for (vector<Interval*>::iterator it = intsMaxClique.begin(); it != intsMaxClique.end(); ++it)
//	{
//		Interval * in = *it;
//		in->print();
//	}
//
//	for (vector<Interval*>::iterator it = input.begin(); it != input.end(); )
//	{
//		Interval* i = *it;
//		delete i;
//		it = input.erase(it);
//	}
//}

/**
 * global noise should be similar to offset (but not always exactly due to local noise)
 */
//void Interval::test2()
//{
//	vector<double> ref(12,0);
//	vector<double> input(12,0);
//	srand ( time(NULL) );
//	int global = (rand() % 1000); // 0-999
//	double globalNoise = double(global)/10000.0; // 0-0.0999
//
//	ref[0] = 8.91;
//	ref[1] = 4.36;
//	ref[2] = 2.12;
//	ref[3] = 0.61;
//	ref[4] = 0.79;
//	ref[5] = -0.13;
//	ref[6] = -0.07;
//	ref[7] = 1.30;
//	ref[8] = 4.36;
//	ref[9] = 0.78;
//	ref[10] = 1.99;
//	ref[11] = 8.46;
//	printf("GlobalNoise %6.3f\n",globalNoise);
//	for (int i = 0; i < 12; i++)
//	{
//		int local = (rand() % 1000);
//		double localNoise = min(0.04,globalNoise/4.0)*(double(local)/1000.0); // 0-0.999
//		input[i] = ref[i]+globalNoise+localNoise;
//	}
//
//	double matchH = 0.04;
//	double searchH = 0.1;
//	bool debug = true;
//	double offsetH = -1;
//
//    Interval::calibrate(ref, input, matchH, searchH, offsetH, debug);
//	printf("The offset is %6.3f\n",offsetH);
//}
