/*
 * Rectangle.cpp
 *
 *  Created on: 2012-07-16
 *      Author: e4k2
 */

#include <algorithm>
#include <list>
#include <tr1/unordered_map>
#include <fstream>
#include <sstream>
#include <cstdlib>
// #include <random>
#include <cmath>
#include "Rectangle.h"
#include "Event.h"
#include "ActivePoint.h"
#include "Utilities.h"

using namespace std;

Rectangle::Rectangle() : xMin(0), xMax(0), yMin(0), yMax(0), width(0), height(0), dx(0), dy(0)
{
}

Rectangle::Rectangle(double xMin, double xMax, double yMin, double yMax) : xMin(xMin), xMax(xMax),
		yMin(yMin), yMax(yMax), dx(0), dy(0)
{
	width = xMax - xMin;
	height = yMax - yMin;
}

Rectangle::Rectangle(const Rectangle& r) : xMin(r.xMin), xMax(r.xMax), yMin(r.yMin), yMax(r.yMax),
		width(r.width), height(r.height), dx(r.dx), dy(r.dy)
{
}

Rectangle::~Rectangle()
{
}

void swap(Rectangle& first, Rectangle& second)
{
	using std::swap;
	swap(first.xMin,second.xMin);
	swap(first.xMax,second.xMax);
	swap(first.yMin,second.yMin);
	swap(first.yMax,second.yMax);
	swap(first.width,second.width);
	swap(first.height,second.height);
	swap(first.dx,second.dx);
	swap(first.dy,second.dy);
}

//Rectangle::Rectangle(Rectangle&& r) : xMin(0), xMax(0), yMin(0), yMax(0), width(0), height(0), dx(0), dy(0)
//{
//	swap(*this, r);
//}

Rectangle& Rectangle::operator=(Rectangle r)
{
	swap(*this,r);
	return *this;
}

bool Rectangle::contains(double x, double y)
{
	if (x >= xMin && x <= xMax &&
			y >= yMin && y <= yMax)
		return true;
	else
		return false;
}

void Rectangle::center()
{
	double x = (xMin+xMax)/2.0;
	double y = (yMin+yMax)/2.0;

	xMin = xMin-x;
	xMax = xMax-x;
	yMin = yMin-y;
	yMax = yMax-y;
}

void Rectangle::translate(double x, double y)
{
	xMin = xMin+x;
	xMax = xMax+x;
	yMin = yMin+y;
	yMax = yMax+y;
}

bool Rectangle::intersectBool(const Rectangle& r)
{
	if (xMax < r.xMin || r.xMax < xMin) // left, right
		return false;
	if (yMax < r.yMin || r.yMax < yMin)  // bottom, top
		return false;
	return true;
}

bool Rectangle::intersect(const Rectangle& r, Rectangle& ret)
{
	// check if intersect
	if (!intersectBool(r))
		return false;
	// this and ret can be the same rectangle
	ret.xMin = max(xMin, r.xMin);
	ret.xMax = min(xMax,r.xMax);
	ret.yMin = max(yMin, r.yMin);
	ret.yMax = min(yMax,r.yMax);
	ret.width = ret.xMax - ret.xMin;
	ret.height = ret.yMax - ret.yMin;
	return true;
}

void Rectangle::print()
{
	printf("%6.2f %6.2f %6.2f %6.2f\n",xMin,xMax,yMin,yMax);
}

void Rectangle::getMidPoint(double& x, double& y)
{
	x = (xMin+xMax)/2.0;
	y = (yMin+yMax)/2.0;
}

void Rectangle::getClosestPointToOrigin(double& x, double &y)
{
	double minDist = 9999999;
	double origin[2] = {0,0};
	x = 0;
	y = 0;
	// do for each side of the rectangle
	double p[4][2] = { {xMin,yMax},{xMax,yMax},{xMax,yMin},{xMin,yMin} };
	for (int i = 0; i < 4; i++)
	{
		double a[2] = {p[i][0],p[i][1]};
		double b[2] = {p[(i+1)%4][0],p[(i+1)%4][1]};

		if (fabs(a[0]-b[0]) < 0.000001 && fabs(a[1]-b[1]) < 0.000001)
		{
			double dist = distance(a,origin,2);
			if (dist < minDist)
			{
				minDist = dist;
				x = a[0];
				y = a[1];
			}
			continue;
		}
		double d[2] = {b[0]-a[0],b[1]-a[1]};
		double n = norm(d,2);
		if (n < 0.000001 || isnan(n))
		{
			double dist = distance(a,origin,2);
			if (dist < minDist)
			{
				minDist = dist;
				x = a[0];
				y = a[1];
			}
			continue;
		}
		double t = -dot(d,a,2)/n; // projection of origin to b-a
		if (t < 0)
		{
			double dist = distance(a,origin,2);
			if (dist < minDist)
			{
				minDist = dist;
				x = a[0];
				y = a[1];
			}
		}
		else if (t > 1)
		{
			double dist = distance(b,origin,2);
			if (dist < minDist)
			{
				minDist = dist;
				x = b[0];
				y = b[1];
			}
		}
		else
		{
			double q[2] = {a[0]+t*d[0], a[1]+t*d[1]};
			double dist = distance(q,origin,2);
			if (dist < minDist)
			{
				minDist = dist;
				x = q[0];
				y = q[1];
			}
		}
	}
}

void Rectangle::getMaximumClique(const vector<Rectangle*>& input, Rectangle& output,
		int& maxCliqueSize,vector<Rectangle*>& recsMaxClique)
{
	maxCliqueSize = 0; // size of current max clique
	ActivePoint cliquePtY;
	Event cliquePtX;
	// cliquePtX and cliquePtY define a point inside the max clique
	vector<Event> events;
	for (vector<Rectangle*>::const_iterator it = input.begin(); it != input.end(); ++it)
	{
		Rectangle* r = *it;
		events.push_back(Event(r,true));
		events.push_back(Event(r,false));
	}
	sort(events.begin(),events.end(),EventSortByX()); // by x coords of edge

	// sweep events from left to right
	list<ActivePoint> activePts; // the horizontal edges
	for (vector<Event>::iterator it = events.begin(); it != events.end(); ++it)
	{
		Event& e = *it;
		Rectangle* r = e.r;
		if (e.left)
		{
			// check other rectangles in activePts to see if r contains it;
			// update size of clique and max size clique
			// also check if r is contained by rectangles in activePts
			ActivePoint apBottom(r,false);
			ActivePoint apTop(r,true);
			for (list<ActivePoint>::iterator itAP = activePts.begin(); itAP != activePts.end(); ++itAP)
			{
				ActivePoint& ap = *itAP;

				if (ap.r == r)
					continue;  // ignore interection to self
				// check if activePts rects is contained by r
				if (r->yMin <= ap.y && ap.y <= r->yMax)
				{
					ap.cliqueSize++;
					if (ap.cliqueSize > maxCliqueSize)
					{
						maxCliqueSize = ap.cliqueSize;
						cliquePtY = ap;
						cliquePtX = e;
					}
				}
				// check r's top/bottom edges if they are contained by activePts rectangles
				// to avoid double counting only check r when encounter top edge of ap
				if (ap.top && ap.r->yMin <= apBottom.y && apBottom.y <= ap.r->yMax)
				{
					apBottom.cliqueSize++;
					if (apBottom.cliqueSize > maxCliqueSize)
					{
						maxCliqueSize = apBottom.cliqueSize;
						cliquePtY = apBottom;
						cliquePtX = e;
					}
				}
				if (ap.top && ap.r->yMin <= apTop.y && apTop.y <= ap.r->yMax)
				{
					apTop.cliqueSize++;
					if (apTop.cliqueSize > maxCliqueSize)
					{
						maxCliqueSize = apTop.cliqueSize;
						cliquePtY = apTop;
						cliquePtX = e;
					}
				}
			}
			// insert into active points
			activePts.push_back(apBottom);
			activePts.push_back(apTop);
		}
		else
		{
			// remove from active points
			list<ActivePoint>::iterator it = activePts.begin();
			while (it != activePts.end())
			{
				ActivePoint& ap = *it;
				if (ap.r == r)
					it = activePts.erase(it);
				else
				{
					// update size of cliques if are contained by r
					if (r->yMin <= ap.y && ap.y <= r->yMax)
					{
						ap.cliqueSize--; // current clique size of ap shrinks
					}
					++it;
				}
			}
		}
	} // end for each event

	// only have cliques of size 1
	if (cliquePtX.r == NULL)
	{
		// pick one rectangle; pick randomly
		int pick = rand()%input.size();
		Rectangle* r = input[pick];
		cliquePtY = ActivePoint(r,true);
		cliquePtX = Event(r,true);
		maxCliqueSize = 1;
	}

	// get the intersecting rectangles
	int count = 0;
	for (vector<Rectangle*>::const_iterator it = input.begin(); it != input.end(); ++it)
	{
		Rectangle* r = *it;
		if (r->yMin <= cliquePtY.y && cliquePtY.y <= r->yMax &&
			r->xMin <= cliquePtX.x && cliquePtX.x <= r->xMax)
		{
			recsMaxClique.push_back(r);
			if (count > 0)
			{
				if (!output.intersect(*r,output))
					printf("ERROR: No intersection. Bug in clique finding\n");
			}
			else
				output = *r;

			count++;
			// r->print();
		}
	}
	if (maxCliqueSize != count)
	{
		printf("ERROR: clique size mismatch. Bug in clique finding\n");
	}
	if (count == 0)
	{
		printf("ERROR: output is null. Bug in clique finding\n");
	}
}

//void Rectangle::calibrate(const vector< tr1::array<double,2> >& ref, const vector< tr1::array<double,2> >& input,
//		   double matchTolN, double matchTolHN, double searchN, double searchHN,
//		   double& offsetN, double& offsetHN, bool debug)
//{
//	vector<Rectangle*> rectangleMatches; // input rectangles for finding max clique
//	offsetN = 0;
//	offsetHN = 0;
//	for (vector< tr1::array<double,2> >::const_iterator it = input.begin(); it != input.end(); ++it) // NOEs
//	{
//		const tr1::array<double,2>& inputCS = *it;
//		Rectangle cR(-searchN,searchN,-searchHN,searchHN);
//		double n = inputCS[0];
//		double hn = inputCS[1];
//		for (vector< tr1::array<double,2> >::const_iterator itR = ref.begin(); itR != ref.end(); ++itR) // chem shift assign
//		{
//			const tr1::array<double,2>& refCS = *itR;
//			double nR = refCS[0];
//			double hnR = refCS[1];
//			double dN = fabs(n-nR);
//			double dHN = fabs(hn-hnR);
//			if (dN <= matchTolN+searchN && dHN <= matchTolHN+searchHN)
//			{
//				Rectangle rR(nR-n-matchTolN, nR-n+matchTolN,
//					hnR-hn-matchTolHN, hnR-hn+matchTolHN); // -n,-hn  since rectangle is centered at origin
//			                                // matchTol is furthest any peak can be from it
//				Rectangle* intersect = new Rectangle;
//				if (cR.intersect(rR,*intersect))
//					rectangleMatches.push_back(intersect);
//				else
//					delete intersect;
//			}
//		}
//	}
//	if (rectangleMatches.size() > 0)
//	{
//		Rectangle output;
//		int maxCliqueSize = 0;
//		vector<Rectangle*> recsMaxClique; // rectangles in maximum clique
//		getMaximumClique(rectangleMatches, output, maxCliqueSize, recsMaxClique);
//		if (maxCliqueSize > 0)
//		{
//			if (!output.contains(0,0))
//			{
//				output.getMidPoint(offsetN,offsetHN);
//				if (debug)
//					printf("calibration with %6.2f %6.2f\n",offsetN,offsetHN);
//			}
//			else
//			{
//				if (debug)
//					printf("calibration not necessary\n");
//			}
//		}
//		else
//		{
//			if (debug)
//				printf("Unable to calibrate. No clique\n");
//		}
//	}
//	else
//	{
//		if (debug)
//			printf("Unable to calibrate. Peak lists do not match\n");
//	}
//
//	for (vector<Rectangle*>::iterator it = rectangleMatches.begin(); it != rectangleMatches.end(); ++it)
//	{
//		delete *it;
//	}
//}

/*
 * Use this version instead of the other versions of calibrate
 * ref should be the assigned chemical shifts
 * input should be the NOEs
 */
void Rectangle::calibrate(const vector< tr1::array<double,2> >& ref, const vector< tr1::array<double,2> >& input,
		   double matchTolN, double matchTolHN, double searchN, double searchHN,
		   double& offsetN, double& offsetHN, bool debug)
{
	vector<Rectangle*> rectangleMatches;
	offsetN = 0;
	offsetHN = 0;
	for (vector< tr1::array<double,2> >::const_iterator it = ref.begin(); it != ref.end(); ++it) // chem shift assign
	{
		const tr1::array<double,2>& refCS = *it;
		Rectangle refRec(-searchN,searchN,-searchHN,searchHN);
		double n = refCS[0];
		double hn = refCS[1];
		for (vector< tr1::array<double,2> >::const_iterator itI = input.begin(); itI != input.end(); ++itI) // NOEs
		{
			const tr1::array<double,2>& inputCS = *itI;
			double nI = inputCS[0];
			double hnI = inputCS[1];
			double dN = fabs(n-nI);
			double dHN = fabs(hn-hnI);
			if (dN <= matchTolN+searchN && dHN <= matchTolHN+searchHN)
			{
				Rectangle recInput(nI-n-matchTolN, nI-n+matchTolN,
					hnI-hn-matchTolHN, hnI-hn+matchTolHN); // -n,-hn  since rectangle is centered at origin
			                                // matchTol is furthest any peak can be from it
				Rectangle* intersect = new Rectangle;
				if (refRec.intersect(recInput,*intersect))
				{
					intersect->dx = n-nI; // offset to add to NOEs; store in dx rectangle annotation
					intersect->dy = hn-hnI;
					rectangleMatches.push_back(intersect);
				}
				else
					delete intersect;
			}
		}
	}
	if (rectangleMatches.size() > 0)
	{
		Rectangle output;
		int maxCliqueSize = 0;
		vector<Rectangle*> recsMaxClique;
		getMaximumClique(rectangleMatches, output, maxCliqueSize,recsMaxClique);
		if (maxCliqueSize > 0)
		{
			// filter outliers; take first quartile value
			sort(recsMaxClique.begin(),recsMaxClique.end(),RectangleSortByDx());
			double midpt = recsMaxClique[maxCliqueSize/2]->dx;
			if (midpt > 0)
				offsetN = recsMaxClique[maxCliqueSize/4]->dx;
			else
				offsetN = recsMaxClique[3*maxCliqueSize/4]->dx;
			sort(recsMaxClique.begin(),recsMaxClique.end(),RectangleSortByDy());
			midpt = recsMaxClique[maxCliqueSize/2]->dy;
			if (midpt > 0)
				offsetHN = recsMaxClique[maxCliqueSize/4]->dy;
			else
				offsetHN = recsMaxClique[3*maxCliqueSize/4]->dy;
			if (debug)
				printf("calibration with %6.3f %6.3f\n",offsetN,offsetHN);
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
			printf("Unable to calibrate. Peak lists do not match\n");
	}

	for (vector<Rectangle*>::iterator it = rectangleMatches.begin(); it != rectangleMatches.end(); ++it)
	{
		delete *it;
	}
}

bool Rectangle::calibrateUnique(const vector< tr1::array<double,2> >& ref, const vector< tr1::array<double,2> >& input,
		double searchN, double searchHN,
		   double& offsetN, double& offsetHN, bool debug)
{
	// look for unique match within matchTolN, matchTolHN
	int numUniqueMatches=0;
	offsetN = 0;
	offsetHN = 0;
	for (vector< tr1::array<double,2> >::const_iterator it = input.begin(); it != input.end(); ++it)
	{
		const tr1::array<double,2>& inputCS = *it;
		double nI=inputCS[0];
		double hnI=inputCS[1];
		double dN = 0;
		double dHN = 0;
		int numMatches = 0;
		for (vector< tr1::array<double,2> >::const_iterator itR = ref.begin(); itR != ref.end(); ++itR)
		{
			const tr1::array<double,2>& refCS = *itR;
			double nR = refCS[0];
			double hnR = refCS[1];
			if (fabs(nR-nI) <= searchN && fabs(hnR-hnI) <= searchHN)
			{
				dN = nR-nI;
				dHN = hnR-hnI;
				numMatches++;
				if (numMatches > 1)
					break;
			}
		}
		if (numMatches == 1)
		{
			offsetN += dN;
			offsetHN += dHN;
			numUniqueMatches++;
		}
	}
	if (numUniqueMatches > 0)
	{
		offsetN = offsetN/double(numUniqueMatches);
		offsetHN = offsetHN/double(numUniqueMatches);
		if (debug)
			printf("CalibrateUnique[%d]: %f %f\n",numUniqueMatches,offsetN,offsetHN);
		return true;
	}
	else
	{
		printf("WARNING: No unique matches for calibration\n");
		return false;
	}
}

///*
// * Expected output:
// *  13.00  16.00   3.00  15.00
// *  5.00  15.00  10.00  14.00
// *  15.00  15.00  14.00  16.00
// *  15.00  15.00  14.00  16.00
// *  clique size: 4
// *  15.00  15.00  14.00  14.00
// */
//void Rectangle::test1()
//{
//	Rectangle r1(2,3,1,13);
//	Rectangle r2(1,14,2,5);
//	Rectangle r3(13,16,3,15);
//	Rectangle r4(9,10,4,11);
//	Rectangle r5(4,8,6,12);
//	Rectangle r6(6,7,8,13);
//	Rectangle r7(11,12,7,9);
//	Rectangle r8(5,15,10,14);
//	Rectangle r9(15,15,14,16);
//	Rectangle r10(15,15,14,16);
//
//	vector<Rectangle*> input;
//	input.reserve(10);
//	input.push_back(&r1);
//	input.push_back(&r2);
//	input.push_back(&r3);
//	input.push_back(&r4);
//	input.push_back(&r5);
//	input.push_back(&r6);
//	input.push_back(&r7);
//	input.push_back(&r8);
//	input.push_back(&r9);
//	input.push_back(&r10);
//
//	Rectangle ret;
//	int cliqueSize = 0;
//	vector<Rectangle*> recsMaxClique;
//	Rectangle::getMaximumClique(input,ret,cliqueSize,recsMaxClique);
//	printf("clique size: %d\n",cliqueSize);
//	ret.print();
//}
//
//void Rectangle::test2()
//{
//	ifstream file;
//	string line;
//	file.open("/home/e4k2/Documents/nmr_data/CAM13/15624.list");
//	if (!file.is_open())
//	{
//		printf("Cannot open file 15624.list\n");
//		return;
//	}
//	vector< tr1::array<double,2> > ref;
//	double searchN = 1.0;
//	double searchHN = 0.1;
//	double matchTolN = 0.5;
//	double matchTolHN = 0.05;
//
//	while (file.good())
//	{
//		getline(file,line);
//		trim(line);
//		if (line.empty())
//			continue;
//
//		stringstream tok(line);
//		string temp;
//		tok >> temp; // res num
//		tok >> temp; // n shift
//		const double n = atof(temp.c_str());
//		tok >> temp;
//		const double hn = atof(temp.c_str());
//		tr1::array<double,2> a;
//		a[0] = n;
//		a[1] = hn;
//		ref.push_back(a);
//	}
//
//	std::random_device rd;
//	std::mt19937 gen(rd());
//	std::uniform_real_distribution<> dis(0, 1);
//	double globalNoiseN = 0.9;
//	double globalNoiseHN = -0.09;
//	double localNoiseN = 0; // 0.15;
//	double localNoiseHN = 0; // 0.015;
//	const int NUMTESTS = 100;
//	double avgError = 0;
//	double avgNumMissed = 0;
//	int maxNumMissed = 0;
//	for (int i = 0; i < NUMTESTS; i++)
//	{
//		vector< tr1::array<double,2> > input;
//		input.reserve(ref.size());
//
//		for (vector< tr1::array<double,2> >::iterator it = ref.begin(); it != ref.end(); ++it)
//		{
//			tr1::array<double,2>& a = *it;
//			double n=a[0];
//			double hn=a[1];
//			double dn=0;
//			double dhn=0;
//
//			if (dis(gen) > 0.5)
//				dn = dis(gen)*localNoiseN;
//			else
//				dn = -dis(gen)*localNoiseN;
//
//			n = n+dn+globalNoiseN;
//
//			if (dis(gen) > 0.5)
//				dhn = dis(gen)*localNoiseHN;
//			else
//				dhn = -dis(gen)*localNoiseHN;
//
//			hn = hn+dhn+globalNoiseHN;
//
//			tr1::array<double,2> inA;
//			inA[0] = n;
//			inA[1] = hn;
//			input.push_back(inA);
//		}
//
//		double offsetN = 0;
//		double offsetHN = 0;
//		calibrate(ref, input,matchTolN,matchTolHN,searchN,searchHN, offsetN, offsetHN,false);
//		// calibrateUnique(ref, input, searchN+matchTolN, searchHN+matchTolHN, offsetN, offsetHN, false);
////		for (vector< array<double,2> >::iterator it = input.begin(); it != input.end(); ++it)
////		{
////			array<double,2>& a = *it;
////			a[0] = a[0]+offsetN;
////			a[1] = a[1]+offsetHN;
////		}
////		if (!calibrateUnique(ref, input, 10.0, 1.0, offsetN, offsetHN, true))
////			printf("No unique matches\n");
//
//		// compare
//		vector< tr1::array<double,2> >::iterator itI = input.begin();
//		int numMissed = 0;
//		double error = 0;
//		for (vector< tr1::array<double,2> >::iterator it = ref.begin(); it != ref.end(); ++it)
//		{
//			tr1::array<double,2>& refA = *it;
//			tr1::array<double,2>& inA = *itI;
//			double n = refA[0];
//			double nI = inA[0]+offsetN;
//			double dN = fabs(n-nI);
//			double hn = refA[1];
//			double hnI = inA[1]+offsetHN;
//			double dHN = fabs(hn-hnI);
//			error += dN+10.0*dHN;
//			if (dN > matchTolN || dHN > matchTolHN)
//				numMissed++;
//			++itI;
//		}
//		error = error/double(ref.size());
//		avgError += error;
//		avgNumMissed += numMissed;
//		if (numMissed > maxNumMissed)
//			maxNumMissed = numMissed;
//	}
//	avgError = avgError/double(NUMTESTS);
//	avgNumMissed = avgNumMissed/double(NUMTESTS);
//	printf("avgError: %f\n",avgError);
//	printf("avgNumMissed: %f\n",avgNumMissed);
//	printf("maxNumMissed: %d\n",maxNumMissed);
//	file.close();
//}
//
//void Rectangle::test3()
//{
//	ifstream file;
//	string line;
//	string filename // ("temp.bmrb");
//	      ("2KNR.bmrb");
//	    // ("/home/e4k2/Documents/nmr_data/ubiquitin/5387.shift");
//		//	 ("/home/e4k2/Documents/nmr_data/CAM13/15624.bmrb");
//	file.open(filename.c_str());
//	if (!file.is_open())
//	{
//		printf("Cannot open file 15624.bmrb\n");
//		return;
//	}
//
//	double searchN = 1.0;
//	double searchHN = 0.1;
//	double matchTolN = 0.3;
//	double matchTolHN = 0.035;
//
//	tr1::unordered_map<int, tr1::unordered_map<string, list<double> > > refMap;
//	// [resnum][carbon_atomname, C and H chem shifts] (first entry is C)
//	double minC = 999999;
//	double maxC = 0;
//	double minH = 999999;
//	double maxH = 0;
//	while (file.good())
//	{
//		getline(file,line);
//		trim(line);
//		if (line.empty())
//			continue;
//
//		stringstream tok(line);
//		string temp;
//		tok >> temp; // serial
//		tok >> temp; // resnum
//		int resnum = atoi(temp.c_str());
//		tok >> temp; // resname
//		string atomname;
//		tok >> atomname;
//		tok >> temp; // atom type
//		if (temp[0] == 'C' && atomname.size() > 1)
//		{
//			tok >> temp; // c shift
//			const double c = atof(temp.c_str());
//			if (c < minC)
//				minC = c;
//			if (c > maxC)
//				maxC = c;
//
//			tr1::unordered_map<int, tr1::unordered_map<string, list<double> > >::iterator it = refMap.find(resnum);
//			if (it != refMap.end())
//			{
//				tr1::unordered_map<string, list<double> >& cMap = it->second;
//				list<double> shifts;
//				shifts.push_back(c);
//				cMap[atomname] = shifts;
//			}
//			else
//			{
//				tr1::unordered_map<string, list<double> > cMap;
//				list<double> shifts;
//				shifts.push_back(c);
//				cMap[atomname] = shifts;
//				refMap[resnum] = cMap;
//			}
//		}
//	}
//	file.close();
//	file.open(filename.c_str());
//	while (file.good())
//	{
//		getline(file,line);
//		trim(line);
//		if (line.empty())
//			continue;
//		stringstream tok(line);
//		string temp;
//		tok >> temp; // serial
//		tok >> temp; // resnum
//		int resnum = atoi(temp.c_str());
//		tok >> temp; // resname
//		string atomname;
//		tok >> atomname;
//		tok >> temp; // atom type
//		if (temp[0] == 'H' && atomname.size() > 1)
//		{
//			tok >> temp; // h shift
//			const double h = atof(temp.c_str());
//			if (h < minH)
//				minH = h;
//			if (h > maxH)
//				maxH = h;
//			tr1::unordered_map<int, tr1::unordered_map<string, list<double> > >::iterator it = refMap.find(resnum);
//			if (it != refMap.end())
//			{
//				tr1::unordered_map<string, list<double> >& cMap = it->second;
//				string cname = "C"+atomname.substr(1);
//				tr1::unordered_map<string, list<double> >::iterator itC = cMap.find(cname);
//				if (itC == cMap.end() && cname.size() > 2)
//				{
//					int len = atomname.size();
//					cname = "C"+atomname.substr(1,len-2);
//					itC = cMap.find(cname);
//				}
//				if (itC != cMap.end())
//				{
//					list<double>& shifts = itC->second;
//					shifts.push_back(h);
//				}
//				// else no c shift for this proton
//			}
//			// else no c shifts for this resnum
//		}
//	}
//	vector< tr1::array<double,2> > refTemp;
//
//	for (tr1::unordered_map<int, tr1::unordered_map<string, list<double> > >::iterator it = refMap.begin(); it != refMap.end(); ++it)
//	{
//		tr1::unordered_map<string, list<double> >& cMap = it->second;
//		for (tr1::unordered_map<string, list<double> >::iterator itC = cMap.begin(); itC != cMap.end(); ++itC)
//		{
//			list<double>& shifts = itC->second;
//			double c = 0;
//			int i = 0;
//			for (double s : shifts)
//			{
//				if (i == 0)
//				{
//					c = s;
//				}
//				else
//				{
//					tr1::array<double,2> a;
//					a[0] = c;
//					a[1] = s;
//					refTemp.push_back(a);
//				}
//				i++;
//			}
//		}
//	}
//	double rangeC = maxC-minC;
//	double rangeH = maxH-minH;
//	double globalNoiseN = 0.8; // 0.9;
//	double globalNoiseHN = 0.08; // -0.09;
//	double localNoiseN = 0.25; // 0.3; // 0.2; // 0.15;
//	double localNoiseHN = 0.03; // 0.03; // 0.02; // 0.015;
//	const int NUMTESTS = 100; // 100; // 1000;
//	double avgError = 0;
//	double avgNumMissed = 0;
//	int maxNumMissed = 0;
//
//	double avgError2 = 0;
//	double avgNumMissed2 = 0;
//	int maxNumMissed2 = 0;
//
//	double avgMissedDev = 0;
//	double avgMissedDev2 = 0;
//
//	int numTimesMissed = 0;
//	int numTimesMissed2 = 0;
//
//	int sizeCut = refTemp.size()/4;
//	int size = refTemp.size()-sizeCut;
//	for (int i = 0; i < NUMTESTS; i++)
//	{
//		random_shuffle(refTemp.begin(),refTemp.end());
//
//		// missing peaks
//		vector< tr1::array<double,2> > ref(refTemp.begin(),refTemp.begin()+size);
//
//		std::random_device rd;
//		std::mt19937 gen(rd());
//		std::uniform_real_distribution<> dis(0, 1);
//
//		int sizeAdd = ref.size()/10; // for noise peaks
//		vector< tr1::array<double,2> > input;
//		input.reserve(ref.size()+sizeAdd);
//
//		for (vector< tr1::array<double,2> >::iterator it = ref.begin(); it != ref.end(); ++it)
//		{
//			tr1::array<double,2>& a = *it;
//			double n=a[0];
//			double hn=a[1];
//
//			double dn=0;
//			double dhn=0;
//			if (dis(gen) > 0.5)
//				dn = dis(gen)*localNoiseN;
//			else
//				dn = -dis(gen)*localNoiseN;
//
//			n = n+dn+globalNoiseN;
//
//			if (dis(gen) > 0.5)
//				dhn = dis(gen)*localNoiseHN;
//			else
//				dhn = -dis(gen)*localNoiseHN;
//
//			hn = hn+dhn+globalNoiseHN;
//
//			tr1::array<double,2> inA;
//			inA[0] = n;
//			inA[1] = hn;
//			input.push_back(inA);
//		}
//		// noise peaks
//		for (int i = 0; i < sizeAdd; i++)
//		{
//			double c = dis(gen)*rangeC+minC;
//			double h = dis(gen)*rangeH+minH;
//			tr1::array<double,2> inA;
//			inA[0] = c;
//			inA[1] = h;
//			input.push_back(inA);
//		}
//
//		double offsetN = 0;
//		double offsetHN = 0;
//		// calibrate(ref, input,matchTolN-0.000001,matchTolHN-0.000001,searchN,searchHN, offsetN, offsetHN,false);
//		calibrate(ref, input,matchTolN,matchTolHN,searchN,searchHN, offsetN, offsetHN,false);
//		//calibrate(ref, input,0.25,0.035,searchN,searchHN, offsetN, offsetHN,false);
//		double offsetN2 = 0;
//		double offsetHN2 = 0;
//		//calibrateUnique(ref, input, searchN, searchHN, offsetN2, offsetHN2, false);
//		calibrateUnique(ref, input, searchN+matchTolN, searchHN+matchTolHN, offsetN2, offsetHN2, false);
////		for (vector< array<double,2> >::iterator it = input.begin(); it != input.end(); ++it)
////		{
////			array<double,2>& a = *it;
////			a[0] = a[0]+offsetN;
////			a[1] = a[1]+offsetHN;
////		}
////		if (!calibrateUnique(ref, input, 10.0, 1.0, offsetN, offsetHN, true))
////			printf("No unique matches\n");
//
//		// compare
//		vector< tr1::array<double,2> >::iterator itI = input.begin();
//		int numMissed = 0;
//		double error = 0;
//		int numMissed2 = 0;
//		double error2 = 0;
//		double avgMD = 0;
//		double avgMD2 = 0;
//		for (vector< tr1::array<double,2> >::iterator it = ref.begin(); it != ref.end(); ++it)
//		{
//			tr1::array<double,2>& refA = *it;
//			tr1::array<double,2>& inA = *itI;
//			double n = refA[0];
//			double nI = inA[0]+offsetN;
//			double dN = fabs(n-nI);
//			double hn = refA[1];
//			double hnI = inA[1]+offsetHN;
//			double dHN = fabs(hn-hnI);
//			error += dN+10.0*dHN;
//			if (dN > matchTolN || dHN > matchTolHN)
//			{
//				numMissed++;
//				if (dN > matchTolN)
//					avgMD += (dN-matchTolN);
//				if (dHN > matchTolHN)
//					avgMD += 10.0*(dHN-matchTolHN);
//			}
//
//			nI = inA[0]+offsetN2;
//			dN = fabs(n-nI);
//			hnI = inA[1]+offsetHN2;
//			dHN = fabs(hn-hnI);
//		    error2 += dN+10.0*dHN;
//			if (dN > matchTolN || dHN > matchTolHN)
//			{
//				numMissed2++;
//				if (dN > matchTolN)
//					avgMD2 += (dN-matchTolN);
//				if (dHN > matchTolHN)
//					avgMD2 += 10.0*(dHN-matchTolHN);
//			}
//
//			++itI;
//		}
//		error = error/double(ref.size());
//		avgError += error;
//		avgNumMissed += numMissed;
//		if (numMissed > maxNumMissed)
//			maxNumMissed = numMissed;
//
//		error2 = error2/double(ref.size());
//		avgError2 += error2;
//		avgNumMissed2 += numMissed2;
//		if (numMissed2 > maxNumMissed2)
//			maxNumMissed2 = numMissed2;
//
//		if (numMissed > 0)
//		{
//			avgMissedDev += avgMD/double(numMissed);
//			numTimesMissed++;
//		}
//
//		if (numMissed2 > 0)
//		{
//			avgMissedDev2 += avgMD2/double(numMissed2);
//			numTimesMissed2++;
//		}
//	}
//	avgError = avgError/double(NUMTESTS);
//	avgNumMissed = avgNumMissed/double(NUMTESTS);
//	printf("avgError: %f\n",avgError);
//	printf("avgNumMissed: %f\n",avgNumMissed);
//	printf("maxNumMissed: %d out of %d\n",maxNumMissed,size);
//	if (numTimesMissed > 0)
//	{
//		avgMissedDev = avgMissedDev/double(numTimesMissed);
//		printf("avgMissedDev: %f\n",avgMissedDev);
//	}
//	else
//		printf("avgMissedDev: 0\n");
//
//	avgError2 = avgError2/double(NUMTESTS);
//	avgNumMissed2 = avgNumMissed2/double(NUMTESTS);
//	printf("avgError2: %f\n",avgError2);
//	printf("avgNumMissed2: %f\n",avgNumMissed2);
//	printf("maxNumMissed2: %d out of %d\n",maxNumMissed2,size);
//	if (numTimesMissed2 > 0)
//	{
//		avgMissedDev2 = avgMissedDev2/double(numTimesMissed2);
//		printf("avgMissedDev2: %f\n",avgMissedDev2);
//	}
//	else
//		printf("avgMissedDev2: 0\n");
//	file.close();
//}
//

//void Rectangle::test4()
//{
//	ifstream file;
//	string line;
//	string filename // ("temp.bmrb");
//	     // ("2KNR.bmrb");
//	     ("/home/e4k2/Documents/nmr_data/ubiquitin/5387.shift");
//		//	 ("/home/e4k2/Documents/nmr_data/CAM13/15624.bmrb");
//	file.open(filename.c_str());
//	if (!file.is_open())
//	{
//		printf("Cannot open file 15624.bmrb\n");
//		return;
//	}
//
//	double searchN = 1.0;
//	double searchHN = 0.1;
//	double matchTolN = 0.3;
//	double matchTolHN = 0.04;
//
//	tr1::unordered_map<int, tr1::unordered_map<string, list<double> > > refMap;
//	// [resnum][carbon_atomname, C and H chem shifts] (first entry is C)
//	double minC = 999999;
//	double maxC = 0;
//	double minH = 999999;
//	double maxH = 0;
//	while (file.good())
//	{
//		getline(file,line);
//		trim(line);
//		if (line.empty())
//			continue;
//
//		stringstream tok(line);
//		string temp;
//		tok >> temp; // serial
//		tok >> temp; // resnum
//		int resnum = atoi(temp.c_str());
//		tok >> temp; // resname
//		string atomname;
//		tok >> atomname;
//		tok >> temp; // atom type
//		if (temp[0] == 'C' && atomname.size() > 1)
//		{
//			tok >> temp; // c shift
//			const double c = atof(temp.c_str());
//			if (c < minC)
//				minC = c;
//			if (c > maxC)
//				maxC = c;
//
//			tr1::unordered_map<int, tr1::unordered_map<string, list<double> > >::iterator it = refMap.find(resnum);
//			if (it != refMap.end())
//			{
//				tr1::unordered_map<string, list<double> >& cMap = it->second;
//				list<double> shifts;
//				shifts.push_back(c);
//				cMap[atomname] = shifts;
//			}
//			else
//			{
//				tr1::unordered_map<string, list<double> > cMap;
//				list<double> shifts;
//				shifts.push_back(c);
//				cMap[atomname] = shifts;
//				refMap[resnum] = cMap;
//			}
//		}
//	}
//	file.close();
//	file.open(filename.c_str());
//	while (file.good())
//	{
//		getline(file,line);
//		trim(line);
//		if (line.empty())
//			continue;
//		stringstream tok(line);
//		string temp;
//		tok >> temp; // serial
//		tok >> temp; // resnum
//		int resnum = atoi(temp.c_str());
//		tok >> temp; // resname
//		string atomname;
//		tok >> atomname;
//		tok >> temp; // atom type
//		if (temp[0] == 'H' && atomname.size() > 1)
//		{
//			tok >> temp; // h shift
//			const double h = atof(temp.c_str());
//			if (h < minH)
//				minH = h;
//			if (h > maxH)
//				maxH = h;
//			tr1::unordered_map<int, tr1::unordered_map<string, list<double> > >::iterator it = refMap.find(resnum);
//			if (it != refMap.end())
//			{
//				tr1::unordered_map<string, list<double> >& cMap = it->second;
//				string cname = "C"+atomname.substr(1);
//				tr1::unordered_map<string, list<double> >::iterator itC = cMap.find(cname);
//				if (itC == cMap.end() && cname.size() > 2)
//				{
//					int len = atomname.size();
//					cname = "C"+atomname.substr(1,len-2);
//					itC = cMap.find(cname);
//				}
//				if (itC != cMap.end())
//				{
//					list<double>& shifts = itC->second;
//					shifts.push_back(h);
//				}
//				// else no c shift for this proton
//			}
//			// else no c shifts for this resnum
//		}
//	}
//	vector< tr1::array<double,2> > refTemp;
//
//	for (tr1::unordered_map<int, tr1::unordered_map<string, list<double> > >::iterator it = refMap.begin(); it != refMap.end(); ++it)
//	{
//		tr1::unordered_map<string, list<double> >& cMap = it->second;
//		for (tr1::unordered_map<string, list<double> >::iterator itC = cMap.begin(); itC != cMap.end(); ++itC)
//		{
//			list<double>& shifts = itC->second;
//			double c = 0;
//			int i = 0;
//			for (double s : shifts)
//			{
//				if (i == 0)
//				{
//					c = s;
//				}
//				else
//				{
//					tr1::array<double,2> a;
//					a[0] = c;
//					a[1] = s;
//					refTemp.push_back(a);
//				}
//				i++;
//			}
//		}
//	}
//	double rangeC = maxC-minC;
//	double rangeH = maxH-minH;
//	double globalNoiseN = 0.8; // 0.9;
//	double globalNoiseHN = 0.08; // -0.09;
//	double localNoiseN = 0; // 0.29; // 0.25; // 0.3; // 0.2; // 0.15;
//	double localNoiseHN = 0; // 0.039; // 0.035; // 0.03; // 0.02; // 0.015;
//	int NUMCOPIES = 10; // num NOEs belonging to same C,HC pair
//	const int NUMTESTS = 20; // 100; // 1000;
//	double avgError = 0;
//	double avgNumMissed = 0;
//	int maxNumMissed = 0;
//	int numTimesMissed = 0;
//	double avgMissedDev = 0;
//
//	int sizeCut = refTemp.size()/4;
//	int size = refTemp.size()-sizeCut;
//	std::random_device rd;
//	std::mt19937 gen(rd());
//	std::uniform_real_distribution<> dis(0, 1);
//	std::default_random_engine generator;
//	std::normal_distribution<double> ndis(0,1); // mean 0, std 1
//
//	for (int i = 0; i < NUMTESTS; i++)
//	{
//		random_shuffle(refTemp.begin(),refTemp.end());
//
//		// missing peaks
//		vector< tr1::array<double,2> > ref(refTemp.begin(),refTemp.begin()+size);
//
//		int sizeAdd = ref.size()/10; // for noise peaks
//		vector< tr1::array<double,2> > input;
//		input.reserve(NUMCOPIES*ref.size()+sizeAdd);
//
//		for (vector< tr1::array<double,2> >::iterator it = ref.begin(); it != ref.end(); ++it)
//		{
//			tr1::array<double,2>& a = *it;
//			for (int i = 0; i < NUMCOPIES; i++)
//			{
//				double n = a[0];
//				double hn = a[1];
//				double dn=0;
//				double dhn=0;
//				double nn = dis(gen); // ndis(generator)
//				double hnhn = dis(gen);
//				if (nn > 0.5)
//					dn = nn*localNoiseN;
//				else
//					dn = -nn*localNoiseN;
//
//				n = n+dn+globalNoiseN;
//
//				if (hnhn > 0.5)
//					dhn = hnhn*localNoiseHN;
//				else
//					dhn = -hnhn*localNoiseHN;
//
//				hn = hn+dhn+globalNoiseHN;
//
//				tr1::array<double,2> inA;
//				inA[0] = n;
//				inA[1] = hn;
//				input.push_back(inA);
//			}
//		}
//		// noise peaks
//		for (int i = 0; i < sizeAdd; i++)
//		{
//			double c = dis(gen)*rangeC+minC;
//			double h = dis(gen)*rangeH+minH;
//			tr1::array<double,2> inA;
//			inA[0] = c;
//			inA[1] = h;
//			input.push_back(inA);
//		}
//
//		double offsetN = 0;
//		double offsetHN = 0;
//		calibrate(ref, input,matchTolN,matchTolHN,searchN,searchHN, offsetN, offsetHN,true);
//
//		// compare
//		vector< tr1::array<double,2> >::iterator itI = input.begin();
//		int numMissed = 0;
//		double error = 0;
//		double avgMD = 0;
//		for (vector< tr1::array<double,2> >::iterator it = ref.begin(); it != ref.end(); ++it)
//		{
//			tr1::array<double,2>& refA = *it;
//			double n = refA[0];
//			double hn = refA[1];
//			for (int i = 0; i < NUMCOPIES; i++)
//			{
//				tr1::array<double,2>& inA = *itI;
//				double nI = inA[0]+offsetN;
//				double dN = fabs(n-nI);
//				double hnI = inA[1]+offsetHN;
//				double dHN = fabs(hn-hnI);
//				error += dN+10.0*dHN;
//				if (dN > matchTolN || dHN > matchTolHN)
//				{
//					numMissed++;
//					if (dN > matchTolN)
//						avgMD += (dN-matchTolN);
//					if (dHN > matchTolHN)
//						avgMD += 10.0*(dHN-matchTolHN);
//				}
//				++itI;
//			}
//		}
//		error = error/double(NUMCOPIES*ref.size());
//		avgError += error;
//		double numMissedPerTest = double(numMissed)/double(NUMCOPIES);
//		avgNumMissed += numMissedPerTest;
//		if (numMissedPerTest > maxNumMissed)
//			maxNumMissed = numMissedPerTest;
//		if (numMissedPerTest > 0)
//		{
//			avgMissedDev += avgMD/double(numMissedPerTest);
//			numTimesMissed++;
//		}
//	}
//	avgError = avgError/double(NUMTESTS);
//	avgNumMissed = avgNumMissed/double(NUMTESTS);
//	printf("avgError: %f\n",avgError);
//	printf("avgNumMissed: %f\n",avgNumMissed);
//	printf("maxNumMissed: %d out of %d\n",maxNumMissed,size);
//	if (numTimesMissed > 0)
//	{
//		avgMissedDev = avgMissedDev/double(numTimesMissed);
//		printf("avgMissedDev: %f\n",avgMissedDev);
//	}
//	else
//		printf("avgMissedDev: 0\n");
//	file.close();
//}
