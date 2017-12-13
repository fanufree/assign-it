/*
 * ActivePoint.h
 *
 *  Created on: 2012-07-16
 *      Author: e4k2
 */

#ifndef ACTIVEPOINT_H_
#define ACTIVEPOINT_H_

#include "Rectangle.h"

/**
 * Known as the status. Represents a horizontal edge of some rectangle
 * Reference:
 * Hiroshi Imai and Takao Asano. Finding the connected components and a maximum clique of an
 * intersection graph of rectangles in the plane. Journal of Algorithms, 4(4):310 – 323, 1983.
 */
class ActivePoint
{
public:
	ActivePoint();
	ActivePoint(Rectangle* r, bool top);
	ActivePoint(const ActivePoint& a);
	virtual ~ActivePoint();
	friend void swap(ActivePoint& first, ActivePoint& second);
	// ActivePoint(ActivePoint&& a);
	ActivePoint& operator=(ActivePoint a);
	void print();

	Rectangle* r;
	bool top; // true if top edge of rectangle; else bottom edge
	double y; // coordinate edge
	int cliqueSize; // to be annotated to hold the num rectangles in clique
};

#endif /* ACTIVEPOINT_H_ */
