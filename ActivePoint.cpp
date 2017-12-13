/*
 * ActivePoint.cpp
 *
 *  Created on: 2012-07-16
 *      Author: e4k2
 */
#include <cstdio>
#include <algorithm>
#include "ActivePoint.h"

ActivePoint::ActivePoint() : r(NULL), top(false),y(0),cliqueSize(0)
{
}

ActivePoint::ActivePoint(Rectangle* r, bool top) : r(r), top(top)
{
	if (top)
	{
		y = r->yMax;
	}
	else
	{
		y = r->yMin;
	}
	cliqueSize = 1; // includes itself
}

ActivePoint::ActivePoint(const ActivePoint& a) : r(a.r), top(a.top), y(a.y), cliqueSize(a.cliqueSize)
{
}

ActivePoint::~ActivePoint()
{
}

void swap(ActivePoint& first, ActivePoint& second)
{
	using std::swap;
	swap(first.r,second.r);
	swap(first.top,second.top);
	swap(first.y,second.y);
	swap(first.cliqueSize,second.cliqueSize);
}

//ActivePoint::ActivePoint(ActivePoint&& a) : r(NULL), top(false),y(0),cliqueSize(0)
//{
//	swap(*this,a);
//}

ActivePoint& ActivePoint::operator=(ActivePoint a)
{
	swap(*this,a);
	return *this;
}

void ActivePoint::print()
{
	r->print();
	if (top)
		printf("yMax: %6.2f cliqueSize: %d\n",y,cliqueSize);
	else
		printf("yMin: %6.2f cliqueSize: %d\n",y,cliqueSize);
}
