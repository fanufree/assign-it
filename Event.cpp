/*
 * Event.cpp
 *
 *  Created on: 2012-07-16
 *      Author: e4k2
 */
#include <cstdio>
#include <algorithm>
#include "Event.h"

Event::Event() : r(NULL), left(false), x(0)
{
}

Event::Event(Rectangle* r, bool left) : r(r), left(left)
{
	if (left)
		x = r->xMin;
	else
		x = r->xMax;
}

Event::Event(const Event& e) : r(e.r), left(e.left), x(e.x)
{
}

Event::~Event()
{
}

void swap(Event& first, Event& second)
{
	using std::swap;
	swap(first.r,second.r);
	swap(first.left,second.left);
	swap(first.x,second.x);
}

//Event::Event(Event&& e) : r(NULL), left(false), x(0)
//{
//	swap(*this,e);
//}

Event& Event::operator=(Event e)
{
	swap(*this,e);
	return *this;
}

void Event::print()
{
	r->print();
	if (left)
		printf("xMin: %6.2f\n",x);
	else
		printf("xMax: %6.2f\n",x);
}
