/*
 * Event.h
 *
 *  Created on: 2012-07-16
 *      Author: e4k2
 */

#ifndef EVENT_H_
#define EVENT_H_

#include "Rectangle.h"

/**
 * Represents a vertical edge of some rectangle
 * Reference:
 * Hiroshi Imai and Takao Asano. Finding the connected components and a maximum clique of an
 * intersection graph of rectangles in the plane. Journal of Algorithms, 4(4):310 – 323, 1983.
 */
class Event
{
public:
	Event();
	Event(Rectangle* r, bool left);
	Event(const Event& e);
	virtual ~Event();
	friend void swap(Event& first, Event& second);
	// Event(Event&& e);
	Event& operator=(Event e);
	void print();

	Rectangle* r;
	bool left;
	double x;
};

struct EventSortByX
{
	/**
	 * sort by x. If same x, then the one with left=true gets priority
	 * (so intersection can be checked before it is removed)
	 */
	bool operator()(const Event& a, const Event& b)
	{
		if (a.x != b.x)
		{
			return a.x < b.x;
		}
		else
		{
			if (a.left && !b.left)
				return true; // a < b
			else // if (!a.left && b.left)
				// return false; // a > b
			// else
				return false; // a=b
		}
	}
};

#endif /* EVENT_H_ */
