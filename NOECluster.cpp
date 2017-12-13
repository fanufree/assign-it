/*
 * NOECluster.cpp
 *
 *  Created on: 2012-10-08
 *      Author: e4k2
 */

#include <algorithm>
#include <cstdio>
#include <cmath>
#include "NOECluster.h"

NOECluster::NOECluster() : noes()
{
}

NOECluster::NOECluster(NOE* n) : noes()
{
	noes.push_back(n);
}

NOECluster::NOECluster(const NOECluster& n) : noes()
{
	noes.insert(noes.begin(), n.noes.begin(),n.noes.end());
}

//NOECluster::NOECluster(NOECluster&& n) : noes()
//{
//	swap(*this,n);
//}

NOECluster::~NOECluster()
{
}

void swap(NOECluster& c1, NOECluster& c2)
{
	using std::swap;
	swap(c1.noes,c2.noes);
}

NOECluster& NOECluster::operator=(NOECluster n)
{
	swap(*this,n);
	return *this;
}

void NOECluster::print()
{
	printf("Begin NOECluster\n");
	for (list<NOE*>::iterator it = noes.begin(); it != noes.end(); ++it)
	{
		NOE* np = *it;
		np->print();
	}
	printf("End NOECluster\n");
}

bool NOECluster::testAndAdd(NOE* n, double tolx, double tolhx)
{
	bool matchAll = true;
	for (list<NOE*>::iterator it = noes.begin(); it != noes.end(); ++it)
	{
		NOE* np = *it;
		if (fabs(np->x-n->x) > tolx || fabs(np->hx-n->hx) > tolhx)
		{
			matchAll = false;
			break;
		}
	}
	if (matchAll)
		noes.push_back(n);
	return matchAll;
}

bool NOECluster::test(NOECluster* c, double dist) const
{
	bool matchAll = true;
	for (list<NOE*>::const_iterator it1 = noes.begin(); it1 != noes.end(); ++it1)
	{
		NOE* np1 = *it1;
		for (list<NOE*>::const_iterator it2 = c->noes.begin(); it2 != c->noes.end(); ++it2)
		{
			NOE* np2 = *it2;
			double d = fabs(np1->x-np2->x)+10.0*fabs(np1->hx-np2->hx);
			if (d > dist)
			{
				return false;
			}
		}
	}
	if (matchAll)
		return true;
	return false;
}

bool NOECluster::test(NOE* n, double dist) const
{
	bool matchAll = true;
	for (list<NOE*>::const_iterator it1 = noes.begin(); it1 != noes.end(); ++it1)
	{
		NOE* np1 = *it1;
		double d = fabs(np1->x-n->x)+10.0*fabs(np1->hx-n->hx);
		if (d > dist)
		{
			return false;
		}
	}
	if (matchAll)
		return true;
	return false;
}

list<NOECluster*>::iterator NOECluster::split(list<NOECluster*>::iterator itAll, list<NOECluster*> allClusters, double tolx, double tolhx)
{
	double distCut = 3.0*min(tolx,10.0*tolhx)/4.0; // for testing noe similarity; peaks considered similar if < distCut
	NOECluster* c = *itAll;
	NOE* startNOE = c->front();
	// at most 4 possible ways to split cluster c
	NOE* minX = startNOE; // noe with min x chemical shift in cluster c
	NOE* maxX = startNOE;
	NOE* minHX = startNOE;
	NOE* maxHX = startNOE;

	// find the above noes
	list<NOE*>::iterator startIt = c->noes.begin();
	if (startIt != c->noes.end())
	{
		++startIt;
		for (list<NOE*>::iterator it = startIt; it != c->noes.end(); ++it)
		{
			NOE* np = *it;
			if (np->x < minX->x)
				minX = np;
			if (np->x > maxX->x)
				maxX = np;
			if (np->hx < minHX->hx)
				minHX = np;
			if (np->hx > maxHX->hx)
				maxHX = np;
		}
	}

	// test if can split cluster c
	list<NOECluster*> clusters; // putative new cluster centers from noes in cluster c
	clusters.push_back(new NOECluster(minX));
	// ignore duplicates; and make sure centers are different
	if (minX != maxX) // test if can add maxX
	{
		bool allDiff = true;
		for (list<NOECluster*>::iterator it = clusters.begin(); it != clusters.end(); ++it)
		{
			NOECluster* nc = *it;
			if (nc->test(maxX,distCut))
			{
				allDiff = false;
				break;
			}
		}
		if (allDiff)
			clusters.push_back(new NOECluster(maxX));
	}
	if (minHX != minX && minHX != maxX) // test if can add minHX
	{
		bool allDiff = true;
		for (list<NOECluster*>::iterator it = clusters.begin(); it != clusters.end(); ++it)
		{
			NOECluster* nc = *it;
			if (nc->test(minHX,distCut))
			{
				allDiff = false;
				break;
			}
		}
		if (allDiff)
			clusters.push_back(new NOECluster(minHX));
	}
	if (maxHX != minX && maxHX != maxX && maxHX != minHX) // test if can add maxHX
	{
		bool allDiff = true;
		for (list<NOECluster*>::iterator it = clusters.begin(); it != clusters.end(); ++it)
		{
			NOECluster* nc = *it;
			if (nc->test(maxHX,distCut))
			{
				allDiff = false;
				break;
			}
		}
		if (allDiff)
			clusters.push_back(new NOECluster(maxHX));
	}

	if (clusters.size() == 1)
	{   // stayed as one cluster; no split
		// clean up
		for (list<NOECluster*>::iterator it = clusters.begin(); it != clusters.end(); ++it)
		{
			NOECluster* nc = *it;
			delete nc;
			it = clusters.erase(it);
		}
		return ++itAll; // no split; update the iterator
	}

	// assign the other noes in cluster c to closest cluster; ignore the new cluster centers
	for (list<NOE*>::iterator it = c->noes.begin(); it != c->noes.end(); ++it)
	{
		NOE* n = *it;
		// check if cluster center
		bool match = false; // true if n is one of the cluster centers
		// add to cluster whose center is closest to n
		double minD = 10000000;
		int minClusterIndex = -1; // the index of the closest cluster center
		int clusterIndex = 0;
		for (list<NOECluster*>::iterator itC = clusters.begin(); itC != clusters.end(); ++itC)
		{
			NOECluster* nc = *itC;
			NOE* cen = nc->front();
			if (n == cen)
			{
				match = true;
				break;
			}
			else
			{
				double d = fabs(n->x-cen->x)+10.0*fabs(n->hx-cen->hx);
				if (d < minD)
				{
					minD = d;
					minClusterIndex = clusterIndex;
				}
			}
			clusterIndex++;
		}
		if (match)
		{
			continue;
		}
		else if (minClusterIndex != -1) // add n to it's new cluster
		{
			list<NOECluster*>::iterator itC = clusters.begin(); // jump to cluster at minClusterIndex
			for (int i = 0; i < minClusterIndex; i++)
				++itC;
			NOECluster* nc = *itC;
			nc->add(n);
		}
	}

	// decide if can merge the new clusters
	for (list<NOECluster*>::iterator itC1 = clusters.begin(); itC1 != clusters.end();)
	{
		NOECluster* nc1 = *itC1; // nc1 gets added to nc2 and then nc1 deleted if can merge
		double x1 = 0;
		double hx1 = 0;
		nc1->getAvg(x1,hx1);
		bool merged = false;
		for (list<NOECluster*>::iterator itC2 = clusters.begin(); itC2 != clusters.end(); ++itC2)
		{
			NOECluster* nc2 = *itC2;
			if (nc1 < nc2)
			{
				double x2 = 0;
				double hx2 = 0;
				nc2->getAvg(x2,hx2);
				double d = fabs(x1-x2)+10.0*fabs(hx1-hx2);
				if (d < distCut)
				{
					merged = true;
					nc2->merge(nc1);
					break;
				}
			}
		}
		if (merged)
		{
			delete nc1;
			itC1 = clusters.erase(itC1);
		}
		else
			++itC1;
	}

	// filter the clusters by size; if >= 3, keep as separate cluster
	// if < 3, delete
	for (list<NOECluster*>::iterator itC = clusters.begin(); itC != clusters.end();)
	{
		NOECluster* nc = *itC;
		if (nc->size() < 3)
		{
			delete nc;
			itC = clusters.erase(itC);
		}
		else
			++itC;
	}

	// replace cluster c with the above clusters if number of clusters > 1
	// delete old cluster
	if (clusters.size() > 1)
	{
		delete c;
		itAll = allClusters.erase(itAll); // itAll points to item X, which is after the deleted item
		allClusters.insert(itAll, clusters.begin(),clusters.end()); // inserts just before item X; itAll still points to X
		// do not delete elements (NOECluster*) in clusters since they are now in allClusters
		return itAll;
	}
	else
	{
		// no split
		// but size might have shrank, so replace c if shrank in size as long as size > 3; else delete c
		if (clusters.size() == 1)
		{
			NOECluster* nc = clusters.front();
			if (nc->size() < c->size())
			{
				if (nc->size() >= 3)
				{
					c->replace(nc);
				}
				else
				{
					// delete c
					delete c;
					itAll = allClusters.erase(itAll);
					delete nc;
					return itAll;
				}
			}
			// clean up
			delete nc;
			return ++itAll;
		}
		else
		{
			// the cluster filter step above deleted all clusters because each are too small
			delete c;
			itAll = allClusters.erase(itAll);
			return itAll;
		}
	}
}

int NOECluster::size() const
{
	return noes.size();
}

bool NOECluster::contains(NOE* n) const
{
	if (std::find(noes.begin(),noes.end(),n) != noes.end())
		return true;
	else
		return false;
}

void NOECluster::merge(NOECluster* c)
{
	for (list<NOE*>::iterator it = c->noes.begin(); it != c->noes.end(); ++it)
	{
		NOE* np = *it;
		if (!contains(np))
			noes.push_back(np);
	}
}

void NOECluster::add(NOE* noe)
{
	noes.push_back(noe);
}

NOE* NOECluster::front()
{
	return noes.front();
}

void NOECluster::replace(NOECluster* c)
{
	noes.clear();
	noes.insert(noes.begin(),c->noes.begin(),c->noes.end());
}

bool NOECluster::remOverlap(NOECluster* c)
{
	bool remSomething = false;
	for (list<NOE*>::iterator it1 = noes.begin(); it1 != noes.end();)
	{
		NOE* np1 = *it1;
		bool rem = false;
		for (list<NOE*>::iterator it2 = c->noes.begin(); it2 != c->noes.end();)
		{
			NOE* np2 = *it2;
			if (np1 == np2)
			{
				it2 = c->noes.erase(it2);
				rem = true;
				break;
			}
			else
				++it2;
		}
		if (rem)
		{
			it1 = noes.erase(it1);
			remSomething = true;
		}
		else
			++it1;
	}
	return remSomething;
}

void NOECluster::getAvg(double& x, double& hx) const
{
	x = 0;
	hx = 0;
	for (list<NOE*>::const_iterator it = noes.begin(); it != noes.end(); ++it)
	{
		NOE* np = *it;
		x += np->x;
		hx += np->hx;
	}
	int size = noes.size();
	x = x/double(size);
	hx = hx/double(size);
}
