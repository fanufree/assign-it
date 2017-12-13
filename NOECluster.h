/*
 * NOECluster.h
 *
 *  Created on: 2012-10-08
 *      Author: rjang
 */

#ifndef NOECLUSTER_H_
#define NOECLUSTER_H_

#include <list>
#include "NOE.h"

using namespace std;

/*
 * Used for calibrating the x, hx chemical shifts of NOE peak lists (for 3D or 4D)
 * Used for grouping NOEs with similar x,hx chemical shifts (represents one NOECluster)
 * Clusters can get split
 */
class NOECluster
{
public:
	NOECluster();
	NOECluster(NOE* n);
	NOECluster(const NOECluster& n);
	virtual ~NOECluster(); // does not delete NOE* in list noes because they are not created by this class
	friend void swap(NOECluster& c1, NOECluster& c2);
	// NOECluster(NOECluster&& n);
	NOECluster& operator=(NOECluster n);
	void print();
	bool testAndAdd(NOE* n, double tolx, double tolhx); // returns true if n was added; to be added, n must match all noes within the given tolerances
	bool test(NOECluster* c, double dist) const; // true if noes in this and c all overlap by dist = dx+10*dhx
	bool test(NOE* n, double dist) const; // test for dist = dx+10*dhx overlap

	/**
	 * Test if cluster *itC can be split
	 *  Deletes *itC from allClusters if it can be split, and replaces with the new groups which get added to allClusters
	 *  The correct iterator location is returned
	 *  If a cluster has < 3 elements but there exists at least another cluster >= 3 elements, the < 3 groups will get deleted
	 *  If splitting results in each cluster < 3, then they are all deleted
	 *  Size of input cluster *itC can shrink. If *itC shrinks to < 3, then it is deleted
	 *  tolx and tolhx are used to determine if the split clusters should stay split; that is, the peaks must all be tolx+10*tolhx apart
	 */
	static list<NOECluster*>::iterator split(list<NOECluster*>::iterator itC, list<NOECluster*> allClusters, double tolx, double tolhx);
	void merge(NOECluster* c); // used by split; merge c with this cluster; ignores noes already in this cluster; does not delete c
	bool contains(NOE* noe) const; // used by merge; true if this contains noe
	void add(NOE* noe); // used by split; adds without checking if noe already in noes
	NOE* front(); // used by split; returns first peak in cluster
	void replace(NOECluster* c); // used by split; replaces the contents of noes with c.noes
	void getAvg(double& x, double& hx) const; // used by split; values are undefined if size() == 0

	int size() const;
	bool remOverlap(NOECluster* c); // used in Main.cpp; removes peaks common to both clusters from both clusters; returns true if something was removed
	list<NOE*> noes;
};

#endif /* NOECLUSTER_H_ */
