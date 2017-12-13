/*
 * Rectangle.h
 *
 *  Created on: 2012-07-16
 *      Author: e4k2
 */

#ifndef RECTANGLE_H_
#define RECTANGLE_H_

#include <vector>
#include <tr1/array>

using namespace std;

/**
 * Axis Aligned 2D rectangle
 */
class Rectangle
{
public:
	double xMin;
	double xMax;
	double yMin;
	double yMax;
	double width; // xMax-xMin
	double height; // yMax-yMin
	double dx; // annotated by getMaximumClique
	double dy;

	Rectangle();
	Rectangle(double xMin, double xMax, double yMin, double yMax);
	Rectangle(const Rectangle& r);
	virtual ~Rectangle();
	friend void swap(Rectangle& first, Rectangle& second);
	// Rectangle(Rectangle&& r);
	Rectangle& operator=(Rectangle r);

	bool contains(double x, double y);

	/**
     * Moves this rectangle to the origin based on its centroid
	 */
    void center();

    void translate(double x, double y);

    /**
     * True if intersect
     */
   bool intersectBool(const Rectangle& r);

   /**
	* Returns the intersection of this and r as the rectangle ret.
	* Returns false if there's no intersection (in which case ret is not modified)
	* this and ret can be the same rectangle object
	*/
   bool intersect(const Rectangle& r, Rectangle& ret);

   void print();

   void getMidPoint(double& x, double& y);

   /**
    * Returns closest point on perimeter (does not include interior) of rectangle to origin
    */
   void getClosestPointToOrigin(double& x, double &y);

	/**
	 * Not the most efficient implementation since we don't use any spatial index (e.g. augmented tree)
	 * O(n^2). Possible to improve to O(nlogn) if use spatial index.
     * Returns the intersection in output and the clique size in cliqueSize
     * @param output
     * A rectangle representing the intersection region
     * @param recsMaxClique
     * The rectangles in input in the maximum clique
	 */
   static void getMaximumClique(const vector<Rectangle*>& input, Rectangle& output, int& maxCliqueSize,
		   vector<Rectangle*>& recsMaxClique);

   /**
    * Returns the calibration offsets offsetN, offsetHN
    * Function can also be used to calibration C, HC NOEs
    * @param ref
    * Reference peaks. Usually from BMRB
    * @param input
    * The peaks whose chemical shifts will have offsetN, offsetHN added to them (to be added by caller).
    * Usually the NOEs
    * @param matchN
    * @param matchHN
    * Chemical shifts a, b are said to match if their differences are within these tolerances
    * @param searcN
    * @param searchHN
    * maximum values for offsetN, offsetH
    * @param debug
    */
   static void calibrate(const vector< tr1::array<double,2> >& ref, const vector< tr1::array<double,2> >& input,
		   double matchTolN, double matchTolHN, double searchN, double searchHN,
		   double& offsetN, double& offsetHN, bool debug);

   /**
    * For testing purposes. Unique peak match calibration. True if exist unique peak match.
    */
   static bool calibrateUnique(const vector< tr1::array<double,2> >& ref, const vector< tr1::array<double,2> >& input,
		   double searchN, double searchHN,
		   double& offsetN, double& offsetHN, bool debug);
};

struct RectangleSortByDx
{
	bool operator()(Rectangle* a, Rectangle* b)
	{
		return a->dx < b->dx;
	}
};

struct RectangleSortByDy
{
	bool operator()(Rectangle* a, Rectangle* b)
	{
		return a->dy < b->dy;
	}
};

#endif /* RECTANGLE_H_ */
