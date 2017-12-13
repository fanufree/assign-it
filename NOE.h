/*
 * NOE.h
 *
 *  Created on: 2012-05-03
 *      Author: e4k2
 */

#ifndef NOE_H_
#define NOE_H_

#include <list>
#include <string>

using namespace std;

/**
 * Can support 2D, 3D, and 4D NOEs.
 * Need to write readNOE2D and readNOE4D functions
 */
class NOE
{
public:
	static const int NUMTYPES=7; // num peak list types; excludes UNKNOWN
	enum PEAKLISTTYPE { N15_3D, CALI_3D, CARO_3D, CHCH_4D, NHNH_4D, CHNH_4D, H2D, UNKNOWN }; // UNKNOWN only used for default value
	static const string PEAKTYPENAMES[]; // = {"N15_3D", "CALI_3D", "CARO_3D", "CHCH_4D", "NHNH_4D", "CHNH_4D", "H2D"}; // in same order as PEAKLISTTYPE
	static const char HEAVYATOM1TYPE[]; // = {'N','C','C','C','N','C','X'}; // for iterating over the chemical shifts; indirect dimension
	                                                    // (see computeScore); X = C or N; in same order as PEAKLISTTYPE
	static const char HEAVYATOM2TYPE[]; // = {'X','X','X','C','N','N','X'}; // for iterating over the chemical shifts; direct dimension

	double x;  // C or N; (set to INVALIDSHIFTVAL if 2D NOE)
	double hx;  // HC or HN attached to x
	double x2; // C or N attached to h (for 4D NOEs, otherwise set to INVALIDSHIFTVAL)
	double h; // indirect H
	double volume; // < 0 if not available
	PEAKLISTTYPE type;

	NOE();
	NOE(double xx, double hxx, double hh, double vol, PEAKLISTTYPE t);
	NOE(double hxx, double hh, double vol, PEAKLISTTYPE t);
	NOE(double xx, double hxx, double xx2, double hh, double vol, PEAKLISTTYPE t);

	~NOE();
	NOE(const NOE& n);
	// NOE(NOE&& n);
	friend void swap(NOE& first, NOE& second);
	NOE& operator=(NOE n);
	int getNumDimensions() const;
	void print();
	static void print(list<NOE*>& noes);
	void reference(double offset_x, double offset_hx); // adds the given values to x and hx respectively
	void reference2(double offset_x2, double offset_h); // adds the given values to x2 and h respectively
	void reference(double offset_h); // adds offset to h
	/**
	 * Remove and deletes duplicate noes using the given chemical shift tolerances
	 * Keeps the NOE with the larger intensity
	 */
	static void removeDuplicates3D(list<NOE*>& noes, double cTol, double hcTol, double hTol);
	/**
	 * NOEs get read into noes
	 * If noes not empty, new NOEs get appended
	 * Caller is responsible for deleting the pointers in list noes
	 * If volumneIndex is < 0, then volume set to -1
	 */
	static void readNOE3D(const string& filename, int ncIndex, int hncIndex,
			int indirectIndex, int volumeIndex,
			list<NOE*>& noes, PEAKLISTTYPE t); // Indices must be 1-based

	// TODO: readNOE2D, readNOE4D

	/*
	 * proton 1 is attached to heavy 1, with chemical shifts proton1 and heavy1 respectively
	 * Ignores x, hx, h, x2 with INVALIDSHIFTVAL chemical shift value;
	 * To ignore specific atoms, use INVALIDSHIFTVAL e.g. to ignore heavy1 column, pass in heavy1=INVALIDSHIFTVAL
	 * A match is true, if the given chemical shifts are <= the given tolerances for the values that aren't ignored
	 */
	bool match(double heavy1, double proton1, double heavy2, double proton2,
			double tolHeavy1, double tolProton1, double tolHeavy2, double tolProton2);

	// returns abs(CS-this.cs) or INVALIDSHIFTVAL if cs does not exist or the argument passed in is INVALIDSHIFTVAL
	void getShiftDiff(double heavy1CS, double proton1CS, double heavy2CS, double proton2CS,
			          double& diffx1, double& diffh1, double& diffx2, double& diffh2);
};

struct NOESortByVolume
{
	bool operator()(NOE* a, NOE* b)
	{
		return a->volume < b->volume;
	}
};

#endif /* NOE_H_ */
