/*
 * Atom.h
 *
 *  Created on: 2012-04-29
 *      Author: e4k2
 */

#ifndef ATOM_H_
#define ATOM_H_

#define INVALIDCOORD -999999
#define INVALIDSHIFTVAL -300.0

#include <list>
#include <string>

using namespace::std;

class Atom;

typedef Atom** HIterator;

/**
 * TODO: This should really be 3 classes - one for heavy atoms, one for protons, one base atom class
 *
 * Heavy atoms may have attached protons
 * Proton atoms do not have attached protons (NULL)
 */
class Atom
{
public:
	string name;
    double x;
    double y;
    double z;
    double cs; // chemical shift
	Atom** protons; // NULL if numProtons is 0
	int numProtons;
    Atom();
	Atom(const string& atomName, int numP); // numP = # of protons if this is a heavy atom
	Atom(const Atom& atom);
	virtual ~Atom();
	friend void swap(Atom& first, Atom& second);
	// Atom(Atom&& atom);
	Atom& operator=(Atom atom);
	Atom* operator[](const int index); // return proton at 0-based index or NULL
	void print();
	void addH3(); // if c has name CX or CXY, makes HX1-3 or HXY1-3
	void addH2(); // if CX or CXY, makes HX2,HX3 or HXY2,HXY3
	void addH1(); // if CX or NX, adds HX
	void addH2_N(); // for NX or NXY, makes HX1-2 or HXY1-2
	Atom* getH(const string& name) const; // NULL if not found
	int getNumProtons() const;
	void setShift(double shift);
	void setCoords(double cx, double cy, double cz);
	bool isHeavy() const; // true if this is a heavy atom (a non-proton)
	bool isNMR_SideChainHeavy() const; // true if this a carbon or nitrogen, and non-backbone, and not CA and CB
	                        // ie. a sidechain heavy atom that is not CA or CB or S or O
	bool isBackbone() const; // true if this atom is C,N,H,HA,CA,terminal H1 H2 H3,O,OXT
	bool isProton() const;
	bool isMethyl() const; // true if this is a heavy atom containing 3 protons; assumes this is the heavy atom; do not call this if this Atom is a proton
	bool isCB() const;   // true if this is CB
	bool isGeminal() const; // assumes this is the heavy atom containing protons
	bool isProtonGroupChemShift(double tol) const; // true if this is a heavy atom and the pairwise chemical shifts of
	                                       // the protons are all within the given tolerance
	                                     // e.g. heavy atoms with either methyl and geminal protons tend to return true here
	char getType() const; // returns atom type: 'H', 'C' or 'N'
	double getDistance(const Atom* a) const; // assumes the coodinates exists
	bool hasCoordinates() const; // true if x,y,z have been initialized
	bool hasChemShift() const; // true if cs has been initialized
	/*
	 * Iterator over the protons (if this is a heavy atom)
	 * e.g.
	 * 	for (HIterator it = cb->beginH(); it != cb->endH(); it++)
		{
			Atom* h = *it;
		}
	 */
	HIterator beginH() const;
	HIterator endH() const;
	friend class Residue;
};




#endif /* ATOM_H_ */
