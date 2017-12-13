/*
 * Residue.h
 *
 *  Created on: 2012-04-27
 *      Author: e4k2
 */

#ifndef RESIDUE_H_
#define RESIDUE_H_

#include <string>
#include <map>
#include <vector>
#include "Atom.h"

class Residue;

/**
 * Used by Residue.
 * Iterates over the heavy atoms starting with the carbons then the nitrogens
 *  Or iterate over only the carbons or only the nitrogens
 *  See Residue::begin(char type)
 *  See Atom::begin() to iterate over the protons of the heavy atom
 */
class AtomIterator
{
public:
	AtomIterator(map<string,Atom*>::const_iterator iterC, map<string,Atom*>::const_iterator iterN,
			const Residue* res);
	~AtomIterator();
	AtomIterator(const AtomIterator& iter);
	friend void swap(AtomIterator& first, AtomIterator& second);
	// AtomIterator(AtomIterator&& iter);
	AtomIterator& operator=(AtomIterator iter);
	AtomIterator& operator++();
	AtomIterator operator++(int); // postfix operator
	bool operator==(const AtomIterator& iter) const;
	bool operator!=(const AtomIterator& iter) const;
	Atom* operator*(); // returns NULL if no more atoms
private:
	map<string,Atom*>::const_iterator itC;
	map<string,Atom*>::const_iterator itN;
	const Residue* r;
};

class Residue
{
public:
	int num;  // residue number
	int type; // one of AA3
	map<string,Atom*> cMap; // carbon name to it's atom object
	map<string,Atom*> nMap;  // nitrogen name to it's atom object
	Atom** hList; // Array of Atom* for the protons; NULL if numProtons 0
	int numProtons; // size of hList
	Residue(int resnum, int aaType);
	Residue(const Residue& res);
	~Residue();
	friend void swap(Residue& first, Residue& second);
	// Residue(Residue&& res);
	Residue& operator=(Residue res);
	void print();
	const char* getAAType1() const; // 1 letter name
	const char* getAAType3() const; // 3 letter name

	void addShift(const string& atomname, double shift);
	void addCoords(const string& atomname, double x, double y, double z);
	Atom* getC(const string& name) const; // null if not present
	Atom* getN(const string& name) const;
	Atom* getX(const string& name) const; // returns the given heavy atom; name must be of type C or N; null returned if not found
	int getNumProtons() const;
	int getNumProtonsMethylOnce() const; // returns num protons while treating each methyl group as only 1 proton instead of 3
	AtomIterator begin(char type) const; // iterate over  type='C', or type='N' or type'X'='C' and 'N'
	AtomIterator end(char type) const;
	void getSideChainCentroid(double& x, double& y, double& z) const; // heavy atoms only; no protons
	                                                  // excludes CA, except for GLY where it is used

	bool isCaro(const string& catomname) const; // true if this residue is an aromatic residue and the given
	                                           // carbon atom name is in the aromatic ring
	bool isCali(const string& catomname) const; // true if the given carbon atom is non-aromatic
	bool is15N(const string& natomname) const; // true if the given atom is a nitrogen
	bool isHydrophobic() const;
	bool isPositive() const; // is basic
	bool isNegative() const; // is acid
	bool isCYS() const;

	friend class AtomIterator;
	enum AA3 {ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL};
	enum AA1 {A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,  F,  P,  S,  T,  W,  Y,  V};
	static const char* AA3S[];
	static const char* AA1S[];
	static int index1(const char* aaName1); // returns AA1 with same name
	static int index1(char aaName1);
	static int index3(const char* aaName3); // returns AA3 with same name
	static bool isStandard3(const char* aaName3); // true if aaName3 is one of the 20 standard amino acids
	static bool isStandard1(const char* aaName1);
	static bool isStandard1(char aaName1);
private:
	enum ADDTYPEH {H0, H1, H2, H3}; // number of H's to add, used by addCN
	void addCN(bool cFlag, const string& name, ADDTYPEH addH, int& hIndex); // for creating and adding C or N atom
	   // to this residue and their corresponding h's, hIndex keeps track of next free position in hList for adding hs
	void makeBackbone(int& hIndex); // makes backbone atoms
	void makeBackboneHB2(int& hIndex); // makeBackbone plus CB,HB2,HB3
	/*
	 * Num protons considered
	ALA 5
	ARG 9
	ASN 6
	ASP 4
	CYS 4
	GLN 8
	GLU 6
	GLY 3
	HIS 8
	ILE 11
	LEU 11
	LYS 13
	MET 9
	PHE 9
	PRO 7
	SER 4
	THR 7
	TRP 10
	TYR 8
	VAL 9
	 */
	void makeALA();
	void makeARG();
	void makeASN();
	void makeASP();
	void makeCYS();
	void makeGLN();
	void makeGLU();
	void makeGLY();
	void makeHIS();
	void makeILE();
	void makeLEU();
	void makeLYS();
	void makeMET();
	void makePHE();
	void makePRO();
	void makeSER();
	void makeTHR();
	void makeTRP();
	void makeTYR();
	void makeVAL();
};

#endif /* RESIDUE_H_ */
