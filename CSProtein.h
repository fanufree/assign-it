/*
 * CSProtein.h
 *
 *  Created on: 2012-04-29
 *      Author: e4k2
 */

#ifndef CSPROTEIN_H_
#define CSPROTEIN_H_

#define MAXPROSIZE 250   // Used in Main.cpp for max size of isDynamic, ssTypes

#include <string>
#include <list>
#include "Contact.h"

class Residue;
class PseudoContact;

using namespace std;

class CSProtein
{
public:
	int startResNum; // the first residue number, usually 1 unless specified otherwise in constructor
	        // methods in this class access residues by their 1-based residue number, so access to the
	        // residue array is achieved by residues[resNum-startResNum]
	int size; //  num residues

	CSProtein(const string& seq, int startingResNum=1); // residues labeled X are ignored;
	                                        // startingResNum: first residue number; normally 1
	CSProtein(const CSProtein& other);
	~CSProtein();
	friend void swap(CSProtein& first, CSProtein& second);
	// CSProtein(CSProtein&& other);
	CSProtein& operator=(CSProtein other);
	void getShifts(const string& bmrbFile, int offset=0); // reads in chemical shifts from BMRB file; offset is added to
	                                            // the residue numbers in the BMRB file to give the actual (CSProtein) residue numbers
	                                           // (used in case the BMRB numbers do not match CSProtein)
	void getCoords(const string& pdbFile, int offset=0); // offset is added to the residue number in pdbFile to give
	                                                  // the residue numbers used by this object
	                                                 // only the first model is read if there are multiple models

	int getNumProtons() const; // total number of protons. Note methyls count as 3 protons

	/**
	 * Only considers chain A or nameless chain among all the models
	 * Similar to getCoords except all models are read. The coordinates are stored in models
	 */
	static void getModels(const string& pdbFile, list<CSProtein*>& models, const string& seq, int startingResNum, int offset=0);

	Residue* operator[](int resnum); // resnum (1-based) must be >= startResNum, returns the residue at index = resnum-startResNum
	                                // To get residue with resnum i, pass in i, NOT i-1
	void print();

	/**
	 * Returns true if contact within the distCutoff
	 * All protons in methyl groups are used to determine if there is a contact. The minimum distance is returned in retDist.
	 * retDist is -1 if coordinates do not exist (and false is returned by inContact).
	 * retDist can be bigger than distCutoff in which case false is returned
	 * This method is used by all inContact functions below (either directly or indirectly)
	 */
	bool inContact(int resNum1, const string& parent1, const string& proton1,
			int resNum2, const string& parent2, const string& proton2,
			double distCutoff, double& retDist) const;

	bool inContact(int resNum1, const string& parent1, const string& proton1,
			int resNum2, const string& parent2, const string& proton2,
			double distCutoff) const;

	bool inContact(PseudoContact& pc, double distCutoff) const; // true if at least one contact is within the cutoff
	bool inContact(Contact& c, double distCutoff) const;

	/**
	 * retDist stores the distance or min distance for PseudoContact.
	 * It is set to -1 if coordinates do not exist
	 */
	bool inContact(Contact& c, double distCutoff, double& retDist) const;
	bool inContact(const Contact& c, double distCutoff, double& retDist) const;
	bool inContact(PseudoContact& pc, double distCutoff, double& retDist) const;

	/**
	 * True if at least one model has the contact
	 * retDist is negative if coordinates do not exist (and false returned); otherwise it is the min dist among the models
	 */
	static bool inContact(list<CSProtein*>& models, PseudoContact& pc, double distCutoff);
	static bool inContact(list<CSProtein*>& models, PseudoContact& pc, double distCutoff, double& retDist); // min dist among models returned in retDist
	static bool inContact(list<CSProtein*>& models, Contact& c, double distCutoff);
	static bool inContact(list<CSProtein*>& models, PseudoContact& pc, double distCutoff, double& retDist, int& numModels); // min dist among models returned in retDist
                                                        // returns num models that have this contact within distCutoff
	static bool inContact(list<CSProtein*>& models, Contact& c, double distCutoff, double& retDist);
	static bool inContact(list<CSProtein*>& models, const Contact& c, double distCutoff, double& retDist);
	static bool inContact(list<CSProtein*>& models, Contact& c, double distCutoff, double& retDist, int& numModels);

	/**
	 * The min proton-proton dist among all the contacts in pc. All protons in methyl groups are tested.
	 * -1 is returned if coordinates do not exist
	 */
	double getMinDistance(PseudoContact& pc);

	/**
	 * gets the min of the side chain centroid distances between the given residues in the given models
	 * returns -1 if residue or coordinates do not exist
	 */
	static double getSideChainCentroidDistance(list<CSProtein*>& models, int res1, int res2);
	double getSideChainCentroidDistance(int res1, int res2) const; // returns -1 if residue or coodinates do not exist
	double getSideChainCentroidCADistance(int res1, int res2) const; // returns -1 if residue or coodinates do not exist; CA atom for res2
	double getSideChainCentroidCBDistance(int res1, int res2) const; // returns -1 if residue or coodinates do not exist; CB atom for res2
	double getSideChainCentroidNDistance(int res1, int res2) const; // returns -1 if residue or coodinates do not exist; N atom for res2
	double getDistance(int res1, int res2, Contact::ResidueContactType type) const; // return -1 if residue or coodinates do not exist;


	/**
	 * returns -1 if coordinates do not exist for one of the protons
	 * Does NOT take average distance between methyl groups; each methyl proton is treated as unique
	 */
	double getDistance(int resNum1, const string& parent1, const string& proton1,
	                   int resNum2, const string& parent2, const string& proton2) const;

	/**
	 * Writes N, HN (including side chain ASN, GLN, TRP N,HN's) chemical shifts in BMRB 2.0 format
	 * If ha=true, then ha and ca chem shifts are also written. Similarly for hb=true
	 * offset is added to res num
	 */
	void write15NChemShift(const string& filename, int offset, bool ha, bool hb) const;
private:
	Residue** residues;
};

#endif /* CSPROTEIN_H_ */
