/*
 * SSType.h
 *
 *  Created on: 2013-08-01
 *      Author: e4k2
 */

#ifndef SSTYPE_H_
#define SSTYPE_H_

#include <vector>
#include <string>
#include <list>
#include <tr1/array>
#include "CSProtein.h"

using namespace std;

/*
 * Functions for parsing stride file and accessing the secondary structure information
 */
enum SSTYPE {HELIX, STRAND, LOOP, SHEET, UNKNOWNSS}; // SHEET is reserved for a group of strands; UNKNOWN for unspecified or nonconfident predictions

const int INVALIDSS = 16; // for marking for deletion
const int TYPEINDEX = 0;
const int STARTRESINDEX = 1;
const int ENDRESINDEX = 2;
const int SSIDINDEX = 3;
const int SHEETIDINDEX = 4;     // for strands only
const int NUMSTRANDSINDEX = 1;  // for sheets only
const int SSIDINDEX_SHEET = 2; // for sheets only
const int STRANDSTARTINDEX = 3; // for sheets only
const int DEFAULTID = -1; // indicates no id available

const int SHEETLOWERBOUNDLEN = 2; // sheet len must be greater than this to be considered as interacting
const int HELIXLOWERBOUNDLEN = 3; // helices must be greater than this
const int MINHBONDS = 3; // min num hbonds between a pair of ss elements before can consider as interacting
/*
 * @param sse
 * list of secondary structure elements
 * Each ss element/row has the following entries
 * [0]= ss type HELIX, STRAND, LOOP, SHEET
 * (strands will be contained inside some sheet)
 * [1]=start res 1, [2] = end res 2 (inclusive)
 * [3] = a unique id for the ss element. Equal to the index of this element in the vector (0-based)
 * [4] = for strands only, the sheet index into sse
 * The entry for sheet is slightly different:
 * If [0] == SHEET, [1] = number of strands K, [2] sheet index into sse,
 * [3]...[K+2] give the strand indices of all K strands in this sheet
 */
void parseStride(string stride, vector< vector<int> >& sse, int pdbOffset);

/**
 * Similar to parseStride except for ssTypes[resnum-1]
 * Assumes ssTypes[i] = UNKNOWNSS if ss type unknown or residue does not exist
 */
void parseSSTypes(tr1::array<SSTYPE,MAXPROSIZE>& ssTypes, vector< vector<int> >& sse);

/**
 * Returns the secondary structure element of res1 and res2 sse[ssid1][0], sse[ssid2][0]
 * If both res1 and res2 are in the same sheet, then the returned ids (index into sse) of their strands
 * are returned; else the id is set to their sheets
 * Returns -1 for the id if res is in a loop or short segment
 */
void getContactSSType(vector< vector<int> >& sse, int res1, int res2, int& ssid1, int& ssid2);

/**
 * Returns the residues in the given secondary structure element with id ssid
 * If ssid is a sheet, then the residues in all its strands are turned
 */
void getResidueList(vector< vector<int> >& sse, int ssid, list<int>& residues);

int getSSLength(vector< vector<int> >& sse, int ssid);

void printSSE(vector< vector<int> >& sse);

#endif /* SSTYPE_H_ */
