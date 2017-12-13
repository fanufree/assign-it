/*
 * XplorFormat.h
 *
 *  Created on: 2013-08-06
 *      Author: e4k2
 */

#ifndef XPLORFORMAT_H_
#define XPLORFORMAT_H_

#include <tr1/unordered_map>
#include <string>

using namespace std;

class Residue;
class Atom;

class XplorFormat
{
public:
	static tr1::unordered_map< int, tr1::unordered_map<string,string> > pdb2psf; // key=Residue::AA1 or AA3
	                           // pdb atom name -> xplor atom name; for psf2pdb, can simply reverse this map
	static void init(); // sets up pdb2psf
	static string translateProton(Residue* res, int resNumOffset, bool isMethyl, Atom* hx); // non-stereospecific version, e.g. LEU HE2# instead of HE21, HE22
	                                    // returns the new name for atom hx, which is assumed to be in res
	                                // if res->num+resNumOffset == 1, terminal amide proton will be returned if hx backbone amide proton
private:
	static bool initFlag;
};

#endif /* XPLORFORMAT_H_ */
