/*
 * Contact.h
 *
 *  Created on: 2012-05-03
 *      Author: e4k2
 */

#ifndef CONTACT_H_
#define CONTACT_H_

#include <tr1/unordered_map>
#include "Residue.h"
#include "Atom.h"
#include <string>

using namespace std;

/**
 * A directed contact from hx1 to hx2
 */
class Contact
{
public:
	static const int NUMCONTACTTYPES = 17; // includes NONE
	enum ResidueContactType {NN, CACA, CBCB, SCSC, NCA, CAN,
		NCB, CBN, NSC, SCN, CACB, CBCA, CASC, SCCA, CBSC, SCCB, NONE};  // from hx1 -> hx2
	static const char* CONTACTTYPESTR[]; // same order as ResidueContactType
	Residue* r1;
	Atom* x1;
	Atom* hx1;
	Residue* r2;
	Atom* x2;
	Atom* hx2;
	ResidueContactType type;
	Contact();
	Contact(Residue* r1, Atom* x1, Atom* hx1, Residue* r2, Atom* x2, Atom* hx2);
	Contact(const Contact& c);
	// Contact(Contact&& c);
	friend void swap(Contact& first, Contact& second);
	Contact& operator=(Contact c);
	void print() const;
	virtual ~Contact();
	bool operator==(const Contact& c) const; // coordinates and chemical shifts are not compared, only the resnum and atomname
	bool operator!=(const Contact& c) const;
	Contact reverse() const; // returns hx2 -> hx1
	void setReverse(); // reverses this contact
	string typeToString() const; // if type is NONE, the empty string "" is returned
	ResidueContactType getReverseType() const;
	static ResidueContactType getReverseType(ResidueContactType type);
	static ResidueContactType getType(string& atom1, string& atom2); // for SC, atom1, atom2 should have name "SC"
private:
	void setType();
};

namespace std
{
	namespace tr1
	{
		template <>
		class hash<Contact>
		{
			public:
			size_t operator()(const Contact& c) const
			{
				return hash<int>()(c.r1->num) ^ hash<string>()(c.x1->name) ^ hash<string>()(c.hx1->name) ^
						hash<int>()(c.r2->num) ^ hash<string>()(c.x2->name) ^ hash<string>()(c.hx2->name);
			}
		};
	}
}

#endif /* CONTACT_H_ */
