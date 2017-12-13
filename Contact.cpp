/*
 * Contact.cpp
 *
 *  Created on: 2012-05-03
 *      Author: e4k2
 */

#include <stdio.h>
#include "Contact.h"
#include "Residue.h"
#include "Atom.h"
#include "Utilities.h"

using namespace std;

const char* Contact::CONTACTTYPESTR[] = {"NN","CACA","CBCB","SCSC","NCA","CAN",
                    "NCB","CBN","NSC","SCN","CACB","CBCA","CASC","SCCA","CBSC","SCCB","NONE"};

Contact::Contact() : r1(NULL), x1(NULL), hx1(NULL), r2(NULL), x2(NULL), hx2(NULL), type(NONE)
{
}

Contact::~Contact()
{
	// pointers to be deleted externally
}

Contact::Contact(Residue* r1, Atom* x1, Atom* hx1, Residue* r2, Atom* x2, Atom* hx2) :
		r1(r1), x1(x1), hx1(hx1), r2(r2), x2(x2), hx2(hx2)
{
	setType();
}

Contact::Contact(const Contact& c) : r1(c.r1), x1(c.x1), hx1(c.hx1), r2(c.r2), x2(c.x2), hx2(c.hx2)
{
	type = c.type;
}

void swap(Contact& first, Contact& second)
{
	using std::swap;
	swap(first.r1,second.r1);
	swap(first.x1,second.x1);
	swap(first.hx1,second.hx1);
	swap(first.r2,second.r2);
	swap(first.x2,second.x2);
	swap(first.hx2,second.hx2);
	swap(first.type,second.type);
}

//Contact::Contact(Contact&& c) : r1(NULL), x1(NULL), hx1(NULL), r2(NULL), x2(NULL), hx2(NULL)
//{
//	swap(*this,c);
//}

Contact& Contact::operator=(Contact c)
{
	swap(*this,c);
	return *this;
}

bool Contact::operator==(const Contact& c) const
{
	//return hx1==c.hx1 && hx2==c.hx2;
	return r1->num==c.r1->num && x1->name==c.x1->name && hx1->name==c.hx1->name &&
		   r2->num==c.r2->num && x2->name==c.x2->name && hx2->name==c.hx2->name;
}

bool Contact::operator!=(const Contact& c) const
{
	return !(*this == c);
}

void Contact::print() const
{
	printf("%d %s %s -> %d %s %s\n",r1->num,x1->name.c_str(),hx1->name.c_str(),
			r2->num,x2->name.c_str(),hx2->name.c_str());
}

Contact Contact::reverse() const
{
	return Contact(r2,x2,hx2,r1,x1,hx1);
}

void Contact::setReverse()
{
	Residue* tempR = r1;
	r1 = r2;
	r2 = tempR;
	Atom* tempX = x1;
	x1 = x2;
	x2 = tempX;
	Atom* tempHX = hx1;
	hx1 = hx2;
	hx2 = tempHX;
}

string Contact::typeToString() const
{
	string& h1 = hx1->name;
	string& h2 = hx2->name;
	if (h1 == "H")
	{
		if (h2 == "H")
			return "NN";
		else if (starts_with(h2,"HA"))
			return "NCA";
		else if (starts_with(h2,"HB"))
			return "NCB";
		else
			return "NSC";
	}
	else if (starts_with(h1,"HA"))
	{
		if (h2 == "H")
			return "CAN";
		else if (starts_with(h2,"HA"))
			return "CACA";
		else if (starts_with(h2,"HB"))
			return "CACB";
		else
			return "CASC";
	}
	else if (starts_with(h1,"HB"))
	{
		if (h2 == "H")
			return "CBN";
		else if (starts_with(h2,"HA"))
			return "CBCA";
		else if (starts_with(h2,"HB"))
			return "CBCB";
		else
			return "CBSC";
	}
	else
	{
		if (h2 == "H")
			return "SCN";
		else if (starts_with(h2,"HA"))
			return "SCCA";
		else if (starts_with(h2,"HB"))
			return "SCCB";
		else
			return "SCSC";
	}
	return "";
}

Contact::ResidueContactType Contact::getReverseType() const
{
	if (type == NCA)
		return CAN;
	else if (type == CAN)
		return NCA;
	else if (type == NCB)
		return CBN;
	else if (type == CBN)
		return NCB;
	else if (type == NSC)
		return SCN;
	else if (type == SCN)
		return NSC;
	else if (type == CACB)
		return CBCA;
	else if (type == CBCA)
		return CACB;
	else if (type == CASC)
		return SCCA;
	else if (type == SCCA)
		return CASC;
	else if (type == CBSC)
		return SCCB;
	else if (type == SCCB)
		return CBSC;
	return type;
}

Contact::ResidueContactType Contact::getReverseType(ResidueContactType type)
{
	if (type == NCA)
		return CAN;
	else if (type == CAN)
		return NCA;
	else if (type == NCB)
		return CBN;
	else if (type == CBN)
		return NCB;
	else if (type == NSC)
		return SCN;
	else if (type == SCN)
		return NSC;
	else if (type == CACB)
		return CBCA;
	else if (type == CBCA)
		return CACB;
	else if (type == CASC)
		return SCCA;
	else if (type == SCCA)
		return CASC;
	else if (type == CBSC)
		return SCCB;
	else if (type == SCCB)
		return CBSC;
	return type;
}

Contact::ResidueContactType Contact::getType(string& atom1, string& atom2)
{
	if (atom1 == "SC")
	{
		if (atom2 == "SC")
		{
			return SCSC;
		}
		else if (atom2 == "CB")
		{
			return SCCB;
		}
		else if (atom2 == "CA")
		{
			return SCCA;
		}
		else if (atom2 == "N")
		{
			return SCN;
		}
		else
			return NONE;
	}
	else if (atom1 == "CB")
	{
		if (atom2 == "SC")
		{
			return CBSC;
		}
		else if (atom2 == "CB")
		{
			return CBCB;
		}
		else if (atom2 == "CA")
		{
			return CBCA;
		}
		else if (atom2 == "N")
		{
			return CBN;
		}
		else
			return NONE;
	}
	else if (atom1 == "CA")
	{
		if (atom2 == "SC")
		{
			return CASC;
		}
		else if (atom2 == "CB")
		{
			return CACB;
		}
		else if (atom2 == "CA")
		{
			return CACA;
		}
		else if (atom2 == "N")
		{
			return CAN;
		}
		else
			return NONE;
	}
	else if (atom1 == "N")
	{
		if (atom2 == "SC")
		{
			return NSC;
		}
		else if (atom2 == "CB")
		{
			return NCB;
		}
		else if (atom2 == "CA")
		{
			return NCA;
		}
		else if (atom2 == "N")
		{
			return NN;
		}
		else
			return NONE;
	}
	else
		return NONE;
}

void Contact::setType()
{
	if (hx1 != NULL && hx2 != NULL)
	{
		string& h1 = hx1->name;
		string& h2 = hx2->name;
		if (h1 == "H")
		{
			if (h2 == "H")
				type = NN;
			else if (starts_with(h2,"HA"))
				type = NCA;
			else if (starts_with(h2,"HB"))
				type = NCB;
			else
				type = NSC;
		}
		else if (starts_with(h1,"HA"))
		{
			if (h2 == "H")
				type = CAN;
			else if (starts_with(h2,"HA"))
				type = CACA;
			else if (starts_with(h2,"HB"))
				type = CACB;
			else
				type = CASC;
		}
		else if (starts_with(h1,"HB"))
		{
			if (h2 == "H")
				type = CBN;
			else if (starts_with(h2,"HA"))
				type = CBCA;
			else if (starts_with(h2,"HB"))
				type = CBCB;
			else
				type = CBSC;
		}
		else
		{
			if (h2 == "H")
				type = SCN;
			else if (starts_with(h2,"HA"))
				type = SCCA;
			else if (starts_with(h2,"HB"))
				type = SCCB;
			else
				type = SCSC;
		}
	}
}
