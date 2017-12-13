/*
 * Atom.cpp
 *
 *  Created on: 2012-04-29
 *      Author: e4k2
 */

#include <cstdio>
#include <cmath>
#include "Atom.h"
#include "Utilities.h"

Atom::Atom() : name(), x(INVALIDCOORD), y(INVALIDCOORD), z(INVALIDCOORD), cs(INVALIDSHIFTVAL), protons(NULL), numProtons(0)
{
}

Atom::Atom(const string& atomName,int numP) : name(atomName), x(INVALIDCOORD), y(INVALIDCOORD), z(INVALIDCOORD), cs(INVALIDSHIFTVAL), numProtons(numP)
{
	if (numP > 0)
		protons = new Atom*[numP];
	else
		protons = NULL;
}

// shallow copy
Atom::Atom(const Atom& atom) : name(atom.name), x(atom.x), y(atom.y), z(atom.z), cs(atom.cs), protons(atom.protons), numProtons(atom.numProtons)
{
}

void swap(Atom& first, Atom& second)
{
	using std::swap;
	swap(first.name,second.name);
	swap(first.x,second.x);
	swap(first.y,second.y);
	swap(first.z,second.z);
	swap(first.cs,second.cs);
	swap(first.protons,second.protons);
	swap(first.numProtons,second.numProtons);
}

//Atom::Atom(Atom&& atom) : name(), x(INVALIDCOORD), y(INVALIDCOORD), z(INVALIDCOORD), cs(INVALIDSHIFTVAL), protons(NULL), numProtons(0)
//{
//	swap(*this, atom);
//}

Atom::~Atom()
{
	if (protons)
	{
		for (int i = 0; i < numProtons; i++)
			delete protons[i];
		delete [] protons;
	}
}

Atom& Atom::operator=(Atom atom)
{
	swap(*this,atom);
	return *this;
}

Atom* Atom::operator[](const int index)
{
	if (index < numProtons)
		return protons[index];
	else
		return NULL;
}

void Atom::print()
{
	printf("%s %8.3f %8.3f %8.3f %6.2f\t",name.c_str(),x,y,z,cs);
	for (int i = 0; i < numProtons; i++)
		printf("%s ",protons[i]->name.c_str());
	printf("\n");
}

void Atom::addH3()
{
	string prefix = "H"+name.substr(1,name.length()-1);
	protons[0] = new Atom(prefix+"1",0);
	protons[1] = new Atom(prefix+"2",0);
	protons[2] = new Atom(prefix+"3",0);
}

void Atom::addH2()
{
	string prefix = "H"+name.substr(1,name.length()-1);
	protons[0] = new Atom(prefix+"2",0);
	protons[1] = new Atom(prefix+"3",0);
}

void Atom::addH1()
{
	const int len = name.length();
	string hname = "H";
	if (len > 1)
		hname = hname+name.substr(1,len-1);
	protons[0] = new Atom(hname,0);
}

void Atom::addH2_N()
{
	string prefix = "H"+name.substr(1,name.length()-1);
	protons[0] = new Atom(prefix+"1",0);
	protons[1] = new Atom(prefix+"2",0);
}

Atom* Atom::getH(const string& name) const
{
	for (int i = 0; i < numProtons; i++)
	{
		if (protons[i]->name == name)
			return protons[i];
	}
	return NULL;
}

int Atom::getNumProtons() const
{
	return numProtons;
}

void Atom::setShift(double shift)
{
	cs = shift;
}

void Atom::setCoords(double cx, double cy, double cz)
{
	x = cx;
	y = cy;
	z = cz;
}

bool Atom::isHeavy() const
{
	if (starts_with(name,"H"))
		return false;
	else
		return true;
}

bool Atom::isNMR_SideChainHeavy() const
{
	if (starts_with(name,"C") || starts_with(name,"N"))
		if (name != "N" && name != "C" && name != "CA" && name != "CB")
			return true;
		else
			return false;
	else
		return false;
}

bool Atom::isBackbone() const
{
	if (name == "H" || starts_with(name,"HA") || name == "N" || name == "CA" ||
		name == "HN" || name == "C" || name == "H1" || name == "H2" || name == "H3" ||
		name == "O" || name == "OXT")
		return true;
	else
		return false;
}

bool Atom::isProton() const
{
	if (starts_with(name,"H"))
		return true;
	else
		return false;
}

bool Atom::isMethyl() const
{
	return (numProtons == 3); // allow NH3 to be considered a methyl
}

bool Atom::isCB() const
{
	return (name == "CB");
}

bool Atom::isGeminal() const
{
	return (numProtons == 2);
}

char Atom::getType() const
{
	return name[0];
}

double Atom::getDistance(const Atom* a) const
{
	return sqrt(pow(x-a->x,2.0)+pow(y-a->y,2.0)+pow(z-a->z,2.0));
}

HIterator Atom::beginH() const
{
	return protons;
}

HIterator Atom::endH() const
{
	return protons+numProtons;
}

bool Atom::hasCoordinates() const
{
	if (x != INVALIDCOORD)
		return true;
	else
		return false;
}

bool Atom::hasChemShift() const
{
	if (cs != INVALIDSHIFTVAL)
		return true;
	else
		return false;
}

bool Atom::isProtonGroupChemShift(double tol) const
{
	if (numProtons > 0)
	{
		for (int i = 0; i < numProtons-1; i++)
		{
			Atom* a1 = protons[i];
			for (int j = i+1; j < numProtons; j++)
			{
				Atom* a2 = protons[j];
				if (abs(a1->cs-a2->cs) > tol)
					return false;
			}
		}
		return true;
	}
	return false;
}
