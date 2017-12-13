/*
 * Residue.cpp
 *
 *  Created on: 2012-04-27
 *      Author: e4k2
 */

#include <stdio.h>
#include <string>
#include <string.h>
#include <algorithm>
#include "Residue.h"
#include "Atom.h"
#include "Utilities.h"

Residue::Residue(int resnum, int aaType) : num(resnum), type(aaType),
												   cMap(), nMap()
{
	// hList, numProtons set below
	switch(aaType)
	{
	case ALA:
		makeALA();
		break;
	case ARG:
		makeARG();
		break;
	case ASN:
		makeASN();
		break;
	case ASP:
		makeASP();
		break;
	case CYS:
		makeCYS();
		break;
	case GLN:
		makeGLN();
		break;
	case GLU:
		makeGLU();
		break;
	case GLY:
		makeGLY();
		break;
	case HIS:
		makeHIS();
		break;
	case ILE:
		makeILE();
		break;
	case LEU:
		makeLEU();
		break;
	case LYS:
		makeLYS();
		break;
	case MET:
		makeMET();
		break;
	case PHE:
		makePHE();
		break;
	case PRO:
		makePRO();
		break;
	case SER:
		makeSER();
		break;
	case THR:
		makeTHR();
		break;
	case TRP:
		makeTRP();
		break;
	case TYR:
		makeTYR();
		break;
	case VAL:
		makeVAL();
		break;
	default:
		printf("Invalid Amino Acid Type\n");
		break;
	}
}

// shallow copy
Residue::Residue(const Residue& res) : num(res.num), type(res.type),
		 cMap(res.cMap), nMap(res.nMap), hList(res.hList),numProtons(res.numProtons)
{
}

void swap(Residue& first, Residue& second)
{
	using std::swap;
	swap(first.num,second.num);
	swap(first.type,second.type);
	swap(first.cMap,second.cMap);
	swap(first.nMap,second.nMap);
	swap(first.hList,second.hList);
	swap(first.numProtons,second.numProtons);
}

//Residue::Residue(Residue&& res) : num(0), type(res.type),
//                    cMap(), nMap(), hList(NULL), numProtons(0)
//{
//	swap(*this,res);
//}

Residue::~Residue()
{
	int numP = 0;
	for (map<string,Atom*>::iterator it = cMap.begin(); it != cMap.end(); it++)
	{
		Atom* a = it->second;
		numP += a->numProtons;
		delete a; // child protons will get deleted by Atom destructor
	}
	for (map<string,Atom*>::iterator it = nMap.begin(); it != nMap.end(); it++)
	{
		Atom* a = it->second;
		numP += a->numProtons;
		delete a;
	}
	// might have protons at the end of hList, where the protons don't have
	// their heavy atom in cMap or nMap
	int diff = numProtons-numP;
	if (diff > 0)
	{
		int index = numProtons-1;
		while (diff > 0)
		{
			delete hList[index];
			index--;
			diff--;
		}
	}
}

Residue& Residue::operator=(Residue res)
{
	swap(*this,res);
	return *this;
}

void Residue::print()
{
	using namespace std;
	printf("%d %s\n",num,AA3S[type]);
	printf("cMap:\n");
	for (map<string,Atom*>::iterator it = cMap.begin(); it != cMap.end(); it++)
	{
		Atom* ap = it->second;
		ap->print();
	}
	printf("nMap:\n");
	for (map<string,Atom*>::iterator it = nMap.begin(); it != nMap.end(); it++)
	{
		Atom* ap = it->second;
		ap->print();
	}
	printf("hList:\n");
	for (int i = 0; i < numProtons; i++)
		hList[i]->print();
}

const char* Residue::getAAType1() const
{
	return AA1S[type];
}

const char* Residue::getAAType3() const
{
	return AA3S[type];
}

void Residue::addShift(const string& atomname, double shift)
{
	if (starts_with(atomname,"H"))
	{
		bool prefixMatch = false;
		bool match = false;
		for (int i = 0; i < numProtons; i++)
		{
			Atom* ap = hList[i];
			if (ap->name == atomname)
			{
				ap->setShift(shift);
				match = true;
				break;
			}
			else if (ap->name.substr(0,ap->name.length()-1) == atomname)
			{
				prefixMatch = true;
			}
		}
		if (!match) // H might be methyl
		{
			if (prefixMatch)
			{
				string m1 = atomname+"1";
				string m2 = atomname+"2";
				string m3 = atomname+"3";
				for (int i = 0; i < numProtons; i++)
				{
					Atom* ap = hList[i];
					if (ap->name == m1 || ap->name == m2 || ap->name == m3)
					{
						ap->setShift(shift);
					}
				}
			}
		}
	}
	else if (starts_with(atomname,"C"))
	{
		map<string,Atom*>::iterator it = cMap.find(atomname);
		if (it != cMap.end())
		{
			Atom* ap = it->second;
			ap->setShift(shift);
		}
		else
			printf("addShift: Atom %s not found\n",atomname.c_str());
	}
	else if (starts_with(atomname,"N"))
	{
		map<string,Atom*>::iterator it = nMap.find(atomname);
		if (it != nMap.end())
		{
			Atom* ap = it->second;
			ap->setShift(shift);
		}
		else
			printf("addShift: Atom %s not found\n",atomname.c_str());
	}
	else
	{
		printf("addShift: Skipping non-standard atom type %s\n",atomname.c_str());
	}
}

void Residue::addCoords(const string& atomname, double x, double y, double z)
{
	if (starts_with(atomname,"H"))
	{
		bool match = false;
		for (int i = 0; i < numProtons; i++)
		{
			Atom* ap = hList[i];
			if (ap->name == atomname)
			{
				ap->setCoords(x,y,z);
				match = true;
				break;
			}
		}
		if (!match)
		{
			// printf("addCoords: Atom %s not found for res %s\n",atomname.c_str(),AA3S[type]);
		}
	}
	else if (starts_with(atomname,"C"))
	{
		map<string,Atom*>::iterator it = cMap.find(atomname);
		if (it != cMap.end())
		{
			Atom* ap = it->second;
			ap->setCoords(x,y,z);
		}
		// else
		//	printf("addCoords: Atom %s not found for res %s\n",atomname.c_str(),AA3S[type]);
	}
	else if (starts_with(atomname,"N"))
	{
		map<string,Atom*>::iterator it = nMap.find(atomname);
		if (it != nMap.end())
		{
			Atom* ap = it->second;
			ap->setCoords(x,y,z);
		}
		// else
		//	printf("addCoords: Atom %s not found for res %s\n",atomname.c_str(),AA3S[type]);
	}
}

Atom* Residue::getC(const string& name) const
{
	map<string,Atom*>::const_iterator it = cMap.find(name);
	if (it != cMap.end())
		return it->second;
	else
		return NULL;
}

Atom* Residue::getN(const string& name) const
{
	map<string,Atom*>::const_iterator it = nMap.find(name);
	if (it != nMap.end())
		return it->second;
	else
		return NULL;
}

Atom* Residue::getX(const string& name) const
{
	if (starts_with(name, "C"))
	{
		return getC(name);
	}
	else if (starts_with(name,"N"))
	{
		return getN(name);
	}
	else
		return NULL;
}

int Residue::getNumProtons() const
{
	return numProtons;
}

int Residue::getNumProtonsMethylOnce() const
{
	if (type == LEU || type == ILE || type == VAL)
		return max<int>(numProtons-4, 0); // to prevent negative return values
	else if (type == ALA || type == LYS || type == THR || type == MET)
		return max<int>(numProtons-2, 0);
	else
		return numProtons;
}

AtomIterator::AtomIterator(map<string,Atom*>::const_iterator iterC,
		map<string,Atom*>::const_iterator iterN, const Residue* res) : itC(iterC), itN(iterN), r(res)
{
}

AtomIterator::~AtomIterator()
{
}

AtomIterator::AtomIterator(const AtomIterator& iter) : itC(iter.itC), itN(iter.itN), r(iter.r)
{
}

void swap(AtomIterator& first, AtomIterator& second)
{
	using std::swap;
	swap(first.itC,second.itC);
	swap(first.itN,second.itN);
	swap(first.r,second.r);
}

//AtomIterator::AtomIterator(AtomIterator&& iter) : itC(iter.r->cMap.end()), itN(iter.r->nMap.end()), r(iter.r)
//{
//	swap(*this,iter);
//}

AtomIterator& AtomIterator::operator=(AtomIterator iter)
{
	swap(*this,iter);
	return *this;
}

AtomIterator& AtomIterator::operator++()
{
	if (itC != r->cMap.end())
		itC++;
	else if (itN != r->nMap.end())
		itN++;
	return *this;
}

AtomIterator AtomIterator::operator++(int)
{
	AtomIterator ait = *this;
	++*this;
	return ait;
}

bool AtomIterator::operator==(const AtomIterator& iter) const
{
	return (itC==iter.itC) && (itN==iter.itN);
}

bool AtomIterator::operator!=(const AtomIterator& iter) const
{
	return !(*this==iter);
}

Atom* AtomIterator::operator*()
{
	if (itC != r->cMap.end())
		return itC->second;
	else if (itN != r->nMap.end())
		return itN->second;
	else
		return NULL;
}

AtomIterator Residue::begin(char type) const
{
	if (type == 'C')
		return AtomIterator(cMap.begin(),nMap.end(),this);
	else if (type == 'N')
		return AtomIterator(cMap.end(),nMap.begin(),this);
	else
		return AtomIterator(cMap.begin(),nMap.begin(),this);
}

AtomIterator Residue::end(char type) const
{
	if (type == 'C')
		return AtomIterator(cMap.end(),nMap.end(),this);
	else if (type == 'N')
		return AtomIterator(cMap.end(),nMap.end(),this);
	else
		return AtomIterator(cMap.end(),nMap.end(),this);
}

void Residue::getSideChainCentroid(double& x, double& y, double& z) const
{
	x = 0;
	y = 0;
	z = 0;
	int count = 0;
	for (AtomIterator it = begin('X'); it != end('X'); ++it)
	{
		Atom* atom = *it;
		string& name = atom->name;
		if (name != "N" && name != "C" && !starts_with(name,"CA") && name != "O")
		{
			if (atom->x != INVALIDCOORD)
			{
				x += atom->x;
				y += atom->y;
				z += atom->z;
				count++;
			}
		}
		else if (type == GLY && starts_with(name,"CA"))// for GLY count the CA
		{
			if (atom->x != INVALIDCOORD)
			{
				x += atom->x;
				y += atom->y;
				z += atom->z;
				count++;
			}
		}
	}
//	for (int i = 0; i < numProtons; i++)
//	{
//		Atom* atom = hList[i];
//		string& name = atom->name;
//		if (name != "H" && !starts_with(name,"HA"))
//		{
//			if (atom->x != INVALIDCOORD)
//		    {
//				x += atom->x;
//				y += atom->y;
//				z += atom->z;
//				count++;
//		    }
//		}
//	}
	if (count > 0)
	{
		x = x/double(count);
		y = y/double(count);
		z = z/double(count);
	}
	else
	{
		x = y = z = INVALIDCOORD;
	}
}

bool Residue::isCaro(const string& catomname) const
{
	if (type == Residue::HIS)
	{
		if (catomname == "CD2")
			return true;
		if (catomname == "CE1")
			return true;
	}
	else if (type == Residue::PHE)
	{
		if (catomname == "CD1")
			return true;
		if (catomname == "CD2")
			return true;
		if (catomname == "CE1")
			return true;
		if (catomname == "CE2")
			return true;
		if (catomname == "CZ")
			return true;
	}
	else if (type == Residue::TRP)
	{
		if (catomname == "CD1")
			return true;
		if (catomname == "CE3")
			return true;
		if (catomname == "CZ2")
			return true;
		if (catomname == "CZ3")
			return true;
		if (catomname == "CH2")
			return true;
	}
	else if (type == Residue::TYR)
	{
		if (catomname == "CD1")
			return true;
		if (catomname == "CD2")
			return true;
		if (catomname == "CE1")
			return true;
		if (catomname == "CE2")
			return true;
	}
	return false;
}

bool Residue::isCali(const string& catomname) const
{
	if (starts_with(catomname,"C"))
	{
		if (!isCaro(catomname))
			return true;
	}
	return false;
}

bool Residue::is15N(const string& natomname) const
{
	if (starts_with(natomname,"N"))
		return true;
	return false;
}

const char* Residue::AA3S[] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
const char* Residue::AA1S[] = {"A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"};

bool Residue::isHydrophobic() const
{
	if (type == ALA || type == LEU || type == ILE || type == MET || type == PHE ||
		type == VAL || type == PRO || type == GLY)
		return true;
	else
		return false;
}

bool Residue::isPositive() const
{
	if (type == ARG || type == HIS || type == LYS)
		return true;
	else
		return false;
}

bool Residue::isNegative() const
{
	if (type == ASP || type == GLU)
		return true;
	else
		return false;
}

bool Residue::isCYS() const
{
	if (type == CYS)
		return true;
	else
		return false;
}

// ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE ,LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL
int Residue::index3(const char* aaName3)
{
	if (aaName3[0] < 'L')
	{
		if (aaName3[0] < 'C')
		{
			if (aaName3[1] < 'S')
			{
				if (aaName3[1] == 'L')
					return ALA;
				else
					return ARG;
			}
			else if (aaName3[2] == 'N')
				return ASN;
			else
				return ASP;
		}
		else if (aaName3[0] > 'C')
		{
			if (aaName3[0] < 'H')
			{
				if (aaName3[2] <'Y')
				{
					if (aaName3[2] == 'N')
						return GLN;
					else
						return GLU;
				}
				else
					return GLY;
			}
			else if (aaName3[0] == 'I')
				return ILE;
			else
				return HIS;
		}
		else
			return CYS;
	}
	else if (aaName3[0] > 'L')
	{
		if (aaName3[0] < 'T')
		{
			if (aaName3[0] != 'P')
			{
				if (aaName3[0] == 'S')
					return SER;
				else
					return MET;
			}
			else if (aaName3[1] == 'R')
					return PRO;
			else
				return PHE;
		}
		else
		{
			if (aaName3[0] == 'V')
				return VAL;
			if (aaName3[1] < 'R')
				return THR;
			else if (aaName3[1] == 'Y')
				return TYR;
			else
				return TRP;
		}
	}
	else if (aaName3[1] == 'E')
		return LEU;
	else
		return LYS;

}

// A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
int Residue::index1(const char* aaName1)
{
	if (aaName1[0] < 'L')
	{
		if (aaName1[0] < 'F')
		{
			if (aaName1[0] < 'D')
			{
				if (aaName1[0] == 'A')
					return A;
				else
					return C;
			}
			else if (aaName1[0] > 'D')
			{
				return E;
			}
			else
				return D;
		}
		else if (aaName1[0] > 'F')
		{
			if (aaName1[0] < 'I')
			{
				if (aaName1[0] == 'G')
					return G;
				else
					return H;
			}
			else if (aaName1[0] > 'I')
			{
				return K;
			}
			else
				return I;
		}
		else
			return F;
	}
	else if (aaName1[0] > 'L')
	{
		if (aaName1[0] < 'R')
		{
			if (aaName1[0] < 'N')
			{
				return M;
			}
			else if (aaName1[0] > 'N')
			{
				if (aaName1[0] == 'P')
					return P;
				else
					return Q;
			}
			else
				return N;
		}
		else if (aaName1[0] > 'R')
		{
			if (aaName1[0] < 'V')
			{
				if (aaName1[0] == 'T')
					return T;
				else
					return S;
			}
			else if (aaName1[0] > 'V')
			{
				if (aaName1[0] == 'Y')
					return Y;
				else
					return W;
			}
			else
				return V;
		}
		else
			return R;
	}
	else
		return L;
}

// A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
int Residue::index1(char aaName1)
{
	if (aaName1 < 'L')
	{
		if (aaName1 < 'F')
		{
			if (aaName1 < 'D')
			{
				if (aaName1 == 'A')
					return A;
				else
					return C;
			}
			else if (aaName1 > 'D')
			{
				return E;
			}
			else
				return D;
		}
		else if (aaName1 > 'F')
		{
			if (aaName1 < 'I')
			{
				if (aaName1 == 'G')
					return G;
				else
					return H;
			}
			else if (aaName1 > 'I')
			{
				return K;
			}
			else
				return I;
		}
		else
			return F;
	}
	else if (aaName1 > 'L')
	{
		if (aaName1 < 'R')
		{
			if (aaName1 < 'N')
			{
				return M;
			}
			else if (aaName1 > 'N')
			{
				if (aaName1 == 'P')
					return P;
				else
					return Q;
			}
			else
				return N;
		}
		else if (aaName1 > 'R')
		{
			if (aaName1 < 'V')
			{
				if (aaName1 == 'T')
					return T;
				else
					return S;
			}
			else if (aaName1 > 'V')
			{
				if (aaName1 == 'Y')
					return Y;
				else
					return W;
			}
			else
				return V;
		}
		else
			return R;
	}
	else
		return L;
}

bool Residue::isStandard3(const char* aaName3)
{
	for (int i = 0; i < 20; i++)
	{
		if (strcmp(AA3S[i],aaName3) == 0)
			return true;
	}
	return false;
}

bool Residue::isStandard1(const char* aaName1)
{
	for (int i = 0; i < 20; i++)
	{
		if (strcmp(AA1S[i],aaName1) == 0)
			return true;
	}
	return false;
}

bool Residue::isStandard1(char aaName1)
{
	for (int i = 0; i < 20; i++)
	{
		if (AA1S[i][0] == aaName1)
			return true;
	}
	return false;
}

void Residue::makeBackbone(int& hIndex)
{
	addCN(false,"N",H1,hIndex);
	addCN(true,"CA",H1,hIndex);
	addCN(true,"C",H0,hIndex);
	// O has no chem shift
}

void Residue::makeBackboneHB2(int& hIndex)
{
	makeBackbone(hIndex);
	addCN(true,"CB",H2, hIndex);
}

void Residue::addCN(bool cFlag, const string& name, ADDTYPEH addH, int& hIndex)
{
	using namespace std;
	if (cFlag)
	{
		if (addH == H2)
		{
			Atom* cp = new Atom(name,2);
			cp->addH2();
			for (int i = 0; i < 2; i++)
				hList[hIndex++]=cp->protons[i];
			cMap[name] = cp;
		}
		else if (addH == H3)
		{
			Atom* cp = new Atom(name,3);
			cp->addH3();
			for (int i = 0; i < 3; i++)
				hList[hIndex++]=cp->protons[i];
			cMap[name] = cp;
		}
		else if (addH == H1)
		{
			Atom* cp = new Atom(name,1);
			cp->addH1();
		    hList[hIndex++]=cp->protons[0];
			cMap[name] = cp;
		}
		else
		{
			Atom* cp = new Atom(name,0);
			cMap[name] = cp;
		}
	}
	else
	{
		if (addH == H1)
		{
			Atom* np = new Atom(name,1);
			np->addH1();
			hList[hIndex++]=np->protons[0];
			nMap[name] = np;
		}
		else if (addH == H2)
		{
			Atom* np = new Atom(name,2);
			np->addH2_N();
			for (int i = 0; i < 2; i++)
				hList[hIndex++]=np->protons[i];
			nMap[name] = np;
		}
		else if (addH == H0)
		{
			Atom* np = new Atom(name,0);
			nMap[name] = np;
		}
		else
		{
			Atom* np = new Atom(name,3);
			np->addH3();
			for (int i = 0; i < 3; i++)
				hList[hIndex++]=np->protons[i];
			nMap[name] = np;
		}
	}
}

void Residue::makeALA()
{
	type = ALA;
	numProtons = 5;
	hList = new Atom*[5];
	int index = 0;
	makeBackbone(index);
	addCN(true,"CB",H3,index);
}

void Residue::makeARG()
{
	type = ARG;
	numProtons = 9;
	hList = new Atom*[9];
	int index = 0;
	makeBackboneHB2(index);
	addCN(true,"CG",H2,index);
	addCN(true,"CD",H2,index);
	// chem shifts of NE, NH1, NH2, and corresponding H's usually not visible
	// among the 3, HE most visible
	addCN(false,"NE",H1,index);
}

void Residue::makeASN()
{
	type = ASN;
	numProtons = 6;
	hList = new Atom*[6];
	int index = 0;
	makeBackboneHB2(index);
	addCN(true,"CG",H0,index); // no attached H's
	addCN(false,"ND2",H2,index);
}

void Residue::makeASP()
{
	type = ASP;
	numProtons = 4;
	hList = new Atom*[4];
	int index = 0;
	makeBackboneHB2(index);
	// CG rare (<1000 in BMRB stats) & has no attached H's
	// HD2 on O is rare
}

void Residue::makeCYS()
{
	type = CYS;
	numProtons = 4;
	hList = new Atom*[4];
	int index = 0;
	makeBackboneHB2(index);
	// HG on S is rare
}

void Residue::makeGLN()
{
	type = GLN;
	numProtons = 8;
	hList = new Atom*[8];
	int index = 0;
	makeBackboneHB2(index);
	addCN(true,"CG",H2,index);
	addCN(true,"CD",H0,index); // no attached H's
	addCN(false,"NE2",H2,index);
}

void Residue::makeGLU()
{
	type = GLU;
	numProtons = 6;
	hList = new Atom*[6];
	int index = 0;
	makeBackboneHB2(index);
	addCN(true,"CG",H2,index);
	// CD rare, HE2 on OE1 rare
}

void Residue::makeGLY()
{
	type = GLY;
	numProtons = 3;
	hList = new Atom*[3];
	int index = 0;
	addCN(false,"N",H1,index);
	addCN(true,"CA",H2,index);
	addCN(true,"C",H0,index);
}

void Residue::makeHIS()
{
	type = HIS;
	numProtons = 8;
	hList = new Atom*[8];
	int index = 0;
	makeBackboneHB2(index);
	// CG rare
	addCN(true,"CD2",H1,index);
	addCN(true,"CE1",H1,index);
	// ND1, NE2 usually not visible
	// but add attached H's
	addCN(false,"ND1",H1,index);
	addCN(false,"NE2",H1,index);
}

void Residue::makeILE()
{
	type = ILE;
	numProtons = 11;
	hList = new Atom*[11];
	int index = 0;
	makeBackbone(index);
	addCN(true,"CB",H1,index);
	addCN(true,"CG1",H2,index);
	addCN(true,"CG2",H3,index);
	addCN(true,"CD1",H3,index);
}

void Residue::makeLEU()
{
	type = LEU;
	numProtons = 11;
	hList = new Atom*[11];
	int index = 0;
	makeBackboneHB2(index);
	addCN(true,"CG",H1,index);
	addCN(true,"CD1",H3,index);
	addCN(true,"CD2",H3,index);
}

void Residue::makeLYS()
{
	type = LYS;
	numProtons = 13;
	hList = new Atom*[13];
	int index = 0;
	makeBackboneHB2(index);
	addCN(true,"CG",H2,index);
	addCN(true,"CD",H2,index);
	addCN(true,"CE",H2,index);
	// HZ rare but could be useful if present, NZ usually not visible
	addCN(false,"NZ",H3,index);
}

void Residue::makeMET()
{
	type = MET;
	numProtons = 9;
	hList = new Atom*[9];
	int index = 0;
	makeBackboneHB2(index);
	addCN(true,"CG",H2,index);
	addCN(true,"CE",H3,index);
}

void Residue::makePHE()
{
	type = PHE;
	numProtons = 9;
	hList = new Atom*[9];
	int index = 0;
	makeBackboneHB2(index);
	// CG rare, no attached H
	addCN(true,"CD1",H1,index);
	addCN(true,"CD2",H1,index);
	addCN(true,"CE1",H1,index);
	addCN(true,"CE2",H1,index);
	addCN(true,"CZ",H1,index);
}

void Residue::makePRO()
{
	type = PRO;
	numProtons = 7;
	hList = new Atom*[7];
	int index = 0;
	// N chem shift usually not available
	addCN(true,"CA",H1,index);
	addCN(true,"C",H0,index);
	addCN(true,"CB",H2,index);
	addCN(true,"CG",H2,index);
	addCN(true,"CD",H2,index);
}

void Residue::makeSER()
{
	type = SER;
	numProtons = 4;
	hList = new Atom*[4];
	int index = 0;
	makeBackboneHB2(index);
	// HG on OG rare
}

void Residue::makeTHR()
{
	type = THR;
	numProtons = 7;
	hList = new Atom*[7];
	int index = 0;
	makeBackbone(index);
	addCN(true,"CB",H1,index);
	addCN(true,"CG2",H3,index);
	addCN(true,"OG1",H1,index);// treat as O as C; HG1 on OG1 rare, but more common than HG in SER
}

void Residue::makeTRP()
{
	type = TRP;
	numProtons = 10;
	hList = new Atom*[10];
	int index = 0;
	makeBackboneHB2(index);
	// CG rare, no attached H
	addCN(true,"CD1",H1,index);
	// CD2 rare, no attached H
	addCN(false,"NE1",H1,index);
	// CE2 rare, no attached H
	addCN(true,"CE3",H1,index);
	addCN(true,"CZ2",H1,index);
	addCN(true,"CZ3",H1,index);
	addCN(true,"CH2",H1,index);
}

void Residue::makeTYR()
{
	type = TYR;
	numProtons = 8;
	hList = new Atom*[8];
	int index = 0;
	makeBackboneHB2(index);
	// CG rare, no attached H
	addCN(true,"CD1",H1,index);
	addCN(true,"CD2",H1,index);
	addCN(true,"CE1",H1,index);
	addCN(true,"CE2",H1,index);
	// CZ rare, no attached H
	// HH on OH rare
}

void Residue::makeVAL()
{
	type = VAL;
	numProtons = 9;
	hList = new Atom*[9];
	int index = 0;
	makeBackbone(index);
	addCN(true,"CB",H1,index);
	addCN(true,"CG1",H3,index);
	addCN(true,"CG2",H3,index);
}
