/*
 * CSProtein.cpp
 *
 *  Created on: 2012-04-29
 *      Author: e4k2
 */

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <list>
#include <algorithm>
#include <tr1/unordered_set>
#include "CSProtein.h"
#include "Residue.h"
#include "Atom.h"
#include "Utilities.h"
#include "PseudoContact.h"

#define INVALIDDIST 9999999 // used by inContact

using namespace std;


CSProtein::CSProtein(const string& seq, int startingResNum) : startResNum(startingResNum), size(seq.length())
{
	ifstream file;
	string line;
	if (size > MAXPROSIZE)
	{
		printf("Protein too big: size is %d\n",size);
		exit(-1);
	}
	residues = new Residue*[seq.length()];
	for (int i = 0; i < size; i++)
	{
		if (Residue::isStandard1(seq[i]))
		{
			int type = Residue::index1(seq[i]);
			residues[i] = new Residue(startingResNum+i,type);
		}
		else
		{
			printf("Nonstandard amino acid %c\n",line[i]);
			residues[i] = NULL;
		}
	}
}

// shallow copy
CSProtein::CSProtein(const CSProtein& other) : startResNum(other.startResNum), size(other.size)
{
	residues = new Residue*[other.size];
	if (other.residues != NULL)
	{
		for (int i = 0; i < other.size; i++)
			residues[i] = other.residues[i];
	}
}

CSProtein::~CSProtein()
{
	for (int i = 0; i < size; i++)
	{
		if (residues[i] != NULL)
		{
			delete residues[i];
		}
	}
	delete [] residues;

}
void swap(CSProtein& first, CSProtein& second)
{
	using std::swap;
	swap(first.startResNum,second.startResNum);
	swap(first.size,second.size);
	swap(first.residues,second.residues);
}

CSProtein& CSProtein::operator=(CSProtein other)
{
	swap(*this,other);
	return *this;
}

Residue* CSProtein::operator[](int resnum)
{
	return residues[resnum-startResNum];
}

//CSProtein::CSProtein(CSProtein&& other) : startResNum(0), size(0), residues(NULL)
//{
//	swap(*this,other);
//}

void CSProtein::print()
{
	printf("size: %d, start res num: %d\n",size,startResNum);
	for (int i = 0; i < size; i++)
	{
		residues[i]->print();
	}
}

// PDB format
// We consider only up to col 54
// 10-13 Model number (right aligned)
// 0 – 5     Record name e.g. ATOM
// 6 – 10    Integer     serial  Atom serial number.
// 12 – 15   Atom    name    Atom name.
// 16    Character   altLoc  Alternate location indicator.
// 17 – 19   Residue name    resName     Residue name.
// 21    Character   chainID     Chain identifier.
// 22 – 25   Integer     resSeq  Residue sequence number.
// 26    AChar   iCode   Code for insertion of residues.
// 30 – 37   Real(8.3)   x   Orthogonal coordinates for X.
// 38 – 45   Real(8.3)   y   Orthogonal coordinates for Y.
// 46 – 53   Real(8.3)   z   Orthogonal coordinates for Z.
// 54 – 59   Real(6.2)   occupancy   Occupancy.
// 60 – 65   Real(6.2)   tempFactor  Temperature factor.
// 76 – 77   LString(2)  element     Element symbol; right-justified.
// 78 – 79   LString(2)  charge  Charge on the atom.
void CSProtein::getCoords(const string& pdbFile, int offset)
{
	ifstream file;
	string line;
	file.open(pdbFile.c_str());
	double ntermx=0;
	double ntermy=0;
	double ntermz=0;
	int numNtermHs=0; // PRO 2, rest 3 protons
	if (file.is_open())
	{
		while (file.good())
		{
			getline(file,line);
			trim(line);
			if (line.length() < 6)
				continue;
			string recordName = line.substr(0,6);
			trim(recordName);
			if (recordName == "TER" || recordName == "END" || recordName == "ENDMDL") // read only first chain of first model
				break;
			if (line.length() < 54)
				continue;
			if (recordName != "ATOM")
				continue;
			string atomName = line.substr(12,4);
			trim(atomName);
			string resName = line.substr(17,3);
			string temp = line.substr(22,4);
			trim(temp);
			int resSeq = atoi(temp.c_str())+offset;
			if (resSeq-startResNum >= size || resSeq-startResNum < 0)
				continue;
			Residue* res = residues[resSeq-startResNum];
			if (res == NULL)
			{
				printf("getCoords: No such residue in sequence %d\n",resSeq);
				continue;
			}
			if (Residue::AA3S[res->type] != resName)
			{
				printf("WARNING: getCoords: PDB-Seq res type mismatch %s-%s at res %d\n",
						resName.c_str(),Residue::AA3S[res->type],resSeq);
			}
			temp = line.substr(31,8);
			trim(temp);
			double x = atof(temp.c_str());
			temp = line.substr(39,8);
			trim(temp);
			double y = atof(temp.c_str());
			temp = line.substr(47,8);
			trim(temp);
			double z = atof(temp.c_str());
			if (atomName == "H1" || atomName == "H2")
			{   // n-term H, take average of H1,H2,H3
				ntermx += x;
				ntermy += y;
				ntermz += z;
				numNtermHs++;
				continue;
			}
			else if (atomName == "H3")
			{
				x += ntermx;
				y += ntermy;
				z += ntermz;
				numNtermHs++;
				x = x/double(numNtermHs);
				y = y/double(numNtermHs);
				z = z/double(numNtermHs);
				atomName = "H";
				ntermx = 0;
				ntermy = 0;
				ntermz = 0;
				numNtermHs=0;
			}
			res->addCoords(atomName,x,y,z);
		}
	}
	else
	{
		cout << "Unable to open file " << pdbFile << endl;
		file.close();
		exit(-1);
	}
	file.close();
}

int CSProtein::getNumProtons() const
{
	int numProtons = 0;
	for (int r = 1; r <= size; r++)
	{
		Residue* res =  residues[r-startResNum];
		for (AtomIterator itX = res->begin('X'); itX != res->end('X'); itX++)
		{
			Atom* xAtom = *itX;
			numProtons += xAtom->getNumProtons();
		}
	}
	return numProtons;
}

void CSProtein::getModels(const string& pdbFile, list<CSProtein*>& models, const string& seq, int startingResNum, int offset)
{
	ifstream file;
	string line;
	file.open(pdbFile.c_str());
	int size = seq.length();
	double ntermx=0;
	double ntermy=0;
	double ntermz=0;
	int numNtermHs=0; // PRO 2, rest 3 protons
	if (file.is_open())
	{
		CSProtein* current = NULL;
		unsigned int numModels = 1;
		while (file.good())
		{
			getline(file,line);
			trim(line);
			if (line.length() < 6)
				continue;
			string recordName = line.substr(0,6);
			trim(recordName);

			if (recordName == "NUMMDL")
			{
				string temp = line.substr(10,13);
				trim(temp);
				numModels = atoi(temp.c_str());
				continue;
			}
			if (recordName == "MODEL")
			{
				current = new CSProtein(seq,startingResNum);
				models.push_back(current);
				continue;
			}
			if (recordName == "TER")
				continue;
			if (recordName == "END")
				break;
			if (line.length() < 54)
				continue;
			if (recordName != "ATOM")
				continue;

			if (line[21] != 'A' && line[21] != ' ')
				continue;

			if (current == NULL)
			{
				// no MODEL recordName; likely a single model pdb file
				current = new CSProtein(seq,startingResNum);
				models.push_back(current);
			}
			string atomName = line.substr(12,4);
			trim(atomName);
			string resName = line.substr(17,3);
			string temp = line.substr(22,4);
			trim(temp);
			int resSeq = atoi(temp.c_str())+offset;
			if (resSeq-startingResNum >= size || resSeq-startingResNum < 0)
				continue;
			Residue* res = current->residues[resSeq-startingResNum];
			if (res == NULL)
			{
				printf("getCoords: No such residue in sequence %d\n",resSeq);
				continue;
			}
			if (Residue::AA3S[res->type] != resName)
			{
				printf("WARNING: getCoords: PDB-Seq res type mismatch %s-%s at res %d\n",
						resName.c_str(),Residue::AA3S[res->type],resSeq);
			}
			temp = line.substr(31,8);
			trim(temp);
			double x = atof(temp.c_str());
			temp = line.substr(39,8);
			trim(temp);
			double y = atof(temp.c_str());
			temp = line.substr(47,8);
			trim(temp);
			double z = atof(temp.c_str());
			if (atomName == "H1" || atomName == "H2")
			{   // n-term H, take average of H1,H2,H3
				ntermx += x;
				ntermy += y;
				ntermz += z;
				numNtermHs++;
				continue;
			}
			else if (atomName == "H3")
			{
				x += ntermx;
				y += ntermy;
				z += ntermz;
				numNtermHs++;
				x = x/double(numNtermHs);
				y = y/double(numNtermHs);
				z = z/double(numNtermHs);
				atomName = "H";
				ntermx = 0;
				ntermy = 0;
				ntermz = 0;
				numNtermHs=0;
			}
			res->addCoords(atomName,x,y,z);
		}
		if (models.size() != numModels)
			printf("WARNING: Actual num models %zu not equal to expected number %d. Missing NUMMDL record in pdb file?\n",
					models.size(),numModels);
	}
	else
	{
		cout << "Unable to open file " << pdbFile << endl;
		file.close();
		exit(-1);
	}
}

/*
 * @param shifts
 * 	indexed by residue number, so first entry (index 0) is empty
 * version 2.1
 * serial number = _Atom_shift_assign_ID
 * res num =  _Residue_seq_code
 * res name = _Residue_label
 * atom name = _Atom_name
 * shift value = _Chem_shift_value
 * version 3.1
 * serial number =  _Atom_chem_shift.ID
 * res num = _Atom_chem_shift.Seq_ID
 * res name =_Atom_chem_shift.Comp_ID
 * atom name = _Atom_chem_shift.Atom_ID
 * shift value =  _Atom_chem_shift.Val
 *
 * use default 2.1 values if not header found
 * e.g.  1    1    MET HA    H   3.850 0.02 1
 */
void CSProtein::getShifts(const string& bmrbFile, int offset)
{
	ifstream file;
	string line;
	int len;
	int colIndex=1, resNumIndex=2, resLabelIndex=3, atomNameIndex=4, shiftIndex=6; // 1-based
	bool hasHeader = false; // true if have header for chemical shift section
	bool hasSeparator = true; // true if there is a blank line separator between header and body
	int resnum=0;
	string resname;
	string atomname;
	double shift;
	file.open(bmrbFile.c_str());
	if (file.is_open())
	{
		// skip to header
		while (file.good())
		{
			getline(file,line);
			trim(line);
			if (line.empty())
				continue;
			len = line.length();
			if (len > 18)
			{
				if (starts_with(line,"#"))
				{
					continue;
				}
				else if (starts_with(line,"_Atom_shift_assign_ID"))
				{
					hasHeader = true;
					break;
				}
				else if (starts_with(line,"_Atom_chem_shift.ID"))
				{
					hasHeader = true;
					break;
				}
			}
		}
		if (!hasHeader) // use default col index values declared above;
		{
			file.close();
			file.open(bmrbFile.c_str());
			// peak at line to see if format is 1  1  1 MET HA  H 3.850 0.02 1 (2 resnums)
		    int pos = file.tellg(); // get current position
		    getline(file,line);
		    trim(line);
		    stringstream tok(line);
		    string temp;
		    tok >> temp; // serial
		    tok >> temp; // resnum1
		    tok >> temp; // resnum2
		    if (isdigit(temp[0]))
		    {
		    	resLabelIndex=4;
		    	atomNameIndex=5;
		    	shiftIndex=7;
		    }
		    file.seekg(pos ,std::ios_base::beg);// return to pos
		}
		else
		{
			// get the column indices
			while (file.good())
			{
				colIndex++;
				getline(file,line);
				trim(line);
				if (line.empty())
					break; // assume empty line separates header from shifts
				           // in case not the case, look for serial number 1; put back line
				if (starts_with(line,"1"))
				{
					hasSeparator = false;
					break;
				}
				if (starts_with(line,"#"))
				{
					continue;
				}
				else if (line == "_Residue_seq_code" ||
						line == "_Atom_chem_shift.Seq_ID")
				{
					resNumIndex = colIndex;
				}
				else if (line == "_Residue_label" ||
						line == "_Atom_chem_shift.Comp_ID")
				{
					resLabelIndex = colIndex;
				}
				else if (line == "_Atom_name" ||
						line == "_Atom_chem_shift.Atom_ID")
				{
					atomNameIndex = colIndex;
				}
				else if (line == "_Chem_shift_value" ||
						line == "_Atom_chem_shift.Val")
				{
					shiftIndex = colIndex;
				}
			}
		}
		while (file.good())
		{
			if (hasSeparator)
			{
				getline(file,line);
				trim(line);
			}
			else
			{
				hasSeparator = true; // so will fetch line in next iter
			}
			if (line.empty())
				continue;
			if (starts_with(line,"#"))
				continue;
			if (starts_with(line,"stop"))
				break;

			stringstream tok(line);
			string temp;
			colIndex = 1;
			shift = INVALIDSHIFTVAL;
			while (tok >> temp)
			{
				if (colIndex == resNumIndex)
				{
					resnum = atoi(temp.c_str());
				}
				else if (colIndex == resLabelIndex)
				{
					resname = temp;
				}
				else if (colIndex == atomNameIndex)
				{
					atomname = temp;
					if (atomname == "HN")
						atomname = "H"; // H is the standard name rather than HN
				}
				else if (colIndex == shiftIndex)
				{
					shift = atof(temp.c_str());
					if (resnum > 0)
					{
						int index = resnum+offset-startResNum;
						if (index >= 0)
						{
							Residue* res = residues[index];
							if (res == NULL)
							{
								printf("getShifts: No such residue in sequence %d\n",(resnum+offset));
								continue;
							}
							res->addShift(atomname,shift);
							// print warning if AA types in BMRB don't match sequence
							// note: CYSS will match CYS (only first 3 chars compared)
							// if use int type2 = Residue::index3(resname.c_str());
							if (Residue::AA3S[residues[index]->type] != resname)
							{
								printf("WARNING: BMRB type %s does not match seq type %s at res %d\n",
										resname.c_str(),Residue::AA3S[res->type],resnum);
							}
						}
					}
					resnum = -1;
					resname.clear();
					atomname.clear();
				}
				colIndex++;
			}
		}
	}
	else
	{
		cout << "Unable to open file " << bmrbFile << endl;
		file.close();
		exit(-1);
	}
	file.close();
}

bool CSProtein::inContact(int resNum1, const string& parent1, const string& proton1,
		int resNum2, const string& parent2, const string& proton2, double distCutoff) const
{
	double dist;
	return inContact(resNum1,parent1,proton1,resNum2,parent2,proton2,distCutoff,dist);
}

bool CSProtein::inContact(int resNum1, const string& parent1, const string& proton1,
		int resNum2, const string& parent2, const string& proton2,
		double distCutoff, double& retDist) const
{
	Residue* res1 = residues[resNum1-startResNum];
	Residue* res2 = residues[resNum2-startResNum];
	Atom* c1=NULL;
	Atom* c2=NULL;
	retDist = INVALIDDIST; // we will set to the min dist
	if (parent1[0] == 'C')
		c1 = res1->getC(parent1);
	else if (parent1[0] == 'N')
		c1 = res1->getN(parent1);
	else // treat as C
		c1 = res1->getC(parent1);

	if (parent2[0] == 'C')
		c2 = res2->getC(parent2);
	else if (parent2[0] == 'N')
		c2 = res2->getN(parent2);
	else // treat as C
		c2 = res2->getC(parent2);

	// check if methyl
	if (c1->isMethyl())
	{
		if (c2->isMethyl())
		{
			for (HIterator it1 = c1->beginH(); it1 != c1->endH(); ++it1)
			{
				Atom* h1 = *it1;
				for (HIterator it2 = c2->beginH(); it2 != c2->endH(); ++it2)
				{
					Atom* h2 = *it2;
					if (h1->x != INVALIDCOORD && h2->x != INVALIDCOORD)
					{
						double dist = h1->getDistance(h2);
						if (dist < retDist)
							retDist = dist;
					}
				}
			}
		}
		else
		{
			Atom* h2 = c2->getH(proton2);
			for (HIterator it1 = c1->beginH(); it1 != c1->endH(); ++it1)
			{
				Atom* h1 = *it1;
				if (h1->x != INVALIDCOORD && h2->x != INVALIDCOORD)
				{
					double dist = h1->getDistance(h2);
					if (dist < retDist)
						retDist = dist;
				}
			}
		}
	}
	else
	{
		Atom* h1 = c1->getH(proton1);
		if (c2->isMethyl())
		{
			retDist = INVALIDDIST;
			for (HIterator it2 = c2->beginH(); it2 != c2->endH(); ++it2)
			{
				Atom* h2 = *it2;
				if (h1->x != INVALIDCOORD && h2->x != INVALIDCOORD)
				{
					double dist = h1->getDistance(h2);
					if (dist < retDist)
						retDist = dist;
				}
			}
		}
		else
		{
			Atom* h2 = c2->getH(proton2);
			if (h1->x != INVALIDCOORD && h2->x != INVALIDCOORD)
			{
				retDist = h1->getDistance(h2);
			}
		}
	}
	if (retDist <= distCutoff)
		return true;
	else
	{
		if (retDist == INVALIDDIST)
			retDist = -1; // invalid coordinates
		return false;
	}
}

bool CSProtein::inContact(PseudoContact& pc, double distCutoff) const
{
	tr1::unordered_set<Contact>& contacts = pc.getContacts();
	for (tr1::unordered_set<Contact>::const_iterator it = contacts.begin(); it != contacts.end(); ++it)
	{
		const Contact& c = *it;
		if (inContact(c.r1->num, c.x1->name, c.hx1->name, c.r2->num, c.x2->name, c.hx2->name,distCutoff))
			return true;
	}
	return false;
}

bool CSProtein::inContact(PseudoContact& pc, double distCutoff, double& retDist) const
{
	tr1::unordered_set<Contact>& contacts = pc.getContacts();
	bool inContactFlag = false;
	retDist = INVALIDDIST; // we will set to the min dist
	for (tr1::unordered_set<Contact>::const_iterator it = contacts.begin(); it != contacts.end(); ++it)
	{
		const Contact& c = *it;
		double dist = 0;
		if (inContact(c.r1->num, c.x1->name, c.hx1->name, c.r2->num, c.x2->name, c.hx2->name,distCutoff,dist))
			inContactFlag = true;
		if (dist >= 0 && dist < retDist)
			retDist = dist;
	}
	if (retDist == INVALIDDIST)
	{
		retDist = -1; // invalid coordinates
		return false;
	}
	return inContactFlag;
}

bool CSProtein::inContact(list<CSProtein*>& models, PseudoContact& pc, double distCutoff)
{
	for (list<CSProtein*>::iterator it = models.begin(); it != models.end(); ++it)
	{
		CSProtein* p = *it;
		if (p->inContact(pc,distCutoff))
			return true;
	}
	return false;
}

bool CSProtein::inContact(list<CSProtein*>& models, PseudoContact& pc, double distCutoff, double& retDist)
{
	double dist;
	retDist = INVALIDDIST;
	for (list<CSProtein*>::iterator it = models.begin(); it != models.end(); ++it)
	{
		CSProtein* p = *it;
		p->inContact(pc,distCutoff,dist);
		if (dist >= 0 && dist < retDist)
			retDist = dist;
	}
	if (retDist < distCutoff)
		return true;
	if (retDist == INVALIDDIST)
		retDist = -1;
	return false;
}

bool CSProtein::inContact(list<CSProtein*>& models, PseudoContact& pc, double distCutoff,
		double& retDist, int& numModels)
{
	double dist;
	retDist = INVALIDDIST;
	numModels = 0;
	for (list<CSProtein*>::iterator it = models.begin(); it != models.end(); ++it)
	{
		CSProtein* p = *it;
		if (p->inContact(pc,distCutoff,dist))
			numModels++;
		if (dist >= 0 && dist < retDist)
			retDist = dist;
	}
	if (retDist < distCutoff)
		return true;
	if (retDist == INVALIDDIST)
		retDist = -1;
	return false;
}


bool CSProtein::inContact(list<CSProtein*>& models, Contact& c, double distCutoff)
{
	for (list<CSProtein*>::iterator it = models.begin(); it != models.end(); ++it)
	{
		CSProtein* p = *it;
		if (p->inContact(c,distCutoff))
			return true;
	}
	return false;
}

bool CSProtein::inContact(list<CSProtein*>& models, Contact& c, double distCutoff,
		double& retDist)
{
	double dist;
	retDist =  INVALIDDIST;
	for (list<CSProtein*>::iterator it = models.begin(); it != models.end(); ++it)
	{
		CSProtein* p = *it;
		p->inContact(c,distCutoff,dist);
		if (dist >= 0 && dist < retDist)
			retDist = dist;
	}
	if (retDist < distCutoff)
		return true;
	if (retDist == INVALIDDIST)
		retDist = -1;
	return false;
}

bool CSProtein::inContact(list<CSProtein*>& models, const Contact& c, double distCutoff,
		double& retDist)
{
	double dist;
	retDist =  INVALIDDIST;
	for (list<CSProtein*>::iterator it = models.begin(); it != models.end(); ++it)
	{
		CSProtein* p = *it;
		p->inContact(c,distCutoff,dist);
		if (dist >= 0 && dist < retDist)
			retDist = dist;
	}
	if (retDist < distCutoff)
		return true;
	if (retDist == INVALIDDIST)
		retDist = -1;
	return false;
}

bool CSProtein::inContact(list<CSProtein*>& models, Contact& c, double distCutoff,
		double& retDist, int& numModels)
{
	double dist;
	retDist =  INVALIDDIST;
	numModels = 0;
	for (list<CSProtein*>::iterator it = models.begin(); it != models.end(); ++it)
	{
		CSProtein* p = *it;
		if (p->inContact(c,distCutoff,dist))
			numModels++;
		if (dist >= 0 && dist < retDist)
			retDist = dist;
	}
	if (retDist < distCutoff)
		return true;
	if (retDist == INVALIDDIST)
		retDist = -1;
	return false;
}

bool CSProtein::inContact(Contact& c, double distCutoff) const
{
	if (inContact(c.r1->num, c.x1->name, c.hx1->name, c.r2->num, c.x2->name, c.hx2->name,distCutoff))
		return true;
	else
		return false;
}

bool CSProtein::inContact(Contact& c, double distCutoff, double& retDist) const
{
	if (inContact(c.r1->num, c.x1->name, c.hx1->name, c.r2->num, c.x2->name, c.hx2->name,distCutoff, retDist))
		return true;
	else
		return false;
}

bool CSProtein::inContact(const Contact& c, double distCutoff, double& retDist) const
{
	if (inContact(c.r1->num, c.x1->name, c.hx1->name, c.r2->num, c.x2->name, c.hx2->name,distCutoff, retDist))
		return true;
	else
		return false;
}

double CSProtein::getMinDistance(PseudoContact& pc)
{
	double minDist = INVALIDDIST;
	tr1::unordered_set<Contact>& contacts = pc.getContacts();
   // for each contact get min dist
   for (tr1::unordered_set<Contact>::const_iterator itC = contacts.begin(); itC != contacts.end(); ++itC)
   {
	   const Contact& contact = *itC;
	   int resNum1 = contact.r1->num;
	   int resNum2 = contact.r2->num;
	   string& parent1 = contact.x1->name;
	   string& parent2 = contact.x2->name;
	   string& proton1 = contact.hx1->name;
	   string& proton2 = contact.hx2->name;
	   Residue* res1 = residues[resNum1-startResNum];
	   Residue* res2 = residues[resNum2-startResNum];
	   Atom* c1=NULL;
	   Atom* c2=NULL;
	   if (parent1[0] == 'C')
		   c1 = res1->getC(parent1);
	   else if (parent1[0] == 'N')
		   c1 = res1->getN(parent1);
	   else // treat as C
		   c1 = res1->getC(parent1);
	   if (parent2[0] == 'C')
		   c2 = res2->getC(parent2);
	   else if (parent2[0] == 'N')
		   c2 = res2->getN(parent2);
	   else // treat as C
		   c2 = res2->getC(parent2);
	   // check if methyl
	   if (c1->isMethyl())
	   {
		   if (c2->isMethyl())
		   {
			   for (HIterator it1 = c1->beginH(); it1 != c1->endH(); ++it1)
			   {
				   Atom* h1 = *it1;
				   for (HIterator it2 = c2->beginH(); it2 != c2->endH(); ++it2)
				   {
					   Atom* h2 = *it2;
					   if (h1->x != INVALIDCOORD && h2->x != INVALIDCOORD)
					   {
						   double dist = h1->getDistance(h2);
						   if (dist < minDist)
							   minDist = dist;
					   }
				   }
			   }
		   }
		   else
		   {
			   Atom* h2 = c2->getH(proton2);
			   for (HIterator it1 = c1->beginH(); it1 != c1->endH(); ++it1)
			   {
				   Atom* h1 = *it1;
				   if (h1->x != INVALIDCOORD && h2->x != INVALIDCOORD)
				   {
					   double dist = h1->getDistance(h2);
					   if (dist < minDist)
						   minDist = dist;
				   }
			   }
		   }
	   }
	   else
	   {
		   Atom* h1 = c1->getH(proton1);
		   if (c2->isMethyl())
		   {
			   for (HIterator it2 = c2->beginH(); it2 != c2->endH(); ++it2)
			   {
				   Atom* h2 = *it2;
				   if (h1->x != INVALIDCOORD && h2->x != INVALIDCOORD)
				   {
					   double dist = h1->getDistance(h2);
					   if (dist < minDist)
						   minDist = dist;
				   }
			   }
		   }
		   else
		   {
			   Atom* h2 = c2->getH(proton2);
			   if (h1->x != INVALIDCOORD && h2->x != INVALIDCOORD)
			   {
				   double dist = h1->getDistance(h2);
				   if (dist < minDist)
					   minDist = dist;
			   }
		   }
	   }
   }
   if (minDist == INVALIDDIST)
	   return -1;
   else
	   return minDist;
}

double CSProtein::getDistance(int resNum1, const string& parent1, const string& proton1,
		                      int resNum2, const string& parent2, const string& proton2) const
{
	Residue* res1 = residues[resNum1-startResNum];
	Residue* res2 = residues[resNum2-startResNum];
	Atom* c1=NULL;
	Atom* c2=NULL;
	if (parent1[0] == 'C')
		c1 = res1->getC(parent1);
	else if (parent1[0] == 'N')
		c1 = res1->getN(parent1);
	else // treat as C
		c1 = res1->getC(parent1);

	if (parent2[0] == 'C')
		c2 = res2->getC(parent2);
	else if (parent2[0] == 'N')
		c2 = res2->getN(parent2);
	else // treat as C
		c2 = res2->getC(parent2);

	Atom* h1 = c1->getH(proton1);
	Atom* h2 = c2->getH(proton2);
	if (h1->x != INVALIDCOORD && h2->x != INVALIDCOORD)
		return h1->getDistance(h2);
	else
		return -1;
}

double CSProtein::getSideChainCentroidDistance(int res1, int res2) const
{
	Residue* r1 = residues[res1-startResNum];
	Residue* r2 = residues[res2-startResNum];
	double x1,y1,z1,x2,y2,z2;
	r1->getSideChainCentroid(x1,y1,z1);
	r2->getSideChainCentroid(x2,y2,z2);
	if (x1 != INVALIDCOORD && x2 != INVALIDCOORD)
		return distance(x1,y1,z1,x2,y2,z2);
	else
		return -1;
}

double CSProtein::getSideChainCentroidCADistance(int res1, int res2) const
{
	Residue* r1 = residues[res1-startResNum];
	Residue* r2 = residues[res2-startResNum];
	double x1,y1,z1,x2,y2,z2;
	r1->getSideChainCentroid(x1,y1,z1);
	Atom* a2 = r2->getC("CA");
	if (a2 != NULL)
	{
		x2 = a2->x;
		y2 = a2->y;
		z2 = a2->z;
	}
	else
	{
		return -1;
	}
	if (x1 != INVALIDCOORD && x2 != INVALIDCOORD)
		return distance(x1,y1,z1,x2,y2,z2);
	else
		return -1;
}

double CSProtein::getSideChainCentroidCBDistance(int res1, int res2) const
{
	Residue* r1 = residues[res1-startResNum];
	Residue* r2 = residues[res2-startResNum];
	double x1,y1,z1,x2,y2,z2;
	r1->getSideChainCentroid(x1,y1,z1);
	Atom* a2 = r2->getC("CB");
	if (a2 != NULL)
	{
		x2 = a2->x;
		y2 = a2->y;
		z2 = a2->z;
	}
	else
	{
		return -1;
	}
	if (x1 != INVALIDCOORD && x2 != INVALIDCOORD)
		return distance(x1,y1,z1,x2,y2,z2);
	else
		return -1;
}

double CSProtein::getSideChainCentroidNDistance(int res1, int res2) const
{
	Residue* r1 = residues[res1-startResNum];
	Residue* r2 = residues[res2-startResNum];
	double x1,y1,z1,x2,y2,z2;
	r1->getSideChainCentroid(x1,y1,z1);
	Atom* a2 = r2->getN("N");
	if (a2 != NULL)
	{
		x2 = a2->x;
		y2 = a2->y;
		z2 = a2->z;
	}
	else
	{
		return -1;
	}
	if (x1 != INVALIDCOORD && x2 != INVALIDCOORD)
		return distance(x1,y1,z1,x2,y2,z2);
	else
		return -1;
}

double CSProtein::getDistance(int res1, int res2, Contact::ResidueContactType type) const
{
	Residue* r1 = residues[res1-startResNum];
	Residue* r2 = residues[res2-startResNum];
	Atom* a1 = NULL;
	Atom* a2 = NULL;

	if (type == Contact::SCSC)
		return getSideChainCentroidDistance(res1,res2);
	else if (type == Contact::SCCB)
		return getSideChainCentroidCBDistance(res1,res2);
	else if (type == Contact::CBSC)
		return getSideChainCentroidCBDistance(res2,res1);
	else if (type == Contact::SCCA)
		return getSideChainCentroidCADistance(res1,res2);
	else if (type == Contact::CASC)
		return getSideChainCentroidCADistance(res2,res1);
	else if (type == Contact::SCN)
		return getSideChainCentroidNDistance(res1,res2);
	else if (type == Contact::NSC)
		return getSideChainCentroidNDistance(res2,res1);

	if (type == Contact::CBCB || type == Contact::CBCA || type == Contact::CBN)
	{
		a1 = r1->getC("CB");
	}
	if (type == Contact::CACB || type == Contact::CACA || type == Contact::CAN)
	{
		a1 = r1->getC("CA");
	}
	if (type == Contact::NCB || type == Contact::NCA || type == Contact::NN)
	{
		a1 = r1->getC("N");
	}
	if (type == Contact::CBCB || type == Contact::CACB || type == Contact::NCB)
	{
		a2 = r2->getC("CB");
	}
	if (type == Contact::CBCA || type == Contact::CACA || type == Contact::NCA)
	{
		a2 = r2->getC("CA");
	}
	if (type == Contact::CBN || type == Contact::CAN || type == Contact::NN)
	{
		a2 = r2->getC("N");
	}
	if (a1 != NULL && a2 != NULL)
	{
		return a1->getDistance(a2);
	}
	else
	{
		return -1;
	}
}

double CSProtein::getSideChainCentroidDistance(list<CSProtein*>& models, int res1, int res2)
{
	double minDist = INVALIDDIST;
	for (list<CSProtein*>::iterator itS = models.begin(); itS != models.end(); ++itS)
	{
		CSProtein* cp = *itS;
		double dist = cp->getSideChainCentroidDistance(res1,res2);
		if (dist >= 0 && dist < minDist)
			minDist = dist;
	}
	if (minDist != INVALIDDIST)
		return minDist;
	else
		return -1;
}

void CSProtein::write15NChemShift(const string& filename, int offset, bool ha, bool hb) const
{
	FILE* out = fopen(filename.c_str(), "w");
	int serial = 1;
	for (int i = 0; i < size; i++)
	{
		Residue* res = residues[i];
		int resnum = res->num+(startResNum-1)+offset;
		for (AtomIterator it = res->begin('X'); it != res->end('X'); it++) // iterate over the heavy atoms
		{
			Atom* atom = *it;
			if (atom->name[0] != 'N')
			{
				if (ha || hb)
				{
					if (ha && hb)
					{
						if (atom->name.substr(0,2) != "CA" && atom->name.substr(0,2) != "CB")
							continue;
					}
					else if (ha)
					{
						if (atom->name.substr(0,2) != "CA")
							continue;
					}
					else // hb only
					{
						if (atom->name.substr(0,2) != "CB")
							continue;
					}
				}
				else
					continue;
			}
			else
			{
				// exclude sidechain HIS, LYS, ARG N's
				if ( (res->type == Residue::HIS || res->type == Residue::LYS || res->type == Residue::ARG) &&
						atom->name.size() > 1 )
					continue;
			}
			// print the heavy atom
			if (atom->cs != INVALIDSHIFTVAL)
			{
				fprintf(out,"%5d %3d %3d %3s %4s %c %6.2f\n",serial,resnum,resnum,res->getAAType3(),
						atom->name.c_str(),atom->getType(),atom->cs);
				serial++;
			}
			// print the attached protons
			if (!atom->isMethyl())
			{
				for (HIterator itH = atom->beginH(); itH != atom->endH(); itH++)
				{
					Atom* hAtom = *itH;
					if (hAtom->cs != INVALIDSHIFTVAL)
					{
						fprintf(out,"%5d %3d %3d %3s %4s %c %6.2f\n",serial,resnum,resnum,res->getAAType3(),
								hAtom->name.c_str(),'H',hAtom->cs);
						serial++;
					}
				}
			}
			else
			{
				// for methyls, truncate the last digit, and print only one of the protons
				for (HIterator itH = atom->beginH(); itH != atom->endH(); itH++)
				{
					Atom* hAtom = *itH;
					if (hAtom->cs != INVALIDSHIFTVAL)
					{
						string name = hAtom->name.substr(0,hAtom->name.length()-1);
						fprintf(out,"%5d %3d %3d %3s %4s %c %6.2f\n",serial,resnum,resnum,res->getAAType3(),
								name.c_str(),'H',hAtom->cs);
						serial++;
					}
					break;
				}
			}
		}
	}
	fclose(out);
}

//int main(int argc, char** argv)
//{
//	string seq;
//	CSProtein::append("/home/e4k2/Documents/nmr_data/casd-nmr/2kky/seq.fasta",seq);
//	CSProtein p(seq,1);
//	p.getShifts("/home/e4k2/Documents/nmr_data/casd-nmr/2kky/oxidized_ET109A/oxidized_ET109A.bmrb",-89);
//	p.getCoords("/home/e4k2/Documents/nmr_data/casd-nmr/2kky/pdb.pdb",-89);
//	p.print();


//	Residue r(1,Residue::TRP);
//	r.print();
//	Residue r2(1,Residue::TYR);
//	r2.print();
//	Residue r3(1,Residue::VAL);
//	r3.print();
//	list<Atom> list1;
//	Atom a1("CA");
//	a1.x = 1;
//	Atom* a3 = new Atom("HA");
//	a3->x = 3;
//	Atom a2("CB");
//	a2.x = 2;
//	list1.push_back(a1);
//	list1.push_back(a2);
//	list<Atom> list2(list1);
//	list<Atom>::iterator it1 = list1.begin();
//	list<Atom>::iterator it2 = list2.begin();
//	while (it1 != list1.end())
//	{
//		Atom& a1 = *it1;
//		Atom& a2 = *it2;
//		cout << (*it1).x << " " << &a1 << " " << (*it2).x << " " << &a2 << endl;
//		it1++;
//		it2++;
//	}
//	delete a3;
//}

