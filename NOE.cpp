/*
 * NOE.cpp
 *
 *  Created on: 2012-05-03
 *      Author: e4k2
 */

#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <cctype>
#include "NOE.h"
#include "Atom.h"
#include "Utilities.h"

const string NOE::PEAKTYPENAMES[] = {"N15_3D", "CALI_3D", "CARO_3D", "CHCH_4D", "NHNH_4D", "CHNH_4D", "H2D"};
const char NOE::HEAVYATOM1TYPE[] = {'N','C','C','C','N','C','X'};
const char NOE::HEAVYATOM2TYPE[] = {'X','X','X','C','N','N','X'};

NOE::NOE() : x(INVALIDSHIFTVAL), hx(INVALIDSHIFTVAL), x2(INVALIDSHIFTVAL), h(INVALIDSHIFTVAL), volume(-1), type(UNKNOWN)
{
}

NOE::NOE(double xx, double hxx, double hh, double vol, PEAKLISTTYPE t) : x(xx), hx(hxx), x2(INVALIDSHIFTVAL), h(hh), volume(vol), type(t)
{
}

NOE::NOE(double xx, double hxx, double xx2, double hh, double vol, PEAKLISTTYPE t) : x(xx), hx(hxx),  x2(xx2), h(hh), volume(vol), type(t)
{
}

NOE::NOE(double hxx, double hh, double vol, PEAKLISTTYPE t) : x(INVALIDSHIFTVAL), hx(hxx), x2(INVALIDSHIFTVAL), h(hh), volume(vol), type(t)
{
}

NOE::~NOE()
{
}

NOE::NOE(const NOE& n) : x(n.x), hx(n.hx), x2(n.x2), h(n.h), volume(n.volume), type(n.type)
{
}

void swap(NOE& first, NOE& second)
{
	using std::swap;
	swap(first.x,second.x);
	swap(first.hx,second.hx);
	swap(first.x2,second.x2);
	swap(first.h,second.h);
	swap(first.volume,second.volume);
	swap(first.type,second.type);
}

//NOE::NOE(NOE&& n) : x(INVALIDSHIFTVAL), hx(INVALIDSHIFTVAL), x2(INVALIDSHIFTVAL), h(INVALIDSHIFTVAL), volume(-1), type(UNKNOWN)
//{
//	swap(*this,n);
//}

NOE& NOE::operator=(NOE n)
{
	swap(*this,n);
	return *this;
}

int NOE::getNumDimensions() const
{
	if (type == CALI_3D || type == N15_3D ||type == CARO_3D)
		return 3;
	else if (type == CHCH_4D || type == NHNH_4D || type == CHNH_4D)
		return 4;
	else
		return 2;
}

void NOE::print()
{
	if (type == CALI_3D || type == N15_3D ||type == CARO_3D)
		printf("%s %5.1f %5.2f %5.2f %10.1f\n",PEAKTYPENAMES[type].c_str(),x,hx,h,volume);
	else if (type == CHCH_4D || type == NHNH_4D || type == CHNH_4D)
		printf("%5.1f %5.2f %5.1f %5.2f %10.1f\n",x,hx,x2,h,volume);
	else
		printf("%5.2f %5.2f %10.1f\n",hx,h,volume);
}

void NOE::print(list<NOE*>& noes)
{
	for (list<NOE*>::iterator it = noes.begin(); it != noes.end(); it++)
	{
		NOE* np = *it;
		np->print();
	}
}

void NOE::reference(double offset_x, double offset_hx)
{
	x += offset_x;
	hx += offset_hx;
}

void NOE::reference2(double offset_x2, double offset_h)
{
	x2 += offset_x2;
	h += offset_h;
}

void NOE::reference(double offset_h)
{
	h += offset_h;
}

void NOE::removeDuplicates3D(list<NOE*>& noes, double cTol, double hcTol, double hTol)
{
	for (list<NOE*>::iterator it1 = noes.begin(); it1 != noes.end(); )
	{
		NOE* np1 = *it1;
		bool np1Deleted = false;
		for (list<NOE*>::iterator it2 = noes.begin(); it2 != noes.end(); )
		{
			NOE* np2 = *it2;
			bool np2Deleted = false;
			if (np1 < np2) // to avoid duplicate comparisions like A with B and B with A
			{
				if (fabs(np1->x-np2->x) <= cTol &&
						fabs(np1->hx-np2->hx) <= hcTol &&
						fabs(np1->h-np2->h) <= hTol &&
						fabs(np1->x2-np2->x2) <= cTol) // works for 2D, 3D, 4D
				{
					if (np1->volume > 0 && np1->volume < np2->volume)
					{
						np1Deleted = true;
						delete np1;
						it1 = noes.erase(it1);
						break;
					}
					else
					{
						delete np2;
						it2 = noes.erase(it2);
						np2Deleted = true;
					}
				}
			}

			if (!np2Deleted)
				++it2;
		}
		if (!np1Deleted)
			++it1;
	}
}

void NOE::readNOE3D(const string& filename, int ncIndex, int hncIndex, int indirectIndex,
		int volumeIndex, list<NOE*>& noes, PEAKLISTTYPE t)
{
	// "Assignment" or "w1" in header = Sparky, # in header = Cyana
	ifstream file;
	string line;
	file.open(filename.c_str());
	if (file.is_open() && file.good())
	{
		do
		{
			double nc=INVALIDSHIFTVAL; // N or C shift
			double hnc=INVALIDSHIFTVAL;
			double h=INVALIDSHIFTVAL;
			double vol=-1; // if volumeIndex < 0 then vol stays -1
			int colIndex = 1; // 1-based
			getline(file,line);
			trim(line);

			if (line.empty())
				continue;

			stringstream tok(line);
			string temp;
			bool isValid = true;
			while (tok >> temp)
			{
				if (colIndex == 1 && temp[0] == '#')
				{
					isValid = false;
					break; // noe comment starts with #
				}

				if (colIndex == ncIndex)
				{
					if (isdigit(temp[0]))
					{
						nc = atof(temp.c_str());
					}
					else
					{
						isValid = false;
						break;
					}
				}
				else if (colIndex == hncIndex)
				{
					if (isdigit(temp[0]) || temp[0] == '-')
					{
						hnc = atof(temp.c_str());
					}
					else
					{
						isValid = false;
						break;
					}
				}
				else if (colIndex == indirectIndex)
				{
					if (isdigit(temp[0]) || temp[0] == '-')
					{
						h = atof(temp.c_str());
					}
					else
					{
						isValid = false;
						break;
					}
				}
				else if (colIndex == volumeIndex)
				{
					if (isdigit(temp[0]) || (temp[0] == '-' && isdigit(temp[1])))
					{
						vol = atof(temp.c_str());
					} // else don't care since might not have vol
					else
						vol = -1;
				}
				colIndex++;
			}
			if (isValid && nc != INVALIDSHIFTVAL && hnc != INVALIDSHIFTVAL && h != INVALIDSHIFTVAL)
			{
//				// check for peak in wrong peak list; for now we will ignore such peaks
//				if (t == CALI_3D && nc > 95.0)
//				{
//					printf("Skipping CALI_3D peak %f>95.0 %f %f %f\n",nc,hnc,h,vol);
//					continue;
//				} else if (t == CARO_3D && nc < 90.0)
//				{
//					printf("Skipping CARO_3D peak %f<90.0 %f %f %f\n",nc,hnc,h,vol);
//					continue;
//				}
				NOE* np = new NOE(nc,hnc,h,vol,t);
				noes.push_back(np);
			}
		} while (file.good());
	}
	else
	{
		printf("Unable to open file %s\n",filename.c_str());
		file.close();
		exit(-1);
	}
	file.close();
}

bool NOE::match(double heavy1, double proton1, double heavy2, double proton2,
		double tolHeavy1, double tolProton1, double tolHeavy2, double tolProton2)
{
	if (heavy1 != INVALIDSHIFTVAL && x != INVALIDSHIFTVAL && fabs(heavy1-x) > tolHeavy1)
		return false;
	if (proton1 != INVALIDSHIFTVAL && hx != INVALIDSHIFTVAL && fabs(proton1-hx) > tolProton1)
		return false;
	if (proton2 != INVALIDSHIFTVAL && h != INVALIDSHIFTVAL && fabs(proton2-h) > tolProton2)
		return false;
	if (heavy2 != INVALIDSHIFTVAL && x2 != INVALIDSHIFTVAL && fabs(heavy2-2) > tolHeavy2)
		return false;

	return true;
}

void NOE::getShiftDiff(double heavy1CS, double proton1CS, double heavy2CS, double proton2CS,
		 double& diffx1, double& diffh1, double& diffx2, double& diffh2)
{
	if (heavy1CS != INVALIDSHIFTVAL && x != INVALIDSHIFTVAL)
		diffx1 = fabs(heavy1CS-x);
	else
		diffx1 = INVALIDSHIFTVAL;
	if (proton1CS != INVALIDSHIFTVAL && hx != INVALIDSHIFTVAL)
		diffh1 = fabs(proton1CS-hx);
	else
		diffh1 = INVALIDSHIFTVAL;
	if (proton2CS != INVALIDSHIFTVAL && h != INVALIDSHIFTVAL)
		diffh2 = fabs(proton2CS-h);
	else
		diffh2 = INVALIDSHIFTVAL;
	if (heavy2CS != INVALIDSHIFTVAL && x2 != INVALIDSHIFTVAL)
		diffx2 = fabs(heavy2CS-x2);
	else
		diffx2 = INVALIDSHIFTVAL;
}


