/*
 * Score.cpp
 *
 *  Created on: 2013-07-31
 *      Author: e4k2
 */

#include <cstdio>
#include <algorithm>
#include "Score.h"

double Score::CS_WT = 0.111;
double Score::STR_WT = 0.111;
double Score::INTENSITY_WT = 0.111;
double Score::SYM_WT = 0.111;
double Score::INTERRES_WT = 0.111;
double Score::NET_WT = 0.111;
double Score::NETSTR_WT = 0.111;
double Score::AMBIG_WT = 0.111;
double Score::DB_WT = 0.111;

Score::Score() : cs(0), str(0), intensity(0), sym(0), interres(0), net(0), netStr(0), ambig(0), db(0), total(0)
{
}

Score::Score(double cs_, double str_, double intensity_, double sym_, double interres_, double net_, double netStr_, double ambig_, double db_) : cs(cs_), str(str_),
		intensity(intensity_), sym(sym_), interres(interres_), net(net_), netStr(netStr_), ambig(ambig_), db(db_),
		total(CS_WT*cs+STR_WT*str+INTENSITY_WT*intensity+SYM_WT*sym+INTERRES_WT*interres+NET_WT*net+NETSTR_WT*netStr+AMBIG_WT*ambig+DB_WT*db)
{
}

Score::~Score()
{
}

Score::Score(const Score& s) : cs(s.cs), str(s.str),
		intensity(s.intensity), sym(s.sym), interres(s.interres), net(s.net), netStr(s.netStr), ambig(s.ambig), db(s.db), total(s.total)
{
}

void swap(Score& s1, Score& s2)
{
	using std::swap;
	swap(s1.cs, s2.cs);
	swap(s1.str, s2.str);
	swap(s1.intensity, s2.intensity);
	swap(s1.sym, s2.sym);
	swap(s1.interres,s2.interres);
	swap(s1.net, s2.net);
	swap(s1.netStr, s2.netStr);
	swap(s1.ambig, s2.ambig);
	swap(s1.db, s2.db);
	swap(s1.total, s2.total);
}

//Score::Score(Score&& s) : cs(0), str(0), intensity(0), sym(0), net(0), netStr(0), ambig(0), db(0), total(0)
//{
//	swap(*this,s);
//}

Score& Score::operator=(Score s)
{
	swap(*this,s);
	return *this;
}

void Score::print() const
{
	printf("cs:%8.6f  str:%8.6f  intensity:%8.6f  sym:%8.6f  interres:%8.6f  net:%8.6f  netStr:%8.6f  ambig:%8.6f  db:%8.6f  total:%8.6f\n",
			cs,str,intensity,sym,interres,net,netStr,ambig,db,total);
}

std::string Score::toString() const
{
	char buffer[160]; // 155
	snprintf(buffer,160,"cs:%8.6f  str:%8.6f  intensity:%8.6f  sym:%8.6f  interres:%8.6f  net:%8.6f  netStr:%8.6f  ambig:%8.6f  db:%8.6f  total:%8.6f",
				cs,str,intensity,sym,interres,net,netStr,ambig,db,total);
	std::string s(buffer);
	return s;
}

void Score::setTotal()
{
	total = CS_WT*cs+STR_WT*str+INTENSITY_WT*intensity+SYM_WT*sym+INTERRES_WT*interres+NET_WT*net+NETSTR_WT*netStr+AMBIG_WT*ambig+DB_WT*db;
}

double Score::getUnweightedTotal() const
{
	return cs+str+intensity+sym+interres+net+netStr+ambig+db;
}

// 1=cs, 2=str, 4=int, 8=sym, 16=interes, 32=net, 64=netstr, 128=ambig, bias=256
double Score::getScore(unsigned int flag) const
{
	double score = 0;
	// commonly used weights
	static double wt505 = CS_WT + SYM_WT + INTERRES_WT + NET_WT + NETSTR_WT + AMBIG_WT + DB_WT;
	static double wt441 = CS_WT + SYM_WT + INTERRES_WT + NET_WT + AMBIG_WT + DB_WT;
	double wt = 0;

	if (flag == 505)
	{
		score = CS_WT*cs + SYM_WT*sym + INTERRES_WT*interres + NET_WT*net + NETSTR_WT*netStr + AMBIG_WT*ambig + DB_WT*db;
		return score/wt505;
	}
	else if (flag == 441)
	{
		score = CS_WT*cs + SYM_WT*sym + INTERRES_WT*interres + NET_WT*net + AMBIG_WT*ambig + DB_WT*db;
		return score/wt441;
	}
	// default case
	if (flag & 1)
	{
		wt += CS_WT;
		score += CS_WT*cs;
	}
	if (flag & 2)
	{
		wt += STR_WT;
		score += STR_WT*str;
	}
	if (flag & 4)
	{
		wt += INTENSITY_WT;
		score += INTENSITY_WT*intensity;
	}
	if (flag & 8)
	{
		wt += SYM_WT;
		score += SYM_WT*sym;
	}
	if (flag & 16)
	{
		wt += INTERRES_WT;
		score += INTERRES_WT*interres;
	}
	if (flag & 32)
	{
		wt += NET_WT;
		score += NET_WT*net;
	}
	if (flag & 64)
	{
		wt += NETSTR_WT;
		score += NETSTR_WT*netStr;
	}
	if (flag & 128)
	{
		wt += AMBIG_WT;
		score += AMBIG_WT*ambig;
	}
	if (flag & 256)
	{
		wt += DB_WT;
		score += DB_WT*db;
	}
	return score/wt;
}

void Score::setMax(Score& s)
{
	using std::max;
	cs = max(cs, s.cs);
	str = max(str, s.str);
	intensity = max(intensity, s.intensity);
	sym = max(sym, s.sym);
	interres = max(interres,s.interres);
	net = max(net, s.net);
	netStr = max(netStr, s.netStr);
	ambig = max(ambig, s.ambig);
	db = max(db, s.db);
	setTotal();
}

void Score::setAdd(Score& s)
{
	cs = cs+s.cs;
	str = str+s.str;
	intensity = intensity+s.intensity;
	sym = sym+s.sym;
	interres = interres + s.interres;
	net = net+s.net;
	netStr = netStr+s.netStr;
	ambig = ambig + s.ambig;
	db = db+s.db;
	setTotal();
}
