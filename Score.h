/*
 * Score.h
 *
 *  Created on: 2013-07-31
 *      Author: e4k2
 */

#ifndef SCORE_H_
#define SCORE_H_

#include <string>

class Score
{
public:
	Score();
	Score(double cs_, double str_, double intensity_, double sym_, double interres_, double net_, double netStr_, double ambig_, double db_);
	virtual ~Score();
	// Score(Score&& s);
	Score(const Score& s);
	friend void swap(Score& s1, Score& s2);
	Score& operator=(Score s);
	void print() const;
	std::string toString() const;
	void setTotal(); // sets total score (using the weights)
	double getUnweightedTotal() const;

	/*
	 * 1=cs, 2=str, 4=int, 8=sym, 16=interes, 32=net, 64=netstr, 128=ambig, bias=256
	 * To turn on all except str, int use 505 (=1+8+16+32+64+128+256)
	 * All except str, int, netstr use 441
	 * Turn on flags to select the given score term.
	 * The score is weighted such that the weights sum to 1
	 */
	double getScore(unsigned int flag) const;
	void setMax(Score& s); // sets the scoring terms to the maximum of this and s, and then setTotal() is called
	void setAdd(Score& s); // sets this to the sum of this and s. Each score term is updatd including total

	double cs;
	double str;
	double intensity;
	double sym;
	double interres;
	double net;
	double netStr;
	double ambig;
	double db;
	double total; // this variable needs to appear last among the energy terms for correct initialization (the total is weighted)

	// energy term weights
	static double CS_WT;
	static double STR_WT;
	static double INTENSITY_WT;
	static double SYM_WT;
	static double INTERRES_WT;
	static double NET_WT;
	static double NETSTR_WT;
	static double AMBIG_WT;
	static double DB_WT;



};

#endif /* SCORE_H_ */
