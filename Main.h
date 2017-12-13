/*
 * Main.h
 *
 * Contains global constants
 *
 *  Created on: Jun 20, 2014
 *      Author: e4k2
 */

#ifndef MAIN_H_
#define MAIN_H_

#include <string>

extern const double INVALIDDISTANCE; // initial distance for finding the minimum distance; initial val should be some big number
                                       // defined in Main.cpp
extern const int LRSEQSEP; //  minimum sequence separation to be considered as a long range contact (inclusive >= LRSEQSEP);
                           // used for filtering assignments; defined in Main.cpp
#endif /* MAIN_H_ */
