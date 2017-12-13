/*
 * Utilities.hpp
 *
 *  Created on: 2012-04-27
 *      Author: e4k2
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <tr1/unordered_map>
#include <tr1/utility>
#include <string>

using namespace std;

bool starts_with(const std::string& text, const std::string& toMatch); // toMatch is usually the shorter string
bool ends_with(const std::string& text, const std::string& toMatch);
string &ltrim(std::string &s);
string &rtrim(std::string &s);
string &trim(std::string &s);
double erfcc(double x);
double distance(double x1, double y1, double z1, double x2, double y2, double z2);
double distance(double a[], double b[], int size);
double dot(double a[], double b[], int size);
double norm(double a[], int size);
double norm2(double a[], int size);
void avg_stdev(double a[], int size, double& avg, double& stdev); // Note: population variance = biased sample variance
                                                       // size can be less than the capacity of the array, but not more

void readFasta(const string& fastaFile, string& seq); // fetches the sequence in fastaFile to seq
	                              // if seq is not empty, the sequence gets appended

using namespace std;

namespace std
{
	namespace tr1
	{
		template<>
		class hash< pair<int,int> >
		{
			public:
			size_t operator()(const pair<int,int> & p) const
			{
				return hash<int>()(p.first) ^ hash<int>()(p.second);
			}
		};
	}
}


#endif /* UTILITIES_H_ */
