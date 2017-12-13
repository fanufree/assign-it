/*
 * Utilities.cpp
 *
 *  Created on: 2012-04-27
 *      Author: e4k2
 */

#include <cstdio>
#include <fstream>
#include <string>
#include <string.h>
#include <algorithm>
#include <cmath>
#include "Utilities.h"


bool starts_with(const std::string& text, const std::string& toMatch)
{
	if (strncmp(text.c_str(),toMatch.c_str(),toMatch.length()) != 0)
		return false;
	else
		return true;
}

bool ends_with(const std::string& text, const std::string& toMatch)
{
	size_t lenT = text.length();
	size_t lenM = toMatch.length();
	if (lenT >= lenM)
	{
		if (text.substr(lenT-lenM,lenM) == toMatch)
			return true;
	}
	return false;
}

// trim from start
std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}

/**
 * Complementary error function for normal distribution
 * For a normally distributed var with std stdev and expected value 0
 * erfcc(a/(stdev*sqrt(2))) = prob var > a union var < -a
 * Ported from python stat library
 */
double erfcc(double x)
{
	double t,z,ans;
	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
			t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
					t*(-0.82215223+t*0.17087277)))))))));
	return  x >= 0.0 ? ans : 2.0-ans;
}

double distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double dx = x1-x2;
	double dy = y1-y2;
	double dz = z1-z2;
	return sqrt(dx*dx+dy*dy+dz*dz);
}

double distance(double a[], double b[], int size)
{
	double d = 0;
	for (int i = 0; i < size; i++)
	{
		d += (a[i]-b[i])*(a[i]-b[i]);
	}
	return sqrt(d);
}

double dot(double a[], double b[], int size)
{
	double result = 0;
	for (int i = 0; i < size; i++)
	{
		result += a[i]*b[i];
	}
	return result;
}

double norm(double a[], int size)
{
	double result = 0;
	for (int i = 0; i < size; i++)
	{
		result += a[i]*a[i];
	}
	return sqrt(result);
}

double norm2(double a[], int size)
{
	double result = 0;
	for (int i = 0; i < size; i++)
	{
		result += a[i]*a[i];
	}
	return result;
}

void avg_stdev(double a[], int size, double& avg, double& stdev)
{
	avg = 0;
	for (int i = 0; i < size; i++)
	{
		avg += a[i];
	}
	avg = avg/(double)size;

	stdev = 0;
	for (int i = 0; i < size; i++)
	{
		stdev += (a[i]-avg)*(a[i]-avg);
	}
	stdev = sqrt(stdev/(double)size);
}

void readFasta(const string& fastaFile, string& seq)
{
	ifstream file;
	string line;
	file.open(fastaFile.c_str());
	if (file.is_open())
	{
		while (file.good())
		{
			getline(file,line);
			trim(line);
			if (line.empty())
				continue;
			if (starts_with(line,">"))
				continue;
			seq.append(line);
		}
	}
	else
	{
		printf("Unable to open file %s\n", fastaFile.c_str());
		file.close();
		exit(-1);
	}
	file.close();
}

