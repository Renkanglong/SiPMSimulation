#ifndef TELPIXEL_H
#define TELPIXEL_H
#include <utility>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

using postype = std::pair<double, double>;
using idxtype = std::pair<int, int>;

class TELPIXEL
{
public:
	TELPIXEL();
	TELPIXEL(const double& d1, const double& d2, const double& d3);
	~TELPIXEL();

	void partIndex(const postype& pos);
	void getIndex();
	bool Boundary();
	std::map<postype, std::vector<postype>> partOut(const std::vector<postype>& pos);

private:
	static double distance1;
	static double distance2;
	static double distance3;
	idxtype number1;
	idxtype number2;
	idxtype number3;
	postype pixel_index;
};

#endif