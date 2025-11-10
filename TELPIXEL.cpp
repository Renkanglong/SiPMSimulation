#include "TELPIXEL.h"

double TELPIXEL::distance1 = 320.4;
double TELPIXEL::distance2 = 104.8;
double TELPIXEL::distance3 = 25.8;



TELPIXEL::TELPIXEL()
{

}

TELPIXEL::TELPIXEL(const double& d1, const double& d2, const double& d3)
{
	distance1 = d1;
	distance2 = d2;
	distance3 = d3;
}

TELPIXEL::~TELPIXEL()
{
}

void TELPIXEL::partIndex(const postype& pos)
{
	const double x_raw = pos.first;
	const double y_raw = pos.second;
	double ratio1 = 1.0 / distance1;
	double ratio2 = 1.0 / distance2;
	double ratio3 = 1.0 / distance3;
	number1 = { static_cast<int>(std::round(x_raw * ratio1)) ,static_cast<int>(std::round(y_raw * ratio1)) };
	if (std::abs(number1.first) > 2 || std::abs(number1.second) > 2)
	{
		return;
	}
	if (std::abs(number1.first) == 2 && std::abs(number1.second) == 2)
	{
		return;
	}
	double x_new = x_raw - number1.first * distance1;
	double y_new = y_raw - number1.second * distance1;
	number2 = { static_cast<int>(std::round(x_new * ratio2)) ,static_cast<int>(std::round(y_new * ratio2)) };
	if (std::abs(number2.first) > 1 || std::abs(number2.second) > 1)
	{
		return;
	}
	if (!Boundary())
		return;
	x_new = x_new - number2.first * distance2 ;
	y_new = y_new - number2.second * distance2 ;
	number3 = { floor(x_new * ratio3) ,floor(y_new * ratio3) };
	if (number3.first > 1 || number3.first < -2 || number3.second > 1 || number3.second < -2)
	{
		return;
	}
	pixel_index = { number1.first * distance1 + number2.first * distance2 + number3.first * distance3,
	number1.second * distance1 + number2.second * distance2 + number3.second * distance3 };
}
void TELPIXEL::getIndex()
{
	std::cout << "Ïà»úÏñËØ±àºÅ1:£¨" << number1.first << "," << number1.second << ")" << "\n";
	std::cout << "Ïà»úÏñËØ±àºÅ2:£¨" << number2.first << "," << number2.second << ")" << "\n";
	std::cout << "Ïà»úÏñËØ±àºÅ3:£¨" << number3.first << "," << number3.second << ")" << "\n";
	std::cout << "Ïà»úÏñËØË÷Òý:£¨" << pixel_index.first << "," << pixel_index.second << ")" << "\n";
}

bool TELPIXEL::Boundary()
{
	if (!(std::abs(number1.first) == 2 || std::abs(number1.second) == 2))
		return true;
	if (number1.first == -2 && number1.second == 0)
	{
		if (number2.first != 1)
		{
			return false;
		}
	}
	if (number1.first == 2 && number1.second == 0)
	{
		if (number2.first != -1)
		{
			return false;
		}
	}
	if (number1.first == 0 && number1.second == 2)
	{
		if (number2.second != -1)
		{
			return false;
		}
	}
	if (number1.first == 0 && number1.second == -2)
	{
		if (number2.second != 1)
		{
			return false;
		}
	}
	if (number1.first == -2 && number1.second == 1)
	{
		if (number2.first != 1 || number2.second != -1)
		{
			return false;
		}
	}
	if (number1.first == -2 && number1.second == -1)
	{
		if (number2.first != 1 || number2.second != 1)
		{
			return false;
		}
	}
	if (number1.first == 2 && number1.second == 1)
	{
		if (number2.first != -1 || number2.second != -1)
		{
			return false;
		}
	}
	if (number1.first == 2 && number1.second == -1)
	{
		if (number2.first != -1 || number2.second != 1)
		{
			return false;
		}
	}
	if (number1.first == -1 && number1.second == 2)
	{
		if (number2.first != 1 || number2.second != -1)
		{
			return false;
		}
	}
	if (number1.first == 1 && number1.second == 2)
	{
		if (number2.first != -1 || number2.second != -1)
		{
			return false;
		}
	}
	if (number1.first == -1 && number1.second == -2)
	{
		if (number2.first != 1 || number2.second != 1)
		{
			return false;
		}
	}
	if (number1.first == 1 && number1.second == -2)
	{
		if (number2.first != -1 || number2.second != 1)
		{
			return false;
		}
	}
}

std::map<postype, std::vector<postype>> TELPIXEL::partOut(const std::vector<postype>& allpos)
{
	std::map<postype, std::vector<postype>> partOut;
	for (const postype& pos : allpos)
	{
		partIndex(pos);
		partOut[pixel_index].emplace_back(pos);
	}
	return partOut;
}



