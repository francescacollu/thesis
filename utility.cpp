#include "utility.hpp"
#include <iostream>

void corner::check(const bool& condition, const char* class_name, const char* msg)
{
	if (!condition)
	{
		std::cout << "ERROR in " << class_name << " : " << msg << std::endl;
		exit(1);
	}
}

void corner::info(const bool& condition, const char* class_name, const char* msg)
{
	if (!condition)
	{
		std::cout << "INFO in " << class_name << " : " << msg << std::endl;
	}
}

void corner::warning(const bool& condition, const char* class_name, const char* msg)
{
	if (!condition)
	{
		std::cout << "WARNING in " << class_name << " : " << msg << std::endl;
	}
}

bool corner::approx_equal(arma::cx_double val1, arma::cx_double val2, double tol )
{
    if(std::abs(val1-val2)<tol)
    {
        return true;
    }
    return false;
}
