#pragma once

#include "stdafx.h"
typedef std::complex<double> cmplex;

double LaplaceInversion(cmplex (*F)(const cmplex &s),const double &t, const double tolerance);


double testf(const double &t);
cmplex testF(const cmplex &s);
void TestLaplaceInversion();
