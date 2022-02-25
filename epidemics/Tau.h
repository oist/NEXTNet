#pragma once

#include "stdafx.h"

class Tau {
public:
    Tau();
    Tau(double mean, double variance, double r0);
    double mean;
    double variance;
    double r0;
    double beta_normalised();
};
