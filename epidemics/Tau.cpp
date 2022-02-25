#include "stdafx.h"
#include "Tau.h"

Tau::Tau():
    mean(5.0),
    variance(3.0),
    r0(2.0){
}
Tau::Tau(
    double m,
    double v,
    double p):
    mean(m),
    variance(v),
    r0(p){
}
