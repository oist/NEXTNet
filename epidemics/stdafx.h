//
//  stdafx.h
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//
// Added to follow apparent conventions of including all libraries in a stdafx file.
#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include <deque>
#include <utility>
#include <numeric>
#include <chrono>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#include <boost/math/distributions.hpp>
#include <boost/math/tools/roots.hpp>
#pragma clang diagnostic pop

/* Make boost::math available as simply bm */
namespace bm = boost::math;
