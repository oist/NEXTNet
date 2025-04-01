#pragma once

#include <complex>
#include <algorithm>
#include <filesystem>

#if defined(_GLIBCXX_RELEASE) && (_GLIBCXX_RELEASE >= 9)
#    define HAVE_STD_EXECUTION 1
#    if PARALLELIZE
#        include <execution>
#    endif
#else
#    define HAVE_STD_EXECUTION 0
#endif

#include <boost/iterator/counting_iterator.hpp>
#include <boost/math/distributions/kolmogorov_smirnov.hpp>

#include "catch/catch.hpp"
#include "../stdafx.h"
