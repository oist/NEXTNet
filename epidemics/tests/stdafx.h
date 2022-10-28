#pragma once

#include <complex>
#include <algorithm>
#include <filesystem>

#if defined(_GLIBCXX_RELEASE) && (_GLIBCXX_RELEASE >= 9)
#define HAVE_STD_EXECUTION 1
#include <execution>
#else
#define HAVE_STD_EXECUTION 0
#endif

#include <boost/iterator/counting_iterator.hpp>

#include "catch/catch.hpp"
#include "../stdafx.h"
