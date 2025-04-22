// PStreams - POSIX Process I/O for C++

//        Copyright (C) 2001 - 2024 Jonathan Wakely
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
//
// SPDX-License-Identifier: BSL-1.0

/**
 * @file pstream.cc
 * @brief Define static members of non-templates
 * @author Florian Pflug
 */

#include "nextnet/pstream/pstream.h"

namespace redi
{

const pstreams::pmode pstreams::pstdin;
const pstreams::pmode pstreams::pstdout;
const pstreams::pmode pstreams::pstderr;

}
