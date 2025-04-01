//
//  main.cpp
//  Epidemics-Mean-Field
//
//  Created by Samuel Curé on 31/01/2022.
//

#include "stdafx.h"

using namespace std::string_literals;

struct dispatcher
{
    typedef int program_t(int argc, const char *argv[]);

    dispatcher(const std::string &name, const std::function<program_t> &program)
    {
        table.insert({ name, program });
    }

    static int dispatch(const std::string &name, int argc, const char *argv[])
    {
        const auto i = table.find(name);
        if (i == table.end())
            throw std::runtime_error("no program named "s + name);
        return i->second(argc, argv);
    }

    static std::unordered_map<std::string, std::function<program_t>> table;
};

std::unordered_map<std::string, std::function<dispatcher::program_t>> dispatcher::table;

#define STRINGIFY(v) STRINGIFY_(v)
#define STRINGIFY_(v) #v

#define DECLARE_PROGRAM(name)                        \
    int program_##name(int argc, const char **argv); \
    dispatcher dispatch_##name(STRINGIFY(name), program_##name)

/*
 * The main function decides which of a set of named "programs" to call
 * based on the first argument. All arguments (including the first) are passed to
 * the program.
 *
 * To add a new type of program named NAME, do the following
 *   (1) Create a new file epidemics/programs/NAME.cpp
 *   (2) Add a function int program_NAME(int argc, const char* argv[]) to that file
 *   (3) Add a line DECLARE_PROGRAM(NAME) below
 *   (4) Add epidemics/programs/NAME.CPP to SIMULATOR_SOURCES in CMakeLists.txt
 */

DECLARE_PROGRAM(simulate);

int main(int argc, const char *argv[])
{

    if (argc < 2)
        throw std::runtime_error("no program name specified");

    const std::string name = argv[1];
    return dispatcher::dispatch(name, argc, argv);
}
