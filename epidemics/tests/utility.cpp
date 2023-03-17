#include "tests/stdafx.h"

#include "utility.h"


TEST_CASE("integer_set", "[utility]") {
	typedef integer_set<int> set_t;
	
	set_t s1;
	s1.insert(1);
	REQUIRE(s1.size() == 1);
	REQUIRE(*s1.begin() == 1);
	auto i = s1.begin();
	REQUIRE(++i == s1.end());
}
