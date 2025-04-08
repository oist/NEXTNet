//
//  factories.cpp
//  NEXTNet
//
//  Created by Florian G. Pflug on 08.04.25.
//

#include "nextnet/stdafx.h"
#include "nextnet/factories/factories.h"

namespace factories {

parsed_expression_t parse_expression(std::string s)
{
	std::string name;
	std::vector<char> bstack;
	bool qinclude = false;
	std::vector<std::string> args;
\
	enum {
		WS1,
		NAME,
		WS2,
		BRA,
		WS3,
		ARG,
		ARG_DQ,
		ARG_DQE,
		WS4,
		KET_OR_COMMA,
		WS5,
		DONE
	} state = WS1;
	for (std::ptrdiff_t i = 0; i < (ptrdiff_t)s.size(); ++i) {
		const char c = s[i];
		switch (state) {
			case WS1:
			case WS2:
			case WS3:
			case WS4:
			case WS5:
				/* scan until non-whitespace, then update state */
				if (isspace(c)) continue;
				switch (state) {
					case WS1:
						state = NAME;
						--i;
						break;
					case WS2:
						state = BRA;
						--i;
						break;
					case WS3:
						state = ARG;
						--i;
						args.push_back(std::string());
						break;
					case WS4:
						state = ARG;
						--i;
						break;
					case WS5:
						state = DONE;
						--i;
						break;
					default:
						throw std::logic_error("invalid state");
				}
				break;
			case NAME:
				/* collect until non-alphanumeric character */
				if (isalnum(c) || (c == '_') || (c == '-')) {
					name.push_back(c);
					continue;
				}
				state = WS2;
				--i;
				break;
			case BRA:
				/* error if not opening bracket */
				if (c != '(')
					throw factory_error("unable to parse '" + s + "', " +
										"invalid symbol '" + c + "''");
				state = WS3;
				break;
			case ARG:
				/* collect until space, comma or bracket */
				if (isspace(c)) {
					state = WS4;
					--i;
					break;
				}
				switch (c) {
					case '(': bstack.push_back(')'); goto bopen;
					case '{': bstack.push_back('}'); goto bopen;
					bopen:
						args.back().push_back(c);
						break;
					case ',':
					case ')':
					case '}':
						if (bstack.empty()) {
							state = KET_OR_COMMA;
							--i;
							break;
						}
						/* in a braced subexpression, braces must match, argument continues */
						if (c != ',') {
							if (c != bstack.back())
								throw factory_error("expected '"s + bstack.back() + "' but found '" + c + "'");
							bstack.pop_back();
						}
						args.back().push_back(c);
						continue;
					case '"':
						qinclude = !args.back().empty();
						if (qinclude)
							args.back().push_back(c);
						state = ARG_DQ;
						break;
					default:
						args.back().push_back(c);
						continue;
				}
				break;
			case ARG_DQ:
				/* scan until closing quote, then update state */
				switch (c) {
					case '"':
						if (qinclude)
							args.back().push_back(c);
						state = ARG;
						break;
					case '\\':
						state = ARG_DQE;
						break;
					default:
						args.back().push_back(c);
				}
				break;
			case ARG_DQE:
				/* escaped character, collect unconditionally */
				args.back().push_back(c);
				state = ARG_DQ;
				break;
			case KET_OR_COMMA:
				switch (c) {
					case ')':
						state = WS5;
						break;
					case ',':
						state = WS3;
						break;
					default:
						throw factory_error("unable to parse '" + s + "', " +
											"invalid symbol '" + c + "''");
				}
				break;
			case DONE:
				throw factory_error("unable to parse '" + s + "', " +
									"invalid symbol '" + c + "''");
		}
	}

	if (state != WS5)
		throw factory_error("unable to parse '" + s + "', " +
							"incomplete expression");

	return { name, args };
}

std::ostream& operator<<(std::ostream& o, const parsed_expression_t& expr)
{
	o << expr.first << "(";
	std::size_t i=0;
	for(const std::string& s: expr.second)
		o << (i++ ? "," : "") << s;
	o << ")";
	return o;
}

/**
 * rng. Implicit argument used to pass the RNG to network constructors
 */

rng_t *random_engine = nullptr;

const char *rng::name = nullptr;

const bool rng::implicit = 0;

const std::function<std::string (const rng::value_type&)> rng::renderer =
	[](const std::reference_wrapper<rng_t>&) { return ""; };

template <>
std::pair<std::reference_wrapper<rng_t>, int> argument<rng>(const std::vector<std::string> &vs, size_t i)
{
   return { *random_engine, 0 };
}

template <>
std::string description<rng>()
{
	return "";
}


} /* namespace factories */
