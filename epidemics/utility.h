////  utility.h
//  epidemics_update
//
//  Created by Cure Samuel Cyrus on 2022/02/19.
//
#pragma once

#include "stdafx.h"
#include "types.h"

/**
 * @brief Hashing support for std::pair, required for use in std::unorderedet_map
 * 
 * The C++ standard library weirdly enough does not provide this. Simply
 * XOR-ing the two hashes is pretty silly, for it seems sufficient for the
 * std::unordered_map declared below.
 */
struct pair_hash {
    template<typename T, typename U>
    std::size_t operator () (const std::pair<T,U> &p) const {
        return std::hash<T>{}(p.first) ^ std::hash<U>{}(p.second);
    }
};  

/**
 * @brief Approximates m with p/q under the constraint that q <= L
 * @param x value to approximate
 * @param L largest allowed denominator
 * @return the pair (p, q)
 */
std::pair<unsigned int, unsigned int> fraction(double x, unsigned int L);

template<typename T, typename ...Args>
double inverse_survival_function(double u, double precision, T f, Args... args) {
    if ((u > 1) || (u < 0))
        return NAN;
    /* Use bisection to invert f. We assume that f(0) = 1 >= u,
     * so we canstart with the left bound l = 0, but we have to find
     * a right bound r such that  phi(r) <= u */
    double l = 0;
    double f_l = 1;
    double r = 1;
    double f_r = f(r, std::forward(args)...);
    while (std::isfinite(r) && (f_r > u)) {
        r *= 2;
        f_r = f(r, std::forward(args)...);
    }
    /* Now we split the interval and pick the left or right subinterval
     * until we reach the desired precision
     */
    while ((f_l != f_r) && ((r - l) > precision)) {
        const double m = std::isfinite(r) ? (l + r) / 2 : l*2 ;
        const double f_m = f(m, std::forward(args)...);
        if (f_m >= u) {
            l = m;
            f_l = f_m;
        }
        else {
            r = m;
            f_r = f_m;
        }
    }
    /* Return mid-point of the final interval */
    return (l + r) / 2;
}

/**
 * @brief Represents a set of integers the allows random elements and non-elements to be sampled
 *
 * Operations `insert()` and `erase()` take time O(log n) where n is the number of elements in the set.
 * Operations `draw_present()` and `draw_absent()` take time O(k) where k is the number of
 * consecutive ranges of elements in the set.
 *
 * NOTE: The performance of `draw_present()` and `draw_absent()` could be improved to
 * almost O(1) by using a dynamic distribution to draw the range. This is a possible future improvement.
 */
template<typename T>
struct integer_set {
	/**
	 * @brief The type of the set elements. Must be an integer type.
	 */
	typedef T element_type;
	typedef T value_type;
	
private:
	/**
	 * @brief A range of consecutive elements
	 */
	struct range {
		range() :first(0), last(0) {}
		
		range(element_type f, element_type l)
			:first(f), last(l)
		{}

		element_type first;
		element_type last;		
	};
	
	/**
	 * @brief Comparisom functor for `range`
	 *
	 * Sorts lexicographically by (`last`, `first`). Also allows comparison with single integers, in
	 * which case they are compared against `last`.
	 */
	struct range_cmp {
		typedef void is_transparent;
		
		bool operator()(const range& a, const range& b) const {
			return (a.last < b.last) || ((a.last == b.last) && (a.first < b.first));
		}
		
		bool operator()(const range& a, const element_type b) const {
			return (a.last < b);
		}
	};

	std::size_t ntotal = 0;
	std::set<range, range_cmp> ranges;
	
public:
	/**
	 * @brief Interator for interating over the members of a `integer_set`
	 *
	 * Conforms to BidirectionalIterator. Only a const iterator is provided since
	 * the semantics of iterating while replacing individual elements of a
	 * `integer_set` are not well-defined.
	 */
	struct const_iterator {
		typedef typename std::set<range>::const_iterator range_iterator;
		
		const_iterator(range_iterator rng, std::size_t idx)
			:current_range(rng), current_index(idx)
		{}
		
		T operator*() const { return current_range->first + (element_type)current_index; }
		
		const_iterator& operator++() {
			const std::size_t last_index = (current_range->last - current_range->first);
			if (current_index == last_index) {
				++current_range;
				current_index = 0;
			} else
				++current_index;
			return *this;
		}
		const_iterator operator++(int) { const_iterator r = *this; *this++; return r; }

		const_iterator& operator--() {
			if (current_index == 0) {
				--current_range;
				current_index = (current_range->last - current_range->first);
			} else
				--current_index;
			return *this;
		}
		const_iterator operator--(int) { const_iterator r = *this; *this--; return r; }
		
		bool operator==(const const_iterator& other) const {
			return (this->current_range == other.current_range) &&
				   (this->current_index == other.current_index);
		}
		
		bool operator!=(const const_iterator& other) const {
			return !(*this == other);
		}
		
	private:
		range_iterator current_range;
		std::size_t current_index;
	};
	
	/**
	 * @brief For compatibility with STL algorithms
	 */
	using iterator = const_iterator;
	
	/**
	 * @brief Returns an iterator to the first element.
	 */
	const_iterator begin() const { return const_iterator(ranges.begin(), 0); }

	/**
	 * @brief Returns an iterator past the last element.
	 */
	const_iterator end() const { return const_iterator(ranges.end(), 0); }

	/**
	 * @brief Creates an empty `integer_set`
	 */
	integer_set() {};

	/**
	 * @brief Creates a set containg the elements between iterators `first` and `last`.
	 */
	template<typename It>
	integer_set(It first, It last) {
		std::copy(first, last, std::inserter(*this, end()));
	}

	/**
	 * @brief Creates a set containg the specified elements
	 */
	integer_set(std::initializer_list<T> init)
		:integer_set(init.begin(), init.end())
	{}

	/**
	 * @brief Number of elements of the set
	 */
    std::size_t size() const {
        return ntotal;
    }

	/**
	 * @brief Inserts the interger into the set unless it is already a member
	 * @return `true` if the integer was inserted, `false` if it was already a member.
	 * @throw `std::range_error` if the element is too large or small to be inserted
	 */
    std::pair<const_iterator, bool> insert(element_type e) {
        /* The possible cases of existing ranges are:
         *    I: [a,     e-1]  [c > e+1, d] ==> extend [a, e-1] to [a, e]
         *   II: [a, b < e-1]  [e+1    , d] ==> extend [e+1, d] to [e, d]
         *  III: [a,     e-1]  [e+1    , d] ==> merge [a, e-1] and [e+1, d] to [a, d]
         *  IVa: [a, b < e-1]               ==> insert [e, e]
         *  IVb: [a, b < e-1]  [c > e+1, d] ==> insert [e, e]
         *   Va: [a <= e-1, b >= e]         ==> do nothing
         *   Vb: [e       , b >= e]         ==> do nothing
         */

        if ((e == std::numeric_limits<T>::min()) || (e == std::numeric_limits<T>::max()))
            throw std::range_error("element in a integer_set cannot be maximum or minimum of type's domain");

        /* Find the first range [x, y] that does not wholly lie left of e - 1, i.e. where y >= e - 1*/
        const auto i = ranges.lower_bound(e - 1);

        /* IVa or IVb. Just insert [e, e] */
        if ((i == ranges.end()) || (i->first > e + 1)) {
            const range r(e, e);
            const auto i2 = ranges.insert(r).first;
			ntotal++;
			return { const_iterator(i2, 0), true };
        }

        /* I, III or Va. If the range contains e - 1, extend or merge it */
        if (i->first <= e - 1) {
            /* Va. If it already extends past e, we're done */
            if (i->last >= e)
                return { const_iterator(i, (std::size_t)(e - i->first)), false };
            assert(i->last == e - 1);

            /* III. Before we extend the range, check whether we would create two adjacent ranges */
            const auto s = std::next(i);
            if ((s != ranges.end()) && (s->first == e + 1)) {
                /* Next range is adjacent, i.e. we have [a, e - 1], [ e, b ]. Merge ranges. */
                const range r(i->first, s->last);
                ranges.erase(i);
                ranges.erase(s);
                const auto i2 = ranges.insert(r).first;
				ntotal++;
                return { const_iterator(i2, (std::size_t)(e - i2->first)), true };
            }

            /* I. Finally, extend the range */
            const range r(i->first, e);
            ranges.erase(i);
            const auto i2 = ranges.insert(r).first;
			ntotal++;
			return { const_iterator(i2, (std::size_t)(e - i2->first)), true };
        }

        /* Vb. The range contains e, we're done */
        if (i->first <= e) {
            assert(i->last >= e);
			return { const_iterator(i, (std::size_t)(e - i->first)), false };
        }

        /* II. Extend the range. */
        assert(i->first == e + 1);
        const range r(e, i->last);
        ranges.erase(i);
        const auto i2 = ranges.insert(r).first;
		ntotal++;
		return { const_iterator(i2, (std::size_t)(e - i2->first)), true };
    }
	
	/**
	 * @brief Inserts the interger into the set unless it is already a member. The hint is ignored.
	 * @return `true` if the integer was inserted, `false` if it was already a member.
	 * @throw `std::range_error` if the element is too large or small to be inserted
	 *
	 * This exists only to be compatible with `std::set`. The hint for the insertion position is currently ignored.
	 */
	const_iterator insert(const_iterator hint, element_type e) {
		return insert(e).first;
	}
	
	/**
	 * @brief Removes the interger from the set if it is currently  a member
	 * @return Number of elements removed (0 or 1)
	 * @throw `std::range_error` if the element is too large or small to be contained in the set.
	 */
    std::size_t erase(element_type e) {
        /* The possible cases are:
         *    Ia: [a < e, e]      ==> restrict to [a, e-1]
         *    Ib: [e, b > e]      ==> restrict to [e+1, b]
         *    Ic: [e, e]          ==> remove
         *    II: [a < e, b > e]  ==> split into [a, e-1], [e+1, b]
         *  IIIa: [a > e, b]      ==> do nothing
         *  IIIb: [a, b < e]      ==> do nothing
         */

        if ((e == std::numeric_limits<T>::min()) || (e == std::numeric_limits<T>::max()))
            throw std::range_error("element in a integer_set cannot be maximum or minimum of type's domain");

        /* Find the first range [x, y] that does not wholly lie left of e, i.e. where y >= e*/
        const auto i = ranges.lower_bound(e);

        /* IIIa or IIIb. Do nothing */
        if ((i == ranges.end()) || (i->first > e))
            return 0;
        assert(i->last >= e);

        /* Ia or Ic. Restrict range */
        if (i->last == e) {
			const element_type f = i->first;
            ranges.erase(i);
            /* Ia. Insert restricted range if non-empty */
            if (f < e) {
                const range r(f, e-1);
                ranges.insert(r);
            }
			ntotal--;
            return 1;
        }

        /* Ib. Restrict range */
        if (i->first == e) {
            const range r(e+1, i->last);
            ranges.erase(i);
            ranges.insert(r);
			ntotal--;
            return 1;
        }

        /* II. Split range */
        assert(i->first <= e-1);
        assert(i->last >= e+1);
        const range r1(i->first, e-1);
        const range r2(e+1, i->last);
        ranges.erase(i);
        ranges.insert(r1);
        ranges.insert(r2);
		ntotal--;
        return 1;
    }

	/**
	 * @brief Removes the element pointed to by the iterator
	 * @return Iterator to the first element past the removed one (can be end()).
	 */
	const_iterator erase(const const_iterator& i) {
		const const_iterator j = i;
		++j;
		return erase(i, j);
	}

	/**
	 * @brief Removes the elements between the two iterators ([first, last)]
	 * @return Iterator to the first element past the last removed one (can be end()).
	 */
	const_iterator erase(const const_iterator& first, const const_iterator& last) {
		throw std::logic_error("implement me");
	}
	
	/**
	 * @brief Draws a random element from the set
	 * @return The element selected
	 * @throw `std::runtime_error` if the set is empty
	 */
    T draw_present(rng_t& engine) {
		if (ntotal == 0)
			throw std::runtime_error("cannot draw a present element from an empty integer_set");
		
        /* Draw random range, with probabilities prop. to the ranges' lengths */
        const double p = std::uniform_real_distribution<>(0, 1)(engine);
        std::size_t l = 0;
        for(const auto& r: ranges) {
            /* Compute fraction q covered by ranges up to the current one */
            l += (r.last - r.first + 1);
            const double q = (double)l / ntotal;
            /* If that fraction exceeds p, we found our range */
            if (q > p) {
                /* Pick uniformly random from the range */
                return std::uniform_int_distribution<T>(r.first, r.last)(engine);
            }
        }

        throw std::logic_error("draw_present() failed due to internal inconsistency");
    }

	/**
	 * @brief Draws a random integer [lb, ub] that is not an element of the set
	 * @return The integer selected
	 * @throw `std::runtime_error` if the set contains all integers [lb, ub].
	 *        `std::range_error` if the lower or upper bound is invalid.
	 */
    T draw_absent(element_type lb, element_type ub, rng_t& engine) {
		if (lb > ub)
			throw std::range_error("lower bound exceeds upper bound");
		if (!ranges.empty() && (ranges.begin()->first < lb))
            throw std::range_error("lower bound is larger than smallest element of integer_set");
        if (!ranges.empty() && (ranges.rbegin()->last > ub))
            throw std::range_error("upper bound is smaller than largest element of integer_set");
		if (ntotal == (ub - lb + 1))
			throw std::runtime_error("cannot draw from empty integer_set");

        /* Draw random gap between ranges, with probabilities prop. to the gaps' lengths */
        const double p = std::uniform_real_distribution<>(0, 1)(engine);
        const std::size_t size_gaps = ((std::size_t)(ub - lb) + 1 - ntotal);
		auto i = ranges.begin();
		const auto e = ranges.end();
        std::size_t l = 0;
        T a = lb;
		bool done = false;
        do {
            /* Gap extends to beginning of next range or upper bound for last range */
            range gap;
            if (i != e) {
                gap = range(a, i->first - 1);
                a = i->last + 1;
                i++;
            } else {
                gap = range(a, ub);
                done = true;
            }
            /* Compute fraction q covered by gaps up to the current one */
            l += (gap.last - gap.first + 1);
            const double q = (double)l / size_gaps;
            /* If that fraction exceeds p, we found our range */
            if (q > p) {
                /* Pick uniformly random from the gap */
                return std::uniform_int_distribution<T>(gap.first, gap.last)(engine);
            }
        } while (!done);

        throw std::logic_error("draw_present() failed due to internal inconsistency");
    }

	/**
	 * @brief Finds the specified element in the set.
	 * @return An iterator pointing to the element, or `end()` if the element is not in the set.
	 */
	const_iterator find(element_type e) const {
		/* The possible cases are:
		 *   Ia:                  ==> not found
		 *   Ib: [a > e, b > e]   ==> not found
		 *   II: [a <= e, b >= e] ==> found
		 */
		
		/* Find the first range [x, y] that does not wholly lie left of e, i.e. where y >= e*/
		const auto i = ranges.lower_bound(e);
		
		/* Ia or Ib. Not found */
		if ((i == ranges.end()) || (i->first > e))
			return end();
		
		/* II. Found */
		return const_iterator(i, (std::size_t)(e - i->first));
	}
};
