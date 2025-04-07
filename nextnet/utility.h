////  utility.h
//  epidemics_update
//
//  Created by Cure Samuel Cyrus on 2022/02/19.
//
#pragma once

#include "nextnet/stdafx.h"
#include "nextnet/types.h"

inline std::size_t hash_combine(std::size_t seed) { return seed; }

/**
 * @brief Updates the `seed` with the hash values of the other parameters
 *
 * All parameter types must have a working specialization of `std::hash`.
 *
 * Based on: https://stackoverflow.com/questions/2590677/how-do-i-combine-hash-values-in-c0x
 *
 * @param seed seed value to update
 * @param v first object whose hash to update the seed with
 * @param rest next object(s) whose hash to update the seed with
 */
template <typename T, typename... Rest>
inline std::size_t hash_combine(std::size_t seed, const T &v, Rest... rest)
{
    return hash_combine(seed ^ (std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2)), rest...);
}

/**
 * @brief Hashing support for std::pair, required for use in std::unorderedet_map
 *
 * The C++ standard library weirdly enough does not provide this.
 */
struct pair_hash
{
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &p) const
    {
        return hash_combine(0, p.first, p.second);
    }
};

/**
 * @brief Approximates m with p/q under the constraint that q <= L
 * @param x value to approximate
 * @param L largest allowed denominator
 * @return the pair (p, q)
 */
std::pair<unsigned int, unsigned int> fraction(double x, unsigned int L);

/**
 * @brief Inverts a survivial function on [0, infinity), i.e. a decreasing function with f(0)=1, f(infinity)=c > 0.
 *
 * Uses bisection to invert a survival function of a probability distribution. Bisection is numerically stable
 * and accurate, but slower than gradient-based method such as Newton-Raphson.
 *
 * @param u the at wich to evalutate f^-1, i.e. for which to solve u = f(x).
 * @param precision the absolute precision up to which to determine u
 * @param f the function (or other callable object) to invert
 * @param args additional arguments to pass to f
 * @return the value f^-1(u), i.e. x such that f(x, args...) = u.
 */
template <typename T, typename... Args>
double inverse_survival_function(double u, double precision, T f, Args... args)
{
    if ((u > 1) || (u < 0))
        return NAN;
    /* Use bisection to invert f. We assume that f(0) = 1 >= u,
     * so we canstart with the left bound l = 0, but we have to find
     * a right bound r such that  phi(r) <= u */
    double l   = 0;
    double f_l = 1;
    double r   = 1;
    double f_r = f(r, std::forward(args)...);
    while (std::isfinite(r) && (f_r > u)) {
        r *= 2;
        f_r = f(r, std::forward(args)...);
    }
    /* Now we split the interval and pick the left or right subinterval
     * until we reach the desired precision
     */
    while ((f_l != f_r) && ((r - l) > precision)) {
        const double m   = std::isfinite(r) ? (l + r) / 2 : l * 2;
        const double f_m = f(m, std::forward(args)...);
        if (f_m >= u) {
            l   = m;
            f_l = f_m;
        } else {
            r   = m;
            f_r = f_m;
        }
    }
    /* Return mid-point of the final interval */
    return (l + r) / 2;
}

/**
 * @brief edge_index_undirected One-to-one map of edges (i,j) of an undirected graph to {0,...,|E|-1}
 *
 * Translate (i, j) where 0 <= i,j < N into an index 0 <= i < N*(N - 1)/2.
 * We first translate (i, j) into (n1, n2) where n1 > n2. Since there
 * are n1 valid edge descriptors (n1, x) for a node n1, there are
 * 0 + 1 + ... + (n1 - 1) = n * (n1 - 1) / 2 edge descriptors (n1p, x)
 * with n1p < n1. The mapping of (n1, n2) to (n1 * (n1 - 1) / 2) + n2
 * is thus bijective onto [0, E) where E is the number of possible edges
 * in an undirected graph with N nodes.
 *
 * @param i first incident node
 * @param j second incident node
 * @return an integer from {0,...,|E|-1} where |E| is the number of edges
 */
inline unsigned long edge_index_undirected(node_t i, node_t j)
{
    if ((i < 0) || (j < 0) || (i == j))
        throw std::range_error("node (i,j) indices must be non-negative and different");
    const unsigned long n1 = (unsigned long)std::max(i, j);
    const unsigned long n2 = (unsigned long)std::min(i, j);
    assert(n1 > n2);
    const unsigned long k = (n1 * (n1 - 1) / 2) + n2;
    assert(k >= 0);
    return k;
}

/**
 * @brief edge_index One-to-on map of {0,...,|E|-1} onto edges (i,j) of an undirected graph
 *
 * @param e edge index from {0,...,|E|-1} where |E| is the number of edges
 * @return an undirected edge (i,j), normalized such that i > j.
 */

inline std::pair<node_t, node_t> edge_undirected(unsigned long k)
{
    // Check that double precision suffices to distinguish edge indices
    if (8.0 * (double)k == (8.0 * (double)k + 1.0))
        throw std::range_error("edge index out of range");
    // Find n1 by inverting edge_index_undirected()
    const double n1p = std::floor(0.5 + 0.5 * std::pow(1 + 8 * k, 0.5));
    if (n1p > (double)std::numeric_limits<node_t>::max())
        throw std::range_error("edge index out of range");
    const node_t n1 = (node_t)n1p;
    const node_t n2 = (node_t)(k - (unsigned long)n1p);
    assert(n1 > n2);
    return std::make_pair(n1, n2);
}

/**
 * @brief Represents a set of integers the allows random elements and non-elements to be sampled
 *
 * Operations `insert()` and `erase()` take time O(log n) where n is the number of elements in the set.
 * Operations `draw_present()` and `draw_absent()` take time O(k) where k is the number of
 * consecutive ranges of elements in the set.
 *
 * NOTE: The performance of `draw_element()` and `draw_complement()` could be improved to
 * almost O(1) by using a dynamic distribution to draw the range. This is a possible future improvement.
 */
template <typename T>
struct integer_set
{
    /**
     * @brief The type of the set elements. Must be an integer type.
     */
    typedef T element_type;
    typedef T value_type;

private:
    /**
     * @brief A range of consecutive elements
     */
    struct range
    {
        range()
            : first(0)
            , last(0)
        {
        }

        range(element_type f, element_type l)
            : first(f)
            , last(l)
        {
        }

        element_type first;
        element_type last;
    };

    /**
     * @brief Comparisom functor for `range`
     *
     * Sorts lexicographically by (`last`, `first`). Also allows comparison with single integers, in
     * which case they are compared against `last`.
     */
    struct range_cmp
    {
        typedef void is_transparent;

        bool operator()(const range &a, const range &b) const
        {
            return (a.last < b.last) || ((a.last == b.last) && (a.first < b.first));
        }

        bool operator()(const range &a, const element_type b) const
        {
            return (a.last < b);
        }
    };

    element_type lowerbound = std::numeric_limits<element_type>::min() + 1;
    element_type upperbound = std::numeric_limits<element_type>::max() - 1;
    std::size_t ntotal      = 0;
    std::set<range, range_cmp> ranges;

public:
    /**
     * @brief Interator for interating over the members of a `integer_set`
     *
     * Conforms to BidirectionalIterator. Only a const iterator is provided since
     * the semantics of iterating while replacing individual elements of a
     * `integer_set` are not well-defined.
     */
    struct const_iterator
    {
        typedef element_type value_type;
        typedef const element_type &reference;
        typedef void difference_type;
        typedef void pointer;
        typedef std::bidirectional_iterator_tag iterator_category;

        typedef typename std::set<range>::const_iterator range_iterator;

        const_iterator(range_iterator rng, std::size_t idx)
            : current_range(rng)
            , current_index(idx)
        {
        }

        value_type operator*() const { return current_range->first + (element_type)current_index; }

        const_iterator &operator++()
        {
            const std::size_t last_index = (current_range->last - current_range->first);
            if (current_index == last_index) {
                ++current_range;
                current_index = 0;
            } else
                ++current_index;
            return *this;
        }
        const_iterator operator++(int)
        {
            const_iterator r = *this;
            *this ++;
            return r;
        }

        const_iterator &operator--()
        {
            if (current_index == 0) {
                --current_range;
                current_index = (current_range->last - current_range->first);
            } else
                --current_index;
            return *this;
        }
        const_iterator operator--(int)
        {
            const_iterator r = *this;
            *this --;
            return r;
        }

        bool operator==(const const_iterator &other) const
        {
            return (this->current_range == other.current_range) &&
                   (this->current_index == other.current_index);
        }

        bool operator!=(const const_iterator &other) const
        {
            return !(*this == other);
        }

    private:
        friend struct integer_set;

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
     * @brief Creates an empty `integer_set` whose domain is maximal
     *
     * Note that the maximal domain is [min+1, max-1] where min and max
     * are the minimal and maximal finite value that the underlying integer
     * type can hold.
     */
    integer_set() {}

    /**
     * @brief Creates an empty `integer_set` over the domain [`lb`, `ub`]
     *
     * Note that lb and ub are restricted to lb >= max+1 and ub <= min-1
     * where min and max are the minimal and maximal finite value that the
     * underlying integer type can hold.
     *
     */
    integer_set(element_type lb, element_type ub)
        : lowerbound(lb)
        , upperbound(ub)
    {
        if (lowerbound > upperbound)
            throw std::range_error("lower bound exceeds upper bound");
        if (lowerbound == std::numeric_limits<element_type>::min())
            throw std::range_error("lower bound must be larget than minimal value of underlying integer type");
        if (upperbound == std::numeric_limits<element_type>::max())
            throw std::range_error("upper bound must be smaller than maximal value of underlying integer type");
    }

    /**
     * @brief Creates a set over maximal domain containg the elements between iterators `first` and `last`.
     */
    template <typename It>
    integer_set(It first, It last)
    {
        std::copy(first, last, std::inserter(*this, end()));
    }

    /**
     * @brief Creates a set with domain [`lb`, `ub`] containg the elements between iterators `first` and `last`.
     */
    template <typename It>
    integer_set(element_type lb, element_type ub, It first, It last)
        : integer_set(lb, ub)
    {
        std::copy(first, last, std::inserter(*this, end()));
    }

    /**
     * @brief Creates a set with maximal domain containg the specified elements
     */
    integer_set(std::initializer_list<T> init)
        : integer_set(init.begin(), init.end())
    {
    }

    /**
     * @brief Creates a set with domain [`lb`, `ub`] containg the specified elements
     */
    integer_set(element_type lb, element_type ub, std::initializer_list<T> init)
        : integer_set(lb, ub, init.begin(), init.end())
    {
    }

    /**
     * @brief Number of elements of the set
     */
    std::size_t size() const
    {
        return ntotal;
    }

    /**
     * @brief Inserts the interger into the set unless it is already a member
     * @return `true` if the integer was inserted, `false` if it was already a member.
     * @throw `std::range_error` if the element is too large or small to be inserted
     */
    std::pair<const_iterator, bool> insert(element_type e)
    {
        /* The possible cases of existing ranges are:
         *    I: [a,     e-1]  [c > e+1, d] ==> extend [a, e-1] to [a, e]
         *   II: [a, b < e-1]  [e+1    , d] ==> extend [e+1, d] to [e, d]
         *  III: [a,     e-1]  [e+1    , d] ==> merge [a, e-1] and [e+1, d] to [a, d]
         *  IVa: [a, b < e-1]               ==> insert [e, e]
         *  IVb: [a, b < e-1]  [c > e+1, d] ==> insert [e, e]
         *   Va: [a <= e-1, b >= e]         ==> do nothing
         *   Vb: [e       , b >= e]         ==> do nothing
         */

        if ((e < lowerbound) || (e > upperbound))
            throw std::range_error("element lies outside of the set's domain");

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
    const_iterator insert(const_iterator hint, element_type e)
    {
        return insert(e).first;
    }

    /**
     * @brief Removes the interger from the set if it is currently  a member
     * @return Number of elements removed (0 or 1)
     * @throw `std::range_error` if the element is too large or small to be contained in the set.
     */
    std::size_t erase(element_type e)
    {
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
                const range r(f, e - 1);
                ranges.insert(r);
            }
            ntotal--;
            return 1;
        }

        /* Ib. Restrict range */
        if (i->first == e) {
            const range r(e + 1, i->last);
            ranges.erase(i);
            ranges.insert(r);
            ntotal--;
            return 1;
        }

        /* II. Split range */
        assert(i->first <= e - 1);
        assert(i->last >= e + 1);
        const range r1(i->first, e - 1);
        const range r2(e + 1, i->last);
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
    const_iterator erase(const const_iterator &i)
    {
        return erase(i, std::next(i));
    }

    /**
     * @brief Removes the elements between the two iterators ([first, last)]
     * @return Iterator to the first element past the last removed one (can be end()).
     */
    const_iterator erase(const_iterator a, const const_iterator &b)
    {
        /* Nothing to do for empty range of elements */
        if (a == b)
            return b;

        while (a.current_range != b.current_range) {
            /* Remove remaining part of current range of a */
            const range r(a.current_range->first, *a - 1);
            const auto i2 = ranges.erase(a.current_range);
            if (r.first <= r.last)
                ranges.insert(r);
            /* Let iterator point to first element of next range */
            a = const_iterator(i2, 0);
        }

        /* Sanity check */
        if (a.current_index >= b.current_index)
            std::logic_error("iterators don't form a valid range of elements");

        /* If b point to the first element of a range or end(), we're done */
        if (b.current_index == 0)
            return b;

        /* Remove part between a (inclusive) and b (exclusive) by splitting */
        const range r1(a.current_range->first, *a - 1);
        const range r2(*b, b.current_range->last);
        assert(r1.last + 1 < r2.first);
        ranges.erase(a.current_range);
        if (r1.first <= r1.last)
            ranges.insert(r1);
        const auto i2 = ranges.insert(r2).first;
        return const_iterator(i2, 0);
    }

    /**
     * @brief Draws a random element from the set
     * @return The element selected
     * @throw `std::runtime_error` if the set is empty
     */
    T draw_element(rng_t &engine) const
    {
        if (ntotal == 0)
            throw std::runtime_error("cannot draw from an empty set");

        /* Draw random range, with probabilities prop. to the ranges' lengths */
        const double p = std::uniform_real_distribution<>(0, 1)(engine);
        std::size_t l  = 0;
        for (const auto &r : ranges) {
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
     * @brief Draws a random integer that is not an element of the set
     * @return The integer selected
     * @throw `std::runtime_error` if the set contains all integers [lb, ub].
     *        `std::range_error` if the lower or upper bound is invalid.
     */
    T draw_complement(rng_t &engine) const
    {
        if (ntotal == (std::size_t)(upperbound - lowerbound + 1))
            throw std::runtime_error("cannot draw from an empty complement");

        /* Draw random gap between ranges, with probabilities prop. to the gaps' lengths */
        const double p              = std::uniform_real_distribution<>(0, 1)(engine);
        const std::size_t size_gaps = ((std::size_t)(upperbound - lowerbound) + 1 - ntotal);
        auto i                      = ranges.begin();
        const auto e                = ranges.end();
        std::size_t l               = 0;
        T a                         = lowerbound;
        bool done                   = false;
        do {
            /* Gap extends to beginning of next range or upper bound for last range */
            range gap;
            if (i != e) {
                gap = range(a, i->first - 1);
                a   = i->last + 1;
                i++;
            } else {
                gap  = range(a, upperbound);
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
    const_iterator find(element_type e) const
    {
        /* The possible cases are:
         *   Ia:                  ==> not found
         *   Ib: [a > e, b > e]   ==> not found
         *   II: [a <= e, b >= e] ==> found
         */

        /* Outside the set's range */
        if ((e < lowerbound) || (e > upperbound))
            return end();

        /* Find the first range [x, y] that does not wholly lie left of e, i.e. where y >= e*/
        const auto i = ranges.lower_bound(e);

        /* Ia or Ib. Not found */
        if ((i == ranges.end()) || (i->first > e))
            return end();

        /* II. Found */
        return const_iterator(i, (std::size_t)(e - i->first));
    }
};

/**
 * @brief Sets which elemenets to be accessed by index and randomly drawn
 *
 * Combined a std::vector with an std::unordered_set, i.e. elements can be accessed by
 * index but also efficiently be located by value. Elements can also be drawn randomly
 * and uniformly. his class provides a drop-in replacement for std::unordered_set that
 * provides operator[](int) which returns the i-th object, and operator()(rng) method which
 * returns an interator pointing to a randomly and unformly chosen element.
 *
 * Internally, elements are stored in a vector, and an unordered_map is used as an index
 * to efficiently find elements by value and to guarantee that uniqueness.
 */
template <typename T, class Hash = std::hash<T>, class KeyEqual = std::equal_to<T>>
class indexed_set
{
public:
    typedef T key_type;
    typedef T value_type;
    typedef std::size_t size_type;
    typedef Hash hasher;
    typedef KeyEqual key_equal;
    typedef value_type &reference;
    typedef const value_type &const_reference;
    typedef value_type *pointer;
    typedef const value_type *const_pointer;

    typedef std::vector<T> vector_type;
    typedef typename vector_type::iterator iterator;
    typedef typename vector_type::const_iterator const_iterator;

    indexed_set(){};

    template <class InputIt>
    indexed_set(InputIt first, const InputIt last)
    {
        for (; first != last; ++first) {
            insert(*first);
        }
    }

    size_type size() const { return elements.size(); }

    value_type &operator[](size_type i) { return elements[i]; }
    const value_type &operator[](size_type i) const { return elements[i]; }

    iterator begin() { return elements.begin(); }
    const_iterator begin() const { return elements.begin(); }

    iterator end() { return elements.end(); }
    const_iterator end() const { return elements.end(); }

    iterator find(const value_type &value)
    {
        const auto i = index.find(value);
        if (i == index.end())
            return end();
        return elements.begin() + i->second;
    }

    const_iterator find(const value_type &value) const
    {
        const auto i = index.find(value);
        if (i == index.end())
            return end();
        return elements.begin() + i->second;
    }

    std::pair<iterator, bool> insert(const value_type &value)
    {
        iterator i = find(value);
        if (i != end())
            return { i, false };
        elements.push_back(value);
        index[value] = elements.size() - 1;
        return { elements.end() - 1, true };
    }

    iterator erase(iterator it)
    {
        const std::size_t idx = it - elements.begin();
        // Remove index entry
        index.erase(*it);
        // Swap element with last element in vector
        if (idx != elements.size() - 1) {
            using std::swap;
            swap(elements[idx], elements.back());
            // Update index entry for former last element
            index[elements[idx]] = idx;
        }
        // Remove last element
        elements.pop_back();
        // Return iterator pointing to element after removed one
        return elements.begin() + idx;
    }

    size_type erase(const key_type &key)
    {
        iterator it = find(key);
        if (it == end())
            return 0;
        erase(it);
        return 1;
    }

    void clear()
    {
        elements.clear();
        index.clear();
    }

    template <typename Rng>
    iterator operator()(Rng &rng)
    {
        return elements.begin() + std::uniform_int_distribution<std::size_t>(0, elements.size() - 1)(rng);
    }

private:
    vector_type elements;
    std::unordered_map<T, std::size_t, hasher, key_equal> index;
};

/**
 * @brief Maps whose elemenets to be accessed by index and randomly drawn
 *
 * Combined a std::vector with an std::unordered_map, i.e. elements can be accessed by
 * index but also efficiently be located by key. Elements can also be drawn randomly
 * and uniformly. his class provides a drop-in replacement for std::unordered_map that
 * provides operator[](int) which returns the i-th object, and operator()(rng) method which
 * returns an interator pointing to a randomly and unformly chosen element.
 *
 * Internally, elements are stored in a vector, and an unordered_map is used as an index
 * to efficiently find elements by value and to guarantee that uniqueness.
 */
template <typename K, typename M, class Hash = std::hash<K>, class KeyEqual = std::equal_to<K>>
class indexed_map
{
public:
    typedef K key_type;
    typedef M mapped_type;
    typedef std::pair<const K, M> value_type;
    typedef std::size_t size_type;
    typedef Hash hasher;
    typedef KeyEqual key_equal;
    typedef value_type &reference;
    typedef const value_type &const_reference;
    typedef value_type *pointer;
    typedef const value_type *const_pointer;

    typedef std::vector<value_type> vector_type;
    typedef typename vector_type::iterator iterator;
    typedef typename vector_type::const_iterator const_iterator;

    indexed_map(){};

    template <class InputIt>
    indexed_map(InputIt first, const InputIt last)
    {
        for (; first != last; ++first) {
            insert(*first);
        }
    }

    size_type size() const { return elements.size(); }

    const value_type &operator[](size_type i) const { return elements[i]; }

    mapped_type &operator[](key_type k)
    {
        auto i = find(k);
        if (i != end())
            return i->second;
        return insert({ k, mapped_type() }).first->second;
    }

    iterator begin() { return elements.begin(); }
    iterator end() { return elements.end(); }
    const_iterator begin() const { return elements.begin(); }
    const_iterator end() const { return elements.end(); }

    iterator find(const key_type &value)
    {
        const auto i = index.find(value);
        if (i == index.end())
            return end();
        return elements.begin() + i->second;
    }

    const iterator find(const key_type &value) const
    {
        const auto i = index.find(value);
        if (i == index.end())
            return end();
        return elements.begin() + i->second;
    }

    std::pair<iterator, bool> insert(const value_type &value)
    {
        iterator i = find(value.first);
        if (i != end())
            return { i, false };
        elements.push_back(value);
        index[value.first] = elements.size() - 1;
        return { elements.end() - 1, true };
    }

    iterator erase(iterator it)
    {
        const std::size_t idx = it - elements.begin();
        // Remove index entry
        index.erase(it->first);
        // Swap element with last element in vector
        if (idx != elements.size() - 1) {
            using std::swap;
            // The first element in value_type is const since
            // we the key is immutable. Here, we have to override
            // that though so we const_cast.
            swap(const_cast<key_type &>(elements[idx].first),
                 const_cast<key_type &>(elements.back().first));
            swap(elements[idx].second, elements.back().second);
            // Update index entry for former last element
            index[elements[idx].first] = idx;
        }
        // Remove last element
        elements.pop_back();
        // Return iterator pointing to element after removed one
        return elements.begin() + idx;
    }

    size_type erase(const key_type &key)
    {
        iterator it = find(key);
        if (it == end())
            return 0;
        erase(it);
        return 1;
    }

    void clear()
    {
        elements.clear();
        index.clear();
    }

    template <typename Rng>
    iterator operator()(Rng &rng)
    {
        return elements.begin() + std::uniform_int_distribution<std::size_t>(0, elements.size() - 1)(rng);
    }

private:
    vector_type elements;
    std::unordered_map<key_type, std::size_t, hasher, key_equal> index;
};
