//
// Created by etien on 31/01/2022.
//

#include "miller-rabin-gmp.h"

gmp_randclass* generate_prng()
{
    auto *prng = new gmp_randclass(gmp_randinit_default);
    prng->seed(time(NULL));
    return prng;
}

mpz_class randint(const mpz_class& lowest, const mpz_class& highest, gmp_randclass* prng)
{
    if ( lowest == highest )
        return lowest;

    return prng->get_z_range(highest - lowest + 1) + lowest;
}

mpz_class pow_mod(mpz_class a, mpz_class x, const mpz_class& n)
{
    mpz_class r = 1;

    while (x > 0) {
        if ((x & 1) == 1) {
            r = a*r % n;
        }
        x >>= 1;
        a = a*a % n;
    }

    return r;
}

bool miller_rabin_backend(const mpz_class& n, const size_t rounds, gmp_randclass* prng)
{
    // Treat n==1, 2, 3 as a primes
    if (n == 1 || n == 2 || n == 3)
        return true;

    // Treat negative numbers in the frontend
    if (n <= 0)
        return false;

    // Even numbers larger than two cannot be prime
    if ((n & 1) == 0)
        return false;

    // Write n-1 as d*2^s by factoring powers of 2 from n-1
    size_t s = 0;
    {
        mpz_class m = n - 1;
        while ((m & 1) == 0) {
            ++s;
            m >>= 1;
        }
    }
    const mpz_class d = (n - 1) / (mpz_class(1) << s);

    for (size_t i = 0; i < rounds; ++i) {
        const mpz_class a = randint(2, n - 2, prng);
        mpz_class x = pow_mod(a, d, n);

        if (x ==1 || x == (n - 1))
            continue;

        for (size_t r = 0; r < (s-1); ++r) {
            x = pow_mod(x, 2, n);
            if (x == 1) {
                // Definitely not a prime
                return false;
            }
            if (x == n - 1)
                break;
        }

        if (x != (n - 1)) {
            // Definitely not a prime
            return false;
        }
    }

    // Might be prime
    return true;
}

/*
 * The Miller-Rabin front end.
 */
bool prob_prime(const mpz_class& n, const size_t rounds, gmp_randclass* prng)
{
    return miller_rabin_backend(n > 0 ? n : -n, rounds, prng);
}