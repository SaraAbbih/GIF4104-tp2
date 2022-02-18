//
// Created by etien on 31/01/2022.
//

#ifndef TP1_V2_MILLER_RABIN_GMP_H
#define TP1_V2_MILLER_RABIN_GMP_H

#include <gmpxx.h>

gmp_randclass* generate_prng();
mpz_class randint(const mpz_class& lowest, const mpz_class& highest, gmp_randclass* prng);
mpz_class pow_mod(mpz_class a, mpz_class x, const mpz_class& n);
bool prob_prime(const mpz_class& n, const size_t rounds, gmp_randclass* prng);

#endif //TP1_V2_MILLER_RABIN_GMP_H
