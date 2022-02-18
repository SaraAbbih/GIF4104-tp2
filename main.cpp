#include <iostream>
#include <math.h>
#include <pthread.h>
#include <fstream>
#include <utility>
#include <vector>
#include <gmpxx.h>
#include <filesystem>

#include "miller-rabin-gmp.h"
#include "Chrono.hpp"

typedef struct {
    mpz_class min;
    mpz_class max;
} Range;

pthread_mutex_t gGetMutex = PTHREAD_MUTEX_INITIALIZER;

int gRangeCounter = 0;
int gInRangeCounter = 0;
std::vector<Range> gRanges;

bool getNextNumber(mpz_class& oNextNumber)
{
    pthread_mutex_lock(&gGetMutex);

    while (true)
    {
        if (gRangeCounter >= gRanges.size()) {
            pthread_mutex_unlock(&gGetMutex);
            return false;
        }

        mpz_class lNextNumber = gRanges[gRangeCounter].min + gInRangeCounter;

        if (lNextNumber > gRanges[gRangeCounter].max) {
            gRangeCounter++;
            gInRangeCounter = 0;
            continue;
        }

        if ((lNextNumber & 1) == 0) {
            gInRangeCounter++;
            continue;
        }

        oNextNumber = lNextNumber;
        gInRangeCounter += 2;

        pthread_mutex_unlock(&gGetMutex);
        return true;
    }
}

void* computePrimes(void* iArg)
{
    auto lPrng = generate_prng();

    auto lThFoundPrimes = new std::vector<mpz_class>;

    mpz_class lEvaluatedNumber;
    while (getNextNumber(lEvaluatedNumber))
        if (prob_prime(lEvaluatedNumber, 3, lPrng))
            lThFoundPrimes->push_back(lEvaluatedNumber);

    pthread_exit((void*)lThFoundPrimes);
}

std::vector<Range> getRangesFromFile(const std::string& iFilename)
{
    std::vector<Range> lRanges;

    std::ifstream lFile(iFilename);
    std::string lLine;
    while (std::getline(lFile, lLine))
    {
        std::size_t lDelimiterPos = lLine.find(' ');

        lRanges.push_back({
            .min = mpz_class(lLine.substr(0, lDelimiterPos)),
            .max = mpz_class(lLine.substr(lDelimiterPos))
        });
    }

    return lRanges;
}

std::vector<Range> collapseRangesOverlaps(std::vector<Range> iRanges)
{
    std::vector<Range> lMergedRanges = std::move(iRanges);

    bool lHasOverlaps = true;
    while (lHasOverlaps)
    {
        lHasOverlaps = false;
        std::vector<Range> lTestedRanges;

        for (const auto &lBaseRange : lMergedRanges)
        {
            bool lCollapsed = false;
            for (auto &lTestedRange : lTestedRanges)
            {
                // base range in tested range, ignoring it
                if (lBaseRange.min >= lTestedRange.min && lBaseRange.max <= lTestedRange.max)
                    lCollapsed = true;

                if (lBaseRange.min <= lTestedRange.min && lBaseRange.max >= lTestedRange.min) {
                    lTestedRange.min = lBaseRange.min;
                    lHasOverlaps = true;
                    lCollapsed = true;
                }

                if (lBaseRange.max >= lTestedRange.max && lBaseRange.min <= lTestedRange.max) {
                    lTestedRange.max = lBaseRange.max;
                    lHasOverlaps = true;
                    lCollapsed = true;
                }

                if (lCollapsed) break;
            }

            if (!lCollapsed)
                lTestedRanges.push_back(lBaseRange);
        }

        lMergedRanges = std::move(lTestedRanges); // todo: vérifier le fonctionnement du std::move
        lTestedRanges.clear();
    }

    return lMergedRanges;
}

void testExecutionTimes(int iThreadAmount, const std::string& iFilename)
{
    auto lChrono = new Chrono();

    std::ofstream lExecTimesFile;
    lExecTimesFile.open("exec_times.txt", std::ios_base::app);
    lExecTimesFile << iFilename << " with " << iThreadAmount << " thread(s)" << std::endl;

    for (int i = 1; i <= iThreadAmount; i++)
    {
        lChrono->resume();

        gRangeCounter = 0;
        gInRangeCounter = 0;

        std::cout << "Working with " << i << " thread(s)..." << std::endl;

        pthread_t lThreads[i];

        for (int j = 0; j < i; j++)
        {
            pthread_create(&lThreads[j], NULL, computePrimes, NULL);
        }

        std::vector<mpz_class> lFoundPrimes;

        for (int j = 0; j < i; j++)
        {
            void* lThreadRetVal = nullptr;
            pthread_join(lThreads[j], &lThreadRetVal);

            std::vector<mpz_class> lThreadFinds = *(std::vector<mpz_class>*)lThreadRetVal;
            lFoundPrimes.insert(lFoundPrimes.end(), lThreadFinds.begin(), lThreadFinds.end());
        }

        lChrono->pause();

        std::cout << "It took " << lChrono->get() << "s to complete." << std::endl;


        lExecTimesFile << i << ";" << lChrono->get() << std::endl;


        lChrono->reset();
    }

    lExecTimesFile.close();
}

bool isPrime(const mpz_class& n)
{
    mpz_class i;
    bool isPrime = true;

    // 0 and 1 are not prime numbers
    if (n == 0 || n == 1) {
        isPrime = false;
    }
    else {
        for (i = 2; i <= n / 2; ++i) {
            if (n % i == 0) {
                isPrime = false;
                break;
            }
        }
    }
    return isPrime;
}

int main(int argc, char **argv) {

    if (argc < 3) {
        std::cout << "Please enter thread amount and filename !" << std::endl;
        return 0;
    }

    int lThreadAmount = atoi(argv[1]);
    std::string lFilename = argv[2];

    if (lThreadAmount < 1) {
        std::cout << "Thread amount must be greater than 0 !" << std::endl;
        return 0;
    }

    if (!std::filesystem::exists(lFilename)) {
        std::cout << "File with given filename does not exists !" << std::endl;
        return 0;
    }

    std::cout << "Reading ranges from given file..." << std::endl;
    std::vector<Range> lFileRanges = getRangesFromFile(lFilename);
    std::cout << "OK!" << std::endl;

    std::cout << "Collapsing overlapping ranges..." << std::endl;
    gRanges = collapseRangesOverlaps(lFileRanges);
    std::cout << "OK!" << std::endl;

    std::cout << "Beginning prime number check..." << std::endl;

    if (argc > 3 && argv[3] == std::string("test"))
    {
        testExecutionTimes(lThreadAmount, lFilename);
        return 0;
    }

    auto lChrono = new Chrono();

    pthread_t lThreads[lThreadAmount];

    for (int i = 0; i < lThreadAmount; i++)
    {
        pthread_create(&lThreads[i], NULL, computePrimes, NULL);
    }

    std::vector<mpz_class> oAllThreadFinds;

    for (int i = 0; i < lThreadAmount; i++)
    {
        void* lThreadRetVal = nullptr;
        pthread_join(lThreads[i], &lThreadRetVal);

        std::vector<mpz_class> lThFoundPrimes = *(std::vector<mpz_class>*)lThreadRetVal;
        oAllThreadFinds.insert(oAllThreadFinds.end(), lThFoundPrimes.begin(), lThFoundPrimes.end());
    }

    lChrono->pause();

    std::ofstream myfile;
    myfile.open ("temps.txt");
    myfile << oAllThreadFinds.size() << " prime numbers found in " << lChrono->get() << "s" << std::endl;
    myfile.close();

    std::cerr << oAllThreadFinds.size() << " prime numbers found in " << lChrono->get() << "s" << std::endl;

    std::sort(oAllThreadFinds.begin(), oAllThreadFinds.end());
    std::cout << "Found primes are:" << std::endl;
    for (const auto &item : oAllThreadFinds)
    {
        /*if (!isPrime(item))
        {
            std::cout << std::endl;
            std::cout << item << " is not prime!";
            std::cout << std::endl;
            return 0;
        }*/
        std::cout << item << ", ";
    }

    std::cout << std::endl;
    return 0;
}