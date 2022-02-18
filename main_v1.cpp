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
pthread_mutex_t gSetMutex = PTHREAD_MUTEX_INITIALIZER;

int gRangeCounter = 0;
int gInRangeCounter = 0;
std::vector<Range> gRanges;
std::vector<mpz_class> gFoundPrimes;

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

void addPrimeNumber(const mpz_class& iPrime)
{
    pthread_mutex_lock(&gSetMutex);
    gFoundPrimes.push_back(iPrime);
    pthread_mutex_unlock(&gSetMutex);
}

void* computePrimes(void* arg)
{
    auto lPrng = generate_prng();

    mpz_class lEvaluatedNumber;
    while (getNextNumber(lEvaluatedNumber))
        if (prob_prime(lEvaluatedNumber, 3, lPrng))
            addPrimeNumber(lEvaluatedNumber);

    pthread_exit(NULL);
    return NULL;
}

std::vector<Range> getRangesFromFile(const std::string& filename)
{
    std::vector<Range> ranges;

    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line))
    {
        std::size_t delimiterPos = line.find(' ');

        ranges.push_back({
            .min = mpz_class(line.substr(0, delimiterPos)),
            .max = mpz_class(line.substr(delimiterPos))
        });
    }

    return ranges;
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

    if (!std::filesystem::exists("../files/" + lFilename)) {
        std::cout << "File with given filename does not exists !" << std::endl;
        return 0;
    }

    std::cout << "Reading ranges from given file..." << std::endl;
    std::vector<Range> lFileRanges = getRangesFromFile("../files/" + lFilename);
    std::cout << "OK!" << std::endl;

    std::cout << "Collapsing overlapping ranges..." << std::endl;
    gRanges = collapseRangesOverlaps(lFileRanges);
    std::cout << "OK!" << std::endl;


    /*std::vector<Range> fileRanges = getRangesFromFile("../files/" + lFilename);
    for (const auto &item : fileRanges)
    {
        std::cout << "min: " << item.min << " max: " << item.max << std::endl;
    }
    fileRanges = collapseRangesOverlaps(fileRanges);
    std::cout << "---" << std::endl;
    for (const auto &item : fileRanges)
    {
        std::cout << "min: " << item.min << " max: " << item.max << std::endl;
    }
    return 0;*/

    std::cout << "Beginning prime number check..." << std::endl;

    auto lChrono = new Chrono();

    pthread_t lThreads[lThreadAmount];

    for (int i = 0; i < lThreadAmount; i++)
    {
        pthread_create(&lThreads[i], NULL, computePrimes, NULL);
    }

    for (int i = 0; i < lThreadAmount; i++)
    {
        pthread_join(lThreads[i], NULL);
    }

    lChrono->pause();

    std::ofstream myfile;
    myfile.open ("temps.txt");
    myfile << gFoundPrimes.size() << " prime numbers found in " << lChrono->get() << "s" << std::endl;
    myfile.close();

    std::cerr << gFoundPrimes.size() << " prime numbers found in " << lChrono->get() << "s" << std::endl;

    std::sort(gFoundPrimes.begin(), gFoundPrimes.end());
    std::cout << "Found primes are:" << std::endl;
    for (const auto &item : gFoundPrimes)
    {
        std::cout << item << ", ";
    }

    return 0;
}