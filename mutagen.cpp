#include <iostream>
#include <array>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <thread>
#include <atomic>
#include <chrono>
#include <queue>
#include <mutex>
#include <cstring>
#include <unordered_map>
#include <cmath>
#include <immintrin.h>
#include <omp.h>
#include <csignal>
#include <random>
#include <algorithm>
#include <getopt.h>

#ifdef _WIN32
#include <windows.h>
#endif

// Include the required headers
#include "sha256_avx2.h"
#include "ripemd160_avx2.h"
#include "SECP256K1.h"
#include "Point.h"
#include "Int.h"
#include "IntGroup.h"

using namespace std;

// Cross-platform terminal functions
void initConsole()
{
#ifdef _WIN32
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD mode = 0;
    GetConsoleMode(hConsole, &mode);
    SetConsoleMode(hConsole, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
#endif
}

void clearTerminal()
{
#ifdef _WIN32
    HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
    COORD coord = {0, 0};
    DWORD count;
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(hStdOut, &csbi);
    FillConsoleOutputCharacter(hStdOut, ' ', csbi.dwSize.X * csbi.dwSize.Y, coord, &count);
    SetConsoleCursorPosition(hStdOut, coord);
#else
    std::cout << "\033[2J\033[H";
#endif
    std::cout.flush();
}

void moveCursorTo(int x, int y)
{
#ifdef _WIN32
    HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
    COORD coord = {(SHORT)x, (SHORT)y};
    SetConsoleCursorPosition(hStdOut, coord);
#else
    std::cout << "\033[" << y << ";" << x << "H";
#endif
    std::cout.flush();
}

// Configuration defaults
int PUZZLE_NUM = 20;
int WORKERS = omp_get_num_procs();
int FLIP_COUNT = -1;
static constexpr int POINTS_BATCH_SIZE = 256;
static constexpr int HASH_BATCH_SIZE = 8;

// Historical puzzle data (puzzle number: {flip count, target hash, private key decimal})
const unordered_map<int, tuple<int, string, string>> PUZZLE_DATA = {
    {20, {8, "b907c3a2a3b27789dfb509b730dd47703c272868", "357535"}},
    {21, {9, "29a78213caa9eea824acf08022ab9dfc83414f56", "863317"}},
    {22, {11, "7ff45303774ef7a52fffd8011981034b258cb86b", "1811764"}},
    {23, {12, "d0a79df189fe1ad5c306cc70497b358415da579e", "3007503"}},
    {24, {9, "0959e80121f36aea13b3bad361c15dac26189e2f", "5598802"}},
    {25, {12, "2f396b29b27324300d0c59b17c3abc1835bd3dbb", "14428676"}},
    {26, {14, "bfebb73562d4541b32a02ba664d140b5a574792f", "33185509"}},
    {27, {13, "0c7aaf6caa7e5424b63d317f0f8f1f9fa40d5560", "54538862"}},
    {28, {16, "1306b9e4ff56513a476841bac7ba48d69516b1da", "111949941"}},
    {29, {18, "5a416cc9148f4a377b672c8ae5d3287adaafadec", "227634408"}},
    {30, {16, "d39c4704664e1deb76c9331e637564c257d68a08", "400708894"}},
    {31, {13, "d805f6f251f7479ebd853b3d0f4b9b2656d92f1d", "1033162084"}},
    {32, {14, "9e42601eeaedc244e15f17375adb0e2cd08efdc9", "2102388551"}},
    {33, {15, "4e15e5189752d1eaf444dfd6bff399feb0443977", "3093472814"}},
    {34, {16, "f6d67d7983bf70450f295c9cb828daab265f1bfa", "7137437912"}},
    {35, {19, "f6d8ce225ffbdecec170f8298c3fc28ae686df25", "14133072157"}},
    {36, {14, "74b1e012be1521e5d8d75e745a26ced845ea3d37", "20112871792"}},
    {37, {23, "28c30fb11ed1da72e7c4f89c0164756e8a021d", "42387769980"}},
    {38, {21, "b190e2d40cfdeee2cee072954a2be89e7ba39364", "100251560595"}},
    {39, {23, "0b304f2a79a027270276533fe1ed4eff30910876", "146971536592"}},
    {40, {20, "95a156cd21b4a69de969eb6716864f4c8b82a82a", "323724968937"}},
    {41, {25, "d1562eb37357f9e6fc41cb2359f4d3eda4032329", "1003651412950"}},
    {42, {24, "8efb85f9c5b5db2d55973a04128dc7510075ae23", "1458252205147"}},
    {43, {19, "f92044c7924e5525c61207972c253c9fc9f086f7", "2895374552463"}},
    {44, {24, "80df54e1f612f2fc5bdc05c9d21a83aa8d20791e", "7409811047825"}},
    {45, {21, "f0225bfc68a6e17e87cd8b5e60ae3be18f120753", "15404761757071"}},
    {46, {24, "9a012260d01c5113df66c8a8438c9f7a1e3d5dac", "19996463086597"}},
    {47, {27, "f828005d41b0f4fed4c8dca3b06011072cfb07d4", "51408670348612"}},
    {48, {21, "8661cb56d9df0a61f01328b55af7e56a3fe7a2b2", "119666659114170"}},
    {49, {30, "0d2f533966c6578e1111978ca698f8add7fffdf3", "191206974700443"}},
    {50, {29, "de081b76f840e462fa2cdf360173dfaf4a976a47", "409118905032525"}},
    {51, {25, "ef6419cffd7fad7027994354eb8efae223c2dbe7", "611140496167764"}},
    {52, {27, "36af659edbe94453f6344e920d143f1778653ae7", "2058769515153876"}},
    {53, {26, "2f4870ef54fa4b048c1365d42594cc7d3d269551", "4216495639600700"}},
    {54, {30, "cb66763cf7fde659869ae7f06884d9a0f879a092", "6763683971478124"}},
    {55, {31, "db53d9bbd1f3a83b094eeca7dd970bd85b492fa2", "9974455244496707"}},
    {56, {31, "48214c5969ae9f43f75070cea1e2cb41d5bdcccd", "30045390491869460"}},
    {57, {33, "328660ef43f66abe2653fa178452a5dfc594c2a1", "44218742292676575"}},
    {58, {28, "8c2a6071f89c90c4dab5ab295d7729d1b54ea60f", "138245758910846492"}},
    {59, {30, "b14ed3146f5b2c9bde1703deae9ef33af8110210", "199976667976342049"}},
    {60, {31, "cdf8e5c7503a9d22642e3ecfc87817672787b9c5", "525070384258266191"}},
    {61, {25, "68133e19b2dfb9034edf9830a200cfdf38c90cbd", "1135041350219496382"}},
    {62, {35, "e26646db84b0602f32b34b5a62ca3cae1f91b779", "1425787542618654982"}},
    {63, {34, "ef58afb697b094423ce90721fbb19a359ef7c50e", "3908372542507822062"}},
    {64, {34, "3ee4133d991f52fdf6a25c9834e0745ac74248a4", "8993229949524469768"}},
    {65, {37, "52e763a7ddc1aa4fa811578c491c1bc7fd570137", "17799667357578236628"}},
    {66, {35, "20d45a6a762535700ce9e0b216e31994335db8a5", "30568377312064202855"}},
    {67, {31, "739437bb3dd6d1983e66629c5f08c70e52769371", "46346217550346335726"}},
    {68, {34, "e0b8a2baee1b77fc703455f39d51477451fc8cfc", "132656943602386256302"}}};

// Global variables
alignas(32) vector<unsigned char> TARGET_HASH160_RAW(20);
string TARGET_HASH160;
Int BASE_KEY;
atomic<bool> stop_event(false);
mutex result_mutex;
queue<tuple<string, __uint128_t, int>> results;
atomic<int> completed_threads(0);
int num_threads = WORKERS;

// AVXCounter implementation
union AVXCounter
{
    __m256i vec;
    uint64_t u64[4];
    __uint128_t u128[2];

    AVXCounter() : vec(_mm256_setzero_si256()) {}
    AVXCounter(__uint128_t value) { store(value); }

    void increment()
    {
        __m256i one = _mm256_set_epi64x(0, 0, 0, 1);
        vec = _mm256_add_epi64(vec, one);
        if (u64[0] == 0)
        {
            __m256i carry = _mm256_set_epi64x(0, 0, 1, 0);
            vec = _mm256_add_epi64(vec, carry);
        }
    }

    void add(__uint128_t value)
    {
        __m256i add_val = _mm256_set_epi64x(0, 0, value >> 64, value);
        vec = _mm256_add_epi64(vec, add_val);
        if (u64[0] < (value & 0xFFFFFFFFFFFFFFFFULL))
        {
            __m256i carry = _mm256_set_epi64x(0, 0, 1, 0);
            vec = _mm256_add_epi64(vec, carry);
        }
    }

    __uint128_t load() const
    {
        return (static_cast<__uint128_t>(u64[1]) << 64) | u64[0];
    }

    void store(__uint128_t value)
    {
        u64[0] = static_cast<uint64_t>(value);
        u64[1] = static_cast<uint64_t>(value >> 64);
        u64[2] = 0;
        u64[3] = 0;
    }

    static AVXCounter div(const AVXCounter &num, uint64_t denom)
    {
        __uint128_t n = num.load();
        __uint128_t q = n / denom;
        return AVXCounter(q);
    }

    static uint64_t mod(const AVXCounter &num, uint64_t denom)
    {
        __uint128_t n = num.load();
        return n % denom;
    }

    static AVXCounter mul(uint64_t a, uint64_t b)
    {
        __uint128_t result = static_cast<__uint128_t>(a) * b;
        return AVXCounter(result);
    }
};

static AVXCounter total_checked_avx;
__uint128_t total_combinations = 0;
vector<string> g_threadPrivateKeys;
mutex progress_mutex;

// Performance tracking
atomic<uint64_t> globalComparedCount(0);
atomic<uint64_t> localComparedCount(0);
double globalElapsedTime = 0.0;
double mkeysPerSec = 0.0;
chrono::time_point<chrono::high_resolution_clock> tStart;

// Helper functions
string formatElapsedTime(double seconds)
{
    int hrs = static_cast<int>(seconds) / 3600;
    int mins = (static_cast<int>(seconds) % 3600) / 60;
    int secs = static_cast<int>(seconds) % 60;
    ostringstream oss;
    oss << setw(2) << setfill('0') << hrs << ":"
        << setw(2) << setfill('0') << mins << ":"
        << setw(2) << setfill('0') << secs;
    return oss.str();
}

string to_string_128(__uint128_t value)
{
    if (value == 0)
        return "0";
    char buffer[50];
    char *p = buffer + sizeof(buffer);
    *--p = '\0';
    while (value != 0)
    {
        *--p = "0123456789"[value % 10];
        value /= 10;
    }
    return string(p);
}

void signalHandler(int signum)
{
    stop_event.store(true);
    cout << "\nInterrupt received, shutting down...\n";
}

class FastCombinationGenerator
{
    const int n;                       // Total bits (e.g., 20 for puzzle 20)
    const int k;                       // Number of bits to flip (e.g., 8)
    __uint128_t next_rank;             // Current combination rank
    std::vector<int> last_combination; // Last generated combination

public:
    FastCombinationGenerator(int n, int k)
        : n(n), k(k), next_rank(0), last_combination(k, -1) {}

    // Generate a batch of 512 combinations
    std::vector<int> generate512Batch()
    {
        std::vector<int> batch;
        batch.reserve(512 * k);

        // check for completion of all combinations
        if (next_rank >= combinations_count(n, k))
        {
            return batch;
        }

        for (int i = 0; i < 512 && next_rank < combinations_count(n, k); ++i)
        {
            generate_next(batch);
        }
        return batch;
    }

    // Get current combination rank
    __uint128_t currentRank() const
    {
        return next_rank;
    }

    // Calculate total possible combinations
    static __uint128_t combinations_count(int n, int k)
    {
        if (k < 0 || k > n)
            return 0;

        __uint128_t res = 1;
        for (int i = 1; i <= k; ++i)
        {
            res = res * (n - k + i) / i;
        }
        return res;
    }

private:
    // Generate next combination and append to batch
    void generate_next(std::vector<int> &batch)
    {
        if (next_rank == 0)
        {
            // First combination: 0,1,2,...,k-1
            for (int i = 0; i < k; ++i)
            {
                last_combination[i] = i;
                batch.push_back(i);
            }
        }
        else
        {
            // Generate next combination in lexicographic order
            int i = k - 1;
            while (i >= 0 && last_combination[i] == n - k + i)
            {
                --i;
            }

            if (i < 0)
                return; // All combinations generated

            ++last_combination[i];
            for (int j = i + 1; j < k; ++j)
            {
                last_combination[j] = last_combination[j - 1] + 1;
            }

            // Append to batch
            batch.insert(batch.end(), last_combination.begin(), last_combination.end());
        }

        ++next_rank;
    }

    // Helper function for unranking (kept for compatibility)
    void unrank(__uint128_t rank, int *output)
    {
        __uint128_t x = combinations_count(n, k) - 1 - rank;
        int a = n, b = k;

        for (int i = 0; i < k; ++i)
        {
            a = largest_a(a, b, x);
            output[i] = (n - 1) - a;
            x -= combinations_count(a, b);
            b--;
        }
    }

    int largest_a(int a, int b, __uint128_t x) const
    {
        while (a >= b && combinations_count(a, b) > x)
        {
            a--;
        }
        return a;
    }
};

void prepareShaBlock(const uint8_t *dataSrc, __uint128_t dataLen, uint8_t *outBlock)
{
    fill_n(outBlock, 64, 0);
    memcpy(outBlock, dataSrc, dataLen);
    outBlock[dataLen] = 0x80;
    uint32_t bitLen = (uint32_t)(dataLen * 8);
    outBlock[60] = (uint8_t)((bitLen >> 24) & 0xFF);
    outBlock[61] = (uint8_t)((bitLen >> 16) & 0xFF);
    outBlock[62] = (uint8_t)((bitLen >> 8) & 0xFF);
    outBlock[63] = (uint8_t)(bitLen & 0xFF);
}

void prepareRipemdBlock(const uint8_t *dataSrc, uint8_t *outBlock)
{
    fill_n(outBlock, 64, 0);
    memcpy(outBlock, dataSrc, 32);
    outBlock[32] = 0x80;
    uint32_t bitLen = 256;
    outBlock[60] = (uint8_t)((bitLen >> 24) & 0xFF);
    outBlock[61] = (uint8_t)((bitLen >> 16) & 0xFF);
    outBlock[62] = (uint8_t)((bitLen >> 8) & 0xFF);
    outBlock[63] = (uint8_t)(bitLen & 0xFF);
}

void computeHash160BatchBinSingle(int numKeys, uint8_t pubKeys[][33], uint8_t hashResults[][20])
{
    alignas(32) array<array<uint8_t, 64>, HASH_BATCH_SIZE> shaInputs;
    alignas(32) array<array<uint8_t, 32>, HASH_BATCH_SIZE> shaOutputs;
    alignas(32) array<array<uint8_t, 64>, HASH_BATCH_SIZE> ripemdInputs;
    alignas(32) array<array<uint8_t, 20>, HASH_BATCH_SIZE> ripemdOutputs;

    __uint128_t totalBatches = (numKeys + (HASH_BATCH_SIZE - 1)) / HASH_BATCH_SIZE;

    for (__uint128_t batch = 0; batch < totalBatches; batch++)
    {
        __uint128_t batchCount = min<__uint128_t>(HASH_BATCH_SIZE, numKeys - batch * HASH_BATCH_SIZE);

        // Prepare SHA-256 input blocks
        for (__uint128_t i = 0; i < batchCount; i++)
        {
            prepareShaBlock(pubKeys[batch * HASH_BATCH_SIZE + i], 33, shaInputs[i].data());
        }

        if (batchCount < HASH_BATCH_SIZE)
        {
            static array<uint8_t, 64> shaPadding = {};
            prepareShaBlock(pubKeys[0], 33, shaPadding.data());
            for (size_t i = batchCount; i < HASH_BATCH_SIZE; i++)
            {
                memcpy(shaInputs[i].data(), shaPadding.data(), 64);
            }
        }

        const uint8_t *inPtr[HASH_BATCH_SIZE];
        uint8_t *outPtr[HASH_BATCH_SIZE];
        for (int i = 0; i < HASH_BATCH_SIZE; i++)
        {
            inPtr[i] = shaInputs[i].data();
            outPtr[i] = shaOutputs[i].data();
        }

        sha256avx2_8B(inPtr[0], inPtr[1], inPtr[2], inPtr[3],
                      inPtr[4], inPtr[5], inPtr[6], inPtr[7],
                      outPtr[0], outPtr[1], outPtr[2], outPtr[3],
                      outPtr[4], outPtr[5], outPtr[6], outPtr[7]);

        for (__uint128_t i = 0; i < batchCount; i++)
        {
            prepareRipemdBlock(shaOutputs[i].data(), ripemdInputs[i].data());
        }

        if (batchCount < HASH_BATCH_SIZE)
        {
            static array<uint8_t, 64> ripemdPadding = {};
            prepareRipemdBlock(shaOutputs[0].data(), ripemdPadding.data());
            for (size_t i = batchCount; i < HASH_BATCH_SIZE; i++)
            {
                memcpy(ripemdInputs[i].data(), ripemdPadding.data(), 64);
            }
        }

        for (int i = 0; i < HASH_BATCH_SIZE; i++)
        {
            inPtr[i] = ripemdInputs[i].data();
            outPtr[i] = ripemdOutputs[i].data();
        }

        ripemd160avx2::ripemd160avx2_32(
            (unsigned char *)inPtr[0], (unsigned char *)inPtr[1],
            (unsigned char *)inPtr[2], (unsigned char *)inPtr[3],
            (unsigned char *)inPtr[4], (unsigned char *)inPtr[5],
            (unsigned char *)inPtr[6], (unsigned char *)inPtr[7],
            outPtr[0], outPtr[1], outPtr[2], outPtr[3],
            outPtr[4], outPtr[5], outPtr[6], outPtr[7]);

        for (__uint128_t i = 0; i < batchCount; i++)
        {
            memcpy(hashResults[batch * HASH_BATCH_SIZE + i], ripemdOutputs[i].data(), 20);
        }
    }
}

void pointToCompressedBin(Point &p, uint8_t out[33])
{
    out[0] = p.y.IsEven() ? 0x02 : 0x03;
    p.x.Get32Bytes(out + 1);
}

void worker(Secp256K1 *secp, int bit_length, int flip_count, int threadId,
            __uint128_t start_rank, __uint128_t end_rank)
{

    if (start_rank >= end_rank)
    {
        completed_threads++;
        return;
    }
    FastCombinationGenerator gen(bit_length, flip_count);

    // Precomputed points
    alignas(32) Point plusPoints[POINTS_BATCH_SIZE];
    alignas(32) Point minusPoints[POINTS_BATCH_SIZE];
    for (int i = 0; i < POINTS_BATCH_SIZE; i++)
    {
        Int tmp;
        tmp.SetInt32(i);
        plusPoints[i] = secp->ComputePublicKey(&tmp);
        minusPoints[i] = plusPoints[i];
        minusPoints[i].y.ModNeg();
    }

    const int BATCH_SIZE = 512;
    alignas(32) Int privKeys[BATCH_SIZE];
    alignas(32) Point pubKeys[BATCH_SIZE];
    alignas(32) uint8_t hashBatch[BATCH_SIZE][20];
    alignas(32) uint8_t pubKeyBatch[HASH_BATCH_SIZE][33];

    while (!stop_event.load() && gen.currentRank() < end_rank)
    {
        auto masks = gen.generate512Batch();
        if (masks.empty())
        {
            break;
        }

        int actual_batch_size = masks.size() / flip_count;
        if (actual_batch_size == 0)
        {
            break;
        }

        for (int i = 0; i < actual_batch_size; i++)
        {
            privKeys[i].Set(&BASE_KEY);
            for (int j = 0; j < flip_count; j++)
            {
                Int bit;
                bit.SetInt32(1);
                bit.ShiftL(masks[i * flip_count + j]);
                privKeys[i].Xor(&bit);
            }

            pubKeys[i] = secp->ComputePublicKey(&privKeys[i]);

            // std::cout << "Private key " << i << ": " << privKeys[i].GetBase16() << std::endl;
            // std::cout << "Public key:  " << pubKeys[i].x.GetBase16() << std::endl;
        }

        for (int batch_offset = 0; batch_offset < actual_batch_size; batch_offset++)
        {
            Point startPoint = pubKeys[batch_offset];
            Int startPointX, startPointY, startPointXNeg;
            startPointX.Set(&startPoint.x);
            startPointY.Set(&startPoint.y);
            startPointXNeg.Set(&startPointX);
            startPointXNeg.ModNeg();

            alignas(32) Int deltaX[POINTS_BATCH_SIZE];
            for (int i = 0; i < POINTS_BATCH_SIZE; i++)
            {
                deltaX[i].ModSub(&plusPoints[i].x, &startPointX);
            }

            IntGroup modGroup(POINTS_BATCH_SIZE);
            modGroup.Set(deltaX);
            modGroup.ModInv();

            // 10. Generate derived keys
            alignas(32) Int pointBatchX[2 * POINTS_BATCH_SIZE];
            alignas(32) Int pointBatchY[2 * POINTS_BATCH_SIZE];

            // Plus points
            for (int i = 0; i < POINTS_BATCH_SIZE; i++)
            {
                Int deltaY;
                deltaY.ModSub(&plusPoints[i].y, &startPointY);
                Int slope;
                slope.ModMulK1(&deltaY, &deltaX[i]);
                Int slopeSq;
                slopeSq.ModSquareK1(&slope);

                pointBatchX[i].Set(&startPointXNeg);
                pointBatchX[i].ModAdd(&slopeSq);
                pointBatchX[i].ModSub(&plusPoints[i].x);

                Int diffX;
                diffX.ModSub(&startPointX, &pointBatchX[i]);
                diffX.ModMulK1(&slope);

                pointBatchY[i].Set(&startPointY);
                pointBatchY[i].ModNeg();
                pointBatchY[i].ModAdd(&diffX);
            }

            // Minus points
            for (int i = 0; i < POINTS_BATCH_SIZE; i++)
            {
                Int deltaY;
                deltaY.ModSub(&minusPoints[i].y, &startPointY);
                Int slope;
                slope.ModMulK1(&deltaY, &deltaX[i]);
                Int slopeSq;
                slopeSq.ModSquareK1(&slope);

                pointBatchX[POINTS_BATCH_SIZE + i].Set(&startPointXNeg);
                pointBatchX[POINTS_BATCH_SIZE + i].ModAdd(&slopeSq);
                pointBatchX[POINTS_BATCH_SIZE + i].ModSub(&minusPoints[i].x);

                Int diffX;
                diffX.ModSub(&startPointX, &pointBatchX[POINTS_BATCH_SIZE + i]);
                diffX.ModMulK1(&slope);

                pointBatchY[POINTS_BATCH_SIZE + i].Set(&startPointY);
                pointBatchY[POINTS_BATCH_SIZE + i].ModNeg();
                pointBatchY[POINTS_BATCH_SIZE + i].ModAdd(&diffX);
            }

            // Hash the batch
            int keys_processed = 0;
            while (keys_processed < 2 * POINTS_BATCH_SIZE)
            {
                int batch_size = min(HASH_BATCH_SIZE, 2 * POINTS_BATCH_SIZE - keys_processed);

                for (int i = 0; i < batch_size; i++)
                {
                    Point p;
                    p.x.Set(&pointBatchX[keys_processed + i]);
                    p.y.Set(&pointBatchY[keys_processed + i]);
                    pointToCompressedBin(p, pubKeyBatch[i]);
                }

                computeHash160BatchBinSingle(batch_size, pubKeyBatch, &hashBatch[keys_processed]);

                for (int i = 0; i < batch_size; i++)
                {
                    printf("Checking key %d:\n", keys_processed + i);

                    printf("  Compressed pubkey: ");
                    for (int j = 0; j < 33; j++)
                    {
                        printf("%02x", pubKeyBatch[i][j]);
                    }
                    printf("\n");

                    printf("  Computed hash160: ");
                    for (int j = 0; j < 20; j++)
                    {
                        printf("%02x", hashBatch[keys_processed + i][j]);
                    }
                    printf("\n");

                    printf("  Target hash160:   ");
                    for (int j = 0; j < 20; j++)
                    {
                        printf("%02x", TARGET_HASH160_RAW.data()[j]);
                    }
                    printf("\n");

                    bool match = true;
                    printf("  Byte comparison:  ");
                    for (int j = 0; j < 20; j++)
                    {
                        if (hashBatch[keys_processed + i][j] != TARGET_HASH160_RAW.data()[j])
                        {
                            printf("XX ");
                            match = false;
                        }
                        else
                        {
                            printf("%02x ", hashBatch[keys_processed + i][j]);
                        }
                    }
                    printf("\n");

                    if (match)
                    {
                        printf("!!! EXACT MATCH FOUND !!!\n");

                        Int foundKey;
                        foundKey.Set(&privKeys[batch_offset]);
                        if (keys_processed + i < POINTS_BATCH_SIZE)
                        {
                            printf("Adding offset: %d\n", keys_processed + i);
                            foundKey.Add(keys_processed + i);
                        }
                        else
                        {
                            printf("Subtracting offset: %d\n", keys_processed + i - POINTS_BATCH_SIZE);
                            foundKey.Sub(keys_processed + i - POINTS_BATCH_SIZE);
                        }

                        string hexKey = foundKey.GetBase16();
                        hexKey = string(64 - hexKey.length(), '0') + hexKey;
                        printf("Found private key: %s\n", hexKey.c_str());

                        lock_guard<mutex> lock(result_mutex);
                        results.push(make_tuple(hexKey, total_checked_avx.load(), flip_count));
                        stop_event.store(true);
                        return;
                    }
                    else
                    {
                        printf("  No match\n\n");
                    }
                }
                keys_processed += batch_size;
            }
        }
        total_checked_avx.add(actual_batch_size * 2 * POINTS_BATCH_SIZE);

        if (threadId == 0)
        {
            static auto last_update = chrono::high_resolution_clock::now();
            auto now = chrono::high_resolution_clock::now();

            if (chrono::duration<double>(now - last_update).count() >= 5.0) // 5 seconds
            {
                last_update = now;
                globalElapsedTime = chrono::duration<double>(now - tStart).count();
                mkeysPerSec = (double)total_checked_avx.load() / globalElapsedTime / 1e6;

                double progress = min(100.0, (double)total_checked_avx.load() / total_combinations * 100.0);

                lock_guard<mutex> lock(progress_mutex);
                moveCursorTo(0, 10);
                cout << "Progress: " << fixed << setprecision(6) << progress << "%\n";
                cout << "Processed: " << to_string_128(total_checked_avx.load()) << "\n";
                cout << "Speed: " << fixed << setprecision(2) << mkeysPerSec << " Mkeys/s\n";
                cout << "Elapsed Time: " << formatElapsedTime(globalElapsedTime) << "\n";
                cout.flush();
            }
        }

        if (gen.currentRank() >= end_rank)
        {
            completed_threads++;
            if (completed_threads >= WORKERS || total_checked_avx.load() >= total_combinations)
            {
                stop_event.store(true);
            }
        }
    }
}

void printUsage(const char *programName)
{
    cout << "Usage: " << programName << " [options]\n";
    cout << "Options:\n";
    cout << "  -p, --puzzle NUM    Puzzle number to solve (default: 20)\n";
    cout << "  -t, --threads NUM   Number of CPU cores to use (default: all)\n";
    cout << "  -f, --flips NUM     Override default flip count for puzzle\n";
    cout << "  -h, --help          Show this help message\n";
    cout << "\nExample:\n";
    cout << "  " << programName << " -p 38 -t 8 -f 21\n";
}

int main(int argc, char *argv[])
{
    signal(SIGINT, signalHandler);
    initConsole();

    // Parse command line arguments
    int opt;
    int option_index = 0;
    static struct option long_options[] = {
        {"puzzle", required_argument, 0, 'p'},
        {"threads", required_argument, 0, 't'},
        {"flips", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "p:t:f:h", long_options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 'p':
            PUZZLE_NUM = atoi(optarg);
            if (PUZZLE_NUM < 20 || PUZZLE_NUM > 68)
            {
                cerr << "Error: Puzzle number must be between 20 and 68\n";
                return 1;
            }
            break;
        case 't':
            WORKERS = atoi(optarg);
            if (WORKERS < 1)
            {
                cerr << "Error: Thread count must be at least 1\n";
                return 1;
            }
            break;
        case 'f':
            FLIP_COUNT = atoi(optarg);
            if (FLIP_COUNT < 1)
            {
                cerr << "Error: Flip count must be at least 1\n";
                return 1;
            }
            break;
        case 'h':
            printUsage(argv[0]);
            return 0;
        default:
            printUsage(argv[0]);
            return 1;
        }
    }

    tStart = chrono::high_resolution_clock::now();

    Secp256K1 secp;
    secp.Init();

    auto puzzle_it = PUZZLE_DATA.find(PUZZLE_NUM);
    if (puzzle_it == PUZZLE_DATA.end())
    {
        cerr << "Error: Invalid puzzle number\n";
        return 1;
    }

    auto [DEFAULT_FLIP_COUNT, TARGET_HASH160_HEX, PRIVATE_KEY_DECIMAL] = puzzle_it->second;
    if (FLIP_COUNT == -1)
    {
        FLIP_COUNT = DEFAULT_FLIP_COUNT;
    }

    TARGET_HASH160 = TARGET_HASH160_HEX;
    for (int i = 0; i < 20; i++)
    {
        TARGET_HASH160_RAW[i] = stoul(TARGET_HASH160.substr(i * 2, 2), nullptr, 16);
    }

    BASE_KEY.SetBase10(const_cast<char *>(PRIVATE_KEY_DECIMAL.c_str()));

    // Verify base key
    Int testKey;
    testKey.SetBase10(const_cast<char *>(PRIVATE_KEY_DECIMAL.c_str()));
    if (!testKey.IsEqual(&BASE_KEY))
    {
        cerr << "Base key initialization failed!\n";
        return 1;
    }

    if (BASE_KEY.GetBitLength() > PUZZLE_NUM)
    {
        cerr << "Base key exceeds puzzle bit length!\n";
        return 1;
    }

    // Calculate total combinations
    total_combinations = FastCombinationGenerator::combinations_count(PUZZLE_NUM, FLIP_COUNT);
    num_threads = WORKERS;

    // Format base key for display
    string paddedKey = BASE_KEY.GetBase16();
    size_t firstNonZero = paddedKey.find_first_not_of('0');
    paddedKey = paddedKey.substr(firstNonZero);
    paddedKey = "0x" + paddedKey;

    // Print initial header
    clearTerminal();
    cout << "=======================================\n";
    cout << "== Mutagen Puzzle Solver by Denevron ==\n";
    cout << "=======================================\n";
    cout << "Starting puzzle: " << PUZZLE_NUM << " (" << PUZZLE_NUM << "-bit)\n";
    cout << "Target HASH160: " << TARGET_HASH160.substr(0, 10) << "..."
         << TARGET_HASH160.substr(TARGET_HASH160.length() - 10) << "\n";
    cout << "Base Key: " << paddedKey << "\n";
    cout << "Flip count: " << FLIP_COUNT << " ";
    if (FLIP_COUNT != DEFAULT_FLIP_COUNT)
    {
        cout << "(override, default was " << DEFAULT_FLIP_COUNT << ")";
    }
    cout << "\n";
    cout << "Total Flips: " << to_string_128(total_combinations) << "\n";
    cout << "Using: " << WORKERS << " threads\n";
    cout << "\n";

    g_threadPrivateKeys.resize(WORKERS, "0");
    vector<thread> threads;

    // Distribute work evenly across threads
    AVXCounter total_combinations_avx;
    total_combinations_avx.store(total_combinations);

    AVXCounter comb_per_thread = AVXCounter::div(total_combinations_avx, WORKERS);
    uint64_t remainder = AVXCounter::mod(total_combinations_avx, WORKERS);

    for (int i = 0; i < WORKERS; i++)
    {
        __uint128_t start = i * comb_per_thread.load() + min((uint64_t)i, remainder);
        __uint128_t end = start + comb_per_thread.load() + (i < remainder ? 1 : 0);
        threads.emplace_back(worker, &secp, PUZZLE_NUM, FLIP_COUNT, i, start, end);
    }

    {
        lock_guard<mutex> lock(progress_mutex);
        if (total_checked_avx.load() >= total_combinations && !stop_event.load())
        {
            stop_event.store(true);
        }
    }

    // Wait for all threads to finish
    for (auto &t : threads)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    // Get final statistics
    __uint128_t final_count = total_checked_avx.load();
    globalElapsedTime = chrono::duration<double>(chrono::high_resolution_clock::now() - tStart).count();
    mkeysPerSec = (double)final_count / globalElapsedTime / 1e6;

    // Display results
    if (!results.empty())
    {
        auto [hex_key, checked, flips] = results.front();

        // Format the found key
        string compactHex = hex_key;
        size_t firstNonZeroHex = compactHex.find_first_not_of('0');
        if (firstNonZeroHex != string::npos)
        {
            compactHex = "0x" + compactHex.substr(firstNonZeroHex);
        }
        else
        {
            compactHex = "0x0";
        }

        // Print solution banner
        cout << "\n";
        cout << "=======================================\n";
        cout << "=========== SOLUTION FOUND ============\n";
        cout << "=======================================\n";
        cout << "Private key: " << compactHex << "\n";
        cout << "Checked " << to_string_128(checked) << " combinations\n";
        cout << "Bit flips: " << flips << endl;
        cout << "Time: " << fixed << setprecision(2) << globalElapsedTime << " seconds ("
             << formatElapsedTime(globalElapsedTime) << ")\n";
        cout << "Speed: " << fixed << setprecision(2) << mkeysPerSec << " Mkeys/s\n";
        cout << "Completion: " << fixed << setprecision(6)
             << (double)checked / total_combinations * 100.0 << "%\n";

        // Save solution to file
        ofstream out("puzzle_" + to_string(PUZZLE_NUM) + "_solution.txt");
        if (out)
        {
            out << hex_key;
            out.close();
            cout << "Solution saved to puzzle_" << PUZZLE_NUM << "_solution.txt\n";
        }
        else
        {
            cerr << "Warning: Failed to save solution to file!\n";
        }
    }
    else
    {
        // No solution found case
        cout << "\n\n=======================================\n";
        cout << "========= SEARCH COMPLETED ===========\n";
        cout << "=======================================\n";
        cout << "No solution found.\n";
        cout << "Checked " << to_string_128(final_count) << " combinations\n";
        cout << "Total possible: " << to_string_128(total_combinations) << "\n";
        cout << "Completion: " << fixed << setprecision(6)
             << (double)final_count / total_combinations * 100.0 << "%\n";
        cout << "Time: " << fixed << setprecision(2) << globalElapsedTime << " seconds ("
             << formatElapsedTime(globalElapsedTime) << ")\n";
        cout << "Speed: " << fixed << setprecision(2) << mkeysPerSec << " Mkeys/s\n";

        // Additional warning if not all combinations were checked
        if (final_count < total_combinations)
        {
            cout << "Warning: " << to_string_128(total_combinations - final_count)
                 << " combinations were not checked!\n";
        }

        return 0;
    }
}