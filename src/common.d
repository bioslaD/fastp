__gshared string FASTP_VER="0.13.2";

__gshared bool _DEBUG = false;

// the limit of the queue to store the packs
// error may happen if it generates more packs than this number
static const int PACK_NUM_LIMIT  = 1000000;

// how many reads one pack has
static const int PACK_SIZE = 1000;

// if one pack is produced, but not consumed, it will be kept in the memory
// this number limit the number of in memory packs
// if the number of in memory packs is full, the producer thread should sleep
static const int PACK_IN_MEM_LIMIT = 500;

// if read number is more than this, warn it
static const int WARN_STANDALONE_READ_LIMIT = 10000;

// different filtering results, bigger number means worse
// if r1 and r2 are both failed, then the bigger one of the two results will be recorded
// we reserve some gaps for future types to be added
static const int PASS_FILTER = 0;
static const int FAIL_POLY_X = 4;
static const int FAIL_OVERLAP = 8;
static const int FAIL_N_BASE = 12;
static const int FAIL_LENGTH = 16;
static const int FAIL_QUALITY = 20;
static const int FAIL_COMPLEXITY = 24;

// how many types in total we support
static const int FILTER_RESULT_TYPES = 32;