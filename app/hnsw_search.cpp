//
// Created by Zhen Peng on 01/24/2021.
//

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <unordered_set>
#include <sys/resource.h>
#include "../hnswlib/hnswlib.h"

using idi = unsigned;

class StopW {
    std::chrono::steady_clock::time_point time_begin;
public:
    StopW() {
        time_begin = std::chrono::steady_clock::now();
    }

    float getElapsedTimeMicro() {
        std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
        return (std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_begin).count());
    }

    void reset() {
        time_begin = std::chrono::steady_clock::now();
    }

};

void load_true_NN(
        const char *filename,
        const unsigned num_queries,
        unsigned int *massQA)
//        std::vector< std::vector<idi> > &true_nn_list)
//        unsigned &t_K)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        fprintf(stderr, "Error: cannot open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
    idi t_query_num;
    idi t_K;
//    unsigned t_K;
    fin.read(reinterpret_cast<char *>(&t_query_num), sizeof(t_query_num));
    fin.read(reinterpret_cast<char *>(&t_K), sizeof(t_K));
//    if (t_query_num != query_num) {
//        fprintf(stderr, "Error: query_num %u is not equal to the record %u in true-NN file %s\n",
//                query_num, t_query_num, filename);
//        exit(EXIT_FAILURE);
//    }
    if (t_query_num != num_queries) {
        fprintf(stderr, "Error: t_query_num %u is not num_queries_ %u\n", t_query_num, num_queries);
        exit(EXIT_FAILURE);
    }
    if (t_K != 100) {
        fprintf(stderr, "Error: t_K %u is not 100.\n", t_K);
        exit(EXIT_FAILURE);
    }

////    data = new unsigned[(size_t) t_query_num * (size_t) t_K];
//    true_nn_list.resize(t_query_num);
//    for (idi q_i = 0; q_i < t_query_num; ++q_i) {
//        true_nn_list[q_i].resize(t_K);
//    }

    for (unsigned q_i = 0; q_i < t_query_num; ++q_i) {
        size_t offset = q_i * t_K;
        for (unsigned n_i = 0; n_i < t_K; ++n_i) {
            unsigned id;
            float dist;
            fin.read(reinterpret_cast<char *>(&id), sizeof(id));
            fin.read(reinterpret_cast<char *>(&dist), sizeof(dist));
            massQA[offset + n_i] = id;
//            true_nn_list[q_i][n_i] = id;
        }
    }

    fin.close();
}

void load_data(
        const char *filename,
        float *data,
        size_t &num,
        size_t &dim)
{  // load data with sift10K pattern
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        fprintf(stderr, "Error: cannot open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
//    in.read((char*)&dim, 4);
    uint32_t t_d;
    in.read((char*) &t_d, 4);
    dim = (size_t) t_d;
//    std::cout << "data dimension: " << dim << std::endl;
    in.seekg(0, std::ios::end);
    std::ios::pos_type ss = in.tellg();
    size_t fsize = (size_t)ss;
    num = (fsize / (dim + 1) / 4);
//    data = new float[static_cast<uint64_t>(num) * static_cast<uint64_t>(dim)];

    in.seekg(0, std::ios::beg);
    for (size_t i = 0; i < num; i++) {
        in.seekg(4, std::ios::cur);
        in.read((char*)(data + i * dim), dim * 4);
    }
    in.close();
}

/**
* Returns the peak (maximum so far) resident set size (physical
* memory use) measured in bytes, or zero if the value cannot be
* determined on this OS.
*/
static size_t getPeakRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
        return (size_t)0L;      /* Can't open? */
    if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo))
    {
        close(fd);
        return (size_t)0L;      /* Can't read? */
    }
    close(fd);
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t) (rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}


/**
* Returns the current resident set size (physical memory use) measured
* in bytes, or zero if the value cannot be determined on this OS.
*/
static size_t getCurrentRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount) != KERN_SUCCESS)
        return (size_t)0L;      /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE *fp = NULL;
    if ((fp = fopen("/proc/self/statm", "r")) == NULL)
        return (size_t) 0L;      /* Can't open? */
    if (fscanf(fp, "%*s%ld", &rss) != 1) {
        fclose(fp);
        return (size_t) 0L;      /* Can't read? */
    }
    fclose(fp);
    return (size_t) rss * (size_t) sysconf(_SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;          /* Unsupported. */
#endif
}

static void get_gt(
        const unsigned int *massQA,
//        unsigned char *massQ,
//        unsigned char *mass,
//        size_t vecsize,
        size_t qsize,
//        hnswlib::L2SpaceI &l2space,
//        size_t vecdim,
        std::vector<std::priority_queue<std::pair<float, hnswlib::labeltype >>> &answers,
        size_t k)
{
    (std::vector<std::priority_queue<std::pair<float, hnswlib::labeltype >>>(qsize)).swap(answers);
//    hnswlib::DISTFUNC<int> fstdistfunc_ = l2space.get_dist_func();
    size_t t_k = 100;
//    std::cout << qsize << "\n";
    for (int i = 0; i < qsize; i++) {
        for (int j = 0; j < k; j++) {
            answers[i].emplace(0.0f, massQA[t_k * i + j]);
        }
    }
}

static float test_approx(
        const float *massQ,
        size_t vecsize,
        size_t qsize,
        hnswlib::HierarchicalNSW<float> &appr_alg,
        size_t vecdim,
        std::vector<std::priority_queue<std::pair<float, hnswlib::labeltype >>> &answers,
        size_t k) {
    size_t correct = 0;
    size_t total = 0;
    //uncomment to test in parallel mode:
    //#pragma omp parallel for
    for (int i = 0; i < qsize; i++) {

        std::priority_queue<std::pair<float, hnswlib::labeltype >> result = appr_alg.searchKnn(massQ + vecdim * i, k);
        std::priority_queue<std::pair<float, hnswlib::labeltype >> gt(answers[i]);
        std::unordered_set<hnswlib::labeltype> g;
        total += gt.size();

        while (gt.size()) {


            g.insert(gt.top().second);
            gt.pop();
        }

        while (result.size()) {
            if (g.find(result.top().second) != g.end()) {

                correct++;
            } else {
            }
            result.pop();
        }

    }
    return 1.0f * correct / total;
}

static void test_vs_recall(
        const float *massQ,
        size_t vecsize,
        size_t qsize,
        hnswlib::HierarchicalNSW<float> &appr_alg,
        size_t vecdim,
        std::vector<std::priority_queue<std::pair<float, hnswlib::labeltype >>> &answers,
        size_t k)
{
    std::vector<size_t> efs;// = { 10,10,10,10,10 };
    for (int i = k; i < 30; i++) {
        efs.push_back(i);
    }
    for (int i = 30; i < 100; i += 10) {
        efs.push_back(i);
    }
    for (int i = 100; i < 500; i += 40) {
        efs.push_back(i);
    }
    for (size_t ef : efs) {
        appr_alg.setEf(ef);
        StopW stopw = StopW();

        float recall = test_approx(massQ, vecsize, qsize, appr_alg, vecdim, answers, k);
        float time_us_per_query = stopw.getElapsedTimeMicro() / qsize;

        std::cout << ef << "\t" << recall << "\t" << time_us_per_query << " us\n";
        if (recall > 1.0) {
            std::cout << recall << "\t" << time_us_per_query << " us\n";
            break;
        }
    }
}

void search(
        hnswlib::HierarchicalNSW<float> *appr_alg,
        const float *massQ,
        const unsigned int *massQA,
        const int subset_size_milllions,
        const size_t vecdim,
        const size_t qsize)
//        const char *path_data,
//        const char *path_index,
//        const int M,
//        const int efConstruction)
{


//    int subset_size_milllions = 200;
//    int efConstruction = 40; // ?
//    int M = 16; // ?

    size_t vecsize = subset_size_milllions * 1000000;

//    size_t qsize = 10000;
//    size_t vecdim = 128;
//    char path_index[1024];
//    char path_gt[1024];
//    char *path_q = "../bigann/bigann_query.bvecs";
//    char *path_data = "../bigann/bigann_base.bvecs";
//    sprintf(path_index, "sift1b_%dm_ef_%d_M_%d.bin", subset_size_milllions, efConstruction, M);

//    sprintf(path_gt, "../bigann/gnd/idx_%dM.ivecs", subset_size_milllions);


//    unsigned char *massb = new unsigned char[vecdim];
//    float massb = new float[vecdim];
//
//    std::cout << "Loading GT:\n";
//    std::ifstream inputGT(path_gt, std::ios::binary);
//    unsigned int *massQA = new unsigned int[qsize * 1000];
//    for (int i = 0; i < qsize; i++) {
//        int t;
//        inputGT.read((char *) &t, 4);
//        inputGT.read((char *) (massQA + 1000 * i), t * 4);
//        if (t != 1000) {
//            std::cout << "err";
//            return;
//        }
//    }
//    inputGT.close();

//    std::cout << "Loading queries:\n";
//    unsigned char *massQ = new unsigned char[qsize * vecdim];
//    std::ifstream inputQ(path_q, std::ios::binary);
//
//    for (int i = 0; i < qsize; i++) {
//        int in = 0;
//        inputQ.read((char *) &in, 4);
//        if (in != 128) {
//            std::cout << "file error";
//            exit(1);
//        }
//        inputQ.read((char *) massb, in);
//        for (int j = 0; j < vecdim; j++) {
//            massQ[i * vecdim + j] = massb[j];
//        }
//
//    }
//    inputQ.close();


//    unsigned char *mass = new unsigned char[vecdim];
//    std::ifstream input(path_data, std::ios::binary);
//    if (!input.is_open()) {
//        fprintf(stderr, "Error: cannot open file %s .\n", path_data);
//        exit(EXIT_FAILURE);
//    }
//    int in = 0;
//    hnswlib::L2SpaceI l2space(vecdim);
    hnswlib::L2Space l2space(vecdim);

//    hnswlib::HierarchicalNSW<int> *appr_alg;
//    hnswlib::HierarchicalNSW<float> *appr_alg = new hnswlib::HierarchicalNSW<float>(&l2space, vecsize, M, efConstruction);
//    if (exists_test(path_index)) {
////        std::cout << "Loading index from " << path_index << ":\n";
////        appr_alg = new hnswlib::HierarchicalNSW<int>(&l2space, path_index, false);
////        std::cout << "Actual memory usage: " << getCurrentRSS() / 1000000 << " Mb \n";
//        printf("File %s exits.\n", path_index);
//    } else {
//        std::cout << "Building index:\n";
////        appr_alg = new hnswlib::HierarchicalNSW<int>(&l2space, vecsize, M, efConstruction);
//        appr_alg = new hnswlib::HierarchicalNSW<float>(&l2space, vecsize, M, efConstruction);
//
//        {
//            float *mass = (float *) malloc(vecdim * sizeof(float));
//            input.read((char *) &in, 4);
//            if (in != vecdim) {
//                std::cout << "file error";
//                exit(1);
//            }
////        input.read((char *) massb, in);
////
////        for (int j = 0; j < vecdim; j++) {
////            mass[j] = massb[j] * (1.0f);
////        }
////        appr_alg->addPoint((void *) (massb), (size_t) 0); // ??? massb or mass? Strange
//            input.read(reinterpret_cast<char *>(mass), 4 * vecdim);
//            appr_alg->addPoint((void *) mass, (size_t) 0);
//            free(mass);
//        }
//        int j1 = 0;
//        StopW stopw = StopW();
//        StopW stopw_full = StopW();
//        size_t report_every = 100000;
////        float *mass = (float *) malloc(vecdim * sizeof(float));
//#pragma omp parallel for
//        for (int i = 1; i < vecsize; i++) {
////            unsigned char mass[128];
//            float *mass = (float *) malloc(vecdim * sizeof(float));
//            int j2=0;
//#pragma omp critical
//            {
//
//                input.read((char *) &in, 4);
//                if (in != vecdim) {
//                    std::cout << "file error";
//                    exit(1);
//                }
////                input.read((char *) massb, in);
////                for (int j = 0; j < vecdim; j++) {
////                    mass[j] = massb[j];
////                }
//                input.read(reinterpret_cast<char *>(mass), 4 * vecdim);
//                j1++;
//                j2=j1;
//                if (j1 % report_every == 0) {
//                    std::cout << j1 / (0.01 * vecsize) << " %, "
//                              << report_every / (1000.0 * 1e-6 * stopw.getElapsedTimeMicro()) << " kips " << " Mem: "
//                              << getCurrentRSS() / 1000000 << " Mb \n";
//                    stopw.reset();
//                }
//            }
//            appr_alg->addPoint((void *) (mass), (size_t) j2);
//
//            free(mass);
//        }
////        free(mass);
//        input.close();
//        std::cout << "Build time:" << 1e-6 * stopw_full.getElapsedTimeMicro() << "  seconds\n";
//        appr_alg->saveIndex(path_index);
//    }


    std::vector<std::priority_queue<std::pair<float, hnswlib::labeltype >>> answers;
    size_t k = 100; // The K as in K-NN
    std::cout << "Parsing gt:\n";
    get_gt(
            massQA,
            qsize,
            answers,
            k);
    std::cout << "Loaded gt\n";
    for (int i = 0; i < 1; i++)
        test_vs_recall(massQ, vecsize, qsize, *appr_alg, vecdim, answers, k);
    std::cout << "Actual memory usage: " << getCurrentRSS() / 1000000 << " Mb \n";
}

void usage(int argc, char *argv[])
{
    if (argc != 7) {
        fprintf(stderr,
                "Usage: %s <index_file> <query_file> <groundtrue_file> "
                "<size_in_millions> <dimension> <num_queries>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[])
{
    usage(argc, argv);
    setbuf(stdout, nullptr); // Remove stdout buffer.

//    const char *path_data = argv[1]; // data file
    const char *path_index = argv[1]; // index
    const char *path_q = argv[2]; // query
    const char *path_gt = argv[3];
    const int subset_size_milllions = strtoull(argv[4], nullptr, 0); // number of vectors in millions
    const size_t vecdim = strtoull(argv[5], nullptr, 0); // dimention of vector
    const size_t qsize = strtoull(argv[6], nullptr, 0); // number of queries
//    const int M = strtoull(argv[5], nullptr, 0);
//    const int efConstruction = strtoull(argv[6], nullptr, 0);

    // Read ground truth
    printf("loading gt...\n");
    unsigned int *massQA = (unsigned *) malloc(qsize * 100 * sizeof(unsigned ));
    load_true_NN(
            path_gt,
            qsize,
            massQA);

    // Read queries
    printf("loading queries...\n");
    float *massQ = (float *) malloc(qsize * vecdim * sizeof(float));
    {
        size_t num_q;
        size_t dim;
        load_data(
                path_q,
                massQ,
                num_q,
                dim);
        if (num_q != qsize) {
            fprintf(stderr, "Error: query file wrong num_q %lu (qsize %lu).\n", num_q, qsize);
            exit(EXIT_FAILURE);
        }
        if (dim != vecdim) {
            fprintf(stderr, "Error: query file wrong dim %lu (vecdim %lu).\n", dim, vecdim);
            exit(EXIT_FAILURE);
        }
    }

    hnswlib::L2Space l2space(vecdim);
    printf("loading index...\n");
    hnswlib::HierarchicalNSW<float> *appr_alg = new hnswlib::HierarchicalNSW<float>(&l2space, path_index, false);

    search(
            appr_alg,
            massQ,
            massQA,
            subset_size_milllions,
            vecdim,
            qsize);
//            path_data,
//            path_index,
//            M,
//            efConstruction);

    // Cleanup
    {
        free(massQA);
        free(massQ);
        delete appr_alg;
    }

    return EXIT_SUCCESS;
}

