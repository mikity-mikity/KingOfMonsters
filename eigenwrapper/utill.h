
#ifndef IDX2F
#define IDX2F(i,j,lda) ((((j)-1)*((size_t)lda))+((i)-1))
#endif /* IDX2F */

#ifndef IDX1F
#define IDX1F(i) ((i)-1)
#endif /* IDX1F */

static int enablePeerAccess(const int nbGpus, const int* deviceList)
{
    int currentDevice = 0;
    cudaGetDevice(&currentDevice);

    /* Remark: access granted by this cudaDeviceEnablePeerAccess is unidirectional */
    /* Rows and columns represents a connectivity matrix between GPUs in the system */
    for (int row = 0; row < nbGpus; row++) {
        cudaSetDevice(row);
        for (int col = 0; col < nbGpus; col++) {
            if (row != col) {
                cudaError_t cudaStat1 = cudaSuccess;
                cudaError_t cudaStat2 = cudaSuccess;
                int canAccessPeer = 0;
                cudaStat1 = cudaDeviceCanAccessPeer(&canAccessPeer, row, col);
                if (canAccessPeer) {
                    printf("\t Enable peer access from gpu %d to gpu %d\n", row, col);
                    cudaStat2 = cudaDeviceEnablePeerAccess(col, 0);
                }
                assert(cudaStat1 == cudaSuccess);
                assert(cudaStat2 == cudaSuccess);
            }
        }
    }
    cudaSetDevice(currentDevice);
    return 0;
}

static int workspaceFree(
    int num_devices,
    const int* deviceIdA, /* <int> dimension num_devices */
    void** array_d_work  /* <t> num_devices, host array */
                         /* array_d_work[j] points to device workspace of device j */
)
{
    int currentDev = 0; /* record current device ID */
    cudaGetDevice(&currentDev);

    for (int idx = 0; idx < num_devices; idx++) {
        int deviceId = deviceIdA[idx];
        /* WARNING: we need to set device before any runtime API */
        cudaSetDevice(deviceId);

        if (NULL != array_d_work[idx]) {
            cudaFree(array_d_work[idx]);
        }
    }
    cudaSetDevice(currentDev);
    return 0;
}

static int workspaceAlloc(
    int num_devices,
    const int* deviceIdA, /* <int> dimension num_devices */
    size_t sizeInBytes,  /* number of bytes per device */
    void** array_d_work  /* <t> num_devices, host array */
                         /* array_d_work[j] points to device workspace of device j */
)
{
    cudaError_t cudaStat1 = cudaSuccess;

    int currentDev = 0; /* record current device ID */
    cudaGetDevice(&currentDev);

    memset(array_d_work, 0, sizeof(void*) * num_devices);
    for (int idx = 0; idx < num_devices; idx++) {
        int deviceId = deviceIdA[idx];
        /* WARNING: we need to set device before any runtime API */
        cudaSetDevice(deviceId);

        void* d_workspace = NULL;

        cudaStat1 = cudaMalloc(&d_workspace, sizeInBytes);
        assert(cudaSuccess == cudaStat1);
        array_d_work[idx] = d_workspace;
    }
    cudaSetDevice(currentDev);
    return 0;
}

/* create a empty matrix A with A := 0 */
template <typename T_ELEM>
int createMat(
    int num_devices,
    const int* deviceIdA, /* <int> dimension num_devices */
    int N_A,   /* number of columns of global A */
    int T_A,   /* number of columns per column tile */
    int LLD_A, /* leading dimension of local A */
    T_ELEM** array_d_A  /* host pointer array of dimension num_devices */
)
{
    printf("aa-1\n");
    cudaError_t cudaStat1 = cudaSuccess;
    int currentDev = 0; /* record current device id */
    cudaGetDevice(&currentDev);
    cudaDeviceSynchronize();
    printf("aa-2\n");
    printf("%d\n", num_devices);
    printf("%d", T_A);
    const int A_num_blks = (N_A + T_A - 1) / T_A;
    const int max_A_num_blks_per_device = (A_num_blks + num_devices - 1) / num_devices;
    /* Allocate base pointers */
    printf("aa-3\n");
    memset(array_d_A, 0, sizeof(T_ELEM*) * num_devices);
    for (int p = 0; p < num_devices; p++) {
        cudaStat1 = cudaSetDevice(deviceIdA[p]);
        assert(cudaSuccess == cudaStat1);
        /* Allocate max_A_num_blks_per_device blocks per device */
        cudaStat1 = cudaMalloc(&(array_d_A[p]), sizeof(T_ELEM) * LLD_A * T_A * max_A_num_blks_per_device);
        assert(cudaSuccess == cudaStat1);
        /* A := 0 */
        cudaStat1 = cudaMemset(array_d_A[p], 0, sizeof(T_ELEM) * LLD_A * T_A * max_A_num_blks_per_device);
        assert(cudaSuccess == cudaStat1);
    }
    printf("aa-4\n");
    cudaDeviceSynchronize();
    cudaSetDevice(currentDev);
    return 0;
}

static int destroyMat(
    int num_devices,
    const int* deviceIdA, /* <int> dimension num_devices */
    int N_A,  /* number of columns of global A */
    int T_A,  /* number of columns per column tile */
    void** array_d_A)  /* host pointer array of dimension num_devices */
{
    cudaError_t cudaStat = cudaSuccess;

    int currentDev = 0; /* record current device id */
    cudaGetDevice(&currentDev);

    const int num_blocks = (N_A + T_A - 1) / T_A;
    for (int p = 0; p < num_devices; p++) {
        cudaStat = cudaSetDevice(deviceIdA[p]);
        assert(cudaSuccess == cudaStat);

        if (NULL != array_d_A[p]) {
            cudaFree(array_d_A[p]);
        }
    }
    memset(array_d_A, 0, sizeof(void*) * num_devices);
    cudaSetDevice(currentDev);
    return 0;
}

template <typename T_ELEM>
static int mat_pack2unpack(
    int num_devices,
    int N_A,   /* number of columns of global A */
    int T_A,   /* number of columns per column tile */
    int LLD_A, /* leading dimension of local A */
    T_ELEM** array_d_A_packed,  /* host pointer array of dimension num_devices */
/* output */
T_ELEM** array_d_A_unpacked /* host pointer array of dimension num_blks */
)
{
    const int num_blks = (N_A + T_A - 1) / T_A;

    for (int p_a = 0; p_a < num_devices; p_a++) {
        T_ELEM* d_A = array_d_A_packed[p_a];
        int nz_blks = 0;
        for (int JA_blk_id = p_a; JA_blk_id < num_blks; JA_blk_id += num_devices) {
            array_d_A_unpacked[JA_blk_id] = d_A + (size_t)LLD_A * T_A * nz_blks;
            nz_blks++;
        }
    }
    return 0;
}


template <typename T_ELEM>
static int memcpyH2D(
    int num_devices,
    const int* deviceIdA, /* <int> dimension num_devices */
    int M,  /* number of rows in local A, B */
    int N,  /* number of columns in local A, B */
/* input */
const T_ELEM* h_B,  /* host array, h_B is M-by-N with leading dimension ldb  */
int ldb,
/* output */
int N_A,   /* number of columns of global A */
int T_A,  /* number of columns per column tile */
int LLD_A, /* leading dimension of local A */
T_ELEM** array_d_A_packed, /* host pointer array of dimension num_devices */
int IA,  /* base-1 */
int JA   /* base-1 */
)
{
    cudaError_t cudaStat1 = cudaSuccess;

    int currentDev = 0; /* record current device id */

/*  Quick return if possible */
    if ((0 >= M) || (0 >= N)) {
        return 0;
    }

    /* consistent checking */
    if (ldb < M) {
        return 1;
    }

    cudaGetDevice(&currentDev);
    cudaDeviceSynchronize();

    const int num_blks = (N_A + T_A - 1) / T_A;

    T_ELEM** array_d_A_unpacked = (T_ELEM**)malloc(sizeof(T_ELEM*) * num_blks);
    assert(NULL != array_d_A_unpacked);

    mat_pack2unpack<T_ELEM>(
        num_devices,
        N_A,   /* number of columns of global A */
        T_A,   /* number of columns per column tile */
        LLD_A, /* leading dimension of local A */
        array_d_A_packed,  /* host pointer array of size num_devices */
/* output */
array_d_A_unpacked /* host pointer arrya of size num_blks */
);

    /* region of interest is A(IA:IA+N-1, JA:JA+N-1) */
    const int N_hat = (JA - 1) + N; /* JA is base-1 */

    const int JA_start_blk_id = (JA - 1) / T_A;
    const int JA_end_blk_id = (N_hat - 1) / T_A;

    for (int p_a = 0; p_a < num_devices; p_a++) {
        /* region of interest: JA_start_blk_id:1:JA_end_blk_id */
        for (int JA_blk_id = p_a; JA_blk_id <= JA_end_blk_id; JA_blk_id += num_devices) {
            if (JA_blk_id < JA_start_blk_id) { continue; }
            /*
             * process column block of A
             *       A(A_start_row:M_A, A_start_col : (A_start_col + IT_A-1) )
             */
            const int IBX_A = (1 + JA_blk_id * T_A); /* base-1 */
            const int A_start_col = std::max(JA, IBX_A);   /* base-1 */
            const int A_start_row = IA;  /* base-1 */

            const int bdd = std::min(N_hat, (IBX_A + T_A - 1));
            const int IT_A = std::min(T_A, (bdd - A_start_col + 1));

            const int loc_A_start_row = A_start_row;   /* base-1 */
            const int loc_A_start_col = (A_start_col - IBX_A) + 1;  /* base-1 */

            T_ELEM* d_A = array_d_A_unpacked[JA_blk_id] + IDX2F(loc_A_start_row, loc_A_start_col, LLD_A);
            const T_ELEM* h_A = h_B + IDX2F(A_start_row - IA + 1, A_start_col - JA + 1, ldb);

            cudaStat1 = cudaMemcpy2D(
                d_A,  /* dst */
                (size_t)LLD_A * sizeof(T_ELEM),
                h_A,  /* src */
                (size_t)ldb * sizeof(T_ELEM),
                (size_t)M * sizeof(T_ELEM),
                (size_t)IT_A,
                cudaMemcpyHostToDevice
            );
            assert(cudaSuccess == cudaStat1);
        }/* for each tile per device */
    }/* for each device */
    cudaDeviceSynchronize();
    cudaSetDevice(currentDev);

    if (NULL != array_d_A_unpacked) { free(array_d_A_unpacked); }
    return 0;
}

template <typename T_ELEM>
static int memcpyD2H(
    int num_devices,
    const int* deviceIdA, /* <int> dimension num_devices */
    int M,  /* number of rows in local A, B */
    int N,  /* number of columns in local A, B */
 /* input */
    int N_A,  /* number of columns of global A */
    int T_A,  /* number of columns per column tile */
    int LLD_A, /* leading dimension of local A */
    T_ELEM** array_d_A_packed, /* host pointer array of dimension num_devices */
    int IA,  /* base-1 */
    int JA,   /* base-1 */
/* output */
T_ELEM* h_B,  /* host array, h_B is M-by-N with leading dimension ldb  */
int ldb
)
{
    cudaError_t cudaStat1 = cudaSuccess;
    int currentDev = 0; /* record current device id */

/*  Quick return if possible */
    if ((0 >= M) || (0 >= N)) {
        return 0;
    }
    /* consistent checking */
    if (ldb < M) {
        return 1;
    }
    cudaGetDevice(&currentDev);
    cudaDeviceSynchronize();

    const int num_blks = (N_A + T_A - 1) / T_A;
    T_ELEM** array_d_A_unpacked = (T_ELEM**)malloc(sizeof(T_ELEM*) * num_blks);
    assert(NULL != array_d_A_unpacked);

    mat_pack2unpack<T_ELEM>(
        num_devices,
        N_A,   /* number of columns of global A */
        T_A,   /* number of columns per column tile */
        LLD_A, /* leading dimension of local A */
        array_d_A_packed,  /* host pointer array of size num_devices */
        array_d_A_unpacked /* host pointer arrya of size num_blks */
        );
    /* region of interest is A(IA:IA+N-1, JA:JA+N-1) */
    const int N_hat = (JA - 1) + N; /* JA is base-1 */
    const int JA_start_blk_id = (JA - 1) / T_A;
    const int JA_end_blk_id = (N_hat - 1) / T_A;
    for (int p_a = 0; p_a < num_devices; p_a++) {
        /* region of interest: JA_start_blk_id:1:JA_end_blk_id */
        for (int JA_blk_id = p_a; JA_blk_id <= JA_end_blk_id; JA_blk_id += num_devices) {
            if (JA_blk_id < JA_start_blk_id) { continue; }
            /* process column block, A(A_start_row:M_A, A_start_col : (A_start_col + IT_A-1) ) */
            const int IBX_A = (1 + JA_blk_id * T_A); /* base-1 */
            const int A_start_col = std::max(JA, IBX_A);   /* base-1 */
            const int A_start_row = IA;  /* base-1 */
            const int bdd = std::min(N_hat, (IBX_A + T_A - 1));
            const int IT_A = std::min(T_A, (bdd - A_start_col + 1));
            const int loc_A_start_row = A_start_row;   /* base-1 */
            const int loc_A_start_col = (A_start_col - IBX_A) + 1;  /* base-1 */
            const T_ELEM* d_A = array_d_A_unpacked[JA_blk_id] + IDX2F(loc_A_start_row, loc_A_start_col, LLD_A);
            T_ELEM* h_A = h_B + IDX2F(A_start_row - IA + 1, A_start_col - JA + 1, ldb);
            cudaStat1 = cudaMemcpy2D(
                h_A,  /* dst */
                (size_t)ldb * sizeof(T_ELEM),
                d_A,  /* src */
                (size_t)LLD_A * sizeof(T_ELEM),
                (size_t)M * sizeof(T_ELEM),
                (size_t)IT_A,
                cudaMemcpyDeviceToHost
            );
            assert(cudaSuccess == cudaStat1);
        }/* for each tile per device */
    }/* for each device */
    cudaDeviceSynchronize();
    cudaSetDevice(currentDev);
    if (NULL != array_d_A_unpacked) { free(array_d_A_unpacked); }
    return 0;
}
