import numpy as np
cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def dist(double [:, :] A, double [:, :] A2, double [:,:] name_A, double [:,:] name_A2):
    cdef:
        Py_ssize_t nrow = A.shape[0]
        Py_ssize_t ncol = A.shape[1]
        Py_ssize_t nrow2 = A2.shape[0]
        Py_ssize_t row_A1, row_A2, seg
        np.ndarray[np.float64_t, ndim=2] D = np.zeros((nrow*nrow2, 7), np.double)
        double tmp_SD, tmp_DD, tmp_AD, dn, tmp1_root, tmp2_root

    cdef int idx = 0

    for row_A1 in range(nrow):
        for row_A2 in range(nrow2):
            tmp_SD = 0
            tmp_DD = 0
            tmp_AD = 0
            # tmp1_root = 0
            # tmp2_root = 0
            dn = 0
            for seg in range(ncol):
                if A[row_A1, seg] == A2[row_A2, seg]:
                    tmp_SD += 1
                elif (A[row_A1, seg] + A2[row_A2, seg]) % 2 == 0:
                    tmp_AD += 1
                else:
                    tmp_DD += abs(A[row_A1, seg] - A2[row_A2, seg])
                    dn += 1
                # tmp1_root += abs(A[row_A1, seg] - 2)
                # tmp2_root += abs(A2[row_A2, seg] - 2)

            D[idx, 0] = name_A[0,row_A1]
            D[idx, 1] = name_A2[0,row_A2]
            D[idx, 2] = tmp_SD
            D[idx, 3] = tmp_AD
            #D[idx, 4] = tmp_DD
            if dn == 0:
                D[idx, 4] = 0
            else:
                D[idx, 4] = tmp_DD / dn
            # D[idx, 5] = tmp1_root
            # D[idx, 6] = tmp2_root
            idx += 1
        #print(aa)
        #aa += 1
    return D



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def dist2(double [:, :] A, double [:, :] A2,
          double [:,:] name_A, double [:,:] name_A2):
    cdef:
        Py_ssize_t nrow = A.shape[0]
        Py_ssize_t ncol = A.shape[1]
        Py_ssize_t row_A1, row_A2, seg
        np.ndarray[np.float64_t, ndim=2] D = np.zeros((nrow, 7), np.double)
        double tmp_SD, tmp_DD, tmp_AD, dn,  tmp1_root, tmp2_root

    cdef int idx = 0

    for i in range(nrow):
        tmp_SD = 0
        tmp_DD = 0
        tmp_AD = 0
        # tmp1_root = 0
        # tmp2_root = 0
        dn = 0
        for seg in range(ncol):
            if A[i, seg] == A2[i, seg]:
                tmp_SD += 1
            elif (A[i, seg] + A2[i, seg]) % 2 == 0:
                tmp_AD += 1
            else:
                tmp_DD += abs(A[i, seg] - A2[i, seg])
                dn += 1
            # tmp1_root += abs(A[i, seg] - 2)
            # tmp2_root += abs(A2[i, seg] - 2)

        D[idx, 0] = name_A[0,i]
        D[idx, 1] = name_A2[0,i]
        D[idx, 2] = tmp_SD
        D[idx, 3] = tmp_AD
        # D[idx, 4] = tmp_DD
        if dn == 0:
            D[idx, 4] = 0
        else:
            D[idx, 4] = tmp_DD / dn
        # D[idx, 5] = tmp1_root
        # D[idx, 6] = tmp2_root
        idx += 1
        #print(aa)
        #aa += 1
    return D