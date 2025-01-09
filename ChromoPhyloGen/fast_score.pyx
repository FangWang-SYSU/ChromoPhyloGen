#import numpy as np
import numpy as pynp
cimport numpy as cnp
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef double skew(double[:] data):
    cdef int n = data.shape[0]
    cdef double mean = pynp.mean(data)
    cdef double std = pynp.std(data, ddof=0)
    cdef double skewness = 0.0

    # skewness = np.mean(np.sum(((np.array(data) - np.array(mean))/std)**3))
    for i in range(n):
        skewness += ((data[i] - mean) / std) ** 3

    skewness *= (1.0 / n)

    return skewness


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef dict fast_score(double [:]  cnv_data,
                      double [:, :]  chr_pos):
    cdef int min_length = 0  #
    cdef int i
    cdef double breakpoint_prop, uniformity_score, avg_amplitude, rearrange_score
    cdef cnp.ndarray[int, ndim=1] bkps
    cdef cnp.ndarray[cnp.double_t, ndim=2] segments

    # breakpoint
    bkps = pynp.where(pynp.diff(cnv_data) != 0)[0]
    # segment
    segments = pynp.zeros((len(bkps) + 1, 3))
    segments[:, 0] = pynp.concatenate(([0], bkps))#chr_pos[np.concatenate(([0], bkps)), 1]# np.concatenate(([0], bkps)) # start
    segments[:, 1] = pynp.concatenate((bkps, [len(cnv_data) - 1]))#chr_pos[np.concatenate((bkps, [len(cnv_data) - 1])), 2]#np.concatenate((bkps, [len(cnv_data) - 1])) # end
    segments[:, 2] = segments[:, 1] - segments[:, 0] + 1 # len

    # segments$chr_start = as.numeric(chr_pos[segments$start, 'start'])
    # segments$chr_end = as.numeric(chr_pos[segments$end, 'end'])
    # segments$chr_len = segments$chr_end - segments$chr_start

    #segments[:, 3] = cnv_data[segments[:,0]] # cnv
    # estimate
    try:
        if segments.shape[0] <=1:
            rearrange_score = 0
        else:
            breakpoint_count = segments.shape[0]
            breakpoint_prop = breakpoint_count / pynp.sum(segments[:, 2])
            uniformity_score = -skew(segments[:, 2])
            # sigmoid
            uniformity_score = 1 / (1 + pynp.exp(-uniformity_score))
            avg_amplitude = pynp.var(cnv_data)
            rearrange_score = breakpoint_prop * uniformity_score * avg_amplitude
    except:
        rearrange_score = 0

    return {'smooth_segment': segments.tolist(),
            'rearrange_score': rearrange_score}



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cal_limit(double[:] x, double coverage_thr=0.6):
    cdef dict counter
    cdef str result
    cdef double top1_value, top2_value, element,i,count, all_count
    cdef list counter_list
    counter = {}
    # 遍历列表，统计元素出现次数
    for element in x:
        if element in counter:
            counter[element] += 1
        else:
            counter[element] = 1
    # unique_elements = np.unique(np.array(x), return_counts=True)

    if len(counter) <= 1:
        return -1 #'neutral'

    counter_list = list(counter.items())
    counter_list.sort(key=lambda a: -a[1])
    top1_value = counter_list[0][0]
    # top2_value = counter_list[1][1]
    # top2_value = np.array(x)-np.array(top1_value)
    # count = np.sum(np.isin(top2_value, [1, -1])) / np.sum(top2_value!=0)
    all_count = 0
    for i in x:
        i = i - top1_value
        if i != 0:
            all_count += 1
            if i==1 or i==-1:
                count += 1.0

    # # if count / len(x) >= coverage_thr:
    # #     return 'limit'
    return count/all_count #'seismic'



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef permutation_bg(double[:] data,
                    double [:, :]  chr_pos,
                    double lbd, double possion_lbd,  int size, int seed=1):
    cdef int ncol = data.shape[0]
    cdef cnp.ndarray[cnp.double_t, ndim=1] new_seg, new_loc
    # cdef np.ndarray[np.int_t, ndim=1] new_loc2
    cdef list tmp_res
    cdef list new_data, seg_len, last_cn
    cdef int i, s, e, tmp_len
    cdef float mean_value
    cdef cnp.ndarray[cnp.double_t, ndim=2] segments
    cdef double breakpoint_prop, uniformity_score, avg_amplitude, rearrange_score

    pynp.random.seed(seed)
    seg_size = pynp.random.poisson(possion_lbd)
    new_seg = pynp.random.exponential(lbd, size=max(1,seg_size))
    if len(new_seg)<=2:
        return pynp.inf

    new_seg /= pynp.sum(new_seg)
    new_loc = pynp.cumsum(pynp.round(new_seg * ncol))
    new_loc = new_loc[new_loc != 0]
    new_loc = new_loc[new_loc != len(data)]
    new_loc = pynp.concatenate(([0], new_loc, [len(data)]))

    seg_len = []
    last_cn = []
    new_data = []
    for i in range(len(new_loc)-1):
        s, e = int(new_loc[i]), int(new_loc[i+1])
        tmp_len = e-s
        if s<=e :
            mean_value = pynp.round(pynp.mean(data[s:e]))
            new_data.append(pynp.repeat(mean_value, tmp_len))
            if len(last_cn) == 0:
                last_cn.append(mean_value)
                seg_len.append(tmp_len)
            if mean_value == last_cn[len(last_cn)-1]:
                seg_len[len(seg_len)-1] = seg_len[len(seg_len)-1] + tmp_len
            else:
                last_cn.append(mean_value)
                seg_len.append(tmp_len)
    #
    # # segment
    # segments = np.zeros((len(new_loc)-1, 3))
    # new_data = []
    # for i in range(segments.shape[0]):
    #     s, e = int(new_loc[i]), int(new_loc[i+1])
    #     if s<=e :
    #         segments[i, 0] = s
    #         segments[i, 1] = e
    #         segments[i, 2] = e - s
    #         mean_value = np.round(np.mean(data[s:e]))
    #         # print('>>>>', mean_value,s,e,segments[i, 2], '<<<')
    #         new_data.append(np.repeat(mean_value, segments[i, 2]))
    try:
        if len(seg_len) <=1:
            rearrange_score = 0
        else:
            # breakpoint_count = segments.shape[0]
            breakpoint_prop = len(seg_len) / len(data)
            uniformity_score = -skew(seg_len)
            # sigmoid
            uniformity_score = 1 / (1 + pynp.exp(-uniformity_score))
            avg_amplitude = pynp.var(pynp.concatenate(new_data))
            rearrange_score = breakpoint_prop * uniformity_score * avg_amplitude
    except:
        rearrange_score = 0

    return rearrange_score



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef dict permutation_score(int x, double[:] cnv, double[:, :] gene_loc,
                             double[:] exp_lambda,
                             double[:] possion_lambda, int chr_num, int random_num, bint random):
    cdef dict dc = fast_score(cnv, gene_loc)
    cdef float limit_num
    cdef int num_greater = 0
    cdef int rn, count
    cdef double tmp_lambda
    cdef cnp.ndarray[cnp.double_t, ndim=1]  diff_res
    cdef float permuted_data
    cdef dict permuted_res
    # cdef float r, s

    limit_num = cal_limit(cnv, 0.6)

    dc['limit'] = limit_num #/ len(cnv)
    tmp_lambda = exp_lambda[x]
    tmp_possion_lambda = possion_lambda[chr_num]
    # print('>>>>')
    # print(tmp_possion_lambda)

    if random:
        # np.random.poisson()
        if pynp.isnan(tmp_lambda):
            tmp_lambda = pynp.mean(exp_lambda, where=~pynp.isnan(exp_lambda))
        #len(dc['smooth_segment']) <= 2 or \
        if pynp.isnan(tmp_possion_lambda) or tmp_possion_lambda is None:
            # num_greater = random_num
            tmp_possion_lambda = pynp.mean(possion_lambda, where=~pynp.isnan(possion_lambda))
        if  pynp.isnan(tmp_lambda) or pynp.isnan(tmp_possion_lambda):
            num_greater = random_num
        else:

            for rn in range(1, random_num + 1):
                #ss = time.time()
                permuted_data = permutation_bg(cnv,
                                               gene_loc,
                                               tmp_lambda, tmp_possion_lambda, len(dc['smooth_segment']), seed=rn)
                #e1 = time.time()
                # permuted_res = fast_score(permuted_data, gene_loc)
                #e2 = time.time()
                #print(f'time:{e1-ss}, {e2-e1}')
                # r =  permuted_res['rearrange_score']# np.var(permuted_data)#
                # s =  dc['rearrange_score']# np.var(cnv)# dc['rearrange_score']
                if permuted_data >= dc['rearrange_score']:
                    num_greater += 1

    return {'dc': dc, 'num_greater': num_greater, 'x': x}


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef double cal_R_pyx(double[:] cnv_data):
    cdef int len_cnv = cnv_data.shape[0]
    cdef cnp.ndarray[cnp.double_t, ndim=1] bkps = pynp.where(pynp.diff(cnv_data) != 0)[0].astype('double')
    cdef int num_bkps = bkps.shape[0]

    if num_bkps <= 1:
        return 0.0

    # B
    # cdef double B = <double> num_bkps / len_cnv

    # U
    cdef cnp.ndarray[double, ndim=1] seg_len = pynp.abs(pynp.diff(bkps))
    cdef double U = 1.0 / (1.0 + pynp.exp(-skew(seg_len)))
    # A
    cdef double A = pynp.std(cnv_data)
    # A = 1-1/(1+A)
    A = A / pynp.mean(cnv_data)
    if pynp.isnan(A):
        A = 0
    if A > 1:
        A = 1

    return U * A


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef perm_cnv_pyx(double[:] cnv_data, double exp_lam, double possion_lam,int random_num, int seed):
    cdef int len_x = len(cnv_data)
    pynp.random.seed(seed)
    cdef cnp.ndarray[cnp.double_t, ndim=1] seg_size = pynp.random.poisson(possion_lam, size=random_num).astype('double')
    cdef cnp.ndarray[cnp.double_t, ndim=1] new_seg, new_loc
    cdef int e, s, k, idx, i
    cdef list tmp_new_data, res


    res = []
    for idx,i in enumerate(seg_size):
        # print(idx)
        if i == 0:
            res.append(0)
        else:
            pynp.random.seed(idx)
            new_seg = pynp.random.exponential(exp_lam, size=i+1)
            new_seg /= pynp.sum(new_seg)
            new_loc = pynp.cumsum(pynp.round(new_seg * len_x))
            new_loc = new_loc[new_loc != 0]
            new_loc = new_loc[new_loc < len_x]
            new_loc = pynp.concatenate(([0], new_loc, [len_x]))

            tmp_new_data = []
            for k in range(len(new_loc)-1):
                s, e = int(new_loc[k]), int(new_loc[k+1])
                tmp_new_data.append(pynp.repeat(pynp.mean(cnv_data[s:e]), e-s))
            res.append(cal_R_pyx(pynp.concatenate(tmp_new_data)))
    return res


