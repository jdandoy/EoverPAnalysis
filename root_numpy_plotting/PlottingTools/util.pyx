import numpy as np
cimport numpy as np

ctypedef np.int_t INTDTYPE_t
ctypedef np.float_t FLOATDTYPE_t

ctypedef fused float_or_int_array:
    np.ndarray[FLOATDTYPE_t, ndim=1]
    np.ndarray[INTDTYPE_t, ndim=1]

cpdef getWeightsFromBins(float_or_int_array variable_in_histogram, float_or_int_array low_edges, float_or_int_array high_edges, np.ndarray[FLOATDTYPE_t, ndim = 1] normalization):
    cdef int maxbin = low_edges.shape[0]
    cdef int tracks = variable_in_histogram.shape[0]
    cdef np.ndarray[FLOATDTYPE_t, ndim =1] weights = np.ones(tracks, dtype = np.float)

    for i in range(maxbin):
        for j in range(tracks):
            if (variable_in_histogram[j] >= low_edges[i]) and (variable_in_histogram[j] < high_edges[i]):
                weights[j] = normalization[i]
    return weights
