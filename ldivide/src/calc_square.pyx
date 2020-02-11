
cpdef calc_covar(haplos, double [::1] autocovars, int32_t [::1] mids_between_breakpoints, double theta):

    cdef:
        int mid_x, mid1, mid2, i, j, k
        double rowsum
        int8_t[:, :] haps_view

    haps_view = haplos
    thetas = (1-theta)*(1-theta)
    colvec = np.zeros(len_haps)
    colvec_view = colvec

    rowvec = np.zeros(len_haps)
    rowvec_view = rowvec

    for mid_x in range(len(mids_between_breakpoints) - 1):

        mid1 = mids_between_breakpoints[mid_x]
        mid2 = mids_between_breakpoints[mid_x + 1]
        rowsum = 0
        ai = autocovars[i]
        aj = autocovars[j]

        for i in range(mid1, mid2):
            for j in range(i + 1, mid2):

                n11, n01, n10 = 0, 0, 0
                for k in range(len_g1_int):
                    if haps_view[i][k] == 1 and haps_view[j][k] == 1:
                        n11 += 1
                    elif haps_view[i][k] == 0 and haps_view[j][k] == 1:
                        n01 += 1
                    elif haps_view[i][k] == 1 and haps_view[j][k] == 0:
                        n10 += 1

                f11 = n11/len_g1
                f1 = (n11+n10)/len_g1
                f2 = (n11+n01)/len_g1
                D = f11 - f1*f2
                Ds2 = (thetas*D)
                Ds2 = Ds2 * Ds2
                rsq = Ds2 / (aj1 * aj2)

                rowvec_view[i] += rsq
                colvec_view[j] += rsq

            if i > 1:
                rowvec_view[i] = rowvec_view[i] + rowvec_view[i-1] - colvec_view[i-1]

        for i in range(mid1, mid2):
            if i > 0:
                rowvec_view[i] = rowvec_view[i] / (i * (len_haps - i))
            else:
                rowvec_view[i] = rowvec_view[i] / len_haps


