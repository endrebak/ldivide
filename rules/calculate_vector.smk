import numpy as np
from ldivide.src.calc_autocovar import calc_autocovar
from ldivide.src.calc_covar import calc_covar
from ldivide.src.calc_square import calc_rowvec, calc_colvec, calc_final

rule theta:
    input:
        "{prefix}/partitions/{population}/{chromosome}/samples.txt"
    output:
        "{prefix}/theta/{chromosome}/{population}.txt"
    run:
        inds = pd.read_table(input[0], header=None).squeeze().to_list()

        haps = list()
        nind_int = len(inds)
        s = 0

        for _i in range(1, 2*nind_int):
            s = s+ 1.0/float(_i)

        s = 1/s

        theta = s/(2.0*float(nind_int)+s)

        with open(output[0], "w+") as out_handle:
            out_handle.write(str(theta) + "\n")


rule covar_self_vs_self:
    input:
        haplos = "{prefix}/1kg/{population}/{chromosome}.pq",
        theta = "{prefix}/theta/{chromosome}/{population}.txt"
    output:
        "{prefix}/1kg/autocorr/{population}/{chromosome}.pq"
    benchmark:
        "{prefix}/benchmark/1kg/autocorr/{population}/{chromosome}.tsv.gz"
    run:
        theta = float(open(input.theta).readline().strip())
        print(theta)

        df = pd.read_parquet(input.haplos).values
        outvec = calc_autocovar(df, theta)

        df = pd.Series(outvec).to_frame()
        df.columns = df.columns.astype(str)
        df.to_parquet(output[0])



rule vector:
    input:
        autocorr = "{prefix}/1kg/autocorr/{population}/{chromosome}.pq",
        haplos = "{prefix}/1kg/{population}/{chromosome}.pq",
        theta = "{prefix}/theta/{chromosome}/{population}.txt"
    output:
        "{prefix}/1kg/vector/{population}/{window_size}/{chromosome}.pq"
    benchmark:
        "{prefix}/benchmark/1kg/vector/{population}/{window_size}/{chromosome}.tsv.gz"
    run:
        window_size = int(wildcards.window_size)
        theta = float(open(input.theta).readline().strip())

        nrows = None
        print("Reading haplos")
        df = pd.read_parquet(input.haplos).values.copy(order="C")
        print("read haplos")
        autocovars = pd.read_parquet(input.autocorr).squeeze()
        print("read autocovars")
        # print(autocovars)
        autocovars = autocovars.values
        # print(autocovars)
        # raise
        print(df)
        print(autocovars)

        outvec = calc_covar(df, autocovars, theta, window_size)

        df = pd.DataFrame({"Covars": outvec})
        df.columns = df.columns.astype(str)

        df.to_parquet(output[0])


rule compute_row:
    input:
        haplos = "{prefix}/1kg/{population}/{chromosome}.pq",
        theta = "{prefix}/theta/{chromosome}/{population}.txt",
        autocorr = "{prefix}/1kg/autocorr/{population}/{chromosome}.pq"
    output:
        "{prefix}/1kg/rowvectors/{population}/{window_size}/{chromosome}.pq"
    run:
        window_size = int(wildcards.window_size)

        theta = float(open(input.theta).readline().strip())
        haplos = pd.read_parquet(input.haplos)
        autocovars = pd.read_parquet(input.autocorr).squeeze().values

        result = calc_rowvec(haplos.values, autocovars, theta, window_size)
        result = pd.Series(result)
        result = result.to_frame()
        result.columns = result.columns.astype(str)
        result.to_parquet(output[0])


rule compute_col:
    input:
        haplos = "{prefix}/1kg/{population}/{chromosome}.pq",
        theta = "{prefix}/theta/{chromosome}/{population}.txt",
        autocorr = "{prefix}/1kg/autocorr/{population}/{chromosome}.pq"
    output:
        "{prefix}/1kg/colvectors/{population}/{window_size}/{chromosome}.pq"
    run:
        window_size = int(wildcards.window_size)
        theta = float(open(input.theta).readline().strip())
        haplos = pd.read_parquet(input.haplos)
        autocovars = pd.read_parquet(input.autocorr).squeeze().values

        result = calc_colvec(haplos.values, autocovars, theta, window_size)
        result = pd.Series(result)
        print(result)
        result = result.to_frame()
        result.columns = result.columns.astype(str)
        print(result)
        result.to_parquet(output[0])


rule compute_square:
    input:
        col = "{prefix}/1kg/colvectors/{population}/{window_size}/{chromosome}.pq",
        row = "{prefix}/1kg/rowvectors/{population}/{window_size}/{chromosome}.pq",
    output:
        "{prefix}/1kg/square/{population}/{window_size}/{chromosome}.pq"
    run:

        window_size = int(wildcards.window_size)
        col = pd.read_parquet(input.col).squeeze()
        row = pd.read_parquet(input.row).squeeze()

        print(row.isnull().sum())

        final = calc_final(row.values, col.values, window_size)

        df = pd.Series(final).to_frame()
        print(df)
        df.columns = df.columns.astype(str)
        df.to_parquet(output[0])


rule normalize_square:
    input:
        "{prefix}/1kg/square/{population}/{window_size}/{chromosome}.pq"
    output:
        "{prefix}/1kg/normalized_square/{population}/{window_size}/{chromosome}.pq"
    run:
        s = pd.read_parquet(input[0]).squeeze()

        window_size = int(wildcards.window_size)

        half_window_size = (window_size // 2)

        l = len(s)
        n = np.ones(l) * (((half_window_size - 1) * (half_window_size - 2))/2)
        flank_n = (np.arange(2, half_window_size + 1) * (np.arange(2, half_window_size + 1) - 1))/2
        n[:half_window_size - 1] = flank_n
        n[len(n)-(half_window_size-1):] = flank_n[::-1]

        print(s)
        s = s / n
        print(s)
        df = s.to_frame()
        df.columns = df.columns.astype(str)
        df.to_parquet(output[0])















