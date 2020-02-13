import numpy as np
from ldivide.src.calc_autocovar import calc_autocovar
from ldivide.src.calc_covar import calc_covar

rule theta:
    input:
        "{prefix}/partitions/{population}/{chromosome}/samples.txt"
    output:
        "{prefix}/theta/{chromosome}/{population}.txt"
    run:
        inds = pd.read_table(input[0], header=None).squeeze().to_list()

        haps= list()
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

        # autocovars = []
        # df = pd.read_table(input.haplos, sep="\t", dtype=np.int8, header=None).values.copy(order="C")
        df = pd.read_parquet(input.haplos).values
        # print(df.head())
        print("Done reading!")
        outvec = calc_autocovar(df, theta)

        # outvec = np.concatenate(autocovars)

        df = pd.Series(outvec).to_frame()
        df.columns = df.columns.astype(str)
        df.to_parquet(output[0])



rule vector:
    input:
        autocorr = "{prefix}/1kg/autocorr/{population}/{chromosome}.pq",
        haplos = "{prefix}/1kg/{population}/{chromosome}.pq",
        theta = "{prefix}/theta/{chromosome}/{population}.txt"
    output:
        "{prefix}/1kg/vector/{population}/{chromosome}.pq"
    benchmark:
        "{prefix}/benchmark/1kg/vector/{population}/{chromosome}.tsv.gz"
    run:
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
















