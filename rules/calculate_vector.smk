from ldivide.src.calc_autocovar import calc_autocovar 

rule theta:
    input:
        "{prefix}/{population}.txt"
    output:
        "{prefix}/theta/{population}.txt"
    run:
        inds = pd.read_table(indfile, header=None, squeeze=True).to_list()

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
        haplos = "{prefix}/1kg/{population}/{chromosome}.vcf.gz",
        theta = "{prefix}/theta/{population}.txt"
    output:
        "{prefix}/1kg/autocorr/{population}/{chromosome}.vcf.gz"
    run:

        theta = float(open(input.theta).readline().strip())

        autocovars = []
        for chunk, cdf in pd.read_table(input[0], sep="\t", chunksize=50000):
            calc_autocovar(cdf, theta, ) 














