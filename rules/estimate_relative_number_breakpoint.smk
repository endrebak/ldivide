

# rule random_partition:
#     input:
#         "{prefix}/1kg/vector/{population}/{window_size}/{chromosome}.pq"
#     output:
#         "{prefix}/1kg/vector/breakpoints/{population}/{window_size}/{chromosome}.txt"
#     run:
#         s = pd.read_parquet(input[0]).squeeze()
#         a = np.array_split(s, 50)
#         sums = pd.Series([np.sum(a.values) for a in np.array_split(s, 50)])

#         s = pd.DataFrame({"PartitionSum": sums})

#         s.to_csv(output[0], sep="\t", index=False)
