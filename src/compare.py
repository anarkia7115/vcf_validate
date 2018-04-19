chr_pos_file = "/home/shawn/PycharmProjects/vcf_validate/data/chr_pos_to_join.csv"
import pandas as pd

chr_pos_df = pd.read_csv(chr_pos_file, sep=',', header=-1, names=['chr_num', 'position'])

identifier_df = pd.read_pickle("/home/shawn/PycharmProjects/vcf_validate/data/test_identifier_04_02.pkl")

joined_df = chr_pos_df.merge(
    identifier_df,
    how="inner", suffixes=['_chr', '_iden'], left_on=['chr_num', 'position'],
    right_on=['chr_num', 'position'])

pass
