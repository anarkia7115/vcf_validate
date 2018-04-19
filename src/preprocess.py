import pandas as pd


def flatten_features(csv_file):
    csv_col_names = ['chr_num', 'position', 'orig', 'alter', 'label', 'gatk_label',
                 'combined_attrs1', 'attr2_names', 'attr2_values']
    vcf_df = pd.read_csv(csv_file, header=-1, names=csv_col_names, dtype={'label': bool})

    # get combined attribute columns
    combined_attrs1_col = vcf_df.combined_attrs1
    attr2_names_col = vcf_df.attr2_names
    attr2_values_col = vcf_df.attr2_values


def main():
    vcf_file = "/home/shawn/PycharmProjects/vcf_validate/data/B1700_training.csv"
    flatten_features(vcf_file)

if __name__ == "__main__":
    main()