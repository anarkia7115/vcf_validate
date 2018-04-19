import pandas as pd


def load():
    train = pd.read_pickle("/home/shawn/PycharmProjects/vcf_validate/data/train.pkl")
    test1  = pd.read_pickle("/home/shawn/PycharmProjects/vcf_validate/data/test.pkl")
    test2 = pd.read_pickle("/home/shawn/PycharmProjects/vcf_validate/data/test_04_02.pkl")

    identifier_df = pd.read_pickle("/home/shawn/PycharmProjects/vcf_validate/data/test_identifier_04_02.pkl")

    gene_56_file = "/home/shawn/PycharmProjects/vcf_validate/data/56_gene.csv"
    gene_56_df = pd.read_csv(gene_56_file, header=-1, names=['chr_num', 'position'])


    y_train = train.label.astype(int)
    y_test1 = test1.label.astype(int)
    y_test2 = test2.label.astype(int)

    X_train = train.drop('label', axis=1)
    X_test1 = test1.drop('label', axis=1)
    X_test2 = test2.drop('label', axis=1)

    return y_train, y_test1, y_test2, X_train, X_test1, X_test2, gene_56_df, identifier_df