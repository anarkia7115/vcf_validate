import pandas as pd
from sklearn.tree import export_graphviz


def read_vcf(vcf_file):
    buffer = []
    header = ""
    with open(vcf_file, 'r') as f:
        for l in f:
            l = l.rstrip()
            if l.startswith("##"):
                continue
            elif l.startswith("#"):
                l = l[1:]
                header = l
            else:
                buffer.append(l)
    vcf_df = pd.DataFrame.from_records([xx.split("\t") for xx in buffer], columns=header.split("\t"))
    return vcf_df


def annotation_vs_standard(vcf_df, gene_56_df):

    joined_df = gene_56_df.merge(vcf_df, how='inner',
                     on=['chr_num', 'position'])

    return joined_df[['label', 'INFO']]


def eval_clf(y_predicted, identifier_df, gene_56_df):
    y_pred_df = pd.DataFrame(y_predicted, columns=['predict'])
    # gene_56_df.merge(compare_df.join(identifier_df), how='left',
    #                  on=['chr_num', 'position']).to_csv("compare_tree.csv")
    joined_df = gene_56_df.merge(y_pred_df.join(identifier_df), how='inner',
                     on=['chr_num', 'position'])

    df_size = joined_df.count()[0]
    predict_sum = joined_df.predict.sum()
    # compare_df.join(identifier_df).to_csv("compare_tree2.csv")

    print("test accuracy 2: {}".format(predict_sum / df_size))

    print("Done")


def plot_tree(tree_clf, X_train):
    # Plot tree
    export_graphviz(
        tree_clf,
        out_file="vcf_tree.dot",
        feature_names=X_train.columns,
        class_names=['False', 'True'],
        rounded=True,
        filled=True
    )