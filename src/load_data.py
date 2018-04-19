import pandas as pd
import numpy as np
from sklearn import preprocessing

col_names = ['chr_num', 'position', 'orig', 'alter', 'label', 'gatk_label',
             'combined_attrs1', 'attr2_names', 'attr2_values']

# MODE = "train"
SCALE = False
MODE = "test"

if MODE == "test":
    # vcf_file = "/home/shawn/PycharmProjects/vcf_validate/data/B1702_training.csv"
    vcf_file = "/home/shawn/PycharmProjects/vcf_validate/data/" \
               "seq_raw_vcf_file-20180402_biaozhunping_AE1800182FFP_56gene_raw.csv"
    vcf_df = pd.read_csv(vcf_file, sep='\t', header=-1, names=col_names, dtype={'label': bool})
else:
    vcf_file = "/home/shawn/PycharmProjects/vcf_validate/data/B1700_training.csv"
    vcf_df = pd.read_csv(vcf_file, header=-1, names=col_names, dtype={'label': bool})

# chr_indicator_col = pd.get_dummies(vcf_df.chr_num.astype('category'))

# scaler
# min_max_scaler = preprocessing.MinMaxScaler()
# scaled_positions = min_max_scaler.fit_transform(
#  vcf_df.position.astype('float').values.reshape(-1, 1)).flatten()
position_col = vcf_df.position.astype('float')


def gene_base_count(row):
    # init empty base count map
    return_row = {base: 0 for base in set(row)}
    for base in row:
        # count base in rows
        try:
            return_row[base] += 1
        except IndexError:
            print(base)
            raise
    return return_row


def gene_base_count_row_func(row):
    return_row = gene_base_count(row)
    return pd.Series(return_row)


orig_cols = vcf_df.orig.apply(gene_base_count_row_func) \
    .rename(columns=lambda x: 'orig_' + x).fillna(0)
alter_cols = vcf_df.alter.apply(gene_base_count_row_func) \
    .rename(columns=lambda x: 'alter_' + x).fillna(0)


# gatk_label_cols = pd.get_dummies(vcf_df.gatk_label.astype('category'))\
#   .rename(columns=lambda x: 'gatk_label_' + x.replace(';', '_'))

# deal with combined attrs (vcf_df.combined_attrs1)
def split_attrs1(row):
    # init return_row
    return_row = dict()
    # split row string
    for attr in row.split(';'):
        if '=' in attr:  # if key value mode
            [k, v] = attr.split('=')  # get key, value
            # Concider 'RPA', and 'RU' as key
            if k == 'RPA':  # RPA: int,int
                v1, v2 = v.split(',')
                return_row['RPA_1'] = np.int32(v1)
                return_row['RPA_2'] = np.int32(v2)
            elif k == 'RU':  # RU: ATCG
                gb_count = gene_base_count(v)
                for gb in gb_count.copy():  # gb will be updated, make a copy
                    new_attr_name = 'RU_' + gb  # new attr name
                    gb_count[new_attr_name] = gb_count.pop(gb)  # update key name
                return_row.update(gb_count)
            else:  # float values
                if v == '.':  # missing value
                    return_row[k] = 0  # use fill na after
                else:  # has float value
                    return_row[k] = np.float64(v)
        else:  # boolean mode
            return_row[attr] = True  # infer False after
    return pd.Series(return_row)


attrs1_cols = vcf_df.combined_attrs1.apply(split_attrs1)


# DB: True, False
# RPA: int,int
# RU: ATCG
# STR: True, False
# Other: Float


def split_attrs2(row):
    attr_names = row['attr2_names'].split(':')
    attr_values = row['attr2_values'].split(':')

    result_rows = dict()

    for name, value in zip(attr_names, attr_values):
        if name == 'AD' or name == 'QSS':  # split by `,`
            v1, v2 = value.split(',')
            result_rows['{}_1'.format(name)] = np.float64(v1)
            result_rows['{}_2'.format(name)] = np.float64(v2)
        elif name == 'GT' or name == 'PGT':  # if exists, true
            result_rows[name] = True
        elif name == 'PID':  # drop
            pass
        else:  # float
            if value == '.':  # missing value
                result_rows[name] = 0  # use fill na after
            else:  # has float value
                result_rows[name] = np.float64(value)

    return pd.Series(result_rows)


attrs2_cols = vcf_df[['attr2_names', 'attr2_values']].apply(split_attrs2, axis=1)
# AD: int,int
# QSS: int,int
# GT: 0/1 -> True, False
# PGT: 0|1 -> True, False
# PID: Drop
attrs2_cols = attrs2_cols.apply(
    lambda x: x.fillna(0) if x.dtype.kind == 'f' else x.fillna(False))

attrs1_cols = attrs1_cols.apply(
    lambda x: x.fillna(0) if x.dtype.kind == 'f' else x.fillna(False))

scaler = preprocessing.MinMaxScaler()
# new_df = pd.concat([vcf_df.label, position_col, orig_cols, alter_cols, gatk_label_cols,
#                     attrs1_cols, attrs2_cols, chr_indicator_col], axis=1)

# new_df = pd.concat([vcf_df.label, gatk_label_cols, attrs1_cols, attrs2_cols], axis=1)
new_df = pd.concat([vcf_df.label, attrs1_cols, attrs2_cols], axis=1)

# selected_columns = ['label', 'gatk_label_PASS', 'gatk_label_clustered_events',
#                     'gatk_label_clustered_events_homologous_mapping_event', 'DB', 'ECNT',
#                     'HCNT', 'MAX_ED', 'MIN_ED', 'RPA_1', 'RPA_2', 'RU_A', 'RU_C', 'RU_G',
#                     'RU_T', 'STR', 'TLOD', 'AD_1', 'AD_2', 'AF', 'ALT_F1R2', 'ALT_F2R1',
#                     'FOXOG', 'GT', 'PGT', 'QSS_1', 'QSS_2', 'REF_F1R2', 'REF_F2R1']

selected_columns = ['label', 'DB', 'ECNT',
                    'HCNT', 'MAX_ED', 'MIN_ED', 'RPA_1', 'RPA_2', 'RU_A', 'RU_C', 'RU_G',
                    'RU_T', 'STR', 'TLOD', 'AD_1', 'AD_2', 'AF', 'ALT_F1R2', 'ALT_F2R1',
                    'FOXOG', 'GT', 'PGT', 'QSS_1', 'QSS_2', 'REF_F1R2', 'REF_F2R1']

# add empty columns
for name in selected_columns:
    if name not in new_df.columns:
        new_df[name] = np.nan

# normalize data
if SCALE:
    new_df = scaler.fit_transform(new_df[selected_columns])
else:
    # fill na with 0
    new_df = new_df[selected_columns].fillna(0)

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

pd.set_option("display.max_columns", 100)
pd.set_option("display.max_rows", 100)

if MODE == "test":
    new_df[selected_columns].to_pickle("/home/shawn/PycharmProjects/vcf_validate/data/test_04_02.pkl")
    vcf_df[['chr_num', 'position']].to_pickle("/home/shawn/PycharmProjects/vcf_validate/data/test_identifier_04_02.pkl")
else:
    train, test = train_test_split(new_df[selected_columns], test_size=0.2)
    train.to_pickle("../data/train.pkl")
    test.to_pickle("../data/test.pkl")
