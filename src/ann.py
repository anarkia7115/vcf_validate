import tensorflow as tf
from load_mat import load
from utils import eval_clf

# import os
# os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

y_train, y_test1, y_test2, X_train, X_test1, X_test2, gene_56_df, identifier_df = load()

feature_columns = tf.contrib.learn.infer_real_valued_columns_from_input(X_train)
dnn_clf = tf.contrib.learn.DNNClassifier(hidden_units=[150, 150], n_classes=2,
                                         feature_columns=feature_columns, model_dir="vcf.model")

# dnn_clf.fit(x=X_train, y=y_train, batch_size=50, steps=40000)

y_predicted = dnn_clf.predict_classes(X_test2)


print("train accuracy: {}".format(dnn_clf.evaluate(X_train, y_train)['accuracy']))
print("test accuracy 1: {}".format(dnn_clf.evaluate(X_test1, y_test1)['accuracy']))
eval_clf(y_predicted, identifier_df, gene_56_df)
