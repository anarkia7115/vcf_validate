from sklearn.tree import DecisionTreeClassifier

from load_mat import load

y_train, y_test1, y_test2, X_train, X_test1, X_test2, gene_56_df, identifier_df = load()

# Train
tree_clf = DecisionTreeClassifier(max_depth=3)
tree_clf.fit(X_train, y_train)

print("train accuracy: {}".format(tree_clf.score(X_train, y_train)))
print("test accuracy 1: {}".format(tree_clf.score(X_test1, y_test1)))

y_predicted = tree_clf.predict(X_test2)

from utils import eval_clf, plot_tree

# plot tree
plot_tree(tree_clf, X_train)

# print evaluation
eval_clf(y_predicted, identifier_df, gene_56_df)