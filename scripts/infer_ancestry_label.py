import os
from run_vcf_plink_filter import HERE, INPUT, OUTPUT
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
# Using scikit-learn to split data into training and testing sets
from sklearn.model_selection import train_test_split
# Import Random Forest Model
from sklearn.ensemble import RandomForestClassifier
# Import scikit-learn metrics module for accuracy calculation
from sklearn import metrics

# create figure folder if it doesn't exist already
if not os.path.isdir(os.path.join(OUTPUT, 'figures')):
    os.makedirs(os.path.join(OUTPUT, 'figures'))

# create result folder if it doesn't exist already
if not os.path.isdir(os.path.join(OUTPUT, 'result')):
    os.makedirs(os.path.join(OUTPUT, 'result'))


def clean_pca_format(eigenvec, label):
    """

	:param eigenvec: eigenvector file
	:param label: ancestry label file
	:return: file with proper format
	"""
    pca = pd.read_table(eigenvec)
    # rename colun #IID to sample_id to match with the ancestry label file
    pca.rename(columns={'#IID': 'sample_id'}, inplace=True)
    # print(pca.head())
    label = pd.read_table(label)
    # join data frame by sample_id to get the ancestry label
    pca = pd.merge(label, pca, on='sample_id')
    # pca.ancestry = pca.ancestry.astype('category')
    return pca


def calc_variance_explained(val):
    """

	:param val: eigenvalues file
	:return: variance explained by each PC
	"""
    eigenval = pd.read_table(val, names=['Value'])
    per_exp = pd.DataFrame([i for i in range(1, 10 + 1)], columns=['PC'])
    Per_Val = round(eigenval / np.sum(eigenval) * 100, 2)
    per_exp['Value'] = Per_Val
    return per_exp


def elbow_plot(data):
    """

	:param data: variance explained by each PC
	:return: elbow plot
	"""
    fig, ax = plt.subplots(figsize=(5, 2.7))
    ax.plot('PC', 'Value', data=data, marker='o', ls='--')
    ax.set_xlabel('PC')
    ax.set_ylabel('Explained variance')
    # plt.show()
    plt.savefig(os.path.join(OUTPUT, 'figures', 'elbow_plot.png'))


def pca_plot(data, data_, data_anc, per):
    """

	:param data: PC1 column
	:param data_: PC2 column
	:param data_anc: ancestry label column
	:param per: percentage of variance explained
	:return: pca plot
	"""
    fig, ax = plt.subplots(figsize=(5, 2.7))
    ancestrys = {'afr': 'tab:blue', 'amr': 'tab:orange', 'eas': 'tab:green', 'fin': 'tab:red', 'nfe': 'tab:purple',
                 'sas': 'tab:pink',
                 'none': 'tab:gray'}
    plt.scatter(data, data_, c=data_anc.map(ancestrys))
    handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=v, label=k, markersize=8) for k, v in
               ancestrys.items()]
    plt.legend(title='Ancestry', handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.hlines(0.00, -0.05, 0.06, ls='dotted', alpha=0.3, color='black')
    plt.xlim(-0.05, 0.06)
    plt.vlines(0.00, -0.05, 0.06, ls='dotted', alpha=0.3, color='black')
    plt.ylim(-0.05, 0.06)
    plt.xlabel('PC1 ' + '(' + str(per.Value[0]) + ')')
    plt.ylabel('PC2 ' + '(' + str(per.Value[1]) + ')')
    plt.grid(color='lightblue', ls='--', alpha=0.4)
    plt.title('PCA plot of human ancestry with missing labels', fontsize=20)
    ax.spines[:].set_visible(False)
    # plt.show()
    plt.savefig(os.path.join(OUTPUT, 'figures', 'pca_plot_missing_labels.png'))


def rf_implementation(data, data_two):
    """

	:param data: pca data with ancestry label
	:param data_two: pca data without ancestry label
	:return: predicts ancestry label for individuals missing labels
	"""
    # Labels are the values we want to predict
    labels = np.array(data['ancestry'])

    # Remove the labels and sample id from the features
    # axis 1 refers to the columns
    features = data.drop(['sample_id', 'ancestry'], axis=1)

    # Convert to numpy array
    features = np.array(features)

    # Split the data into training and testing sets
    train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size=0.20,
                                                                                random_state=88)
    # create the model with 1000 decision trees
    rf = RandomForestClassifier(n_estimators=1000, random_state=88)

    # Train the model on training set
    rf.fit(train_features, train_labels)

    # Predict on the test set
    pred = rf.predict(test_features)
    # Model Accuracy
    # print("Accuracy: ", metrics.accuracy_score(test_labels,pred))

    # Remove the labels and sample id from the features
    # axis 1 refers to the columns
    features_two = data_two.drop(['sample_id', 'ancestry'], axis=1)
    features_two_np = np.array(features_two)
    pred_label = rf.predict(features_two_np)
    pred_label_df = pd.DataFrame(pred_label, columns=['ancestry'])
    features_two.insert(0, 'sample_id', data_two['sample_id'], True)
    features_two.insert(1, 'ancestry', pred_label_df['ancestry'].values, True)
    features_two['label'] = 'Inferred'
    return features_two


def combine_original_predicted_dataframe(data, data_two):
    """

	:param data: data with predicted ancestry labels
	:param data_two: data with original ancestry labels
	:return: merged dataframe
	"""
    full_data = pd.concat([data, data_two]).reset_index(drop=True)
    return full_data


def pca_plot_with_full_label(data, data_, data_anc, per):
    """

	:param data: PC1 column
	:param data_: PC2 column
	:param data_anc: ancestry label column
	:param per: percentage of variance explained
	:return: pca plot
	"""
    fig, ax = plt.subplots(figsize=(5, 2.7))
    ancestrys = {'afr': 'tab:blue', 'amr': 'tab:orange', 'eas': 'tab:green', 'fin': 'tab:red', 'nfe': 'tab:purple',
                 'sas': 'tab:pink'}
    plt.scatter(data, data_, c=data_anc.map(ancestrys))
    handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=v, label=k, markersize=8) for k, v in
               ancestrys.items()]

    plt.legend(title='Ancestry', handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.hlines(0.00, -0.05, 0.06, ls='dotted', alpha=0.3, color='black')
    plt.xlim(-0.05, 0.06)
    plt.vlines(0.00, -0.05, 0.06, ls='dotted', alpha=0.3, color='black')
    plt.ylim(-0.05, 0.06)
    plt.xlabel('PC1 ' + '(' + str(per.Value[0]) + ')')
    plt.ylabel('PC2 ' + '(' + str(per.Value[1]) + ')')
    plt.grid(color='lightblue', ls='--', alpha=0.4)
    plt.title('PCA plot of human ancestry with inferred labels', fontsize=20)
    ax.spines[:].set_visible(False)
    # plt.show()
    plt.savefig(os.path.join(OUTPUT, 'figures', 'pca_plot_with_full_labels.png'))


if __name__ == '__main__':
    EIGENVEC = [os.path.join(OUTPUT, 'filter', i) for i in os.listdir(os.path.join(OUTPUT, 'filter'))
                if i.endswith('eigenvec')][0]
    EIGENVAL = [os.path.join(OUTPUT, 'filter', i) for i in os.listdir(os.path.join(OUTPUT, 'filter'))
                if i.endswith('eigenval')][0]
    LABEL = [os.path.join(INPUT, i) for i in os.listdir(os.path.join(INPUT))
             if i.endswith('txt')][0]
    pca_df = clean_pca_format(EIGENVEC, LABEL)
    pca_df.ancestry = pca_df.ancestry.fillna('none')
    pca_df.ancestry = pca_df.ancestry.astype('category')
    pca_df_original = clean_pca_format(EIGENVEC, LABEL)
    pca_no_label = pca_df_original[pca_df_original.ancestry.isna()]
    pca_label = pca_df_original.dropna()
    percentage_exp = calc_variance_explained(EIGENVAL)
    elbow_plot(percentage_exp)
    pca_plot(pca_df.PC1, pca_df.PC2, pca_df.ancestry, percentage_exp)
    original_label = pca_label.copy()
    original_label['label'] = 'Known'
    pred_label = rf_implementation(pca_label, pca_no_label)
    full_data = combine_original_predicted_dataframe(pred_label, original_label)
    full_data.to_csv(os.path.join(OUTPUT, 'result', 'inferred_ancestry_label.csv'))
    pca_plot_with_full_label(full_data.PC1, full_data.PC2, full_data.ancestry, percentage_exp)
