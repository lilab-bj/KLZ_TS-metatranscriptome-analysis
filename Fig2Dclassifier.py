import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
import numpy as np, matplotlib.pyplot as plt,  pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import cross_val_score,cross_val_predict,  KFold,  LeaveOneOut, StratifiedKFold
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import plot_roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import OneHotEncoder
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

def dummmize_matrix(df,col):
    """
    dummize category variable
    col:column name to dummized.
    """
    vals = [[i] for i in df.loc[:,col].values]
    enc = OneHotEncoder(handle_unknown='ignore')
    enc.fit(vals)
    enc_df = enc.transform(vals).toarray()
    cols = enc.categories_[0]
    for i in range(len(cols)):
        print(cols[i])
        new_name = col + '_' +str(cols[i])
        new_val = enc_df[:,i]
        df.loc[:,new_name] = new_val
    return df.drop(col,axis=1)

def clean_dataset(df):
    # df needs to be a pd.DataFrame
    assert isinstance(df, pd.DataFrame)
    df.dropna(inplace=True)
    indices_to_keep = ~df.isin([np.nan, np.inf, -np.inf]).any(1)
    return df[indices_to_keep].astype(np.float64)

if __name__ == "__main__":

    ### load dataset
    meta = pd.read_excel("map_A.xlsx",index_col = 0)
    df = pd.read_csv("genus_A.csv",
                    index_col= 0)
    df = df.T

    ### arcsin transform
    scaler = StandardScaler()
    for i in df.columns:
        df.loc[:,i] = np.arcsin(df.loc[:,i])

    ### get subset of dataframe
    #meta = meta.loc[meta['day_2'].isin(['d1','d5'])]
    meta1 = meta.loc[:,'Death']
    meta1 = pd.DataFrame(meta1)
    meta1 = clean_dataset(meta1) ##remove na
    overlap = df.index & meta1.index

    meta = meta.loc[overlap,:]
    df = df.loc[overlap,:]

    ### check consistency
    print(sum(df.index != meta.index))

    ### add more features: Severity(A),Gender_Bin,Treatment_Bin,copy,age
    df.loc[:,'Severity_A'] = meta.loc[:,'Severity_A'].values
    df.loc[:,'Age']  = scaler.fit_transform(meta.Age.values.reshape(-1, 1))
    df.loc[:,'CORTICOSTEROID'] = meta.loc[:,'CORTICOSTEROID'].values
    df.loc[:,'COMORBIDITY'] = meta.loc[:,'COMORBIDITY'].values
    df.loc[:,'ANTB_CLS3_STATUS'] = meta.loc[:,'ANTB_CLS3_STATUS'].values

    ## set dummy variable
    df = dummmize_matrix(df,'Severity_A')
    df = df.drop('Severity_A_3',axis=1) 
    df = df.drop('Severity_A_4',axis=1)


    ##############################################################################
    X = df.values
    Y = meta.Death.values
    n_samples, n_features = X.shape
    random_state = np.random.RandomState(0)
    #############################################################################
    # Classification and ROC analysis

    # Run classifier with cross-validation and plot ROC curves
    cv = StratifiedKFold(n_splits=10)

    clf = LogisticRegression(
                            #random_state=0,
                            penalty='l1',
                            C=0.93, 
                            solver='liblinear',
                            )

    tprs = []
    aucs = []

    mean_fpr = np.linspace(0, 1, 100)

    plt.rc('font', size=20)
    plt.rc('font', family='Arial')
    fig = plt.figure(figsize=[10,8])
    ax = fig.subplots()
    for i, (train, test) in enumerate(cv.split(X, Y)):
        clf.fit(X[train], Y[train])


        fpr, tpr, _ = roc_curve(Y[test],
                                clf.predict_proba(X[test])[:,1],
                                pos_label=1)
        auc_ = roc_auc_score(Y[test],clf.predict_proba(X[test])[:,1])
        print(clf.predict_proba(X[test]).shape)
        print("auc=",auc_)
        interp_tpr = np.interp(mean_fpr, fpr, tpr)

        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(auc_)

    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
            #label='Chance',
            alpha=.8)


    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)

    std_auc = np.std(aucs)


    #00A087FF
    ax.plot(mean_fpr, mean_tpr, color='#3C5488FF',
            label=r'Mean ROC (AUC = %0.3f $\pm$ %0.3f)' % (mean_auc, std_auc),
            lw=2, alpha=.8)


    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='#3C5488FF', alpha=.1,
                # label=r'$\pm$ 1 std. dev.'
                )

    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
        title="")
    #ax.legend(loc="lower right",bbox_to_anchor=(1.7,0))
    ax.legend(loc="lower right",bbox_to_anchor=(1,0))
    ax.grid(linestyle='--')
    plt.xlabel('1-Specificity')
    plt.ylabel('Sensitivity')

    y_pred_proba = clf.predict_proba(X)
    y_pred = clf.predict(X)

    print('confusion_matrix:')
    tn, fp, fn, tp = confusion_matrix(Y, y_pred).ravel()
    print("true positive:",tp)
    print("false positive:",fp)
    print("true neg:",tn)
    print("false neg:",fn)

    print('summary:')
    print(classification_report(Y, y_pred, target_names=['ndeath','death']))


    coef_map = dict(zip(df.columns.values,clf.coef_[0]))
    for i,j in coef_map.items():
        if j != 0:
            print(i,':  ',j)
    plt.show()
    #fig.set_size_inches(9, 8.5)
    #plt.savefig("F:/AAAA/KLZ_revise/withoutCopy/species+metadata.pdf",
                #dpi=600,
                #bbox_inches="tight")