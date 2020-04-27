#limit threads to 1 for cluster
import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"



import json
import random
from sklearn import decomposition, svm
import numpy as np
import sys
from sklearn.feature_selection import SelectPercentile, f_classif
from sklearn.metrics import roc_auc_score
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import Imputer, StandardScaler


# arguments: path to data base json, percentage of training set size, booleans which parts of the data base should be used for predictions
DATA_DIR = sys.argv[1]
TRAINING_P = float(sys.argv[2])
KMER = bool(int(sys.argv[3]))
MDS = bool(int(sys.argv[4]))
C_FEATURES = bool(int(sys.argv[5]))
D_FEATURES = bool(int(sys.argv[6]))
RUNS = int(sys.argv[7])
PCA = False


# count number of discrete states
# input: database
# output: list of discrete states in database
def count_states(data):
    states_list = []
    for i in data:
        if i != 'kmer_list':
            for j in range(len(data[i])):
                index = 0
                for feature in data[i][j]['d_features']:
                    if len(states_list) == index:
                        states_list.append(0)
                    if feature > states_list[index]:
                        states_list[index] = int(feature)
                    index += 1
    return states_list


# convert discrete states to binary states
# input: discrete features, list of discrete states
# output: binary features
def binary(features, state_list):
    bin_features = []
    for i in range(len(features)):
        for j in range(state_list[i]):
            if features[i] == j+1:
                bin_features.append(1)
            else:bin_features.append(0)
    return bin_features




# get data for machine learning
# input: all data, species name, individual number
# output: selected data
def get_features(data, key, i):
    feature_data = []
    if C_FEATURES:
        feature_data.extend(data[key][i]['c_features'])
    if D_FEATURES:
        feature_data.extend(binary(data[key][i]['d_features'], count_states(data)))
    if MDS:
        feature_data.extend((data[key][i]['mds_features']))
    if KMER:
        for kmer in data['kmer_list']:
            if kmer in data[key][i]['kmer_profile']:
                feature_data.append(data[key][i]['kmer_profile'][kmer])
            else:
                feature_data.append(0)
    return feature_data



# compare two list if their items are the same
# input: to be compared lists
# output: number of identical items
def compare(list1, list2):
    i = 0
    match = 0
    for item in list1:
        if item == list2[i]:
            match += 1
        i += 1
    return match


# compares two lists and returns a list with 1 on positions where they are identical and 0 where they different
def evaluate(list1, list2):
    result = []
    for i in range(len(list1)):
        if list1[i] == list2[i]:
            result.append(1)
        else:
            result.append(0)
    return result


# get the highest propabilitys for each prediction
# input: list of propabiltys
# output: list of highest propabilitys
def get_selected_propa(propa_list):
    result = []
    for entry in propa_list:
        result.append(max(entry))
    return result


# read data base file
with open(DATA_DIR) as json_data:
    data = json.load(json_data)
    json_data.close()

pre_features = []
features = []
category = []
features_test = []
test_category = []



# crossvalidate to find the best percentile (p) for feature selection
p_min = 0
best_score = 0

for p in range(1, 101):
    svm_score = 0
    nn_score = 0
    svm_p = 0
    nn_p = 0
    for run in range(RUNS):
        pre_features = []
        features_test = []
        test_category = []
        category = []
        features =[]

        # collect the data needed for svm and nn
        for key in data:
            i = 0
            if key != 'kmer_list':
                random.shuffle(data[key])
                while i < len(data[key]):
                    training_part = round(len(data[key]) * TRAINING_P, 0)
                    if training_part == 0:
                        training_part = 1
                    if i < len(data[key]) - training_part:
                        pre_features.append(get_features(data, key, i))
                        category.append(key)

                    else:
                        features_test.append(get_features(data, key, i))
                        test_category.append(key)

                    i += 1


                features.extend(pre_features)
                pre_features = []


        # transform data to numpy arrays in order to feed them to the predictors
        x = np.array(features)
        y = np.array(category)
        t = np.array(features_test)




        if PCA and KMER:
            gen_data = x[...,len(x[0])-len(data['kmer_list'])+1:]
            gen_data2 = t[...,len(t[0])-len(data['kmer_list'])+1:]
            gen_data = np.append(gen_data, gen_data2, axis=0)
            pca = decomposition.PCA(n_components=10)
            pca.fit(gen_data)
            #print(pca.explained_variance_ratio_)
            gen_data = pca.fit_transform(gen_data)

            x = np.append(x[...,:len(x[0])-len(data['kmer_list'])+1], gen_data[:-len(gen_data2)], axis=1)
            t = np.append(t[...,:len(t[0])-len(data['kmer_list'])+1], gen_data[-len(gen_data2):], axis=1)


        #scale data
        scaler = StandardScaler()

        x = scaler.fit_transform(np.array(x))
        y = np.array(category)
        t = scaler.transform(np.array(t))


        # select features
        select_percentile = SelectPercentile(score_func=f_classif, percentile=p)
        x = select_percentile.fit_transform(x, y)
        t = select_percentile.transform(t)




        # implement predictors
        svc = svm.LinearSVC(C=1, penalty="l1", dual=False).fit(x, y)
        mlp = MLPClassifier(hidden_layer_sizes=(100, 100, 100, 100), max_iter=400).fit(x, y)

        # run predictions
        svm_res = svc.predict(t)
        svm_poba = get_selected_propa(svc.decision_function(t))



        # compare the predicted svm result with the expected results and calculate the auc score using the propability of prediction
        svm_res_eva = evaluate(test_category, svm_res)
        is_not_perfect = False
        is_not_all_wrong = False
        for value in svm_res_eva:
            if value == 0:
                is_not_perfect = True
            else:
                is_not_all_wrong = True
        if is_not_perfect and is_not_all_wrong:
            svm_score += roc_auc_score(svm_res_eva, svm_poba)
        elif is_not_all_wrong:
            svm_score += 1
        elif is_not_perfect:
            svm_score += 0.5

        svm_p += compare(svm_res, test_category)

        # compare the predicted neural network result with the expected results and calculate the auc score using the propability of prediction
        nn_res = mlp.predict(t)
        nn_poba = get_selected_propa(mlp.predict_proba(t))
        nn_res_eva = np.array(evaluate(test_category, nn_res))

        is_not_perfect = False
        is_not_all_wrong = False
        for value in nn_res_eva:
            if value == 0:
                is_not_perfect = True
            else:
                is_not_all_wrong = True
        if is_not_perfect and is_not_all_wrong:
            nn_score += roc_auc_score(nn_res_eva, nn_poba)
        elif is_not_all_wrong:
            nn_score += 1
        elif is_not_perfect:
            nn_score += 0.5
        nn_p += compare(nn_res, test_category)

    svm_score /= RUNS
    nn_score /= RUNS
    svm_p /= RUNS
    nn_p /= RUNS
    if svm_p > best_score:
        p_min = p
        best_score = svm_p
    if svm_p == best_score:
        p_max = p

# write prediction results to json file
results = [abs(0.5 - svm_score), abs(0.5 - nn_score), str(best_score) + '/' + str(len(test_category)),
           str(nn_p) + '/' + str(len(test_category)), float(best_score) / len(test_category), float(nn_p) / len(test_category), p_min, p_max]
with open('svm.json', 'w') as file:
    json.dump(results, file)

