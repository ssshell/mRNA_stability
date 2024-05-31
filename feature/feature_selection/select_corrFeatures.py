
####### select features by correlation and significance level on target class
####### correlation includes Spearman correlation among features and Kendall correlation between feature and target

import sys

####### load correlation & significance level
feature_list = []
featureCorr_pair = []
transcript_region = {}
N_sigTest = {}
N_featureCorr = {}
featureCorr_pair_rs = {}
targetCorr_tau = {}
avg_sigTest = {}

with open(sys.argv[1]) as f_featureCorrSig:
    next(f_featureCorrSig)
    for line in f_featureCorrSig:
        line = line.rstrip().split('\t')
        f1 = line[0]
        rs = float(line[1])
        region_f1 = line[2]
        type_f1 = line[3]
        avgSig_f1 = float(line[6])
        tau_f1 = float(line[7])
        f2 = line[8]
        region_f2 = line[9]
        type_f2 = line[10]
        avgSig_f2 = float(line[13])
        tau_f2 = float(line[14])

        """transcript region"""
        transcript_region[f1] = region_f1
        transcript_region[f2] = region_f2

        """number of significant tests for each feature"""
        sig_f1 = [float(line[4]), float(line[5])]
        n_sig_f1 = sum(i < 0.05 for i in sig_f1)
        sig_f2 = [float(line[11]), float(line[12])]
        n_sig_f2 = sum(i < 0.05 for i in sig_f2)

        if f1 not in N_sigTest.keys():
            N_sigTest[f1] = n_sig_f1
        if f2 not in N_sigTest.keys():
            N_sigTest[f2] = n_sig_f2

        """average FDR adjusted P value"""
        if f1 not in avg_sigTest.keys():
            avg_sigTest[f1] = avgSig_f1
        if f2 not in avg_sigTest.keys():
            avg_sigTest[f2] = avgSig_f2

        """Kendall correlation with target class"""
        if f1 not in targetCorr_tau.keys():
            targetCorr_tau[f1] = tau_f1
        if f2 not in targetCorr_tau.keys():
            targetCorr_tau[f2] = tau_f2

        """number of correlated feature for each feature
        without region restriction"""
        for feature in [f1, f2]:
            N_featureCorr[feature] = N_featureCorr.get(feature, 0) + 1

        """pairwise Spearman correlation"""
        feature_pair = f1 + ':' + f2
        featureCorr_pair.append(feature_pair)
        featureCorr_pair_rs[feature_pair] = rs

        """complete feature list for sanity check at the end"""
        for feature in [f1, f2]:
            if feature not in feature_list:
                feature_list.append(feature)

####### test session
# print(N_sigTest, len(N_sigTest))
# print(featureCorr_pair, len(featureCorr_pair))
# print(feature_list, len(feature_list))
# print(N_featureCorr, len(N_featureCorr))
# print(avg_sigTest, len(avg_sigTest))
# print(featureCorr_pair_rs, len(featureCorr_pair_rs))
# print(targetCorr_tau, len(targetCorr_tau))
####### end of test session


####### create correlated features group for feature in the same region and type
###### for leadered: 5'transcript and 5'UTR are considered the same region
##### limit features to have same region and type
featureCorr_pair_list = {}
for feature_pair in featureCorr_pair:
    f1 = feature_pair.split(':')[0]
    f2 = feature_pair.split(':')[1]
    region_f1 = transcript_region[f1]
    region_f2 = transcript_region[f2]
    leadered_check = [region_f1, region_f2]
    leadered_check.sort()
    if type_f1 == type_f2:
        if region_f1 == region_f2 or leadered_check == ['fpr_UTR', 'fpr_transcript']:
            for feature in [f1, f2]:
                if feature not in featureCorr_pair_list.keys():
                    featureCorr_pair_list[feature] = []
            featureCorr_pair_list[f1].append(f2)
            featureCorr_pair_list[f2].append(f1)

##### create correlated feature group and remove duplicates
featureCorr_features_group = []
for k, v in featureCorr_pair_list.items():
    temp_list = []
    temp_list.append(k)
    for feature in v:
        temp_list.append(feature)
    temp_list.sort()
    if temp_list not in featureCorr_features_group:
        featureCorr_features_group.append(temp_list)

##### first sort by ascending length 
#### this is the order preserved if sort next by kendall correlation is tied
featureCorr_features_group.sort(key = len)

#### then sort feature group by the max kendall correlation in each group
def get_max_kendall(feature_group):
    feature_kendall = []
    for feature in feature_group:
        feature_kendall.append(targetCorr_tau[feature])
    max_kendall = max(feature_kendall)
    # print(max_kendall)
    # max_index = [index for index in range(len(feature_group)) if targetCorr_tau[feature_group[index]] == max_kendall]
    # print(max_index)
    # print(max_kendall, feature_group[max_index[0]])
    return max_kendall

featureCorr_features_group.sort(key = get_max_kendall, reverse = True)

####### test session
# print(featureCorr_pair_list, len(featureCorr_pair_list))
# print(featureCorr_features_group, len(featureCorr_features_group))

# def get_max_kendall(feature_group):
#     feature_kendall = []
#     for feature in feature_group:
#         feature_kendall.append(targetCorr_tau[feature])
#     max_kendall = max(feature_kendall)
#     # print(max_kendall)
#     max_index = [index for index in range(len(feature_group)) if targetCorr_tau[feature_group[index]] == max_kendall]
#     # print(max_index)
#     print(max_kendall, feature_group[max_index[0]])

# featureCorr_features_group.sort(key = get_max_kendall, reverse = True)
# for feature_group in featureCorr_features_group:
#     print(feature_group)
#     get_max_kendall(feature_group)
####### end of test session


####### select features
feature_selected = []
feature_excluded = []
def get_selected():
    for feature_group in featureCorr_features_group:
        N_sigTest_group = []
        N_featureCorr_group = []
        V_targetCorr_group = []
        for feature in feature_group:
            N_sigTest_group.append(N_sigTest[feature])
            N_featureCorr_group.append(N_featureCorr[feature])
            V_targetCorr_group.append(targetCorr_tau[feature])

        ####### test session
        # print(feature_group)
        # print(N_sigTest_group, len(N_sigTest_group))
        # print(N_featureCorr_group, len(N_featureCorr_group))
        # print(V_targetCorr_group)

        max_N_sigTest = max(N_sigTest_group)
        # print(max_N_sigTest)
        max_N_sigTest_list = [index for index in range(len(N_sigTest_group)) if N_sigTest_group[index] == max_N_sigTest]
        # print(max_N_sigTest_list)
        if len(max_N_sigTest_list) == 1:
            feature_keep = feature_group[max_N_sigTest_list[0]]
            if feature_keep not in feature_selected:
                # print(feature_keep)
                feature_selected.append(feature_keep)
                feature_group.remove(feature_keep)
            for feature in feature_group:
                if feature not in feature_excluded and feature not in feature_selected:
                    feature_excluded.append(feature)
        else:
            min_N_featureCorr = min(N_featureCorr_group)
            min_N_featureCorr_list = [index for index in range(len(N_featureCorr_group)) if N_featureCorr_group[index] == min_N_featureCorr]
            # print(min_N_featureCorr_list)
            if len(min_N_featureCorr_list) == 1:
                feature_keep = feature_group[min_N_featureCorr_list[0]]
                if feature_keep not in feature_excluded and feature_keep not in feature_selected:
                    # print(feature_keep)
                    feature_selected.append(feature_keep)
                    feature_group.remove(feature_keep)
                    for feature in feature_group:
                        if feature not in feature_excluded and feature not in feature_selected:
                            feature_excluded.append(feature)
            else:
                # print(feature_group)
                temp_group = []
                for feature in feature_group:
                    if feature not in feature_excluded:
                        temp_group.append(feature)
                if len(temp_group) == 1:
                    # print(temp_group)
                    feature_keep = temp_group[0]
                    if feature_keep not in feature_selected:
                        feature_selected.append(feature_keep)
                        # print(feature_keep)
                else:
                    # print(temp_group)
                    if len(temp_group) != 0:
                        temp_group_rs = []
                        for feature in temp_group:
                            feature_rs = 0
                            for feature_pair in featureCorr_pair_rs.keys():
                                f1 = feature_pair.split(':')[0]
                                f2 = feature_pair.split(':')[1]
                                if feature == f1 or feature == f2:
                                    # print(feature_pair)
                                    # print(featureCorr_pair_rs[feature_pair])
                                    feature_rs += featureCorr_pair_rs[feature_pair]
                            temp_group_rs.append(feature_rs)
                        # print(temp_group_rs)
                        min_group_rs = min(temp_group_rs)
                        min_group_rs_list = [index for index in range(len(temp_group_rs)) if temp_group_rs[index] == min_group_rs]
                        if len(min_group_rs_list) == 1:
                            feature_keep = temp_group[min_group_rs_list[0]]
                            if feature_keep not in feature_excluded and feature_keep not in feature_selected:
                                # print(feature_keep)
                                feature_selected.append(feature_keep)
                                feature_group.remove(feature_keep)
                                for feature in feature_group:
                                    if feature not in feature_excluded and feature not in feature_selected:
                                        feature_excluded.append(feature)
                        else:
                            # print(feature_group)
                            max_targetCorr = max(V_targetCorr_group)
                            # print(max_targetCorr)
                            max_targetCorr_list = [index for index in range(len(V_targetCorr_group)) if V_targetCorr_group[index] == max_targetCorr]
                            # print(max_targetCorr_list)
                            if len(max_targetCorr_list) == 1:
                                feature_keep = feature_group[max_targetCorr_list[0]]
                                if feature_keep not in feature_excluded and feature_keep not in feature_selected:
                                    feature_selected.append(feature_keep)
                                    feature_group.remove(feature_keep)
                                    for feature in feature_group:
                                        if feature not in feature_excluded and feature not in feature_selected:
                                            # print(feature)
                                            feature_excluded.append(feature)
                            else:
                                # print(feature_group)
                                avg_sigTest_group = []
                                for feature in feature_group:
                                    avg_sigTest_group.append(avg_sigTest[feature])
                                # print(avg_sigTest_group)
                                feature_keep = feature_group[avg_sigTest_group.index(min(avg_sigTest_group))]
                                if feature_keep not in feature_excluded and feature_keep not in feature_selected:
                                    # print(feature_keep)
                                    feature_selected.append(feature_keep)
                                    feature_group.remove(feature_keep)
                                    for feature in feature_group:
                                        if feature not in feature_excluded and feature not in feature_selected:
                                            feature_excluded.append(feature)

####### test session
# print(feature_selected, len(feature_selected))
# print(feature_excluded, len(feature_excluded))
# test = []
# for feature in feature_selected:
#     test.append(feature)
# for feature in feature_excluded:
#     test.append(feature)
# print([x for x in feature_list if x not in test])
####### end of test session


####### get feature selected
get_selected()
for feature in feature_excluded:
    print(feature)
