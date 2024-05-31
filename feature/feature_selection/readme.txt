
############## final complete version of feature selection 
####### get files from ../feature_selection_vTest
cp ../feature_selection_vTest/features_CorrSigSelection_halfLifeCls_logPhase_leadered.R ./
cp ../feature_selection_vTest/select_features_vTest.py ./select_corrFeatures.py

####### get correlations of features
###### features_CorrSigSelection_*.R

####### get excluded features for different class
###### leadered genes
python select_corrFeatures.py corrSig_halfLifeCls_hypoxia_leadered.txt > feature_halfLifeCls_hypoxia_excludedByLeadered.txt
python select_corrFeatures.py corrSig_halfLifeCls_logPhase_leadered.txt > feature_halfLifeCls_logPhase_excludedByLeadered.txt
python select_corrFeatures.py corrSig_halfLifeFcCls_hypoxiaToLogPhase_leadered.txt > feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeadered.txt

###### leaderless genes
python select_corrFeatures.py corrSig_halfLifeCls_hypoxia_leaderless.txt > feature_halfLifeCls_hypoxia_excludedByLeaderless.txt
python select_corrFeatures.py corrSig_halfLifeCls_logPhase_leaderless.txt > feature_halfLifeCls_logPhase_excludedByLeaderless.txt
python select_corrFeatures.py corrSig_halfLifeFcCls_hypoxiaToLogPhase_leaderless.txt > feature_halfLifeFcCls_hypoxiaToLogPhase_excludedByLeaderless.txt

