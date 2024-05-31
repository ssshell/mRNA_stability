
import sys
import pandas as pd


file = open(sys.argv[1], 'r')
feature = sys.argv[2]

data = pd.read_csv(file, sep = "\t", header = None)
data.columns = [" ", feature]
# print(data)

data_onehot = pd.concat([data, pd.get_dummies(data.iloc[:, 1], prefix = feature)], axis = 1)
data_onehot.drop([feature], axis = 1, inplace = True)
print(data_onehot.to_csv(index = False, sep = "\t"))
