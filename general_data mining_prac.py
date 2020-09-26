from Orange.data import Table,Domain
from Orange.classification import NNClassificationLearner
from Orange.evaluation import CrossValidation, scoring
from Orange.data import Table
from sklearn.cluster import KMeans
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report,confusion_matrix
import warnings
import Orange
import Orange.preprocess
import pandas as pd

data = Table.from_file('white wine.csv')
CLabel = Orange.preprocess.Discretize()
CLabel.method = Orange.preprocess.discretize.EqualWidth(n=3)
newCLabel = CLabel(data[:,0])
FeatureV = data[:,1:].domain.variables
CLabelSet = newCLabel.domain.variables
WineDomain = Domain(FeatureV, CLabelSet)
data = Table.from_table(domain=WineDomain, source=data)

TLearner = NNClassificationLearner(hidden_layer_sizes=(10,1),max_iter=4000)
evalR = CrossValidation(data, [TLearner], k=10)
print("Accuracy: {:.3f}".format(scoring.CA(evalR)[0]))
print("AUC: {:.3f}".format(scoring.AUC(evalR)[0]))

n = sum(1 for d in data if (d["quality"] == 1.0 or d["quality"] == 2.0 or d["quality"] == 9.0))
SubSet = Table(data.domain,[d for d in data if (d["quality"] < 4.0 or d["quality"] > 8.0)])
for d in SubSet:
    del data[d]
CLabel = Orange.preprocess.Discretize()
CLabel.method = Orange.preprocess.discretize.EqualWidth(n=3)
DataSet = CLabel(data[:,0])
FeatureV = data[:,1:].domain.variables
newCLabelSet = DataSet.domain.variables
wineDomain = Domain(FeatureV, newCLabelSet)
data = Table.from_table(domain=wineDomain, source=data)
TLearner = NNClassificationLearner(hidden_layer_sizes=(10,1),max_iter=4000)
evalR = CrossValidation(data, [TLearner], k=10)
print("Accuracy after deletion: {:.3f}".format(scoring.CA(evalR)[0]))
print("AUC after deletion: {:.3f}".format(scoring.AUC(evalR)[0]))
    
X = data.drop('quality',axis = 1)
Y = data['quality'].astype(float)
XTrain,XTest,YTrain,YTest= train_test_split(X,Y)
Scal = StandardScaler()
Scal.fit(XTrain)
yTr = YTrain.values.reshape(-1,1)
yTe = YTest.values.reshape(-1,1)
Scal.fit(yTr)
M = MLPClassifier(hidden_layer_sizes = (10,3),max_iter = 6000)
M.fit(XTrain,yTr.ravel())
Predics = M.predict(XTest)
print(confusion_matrix(yTe,Predics))
print(classification_report(yTe,Predics))

data = Table('ionosphere')
K = KMeans(n_clusters=2,random_state = 0)
Clus = K.fit(data)

Siz = sum(1 for d in data)
print(Siz)

Dict = {i: np.where(Clus.labels_ == i)[0] for i in range(Clus.n_clusters)}


gO = 0
bO = 0
g1 = 0
b1 = 0
for x in Dict[0]:
    if data[x]["y"] == "g":
        gO = gO + 1
    elif data[x]["y"] == "b":
        bO = bO + 1 
for y in Dict[1]:
    if data[y]["y"] == "g":
        g1 = g1 + 1
    elif data[y]["y"] == "b":
        b1 = b1 + 1
print("{} & {}  : {}".format("cluster 0", "label = g", gO))
print("{} & {}  : {}".format("cluster 0", "label = b", bO))
print("{} & {}  : {}".format("cluster 1", "label = g", g1))
print("{} & {}  : {}".format("cluster 1", "label = b", b1))


