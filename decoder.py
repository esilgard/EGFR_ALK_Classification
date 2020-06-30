'''author@esilgard'''
#
# Copyright (c) 2015-2017 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#

### turn feature vectors into arrays for sci-learn SVM ###
import os
import numpy as np
from scipy.sparse import dok_matrix
import sklearn
from sklearn.externals import joblib
print('The scikit-learn version is {}.'.format(sklearn.__version__))

def main(num_features, model_file, algorithm, train_batch, test_name, label):
    feature_d = dict((x.strip().split('\t')[1],x.split('\t')[0]) for x in \
        open(algorithm + os.sep + train_batch + '_features_mapping.txt','r').readlines())
    sparse_arrays_file =  algorithm + os.sep + test_name + '_sparse_arrays.txt'
    class_map_file = label + '_label_mapping.txt'

    # maintain order of instances for decoding batch classification
    instances = [x.strip().split() for x in \
        open(sparse_arrays_file,'r').readlines()]
    class_map=dict((x.strip().split('\t')[1], int(x.strip().split('\t')[0])) \
        for x in open(class_map_file,'r').readlines())
    reverse_class_map = dict((v,k) for k,v in class_map.items())    
    X = dok_matrix((len(instances), num_features),dtype = np.float64)
    index = 0
    print ('{} instances'.format(len(instances)))
    for index in range(len(instances)):
        #instance_id = instances[index][0]
        features =  instances[index][1:]
        for feat in features:
            if feat in feature_d.keys():
                X[index,int(feat)]=1

    X.tocsc()
    clf = joblib.load(model_file)
    # the insufficient model was pickled with a float for feature nums? temp hack fix
    clf.named_steps.feature_selection.k = int(clf.named_steps.feature_selection.k)

    try:
        output = clf.predict(X)
        print ('instances classified\n')
     
        # output files (depending on pos neg flag) 
        # so they may be used by downstream algorithms
        output_pos =  algorithm + os.sep + test_name + '_pos_instances.txt'
        output_neg =  algorithm + os.sep + test_name + '_neg_instances.txt'
        
        with open(output_pos,'w') as pos_out:
            with open(output_neg,'w') as neg_out:
                for i in range(len(output)):                
                    system_out = reverse_class_map[output[i]]
                    if algorithm == 'positive' : 
                        positive_hit = set(['Positive'])
                    elif algorithm == 'method':
                        positive_hit = set(['MutationalAnalysis','FISH'])
                    elif algorithm == 'insufficient':
                        positive_hit = set(['Insufficient'])
                    else:
                        positive_hit = set(['Negative','Positive', 'Reported'])
                    if system_out in positive_hit: 
                        pos_out.write(instances[i][0] + '\t' + system_out + '\n')
                    else:
                        neg_out.write(instances[i][0] + '\t' + system_out + '\n')
    except:
        print ('ERR: outputs not processed')