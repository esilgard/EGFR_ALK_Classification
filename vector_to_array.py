# -*- coding: utf-8 -*-
'''author@esilgard'''
#
# Copyright (c) 2015-2017 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#

import os

def vector_to_array(sparse_array, feature_d, batch_set, input_vector_file):    
    with open(sparse_array,'w') as out:
        with open(input_vector_file,'r') as f:   
            for lines in f:                
                l=lines.strip().split('\t')
                if l[0].split('\t')[0] in batch_set:
                    out.write(l[0])
                    #sparse vector for binary features (skip over feature counts)
                    for f in range(1,len(l),2):                           
                        if l[f] in feature_d:
                            out.write(' '+feature_d[l[f]])
                    out.write('\n')

def main(test_name, algorithm, train_batch):
    '''
    get appropriate instances for classification depending on the algorithm
    '''
    input_vector_file = test_name + '_feature_vectors.txt'
    
    if algorithm == 'reported':
        batch_set = set([y.split('\t')[0] for y in \
            open('rule_based_reported'  + os.sep + test_name + '_pos_instances.txt','r').readlines()])
    elif algorithm == 'positive' or algorithm == 'method':
        batch_set = set([y.split('\t')[0] for y in\
            open('reported' + os.sep + test_name + '_pos_instances.txt','r').readlines()])
    elif algorithm == 'insufficient':
        test_reported_rules = set([y.split('\t')[0] for y in \
            open('rule_based_reported' + os.sep + test_name + '_neg_instances.txt','r').readlines()])
        test_reported_svm = set([y.split('\t')[0] for y in \
            open('reported' + os.sep + test_name + '_neg_instances.txt','r').readlines()])
        ## combine rule based and svm 'reported' output
        batch_set = test_reported_rules.union(test_reported_svm)
    
    ## make smaller sparse array file from feature vectors
    array_file = algorithm + os.sep + test_name + '_sparse_arrays.txt'
    feature_mapping_file = algorithm + os.sep + train_batch + '_features_mapping.txt'
    feature_d = dict((a.split('\t')[0],a.strip().split('\t')[1]) for a in \
        open(feature_mapping_file,'r').readlines())
    vector_to_array(array_file, feature_d, batch_set, input_vector_file)

    return len(feature_d)

