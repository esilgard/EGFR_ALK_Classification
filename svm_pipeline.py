# -*- coding: utf-8 -*-
'''author@esilgard'''
#
# Copyright (c) 2015-2017 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#

import os
import warnings
import vector_to_array
import decoder
import final_output


TEST_NAME = 'ALK'
TRAIN_BATCH = 'internal_validation'
VECTOR_FILE = TEST_NAME + '_feature_vectors.txt'
ALGORITHM_ORDER = [('reported','result'),('insufficient','result'),
    ('method','method'),('positive','result')]
FINAL_OUTPUT_DIRECTORY = 'final_output'

# initial output files for keyword rule based filter
pos_inst = 'rule_based_reported' + os.path.sep + TEST_NAME + '_pos_instances.txt'
neg_inst = 'rule_based_reported' + os.path.sep + TEST_NAME + '_neg_instances.txt'

def run_pipeline():  
    '''
    pipeline for classification of EGFR and ALK test use, result, and method
     - reported will further classify the reports that passed through 
        the rule based filter as "reported"
     - positive will ascertain the test result (pos or neg) for the reports 
        that pass through the reported SVM as "reported"
     - method will ascertain the test method (FISH, MutAnalysis, Other) for 
        reports that pass through the reported SVM as "reported"
     - insufficient will ascertain insuff vs unknown for the reports that 
         were labeled "not reported" by rule based and svm reported classifiers
     - final class labels will be output to the FINAL_OUTPUT_DIRECTORY 
    '''
    # rule based keyword filter
    with open(pos_inst, 'w') as pos_out:
        with open(neg_inst, 'w') as neg_out:
            for instances in open(VECTOR_FILE,'r').readlines():
                vec = instances.strip().split('\t')
                if 'NO_KEYWORD_IN_TEXT' in vec:
                    neg_out.write(vec[0] + '\tNotReported\n')
                else:
                    pos_out.write(vec[0] + '\tReported\n')
    
    # loop through SVM classifiers; "positive" and "negative" instances
    # will be written to algorithm specific directories s
    # o they can be used by subsequent algorithms
    for algorithm_tuple in ALGORITHM_ORDER:  
        algorithm = algorithm_tuple[0]
        label = algorithm_tuple[1]
        
        model_file = algorithm + os.path.sep + TRAIN_BATCH + '.pkl' 
        num_features = int(open(algorithm + os.path.sep + 'num_features.txt', \
            'r').read().strip())
        print 'classifying with model ',model_file, '-', num_features, 'features'
        # turn text feature vector into integer array
        vector_to_array.main(TEST_NAME, algorithm, TRAIN_BATCH)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fxn()            
            decoder.main(num_features, model_file, algorithm, TRAIN_BATCH, TEST_NAME, label)
    final_output.output_final_class_labels(TEST_NAME, FINAL_OUTPUT_DIRECTORY)
    
def fxn():
    '''
    silence ski-kit learn deprication warning
    '''
    warnings.warn("deprecated", DeprecationWarning)
   
if __name__ == "__main__":
    run_pipeline()
