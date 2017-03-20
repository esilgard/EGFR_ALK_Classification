# -*- coding: utf-8 -*-
'''author@esilgard'''
#
# Copyright (c) 2015-2017 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#
import os

def output_final_class_labels(test_name, FINAL_OUTPUT_DIRECTORY):
    '''
    aggregate the results from the insufficient, positive, and method
    classifications and produce a final label for the report
    '''
    report_label_dict = {}
    for algorithm in ['insufficient','positive','method']:
        system_labels = [x.strip().split('\t') for x in open(algorithm + \
            os.sep + test_name + '_pos_instances.txt','r').readlines()] + \
            [y.strip().split('\t') for y in open(algorithm + os.sep + \
             test_name + '_neg_instances.txt','r').readlines()]
        for each in system_labels:
            report_id = each[0]
            report_label_dict[report_id] = report_label_dict.get(report_id, {})
            report_label_dict[report_id][algorithm] = each[1]
    
    with open(FINAL_OUTPUT_DIRECTORY + os.sep + test_name + \
        '_final_output.txt','w') as out:
        for report, algorithm_labels in report_label_dict.items():
            out.write(report + '\t')
            if 'positive' in algorithm_labels:
                out.write(algorithm_labels['positive'] + ' by ' + \
                    algorithm_labels['method'] + '\n')
            else:
                out.write(algorithm_labels['insufficient'] + '\n')
        
    
