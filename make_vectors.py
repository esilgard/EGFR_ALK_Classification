# -*- coding: utf-8 -*-
'''author@esilgard'''
#
# Copyright (c) 2015-2017 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#
'''
script to make feature vectors for EGFR/ALK (and other) tests from pathology reports
this process consists of iterative text munging with feature vector additions
'''
from datetime import datetime
import os, re, json

# a file that (at minimum) contains the unique id of the report (instance)
# the pathology report accession number 
# and the raw text of the pathology report
# assumed to be a tab delimited file
INSTANCES_FILE = 'Input/IR 10469_Path_Text_2020-06-29.txt'
# identifiey columns for report id, accession, and report text (initial column = 0)
TEXT_COL = 7
REPORT_ID_COL = 2
ACC_NUM_COL = 5

# the test/marker to be classified - note this needs to be the standardized
# key that is present in the "test_patterns.json" resource file
TEST_NAME = 'ALK'

# one output file per batch, per test type
OUTPUT_FILE = 'Input' + os.sep + TEST_NAME + '_feature_vectors.txt'
RESOURCE_DIR = 'Resources'
STOP_LIST = '[\s\^](TO|THE|FOR|A|AN|AS|THIS|THAT|THESE|THEY|IN|OF|ON|OR|BY)( THE|A|AN)?[\s\$]'

# windows dictate the number of tokens on either side of the test name
# that are considered in the feature engineering
PRE_WINDOW = 10
POST_WINDOW = 10

    
def vector_creation():
    '''
    initial/main method for vector creation
    '''  
    def compile_patterns(pattern_dictionary, uppercase_boolean, cushion1, cushion2):
        '''
        helper method to compile patterns for each regex dictionary 
        with appropriate "cushions" before and after primary capture groups
        allows for optomizing pattern matching while 
        allowing for slightly better readibility in json docs
        '''
        compiled_d = {}
        for key, val in pattern_dictionary.items():
            compiled_patterns = []
            k = key.replace('<newline>','\n')
            for each in val:
                ## make sure match pattern is isolated from alphanumeric characters
                compiled_patterns.append(re.compile(r'' + cushion1 + \
                    '('+each+')' + cushion2, re.MULTILINE))
                if uppercase_boolean:
                    compiled_patterns.append(re.compile(r'' + cushion1 + \
                        '('+each.upper()+')' + cushion2, re.MULTILINE))
                
            compiled_d[k] = compiled_patterns
        return compiled_d
    
    # resources needed for text standardization and feature engineering
    test_pattern_file = '{}{}{}'.format(RESOURCE_DIR, os.path.sep, 'test_patterns.json')
    other_pattern_file = '{}{}{}'.format(RESOURCE_DIR, os.path.sep, 'other_keyword_patterns.json')
    section_pattern_file = '{}{}{}'.format(RESOURCE_DIR, os.path.sep, 'section_patterns.json')
    
    # load json objects containing regex patterns        
    test_patterns = json.load(open(test_pattern_file,'r'))
    other_patterns = json.load(open(other_pattern_file,'r'))
    section_patterns = json.load(open(section_pattern_file,'r'))

    # here we only uppercase all patterns (based on boolean flag)
    # as long as no regex character classes are used in pattern 
    # e.g. we don't want [\w] to turn into [\W]
    compiled_test_patterns = compile_patterns(test_patterns, True, '[\W\^]', '[\W$]')
    compiled_other_patterns = compile_patterns(other_patterns, False, '[\W\^]', '[\W$]')
    compiled_section_patterns = compile_patterns(section_patterns, True, '^', '$')

    counter = 0
    # loop through instances
    with open(OUTPUT_FILE, 'w') as out:
        num_processed = 0
        for lines in open(INSTANCES_FILE, 'r', encoding='unicode_escape').readlines()[1:]:   #, encoding='unicode_escape'
            separate_columns = lines.split('\t')
            if True:  # place holder for individual instance debugging
                counter += 1  
                pathnum = separate_columns[ACC_NUM_COL]
                report_id = separate_columns[REPORT_ID_COL]                
                text = separate_columns[TEXT_COL]

                # seems unnecessary, but need to maintain regex behavior
                text = text.replace('<newline>','\n').strip()
                # dictionary of feature counts for models
                vector = {} 
                
                # cytometry check
                vector['CYTO_RELATED_REPORT'] = int(get_cyto(text))
                # standardize test instance
                text, vector = strip_test_name(compiled_test_patterns[TEST_NAME], text, \
                    vector)
                # standardize all mentions of other test names in the text 
                text = make_standardized_text(compiled_test_patterns, text, \
                    compiled_other_patterns, counter, compiled_section_patterns)                
                
                ## get reference to other path accession feature
                other_acc_num, text = get_other_acc_num(text, pathnum)
                vector['OTHER_ACC_NUM_IN_TEXT'] = vector.get('OTHER_ACC_NUM_IN_TEXT', 0) + int(other_acc_num)
                # general check for language about technical difficulties
                vector['INSUFFICIENT'] = vector.get('INSUFFICIENT', 0) + int(get_insufficient(text))
                ## condense some duplicate standardizations    
                for string in ['OTHER_TEST','PUBLICATION','TEST_INSTANCE','IHC',
                    'PATHOLOGIST','BLOCK_ACC', 'SPECIFIC_MUT','MUT_ANALYSIS',
                    'FISH','AUTHOR']:
                    text = re.sub(string + r'[,.\(\):\-;andor' + string + \
                        ' ]{1,}' + string, ' ' + string + ' ', text)
                
                # condense coordinated test instance 
                # (eg test 1 and test 2 are pending) 
                # we dont want to break on test2
                text = re.sub('TEST_INSTANCE[,.\(\):\-;andor ]{1,}OTHER_TEST', 'TEST_INSTANCE', text)
                # strip out punctuation at the last minute for skipgram features 
                text = re.sub('[0-9]{2}[\-\\\/][0-9]{2}[\-\\\/][0-9]{2,5}',' DATE ',text)    
                text = re.sub('($|\s)[A-H][\)]?[.]',' SPECIMEN_LABEL', text)
                text = re.sub('[\"\(\\\)\-\/\']', ' ', text)                
                text = re.sub('[.,;:\?]', ' PUNCTUATION ', text)
                text = text.replace('[', ' ').replace(']', ' ')
                
                # strip out long trails of  ____ - 
                # there ARE single underscores IN the standardized features
                text = re.sub('_{2,}',' ', text)
                text = text.upper()
                # twice through to capture SOME MORE of the coordinations
                text = re.sub(STOP_LIST,' ', text)
                text = re.sub(STOP_LIST,' ', text)  
                text = text.split()  
                vector = make_ngrams(text, vector)                  
               
                vector['COUNT_TEST_INSTANCE'] = text.count('TEST_INSTANCE')
                if not vector['COUNT_TEST_INSTANCE']:
                    vector['NO_KEYWORD_IN_TEXT'] = 1    
                
                # output feature vectors to text specific file
                out.write(report_id)
                try:
                    for k,v in vector.items():
                        out.write('\t{}\t{}'.format(k, v))
                except:
                    print ('{} processed {}{}'.format(num_processed, report_id,k))
                    sys.exit()
                out.write('\n')
                num_processed += 1
   

def strip_test_name(regex, text, vector):
    '''
    make all test instances look the same 
    so that the model can be applied to all
    '''
    for expression in regex:
        # the '+' get lost in the [\W] buffer around tests - pull them out 
        # not pulling out '-', since it's ambiguous; minus or just a dash?
        # strips off the "cushion" from the end of the test instance pattern
        if re.search(expression.pattern[:-5] + '[\s]*[\+]',text):
            vector['post_window=POSITIVE'] = vector.get('post_window=POSITIVE',0) + 1
            vector['post_window=TEST_INSTANCE_POSITIVE'] = \
                vector.get('post_window=TEST_INSTANCE_POSITIVE',0) + 1
        text = re.sub(expression, ' TEST_INSTANCE ',text)       
    return text, vector


def get_other_acc_num(text, pathnum):
    '''
    attempt to acertain whether other/previous pathology reports are mentioned
    '''
    other_acc_bool = False
    pathnum = pathnum.replace('-', '')
    for acc in re.finditer('[\W]([\(]?[A-Z]{1,2}[\- ]?[\d]{2,4}[\- ]{1,3}[\d]{2,8}[\)]?)[\W]', text):
        current_acc = re.sub('[\- ]', '', acc.group(1))
        if current_acc == pathnum:            
            text = re.sub(re.sub('[\(\)]', '' ,acc.group(1)), ' THIS_ACC_NUM ', text)
        else:
            text = re.sub(re.sub('[\(\)]', '', acc.group(1)), ' OTHER_ACC_NUM ', text)
            other_acc_bool = True
    return other_acc_bool, text
    

def get_insufficient(text):
    '''
    attempt to capture insufficient samples EVEN when test name isn't mentioned
    '''
    if re.search('insufficient (tumor|sample)?', text, re.IGNORECASE) or\
        re.search('technical diff', text, re.IGNORECASE) or \
        re.search('(tumor|sample) insufficient', text, re.IGNORECASE):
        return True
    else:
        return False

def get_cyto(text):
    '''
    attempt to capture cytology related reports; unlikely to have reliable tests
    '''
    if re.search('(cytoprep)|(cytolog)',text, re.IGNORECASE):
        return True
    else:
        return False
   
def make_standardized_text(compiled_test_patterns, text, compiled_other_patterns, counter, compiled_section_patterns):
    def regex_sub(regex, standardization, text):        
        for expression in regex:
            text = re.sub(expression,' ' + standardization + ' ',text)
        return text
       
    ## kinda gross hack - this runs through twice to catch overlapping patterns
    ## (because of [\W] buffer in pattern match)    
    for i in range(2):
        ## replace all other tests in the text    
        for test, reg in compiled_test_patterns.items():            
            text = regex_sub(reg, ' OTHER_TEST ' , text)
        for section, pattern in compiled_section_patterns.items():
            text = regex_sub(pattern, section, text)
        ## replace patterns with variable standardizations
        for standardization,regex in compiled_other_patterns.items():
            text = regex_sub(regex, standardization, text)
    return text
    
def make_ngrams(text, vector):
    '''
    loop through tokens to look for windows around test instances
    create unigrams, bigrams, and skipgrams
    '''
    for v in range(len(text)):
        current = text[v]                          
        if current == 'TEST_INSTANCE':
            try:
                ## find the first and closest second section to the test instance
                reversed_snippet = [text[b] for b in range(v-1, -1, -1)]
                if '_SECTION_' in reversed_snippet:
                    section = reversed_snippet[reversed_snippet.index('_SECTION_') - 1]
                    vector['SECTION='+section] = vector.get('SECTION='+section,0) + 1
            except IndexError:
                pass                        
            
            # include standardized test name in feature set
            vector[TEST_NAME] = 1
            # create window around test mention without extending
            # beyohnd beginning or end of full text
            pre_window_index = max(v - PRE_WINDOW,0)                           
            post_window_index = min(len(text), v + POST_WINDOW)                            
            window = text[pre_window_index:post_window_index]
            
            # break window size for new sections or other tests, etc
            breaks = [i for i in range(len(window)) if \
            (window[i] == "_SECTION_" or window[i] == "PUNCTUATION" \
             or  window[i] == "SPECIMEN_LABEL" or window[i] == "OTHER_TEST")]                                        
            breaks.append(0)                            
            pre_break = max([x + 1 for x in breaks if x < len(window)/ 2])
            breaks[-1] = len(window)    
            post_break = min([x for x in breaks if x > len(window)/ 2])                          
            window = text[pre_window_index + pre_break:pre_window_index + post_break]
            post_window_index = post_break + pre_window_index
            pre_window_index = pre_break + pre_window_index
      
            # move from the current token backwards for pre window features
            pre_pointer = v-1      
            if v > 0 and pre_pointer >= pre_window_index :\
                vector['immediately_pre_window=' + text[pre_pointer]]=1
            while pre_pointer >= pre_window_index:                              
                #unigrams 
                vector['pre_window='+text[pre_pointer]] = \
                    vector.get('pre_window=' + text[pre_pointer],0) + 1
                if  v - 1 > pre_window_index:
                    #bigrams
                    vector['pre_window=' + text[pre_pointer-1] + '_' + \
                        text[pre_pointer]] = vector.get('pre_window=' + \
                        text[pre_pointer-1] + '_' + text[pre_pointer],0) + 1
                    if v - 2 > pre_window_index:
                        #one skipgram                                  
                        vector['pre_window=' + text[pre_pointer-2] + \
                            '_' + text[pre_pointer]] = \
                            vector.get('pre_window=' + text[pre_pointer-2] + \
                            '_' + text[pre_pointer],0) +1
                        if v - 3 > pre_window_index:
                            #two skipgram
                            vector['pre_window=' + text[pre_pointer-3] + \
                                '_' + text[pre_pointer]] = vector.get('pre_window=' + \
                                text[pre_pointer-3] +'_' + text[pre_pointer], 0) + 1 
                pre_pointer -= 1
            
          
            ## move from the current token forwards for post window features
            post = v+1
            if post < len(text) and post < post_window_index: vector['immediately_post_window=' + text[post]] = 1
            while post < min(len(text) ,post_window_index):                                   
                #unigrams                                                                    
                vector['post_window=' + text[post]] = vector.get('post_window=' + text[post], 0) + 1
                if  post < post_window_index - 1:
                    #bigrams
                    vector['post_window='+text[post]+'_'+text[post+1]] = vector.get('post_window='+text[post]+'_'+text[post+1],0)+1
                    if post < post_window_index - 2:
                        #one skipgram
                        vector['post_window='+text[post]+'_'+text[post+2]]= vector.get('post_window='+text[post]+'_'+text[post+2],0)+1
                        if post < post_window_index - 3:
                            #two skipgram
                            vector['post_window='+text[post]+'_'+text[post+3]] = vector.get('post_window='+text[post]+'_'+text[post+3],0)+1
                post+=1 
    return vector

if __name__ == '__main__':
    ## timeit variable for performance testing ##
    BEGIN = datetime.today()
    print ('vector creation started at {}'.format(BEGIN))
    vector_creation()
    ## timeit - print out the amount of time it took to process all the reports ##
    print ('{} seconds to create vectors'.format((datetime.today()-BEGIN).days * 86400 + \
        (datetime.today()-BEGIN).seconds))

