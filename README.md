# EGFR_ALK_Classification
===========================
a hybrid model to classify EGFR and ALK test use, results, and methods from free text pathology reports 

### This system classifies EGFR and ALK molecular tests from the free text of pathology reports.
None of the SVM models needed for end to end classification are included in this repository


*There are four principal algorithms*
-------------------------------------

* reported: this is a combination of a keyword filter and an SVM that classifies all reports in the input as either "Reported" or "NotReported" for each test
    (for the purposes of sensitivity/specificity metrics and for consumption by downstream algorithms, "Reported" is considered a positive classification and "NotReported" a negative)

* insufficient: this is an SVM that classifies all reports classified as "NotReported" by the previous reported algorithm as either "Insufficient" or "Unknown" for each test
    (for the purposes of sensitivity/specificity metrics and for consumption by downstream algorithms, "Insufficient" is considered a positive classification and "Unknown" a negative)

* positive: this is an SVM that classifies all reports classified as "Reported" by the previous reported algorithm as either "Positive" or "Negative" for each test
    (for the purposes of sensitivity/specificity metrics and for consumption by downstream algorithms, "Positive" is considered a positive classification and "Negative" a negative)

* method: this is an SVM that classifies all reports classified as "Reported" by the previous reported algorithm as either the standard testing methodology or "OTHER" for each test
    (for the purposes of sensitivity/specificity metrics and for consumption by downstream algorithms, the standard testing methodology ("MutationalAnalysis" for EGFR and "FISH" for ALK) is considered a positive classification and "OTHER" a negative)
    

- make_vectors.py is a standalone script that creates individual feature vectors for each report&test instance
    this  will need to be modified to include the name/location of the input file (now expected to be in the "Input" sub directory)
    also the input file is expected to be a tab delimited text file, the column indices of the text, report identifiers, and accession numbers will need to be manually edited. (see TEXT_COL, REPORT_ID_COL, ACC_NUM_COL). All instances will have one single feature vector per molecular test (although they will have multiple arrays created in the pipeline; one per algorithm)

- svm_pipeline.py is the main script to run the end to end classification pipeline
    - vector_to_array.py creates arrays based on feature sets (not included in the public repository) from the training set; one array   per instance, per test, per algorithm
    - decoder.py uses the svm models learned in training (not included in the public repository) to classify each instance
    - final_output.py aggregates all the individual classifications to produce one label per report instance (e.g. EGFR Negative by MutationalAnalysis)

Currently the vector creation and classification pipeline are run for one test at a time, the "TEST_NAME" will have to be manually edited in both make_vectors.py, as well as svm_pipeline.py


Internal validation performance as well as a project overview is reported in the attached abstract "Validation of Natural Language Processing (NLP) for Automated Ascertainment of EGFR and ALK Tests in SEER Cases of Non-Small Cell Lung Cancer (NSCLC)"
