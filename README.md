# EGFR_ALK_Classification
a hybrid model to classify EGFR and ALK test use, results, and methods from free text pathology reports 

This system classifies EGFR and ALK molecular tests from the free text of pathology reports.
None of the SVM models needed for end to end classification are included in this repository


*There are four principal algorithms:

reported: this is a combination of a keyword filter and an SVM that classifies all reports in the input as either "Reported" or "NotReported" for each test
    (for the purposes of sensitivity/specificity metrics and for consumption by downstream algorithms, "Reported" is considered a positive classification and "NotReported" a negative)

reported: this is an SVM that classifies all reports classified as "NotReported" by the previous reported algorithm as either "Insufficient" or "Unknown" for each test
    (for the purposes of sensitivity/specificity metrics and for consumption by downstream algorithms, "Insufficient" is considered a positive classification and "Unknown" a negative)

positive: this is an SVM that classifies all reports classified as "Reported" by the previous reported algorithm as either "Positive" or "Negative" for each test
    (for the purposes of sensitivity/specificity metrics and for consumption by downstream algorithms, "Positive" is considered a positive classification and "Negative" a negative)

positive: this is an SVM that classifies all reports classified as "Reported" by the previous reported algorithm as either the standard testing methodology or "OTHER" for each test
    (for the purposes of sensitivity/specificity metrics and for consumption by downstream algorithms, the standard testing methodology ("MutationalAnalysis" for EGFR and "FISH" for ALK) is considered a positive classification and "OTHER" a negative)
    

- make_vectors.py is a standalone script that creates individual feature vectors for each report&test instance

- svm_pipeline.py is the main script to run the end to end classification pipeline



Internal validation performance as well as a project overview is reported in the attached abstract "Validation of Natural Language Processing (NLP) for Automated Ascertainment of EGFR and ALK Tests in SEER Cases of Non-Small Cell Lung Cancer (NSCLC)"
