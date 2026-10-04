# Results of the latest run

2026-09-24: 24 children, 5 folds x 1 repeats, quick = TRUE. Balanced accuracy and macro F1 on held-out children.

|method                             | balanced_accuracy|    sd| macro_f1| folds| fit_seconds|
|:----------------------------------|-----------------:|-----:|--------:|-----:|-----------:|
|XGBoost, all features              |             0.963| 0.051|    0.963|     5|         3.7|
|XGBoost, top-27 in-fold            |             0.960| 0.043|    0.959|     5|         6.2|
|SOM features, unit majority        |             0.884| 0.068|    0.882|     5|         0.7|
|DTW-SOM segments, unit majority    |             0.716| 0.088|    0.673|     5|         1.7|
|SOM features, 7 super-clusters     |             0.696| 0.052|    0.640|     5|         0.7|
|DTW-SOM windows, unit majority     |             0.576| 0.027|    0.478|     5|         0.5|
|DTW-SOM segments, 7 super-clusters |             0.531| 0.077|    0.435|     5|         1.7|
|DTW-SOM windows, 7 super-clusters  |             0.494| 0.024|    0.379|     5|         0.5|
|Duration only                      |             0.402| 0.055|    0.361|     5|          NA|
