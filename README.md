Source codes of algorithms and datasets for our paper "Ahpatron: A New Budgeted Online Kernel Learning Machine with Tighter Mistake Bound", accepted in AAAI 2024.

We implement all algorithms with R on a Windows machine with 2.8 GHz Core(TM) i7-1165G7 CPU. execute each experiment 10 times with random permutation of all datasets and average all of the results.

To run the code, you must set the paths following the code, or set new paths. The default path of codes is "D:/experiment/AAAI/AAAI2024/code". The path of datasets is "D:/experiment/AAAI/AAAI2024/dataset". The store path is "D:/experiment/AAAI/AAAI2024/Result/".

The baseline algorithms include: Projectron, Projectron++, FOGD, BOGD++, NOGD and POMDR. Our algorithm is Ahpatron.

The datasets are downloaded from: https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/ and http://archive.ics.uci.edu/ml/datasets.php

binary classification datasets: w8a (Num:49749, Fea:300), cod-rna (Num:271617, Fea:8), phishing (Num:11055, Fea:69), a9a (Num:48842, Fea:123), SUSY (Num:50000, Fea:18), ijcnn1 (Num:141691, Fea:22)
