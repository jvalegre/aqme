"""
In this file, the program will split your data
into train and test, and then it'll determine
how good your model is as a predictor based on the
4 rules that we specified in Jeff's article
(R2 and Q2 of the test, k and R2 diff)


Also, add apply the "X random" method and tje Y random
method (they should be either in Mihais or Bens code). Also,
add new test Rob mentioned in Slack about adding 1s and 0s
with dummy function of numpy. Also, add Q-Q plots of errors to ensure
there are no systematic errors in a particular region.
Add outlier detector from the repository

Add options about PCAs using chunks of variables: the PCA
should be able to find the elbow in the scree plot, and then
grab that number of PCS to reduce the number of variables
(i.e. 20 variables correlated regarding sterics -- converted
to --  4 PCs used by the model, and the same for electronics,
and so on)

Add one hot encoding:
Here is basic one hot encoding.
The data contains several descriptors, and the regressions are first run on a train/test 70:30 split using the descriptors. The best RF model is R2 0.89 and RMSE 0.4 kcal/mol – pretty good!
The y-scramble and then x-scramble are run, and the regressions are much worse – so far so good!
BUT – the final regression uses no descriptors, just one hot encoding. The best RF model is R2 0.88 and RMSE 0.4 kcal/mol !!!
To support a descriptor-based model it is necessary to show the failure of the one hot regression on a test or other external prediction set.

Add filter for duplicates in the database
"""
