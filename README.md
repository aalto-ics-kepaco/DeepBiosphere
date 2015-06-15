# DeepBiosphere

The codes in this repository can be used to perform either sparse canonical correlation analysis (SCCA) or kernel canonical correlation analysis (KCCA) on a set of two-view data. The ClustergramFig, CorrelationPlot and ScorePlot functions can be used to visualise the results. The functions PredictiveCC_GaussianKernel and PredictiveCC perform 3-fold cross validation in order to determine the optimal value of the sparsity hyper parameter in SCCA when a Gaussian or linear kernel are applied respectively. The stability of the obtained value can be assessed by the StabilityAnalysis function. The statistical siginificance of the obtained canonical correlations can be evaluated by randomisation tests that are implemented in Significance_KCCA and Significance_SCCA.


© 09/06/2015 Viivi Uurtio, Aalto University
viivi.uurtio@aalto.fi

Commercial use is not allowed.