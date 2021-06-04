################ MALDI-TOF MS data analysis ################ 

### Spectra pre-processing: feature matrix generation ###

library("MALDIquant")
library("MALDIquantForeign")
library("cluster")
library("factoextra")
library("binda")
library("dplyr")
library("ggplot2")
library("crossval")
library("caret")

Spectra_list <- importBrukerFlex(Espectra_path, verbose=FALSE)

Spectra_list <- trim(Spectra_list)
Spectra_list <- transformIntensity(Spectra_list, method = "sqrt")
Spectra_list <- smoothIntensity(Spectra_list, method = "SavitzkyGolay",
                       halfWindowSize = 40)
Spectra_list <- removeBaseline(Spectra_list, method = "SNIP", iterations = 100)
Spectra_list <- calibrateIntensity(Spectra_list, method = "TIC")
Spectra_list <- alignSpectra(Spectra_list, halfWindowSize = 40,SNR = 4,
                    tolerance = 0.2, warpingMethod = "quadratic")
Spectra_list <- averageMassSpectra(Spectra_list, labels = spot.factor, method = "sum")

# spot.factor are the categories of each spectra.

peaks <- detectPeaks(Spectra_list, SNR = 4, 
                     method = "MAD", halfWindowSize = 40)

peaks <- binPeaks(peaks, tolerance = 0.2)

peaks <- filterPeaks(peaks, minFrequency = c(0.33, 0.33, 0.33),
                     labels = spot.factor,
                     mergeWhitelists=TRUE)

featureMatrix_Int <- intensityMatrix(peaks, Spectra_list)

# spot.factor.2 are the categories of each average spectra.

thr <- optimizeThreshold(featureMatrix, spot.factor.2)
featureMatrix_dicho <- dichotomize(featureMatrix, thr)

### Unsupervised statistical analysis ### 

### Top peaks selected by the binary discriminant analysis (BDA) algorithm

br <- binda.ranking(featureMatrix_dicho, spot.factor.2, verbose = FALSE)
top.b <- br[Number_of_top_peaks]

### (a) Hierarchical k-means clustering (Hkmc)-Principal Component Analysis (PCA) cluster plot 

K.num <- 3 # Number of clusters
rownames(featureMatrix_dicho) <- spot.factor.2
km1.res.t <-hkmeans(featureMatrix_dicho[, top.b], K.num)

fviz_cluster(km1.res.t, 
             frame.type = "norm", 
             frame.level = 0.95,
             repel = TRUE)

### (b) Hierarchical k-means clustering-PCA cluster composition

Results <- data.frame(spot.factor.2)

Results$spot.factor.2 <- factor(Results$spot.factor.2
                    levels=c("Cnt","IS","LPS"))

Results <- Results %>%
    mutate(HKmeans = km1.res.t$cluster)

Results %>% 
    group_by(HKmeans) %>% 
    count(spot.factor.2) %>% 
    mutate(prop = (n/sum(n))*100) %>%
    ggplot(aes(x = factor(HKmeans.9), y = prop,
               label = paste(round(prop),"%"),
               fill = factor(spot.factor.2)))+ 
    geom_bar(stat = "identity") + 
    geom_text(position = "stack",
              aes(ymax = 100),vjust = 1.2) +
    scale_fill_manual(values = c("IS"="#3F9AF4", 
                                 "LPS"="#F44A3F", 
                                 "Cnt"="#229954")) 

### Machine learning analysis ### 

### Train/test splits

rownames(featureMatrix_dicho) <- spot.factor.2
set.seed(seed)
ind <- sample(2, nrow(featureMatrix_dicho), replace = TRUE, 
              prob = c(0.4, 0.6))

Train <- featureMatrix_dicho[ind == 1, ] 
Test <- featureMatrix_dicho[ind == 2, ]

Ytrain11 <- rownames(Train)
Ytest1 <- rownames(Test)
Ytrain11 <- as.factor(Ytrain11)
Ytest1 <- as.factor(Ytest1)

### BDA 

predfun1 <- function(Xtrain, Ytrain, Xtest, Ytest, selPeaks) {
    binda.out = binda(Xtrain[, selPeaks, drop = FALSE], Ytrain, verbose = FALSE)
    ynew = predict.binda(binda.out, Xtest[, selPeaks, drop=FALSE], verbose = FALSE)$class 
    cm = crossval::confusionMatrix(Ytest, ynew, negative = "CNT") 
    return(cm)  
}

br = binda.ranking(featureMatrix_Int, spot.factor.2)

# 20 peaks CV

ourPeaks.20 <- br[_20_top_peaks]
    
cvp.train.20 <- crossval::crossval(predfun1, Train, Ytrain11, 
                                   K = 5, B = 20, selPeaks = ourPeaks.20, verbose = FALSE)
c1 <- diagnosticErrors(cvp.train.20$stat)


cvp.test.20 <- crossval::crossval(predfun1, Test, Ytest1, 
                                  K = 5, B = 20, selPeaks = ourPeaks.20, verbose = FALSE)
c2 <- diagnosticErrors(cvp.test.20$stat)

# The same code but changing the number of top selected peaks is run to evaluate the performances with
# 5, 10 and 15 best peaks. 

### Rf 

# Both the train and the test sets were subsetted with the top 5, 15 and 20 best peaks.
top_five <- br[_5_top_peaks]
Train <- Train[, top_five]

ctrl <- trainControl(method = "cv",
                     number = 5,
                     repeats = 20,
                     summaryFunction = multiClassSummary,
                     search = "grid") 

set.seed(1) 
x <- ncol(Train)
mtry <- seq(4, ncol(x) * 0.8, 2)
hyper_grid <- expand.grid(.mtry = mtry)

Rf_caret_model.CV <- train(Ytrain11 ~ .,
                               data = Train, 
                               method = "rf",
                               metric = "Accuracy",
                               trControl = ctrl,
                               tuneGrid = hyper_grid)

pred.cv <- predict(object = Rf_caret_model.CV, 
                   newdata = Test)

Tabla.CV.Diagnos <- crossval::confusionMatrix(predicted  = pred.cv,       
                                            actual = Ytest1, 
                                            negative = "CTL" )

# Rf best peaks

Picos.CV <- data.frame(Rf_caret_model.CV[["finalModel"]][["importance"]])











