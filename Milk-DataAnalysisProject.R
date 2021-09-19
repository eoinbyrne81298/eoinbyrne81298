#Question 1 -Reading data and removing sample at random
data1 <- read.csv("Grain_data.csv")
set.seed(17419782)
n <- nrow(data1)
data <- data1[-c(sample(1:n,1)),]

#Question 2 - Data inspection/visualization

#Checking variable names
colnames(data)

#Scatterplot of glucose and ethanol values of samples to check for any obvious relationship.
plot(data$Glucose,data$Ethanol, xlab = "Glucose",ylab="Ethanol")
#Notice that some samples seem to be grouped together linearly in the graph.
#There are several groups of samples for which there seems to be a linear relationship between their glucose and ethanol values.  

#Two pairs plots of NIR data corresponding to spectra that are close together
pairs(data[,7:11])
pairs(data[,100:104])
pairs(data[,225:229])
#Notice that each of the plots have a similar structure in both of these pairs plots indicating that the absorbance values for close together spectra are correlated.

#Pairs plot of NIR data corresponding to spectra far from eachother. 
pairs(data[,c(5,95,185,230)])
#More variation in the graphs in this pairs plot. This is an indication of these spectra being less correlated.

#Question 3
#Creating matrix of NIR data of the samples
NIRdata = data[,2:236]

#Calculating standard deviation of each column/spectra.
sds = apply(NIRdata,2,sd)

#Dividing each column by its std dev to remove affect of NIR spectra with greater variance on the principal components.
std.NIR = sweep(NIRdata,2,sds,"/")

#Calculating principal components using prcomp and giving a summary of the principal components
fit = prcomp(std.NIR)
summary(fit)
#Even though there are 235 spectra, we notice that prcomp only returns 165 principal components.
#By inspecting the % of variance explained, we can see that 165 principal components account for ~100% of the variance and therefore the remaining 70 principal components are negligible. 

#Putting the standard deviations of the principal components into a data frame:
sdev = data.frame(fit$sdev)

#Converting the values for the standard deviations into variances.
variance = sdev^2
#This gives us a data frame with the variances of the principal components.

#Summing up the variances of the principal components to calculate the total variance.
sum_variance = sum(variance[1:165,])
sum_variance
#Observe that the summed variance is 235 which is what we should expect since the principal components account for the variance of 235 the standardized spectra.

#Initialising vector of 10 0s which will store the proportion of variance explained by the first 10 PCs.
prop_var = rep(0,10)

#Using for loop to set the kth entry in the vector to the cumulative variance of the first k PCs.
for(k in 1:10)
{
   prop_var[k] = sum(variance[1:k,])/sum_variance
}

#Plotting the cumulative variance for first 10 PCs.
plot(prop_var, xlab = "# of PCs", ylab ="Proportion of Variance")

#From this graph and the summary of the PCs I think 2 principal components are required to represent the NIR data.
#I think this because from the graph we see that there is a relatively large increase in cumulative explained variance between PCs 1,2. However after 2 PCs, the increase in cumulative explained variance is relatively small. 

#Question 4
#Calculating the PC scores for each sample using predict function. 
PCscores = predict(fit)

#Plotting PC scores against glucose values for each sample. 
plot(PCscores[,1],data[,237], xlab="PC1 score", ylab="Glucose Value")
plot(PCscores[,2],data[,237], xlab="PC2 score", ylab="Glucose Value")
#I can't notice any clear relationship between the glucose value and the PC scores.
#However there are some things worth noting.
#For example, from the plot we see that there is a cluster of samples with low glucose level(~3) that all have a high PC1 score of over ~14.
#From the 2nd plot, again there does not seem to be any clear relationship but we notice that all samples with PC2 score less than -10 have glucose levels of greater than 20.

#Plotting the PC scores against the ethanol values.
plot(PCscores[,1],data[,238], xlab="PC1 score", ylab="Ethanol Value")
plot(PCscores[,2],data[,238], xlab="PC2 score", ylab="Ethanol Value")
#From these plots, it seems that there is a trend for samples with very high ethanol values of over 80 to have a high PC1 score of over 12.
#When we plot the ethanol values against the PC2 scores, we notice that the graph has a semi- "upside down v" shape to it.
#The samples with the highest ethanol values have middle of the range PC2 scores of between about -5 and 0..
#Also all samples with ethanol values below ~40 have PC2 scores of greater than 0.  

#Question 5
#Taking the first 2 columns of the rotation matrix, representing the loadings of the PC1 and PC2.  
loadings_mat = fit$rotation[,1:2]

#Plotting the loadings:
plot(loadings_mat[,1],loadings_mat[,2],xlab="PC1 loadings",ylab="PC2 loadings")
#Observe that there is a circular structure to this plot.
#This indicates that there is a relationship between the PC1 loadings and PC2 loadings on the absorbance values for each spectra.  

#Question 7
#Creating test data set with 1/3 the samples from the data set.
#Taking a random sample instead of the 1st third of the samples in case samples are grouped together by grain type or some other variable.
set.seed(17419782)
index <- sample(1:nrow(data), size = (1/3)*nrow(data))
test_data <- data[index,]
train_data <- data[-index,]

#Installing pls packages in order to use pcr function.
install.packages("pls")
library(pls)

#Using pcr function to create a model with response variable as glucose and predictor variables as the NIR data.
#Setting scale=TRUE to ensure the NIR data is standardized. Setting number of PCs = 8 and performing leave out one validation(LOO). 
#LOO validation will run the model n times, where n is the number of samples, omitting a different sample each time and test the model on this omitted sample. This will give us an indication of how the model performs.
model = pcr(train_data$Glucose~., data = train_data[,2:236], scale=TRUE, ncomp = 8, validation = "LOO")

#Using summary function we can get several details about the model.
#For validation, we can see the square root of the mean squared error of the predicted glucose values(RMSEP) for each # of PCs.
#We can see the percentage of variance in the NIR data explained by the PCs.
#We can see the percentage of variance in the glucose data explained by the PCs.
summary(model)

#Plotting the RMSEP for each of the # of PCs to decide how many PCs to use.
validationplot(model)
#From the graph we see that the reduction in RMSEP is relatively large up until 4 PCs, after which there is a levelling off.
#We also can see that 4 components explains ~54% of the variance in the glucose data. 
#There is a relatively large increase in % variance explained between 3 and 4 components.
#However, after 4 components, each subsequent increase in % variance in glucose data explained is relatively small. 
#Therefore it seems 4 principal components is a good choice for our model.

#Calculating the predicted glucose values for our test data set.
glucose_predict <- predict(model,test_data,ncomp = 4)

#Plotting the glucose values predicted by our model and the actual values for the test data set.
plot(glucose_predict,test_data$Glucose, xlab = "Predicted Glucose Values", ylab="True Glucose Values")
#We notice that while the model is not perfect there is a positive correlation between the predicted values and true values.

#Checking the correlation coefficient between the results.
cor(glucose_predict,test_data$Glucose)

#Question 8
#Creating training data and testing NIR data sets.
train_NIR <- NIRdata[-index,]
test_NIR <- NIRdata[index,]

#Using prcomp function to calculate the principal components of the training NIR data.
#Dividing by standard deviation since scale in the pcr function was set to TRUE
sd1 <- apply(train_NIR,2,sd)
train_NIRst <- sweep(train_NIR,2,sd1,"/")
train <- prcomp(train_NIRst)

#Calculating PC scores of training NIR from first principles.
#Using 4 principal components to obtain results that are comparable to results from q.7.
v <- as.matrix(train_NIRst)
m <- as.matrix(train$rotation[,c(1,2,3,4)])
PCtrainscores <- v%*%m
#Creating linear model with response variable as the glucose data and predictor variables as the PC scores for components 1,2,3,4.
linear_model = lm(data[-index,237] ~ PCtrainscores)
linear_model

#Dividing the test data by the standard deviation of the train data variables in order for the test data to be appropriately scaled. 
v1 <- as.matrix(sweep(test_NIR,2,sd1,"/"))
#Calculating PC scores for test data from 1st principles.
PCtestscores = v1%*%m

#Initialising vector to store glucose values of the test samples.
glucose_test <-rep(0,55)
#Storing the coefficients from the linear model as a,b,c,d,e
a = linear_model$coefficients[1]
b = linear_model$coefficients[2]
c = linear_model$coefficients[3]
d = linear_model$coefficients[4]
e = linear_model$coefficients[5]

#Using for loop to calculate predicted glucose value from PC test scores and the linear model coefficients.
for(k in 1:55)
{
   glucose_test[k] = a+b*PCtestscores[k,1] +c*PCtestscores[k,2]+d*PCtestscores[k,3]+e*PCtestscores[k,4]
}

#Plotting predicted glucose from q.7 and predicted glucose from this question to check that they agree for the test samples.
plot(glucose_predict,glucose_test, xlab="Q.7",ylab="Q.8")
glucose_predict-glucose_test
#We notice that the values agree from the plot.
#Small difference in values may be as a result of rounding/truncation errors in the computing process.