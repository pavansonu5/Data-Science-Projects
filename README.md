# Data-Science-Projects
## Projects
1. [Bank Marketing Analysis using PySpark](https://github.com/pavansonu5/Data-Science-Projects/tree/main/Bank%20Marketing%20Analysis%20using%20PySpark)
#### Abstract
The purpose of marketing campaigns is to help organizations improve financial posture as well as to gain a competitive edge over their competitors. Organizations use direct marketing to target specific client segments by reaching out to them through remote communication centers. Direct marketing campaigns can be managed more efficiently where clients are reached through remote communication centers. This project focuses on collecting the marketing campaign data, visualizing the data using Tableau, and identifying and selecting appropriate analytical techniques for big data analysis. It involves the application of Apache Spark, specifically PySpark, which is compatible with the python programming language. Apache spark has the advantages of distributed computed frameworks and libraries for real-time, large-scale data processing. 
#### Implementation and Observations
* The Dataset contains 12870 rows and 17 columns, where each row represents client information. To make predictions realistic, Compared to original dataset[Moro et al., 2014] from a research paper, a variable named "duration time" was removed, as it directly affects the output. It is a binary classification problem; find whether the client subscribed to a term deposit.
* All the visualisations(univariate, bivariate analysis) were done in Tableau.
* Dataset is moderate imabalanced(approx 67%-NO and 33%-YES), additionally over-sampled and under-sampled analysis were perfomed(which are not required, just for comparision)
* Followed by some statistical analysis on the dataset(count, mean) on combination of variables.
* Created some custom functions
1. a pipeline of stringindexer, onehotencoder, standardscaler 
2. evaluation metrics
* 3 different Classification algorithms were used- logistic regression, random forest and gradient boosting algorithm. Suprisingly, after variable transformations, logistic regression with accuracy=75.1% performed better than random forest with accuracy=74.2% and gradient boosting outperformed remaing all algorithms with accuracy=78.1%.
2. [Diamond price predictions](https://github.com/pavansonu5/Data-Science-Projects/tree/main/Diamond-price-predictions)
#### Implementation and Observations
* Dataset contains 11 variables, which includes measurements, cut, color, clarity, depth, table.
* Created custom functions for visualization for both numeirc and categorical variables(these functions can be fit to any dataset)
* Dataset has no null values. From the boxplots, it can be observed some variables have extreme values(might be outliers).
* Categorical variables are converted to numeric data based on data description.
* Finally, applied 7 different machine learning algorithms along with cross-validation. Based on the results(metrics=MSE), xgb REGRESSOR has low MSE valye.
4. [Cognitive Vehicular Networks](https://github.com/pavansonu5/Data-Science-Projects/tree/main/Cognitive-Vehicular-Networks)
5. [MNIST using different models](https://github.com/pavansonu5/Data-Science-Projects/tree/main/MNIST_CNN)
6. [Spaceship Titanic Complete Analysis](https://github.com/pavansonu5/Data-Science-Projects/tree/main/Spaceship-Titanic)
7. [Temparature Time Series Forecasting](https://github.com/pavansonu5/Data-Science-Projects/tree/main/Temperature-Time-Series-Forecasting)
