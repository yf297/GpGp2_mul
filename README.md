# GpGp2_mul
Enter directory and write "make".
The demo.R file in the R directory fits two models as examples.\
Do not load GpGp or GpGpm, since some functions have name overlap.\
The fit_bivariate_matern function assumes mean 0 responses, so y - mean(y) is done before passing in.\
After saving estimated paramaters, use predict function from GpGp to make predictions. Account for the y - mean(y).
