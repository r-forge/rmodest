\name{constrshow}
\Rdversion{1.1}
\alias{constrshow}
\title{
Function to read in and display the log-ratio distributions produced by the function empdist.
}
\description{
This function plots multiple null distributions in the same window, with the actual log-ratio and the chi-squared distribution on the same axes. It also invisibly returns a data.frame for the actual log-ratio comparisons with an empirical p-value column added to it. This function is similar to \code{modelshow} except it works on the output of \code{empdist} to plot log-ratio null distributions for two-group models with a constrained parameter.
}
\usage{
constrshow(d, x = NULL, k = 2, datacol = "LR", idcol = "id", idobs = 1, modelcol = "model", nullmodelcol = "null_model", chicol = "p (chi squared)", figure = T, models = c("wa", "wb", "ga", "gb"), breaks = 200, hcol = "black", vcol = "red", chidf = 1, rc = NULL)
}
\arguments{
  \item{d}{
A data.frame as read from the output file produced by \code{empdist}.
}
  \item{x}{
Legacy/debug argument that will be removed in a future version. Do not use.
}
  \item{k}{
Numeric multiplier for converting the distribution of log-ratios to a chi-squared distribution, for plotting the chi-squared overlay.
}
  \item{datacol}{
Character string giving the name of the column containing the data whose distributions are to be plotted. By default the value of this argument is 'LR'.
}
  \item{idcol}{
Character string giving the name of the column containing the id of the permutation or sample. By default the value of this argument is 'id'.
}
  \item{idobs}{
An integer indicating which id corresponds to the original model fits.
}
  \item{modelcol}{
Character string giving the name of the column containing the abbreviated model designation of each datapoint. By default the value of this argument is 'model'.
}
  \item{nullmodelcol}{
Character string giving the name of the column containing the abbreviated null model designation shared by each datapoint obtained for a particular simulation. By default the value of this argument is 'null_model'.
}
  \item{chicol}{
Character string giving the name of the column containing the p-value returned by a chi-squared test. By default the value of this argument is 'p (chi squared)'.
}
  \item{figure}{
Logical value indicating whether to plot the distributions (as opposed to running this function just to obtain the bootstrapped p-values without plotting anything). Set to TRUE by default.
}
  \item{models}{
A character vector of models to try which can be any combination of: 'wa','wb','ga','gb','gma','gmb','gmc','la','lb','ls','lma','lmb','lmc', or 'lms'. Default value is c('ga','gb').
}
  \item{breaks}{
An integer which is passed to the \code{breaks} argument of the \code{hist} function and controls the size and number of columns in the histograms. Default value is 200.
}
  \item{hcol}{
An integer or character string which is passed to the \code{border} argument of the \code{hist} function and controls the color of the histograms. Default value is 'black'.
}
  \item{vcol}{
An integer or character string which determines the color of the chi-squared overlay and the non-bootstrapped log-ratio (if it is close enough to 0 to be on the same horizontal scale as the distribution). Default value is 'red'.
}
  \item{chidf}{
A number indicating how many degrees of freedom the chi-squared distribution should have. Default value is 1.
}
  \item{rc}{
A vector of two integers, the first setting the number of rows and the second, the number of columns in the figure that is plotted. The default value is NULL, which tells the function to choose reasonable values for both based on how many models are being plotted.
}
}
\details{
This function plots a collection of histograms are produced the data.frame that is passed to the first argument. These histograms allow the user to visually compare each observed log-ratio to the distribution of log-ratios that would be observed if the null hypothesis was correct for that comparison. This function also extracts a data.frame of just the output for the non-simulated model fit and returns it silently.
}
\value{
A data.frame is silently returned that has the same columns as the one output by \code{findpars}, but with one additional column named 'emp.p' which contains the p-values obtained by comparing the observed log-ratio to the bootstrapped distribution rather than to the ideal chi-squared distribution.
}
\references{
Pletcher,S.D., Khazaeli,A.A., and Curtsinger,J.W. (2000). Why do life spans differ? Partitioning mean longevity differences in terms of age-specific mortality parameters. Journals of Gerontology Series A-Biological Sciences and Medical Sciences 55, B381-B389
}
\note{
This function is similar but not to be confused with \code{modelshow}. This one takes a character string as the first argument, the other one takes a data.frame. This one is the display function for \code{empdist} and is used for identifying parameters that differ between two datasets once a model or models have already been chosen; the other one is the display function for \code{simdist} and is used for selecting an unconstrained model that best fits a dataset.
}
\author{
Alex F. Bokov
}
\seealso{
\code{\link{findpars}}, \code{\link{empdist}}
}
\section{Warning}{
These examples may take 10 minutes or more to finish running.
}
\examples{
data(demo1);
data(demo2);
empdist(demo1,demo2,'test',rounds=100);
test.in<-read.table(file='test.txt',sep='\t',header=T,check.names=F)
test.table<-constrshow(test.in);
test.table[1:10,];
names(test.table);
# Invoking constrshow with a different layout and with models explicitly specified.
constrshow(test.in,models=c('ga','lma','lms'),rc=c(2,2));
}
\keyword{ survival }
