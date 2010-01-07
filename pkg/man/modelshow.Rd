\name{modelshow}
\Rdversion{1.1}
\alias{modelshow}
\title{
Function to read in and display the log-ratio distributions produced by the function simdist.
}
\description{
This function plots multiple null distributions in the same window, with the actual log-ratio and the chi-squared distribution on the same axes. It also invisibly returns a data.frame for the actual log-ratio comparisons with an empirical p-value column added to it. This function is similar to \code{constrshow} except it works on the output of \code{simdist} to plot log-ratio null distributions for one-group unconstrained models.
}
\usage{
modelshow(label, x = NULL, reload = F, k = 2, datacol = "LR", modelcol = "model", nullmodelcol = "null_model", chicol = "p (chi squared)", figure = T, models = c("w", "g", "gm", "l", "lm"), breaks = 200, hcol = "black", vcol = "red", chidf = 1, rc = NULL, namesep = "_", fpsuffix = ".fp", header = 1, sep = "\t")
}
\arguments{
  \item{label}{
A character string representing the base name for the input files, such that for each model 'M' the suffix '_M' is appended to the \code{label}. If no local object with the resulting name is found, then an attempt is made to open a file by this name in the working directory and create the corresponding object from it.
}
  \item{x}{
Legacy/debug argument that will be removed in a future version. Do not use.
}
  \item{reload}{
Logical value that, when true, tells the function to look for and read input files even if objects for them already exist. By default set to FALSE.
}
  \item{k}{
Numeric multiplier for converting the distribution of log-ratios to a chi-squared distribution, for plotting the chi-squared overlay. Since the distribution is actually a mixture model between chi-squared and a point value of zero, lowering the value of this argument from the default of 2 may give a better correspondence between the ideal and actual distributions.
}
  \item{datacol}{
Character string giving the name of the column in the input files containing the data whose distributions are to be plotted. By default the value of this argument is 'LR'.
}
  \item{modelcol}{
Character string giving the name of the column in the input files containing the abbreviated model designation of each datapoint. By default the value of this argument is 'model'.
}
  \item{nullmodelcol}{
Character string giving the name of the column in the input files containing the abbreviated null model designation shared by each datapoint obtained for a particular simulation. By default the value of this argument is 'null_model'.
}
  \item{chicol}{
Character string giving the name of the column in the input files containing the p-value returned by a chi-squared test. By default the value of this argument is 'p (chi squared)'.
}
  \item{figure}{
Logical value indicating whether to plot the distributions (as opposed to running this function just to obtain the simulation-based p-values without plotting anything). Set to TRUE by default.
}
  \item{models}{
A character vector of models to try which can be any combination of: 'w','g','gm','l', or 'lm'. Default value is c('g','gm','l','lm'). Invalid models are silently ignored.
}
  \item{breaks}{
An integer which is passed to the \code{breaks} argument of the \code{hist} function and controls the size and number of columns in the histograms. Default value is 200.
}
  \item{hcol}{
An integer or character string which is passed to the \code{border} argument of the \code{hist} function and controls the color of the histograms. Default value is 'black'.
}
  \item{vcol}{
An integer or character string which determines the color of the chi-squared overlay and the non-simulated log-ratio (if it is close enough to 0 to be on the same horizontal scale as the simulated distribution). Default value is 'red'.
}
  \item{chidf}{
A number indicating how many degrees of freedom the chi-squared distribution should have. Default value is 1.
}
  \item{rc}{
A vector of two integers, the first setting the number of rows and the second, the number of columns in the figure that is plotted. The default value is NULL, which tells the function to choose reasonable values for both based on how many models are being plotted.
}
  \item{namesep}{
A character string with which to join the \code{label} to the \code{suffixes}. Default value is '_'.
}
  \item{fpsuffix}{
A character string to append to the data.frame containing the data from the original, non-simulated comparison. Default value is '.fp'.
}
  \item{header}{
Logical value passed to the \code{header} argument of the \code{read.table} function in order to specify what header to use for the data read in from the input files.
}
  \item{sep}{
Character passed to the \code{sep} argument of the \code{read.table} function in order to specify the character that separates the columns in the input file from each other.
}
}
\details{
This function generates a set of names from the \code{label} and \code{models} arguments. If objects by this name already exist in the working environment, histograms are produced from those objects. These histograms allow the user to visually compare each observed log-ratio to the distribution of log-ratios that would be observed if the null hypothesis was correct for that comparison. If objects by this name do not exist in the working environment, the function attempts to open files of the same name in the working directory, and if successful, creates the corresponding objects in the working environment and then plots the histograms. This function also extracts a small data.frame of just the output for the non-simulated model fit, saves it to an object called LABEL.SUFFIX (concatenation of the \code{label} and \code{suffix} arguments) if such an object doesn't already exist in the local environment, and outputs a copy of it silently.
}
\value{
A data.frame is silently returned that has the same columns as the one output by \code{findpars}, but with one additional column named 'emp.p' which contains the p-values obtained by comparing the observed log-ratio to the simulated distribution rather than to the ideal chi-squared distribution.
}
\references{
Pletcher,S.D., Khazaeli,A.A., and Curtsinger,J.W. (2000). Why do life spans differ? Partitioning mean longevity differences in terms of age-specific mortality parameters. Journals of Gerontology Series A-Biological Sciences and Medical Sciences 55, B381-B389
}
\note{
This function is similar but not to be confused with \code{constrshow}. This one takes a character string as the first argument, the other one takes a data.frame. This one is the display function for \code{simdist} and is used for selecting an unconstrained model that best fits a dataset; the other one is the display function for \code{empdist} and is used for identifying parameters that differ between two datasets once a model or models have already been chosen.
}
\author{
Alex F. Bokov
}
\seealso{
\code{\link{findpars}}, \code{\link{simdist}}
}
\keyword{ survival }
