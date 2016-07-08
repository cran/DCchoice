\name{oohbchoice}
\alias{oohbchoice}

\title{Parametric approach to analyze one-and-one-half-bound dichotomous choice contingent valuation data}
\description{This function analyzes one-and-one-half-bound dichotomous choice contingent 
            valuation (CV) data on the basis of the utility difference approach.}
\usage{
oohbchoice(formula, data, subset, na.action = na.omit, dist = "log-logistic",
           par = NULL, \dots)

}

\arguments{
    \item{formula}{an object of S3 class \code{"Formula"} and specifies the model structure.}
    \item{data}{a data frame containing the variables in the model formula.}
    \item{subset}{an optional vector specifying a subset of observations.}
    \item{na.action}{a function which indicates what should happen when the data contains \code{NA}s.}
    \item{dist}{a character string setting the error distribution in the model, which 
                takes one of \code{"logistic"}, \code{"normal"}, \code{"log-logistic"}, 
                \code{"log-normal"} or \code{"weibull"}.}
    \item{par}{a vector of initial parameters over which the optimization is carried out.}
    \item{...}{optional arguments. Currently not in use.}
}

\details{

One-and-one-half-bound dichotomous choice contingent valuation (OOHBDC-CV), which was developed by Cooper et al. (2002), is an intermediate CV format between single-bounded dichotomous choice (SBDC) CV format and double-bounded dichotomous choice (DBDC) CV format. 

On the basis of an example of environmental valuation study, we will explain differences in question format among SBDC-CV, DBDC-CV, and OOHBDC-CV below.

In any of three CV surveys, two situations are firstly explained to the respondents: the current situation and the improved situation where an environmental improvement plan is assumed to be implemented. Question following the explanation of situations differs according to CV types.

In an SBDC-CV survey, after the explanation of situation mentioned above, the respondents are asked whether they would be willing to pay a certain amount of money toward the implementation of the plan. Therefore, there are two possible responses to the SBDC-CV survey: "yes (agree)," and "no (disagree)." The amounts (bids) that respondents are requested to contribute toward the plan are listed in advance. Each respondent is requested to answer a question randomly assigned with one of the listed bids. 

In a DBDC-CV survey, the CV question consists of two stages: after answering the SBDC-CV style question mentioned above (the first stage), the respondents are also asked to answer an additional SBDC-CV style question (the second stage). The bid in the second stage varies according the response in the first stage: a higher bid is displayed in the second stage if the response in the first stage is "yes," whereas a lower bid is displayed when the response in the first stage is "no." Therefore, there are four possible responses to the DBDC-CV survey: "yes-yes" ("yes" in the both stages), "yes-no" ("yes" and "no" in the first and second stages, respectively), "no-yes" ("no" and "yes" in the first and second stages, respectively), and "no-no" ("no" in the both stages).

In the OOHBDC-CV survey, after answering the first SBDC-CV style question (the first stage), only respondents who satisfy certain conditions are requested to answer an additional SBDC-CV style question (the second stage). Details in the OOHBDC-CV survey are as follows:
Step 1) A list of bid ranges [BLj, BHj] (j = 1, 2, ..., J), where BLj < BHj, are decided: i.e., [BL1, BH1], [BL2, BH2], ..., and [BLJ, BHJ]. 
Step 2) One of the bid ranges is randomly presented to respondents (e.g., a bid range of [BL3, BH3] for j = 3).
Step 3) One of the two bids presented to the respondents is selected randomly (i.e., BL3 or BH3 in the case of step 2 example) and then the respondents are asked whether they would be willing to pay the amount of the bid selected (the first stage).
Step 4) The respondents are asked to answer the second stage if they satisfy either condition: (a) their answer in the first stage is "yes" when the lower bid is presented in the first stage, or (b) their answer in the first stage is "no" when the higher bid is presented in the first stage. Therefore, there are six possible responses to the OOHBDC-CV survey: "no", "yes-no", and "yes-yes" when the lower bid is shown in the first stage; and "yes", "no-yes", and "no-no" when the higher bid is shown in the first stage. Refer to Cooper et al. (2002) for detailed explanation of OOHBDC-CV, including the example CV questions.

The function \code{oohbchoice()} implements an analysis of OOHBDC-CV responses (data) on the basis of the utility difference approach (Hanemann, 1984). 

The function returns an object of S3 class \code{oohbchoice} (see below for details), which inherits from an S3 class \code{dbchoice}. The generic functions for the S3 class \code{dbchoice} such as \code{print()}, \code{summary()}, \code{vcov()}, \code{logLik()}, \code{plot()}, and \code{update()}, are available for the S3 class \code{oohbchoice}. In addition, the two functions \code{krCI()} and \code{bootCI()} are available to compute the confidence intervals for the estimates of willingness-to-pays (WTPs): \code{krCI()} implements simulations to construct empirical distributions of the WTP, while \code{bootCI()} carries out nonparametric bootstrapping (see the package \pkg{DCchoice} for details).

Although \code{oohbchoice()} has six arguments, a basic generic call to \code{oohbchoice()} is given as follows:

\code{oohbchoice(formula, data, dist = "log-logistic")}

The argument \code{formula} defines the response variables and covariates (see below for details on the formula). The argument \code{data} specifies the data frame containing the variables in the model. The argument \code{dist} sets the error distribution: one of \code{"logistic"}, \code{"normal"}, \code{"log-logistic"} (default value), \code{"log-normal"}, or \code{"weibull"} is available. The difference between normal and log-normal models or between logistic or log-logistic ones is how the bid variable is incorporated into the model to be estimated. For the Weibull model, the bid variable must be entered in the natural log. Therefore, the user must be careful in defining the model formula that is explained in details below. 

A typical structure of the formula for \code{oohbchoice()} is defined as follows: 

\code{R1 + R2 ~ (the names of the covariates) | BL + BH}

The formula is an object of S3 class \code{Formula} and specifies the model structure. It has to be written in a symbolic expression in R. The formula consists of three parts as follows.

The first part, the left-hand side of the tilde sign (\code{~}), must contain the response variables for the suggested prices in the first and the second stage of CV questions. In the example below, \code{R1} denotes a binary or two-level factor response variable for a bid in the first stage and \code{R2} for a bid in the second stage. \code{R1} contains \code{yes} or \code{no} to the bid in the first stage or \code{1} for \code{yes} and \code{0} for \code{no}. \code{R2} contains \code{yes}, \code{no}, \code{none} to the bid in the second stage or \code{1} for \code{yes}, \code{0} for \code{no}, and \code{-9} for \code{none}. The value of \code{none} (\code{-9}) means that the respondents have no second stage: the respondents are asked to answer the second stage question only if they satisfy either condition: (a) they answer \code{yes} in the first stage when the lower bid is shown in the first stage, or (b) they answer \code{no} in the first stage when the higher bid is shown in the first stage.

The covariates are defined in the second part of the formula in the place of (the names of the covariates). Each covariate is connected with the arithmetic operator \code{+} and (the names of the covariates) in the above syntax should be replaced with \code{var1 + var2} and the like. The plus sign is nothing to do with addition of the two variables in the symbolic expression. When the covariate contains only a constant term, a value of \code{1} is set as the covariate: \code{R1 + R2 ~ 1 | BL + BH} (see the examples section below) 

The last part of the formula starts after the vertical bar (\code{|}). The names of the two variables (\code{BL} and \code{BH}) containing suggested lower and higher prices in OOHBDC-CV question are specified in this part. The two variables are also connected with the arithmetic operator (\code{+}). 

According to the structure of the formula, a data set (data frame) consists of three parts. An example of the data set (first six rows) is as follows (\code{gender} and \code{age} are respondents' characteristics and assumed to be covariates): 

\tabular{rrrrrrr}{
    \code{id} \tab \code{R1}  \tab  \code{R2}  \tab \code{gender} \tab    \code{age} \tab  \code{BL} \tab  \code{BH} \cr
    \code{1}  \tab \code{no}  \tab  \code{none} \tab  \code{male} \tab    \code{51} \tab  \code{2000} \tab  \code{2500} \cr
    \code{2}  \tab \code{yes} \tab  \code{no} \tab  \code{male} \tab    \code{30} \tab  \code{500} \tab  \code{1000} \cr
    \code{3}  \tab \code{yes} \tab  \code{yes} \tab  \code{female} \tab  \code{25} \tab  \code{500} \tab  \code{1000} \cr
    \code{4}  \tab \code{yes} \tab  \code{none} \tab  \code{male} \tab    \code{48} \tab  \code{1000} \tab  \code{1500} \cr
    \code{5}  \tab \code{no}  \tab  \code{yes} \tab  \code{male} \tab    \code{60} \tab  \code{1000} \tab \code{1500} \cr
    \code{6}  \tab \code{no}  \tab  \code{no} \tab  \code{female} \tab  \code{34} \tab  \code{2500} \tab  \code{3000}
}

Respondent 1 faced a bid range [\code{2000}, \code{2500}]; respondents 2 and 3 faced a bid range [\code{500}, \code{1000}]; respondents 4 and 5 faced a bid range [\code{1000}, \code{1500}]; and respondent 6 faced [\code{2500}, \code{3000}]. Respondent 1 answered \code{no} in the first stage of CV question and had no the second stage; respondent 2 answered \code{yes} and \code{no} in the first and second stage, respectively; respondent 3 answered \code{yes} and \code{yes} in the both stages; respondent 4 answered \code{yes} in the first stage and had no the second stage; respondent 5 answered \code{no} and \code{yes} in the first and second stage; and respondent 6 answered \code{no} in the both stages.

Note that \code{BL} and \code{BH} are NOT the first stage bid and the second stage bid, respectively. The function \code{oohbchoice()} understands which bids (\code{BL} and \code{BH}) are presented in the first stage and second stage, respectively, on the basis of values of variables \code{R1} and \code{R2}.


Nonparametric analysis of OOHBDC-CV data can be done by \code{\link{turnbull.oohb}}. 

}

\value{

The function returns an S3 class object \code{oohbchoice}, which inherits from the S3 class \code{dbchoice}. See \code{dbchoice()} for the details on the S3 object \code{dbchoice}.

}

\section{Acknowledgments}{
We would like to thank Dr. Joseph C. Cooper and Dr. Giovanni Signorello for their kindness.
}

\references{

Cooper JC, Hanemann M, Signorello G (2002).
\dQuote{One-and-one-half-bound dichotomous choice contingent valuation},
\emph{The Review of Economics and Statistics},
\bold{84}, 742--750.

Hanemann WM (1984).
\dQuote{Welfare Evaluations in Contingent Valuation Experiments with Discrete Responses},
\emph{American Journal of Agricultural Economics},
\bold{66}(2), 332--341.

}

\seealso{
\code{\link{summary.dbchoice}}, \code{\link{oohbsyn}}, \code{\link{krCI}}, \code{\link{bootCI}}, 
\code{\link{turnbull.oohb}}, \code{\link[Formula]{Formula}}
}

\examples{
## See oohbsyn.
}




\keyword{DCchoice}
\keyword{double-bounded}
\keyword{nonlinear}

