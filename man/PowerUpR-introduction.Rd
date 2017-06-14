
\name{introduction}
\alias{PowerUpR}
\title{Introduction to \code{PowerUpR}}
\description{
  \code{PowerUpR} is an implementation of \emph{PowerUp!}, and \emph{PowerUp!-Moderator} in R environment (R Core Team, 2017). \emph{PowerUp!} series consist of convenient excel based functions to conduct statistical power analysis for various experimental and quasi-experimental designs (Dong & Maynard, 2013).
  It also enables users to conduct statistical power analysis for moderator effects in two- and three-level cluster randomized trials (Dong, Kelcey, & Spybrook, in press; Spybrook, Kelcey, & Dong, in press), and for mediator effects in two-level cluster randomized trials (Kelcey, Dong, Spybrook, & Cox, 2017; Kelcey, Dong, Spybrook, & Shen, in press).

  The \code{PowerUpR} package bases its framework on three fundemental concepts in statistical power analysis; statistical power calculation, minimum detectable effect size calculation, and constrained optimal sample allocation (COSA; Hedges & Borenstein, 2014; Raudenbush, 1997; Raudenbush & Liu, 2000).
  COSA problems can be solved in the following forms,
  (i) under budgetary constraints given marginal costs per unit,
  (ii) under power constraints given marginal costs per unit,
  (iii) under MDES constraints given marginal costs per unit, and
  (iv) under sample size constraints for one or more levels along with any of the i, ii or iii options.

A design parameter (one of the MDES, power, or COSA) can be requested by using approriate function given design characteristics.

Each function begins with an \strong{output} name, followed by a period,  an \strong{effect} (for moderator only), followed by a period, and a \strong{design} name.

There are three types of output; \code{mdes},  \code{power}, and \code{optimal}, and 14 types of design; \code{ira1r1}, \code{bira2r1}, \code{bira2f1}, \code{bira2c1}, \code{cra2r2}, \code{bira3r1}, \code{bcra3r2}, \code{bcra3f2}, \code{cra3r3},
\code{bira4r1}, \code{bcra4r2}, \code{bcra4r3}, \code{bcra4f3}, and \code{cra4r4}. The first three or four letters of the design stands
for the type of assignment, for individual random assignment \code{ira}, for blocked individual random assignment \code{bira},
for cluster random assignment \code{cra}, and for blocked cluster random assignment \code{bcra}.
It is followed by a number indicating number of levels. A single letter followed by a number indicates whether a block is considered
to be \code{r}, random; \code{f}, fixed; or \code{c}, constant and the level at which random assingment takes place.

A function with two keywords (output.design) returns results for main effect. For moderator effect,
there is an additional keyword so the function has three keywords (output.effect.design). There are three moderator effects available for design \code{cra2r2} and five moderator effect available for design \code{cra3r3}. In total there are six moderator effects; \code{mod1n}, \code{mod1r}, \code{mod2}, \code{mod2n}, \code{mod2r}, and \code{mod3}.The number and the single letter at the end stands for the level of moderator variable and whether it varies randomly or non-randomly across higher level unit.


So, to find MDES for main effect in two-level cluster randomized design where random assignment is at level-2, function \code{mdes.cra2r2} is used. Similiarly, to find MDES for a non-randoly varying moderator effect at level-1 for the same design, function \code{mdes.mod1n.cra2r2} is used.

Each function requires slightly different arguments depending on the output it produces and the design. Most of the arguments have default values to provide users a starting point, which can be found in \emph{usage} section of the documentation. For all functions default values are
\itemize{
  \item{\code{mdes} = .25}
  \item{\code{power} (\eqn{1-\beta}) = .80}
  \item{\code{alpha} (\eqn{\alpha}) = .05}
  \item{\code{two.tail} = \code{TRUE}}
  \item{\code{P} = .50}
}
and depending on the effect and design
\itemize{
  \item{any of one of \code{g1, g2, g3, g4} = 0}
  \item{any sequence of \code{R12, R22, R32, R42} = 0}
  \item{any sequence of \code{RT22, RT32, RT42} = 0}
  \item{\code{Q} = \code{NULL}, implies continuous moderator}
}
Users should be aware of default values and change them if necessary. Depending on the function minimum required arguments are
\itemize{
  \item{any sequence of \code{rho2, rho3, rho4}}
  \item{any sequence of \code{omega2, omega3, omega4}}
  \item{any one of, any sequence of, or any combination of \code{n, J, K, L}}
}
For definition of above-mentioned parameters, statistical models and formulas see Dong and Maynard (2013), Dong, Kelcey, and Spybrook (in press), Spybrook, Kelcey, and Dong (in press), Kelcey, Dong, Spybrook, and Cox (2017), and  Kelcey, Dong, Spybrook, and Shen (in press).

For reference intraclass correlation (\code{rho2, rho3})
values see Dong, Reinke, Herman, Bradshaw, and Murray (2016), Hedberg and Hedges (2014), Hedges and Hedberg (2007, 2013), Kelcey, and Phelps (2013), Schochet (2008), Spybrook, Westine, and Taylor (2016). For reference variance (\code{R12, R22, R32}) values see Bloom, Richburg-Hayes, and Black (2007), Deke et al. (2010), Dong et al. (2016), Hedges and Hedberg (2013), Kelcey, and Phelps (2013), Spybrook, Westine,and Taylor (2016), Westine, Spybrook, and Taylor (2013).
Users can also obtain design parameters for various levels using publicly available state or district data.

}


\references{
Bloom, H. S., Richburg- Hayes, L. & Black, A. R. (2007).  Using Covariates to Improve Precision for Studies that Randomize Schools to Evaluate Educational Interventions.  \emph{Educational Evaluation and Policy Analysis, 29(1)}, 0-59.

Deke, John, Dragoset, Lisa, and Moore, Ravaris (2010). Precision Gains from Publically Available School Proficiency Measures Compared to Study-Collected Test Scores in Education Cluster-Randomized Trials (NCEE 2010-4003). Washington, DC: National Center for Education Evaluation and Regional Assistance, Institute of Education Sciences, U.S. Department of Education.  \url{http://ies.ed.gov/ncee/pubs/20104003/}

Dong, N., Kelcey, B., & Spybrook, J. (in press). Power analyses of moderator effects in threelevel
cluster randomized trials. \emph{Journal of Experimental Education}.

Dong, N., & Maynard, R. A. (2013). PowerUp!: A Tool for Calculating Minum Detectable Effect Sizes and Minimum Required Sample Sizes for Experimental and Quasi-Experimental Design Studies,\emph{Journal of Research on Educational Effectiveness, 6(1)}, 24-6.

Dong, N., Reinke, W. M., Herman, K. C., Bradshaw, C. P., & Murray, D. W. (2016). Meaningful effect sizes, intraclass correlations, and proportions of variance explained by covariates for panning two-and three-level cluster randomized trials of social and behavioral outcomes. \emph{Evaluation Review}. doi: 10.1177/0193841X16671283

Hedges, L. V., & Borenstein, M. (2014). Conditional Optimal Design in Three- and Four-Level Experiments. \emph{Journal of Educational and Behavioral Statistics, 39(4)}, 257-281

Hedberg, E., & Hedges, L. V.(2014). Reference Values of Within-District Intraclass Correlations of Academic Achivement by District Characteristics: Results From a Meta-Analysis of District-Specified Values. \emph{Evaluation Review, 38(6)}, 546-582.

Hedges, L. V., & Hedberg, E. (2007). Interclass correlation values for planning group-randomized trials in education. \emph{Educational Evaluation and Policy Analysis, 29(1)}, 60-87.

Hedges, L. V., & Hedberg, E. (2013). Interclass Correlations and Covariate Outcome Correlations for Planning Two- and Three-Level Cluster-Randomized Experiments in Education. \emph{Evaluation Review, 37(6)}, 445-489.

Kelcey, B., Dong, N., Spybrook, J., & Cox, K. (2017). Statistical Power for Causally Defined Indirect Effects in Group-Randomized Trials With Individual-Level Mediators. \emph{Journal of Educational and Behavioral Statistics}, 1076998617695506.

Kelcey, B., Dong, N., Spybrook, J., & Shen, Z. (in press). Statistical power for causally-defined
mediation in group-randomized studies. \emph{Multivariate Behavioral Research}.

Kelcey, B., & Phelps, G. (2013). Strategies for improving power in school randomized studies of professional development. \emph{Evaluation Review, 37(6)}, 520-554.

R Core Team (2017). R: A language and environment for statistical computin . R Foundation for Statistical Computing, Vienna, Austria. \url{https://www.R-project.org}.

Raudenbush, S. W. (1997). Statistical analysis and optimal design for cluster randomized trials. \emph{Psychological Methods, 2}, 173-185.

Raudenbush, S. W., & Liu, X. (2000). Statistical power and optimal design for multisite trials. \emph{Psychological Methods, 5}, 199-213.

Schochet, P. Z. (2008). Statistical Power for Random Assignment Evaluations of Education Programs. \emph{Journal of Educational and Behavioral Statistics, 33(1)}, 62-87

Spybrook, J., Kelcey, B., & Dong, N. (2016). Power for detecting treatment by moderator effects in two-and three-level cluster randomized trials. \emph{Journal of Educational and Behavioral Statistics, 41}(6), 605-627.

Spybrook, J., Westine, C. D., & Taylor, J. A. (2016). Design Parameters for Impact Research in Science Education: A Multisite Anlaysis. \emph{AERA Open, 2(1)}, 1-15.

Westine, C. D., Spybrook, J.,  & Taylor, J. A. (2013). An Empirical Investigation of Variance Design Parameters for Planning Cluster-Randomized Trials of Science Achievement. \emph{Evaluation Review, 37(6)}, 490-519.

}

\keyword{introduction}
