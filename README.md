To install CRAN version
```{r}
install.packages("PowerUpR")
```

Statistical power, minimum detectable effect size (MDES), or minimum required sample size (MRSS) can be requested by using the relevant function given design parameters. Each function begins with an **output** name, following by a period, and a **design** name. There are three types of output; `mdes`,  `power`, and `mrss`, and 14 types of design; `ira1r1`, `bira2r1`, `bira2f1`, `bira2c1`, `cra2r2`, `bira3r1`, `bcra3r2`, `bcra3f2`, `cra3r3`, `bira4r1`, `bcra4r2`, `bcra4r3`, `bcra4f3`, and `cra4r4`. The first three letters of the design stands for the type of assignment, for individual random assignment `ira`, for blocked individual random assignment `bira`, for cluster random assignment `cra`, and for blocked cluster random assignment `bcra`.
It is followed by a number indicating number of levels. A single letter followed by a number indicates whether the top level block is considered to be `r`, random; `f`, fixed; or `c`, constant and the level of the random assingment. For example, to find MDES for three-level blocked cluster-level randomized design where random assignment is at level 2, the function `mdes.bcra3r2` can be called.

Each function requires slightly different arguments depending on the output it produces and the design. Most of the arguments have default values to provide users a starting point. Default values are

- `es` = .25
- `power` = .80
- `alpha` = .05
- `two.tailed` = `TRUE`
- `p` = .50
- `g1`, `g2`, `g3`, `g4` = 0
- any sequence of `r21`, `r22`, `r23`, `r24` = 0
- any sequence of `r2t2`, `r2t3`, `r2t4` = 0

Users should be aware of default values and change them if necessary. 
Minimum required arguments to successfully run a function are

- any sequence of `rho2`, `rho3`, `rho4`
- any sequence of `omega2`, `omega3`, `omega4`
- any sequence of `n`, `J`, `K`, `L`

For definition of above-mentioned parameters see Dong & Maynard (2013) and Hedges & Rhoads (2009), or help files. For reference intraclass correlation (`rho2`, `rho3`) values see Dong, Reinke, Herman, Bradshaw, and Murray (2016), Hedberg and Hedges (2014), Hedges and Hedberg (2007, 2013),  Kelcey, and Phelps (2013), Schochet (2008), Spybrook, Westine, and Taylor (2016). For reference variance (`r21`, `r22`, `r23`) values see Bloom, Richburg-Hayes, and Black (2007), Deke et al. (2010), Dong et al. (2016), Hedges and Hedberg (2013),  Kelcey, and Phelps (2013), Spybrook, Westine,and Taylor (2016), Westine, Spybrook, and Taylor (2013). Users can also obtain design parameters for various levels using publicly available state or district data.

Please email us any issues or suggestions.

Metin Bulus bulus.metin@gmail.com  
Nianbo Dong dong.nianbo@gmail.com  

**Suggested citation:**  

Bulus, M., & Dong, N. (2018).  `PowerUpR`: R version of PowerUp!. R package version 0.2.3.

Dong, N., & Maynard, R. A. (2013). PowerUp!: A Tool for Calculating Minimum Detectable Effect Sizes and Minimum Required Sample Sizes
for Experimental and Quasi-Experimental Design Studies, *Journal of Research on Educational Effectiveness, 6(1)*, 24-6.

## References
Bloom, H. S., Richburg- Hayes, L. & Black, A. R. (2007).
Using Covariates to Improve Precision for Studies that Randomize Schools to Evaluate Educational Interventions.
*Educational Evaluation and Policy Analysis, 29(1)*, 0-59.

Deke, John, Dragoset, Lisa, and Moore, Ravaris (2010). Precision Gains from Publically Available School Proficiency Measures Compared to Study-Collected Test Scores in Education Cluster-Randomized Trials (NCEE 2010-4003). Washington, DC: National Center for Education Evaluation and Regional Assistance, Institute of Education Sciences, U.S. Department of Education. Retrieved from https://ies.ed.gov/ncee/pubs/20104003/

Dong, N., & Maynard, R. A. (2013). PowerUp!: A Tool for Calculating Minimum Detectable Effect Sizes and Minimum Required Sample Sizes
for Experimental and Quasi-Experimental Design Studies, *Journal of Research on Educational Effectiveness, 6(1)*, 24-6.

Dong, N., Reinke, W. M., Herman, K. C., Bradshaw, C. P., & Murray, D. W. (2016). Meaningful effect sizes, intraclass correlations, and proportions of variance explained by covariates for panning two-and three-level cluster randomized trials of social and behavioral outcomes. *Evaluation Review*. doi: 10.1177/0193841X16671283

Hedges, L. V., & Borenstein, M. (2014). Conditional Optimal Design in Three- and Four-Level Experiments.
*Journal of Educational and Behavioral Statistics, 39(4)*, 257-281.

Hedberg, E., & Hedges, L. V.(2014). Reference Values of Within-District Intraclass Correlations of Academic Achivement
by District Characteristics: Results From a Meta-Analysis of District-Specified Values. *Evaluation Review, 38(6)*, 546-582.

Hedges, L. V., & Hedberg, E. (2007). Interclass correlation values for planning group-randomized trials in education.
*Educational Evaluation and Policy Analysis, 29(1)*, 60-87.

Hedges, L. V., & Hedberg, E. (2013). Interclass Correlations and Covariate Outcome Correlations for Planning 
Two- and Three-Level Cluster-Randomized Experiments in Education. *Evaluation Review, 37(6)*, 445-489.

Hedges, L. & Rhoads, C.(2009). Statistical Power Analysis in Education Research (NCSER 2010-3006).
Washington, DC: National Center for Special Education Researc , Institute of Education Sciences, U.S. Department of Education.
Retrieved from https://ies.ed.gov/ncser/.

Kelcey, B., & Phelps, G. (2013). Strategies for improving power in school randomized studies of professional development. *Evaluation Review, 37(6)*, 520-554.

R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.r-project.org/

Raudenbush, S. W. (1997). Statistical analysis and optimal design for cluster randomized trials.
*Psychological Methods, 2*, 173-185.

Raudenbush, S. W., & Liu, X. (2000). Statistical power and optimal design for multisite trials.
*Psychological Methods, 5*, 199-213.

Schochet, P. Z. (2008). Statistical Power for Random Assignment Evaluations of Education Programs.
*Journal of Educational and Behavioral Statistics, 33(1)*, 62-87.

Spybrook, J., Westine, C. D., & Taylor, J. A. (2016). Design Parameters for Impact Research in Science Education:
A Multisite Anlaysis. *AERA Open, 2(1)*, 1-15.

Westine, C. D., Spybrook, J.,  & Taylor, J. A. (2013). An Empirical Investigation of Variance Design Parameters
for Planning Cluster-Randomized Trials 
