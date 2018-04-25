<h2> Power Analysis Tools for Multilevel Randomized Experiments </h2>
<h3> R version of PowerUp! </h3>

To install and load the CRAN version
```{r}
install.packages("PowerUpR")
library(PowerUpR)
```

Statistical power, minimum detectable effect size (MDES), or minimum required sample size (MRSS) can be requested by using the relevant function given design parameters. Each function begins with an **output** name, following by a period, and a **design** name. There are three types of output; `mdes`,  `power`, and `mrss`. Each output can be requested for fourteen types of design for average treatment effect and five types of design for moderation effects; `ira1r1`, `bira2r1`, `bira2f1`, `bira2c1`, `cra2r2`, `bira3r1`, `bcra3r2`, `bcra3f2`, `cra3r3`, `bira4r1`, `bcra4r2`, `bcra4r3`, `bcra4f3`, `cra4r4`, `mod221`, `mod222`, `mod331`, `mod332`, and `mod333`. For mediation effect only `mdes` and `power` can be requested for two types of design; `med211` and `med221`. The first three letters of the design for average treatment effects stands for the type of assignment, for individual random assignment `ira`, for blocked individual random assignment `bira`, for cluster random assignment `cra`, and for blocked cluster random assignment `bcra`. Numbers indicate the total number of levels and the level at which randomization take place correspondingly. Single letter inbetween refers to whether the top level is random or fixed. Naming conventions are slighlty different for moderation and mediation effects. Numbers following `mod` keyword indicates the total number of levels, the level at which randomization takes place and the level at which the moderator resides. As for the mediation effects, numbers following `med` keyword indicates the level at which path `a`, `b` and `cp` resides. 
For example, to find MDES for three-level blocked cluster-level randomized design where random assignment is at level 2, the function `mdes.bcra3r2` can be called.


**Suggested citations**:

Dong, N., Kelcey, B., Spybrook, J., & Maynard, R. A. (2017a). PowerUp!-Moderator: A tool for calculating statistical power and minimum detectable effect size of the moderator effects in cluster randomized trials (Version 1.08) [Software]. Available from [http://www.causalevaluation.org/](http://www.causalevaluation.org/)

Dong, N., Kelcey, B., Spybrook, J., & Maynard, R. A. (2017b). PowerUp!-Mediator: A tool for calculating statistical power for causally-defined mediation in cluster randomized trials. (Beta Version 1.0) [Software]. Available from [http://www.causalevaluation.org/](http://www.causalevaluation.org/)

Dong, N.,  Kelcey, B., & Spybrook, J. (2017). Power analyses of moderator effects in three-level cluster randomized trials. *Journal of Experimental Education*. Advance online publication. doi: 10.1080/00220973.2017.1315714

Dong, N. & Maynard, R. A. (2013). PowerUp!: A tool for calculating minimum detectable effect sizes and sample size requirements for experimental and quasi-experimental designs. *Journal of Research on Educational Effectiveness*, 6(1), 24-67.  doi: 10.1080/19345747.2012.673143

Dong, N., & Maynard, R. A. (2013). PowerUp!: A tool for calculating minimum detectable effect sizes and minimum required sample sizes for experimental and quasi-experimental design studies. [Software]. [http://www.causalevaluation.org/](http://www.causalevaluation.org/)

Kelcey, B., Dong, N., Spybrook, J., & Shen, Z. (2017). Experimental Power for Indirect Effects in Group-randomized Studies with Group-level Mediators. *Multivariate Behavioral Research*. Advance online publication. doi: 10.1080/00273171.2017.1356212

Kelcey, B., Dong, N., Spybrook, J., & Cox, K. (2017). Statistical power for causally-defined individual and contextual indirect effects in group-randomized Trials. *Journal of Educational and Behavioral Statistics*. Advance online publication. doi: 10.3102/1076998617695506

Spybrook, J., Kelcey, B., & Dong, N. (2016). Power for detecting treatment by moderator effects in two and three-level cluster randomized trials. *Journal of Educational and Behavioral Statistics*. doi: 10.3102/1076998616655442


**Acknowledgement**:

This work is supported by National Science Foundation through a collaborative research grant titiled “Power Analyses for Moderator and Mediator Effects in Cluster Randomized Trials” to Benjamin Kelcey (Award Number: 1437679), Jessaca Spybrook (Award Number:1437692). and Nianbo Dong (Award Number: 1437745).


**Bug Reports, Feedback, and Suggestions**:

Please contact individual authors or maintainer for any issues or suggestions. <br>
Metin Bulus [maintainer] bulusmetin@gmail.com  <br>
Nianbo Dong [author] <br>
Benjamin Kelcey [author] <br>
Jessaca Spybrook [author] <br>
