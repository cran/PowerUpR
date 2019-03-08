<h2> Power Analysis Tools for Multilevel Randomized Experiments </h2>
<h3> R version of PowerUp! Excel Series </h3>

To install and load the library
```{r}
install.packages("PowerUpR")
library(PowerUpR)
```

Statistical power, minimum detectable effect size (MDES), MDES difference (MDESD), or minimum required sample size (MRSS) can be requested by using the relevant function given design parameters. Each function begins with an **output** name, follows by a period, and ends with a **design** name in the form `<output>.<design>()`. There are three types of output; `mdes` for main effects (`mdes` or `mdesd` for moderation effects),  `power`, and `mrss`. Each output can be requested for fourteen types of designs to detect main treatment effects; `ira1r1`, `bira2r1`, `bira2f1`, `bira2c1`, `cra2r2`, `bira3r1`, `bcra3r2`, `bcra3f2`, `cra3r3`, `bira4r1`, `bcra4r2`, `bcra4r3`, `bcra4f3`, `cra4r4`, and five types of designs to detect moderator effects; `mod221`, `mod222`, `mod331`, `mod332`, and `mod333`. To detect mediator effects, only `power` can be requested for two types of designs; `med211` and `med221`.

For designs to detect main effects, first three letters stands for the type of assignment; for individual random assignment `ira`, for blocked individual random assignment `bira`, for cluster random assignment `cra`, and for blocked cluster random assignment `bcra`. Numbers indicate total number of levels and the level at which randomization takes place correspondingly. The single letter inbetween refers to whether the top level is random or fixed. Naming conventions are slighlty different for designs to detect moderator and mediator effects. Numbers following `mod` keyword indicate total number of levels, the level at which randomization takes place, and the level at which the moderator resides correspondingly. As for the mediator effects, numbers following `med` keyword indicate the level at which path `a`, `b` and `cp` resides. 

For example, the function `mdes.cra2r2()` can be called to calculate MDES for main treatment effect in a two-level cluster-randomized trial. Similiarly, the function `mdesd.mod222()` can be called to calculate MDESD for moderator effect that resides at level 2 in a two-level cluster-randomized trial. Finally, the function `power.med221()` can be called to calculate statistical power for mediator effect that resides at level 2 in a two-level cluster-randomized trial. 

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



