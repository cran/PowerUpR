# minimum detectable effect size
.mdes.fun <- function(power, alpha, sse, df, two.tailed){
  t1 <- ifelse(two.tailed == TRUE, abs(qt(alpha / 2, df)), abs(qt(alpha, df)))
  t2 <- abs(qt(power, df))
  m <- ifelse(power >= 0.5, t1 + t2, t1 - t2)
  mdes <- m * sse
  lcl <- mdes * (1 - t1 / m)
  ucl <- mdes * (1 + t1 / m)
  mlu <- cbind(mdes, lcl, ucl)
  return(mlu)
}

# statistical power
.power.fun <- function(es, alpha, sse, df, two.tailed){
  lambda <- es/sse
  power <- ifelse(two.tailed == FALSE,
                  1 - pt(qt(alpha, df, lower.tail = FALSE), df, lambda),
                  1 - pt(qt(alpha / 2, df, lower.tail = FALSE), df, lambda) +
                    pt(-qt(alpha / 2, df, lower.tail = FALSE), df, lambda))
  return(power)
}


# summarize mdes output
.summ.mdes <- function(power, alpha, sse, df, two.tailed, mdes) {
  cat("\n Minimum detectable effect size: \n --------------------------------------- \n MDES is ",
      round(mdes[1], 3), " ", 100 * (1 - round(alpha, 2)), "% CI [", round(mdes[2], 3),
      ",", round(mdes[3], 3), "] with ", round(power,3)*100,
      "% power \n ---------------------------------------\n Degrees of freedom:", df,
      "\n Standardized standard error:", round(sse, 3), "\n Type I error rate:", alpha,
      "\n Type II error rate:", round(1 - power, 3), "\n Two-tailed test:", two.tailed, "\n",
      sep = "")
}

# summarize power output
.summ.power <- function(es, alpha, sse, df, two.tailed, power) {
  mlu <- .mdes.fun(power = power, alpha = alpha, sse = sse, df = df, two.tailed = two.tailed)
  cat("\n Statistical power: \n --------------------------------------- \n ",
      round(power, 3) * 100, "% power to detect an ES of ", round(mlu[1], 3), " ",
      100 * (1 - round(alpha, 2)), "% CI [", round(mlu[2], 3), ",", round(mlu[3],3),
      "] \n --------------------------------------- \n Degrees of freedom:", df,
      "\n Standardized standard error:", round(sse, 3), "\n Type I error rate:", alpha,
      "\n Type II error rate:", round(1 - power, 3), "\n Two-tailed test:", two.tailed, "\n",
      sep = "")
}
