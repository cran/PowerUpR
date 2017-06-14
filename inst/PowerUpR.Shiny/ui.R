shinyUI({
  fluidPage(
    fluidRow(style="background-color:rgb(200,200,200)",
      column(width=4,
             helpText(h2("Output"), style="color:white"),
             selectInput('output', 'Request',
                         choices = c(`Statistical Power` = 'power',
                                     `Minimum Detectable Effect Size` = 'mdes',
                                     `Constrained Optimal Sample Allocation` = "optimal"),
                         selected = 'power')
      ),
      column(width=4,
             helpText(h2("Estimand"), style="color:white"),
             uiOutput("effect"),
             conditionalPanel("input.effect == 'mod'",
                              uiOutput("modlev")
                              ),
             conditionalPanel("input.effect == 'mod' & input.totallev > input.modlev",
                              list(
                                selectInput('modtype', 'Moderator Type',
                                          choices =  c(`Non-randomly Varying` = 'n', `Randomly Varying` = 'r'),
                                          selected = 'r'))
                              )

      ),
      column(width=4,
             helpText(h2("Design"), style="color:white"),
             selectInput('totallev', 'Total Levels',
                         choices = c(`Single Level` = 1, `Two-Levels` = 2, `Three-Levels` = 3, `Four-Levels` = 4),
                         selected = 2),
             uiOutput("randlev"),
             uiOutput("block")

      )
    ),
    br(),
    htmlOutput('help')
  )
})
