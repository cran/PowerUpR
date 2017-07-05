shinyServer(function(input, output, session) {

  output$effect <- renderUI({
    req(input$totallev)
    req(input$randlev)
    if(input$totallev == input$randlev & input$totallev %in% c(2,3)){
      selectInput('effect', 'Request',
                  choices = c(`Main Effect` = 'main',
                              `Moderator Effect` = 'mod'),
                  selected = 'main')
    }else{
      selectInput('effect', 'Request',
                  choices = c(`Main Effect` = 'main'),
                  selected = 'main')
    }
  })

  output$block<- renderUI({
    req(input$totallev)
    req(input$randlev)
    if(input$totallev == 2 & input$randlev == 1){
      selectInput('block', 'Block Effect',
                  choices = c(`Random` = 'r', `Fixed` = 'f', `Constant` = 'c'),
                  selected = 'r')
    }else if(input$totallev == input$randlev){
      selectInput('block', 'Block Effect',
                  choices = c(`Random` = 'r'),
                  selected = 'r')
    }else if((input$totallev == 3 & input$randlev == 1) | (input$totallev == 4 & input$randlev%in% c(1,2))){
      selectInput('block', 'Block Effect',
                  choices = c(`Random` = 'r'),
                  selected = 'r')
    }else{
      selectInput('block', 'Block Effect',
                  choices = c(`Random` = 'r', `Fixed` = 'f'),
                  selected = 'r')
    }
  })

  output$randlev <- renderUI({
    req(input$totallev)
    if(input$totallev == 4){
      selectInput('randlev', 'Randomization Level',
                  choices = c(`At Level-1` = 1, `At Level-2` = 2, `At Level-3` = 3, `At Level-4` = 4),
                  selected = 4)
    }else if(input$totallev == 3){
      selectInput('randlev', 'Randomization Level',
                  choices = c(`At Level-1` = 1, `At Level-2` = 2, `At Level-3` = 3),
                  selected = 3)
    }else if(input$totallev == 2){
      selectInput('randlev', 'Randomization Level',
                  choices = c(`At Level-1` = 1, `At Level-2` = 2),
                  selected = 2)
    }else if(input$totallev == 1){
      selectInput('randlev', 'Randomization Level',
                  choices = c(`At Level-1` = 1),
                  selected = 1)
    }
  })

  output$modlev <- renderUI({
    req(input$totallev)
    if(input$totallev == 3){
      selectInput('modlev', 'Moderator Level',
                  choices =  c(`At Level-1` = 1, `At Level-2` = 2, `At Level-3` = 3),
                  selected = 1)
    }else if(input$totallev == 2){
      selectInput('modlev', 'Moderator Level',
                  choices =  c(`At Level-1` = 1, `At Level-2` = 2),
                  selected = 1)
    }
  })

  output.design <- reactive({
    req(input$totallev)
    req(input$randlev)
    req(input$modlev)
    req(input$modtype)
    req(input$output)
    req(input$block)
    if(input$effect == 'mod'){
      if(input$totallev > input$modlev){
        effect <- paste(input$effect, input$modlev, input$modtype, sep="")
      }else{
        effect <- paste(input$effect, input$modlev, sep="")
      }
      paste0(input$output, ".", effect, ".cra", input$totallev, input$block, input$randlev)
    }else if(input$effect == 'main'){
      if(input$randlev == input$totallev & input$totallev != 1){
        paste0(input$output, ".cra", input$totallev, input$block, input$randlev)
      }else if(input$randlev < input$totallev & input$randlev != 1){
        paste0(input$output, ".bcra", input$totallev, input$block, input$randlev)
      }else if(input$randlev < input$totallev & input$randlev == 1){
        paste0(input$output, ".bira", input$totallev, input$block, input$randlev)
      }else if(input$randlev == input$totallev & input$randlev == 1){
        paste0(input$output, ".ira", input$totallev, input$block, input$randlev)
      }
    }

  })

  # output$design <- renderPrint({
  #   req(output.design())
  #   fun.parsed <- scan(text = output.design(), what = "character", sep=".", quiet = TRUE)
  #   design <- fun.parsed[length(fun.parsed)]
  #   output.design()
  # })

  output$help <- renderPrint({
    fun <- output.design()
    help <- utils:::.getHelpFile(help(fun))
    tools:::Rd2HTML(help)
  })

})
