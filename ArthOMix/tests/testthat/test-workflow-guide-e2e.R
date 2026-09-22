## End-to-end: the "Start an analysis" guide drives the existing Transcriptomics pages in a real browser.
## It must use the pages' own inputs and run button, so the stored DE result matches a manual run.

skip_if_not_installed("shinytest2")
skip_if_not_installed("chromote")

wf_bar_text <- function(app) gsub("\\s+", " ", tryCatch(app$get_text("#wf_bar"), error = function(e) ""))
wf_wait_bar <- function(app, pattern, timeout = 300) {
  deadline <- Sys.time() + timeout
  repeat {
    txt <- wf_bar_text(app)
    if (grepl(pattern, txt, fixed = TRUE) || Sys.time() >= deadline) return(txt)
    Sys.sleep(1)
  }
}

test_that("the guide walks Transcriptomics to a locked design using the existing pages", {
  app <- new_app_driver(name = "workflow-guide", height = 1000, width = 1400,
                        timeout = 120 * 1000, load_timeout = 600 * 1000)
  on.exit(app$stop(), add = TRUE)
  try(app$wait_for_idle(timeout = 60 * 1000), silent = TRUE)

  home <- app$get_html(".home-hero")
  expect_match(home, "Start an analysis", fixed = TRUE)
  expect_match(home, "Explore Modules", fixed = TRUE)
  expect_equal(lengths(regmatches(home, gregexpr("home-hero-flow-step", home))), 8L)

  app$click("home_start_analysis")
  Sys.sleep(2)
  ## wait_ = FALSE: picking a radio in the dialog updates no output, so the default wait would just time out.
  app$set_inputs(wf_layer = "transcriptomics", wf_data_choice = "example", wait_ = FALSE)
  app$click("wf_start_confirm")

  ## Step 2 - Check samples: the cohort composition, shown before any analysis design is picked.
  txt <- wf_wait_bar(app, "Step 2 of 8")
  expect_match(txt, "Step 2 of 8", fixed = TRUE)
  expect_match(txt, "Groups (group):", fixed = TRUE)

  ## Step 3 - Choose analysis: preset picker, the comparison (defaulted from the dataset's own
  ## group column, via effective_contrast()'s fallback - it must not wait on the Differential
  ## Expression page's own inputs to mount before showing something) and the sex design, together.
  app$click("wf_next")
  txt <- wf_wait_bar(app, "Step 3 of 8")
  expect_match(txt, "Step 3 of 8", fixed = TRUE)
  expect_match(txt, "Transcriptomics analysis design", fixed = TRUE)
  expect_equal(app$get_value(input = "tx_menu"), "Custom Analysis Design")
  txt <- wf_wait_bar(app, "Contrast column: ")
  expect_match(txt, "Contrast column: group", fixed = TRUE)
  expect_false(grepl("Waiting for the page's comparison controls", txt, fixed = TRUE))

  app$set_inputs(wf_sex_mode = "pooled")
  app$click("wf_next")
  Sys.sleep(3)
  ## Locking the design hands the rest of the pipeline to the existing sidebar-driven navigation
  ## (covered by the second test below) - the guide bar closes and the first chosen module opens.
  expect_equal(app$get_value(input = "tx_menu"), "Differential Expression")
  expect_equal(app$get_value(input = "tx_dge-contrast_col"), "group")
  expect_false(isTRUE(app$get_js("$('#wf_bar').is(':visible')")))
})

test_that("the Transcriptomics page is one step-driven navigation, and tabs follow the pipeline", {
  app <- new_app_driver(name = "workflow-shell", height = 1000, width = 1400,
                        timeout = 120 * 1000, load_timeout = 600 * 1000)
  on.exit(app$stop(), add = TRUE)
  try(app$wait_for_idle(timeout = 60 * 1000), silent = TRUE)
  app$set_inputs(sidebar_tabs = "transcriptomics", wait_ = FALSE)

  ## The sidebar list is a uiOutput, so wait for it rather than assuming it has rendered.
  for (i in 1:60) {
    if (isTRUE(app$get_js("$('.omics-sidebar[data-server-nav] .sidebar-step-header').length") >= 8)) break
    Sys.sleep(0.5)
  }
  expect_equal(app$get_js("$('.omics-sidebar[data-server-nav] .sidebar-step-header').length"), 8L)
  expect_equal(
    app$get_js("$('.omics-sidebar[data-server-nav] .sidebar-step-name').map(function(){return $(this).text()}).get().join('|')"),
    paste(c("Choose data", "Check samples", "Choose analysis", "Run analysis",
            "Explore biomarkers", "Biomarker modeling", "Validate", "Interpret biomarkers"), collapse = "|"))

  ## One navigation per level: the horizontal strip is hidden and the guide bar's own dot strip is gone.
  expect_true(app$get_js("$('#tx_menu').hasClass('nav-hidden') && $('#tx_menu').is(':hidden')"))
  expect_equal(app$get_js("$('.wf-bar-steps:visible').length"), 0L)

  ## Breadcrumb in place of the bare page title.
  expect_equal(app$get_js("$('.page-header-crumbed h2').length"), 0L)
  expect_match(app$get_js("$('.page-breadcrumb:visible').text().replace(/\\s+/g,' ').trim()"),
               "Home", fixed = TRUE)

  ## Added in reverse pipeline order, the tabs still come out in pipeline order.
  for (m in c("candidates", "wgcna", "dge")) { app$click(paste0("sm_toggle_", m)); Sys.sleep(2.5) }
  expect_equal(app$get_js("$('#tx_menu li a').map(function(){return $(this).data('value')}).get().join(' > ')"),
               paste("Dataset > Custom Analysis Design > Differential Expression",
                     "WGCNA Co-expression Network > Candidate Gene Identification > Sub-modules", sep = " > "))

  ## Candidates is locked until WGCNA has run, and says so - but its tab still opens.
  sidebar_row <- function(key) app$get_js(sprintf(
    "(function(){var a=$('.omics-sidebar[data-server-nav] a[onclick*=\"%s\\'\"]');return a.length? a.text().replace(/\\s+/g,' ').trim()+'|'+a.attr('class') : 'ABSENT';})()", key))
  expect_match(sidebar_row("candidates"), "sm-state-locked", fixed = TRUE)
  expect_match(app$get_js("$('#smreason_candidates').text()"), "WGCNA", fixed = TRUE)
  app$run_js("Shiny.setInputValue('tx_nav_go','candidates',{priority:'event'})"); Sys.sleep(2)
  expect_equal(app$get_value(input = "tx_menu"), "Candidate Gene Identification")

  ## Diagnostic Model is one tab with two sidebar rows; the step 7 row opens its External Validation tab.
  app$click("sm_toggle_diagnostic"); Sys.sleep(4)
  expect_match(sidebar_row("diagnostic"), "Diagnostic Model", fixed = TRUE)
  expect_match(sidebar_row("diagnostic_external"), "External Validation", fixed = TRUE)
  app$run_js("Shiny.setInputValue('tx_nav_go','diagnostic_external',{priority:'event'})"); Sys.sleep(5)
  expect_equal(app$get_value(input = "tx_menu"), "Diagnostic Model")
  expect_equal(app$get_value(input = "tx_diagnostic-main_tabs"), "External Validation")
  expect_match(sidebar_row("diagnostic_external"), "active", fixed = TRUE)
  expect_false(grepl("active", sidebar_row("diagnostic"), fixed = TRUE))

  ## Every other jump route still lands, and the breadcrumb follows.
  app$run_js("Shiny.setInputValue('header_search_submit','Immune Deconvolution',{priority:'event'})"); Sys.sleep(3)
  expect_equal(app$get_value(input = "tx_menu"), "Sub-modules")
  expect_match(app$get_js("$('.page-breadcrumb:visible').text()"), "Add analyses", fixed = TRUE)
  app$run_js("Shiny.setInputValue('tx_crumb_home', Math.random(), {priority:'event'})"); Sys.sleep(2)
  expect_equal(app$get_value(input = "sidebar_tabs"), "home")
})
