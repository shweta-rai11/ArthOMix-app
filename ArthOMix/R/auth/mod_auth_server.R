## R/auth/mod_auth_server.R
## Server logic for the authentication screens (R/auth/mod_auth_ui.R):
## validation, the Supabase REST calls (R/auth/auth_api.R), and the

is_valid_email <- function(x) {
  grepl("^[^@[:space:]]+@[^@[:space:]]+\\.[^@[:space:]]+$", x %||% "")
}

is_valid_password <- function(x) nchar(x %||% "") >= 8

auth_status_ui <- function(status) {
  if (is.null(status)) return(NULL)
  div(class = paste0("auth-alert auth-alert-", status$type), status$msg)
}

mod_auth_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    session_info <- reactiveVal(NULL)
    recovery_token <- reactiveVal(NULL)
    pending_email <- reactiveVal(NULL)
    check_email_mode <- reactiveVal("signup")

    login_status_msg <- reactiveVal(NULL)
    signup_status_msg <- reactiveVal(NULL)
    forgot_status_msg <- reactiveVal(NULL)
    reset_status_msg <- reactiveVal(NULL)

    show_screen <- function(target) {
      for (s in c("login", "signup", "check_email", "forgot", "reset_password")) {
        if (identical(s, target)) shinyjs::show(paste0("screen_", s)) else shinyjs::hide(paste0("screen_", s))
      }
    }

    observeEvent(input$login_btn, {
      login_status_msg(NULL)
      email <- trimws(input$login_email %||% "")
      pw <- input$login_password %||% ""
      if (!nzchar(email) || !nzchar(pw)) {
        login_status_msg(list(msg = "Enter your email and password.", type = "error"))
        return()
      }
      shinyjs::disable("login_btn")
      on.exit(shinyjs::enable("login_btn"), add = TRUE)
      res <- supabase_sign_in(email, pw)
      if (res$ok) {
        session_info(list(access_token = res$data$access_token, user = list(email = res$data$user$email %||% email)))
      } else {
        login_status_msg(list(msg = res$error, type = "error"))
      }
    }, ignoreInit = TRUE)

    observeEvent(input$signup_btn, {
      signup_status_msg(NULL)
      email <- trimws(input$signup_email %||% "")
      pw <- input$signup_password %||% ""
      pw2 <- input$signup_password_confirm %||% ""
      if (!is_valid_email(email)) {
        signup_status_msg(list(msg = "Enter a valid email address.", type = "error"))
        return()
      }
      if (!is_valid_password(pw)) {
        signup_status_msg(list(msg = "Password must be at least 8 characters.", type = "error"))
        return()
      }
      if (!identical(pw, pw2)) {
        signup_status_msg(list(msg = "Passwords do not match.", type = "error"))
        return()
      }
      shinyjs::disable("signup_btn")
      on.exit(shinyjs::enable("signup_btn"), add = TRUE)
      supabase_sign_up(email, pw)
      pending_email(email)
      check_email_mode("signup")
      show_screen("check_email")
    }, ignoreInit = TRUE)

    output$check_email_body <- renderUI({
      email <- pending_email() %||% "your email address"
      if (identical(check_email_mode(), "forgot")) {
        p(sprintf("If an account exists for %s, we've sent a link to reset your password.", email))
      } else {
        tagList(
          p(sprintf("We've sent a confirmation link to %s.", email)),
          p(class = "auth-subtext", "Click the link to verify your account, then log in.")
        )
      }
    })

    observeEvent(input$resend_btn, {
      email <- pending_email()
      req(email)
      shinyjs::disable("resend_btn")
      on.exit(shinyjs::enable("resend_btn"), add = TRUE)
      if (identical(check_email_mode(), "forgot")) supabase_recover(email) else supabase_resend(email, "signup")
      showNotification("Email sent.", type = "message", duration = 3)
    }, ignoreInit = TRUE)

    observeEvent(input$forgot_btn, {
      forgot_status_msg(NULL)
      email <- trimws(input$forgot_email %||% "")
      if (!is_valid_email(email)) {
        forgot_status_msg(list(msg = "Enter a valid email address.", type = "error"))
        return()
      }
      shinyjs::disable("forgot_btn")
      on.exit(shinyjs::enable("forgot_btn"), add = TRUE)
      supabase_recover(email)
      pending_email(email)
      check_email_mode("forgot")
      show_screen("check_email")
    }, ignoreInit = TRUE)

    observeEvent(input$reset_btn, {
      reset_status_msg(NULL)
      token <- recovery_token()
      pw <- input$reset_password %||% ""
      pw2 <- input$reset_password_confirm %||% ""
      if (is.null(token)) {
        reset_status_msg(list(msg = "This reset link is no longer valid. Please request a new one.", type = "error"))
        return()
      }
      if (!is_valid_password(pw)) {
        reset_status_msg(list(msg = "Password must be at least 8 characters.", type = "error"))
        return()
      }
      if (!identical(pw, pw2)) {
        reset_status_msg(list(msg = "Passwords do not match.", type = "error"))
        return()
      }
      shinyjs::disable("reset_btn")
      on.exit(shinyjs::enable("reset_btn"), add = TRUE)
      res <- supabase_update_password(token, pw)
      if (res$ok) {
        recovery_token(NULL)
        login_status_msg(list(msg = "Password updated. Please log in with your new password.", type = "success"))
        show_screen("login")
      } else {
        reset_status_msg(list(msg = res$error %||% "Couldn't update your password. Please request a new reset link.", type = "error"))
      }
    }, ignoreInit = TRUE)

    observeEvent(input$go_signup, show_screen("signup"), ignoreInit = TRUE)
    observeEvent(input$go_login_from_signup, show_screen("login"), ignoreInit = TRUE)
    observeEvent(input$go_forgot, show_screen("forgot"), ignoreInit = TRUE)
    observeEvent(input$go_login_from_forgot, show_screen("login"), ignoreInit = TRUE)
    observeEvent(input$go_login_from_check, show_screen("login"), ignoreInit = TRUE)

    observeEvent(session$clientData$url_search, {
      qs <- shiny::parseQueryString(session$clientData$url_search)
      token_hash <- qs$token_hash
      type <- qs$type
      if (is.null(token_hash) || is.null(type)) return()

      if (identical(type, "signup")) {
        res <- supabase_verify_otp(token_hash, "signup")
        if (res$ok) {
          session_info(list(access_token = res$data$access_token, user = list(email = res$data$user$email)))
        } else {
          login_status_msg(list(msg = "This confirmation link is invalid or has expired. Log in, or sign up again to get a new one.", type = "error"))
        }
      } else if (identical(type, "recovery")) {
        res <- supabase_verify_otp(token_hash, "recovery")
        if (res$ok) {
          recovery_token(res$data$access_token)
          show_screen("reset_password")
        } else {
          login_status_msg(list(msg = "This password reset link is invalid or has expired. Please request a new one.", type = "error"))
        }
      }
    }, once = TRUE)

    output$login_status <- renderUI(auth_status_ui(login_status_msg()))
    output$signup_status <- renderUI(auth_status_ui(signup_status_msg()))
    output$forgot_status <- renderUI(auth_status_ui(forgot_status_msg()))
    output$reset_status <- renderUI(auth_status_ui(reset_status_msg()))

    logout <- function() {
      info <- isolate(session_info())
      if (!is.null(info)) supabase_sign_out(info$access_token)
      session$reload()
    }

    list(session_info = session_info, logout = logout)
  })
}
