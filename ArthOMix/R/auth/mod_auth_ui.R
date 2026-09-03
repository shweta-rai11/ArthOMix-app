## R/auth/mod_auth_ui.R
## Authentication screen: sign up / log in / forgot password / reset
## password / "check your email", as sibling <div>s inside one tagList,

auth_password_field <- function(ns, input_id, label, placeholder) {
  fid <- ns(input_id)
  tagList(
    tags$label(class = "control-label", `for` = fid, label),
    tags$div(
      class = "auth-password-wrap",
      tags$input(
        id = fid, type = "password", class = "form-control",
        placeholder = placeholder, autocomplete = "off"
      ),
      tags$button(
        type = "button", class = "auth-password-toggle",
        title = "Show/hide password",
        onclick = sprintf(
          "var f=document.getElementById('%s'); f.type = f.type==='password' ? 'text' : 'password';
           this.querySelector('i').classList.toggle('fa-eye'); this.querySelector('i').classList.toggle('fa-eye-slash');",
          fid
        ),
        icon("eye")
      )
    )
  )
}

AUTH_PASSWORD_STRENGTH_JS <- "
function arthomixPasswordStrength(pw) {
  var score = 0;
  if (pw.length >= 8) score++;
  if (pw.length >= 12) score++;
  if (/[a-z]/.test(pw) && /[A-Z]/.test(pw)) score++;
  if (/[0-9]/.test(pw)) score++;
  if (/[^A-Za-z0-9]/.test(pw)) score++;
  return Math.min(score, 4);
}
function arthomixWirePasswordStrength(inputId, barId, labelId) {
  var input = document.getElementById(inputId);
  if (!input) return;
  input.addEventListener('input', function() {
    var score = arthomixPasswordStrength(input.value);
    var bar = document.getElementById(barId);
    var label = document.getElementById(labelId);
    var labels = ['Too short', 'Weak', 'Okay', 'Good', 'Strong'];
    var pct = input.value.length === 0 ? 0 : (20 + score * 20);
    if (bar) {
      bar.style.width = pct + '%';
      bar.className = 'auth-strength-fill auth-strength-' + score;
    }
    if (label) label.textContent = input.value.length === 0 ? '' : labels[score];
  });
}
"

auth_password_strength_ui <- function(ns, input_id) {
  bar_id <- paste0(ns(input_id), "_strength_bar")
  label_id <- paste0(ns(input_id), "_strength_label")
  tagList(
    tags$div(class = "auth-strength-track", tags$div(id = bar_id, class = "auth-strength-fill")),
    tags$div(id = label_id, class = "auth-strength-label"),
    tags$script(HTML(sprintf(
      "%s\n$(function(){ arthomixWirePasswordStrength('%s', '%s', '%s'); });",
      AUTH_PASSWORD_STRENGTH_JS, ns(input_id), bar_id, label_id
    )))
  )
}

mod_auth_ui <- function(id) {
  ns <- NS(id)

  tags$div(
    class = "auth-shell",
    tags$div(
      class = "auth-card",
      tags$div(class = "auth-brand", "ArthOMix"),
      tags$div(class = "auth-avatar", icon("circle-user", class = "auth-avatar-icon")),

      tags$div(
        id = ns("screen_login"),
        h3("Log in"),
        uiOutput(ns("login_status")),
        tags$label(class = "control-label", `for` = ns("login_email"), "Email address"),
        textInput(ns("login_email"), NULL, placeholder = "Enter your email address"),
        auth_password_field(ns, "login_password", "Password", "Enter your password"),
        actionButton(ns("login_btn"), "Log In", class = "btn btn-primary auth-primary-btn"),
        tags$div(
          class = "auth-links",
          actionLink(ns("go_forgot"), "Forgot password?"),
          actionLink(ns("go_signup"), "Don't have an account? Sign up")
        )
      ),

      tags$div(
        id = ns("screen_signup"), style = "display:none;",
        h3("Create your account"),
        uiOutput(ns("signup_status")),
        tags$label(class = "control-label", `for` = ns("signup_email"), "Email address"),
        textInput(ns("signup_email"), NULL, placeholder = "Enter your email address"),
        auth_password_field(ns, "signup_password", "Password", "Create a password"),
        auth_password_strength_ui(ns, "signup_password"),
        auth_password_field(ns, "signup_password_confirm", "Confirm password", "Confirm your password"),
        actionButton(ns("signup_btn"), "Create Account", class = "btn btn-primary auth-primary-btn"),
        tags$div(class = "auth-links", actionLink(ns("go_login_from_signup"), "Already have an account? Log In"))
      ),

      tags$div(
        id = ns("screen_check_email"), style = "display:none;",
        h3("Check your email"),
        uiOutput(ns("check_email_body")),
        actionButton(ns("resend_btn"), "Resend email", class = "btn btn-default"),
        tags$div(class = "auth-links", actionLink(ns("go_login_from_check"), "Return to Login"))
      ),

      tags$div(
        id = ns("screen_forgot"), style = "display:none;",
        h3("Reset your password"),
        p(class = "auth-subtext", "Enter your email address and we'll send you a link to reset your password."),
        uiOutput(ns("forgot_status")),
        tags$label(class = "control-label", `for` = ns("forgot_email"), "Email address"),
        textInput(ns("forgot_email"), NULL, placeholder = "Enter your email address"),
        actionButton(ns("forgot_btn"), "Send reset link", class = "btn btn-primary auth-primary-btn"),
        tags$div(class = "auth-links", actionLink(ns("go_login_from_forgot"), "Back to Login"))
      ),

      tags$div(
        id = ns("screen_reset_password"), style = "display:none;",
        h3("Choose a new password"),
        uiOutput(ns("reset_status")),
        auth_password_field(ns, "reset_password", "New password", "Create a password"),
        auth_password_strength_ui(ns, "reset_password"),
        auth_password_field(ns, "reset_password_confirm", "Confirm new password", "Confirm your password"),
        actionButton(ns("reset_btn"), "Update password", class = "btn btn-primary auth-primary-btn")
      )
    )
  )
}
