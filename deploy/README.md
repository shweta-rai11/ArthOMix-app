# Remote access

ArthOMix has no application-level login (see the "No authentication" note in
the main [README](../README.md)). Earlier experiments exposed the app through
a public Cloudflare quick-tunnel/Worker and a Caddy reverse proxy on a VPS —
both forwarded every request from the open internet straight to the Shiny
server with no restriction of any kind. Neither is used any more; they were
removed in favor of the approach below.

## SSH port-forward (only supported remote-access method)

Run the app as usual on the host machine:

```sh
cd ArthOMix
shiny::runApp()   # binds to 127.0.0.1:7788 per .Rprofile
```

From a client machine, forward a local port to it over SSH instead of
exposing any port publicly:

```sh
ssh -N -L 7788:localhost:7788 <user>@<host>
```

Then open `http://localhost:7788` in a browser on the client machine. Traffic
never leaves the SSH tunnel, and only accounts with SSH access to the host
(i.e. accounts you've added an authorized key for) can reach the app. There is
no separate credential to configure or rotate.

This does not provide a shareable public demo link — that would require
reintroducing a public listener, which is exactly what this setup avoids.
