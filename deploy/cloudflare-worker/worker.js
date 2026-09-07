// ArthOMix reverse proxy: arthomix.<account>.workers.dev -> Cloudflare tunnel -> localhost:7788
// Forwards every request (including Shiny's WebSocket upgrade) to ORIGIN_HOST unchanged.
export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    url.protocol = "https:";
    url.hostname = env.ORIGIN_HOST;
    url.port = "";

    const headers = new Headers(request.headers);
    headers.set("Host", env.ORIGIN_HOST);
    headers.delete("cf-connecting-ip");

    const init = {
      method: request.method,
      headers,
      body: request.method === "GET" || request.method === "HEAD" ? undefined : request.body,
      redirect: "manual",
    };
    // Shiny's session WebSocket: pass the upgrade through untouched.
    if (request.headers.get("Upgrade")?.toLowerCase() === "websocket") {
      return fetch(url.toString(), init);
    }
    const resp = await fetch(url.toString(), init);
    // Rewrite any absolute redirect back to the public hostname.
    const loc = resp.headers.get("Location");
    if (loc && loc.includes(env.ORIGIN_HOST)) {
      const h = new Headers(resp.headers);
      h.set("Location", loc.replace(env.ORIGIN_HOST, new URL(request.url).hostname));
      return new Response(resp.body, { status: resp.status, headers: h });
    }
    return resp;
  },
};
