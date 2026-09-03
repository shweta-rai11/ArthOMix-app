## ui.R
## ArthOMix Explorer
## Layout follows xOmicsShiny's pattern: no marketing landing page, the app

home_stat_card <- function(value, label, icn, count_to = NULL, detail = NULL) {
  detail_id <- if (!is.null(detail)) paste0("home-stat-detail-", gsub("[^a-z0-9]+", "-", tolower(label))) else NULL
  div(
    class = "home-stat-card",
    tabindex = if (!is.null(detail)) "0" else NULL,
    `aria-describedby` = detail_id,
    div(class = "home-stat-icon", icon(icn)),
    div(
      div(
        class = "home-stat-value",
        if (is.null(count_to)) value else tags$span(class = "home-stat-counter", `data-count-to` = count_to, "0")
      ),
      div(class = "home-stat-label", label)
    ),
    if (!is.null(detail)) div(id = detail_id, class = "home-stat-detail", role = "tooltip", detail)
  )
}

home_feature_card <- function(icn, title, desc, kicker = NULL) {
  div(
    class = "info-card home-feature-card",
    div(class = "info-icon", icon(icn)),
    h4(title),
    if (!is.null(kicker)) div(class = "home-feature-kicker", kicker),
    p(desc)
  )
}

home_hero_title_words <- function(text, step = 0.09) {
  words <- strsplit(text, " ", fixed = TRUE)[[1]]
  tagList(lapply(seq_along(words), function(i) {
    tagList(
      tags$span(
        class = "home-hero-title-word",
        style = sprintf("animation-delay: %.2fs;", (i - 1) * step),
        words[i]
      ),
      if (i < length(words)) " "
    )
  }))
}

home_meta_chip <- function(icn, label, value) {
  div(
    class = "home-hero-meta-chip",
    icon(icn),
    div(
      div(class = "home-hero-meta-label", label),
      div(class = "home-hero-meta-value", value)
    )
  )
}

home_hero_flow_step <- function(n, icn, title, blurb) {
  div(
    class = "home-hero-flow-step",
    div(class = "home-hero-flow-num", n),
    div(
      class = "home-hero-flow-body",
      div(class = "home-hero-flow-title", icon(icn), title),
      p(class = "home-hero-flow-blurb", blurb)
    )
  )
}

home_hero_pill <- function(label, query = label) {
  tags$a(
    label, href = "#", class = "home-hero-pill",
    onclick = sprintf(
      "Shiny.setInputValue('header_search_submit', %s, {priority: 'event'}); return false;",
      jsonlite::toJSON(query, auto_unbox = TRUE)
    )
  )
}

home_faq_item <- function(question, answer) {
  tags$details(
    class = "home-faq-item",
    tags$summary(tags$span(question), tags$span(class = "home-faq-icon", icon("plus"))),
    div(class = "home-faq-answer", p(answer))
  )
}

HOME_HERO_CANVAS_JS <- "
$(function(){
  var canvas = document.getElementById('home_hero_canvas');
  if (!canvas || !canvas.getContext) return;
  var ctx = canvas.getContext('2d');
  var reduceMotion = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  var DPR = Math.min(window.devicePixelRatio || 1, 2);
  var W = 0, H = 0;
  var STAGES = [
    { end: 0.16, color: [100,116,139] },
    { end: 0.36, color: [37,99,235]  },
    { end: 0.60, color: [20,184,166] },
    { end: 0.85, color: [124,58,237] },
    { end: 1.01, color: [37,99,235]  }
  ];
  var N = 170;
  var particles = [];

  function rand(a, b) { return a + Math.random() * (b - a); }

  function spawn(seedX) {
    return {
      x: (seedX === undefined) ? rand(0, 1) : seedX,
      lane: rand(-1, 1),
      depth: rand(0.55, 1.35),
      speed: rand(0.00045, 0.0009),
      wobble: rand(0, Math.PI * 2),
      wobbleSpeed: rand(0.01, 0.025),
      wobbleAmp: rand(1, 3),
      survives1: Math.random() < 0.5,
      survives2: Math.random() < 0.35,
      dropped: false,
      alpha: (seedX === undefined) ? 1 : 0
    };
  }

  for (var i = 0; i < N; i++) particles.push(spawn());

  function stageColor(x) {
    for (var s = 0; s < STAGES.length; s++) if (x <= STAGES[s].end) return STAGES[s].color;
    return STAGES[STAGES.length - 1].color;
  }

  // funnel envelope: how far a lane can stray from the centerline at a
  // given x - near-full height on the left, narrowing toward the
  // centerline on the right, so surviving points visibly converge.
  function spreadAtX(x) {
    return 0.92 * Math.pow(Math.max(0, 1 - x), 1.4) + 0.06;
  }

  function laneY(p, x) {
    return 0.5 * H + p.lane * spreadAtX(x) * 0.47 * H + Math.sin(p.wobble) * p.wobbleAmp;
  }

  function resize() {
    var rect = canvas.getBoundingClientRect();
    W = rect.width; H = rect.height;
    canvas.width = Math.max(1, Math.round(W * DPR));
    canvas.height = Math.max(1, Math.round(H * DPR));
    ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
  }

  function update() {
    for (var i = 0; i < particles.length; i++) {
      var p = particles[i];
      p.wobble += p.wobbleSpeed;
      if (p.dropped) {
        p.alpha -= 0.02;
        if (p.alpha <= 0) particles[i] = spawn(0);
        continue;
      }
      p.x += p.speed;
      if (p.alpha < 1) p.alpha = Math.min(1, p.alpha + 0.02);
      if (p.x > 0.34 && p.x < 0.36 && !p.survives1) p.dropped = true;
      else if (p.x > 0.58 && p.x < 0.60 && p.survives1 && !p.survives2) p.dropped = true;
      else if (p.x > 1.02) particles[i] = spawn(0);
    }
  }

  function lerpColor(c1, c2, t) {
    return [
      Math.round(c1[0] + (c2[0] - c1[0]) * t),
      Math.round(c1[1] + (c2[1] - c1[1]) * t),
      Math.round(c1[2] + (c2[2] - c1[2]) * t)
    ];
  }

  // Soft radial wash, brightest over the convergence zone on the right -
  // pulls the eye toward 'biomarkers' the way a real analytics dashboard's
  // focal highlight would, instead of a flat uniform backdrop.
  function drawSpotlight() {
    var g = ctx.createRadialGradient(W * 0.92, H * 0.5, 0, W * 0.92, H * 0.5, W * 0.42);
    g.addColorStop(0, 'rgba(37,99,235,0.14)');
    g.addColorStop(1, 'rgba(37,99,235,0)');
    ctx.fillStyle = g;
    ctx.fillRect(0, 0, W, H);
  }

  function drawFunnelGuides() {
    ctx.strokeStyle = 'rgba(100,116,139,0.16)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, H * 0.04);
    ctx.quadraticCurveTo(W * 0.55, H * 0.5, W * 0.99, H * 0.485);
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(0, H * 0.96);
    ctx.quadraticCurveTo(W * 0.55, H * 0.5, W * 0.99, H * 0.515);
    ctx.stroke();
  }

  function draw() {
    ctx.clearRect(0, 0, W, H);
    ctx.shadowBlur = 0;
    drawSpotlight();
    drawFunnelGuides();
    ctx.lineWidth = 1;
    for (var i = 0; i < particles.length; i++) {
      var a = particles[i];
      if (a.dropped || a.x < 0.5 || a.x > 0.88) continue;
      for (var j = i + 1; j < particles.length; j++) {
        var b = particles[j];
        if (b.dropped || b.x < 0.5 || b.x > 0.88) continue;
        var ax = a.x * W, ay = laneY(a, a.x), bx = b.x * W, by = laneY(b, b.x);
        var dx = ax - bx, dy = ay - by;
        var dist = Math.sqrt(dx * dx + dy * dy);
        if (dist < 54) {
          var t = 1 - dist / 54;
          var ca = stageColor(a.x), cb = stageColor(b.x);
          var grad = ctx.createLinearGradient(ax, ay, bx, by);
          grad.addColorStop(0, 'rgba(' + ca.join(',') + ',' + (0.26 * t) + ')');
          grad.addColorStop(1, 'rgba(' + cb.join(',') + ',' + (0.26 * t) + ')');
          ctx.strokeStyle = grad;
          ctx.beginPath();
          ctx.moveTo(ax, ay);
          ctx.lineTo(bx, by);
          ctx.stroke();
        }
      }
    }
    for (var k = 0; k < particles.length; k++) {
      var p = particles[k];
      var x = p.x * W;
      var y = laneY(p, p.x);
      var c = stageColor(p.x);
      var isCandidate = p.x > 0.85 && p.survives1 && p.survives2 && !p.dropped;
      var r = (isCandidate ? 5.6 : (p.x > 0.36 ? 3.2 : 2.4)) * (0.75 + p.depth * 0.25);
      var alpha = p.alpha * (isCandidate ? 1 : 0.62) * Math.min(1, p.depth);
      if (isCandidate) {
        var pulse = 7 + Math.sin(p.wobble * 1.4) * 4;
        var glow = ctx.createRadialGradient(x, y, 0, x, y, r + pulse);
        glow.addColorStop(0, 'rgba(' + c.join(',') + ',' + (0.35 * p.alpha) + ')');
        glow.addColorStop(1, 'rgba(' + c.join(',') + ',0)');
        ctx.beginPath();
        ctx.arc(x, y, r + pulse, 0, Math.PI * 2);
        ctx.fillStyle = glow;
        ctx.fill();
        ctx.shadowColor = 'rgba(' + c.join(',') + ',0.75)';
        ctx.shadowBlur = 9;
      } else {
        ctx.shadowColor = 'rgba(' + c.join(',') + ',0.35)';
        ctx.shadowBlur = 2.5;
      }
      ctx.beginPath();
      ctx.arc(x, y, r, 0, Math.PI * 2);
      ctx.fillStyle = 'rgba(' + lerpColor(c, [255,255,255], isCandidate ? 0.35 : 0.1).join(',') + ',' + alpha + ')';
      ctx.fill();
    }
    ctx.shadowBlur = 0;
  }

  function drawStatic() {
    resize();
    for (var i = 0; i < particles.length; i++) {
      particles[i].x = i / particles.length;
      particles[i].alpha = 1;
      particles[i].dropped = false;
    }
    draw();
  }

  var running = false, rafId = null;
  function step() {
    if (!running) return;
    update();
    draw();
    rafId = window.requestAnimationFrame(step);
  }
  function start() {
    if (running || reduceMotion) return;
    running = true;
    step();
  }
  function stop() {
    running = false;
    if (rafId) window.cancelAnimationFrame(rafId);
  }

  resize();
  window.addEventListener('resize', resize);

  if (reduceMotion) {
    drawStatic();
  } else if ('IntersectionObserver' in window) {
    var obs = new IntersectionObserver(function(entries){
      entries.forEach(function(entry){
        if (entry.isIntersecting) start(); else stop();
      });
    }, {threshold: 0.05});
    obs.observe(canvas);
  } else {
    start();
  }
});
"

HOME_STORY_CANVAS_JS <- "
$(function(){
  var wrapper = document.querySelector('.home-page-root');
  var canvas = document.getElementById('home_story_canvas');
  if (!wrapper || !canvas || !canvas.getContext) return;
  var ctx = canvas.getContext('2d');
  var reduceMotion = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  var DPR = Math.min(window.devicePixelRatio || 1, 2);
  var W = 0, H = 0;

  function rand(a, b) { return a + Math.random() * (b - a); }

  function resize() {
    var rect = canvas.getBoundingClientRect();
    W = rect.width; H = rect.height;
    canvas.width = Math.max(1, Math.round(W * DPR));
    canvas.height = Math.max(1, Math.round(H * DPR));
    ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
  }

  var FONT = '10px -apple-system, BlinkMacSystemFont, sans-serif';

  /* ---- Stage 0: raw data (expression/methylation-style matrix) ---- */
  var matrixState = (function(){
    var cols = 26, rows = 13, cells = [];
    for (var r = 0; r < rows; r++) {
      for (var c = 0; c < cols; c++) {
        cells.push({ cx: (c + 0.5) / cols, cy: (r + 0.5) / rows, v: Math.random(), phase: rand(0, 6.28) });
      }
    }
    return { cols: cols, rows: rows, cells: cells };
  })();
  function drawMatrix(t, alpha) {
    if (alpha <= 0.01) return;
    var cw = W / matrixState.cols, ch = H / matrixState.rows;
    for (var i = 0; i < matrixState.cells.length; i++) {
      var cell = matrixState.cells[i];
      var shimmer = 0.15 * Math.sin(t * 0.6 + cell.phase);
      var v = Math.max(0.06, Math.min(0.95, cell.v + shimmer));
      ctx.fillStyle = 'rgba(37,99,235,' + (v * 0.5 * alpha) + ')';
      var w = cw * 0.72, h = ch * 0.72;
      ctx.fillRect(cell.cx * W - w / 2, cell.cy * H - h / 2, w, h);
    }
  }

  /* ---- Stage 1: explore (PCA-/clustering-style scatter) ---- */
  var explorePoints = (function(){
    var centers = [[0.32, 0.38], [0.62, 0.28], [0.5, 0.68]];
    var colors = [[37,99,235], [20,184,166], [124,58,237]];
    var pts = [];
    for (var g = 0; g < centers.length; g++) {
      for (var i = 0; i < 22; i++) {
        pts.push({
          bx: centers[g][0] + rand(-0.13, 0.13), by: centers[g][1] + rand(-0.13, 0.13),
          color: colors[g], phase: rand(0, 6.28), amp: rand(2, 6)
        });
      }
    }
    return pts;
  })();
  function drawExplore(t, alpha) {
    if (alpha <= 0.01) return;
    ctx.strokeStyle = 'rgba(100,116,139,' + (0.18 * alpha) + ')';
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(W * 0.08, H * 0.5); ctx.lineTo(W * 0.92, H * 0.5); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(W * 0.5, H * 0.1); ctx.lineTo(W * 0.5, H * 0.9); ctx.stroke();
    for (var i = 0; i < explorePoints.length; i++) {
      var p = explorePoints[i];
      var x = p.bx * W + Math.sin(t * 0.4 + p.phase) * p.amp;
      var y = p.by * H + Math.cos(t * 0.35 + p.phase) * p.amp;
      ctx.beginPath();
      ctx.arc(x, y, 2.4, 0, Math.PI * 2);
      ctx.fillStyle = 'rgba(' + p.color.join(',') + ',' + (0.55 * alpha) + ')';
      ctx.fill();
    }
  }

  /* ---- Stage 2: feature discovery (volcano-style ranking) ---- */
  var volcanoPoints = (function(){
    var pts = [];
    for (var i = 0; i < 90; i++) {
      var x = rand(-1, 1);
      var y = 0.15 + 0.55 * (x * x) + rand(-0.08, 0.08);
      pts.push({ x: x, y: Math.max(0.02, Math.min(0.98, y)), phase: rand(0, 6.28) });
    }
    pts.sort(function(a, b){ return b.y - a.y; });
    var labels = ['Feature 01', 'Feature 02', 'Feature 03'];
    for (var k = 0; k < labels.length; k++) pts[k].label = labels[k];
    return pts;
  })();
  function drawFeatureDiscovery(t, alpha) {
    if (alpha <= 0.01) return;
    ctx.font = FONT;
    for (var i = 0; i < volcanoPoints.length; i++) {
      var p = volcanoPoints[i];
      var x = W * (0.5 + p.x * 0.38);
      var y = H * (0.92 - p.y * 0.72) + Math.sin(t * 0.5 + p.phase) * 1.5;
      var hi = !!p.label;
      var r = hi ? 3.4 : 1.6;
      ctx.beginPath();
      ctx.arc(x, y, r, 0, Math.PI * 2);
      ctx.fillStyle = hi ? 'rgba(37,99,235,' + (0.8 * alpha) + ')' : 'rgba(100,116,139,' + (0.32 * alpha) + ')';
      ctx.fill();
      if (hi) {
        ctx.fillStyle = 'rgba(23,32,51,' + (0.55 * alpha) + ')';
        ctx.fillText(p.label, x + 6, y - 6);
      }
    }
  }

  /* ---- Stage 3: biological interpretation (pathway/enrichment network) ---- */
  var pathwayState = (function(){
    var n = 12, nodes = [];
    for (var i = 0; i < n; i++) {
      var ang = (i / n) * Math.PI * 2;
      nodes.push({ bx: 0.5 + Math.cos(ang) * rand(0.22, 0.34), by: 0.5 + Math.sin(ang) * rand(0.22, 0.34), phase: rand(0, 6.28) });
    }
    var edges = [];
    for (var e = 0; e < n; e++) edges.push([e, (e + 1) % n]);
    for (var e2 = 0; e2 < 5; e2++) edges.push([Math.floor(rand(0, n)), Math.floor(rand(0, n))]);
    var labelIdx = [0, 4, 8];
    var labels = ['Pathway A', 'Pathway B', 'Process C'];
    for (var li = 0; li < labelIdx.length; li++) nodes[labelIdx[li]].label = labels[li];
    return { nodes: nodes, edges: edges };
  })();
  function drawPathwayNetwork(t, alpha) {
    if (alpha <= 0.01) return;
    var nodes = pathwayState.nodes, edges = pathwayState.edges;
    ctx.lineWidth = 1;
    ctx.strokeStyle = 'rgba(20,184,166,' + (0.22 * alpha) + ')';
    ctx.beginPath();
    for (var e = 0; e < edges.length; e++) {
      var a = nodes[edges[e][0]], b = nodes[edges[e][1]];
      ctx.moveTo(a.bx * W, a.by * H);
      ctx.lineTo(b.bx * W, b.by * H);
    }
    ctx.stroke();
    ctx.font = FONT;
    for (var i = 0; i < nodes.length; i++) {
      var nd = nodes[i];
      var x = nd.bx * W + Math.sin(t * 0.3 + nd.phase) * 3;
      var y = nd.by * H + Math.cos(t * 0.3 + nd.phase) * 3;
      ctx.beginPath();
      ctx.arc(x, y, nd.label ? 4 : 2.2, 0, Math.PI * 2);
      ctx.fillStyle = 'rgba(20,184,166,' + ((nd.label ? 0.75 : 0.4) * alpha) + ')';
      ctx.fill();
      if (nd.label) {
        ctx.fillStyle = 'rgba(23,32,51,' + (0.55 * alpha) + ')';
        ctx.fillText(nd.label, x + 7, y - 7);
      }
    }
  }

  /* ---- Stage 4: multi-omics connected network (two layers + cross-links) ---- */
  var multiOmicsState = (function(){
    var top = [], bottom = [], side = [];
    for (var i = 0; i < 6; i++) top.push({ bx: 0.15 + i * 0.14, by: 0.28, phase: rand(0, 6.28) });
    for (var j = 0; j < 6; j++) bottom.push({ bx: 0.15 + j * 0.14, by: 0.68, phase: rand(0, 6.28) });
    for (var k = 0; k < 5; k++) {
      var ang = (k / 5) * Math.PI * 2;
      side.push({ bx: 0.86 + Math.cos(ang) * 0.08, by: 0.48 + Math.sin(ang) * 0.14, phase: rand(0, 6.28) });
    }
    var cross = [];
    for (var c = 0; c < 7; c++) cross.push([Math.floor(rand(0, 6)), Math.floor(rand(0, 6))]);
    return { top: top, bottom: bottom, side: side, cross: cross };
  })();
  function drawMultiOmics(t, alpha) {
    if (alpha <= 0.01) return;
    var L = multiOmicsState;
    function withinLayer(list, color) {
      ctx.strokeStyle = 'rgba(' + color.join(',') + ',' + (0.18 * alpha) + ')';
      ctx.beginPath();
      for (var i = 0; i < list.length - 1; i++) {
        ctx.moveTo(list[i].bx * W, list[i].by * H);
        ctx.lineTo(list[i + 1].bx * W, list[i + 1].by * H);
      }
      ctx.stroke();
    }
    withinLayer(L.top, [37,99,235]);
    withinLayer(L.bottom, [20,184,166]);
    ctx.strokeStyle = 'rgba(124,58,237,' + (0.16 * alpha) + ')';
    ctx.beginPath();
    for (var c = 0; c < L.cross.length; c++) {
      var a = L.top[L.cross[c][0]], b = L.bottom[L.cross[c][1]];
      ctx.moveTo(a.bx * W, a.by * H);
      ctx.lineTo(b.bx * W, b.by * H);
    }
    ctx.stroke();
    function dots(list, color, r) {
      for (var i = 0; i < list.length; i++) {
        var nd = list[i];
        var x = nd.bx * W + Math.sin(t * 0.3 + nd.phase) * 2.5;
        var y = nd.by * H + Math.cos(t * 0.3 + nd.phase) * 2.5;
        ctx.beginPath();
        ctx.arc(x, y, r, 0, Math.PI * 2);
        ctx.fillStyle = 'rgba(' + color.join(',') + ',' + (0.6 * alpha) + ')';
        ctx.fill();
      }
    }
    dots(L.top, [37,99,235], 2.6);
    dots(L.bottom, [20,184,166], 2.6);
    dots(L.side, [124,58,237], 2.2);
  }

  /* ---- Stage 5: biomarker prioritisation (network thinning to candidates) ---- */
  var biomarkerState = (function(){
    var n = 22, nodes = [];
    for (var i = 0; i < n; i++) nodes.push({ bx: rand(0.1, 0.9), by: rand(0.15, 0.85), phase: rand(0, 6.28), keep: i < 3 });
    var labels = ['Candidate biomarker', 'Feature importance', 'Prioritized signal'];
    for (var k = 0; k < 3; k++) nodes[k].label = labels[k];
    return nodes;
  })();
  function drawBiomarker(t, alpha) {
    if (alpha <= 0.01) return;
    ctx.font = FONT;
    for (var i = 0; i < biomarkerState.length; i++) {
      var nd = biomarkerState[i];
      var x = nd.bx * W + Math.sin(t * 0.3 + nd.phase) * 2.5;
      var y = nd.by * H + Math.cos(t * 0.3 + nd.phase) * 2.5;
      if (nd.keep) {
        var pulse = 5 + Math.sin(t * 1.3 + nd.phase) * 3;
        ctx.beginPath();
        ctx.arc(x, y, 4 + pulse, 0, Math.PI * 2);
        ctx.fillStyle = 'rgba(37,99,235,' + (0.08 * alpha) + ')';
        ctx.fill();
        ctx.beginPath();
        ctx.arc(x, y, 4, 0, Math.PI * 2);
        ctx.fillStyle = 'rgba(37,99,235,' + (0.85 * alpha) + ')';
        ctx.fill();
        ctx.fillStyle = 'rgba(23,32,51,' + (0.55 * alpha) + ')';
        ctx.fillText(nd.label, x + 8, y - 8);
      } else {
        ctx.beginPath();
        ctx.arc(x, y, 1.6, 0, Math.PI * 2);
        ctx.fillStyle = 'rgba(100,116,139,' + (0.22 * alpha) + ')';
        ctx.fill();
      }
    }
  }

  /* ---- Stage 6: modeling & validation (ROC-style curve + importance bars) ---- */
  var validationBars = (function(){
    var bars = [];
    for (var i = 0; i < 6; i++) bars.push({ h: rand(0.2, 0.85), phase: rand(0, 6.28) });
    return bars;
  })();
  function drawValidation(t, alpha) {
    if (alpha <= 0.01) return;
    ctx.strokeStyle = 'rgba(100,116,139,' + (0.18 * alpha) + ')';
    ctx.setLineDash([3, 4]);
    ctx.beginPath();
    ctx.moveTo(W * 0.12, H * 0.82);
    ctx.lineTo(W * 0.5, H * 0.2);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.strokeStyle = 'rgba(37,99,235,' + (0.5 * alpha) + ')';
    ctx.lineWidth = 1.6;
    ctx.beginPath();
    ctx.moveTo(W * 0.12, H * 0.82);
    ctx.bezierCurveTo(W * 0.2, H * 0.3, W * 0.34, H * 0.16, W * 0.5, H * 0.14);
    ctx.stroke();
    ctx.font = FONT;
    ctx.fillStyle = 'rgba(23,32,51,' + (0.5 * alpha) + ')';
    ctx.fillText('Model evaluation', W * 0.14, H * 0.12);
    var bx0 = W * 0.64, bw = 22, gap = 10;
    for (var i = 0; i < validationBars.length; i++) {
      var b = validationBars[i];
      var hh = (b.h + 0.05 * Math.sin(t * 0.6 + b.phase)) * H * 0.5;
      var x = bx0 + i * (bw + gap);
      ctx.fillStyle = 'rgba(124,58,237,' + (0.4 * alpha) + ')';
      ctx.fillRect(x, H * 0.82 - hh, bw, hh);
    }
    ctx.fillStyle = 'rgba(23,32,51,' + (0.5 * alpha) + ')';
    ctx.fillText('Cross-validation \\u00b7 Feature importance', bx0, H * 0.9);
  }

  var STAGE_DRAWERS = [drawMatrix, drawExplore, drawFeatureDiscovery, drawPathwayNetwork, drawMultiOmics, drawBiomarker, drawValidation];

  function scrollFraction() {
    var rect = wrapper.getBoundingClientRect();
    var vh = window.innerHeight || document.documentElement.clientHeight;
    var total = rect.height - vh;
    if (total <= 0) return 0;
    var scrolled = -rect.top;
    return Math.max(0, Math.min(1, scrolled / total));
  }

  var startTime = null;
  function drawFrame(now) {
    if (startTime === null) startTime = now;
    var t = (now - startTime) / 1000;
    var frac = scrollFraction();
    var idx = frac * (STAGE_DRAWERS.length - 1);
    var i0 = Math.floor(idx);
    var i1 = Math.min(i0 + 1, STAGE_DRAWERS.length - 1);
    var f = idx - i0;
    ctx.clearRect(0, 0, W, H);
    STAGE_DRAWERS[i0](t, 1 - f);
    if (i1 !== i0) STAGE_DRAWERS[i1](t, f);
  }

  function drawStatic() {
    resize();
    STAGE_DRAWERS[0](0, 1);
  }

  var running = false, rafId = null;
  function step(now) {
    if (!running) return;
    drawFrame(now);
    rafId = window.requestAnimationFrame(step);
  }
  function start() {
    if (running || reduceMotion) return;
    running = true;
    rafId = window.requestAnimationFrame(step);
  }
  function stop() {
    running = false;
    if (rafId) window.cancelAnimationFrame(rafId);
  }

  resize();
  window.addEventListener('resize', resize);

  if (reduceMotion) {
    drawStatic();
  } else if ('IntersectionObserver' in window) {
    var obs = new IntersectionObserver(function(entries){
      entries.forEach(function(entry){
        if (entry.isIntersecting) start(); else stop();
      });
    }, {threshold: 0});
    obs.observe(wrapper);
  } else {
    start();
  }
});
"

HOME_FEATURES_CANVAS_JS <- "
$(function(){
  var canvas = document.getElementById('home_features_canvas');
  if (!canvas || !canvas.getContext) return;
  var ctx = canvas.getContext('2d');
  var reduceMotion = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  var DPR = Math.min(window.devicePixelRatio || 1, 2);
  var W = 0, H = 0;
  var BLUE = [37,99,235];
  var BLUE_DARK = [29,78,216];

  function rand(a, b) { return a + Math.random() * (b - a); }
  function mix(a, t) { return 'rgba(' + a.join(',') + ',' + t + ')'; }

  function resize() {
    var rect = canvas.getBoundingClientRect();
    W = rect.width; H = rect.height;
    canvas.width = Math.max(1, Math.round(W * DPR));
    canvas.height = Math.max(1, Math.round(H * DPR));
    ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
  }

  // ---- 1: expression/methylation-style heatmap ----
  var heatmapCells = (function(){
    var cols = 22, rows = 11, cells = [];
    for (var r = 0; r < rows; r++) for (var c = 0; c < cols; c++)
      cells.push({ cx: (c + 0.5) / cols, cy: (r + 0.5) / rows, v: Math.random(), phase: rand(0, 6.28) });
    return { cols: cols, rows: rows, cells: cells };
  })();
  function drawHeatmap(t, a) {
    if (a <= 0.01) return;
    var cw = W / heatmapCells.cols, ch = H / heatmapCells.rows;
    for (var i = 0; i < heatmapCells.cells.length; i++) {
      var cell = heatmapCells.cells[i];
      var v = Math.max(0.05, Math.min(0.9, cell.v + 0.12 * Math.sin(t * 0.4 + cell.phase)));
      ctx.fillStyle = mix(BLUE, v * 0.55 * a);
      var w = cw * 0.7, h = ch * 0.7;
      ctx.fillRect(cell.cx * W - w / 2, cell.cy * H - h / 2, w, h);
    }
  }

  // ---- 2: PCA/UMAP-style projection ----
  var projPoints = (function(){
    var centers = [[0.3, 0.4], [0.62, 0.3], [0.5, 0.7]];
    var pts = [];
    for (var g = 0; g < centers.length; g++)
      for (var i = 0; i < 20; i++)
        pts.push({ bx: centers[g][0] + rand(-0.13, 0.13), by: centers[g][1] + rand(-0.13, 0.13), phase: rand(0, 6.28) });
    return pts;
  })();
  function drawProjection(t, a) {
    if (a <= 0.01) return;
    ctx.strokeStyle = mix([100,116,139], 0.16 * a);
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(W * 0.08, H * 0.5); ctx.lineTo(W * 0.92, H * 0.5); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(W * 0.5, H * 0.1); ctx.lineTo(W * 0.5, H * 0.9); ctx.stroke();
    for (var i = 0; i < projPoints.length; i++) {
      var p = projPoints[i];
      var x = p.bx * W + Math.sin(t * 0.35 + p.phase) * 4;
      var y = p.by * H + Math.cos(t * 0.3 + p.phase) * 4;
      ctx.beginPath();
      ctx.arc(x, y, 2.3, 0, Math.PI * 2);
      ctx.fillStyle = mix(BLUE, 0.4 * a);
      ctx.fill();
    }
  }

  // ---- 3: volcano-style ranking ----
  var volcanoPoints = (function(){
    var pts = [];
    for (var i = 0; i < 70; i++) {
      var x = rand(-1, 1);
      var y = 0.15 + 0.55 * (x * x) + rand(-0.08, 0.08);
      pts.push({ x: x, y: Math.max(0.02, Math.min(0.98, y)), phase: rand(0, 6.28), hi: Math.abs(x) > 0.62 && Math.random() < 0.5 });
    }
    return pts;
  })();
  function drawVolcano(t, a) {
    if (a <= 0.01) return;
    for (var i = 0; i < volcanoPoints.length; i++) {
      var p = volcanoPoints[i];
      var x = W * (0.5 + p.x * 0.4);
      var y = H * (0.9 - p.y * 0.72) + Math.sin(t * 0.5 + p.phase) * 1.5;
      ctx.beginPath();
      ctx.arc(x, y, p.hi ? 2.8 : 1.5, 0, Math.PI * 2);
      ctx.fillStyle = mix(p.hi ? BLUE_DARK : BLUE, (p.hi ? 0.55 : 0.28) * a);
      ctx.fill();
    }
  }

  // ---- 4: pathway/enrichment network ----
  var pathwayState = (function(){
    var n = 11, nodes = [];
    for (var i = 0; i < n; i++) {
      var ang = (i / n) * Math.PI * 2;
      nodes.push({ bx: 0.5 + Math.cos(ang) * rand(0.2, 0.32), by: 0.5 + Math.sin(ang) * rand(0.2, 0.32), phase: rand(0, 6.28) });
    }
    var edges = [];
    for (var e = 0; e < n; e++) edges.push([e, (e + 1) % n]);
    for (var e2 = 0; e2 < 4; e2++) edges.push([Math.floor(rand(0, n)), Math.floor(rand(0, n))]);
    return { nodes: nodes, edges: edges };
  })();
  function drawPathwayNetwork(t, a) {
    if (a <= 0.01) return;
    var nodes = pathwayState.nodes, edges = pathwayState.edges;
    ctx.strokeStyle = mix(BLUE, 0.2 * a);
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (var e = 0; e < edges.length; e++) {
      var na = nodes[edges[e][0]], nb = nodes[edges[e][1]];
      ctx.moveTo(na.bx * W, na.by * H);
      ctx.lineTo(nb.bx * W, nb.by * H);
    }
    ctx.stroke();
    for (var i = 0; i < nodes.length; i++) {
      var nd = nodes[i];
      var x = nd.bx * W + Math.sin(t * 0.3 + nd.phase) * 2.5;
      var y = nd.by * H + Math.cos(t * 0.3 + nd.phase) * 2.5;
      ctx.beginPath();
      ctx.arc(x, y, 2.4, 0, Math.PI * 2);
      ctx.fillStyle = mix(BLUE, 0.42 * a);
      ctx.fill();
    }
  }

  // ---- 5: multi-omics integration network (two layers + cross-links) ----
  var multiOmicsState = (function(){
    var top = [], bottom = [];
    for (var i = 0; i < 5; i++) top.push({ bx: 0.18 + i * 0.16, by: 0.3, phase: rand(0, 6.28) });
    for (var j = 0; j < 5; j++) bottom.push({ bx: 0.18 + j * 0.16, by: 0.7, phase: rand(0, 6.28) });
    var cross = [];
    for (var c = 0; c < 6; c++) cross.push([Math.floor(rand(0, 5)), Math.floor(rand(0, 5))]);
    return { top: top, bottom: bottom, cross: cross };
  })();
  function drawMultiOmicsNetwork(t, a) {
    if (a <= 0.01) return;
    var L = multiOmicsState;
    function within(list) {
      ctx.strokeStyle = mix(BLUE, 0.16 * a);
      ctx.beginPath();
      for (var i = 0; i < list.length - 1; i++) {
        ctx.moveTo(list[i].bx * W, list[i].by * H);
        ctx.lineTo(list[i + 1].bx * W, list[i + 1].by * H);
      }
      ctx.stroke();
    }
    within(L.top); within(L.bottom);
    ctx.strokeStyle = mix(BLUE_DARK, 0.14 * a);
    ctx.beginPath();
    for (var c = 0; c < L.cross.length; c++) {
      var na = L.top[L.cross[c][0]], nb = L.bottom[L.cross[c][1]];
      ctx.moveTo(na.bx * W, na.by * H);
      ctx.lineTo(nb.bx * W, nb.by * H);
    }
    ctx.stroke();
    function dots(list) {
      for (var i = 0; i < list.length; i++) {
        var nd = list[i];
        var x = nd.bx * W + Math.sin(t * 0.3 + nd.phase) * 2;
        var y = nd.by * H + Math.cos(t * 0.3 + nd.phase) * 2;
        ctx.beginPath();
        ctx.arc(x, y, 2.6, 0, Math.PI * 2);
        ctx.fillStyle = mix(BLUE, 0.45 * a);
        ctx.fill();
      }
    }
    dots(L.top); dots(L.bottom);
  }

  // ---- 6: predictive-model evaluation (ROC-style curve + bars) ----
  var modelBars = (function(){
    var bars = [];
    for (var i = 0; i < 5; i++) bars.push({ h: rand(0.2, 0.85), phase: rand(0, 6.28) });
    return bars;
  })();
  function drawModelOutput(t, a) {
    if (a <= 0.01) return;
    ctx.strokeStyle = mix([100,116,139], 0.16 * a);
    ctx.setLineDash([3, 4]);
    ctx.beginPath();
    ctx.moveTo(W * 0.1, H * 0.82); ctx.lineTo(W * 0.48, H * 0.18);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.strokeStyle = mix(BLUE, 0.45 * a);
    ctx.lineWidth = 1.6;
    ctx.beginPath();
    ctx.moveTo(W * 0.1, H * 0.82);
    ctx.bezierCurveTo(W * 0.18, H * 0.28, W * 0.32, H * 0.14, W * 0.48, H * 0.12);
    ctx.stroke();
    var bx0 = W * 0.62, bw = 18, gap = 9;
    for (var i = 0; i < modelBars.length; i++) {
      var b = modelBars[i];
      var hh = (b.h + 0.05 * Math.sin(t * 0.6 + b.phase)) * H * 0.5;
      ctx.fillStyle = mix(BLUE_DARK, 0.32 * a);
      ctx.fillRect(bx0 + i * (bw + gap), H * 0.82 - hh, bw, hh);
    }
  }

  var STAGE_DRAWERS = [drawHeatmap, drawProjection, drawVolcano, drawPathwayNetwork, drawMultiOmicsNetwork, drawModelOutput];
  var HOLD = 3.2, FADE = 1.8;
  var CYCLE = HOLD + FADE;
  var TOTAL = STAGE_DRAWERS.length * CYCLE;

  var startTime = null;
  function drawFrame(now) {
    if (startTime === null) startTime = now;
    var t = (now - startTime) / 1000;
    var tMod = t % TOTAL;
    var idx = Math.floor(tMod / CYCLE);
    var within = tMod - idx * CYCLE;
    var i0 = idx, i1 = (idx + 1) % STAGE_DRAWERS.length;
    var f = within <= HOLD ? 0 : Math.min(1, (within - HOLD) / FADE);
    ctx.clearRect(0, 0, W, H);
    STAGE_DRAWERS[i0](t, 1 - f);
    if (f > 0) STAGE_DRAWERS[i1](t, f);
  }

  function drawStatic() {
    resize();
    drawHeatmap(0, 1);
  }

  var running = false, rafId = null;
  function step(now) {
    if (!running) return;
    drawFrame(now);
    rafId = window.requestAnimationFrame(step);
  }
  function start() {
    if (running || reduceMotion) return;
    running = true;
    rafId = window.requestAnimationFrame(step);
  }
  function stop() {
    running = false;
    if (rafId) window.cancelAnimationFrame(rafId);
  }

  resize();
  window.addEventListener('resize', resize);

  if (reduceMotion) {
    drawStatic();
  } else if ('IntersectionObserver' in window) {
    var obs = new IntersectionObserver(function(entries){
      entries.forEach(function(entry){
        if (entry.isIntersecting) start(); else stop();
      });
    }, {threshold: 0});
    obs.observe(canvas);
  } else {
    start();
  }
});
"

homeUI <- function() {
  tagList(
    div(
      class = "home-page-root",
      tags$canvas(id = "home_story_canvas", class = "home-story-canvas", `aria-hidden` = "true"),
      div(
        class = "home-story-content",
    div(
      class = "home-hero",
      div(class = "home-hero-pattern"),
      div(class = "home-hero-orb home-hero-orb-1", `aria-hidden` = "true"),
      div(class = "home-hero-orb home-hero-orb-2", `aria-hidden` = "true"),
      tags$canvas(id = "home_hero_canvas", class = "home-hero-canvas", `aria-hidden` = "true"),
      div(
        class = "home-hero-grid",
      div(
        class = "home-hero-content",
        span(
          class = "home-hero-eyebrow",
          tags$span(class = "home-hero-live-dot"),
          "Omics biomarker discovery platform"
        ),
        h1(class = "home-hero-title", home_hero_title_words("From raw omics data to potential biomarkers.")),
        p(
          class = "home-hero-subtitle",
          "ArthOMix is an Shiny based web application used for transcriptomics, methylomics and integrated platform for identifying potential biomarkers across different omics layers."
        ),
        div(
          class = "home-hero-actions",
          actionButton("home_browse_modules", tagList(icon("table-cells-large"), "Explore Modules"), class = "btn btn-default"),
          tags$a(icon("comments"), " Ask ArthOChat", href = "#", class = "btn btn-default", onclick = ARTHOCHAT_DRAWER_OPEN_JS)
        ),
        div(class = "home-hero-search-label", "OR GO STRAIGHT TO A MODULE"),
        div(
          class = "home-hero-search",
          tags$span(class = "home-hero-search-icon", icon("magnifying-glass")),
          tags$input(
            type = "text", id = "home_hero_search", class = "form-control",
            placeholder = "Try “WGCNA”, “TNF”, “diagnostic model”…", autocomplete = "off"
          ),
          tags$span(class = "home-hero-search-kbd", "Enter")
        ),
        div(
          class = "home-hero-meta-row",
          home_meta_chip("layer-group", "Omics coverage", "Transcriptomics • Methylomics • Cross-Omics • Multi-omics"),
          home_meta_chip("sliders", "Analysis design", "Sex-stratified • Sex-specific • Sex-pooled"),
          home_meta_chip("file-csv", "Accepted data", "CSV, RDS file formats.")
        )
      ),
      div(
        class = "home-hero-flow",
        home_hero_flow_step("1", "layer-group", "Pick a module layer",
          "Transcriptomics, methylomics and integrated."),
        home_hero_flow_step("2", "route", "Run the analysis",
          "Analysis based on modules."),
        home_hero_flow_step("3", "flask-vial", "Identify potential biomarkers",
          "To validated potential biomarker for wet-lab follow-up.")
      )
      ),
      tags$script(HTML(
        "$(function(){
           $('#home_hero_search').on('keydown', function(e){
             if (e.which === 13) {
               Shiny.setInputValue('header_search_submit', $(this).val(), {priority: 'event'});
             }
           });
         });"
      )),
      tags$script(HTML(HOME_HERO_CANVAS_JS))
    ),
    div(
      class = "home-goal-band",
      div(class = "home-goal-icon", icon("bullseye")),
      div(
        class = "home-goal-body",
        div(class = "home-goal-eyebrow", "The goal"),
        h3(class = "home-goal-title", "Choose the analysis."),
        p(
          class = "home-goal-text",
          "ArthOMix is designed as a modular application. The user can either uplaod raw data or normalised data for each of the omcis layers. For each modules (Transcriptomics, methylomics, Integrated) it has multiple sub-modules. See the upload guide for file format details before you upload the data"
        ),
        tags$a(
          class = "home-goal-chat-link", href = "#", onclick = ARTHOCHAT_DRAWER_OPEN_JS,
          icon("comments"), " Stuck at any step? Ask ArthOChat."
        )
      )
    ),
    div(
      class = "home-features-section",
      tags$canvas(id = "home_features_canvas", class = "home-features-canvas", `aria-hidden` = "true"),
      div(
        class = "page-header",
        h2("Why ArthOMix"),
        p("One environment for reproducible biomarker discovery.")
      ),
      div(
        class = "home-features-grid",
        home_feature_card("route", "Multiple analysis",
          "From raw data to potential biomarkers, feature selection, diagnostic modeling, and biological interpretation within one analytical workflow."),
        home_feature_card("sliders", "Flexible, study-aware analysis",
          "Run sex-stratified, sex-specific, or sex-pooled analysis based on the research design."),
        home_feature_card("comments", "Analysis-grounded AI assistant",
          "Struggling with a basic concept? Learn while analysing and interpret biological findings with responses grounded in your analysis.",
          kicker = "ArthOChat, grounded in your data"),
        home_feature_card("flask", "Supports single and multiomics",
          "All in one app, from loading data to preprocessing, statistical analysis, validation, and interpreting results.")
      ),
      tags$script(HTML(HOME_FEATURES_CANVAS_JS))
    ),
    div(
      class = "page-header",
      h2("Frequently asked"),
      p("What people ask before analysis.")
    ),
    div(
      class = "home-faq-list",
      home_faq_item(
        "Do I need to write code to run any analysis?",
        "No. Every step - preprocessing, differential expression or differential methylation analysis, feature selection, model training and evaluation is click and anlyse inside the platform's sub-modules. No R or Python is required."
      ),
      home_faq_item(
        "What types of omics data can ArthOMix analyze?",
        "It can analyse Transcriptomics,methylomics, intergated study like cross-omics or multiomics."
      ),
      home_faq_item(
        "Can I load my own data?",
        "Yes - Use the uplaod data to upload your own data."
      ),
      home_faq_item(
        "Can I perform pooled and sex-specific analyses?",
        "Yes. It can perform sex specific analysis and all the results are downloadable in CSV or PNG file format."
      ),
      home_faq_item(
        "What data does ArthOMix provide?",
        "A platform to analyse omics data. It loads preloaded data as an example, so you can explore the full pipeline before bringing your own data."
      ),
      home_faq_item(
        "Why it is named as ArthOMix?",
        "This shiny based application allows the user to analyse any omics data but since the preloaded data (example) data is confined to Rheumatoid arthritis so the name \"Arth\", \"O\" stands for omics and \"Mix\" because it can perform both single and multiomics analysis based on the data provided."
      )
    ),
    div(
      class = "home-cta-banner",
      div(class = "home-cta-pattern"),
      div(
        class = "home-cta-content",
        h3("Ready to explore your data?"),
        p("Explore the platform's analysis modules to get started."),
        div(
          class = "home-cta-actions",
          actionButton("home_cta_browse_modules", tagList(icon("table-cells-large"), "Explore Modules"), class = "btn btn-primary")
        ),
        p(class = "home-cta-subtext", "An example dataset is already loaded, so you can explore right away - or swap in your own data from the Dataset tab.")
      )
    ),
    div(
      class = "home-footer",
      div(class = "home-footer-credit", "Developed by the Kovalchuk Lab, University of Lethbridge"),
      div(class = "home-footer-version", "ArthOMix v1.0.0")
    )
      ),
      tags$script(HTML(HOME_STORY_CANVAS_JS))
    )
  )
}

moduleCardUI <- function(m, id_prefix = "") {
  is_available <- identical(m$status, "available")
  div(
    class = paste("module-card", if (!is_available) "module-card-disabled"),
    div(class = "module-card-icon", icon(m$icon)),
    div(
      class = "module-card-body",
      div(
        class = "module-card-title-row",
        h4(m$title),
        span(
          class = if (is_available) "badge-available" else "badge-soon",
          if (is_available) "Available" else "Coming soon"
        )
      ),
      p(class = "module-card-tagline", m$tagline),
      if (identical(m$id, "arthochat")) {
        tags$a("Open module", href = "#", class = "btn btn-primary btn-sm", onclick = ARTHOCHAT_DRAWER_OPEN_JS)
      } else if (is_available) {
        actionButton(paste0(id_prefix, "open_", m$id), "Open module", class = "btn-primary btn-sm")
      } else {
        tags$button("Open module", class = "btn btn-sm btn-default", disabled = "disabled")
      }
    )
  )
}

modulesLandingUI <- function() {
  single <- Filter(function(m) m$kind == "Single-omics", MODULE_REGISTRY)
  multi  <- Filter(function(m) m$kind == "Multi-omics", MODULE_REGISTRY)
  assist <- Filter(function(m) m$kind == "Assistant", MODULE_REGISTRY)
  tagList(
    div(
      class = "page-header",
      h2("Modules"),
      p("Each module is one omics layer. Transcriptomics, Methylomics, and Integrated omics. Open one to load a dataset, or run preloadeed data or fetch data from NCBI GEO datasets, add sub-modules and perform analysis")
    ),
    h4(class = "module-group-title", "Single-omics modules"),
    div(class = "module-grid", lapply(single, moduleCardUI)),
    h4(class = "module-group-title", "Multi-omics modules"),
    div(class = "module-grid", lapply(multi, moduleCardUI)),
    if (length(assist) > 0) tagList(
      h4(class = "module-group-title", "Assistant"),
      div(class = "module-grid", lapply(assist, moduleCardUI))
    )
  )
}

SUBMODULE_GROUP_ORDER <- c(
  "Data", "Network", "Genetics", "Biomarker modeling", "Validation", "Interpretation"
)
SUBMODULE_GROUP_BLURB <- c(
  "Data" = "Load, inspect and prepare the expression matrix.",
  "Network" = "Co-expression structure and the candidate genes it points to.",
  "Genetics" = "Causal evidence from GWAS summary statistics.",
  "Biomarker modeling" = "Turn candidate genes into a panel and a diagnostic model.",
  "Validation" = "Check the panel holds up outside the discovery cohort.",
  "Interpretation" = "What the panel means biologically and clinically."
)

submoduleCardUI <- function(cfg, id_prefix = "") {
  div(
    id = paste0(id_prefix, "smcard_wrap_", cfg$id), class = "sm-card-wrap",
    `data-sm-title` = tolower(cfg$title),
    div(
      class = "sm-card", id = paste0(id_prefix, "smcard_", cfg$id),
      div(class = "sm-card-icon", icon(cfg$icon)),
      div(
        class = "sm-card-body",
        div(class = "sm-card-title-row", h4(cfg$title)),
        p(class = "sm-card-desc", cfg$description)
      ),
      tags$button(
        id = paste0(id_prefix, "sm_toggle_", cfg$id), type = "button",
        class = "btn sm-toggle-btn action-button",
        span(id = paste0(id_prefix, "smstate_", cfg$id), "Add")
      )
    )
  )
}

build_submodule_grid <- function(modules = TX_MODULES, group_order = SUBMODULE_GROUP_ORDER, group_blurb = SUBMODULE_GROUP_BLURB, id_prefix = "") {
  by_group <- split(modules, vapply(modules, function(m) m$config$group %||% "Data", character(1)))
  groups <- intersect(group_order, names(by_group))
  tagList(
    lapply(groups, function(g) {
      tagList(
        div(
          class = "sm-group-header",
          h4(g),
          p(class = "sm-group-blurb", group_blurb[[g]])
        ),
        div(class = "sm-grid", lapply(by_group[[g]], function(m) submoduleCardUI(m$config, id_prefix)))
      )
    })
  )
}

TRANSCRIPTOMICS_SIDEBAR_NAV <- list(
  list(id = "dataset", label = "Datasets", icon = "database", match = "Dataset"),
  list(id = "submodules", label = "Sub-modules", icon = "layer-group", match = "Sub-modules")
)

transcriptomicsUI <- function() {
  fluidRow(
    column(3, div(class = "omics-sidebar-col", omics_sidebar(
      "transcriptomics", "Transcriptomics", TRANSCRIPTOMICS_SIDEBAR_NAV,
      extra_sidebar_content = uiOutput("tx_sidebar_arthochat_hint"),
      dynamic_nav_output_id = "tx_sidebar_dynamic_nav"
    ))),
    column(
      9,
      div(
        class = "page-header page-header-tight",
        div(class = "page-header-pattern"),
        h2(icon("dna"), " Transcriptomics"),
        uiOutput("tx_page_subtitle")
      ),
      div(
        class = "tx-menu-wrap",
        tabsetPanel(
          id = "tx_menu", type = "tabs",
          tabPanel("Dataset", br(), mod_dataset_ui("tx_dataset")),
          tabPanel(
            "Sub-modules", br(),
            div(
              class = "sm-toolbar",
              div(class = "sm-toolbar-count", uiOutput("sm_active_count")),
              div(
                class = "sm-toolbar-search",
                tags$span(icon("magnifying-glass")),
                textInput("sm_search", NULL, placeholder = "Filter sub-modules by name...", width = "260px")
              )
            ),
            build_submodule_grid()
          )
        )
      )
    )
  )
}

MX_SUBMODULE_GROUP_ORDER <- c("Data", "Network", "Genetics", "Biomarker modeling", "Interpretation")
MX_SUBMODULE_GROUP_BLURB <- c(
  "Data" = "Load, quality-control, and prepare the methylation dataset.",
  "Network" = "Co-methylation structure and the candidate CpGs it points to.",
  "Genetics" = "Causal evidence from mQTL and GWAS summary statistics.",
  "Biomarker modeling" = "Turn candidate CpGs into a panel and a diagnostic model.",
  "Interpretation" = "Profile a single candidate biomarker in depth."
)

METHYLOMICS_SIDEBAR_NAV <- list(
  list(id = "dataset", label = "Dataset", icon = "database", match = "Dataset"),
  list(id = "submodules", label = "Sub-modules", icon = "layer-group", match = "Sub-modules")
)

methylomicsUI <- function() {
  fluidRow(
    column(3, div(class = "omics-sidebar-col", omics_sidebar(
      "methylomics", "Methylomics", METHYLOMICS_SIDEBAR_NAV,
      extra_sidebar_content = uiOutput("mx_sidebar_arthochat_hint"),
      dynamic_nav_output_id = "mx_sidebar_dynamic_nav"
    ))),
    column(
      9,
      div(
        class = "page-header page-header-tight",
        div(class = "page-header-pattern"),
        h2(icon("circle-nodes"), " Methylomics"),
        uiOutput("mx_page_subtitle")
      ),
      div(
        class = "tx-menu-wrap",
        tabsetPanel(
          id = "mx_menu", type = "tabs",
          tabPanel("Dataset", br(), mod_methyl_dataset_ui("mx_dataset")),
          tabPanel(
            "Sub-modules", br(),
            div(
              class = "sm-toolbar",
              div(class = "sm-toolbar-count", uiOutput("mx_sm_active_count")),
              div(
                class = "sm-toolbar-search",
                tags$span(icon("magnifying-glass")),
                textInput("mx_sm_search", NULL, placeholder = "Filter sub-modules by name...", width = "260px")
              )
            ),
            build_submodule_grid(MX_MODULES, MX_SUBMODULE_GROUP_ORDER, MX_SUBMODULE_GROUP_BLURB, id_prefix = "mx_")
          )
        )
      )
    )
  )
}

CX_SUBMODULE_GROUP_ORDER <- c("Data", "Genetics")
CX_SUBMODULE_GROUP_BLURB <- c(
  "Data" = "Join the eQTL-MR/mQTL-MR/DEG/DMP tables into one biomarker-convergence view.",
  "Genetics" = "Bridge the two causal panels with a follow-up cross-omics MR stage."
)

CROSSOMICS_SIDEBAR_NAV <- list(
  list(id = "dataset", label = "Dataset", icon = "database", match = "Dataset"),
  list(id = "submodules", label = "Sub-modules", icon = "layer-group", match = "Sub-modules")
)

crossomicsUI <- function() {
  fluidRow(
    column(3, div(class = "omics-sidebar-col", omics_sidebar(
      "crossomics", "Cross-Omics", CROSSOMICS_SIDEBAR_NAV,
      extra_sidebar_content = uiOutput("cx_arthochat_hint"),
      dynamic_nav_output_id = "cx_sidebar_dynamic_nav"
    ))),
    column(
      9,
      div(
        class = "page-header page-header-tight",
        div(class = "page-header-pattern"),
        h2(icon("arrows-left-right"), " Cross-Omics"),
        uiOutput("cx_page_subtitle")
      ),
      div(
        class = "tx-menu-wrap",
        tabsetPanel(
          id = "cx_menu", type = "tabs",
          tabPanel("Dataset", br(), mod_cross_dataset_ui("cx_dataset")),
          tabPanel(
            "Sub-modules", br(),
            div(
              class = "sm-toolbar",
              div(class = "sm-toolbar-count", uiOutput("cx_sm_active_count")),
              div(
                class = "sm-toolbar-search",
                tags$span(icon("magnifying-glass")),
                textInput("cx_sm_search", NULL, placeholder = "Filter sub-modules by name...", width = "260px")
              )
            ),
            build_submodule_grid(CX_MODULES, CX_SUBMODULE_GROUP_ORDER, CX_SUBMODULE_GROUP_BLURB, id_prefix = "cx_")
          )
        )
      )
    )
  )
}

MULTI_SUBMODULE_GROUP_ORDER <- c("Data", "Biomarker modeling", "Interpretation")
MULTI_SUBMODULE_GROUP_BLURB <- c(
  "Data" = "Data preparation.",
  "Biomarker modeling" = "Rank candidate multi-omics biomarkers.",
  "Interpretation" = "Biomarker interpretation."
)

MULTIOMICS_SIDEBAR_NAV <- list(
  list(id = "dataset", label = "Dataset Workspace", icon = "database", match = "Dataset"),
  list(id = "submodules", label = "Sub-modules", icon = "layer-group", match = "Sub-modules")
)

multiomicsUI <- function() {
  fluidRow(
    column(3, div(class = "omics-sidebar-col", omics_sidebar(
      "multiomics", "Multi-Omics", MULTIOMICS_SIDEBAR_NAV,
      extra_sidebar_content = arthochat_shortcut_ui("Questions about DIABLO, SNF, or how these panels were built? Ask ArthOChat.", compact = TRUE),
      dynamic_nav_output_id = "mo_sidebar_dynamic_nav"
    ))),
    column(
      9,
      div(
        class = "page-header page-header-tight",
        div(class = "page-header-pattern"),
        h2(icon("layer-group"), " Multi-Omics"),
        uiOutput("mo_page_subtitle")
      ),
      div(
        class = "tx-menu-wrap",
        tabsetPanel(
          id = "mo_menu", type = "tabs",
          tabPanel("Dataset", br(), mod_multi_dataset_ui("mo_dataset")),
          tabPanel(
            "Sub-modules", br(),
            div(
              class = "sm-toolbar",
              div(class = "sm-toolbar-count", uiOutput("mo_sm_active_count")),
              div(
                class = "sm-toolbar-search",
                tags$span(icon("magnifying-glass")),
                textInput("mo_sm_search", NULL, placeholder = "Filter sub-modules by name...", width = "260px")
              )
            ),
            build_submodule_grid(MULTI_MODULES, MULTI_SUBMODULE_GROUP_ORDER, MULTI_SUBMODULE_GROUP_BLURB, id_prefix = "mo_")
          )
        )
      )
    )
  )
}

comingSoonUI <- function(title, blurb) {
  tagList(
    div(class = "page-header", h2(title)),
    box(
      width = 12, status = "primary", solidHeader = FALSE,
      div(
        class = "coming-soon",
        icon("hammer", class = "coming-soon-icon"),
        h4("This module is not built yet"),
        p(blurb)
      )
    )
  )
}

addCssDepsOnly <- function(tag) {
  deps <- htmltools::findDependencies(shinydashboard:::addDeps(tags$div()))
  js_free <- lapply(deps, function(d) { d$script <- character(0); d })
  htmltools::attachDependencies(tag, js_free, append = TRUE)
}

existing_app_ui <<- function(user_email) {
  addCssDepsOnly(
    fluidPage(
      theme = shinythemes::shinytheme("cerulean"),
      app_header(user_email),
      navbarPage(
        title = "", id = "sidebar_tabs", selected = "home",
        tabPanel(tagList(icon("house"), "Home"), value = "home", homeUI()),
        tabPanel(tagList(icon("table-cells-large"), "Modules"), value = "modules", modulesLandingUI()),
        tabPanel(tagList(icon("dna"), "Transcriptomics"), value = "transcriptomics", transcriptomicsUI()),
        tabPanel(tagList(icon("circle-nodes"), "Methylomics"), value = "methylomics", methylomicsUI()),
        tabPanel(tagList(icon("arrows-left-right"), "Cross-Omics"), value = "crossomics", crossomicsUI()),
        tabPanel(tagList(icon("layer-group"), "Multi-Omics"), value = "multiomics", multiomicsUI())
      ),
      div(
        id = "arthochat_drawer", class = "chat-drawer",
        div(
          class = "chat-drawer-header",
          div(icon("comments"), strong(" ArthOChat")),
          tags$button(
            icon("xmark"), class = "chat-drawer-close", title = "Close",
            onclick = "document.getElementById('arthochat_drawer').classList.remove('open')"
          )
        ),
        div(class = "chat-drawer-body", mod_arthochat_ui("arthochat"))
      ),
      tags$div(
        id = "arthochat_drawer_backdrop", class = "chat-drawer-backdrop",
        onclick = "document.getElementById('arthochat_drawer').classList.remove('open')"
      )
    )
  )
}

ui <- function(request) {
  addCssDepsOnly(
    fluidPage(
      theme = shinythemes::shinytheme("cerulean"),
      useShinyjs(),
      tags$head(
        tags$title("ArthOMix"),
        tags$script(HTML(
          "(function(){
             try {
               var t = localStorage.getItem('arthomix-theme') || 'light';
               document.documentElement.setAttribute('data-theme', t);
             } catch (e) {}
           })();"
        )),
        tags$link(rel = "stylesheet", type = "text/css",
                  href = paste0("custom.css?v=", as.integer(file.mtime("www/custom.css")))),
        tags$link(rel = "stylesheet", type = "text/css",
                  href = paste0("menuhex.css?v=", as.integer(file.mtime("www/menuhex.css")))),
        tags$link(rel = "stylesheet", type = "text/css",
                  href = paste0("auth.css?v=", as.integer(file.mtime("www/auth.css")))),
        tags$link(rel = "stylesheet", type = "text/css",
                  href = paste0("dark-theme.css?v=", as.integer(file.mtime("www/dark-theme.css"))))
      ),
      uiOutput("app_shell")
    )
  )
}
