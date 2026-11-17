// Build the rnaseq-flow talk from scripts/slides/deck/.
//
// scripts/slides/deck/ holds the deck exactly as the Slides artifact keeps it:
//   deck.json         title, slide order, outline and typefaces
//   deck.css          the shared look (the deck-design-flow design system)
//   slides/<id>.html  one <section> per slide; its <aside> is the speaker notes
// Edit those files, or copy them out of the Slides artifact, then rebuild:
//
//   cd scripts/slides
//   npm run deck
//
// Outputs, all at the repository root:
//   rnaseq-flow_slide_deck.html          the talk in a browser: arrow keys or click to move,
//                                        N speaker notes, F full screen, P print
//   rnaseq-flow_slide_deck.pptx          editable PowerPoint: every text box, shape and table
//                                        is placed where the browser draws it
//   rnaseq-flow_slide_deck_images.pptx   one picture per slide, identical to the HTML deck
//
// Both PowerPoint decks are measured and photographed from the HTML deck in a headless
// browser (playwright-core driving Microsoft Edge on Windows, Chrome elsewhere, or the
// browser named by CHROMIUM_PATH). The browser must reach Google Fonts, or have Roboto
// and Roboto Mono installed, so that the slides are laid out in Roboto; the build warns
// if it could not. PowerPoint draws the editable deck in Roboto and Roboto Mono, so
// install both where it will be shown, or present the images deck.
//
// The same builder, with its own paths, makes the refdata-flow talk.

const fs = require("fs");
const path = require("path");
const { pathToFileURL } = require("url");
const PptxGenJS = require("pptxgenjs");
const JSZip = require("jszip");
const { chromium } = require("playwright-core");

const HERE = __dirname;
const REPO = path.join(HERE, "..", "..");
const SRC = path.join(HERE, "deck");
const OUT_HTML = path.join(REPO, "rnaseq-flow_slide_deck.html");
const OUT_PPTX = path.join(REPO, "rnaseq-flow_slide_deck.pptx");
const OUT_IMAGES = path.join(REPO, "rnaseq-flow_slide_deck_images.pptx");

const W = 1920, H = 1080;          // the slide canvas, in CSS px
const PX_PER_IN = 144;             // 1920 px across a 13.333 in (16:9) PowerPoint slide
const PT_PER_PX = 0.5;             // so 1 px is half a point
const inch = (px) => px / PX_PER_IN;
const pt = (px) => Math.round(px * PT_PER_PX * 100) / 100;

// =====================================================================================
// The deck's files
// =====================================================================================
function readDeck() {
  const index = JSON.parse(fs.readFileSync(path.join(SRC, "deck.json"), "utf8"));
  const css = fs.readFileSync(path.join(SRC, "deck.css"), "utf8");
  const slides = index.order.map((id) => {
    if (!/^[A-Za-z0-9_-]{1,64}$/.test(id)) throw new Error(`bad slide id in deck.json: ${id}`);
    const file = path.join(SRC, "slides", `${id}.html`);
    if (!fs.existsSync(file)) throw new Error(`deck.json names slide "${id}" but ${path.relative(REPO, file)} does not exist`);
    const html = fs.readFileSync(file, "utf8").trim();
    if (!new RegExp(`<section[^>]*\\bid="${id}"`).test(html)) {
      throw new Error(`${path.relative(REPO, file)}: its <section> must have id="${id}"`);
    }
    const m = html.match(/<aside>([\s\S]*?)<\/aside>/);
    const notes = m ? decode(m[1]).trim() : "";
    return { id, html, notes };
  });
  const fonts = Object.values(index.faces || {}).filter((f) => f.href).map((f) => f.href);
  return { title: index.title || "rnaseq-flow", css, slides, fonts };
}

function decode(s) {
  return s.replace(/<[^>]+>/g, "")
    .replace(/&#(\d+);/g, (_, n) => String.fromCodePoint(+n))
    .replace(/&#x([0-9a-f]+);/gi, (_, n) => String.fromCodePoint(parseInt(n, 16)))
    .replace(/&lt;/g, "<").replace(/&gt;/g, ">").replace(/&quot;/g, '"').replace(/&amp;/g, "&");
}
const esc = (s) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");

// =====================================================================================
// HTML deck
// =====================================================================================

// What the Slides artifact draws before the deck's own CSS: no default margins, its
// heading sizes, and tables with their own hairlines.
const BASE_CSS = `
* { margin: 0; }
h1 { font-size: 96px; font-weight: 600; }
h2 { font-size: 64px; font-weight: 600; }
h3 { font-size: 44px; font-weight: 600; }
#stage > section { width: ${W}px; height: ${H}px; overflow: hidden; box-sizing: border-box; position: relative; isolation: isolate; }
section aside { display: none !important; }
table { width: 100%; border-collapse: separate; border-spacing: 0; border: solid rgba(0,0,0,0.2); border-width: 2px 0 0 2px; }
th, td { padding: .35em .6em; text-align: left; vertical-align: top; border: solid rgba(0,0,0,0.2); border-width: 0 2px 2px 0; }
th { font-weight: 700; }
`;

const PLAYER_CSS = `
html, body { height: 100%; }
body { background: #0b0f17; overflow: hidden; }
#stage { position: absolute; left: 0; top: 0; width: ${W}px; height: ${H}px; transform-origin: 0 0; }
#stage > section { position: absolute; left: 0; top: 0; display: none; }
#stage > section.on { display: flex; }
#notes { position: fixed; left: 0; right: 0; bottom: 0; max-height: 40vh; overflow: auto; padding: 16px 24px;
  background: rgba(15,27,45,0.94); color: #f1f5f9; font: 18px/1.5 Roboto, 'Segoe UI', Arial, sans-serif; }
#notes[hidden] { display: none; }
#notes b { color: #94a3b8; font-weight: 500; letter-spacing: 0.1em; }
@media print {
  @page { size: ${W}px ${H}px; margin: 0; }
  html, body { height: auto; overflow: visible; background: none; }
  #stage { position: static; transform: none !important; height: auto; }
  #stage > section { display: flex !important; position: relative; break-after: page; }
  #notes { display: none !important; }
  .rise { animation: none !important; }
}
`;

const PLAYER_JS = `
(function () {
  var stage = document.getElementById('stage');
  var slides = Array.prototype.slice.call(stage.children);
  var notes = document.getElementById('notes');
  var i = 0;
  function fit() {
    var s = Math.min(innerWidth / ${W}, innerHeight / ${H});
    stage.style.transform = 'translate(' + (innerWidth - ${W} * s) / 2 + 'px,' + (innerHeight - ${H} * s) / 2 + 'px) scale(' + s + ')';
  }
  function show(n) {
    i = Math.max(0, Math.min(slides.length - 1, n));
    slides.forEach(function (el, k) { el.classList.toggle('on', k === i); });
    var a = slides[i].querySelector('aside');
    notes.innerHTML = '<b>' + (i + 1) + ' / ' + slides.length + '</b>&nbsp; ';
    notes.appendChild(document.createTextNode(a ? a.textContent : ''));
    if (location.hash !== '#' + (i + 1)) history.replaceState(null, '', '#' + (i + 1));
  }
  addEventListener('keydown', function (e) {
    var k = e.key;
    if (k === 'ArrowRight' || k === 'PageDown' || k === ' ' || k === 'Enter') show(i + 1);
    else if (k === 'ArrowLeft' || k === 'PageUp' || k === 'Backspace') show(i - 1);
    else if (k === 'Home') show(0);
    else if (k === 'End') show(slides.length - 1);
    else if (k === 'n' || k === 'N') notes.hidden = !notes.hidden;
    else if (k === 'f' || k === 'F') { if (document.fullscreenElement) document.exitFullscreen(); else document.documentElement.requestFullscreen(); }
    else if (k === 'p' || k === 'P') print();
    else return;
    e.preventDefault();
  });
  stage.addEventListener('click', function () { show(i + 1); });
  addEventListener('resize', fit);
  addEventListener('hashchange', function () { show((parseInt(location.hash.slice(1), 10) || 1) - 1); });
  fit();
  show((parseInt(location.hash.slice(1), 10) || 1) - 1);
})();
`;

function buildHtml(deck) {
  // CSS counters skip slides that are not shown, so give each slide its number outright.
  const numbers = deck.slides.map((s, k) => `#stage > section#${s.id} { counter-set: page ${k + 1}; }`).join("\n");
  const links = deck.fonts.map((h) => `<link rel="stylesheet" href="${esc(h)}">`).join("\n");
  return `<!doctype html>
<html lang="en-GB">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>${esc(deck.title)}</title>
<!-- Generated by scripts/slides/build_deck.js from scripts/slides/deck/. Edit there and run: cd scripts/slides && npm run deck -->
${links}
<style>${BASE_CSS}
/* ---- scripts/slides/deck/deck.css ---- */
${deck.css}
/* ---- player ---- */
${PLAYER_CSS}
${numbers}
</style>
</head>
<body>
<div id="stage">
${deck.slides.map((s) => s.html).join("\n")}
</div>
<div id="notes" hidden></div>
<script>${PLAYER_JS}</script>
</body>
</html>
`;
}

// =====================================================================================
// Measuring the HTML deck
// =====================================================================================

// Shows every slide at full size, one under another, with no motion.
const MEASURE_CSS = `
html, body { overflow: visible !important; height: auto !important; background: #fff !important; }
#stage { position: static !important; transform: none !important; height: auto !important; }
#stage > section { display: flex !important; position: relative !important; }
#notes { display: none !important; }
*, *::before, *::after { animation: none !important; transition: none !important; }
`;

// Runs in the browser. Returns, for one slide, its background, the boxes, lines, text
// and tables to draw, in paint order, in px from the slide's top-left corner, and any
// element whose bottom edge passes the slide's bottom margin.
function measureSlide(index) {
  const sec = document.querySelectorAll("#stage > section")[index];
  const origin = sec.getBoundingClientRect();
  const items = [];

  function rgba(c) {
    const m = c && c.match(/rgba?\(([^)]+)\)/);
    if (!m) return null;
    const p = m[1].split(/[ ,/]+/).filter(Boolean).map(Number);
    return { r: p[0], g: p[1], b: p[2], a: p.length > 3 ? p[3] : 1 };
  }
  const hex = (c) => [c.r, c.g, c.b].map((v) => Math.round(v).toString(16).padStart(2, "0")).join("").toUpperCase();
  const slideBg = rgba(getComputedStyle(sec).backgroundColor);
  // A translucent colour, mixed with the slide's background as PowerPoint has no alpha here.
  function solid(c) {
    if (!c || c.a === 0) return null;
    if (c.a >= 1) return hex(c);
    return hex({ r: c.r * c.a + slideBg.r * (1 - c.a), g: c.g * c.a + slideBg.g * (1 - c.a), b: c.b * c.a + slideBg.b * (1 - c.a) });
  }
  function rect(el) {
    const r = el.getBoundingClientRect();
    return { x: r.left - origin.left, y: r.top - origin.top, w: r.width, h: r.height };
  }
  const n = (v) => parseFloat(v) || 0;

  function runStyle(el) {
    const cs = getComputedStyle(el);
    return {
      font: cs.fontFamily.split(",")[0].replace(/["']/g, "").trim(),
      size: n(cs.fontSize),
      weight: n(cs.fontWeight),
      italic: cs.fontStyle === "italic",
      underline: cs.textDecorationLine.includes("underline"),
      color: solid(rgba(cs.color)),
      spacing: cs.letterSpacing === "normal" ? 0 : n(cs.letterSpacing),
    };
  }
  // The words of an element as runs, a <br> as a line break.
  function runs(el) {
    const out = [];
    (function walk(node) {
      for (const c of node.childNodes) {
        if (c.nodeType === 3) {
          const t = c.textContent.replace(/[ \t\n\r\f]+/g, " ");
          if (t) out.push({ text: t, ...runStyle(c.parentElement) });
        } else if (c.nodeType === 1) {
          const tag = c.tagName.toLowerCase();
          if (tag === "br") out.push({ text: "", br: true });
          else if (tag !== "aside" && getComputedStyle(c).display !== "none") walk(c);
        }
      }
    })(el);
    const before = getComputedStyle(el, "::before").content;
    if (/counter\(page\)/.test(before) || (before && before !== "none" && before !== "normal" && !out.length)) {
      const text = /counter/.test(before) ? String(index + 1) : before.replace(/^["']|["']$/g, "");
      out.unshift({ text, ...runStyle(el) });
    }
    if (out.length) {
      out[0].text = out[0].text.replace(/^ +/, "");
      out[out.length - 1].text = out[out.length - 1].text.replace(/ +$/, "");
    }
    return out.filter((r) => r.text || r.br);
  }
  function para(el) {
    const cs = getComputedStyle(el);
    const lh = cs.lineHeight === "normal" ? 1.2 * n(cs.fontSize) : n(cs.lineHeight);
    let align = cs.textAlign;
    if (align === "start" || align === "justify") align = "left";
    if (align === "end") align = "right";
    return { align, lineHeight: lh, nowrap: cs.whiteSpace === "nowrap" || cs.whiteSpace === "pre" };
  }
  function contentBox(el) {
    const cs = getComputedStyle(el), r = rect(el);
    const l = n(cs.paddingLeft) + n(cs.borderLeftWidth), t = n(cs.paddingTop) + n(cs.borderTopWidth);
    const rr = n(cs.paddingRight) + n(cs.borderRightWidth), b = n(cs.paddingBottom) + n(cs.borderBottomWidth);
    return { x: r.x + l, y: r.y + t, w: r.w - l - rr, h: r.h - t - b };
  }

  function paint(el, cs) {
    const r = rect(el);
    const fill = solid(rgba(cs.backgroundColor));
    const sides = ["Top", "Right", "Bottom", "Left"].map((s) => ({
      w: cs[`border${s}Style`] === "none" ? 0 : n(cs[`border${s}Width`]), c: solid(rgba(cs[`border${s}Color`])),
    }));
    const even = sides.every((s) => s.w === sides[0].w && s.c === sides[0].c);
    const radius = Math.min(n(cs.borderTopLeftRadius), r.w / 2, r.h / 2);
    if (fill || (even && sides[0].w > 0 && sides[0].c)) {
      items.push({ kind: "box", ...r, fill, radius, line: even && sides[0].w > 0 && sides[0].c ? { color: sides[0].c, w: sides[0].w } : null });
    }
    if (!even) {
      const [t, rt, b, l] = sides;
      if (t.w && t.c) items.push({ kind: "line", x1: r.x, y1: r.y + t.w / 2, x2: r.x + r.w, y2: r.y + t.w / 2, color: t.c, w: t.w });
      if (b.w && b.c) items.push({ kind: "line", x1: r.x, y1: r.y + r.h - b.w / 2, x2: r.x + r.w, y2: r.y + r.h - b.w / 2, color: b.c, w: b.w });
      if (l.w && l.c) items.push({ kind: "line", x1: r.x + l.w / 2, y1: r.y, x2: r.x + l.w / 2, y2: r.y + r.h, color: l.c, w: l.w });
      if (rt.w && rt.c) items.push({ kind: "line", x1: r.x + r.w - rt.w / 2, y1: r.y, x2: r.x + r.w - rt.w / 2, y2: r.y + r.h, color: rt.c, w: rt.w });
    }
  }

  function svg(el) {
    const r = rect(el);
    const vb = (el.getAttribute("viewBox") || `0 0 ${r.w} ${r.h}`).split(/[ ,]+/).map(Number);
    const sx = r.w / vb[2], sy = r.h / vb[3];
    const X = (v) => r.x + (v - vb[0]) * sx, Y = (v) => r.y + (v - vb[1]) * sy;
    for (const s of el.querySelectorAll("*")) {
      const cs = getComputedStyle(s);
      const stroke = cs.stroke && cs.stroke !== "none" ? solid(rgba(cs.stroke)) : null;
      const sw = n(cs.strokeWidth) * sx;
      const tag = s.tagName.toLowerCase();
      if (tag === "line") {
        items.push({ kind: "route", x1: X(n(s.getAttribute("x1"))), y1: Y(n(s.getAttribute("y1"))),
          x2: X(n(s.getAttribute("x2"))), y2: Y(n(s.getAttribute("y2"))), color: stroke, w: sw, round: cs.strokeLinecap === "round" });
      } else if (tag === "circle") {
        const rad = n(s.getAttribute("r")) * sx;
        items.push({ kind: "circle", cx: X(n(s.getAttribute("cx"))), cy: Y(n(s.getAttribute("cy"))), r: rad,
          fill: cs.fill && cs.fill !== "none" ? solid(rgba(cs.fill)) : null, stroke, w: sw });
      } else if (tag !== "title" && tag !== "desc") {
        items.push({ kind: "unsupported", what: `<${tag}> in an svg` });
      }
    }
  }

  function table(el) {
    const r = rect(el);
    const rows = Array.from(el.rows);
    const first = Array.from(rows[0].cells).map((c) => rect(c).w);
    const td = rows[0].cells[0], cs = getComputedStyle(td);
    items.push({
      kind: "table", x: r.x, y: r.y, w: r.w, colW: first,
      rowH: rows.map((tr) => rect(tr).h),
      pad: [n(cs.paddingTop), n(cs.paddingRight), n(cs.paddingBottom), n(cs.paddingLeft)],
      border: { color: solid(rgba(cs.borderBottomColor)), w: n(cs.borderBottomWidth) },
      lineHeight: para(td).lineHeight,
      rows: rows.map((tr) => ({
        fill: solid(rgba(getComputedStyle(tr).backgroundColor)),
        cells: Array.from(tr.cells).map((c) => ({ runs: runs(c), align: para(c).align })),
      })),
    });
  }

  function list(el, cs) {
    const box = contentBox(el);
    const indent = n(cs.paddingLeft);
    const lis = Array.from(el.children).filter((c) => c.tagName.toLowerCase() === "li");
    items.push({
      kind: "text", x: box.x - indent, y: box.y, w: box.w + indent, h: box.h,
      paras: lis.map((li) => ({ ...para(li), runs: runs(li), bullet: { indent, numbered: el.tagName.toLowerCase() === "ol" } })),
    });
  }

  const TEXT_TAGS = new Set(["p", "h1", "h2", "h3", "h4", "h5", "h6", "li"]);
  function walk(el) {
    const cs = getComputedStyle(el);
    const tag = el.tagName.toLowerCase();
    if (tag === "aside" || cs.display === "none" || cs.visibility === "hidden") return;
    if (tag === "svg") return svg(el);
    if (tag === "table") return table(el);
    if (el !== sec) paint(el, cs);
    if (tag === "ul" || tag === "ol") return list(el, cs);
    const ownText = Array.from(el.childNodes).some((c) => c.nodeType === 3 && c.textContent.trim());
    if (TEXT_TAGS.has(tag) || ownText) {
      const rs = runs(el);
      if (rs.length) items.push({ kind: "text", ...contentBox(el), paras: [{ ...para(el), runs: rs }] });
      return;
    }
    if (tag.startsWith("x-") || tag === "img" || tag === "video") items.push({ kind: "unsupported", what: `<${tag}>` });
    for (const c of el.children) walk(c);
  }
  walk(sec);

  // Content in the flow that runs past the bottom margin, or off the slide.
  const limit = origin.height - n(getComputedStyle(sec).paddingBottom);
  const overflow = [];
  const pinned = (el) => {
    for (let a = el; a && a !== sec; a = a.parentElement) if (getComputedStyle(a).position === "absolute") return true;
    return false;
  };
  for (const el of sec.querySelectorAll("*")) {
    if (el.closest("aside") || el.closest("svg") || pinned(el)) continue;
    const bottom = rect(el).y + rect(el).h;
    if (bottom > limit + 1 && el.children.length === 0) overflow.push({ tag: el.tagName.toLowerCase(), by: Math.round(bottom - limit) });
  }
  return { bg: solid(slideBg) || "FFFFFF", items, overflow };
}

async function launch() {
  if (process.env.CHROMIUM_PATH) return chromium.launch({ executablePath: process.env.CHROMIUM_PATH });
  const channels = process.platform === "win32" ? ["msedge", "chrome"] : ["chrome", "msedge"];
  for (const channel of channels) {
    try { return await chromium.launch({ channel }); } catch (e) { /* try the next */ }
  }
  throw new Error("No browser found: install Microsoft Edge or Google Chrome, or set CHROMIUM_PATH to a Chromium executable.");
}

async function measure(deck) {
  const browser = await launch();
  try {
    const page = await browser.newPage({ viewport: { width: W, height: H }, deviceScaleFactor: 2 });
    await page.goto(pathToFileURL(OUT_HTML).href, { waitUntil: "networkidle" });
    await page.addStyleTag({ content: MEASURE_CSS });
    await page.evaluate(() => document.fonts.ready);
    const missing = await page.evaluate(() =>
      ["400 32px Roboto", "700 32px Roboto", "800 32px Roboto", "400 28px 'Roboto Mono'"].filter((f) => !document.fonts.check(f)));
    if (missing.length) {
      console.warn(`WARNING: these faces did not load, so the slides were measured in a fallback font: ${missing.join(", ")}.`);
      console.warn("         The PowerPoint layout will not match the deck. Rebuild with access to fonts.googleapis.com.");
    }
    const sections = await page.$$("#stage > section");
    const out = [];
    for (let k = 0; k < sections.length; k++) {
      const m = await page.evaluate(measureSlide, k);
      const png = await sections[k].screenshot({ type: "png" });
      out.push({ ...m, png, notes: deck.slides[k].notes, id: deck.slides[k].id });
    }
    return out;
  } finally {
    await browser.close();
  }
}

// =====================================================================================
// PowerPoint
// =====================================================================================
function newPres(deck) {
  const pres = new PptxGenJS();
  pres.layout = "LAYOUT_WIDE";        // 13.333 x 7.5 in
  pres.author = "Nicholas Owen";
  pres.company = "UCL Advanced Research Computing";
  pres.title = deck.title;
  return pres;
}

const FACES = { "Roboto Mono": "Roboto Mono", Roboto: "Roboto" };
function runOpts(r) {
  return {
    fontFace: FACES[r.font] || r.font,
    fontSize: pt(r.size),
    bold: r.weight >= 600,
    italic: r.italic,
    underline: r.underline ? { style: "sng" } : undefined,
    color: r.color || "000000",
    charSpacing: r.spacing ? pt(r.spacing) : undefined,
  };
}
function textRuns(paras) {
  const out = [];
  paras.forEach((p, pi) => {
    const last = pi === paras.length - 1;
    const pOpts = {
      align: p.align, lineSpacing: pt(p.lineHeight), paraSpaceBefore: 0, paraSpaceAfter: 0,
      bullet: p.bullet ? (p.bullet.numbered ? { type: "number", indent: pt(p.bullet.indent) } : { indent: pt(p.bullet.indent) }) : false,
    };
    p.runs.forEach((r, ri) => {
      if (r.br) { if (out.length) out[out.length - 1].options.breakLine = true; return; }
      out.push({ text: r.text, options: { ...pOpts, ...runOpts(r), breakLine: !last && ri === p.runs.length - 1 } });
    });
  });
  return out;
}

const SLACK = 6;   // px added to a text box's width, so a line PowerPoint measures a hair wider still fits
function drawItem(slide, it, warn) {
  switch (it.kind) {
    case "box":
      slide.addShape(it.radius > 0 ? "roundRect" : "rect", {
        x: inch(it.x), y: inch(it.y), w: inch(it.w), h: inch(it.h),
        fill: it.fill ? { color: it.fill } : { type: "none" },
        line: it.line ? { color: it.line.color, width: pt(it.line.w) } : undefined,
        rectRadius: it.radius > 0 ? inch(it.radius) : undefined,
      });
      break;
    case "line":
      slide.addShape("line", {
        x: inch(Math.min(it.x1, it.x2)), y: inch(Math.min(it.y1, it.y2)),
        w: inch(Math.abs(it.x2 - it.x1)), h: inch(Math.abs(it.y2 - it.y1)),
        line: { color: it.color, width: pt(it.w) },
      });
      break;
    case "route": {
      // A route is a horizontal or vertical line with round caps: drawn as a pill.
      const half = it.round ? it.w / 2 : 0;
      const x = Math.min(it.x1, it.x2) - (it.y1 === it.y2 ? half : it.w / 2);
      const y = Math.min(it.y1, it.y2) - (it.x1 === it.x2 ? half : it.w / 2);
      const w = it.y1 === it.y2 ? Math.abs(it.x2 - it.x1) + 2 * half : it.w;
      const h = it.x1 === it.x2 ? Math.abs(it.y2 - it.y1) + 2 * half : it.w;
      if (it.x1 !== it.x2 && it.y1 !== it.y2) warn("a diagonal svg line, drawn as its bounding box");
      slide.addShape(it.round ? "roundRect" : "rect", {
        x: inch(x), y: inch(y), w: inch(w), h: inch(h), fill: { color: it.color },
        rectRadius: it.round ? inch(it.w / 2) : undefined,
      });
      break;
    }
    case "circle":
      slide.addShape("ellipse", {
        x: inch(it.cx - it.r), y: inch(it.cy - it.r), w: inch(2 * it.r), h: inch(2 * it.r),
        fill: it.fill ? { color: it.fill } : { type: "none" },
        line: it.stroke ? { color: it.stroke, width: pt(it.w) } : undefined,
      });
      break;
    case "text": {
      const nowrap = it.paras.every((p) => p.nowrap);
      const right = it.paras.every((p) => p.align === "right");
      slide.addText(textRuns(it.paras), {
        x: inch(right ? it.x - SLACK : it.x), y: inch(it.y), w: inch(it.w + SLACK), h: inch(it.h),
        margin: 0, valign: "top", fit: "none", wrap: !nowrap,
      });
      break;
    }
    case "table": {
      const rows = it.rows.map((row, ri) => row.cells.map((c) => {
        const rs = c.runs.map((r) => ({ text: r.text, options: { ...runOpts(r), bold: ri === 0 || r.weight >= 600 } }));
        return {
          text: rs.length ? rs : "",
          options: { fill: row.fill ? { color: row.fill } : undefined, align: c.align, valign: "top", lineSpacing: pt(it.lineHeight) },
        };
      }));
      slide.addTable(rows, {
        x: inch(it.x), y: inch(it.y), w: inch(it.w),
        colW: it.colW.map(inch), rowH: it.rowH.map(inch),
        margin: it.pad.map(inch),
        border: { type: "solid", pt: pt(it.border.w), color: it.border.color || "CCCCCC" },
        autoPage: false,
      });
      break;
    }
    case "unsupported":
      warn(`${it.what} is not drawn in the editable deck`);
      break;
  }
}

// pptxgenjs writes an entry for every folder in the package. PowerPoint does not mind,
// but a cloud folder (this repository lives in Google Drive) can fail to open or co-edit
// a file that has them. Rewrite the package with the files alone.
async function writeDeck(pres, file) {
  await pres.writeFile({ fileName: file });
  const zip = await JSZip.loadAsync(fs.readFileSync(file));
  const clean = new JSZip();
  for (const [name, entry] of Object.entries(zip.files)) {
    if (!entry.dir) clean.file(name, await entry.async("nodebuffer"), { createFolders: false });
  }
  fs.writeFileSync(file, await clean.generateAsync({ type: "nodebuffer", compression: "DEFLATE" }));
  console.log("written", path.relative(REPO, file));
}

// =====================================================================================
async function main() {
  const deck = readDeck();
  fs.writeFileSync(OUT_HTML, buildHtml(deck));
  console.log("written", path.relative(REPO, OUT_HTML));

  const measured = await measure(deck);

  const editable = newPres(deck);
  const pictures = newPres(deck);
  for (const m of measured) {
    const warn = (msg) => console.warn(`  slide ${m.id}: ${msg}`);
    if (m.overflow.length) warn(`content runs past the bottom margin by up to ${Math.max(...m.overflow.map((o) => o.by))} px`);
    const s = editable.addSlide();
    s.background = { color: m.bg };
    for (const it of m.items) drawItem(s, it, warn);
    if (m.notes) s.addNotes(m.notes);

    const p = pictures.addSlide();
    p.addImage({ data: "image/png;base64," + m.png.toString("base64"), x: 0, y: 0, w: inch(W), h: inch(H) });
    if (m.notes) p.addNotes(m.notes);
  }
  await writeDeck(editable, OUT_PPTX);
  await writeDeck(pictures, OUT_IMAGES);
}

main().catch((e) => { console.error(e.message || e); process.exit(1); });
