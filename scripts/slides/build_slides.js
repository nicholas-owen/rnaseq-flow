// =============================================================================
// rnaseq-flow deck builder
// -----------------------------------------------------------------------------
// Builds the PDF and the PowerPoint rendering of the deck from rnaseq-flow.html,
// which is the authored source. Neither output holds slide content of its own,
// so the formats cannot drift apart the way the May 2026 .pptx pair did: those
// were a 1.1.0 snapshot still being handed out while the HTML deck moved on
// seven releases.
//
//   cd scripts/slides
//   npm install        # one-time: cheerio and pptxgenjs
//   npm run build      # both outputs
//   npm run pdf        # PDF only
//   npm run screens    # PowerPoint only, about a minute
//
// Writes rnaseq-flow_Slides.pdf and rnaseq-flow_Slides_screens.pptx to the
// repository root. Both are gitignored: the deck they come from is the tracked
// file.
//
// There is no build that reconstructs the slides as editable PowerPoint shapes.
// There was, and it was dropped on 2026-11-01: rebuilding a designed slide from
// its text gives you something that looks like the deck rather than the deck,
// the workflow map cannot survive it at all, and a deck that needs design work
// is better taken back to the tool that designed it. rnaseq-flow-editable.pptx
// at the repository root is that editable copy, maintained by hand.
// =============================================================================

const fs = require('fs');
const path = require('path');
const os = require('os');
const { spawnSync } = require('child_process');
const cheerio = require('cheerio');
const PptxGenJS = require('pptxgenjs');

const REPO_ROOT = path.join(__dirname, '..', '..');
const DECK = path.join(REPO_ROOT, 'rnaseq-flow.html');
const OUT_PDF = path.join(REPO_ROOT, 'rnaseq-flow_Slides.pdf');
const OUT_SCREENS = path.join(REPO_ROOT, 'rnaseq-flow_Slides_screens.pptx');

// The deck is authored at 1920x1080. 144 px to the inch makes that 13.333 x 7.5
// inches, the 16:9 slide PowerPoint uses.
const PX = 144;
const DECK_W = 1920;
const DECK_H = 1080;
const inch = (px) => px / PX;

const squash = (s) => (s || '').replace(/\s+/g, ' ').trim();

// ---- reading the deck -------------------------------------------------------

function readDeck() {
  const $ = cheerio.load(fs.readFileSync(DECK, 'utf8'));
  let notes = [];
  const raw = $('#speaker-notes').html();
  if (raw) {
    try {
      notes = JSON.parse(raw.trim());
    } catch (e) {
      console.warn('  speaker notes are not valid JSON, skipping them:', e.message);
    }
  }

  const slides = [];
  $('deck-stage > section').each((i, sec) => {
    const $sec = $(sec);
    slides.push({
      label: $sec.attr('data-label') || '',
      title: squash($sec.find('.slide-title, .title-name').first().text()),
      notes: notes[i] || '',
    });
  });

  if (notes.length && notes.length !== slides.length) {
    console.warn(`  ${slides.length} slides but ${notes.length} speaker notes: ` +
                 'the notes are matched by position, so one of them is against the wrong slide.');
  }
  return slides;
}

// ---- the browser ------------------------------------------------------------

function findBrowser() {
  if (process.env.CHROME_PATH && fs.existsSync(process.env.CHROME_PATH)) {
    return process.env.CHROME_PATH;
  }
  const candidates = process.platform === 'win32'
    ? [
      'C:/Program Files/Google/Chrome/Application/chrome.exe',
      'C:/Program Files (x86)/Google/Chrome/Application/chrome.exe',
      'C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe',
      'C:/Program Files/Microsoft/Edge/Application/msedge.exe',
      path.join(os.homedir(), 'AppData/Local/Google/Chrome/Application/chrome.exe'),
    ]
    : [
      '/usr/bin/google-chrome', '/usr/bin/google-chrome-stable',
      '/usr/bin/chromium', '/usr/bin/chromium-browser',
      '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
      '/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge',
    ];
  return candidates.find((p) => fs.existsSync(p)) || null;
}

// ---- PDF --------------------------------------------------------------------

function buildPdf() {
  const browser = findBrowser();
  if (!browser) {
    console.error('  no Chrome or Edge found. Set CHROME_PATH to the executable, or open');
    console.error('  rnaseq-flow.html in a browser and print to PDF: the page carries the');
    console.error('  @page rules that put one slide on each landscape page.');
    return false;
  }
  const args = [
    '--headless=new',
    '--disable-gpu',
    '--no-pdf-header-footer',
    '--run-all-compositor-stages-before-draw',
    '--virtual-time-budget=5000',
    `--print-to-pdf=${OUT_PDF}`,
    `file:///${DECK.replace(/\\/g, '/')}`,
  ];
  const r = spawnSync(browser, args, { encoding: 'utf8' });
  if (r.error || !fs.existsSync(OUT_PDF)) {
    console.error('  print-to-pdf failed:', r.error ? r.error.message : (r.stderr || '').trim());
    return false;
  }
  console.log(`  ${path.basename(browser)} wrote ${OUT_PDF} ` +
              `(${fs.statSync(OUT_PDF).size} bytes)`);
  return true;
}

// ---- PowerPoint -------------------------------------------------------------
//
// Chrome renders each slide at 1920x1080 and the image goes on the slide full
// bleed, so the result matches the design exactly, SVG diagrams included, and
// nothing on it can be edited. Speaker notes come across, so it presents.
//
// Same idea as the old rnaseq-flow-screens.pptx, except this one is rebuilt
// from the deck in a minute rather than exported by hand once in May.

function shotSlide(browser, index, outPng) {
  const html = fs.readFileSync(DECK, 'utf8').replace(
    '</head>',
    `<style>deck-stage > section { display: none !important; }
      deck-stage > section:nth-of-type(${index + 1}) {
        display: block !important; margin: 0 !important; box-shadow: none !important; }
    </style></head>`);
  const tmp = path.join(os.tmpdir(), `rnaseq-flow-slide-${index + 1}.html`);
  fs.writeFileSync(tmp, html, 'utf8');
  const r = spawnSync(browser, [
    '--headless=new', '--disable-gpu', '--hide-scrollbars',
    '--window-size=1920,1080',
    '--run-all-compositor-stages-before-draw',
    '--virtual-time-budget=4000',
    `--screenshot=${outPng}`,
    `file:///${tmp.replace(/\\/g, '/')}`,
  ], { encoding: 'utf8' });
  fs.unlinkSync(tmp);
  return !r.error && fs.existsSync(outPng);
}

async function buildScreensPptx(slides) {
  const browser = findBrowser();
  if (!browser) {
    console.error('  no Chrome or Edge found, so the PowerPoint build is skipped.');
    return;
  }
  const shots = fs.mkdtempSync(path.join(os.tmpdir(), 'rnaseq-flow-shots-'));
  const pptx = new PptxGenJS();
  pptx.defineLayout({ name: 'DECK', width: inch(DECK_W), height: inch(DECK_H) });
  pptx.layout = 'DECK';
  pptx.title = 'rnaseq-flow: introduction';
  pptx.author = 'rnaseq-flow';
  pptx.company = 'rnaseq-flow';

  let made = 0;
  slides.forEach((s, i) => {
    const png = path.join(shots, `slide-${i + 1}.png`);
    const slide = pptx.addSlide();
    if (s.notes) slide.addNotes(s.notes);
    if (shotSlide(browser, i, png)) {
      slide.addImage({ path: png, x: 0, y: 0, w: inch(DECK_W), h: inch(DECK_H) });
      made += 1;
    } else {
      slide.addText(`Slide ${i + 1} could not be rendered: ${s.title || s.label}`, {
        x: 1, y: 3, w: inch(DECK_W) - 2, h: 1,
        fontFace: 'Segoe UI', fontSize: 16, color: '868E96',
      });
    }
  });

  await pptx.writeFile({ fileName: OUT_SCREENS });
  console.log(`  wrote ${OUT_SCREENS} (${fs.statSync(OUT_SCREENS).size} bytes, ${made} rendered)`);
  fs.rmSync(shots, { recursive: true, force: true });
}

// ---- main -------------------------------------------------------------------

async function main() {
  const args = process.argv.slice(2);
  const wantPdf = args.length === 0 || args.includes('--pdf');
  const wantScreens = args.length === 0 || args.includes('--screens') || args.includes('--images');

  if (!fs.existsSync(DECK)) throw new Error(`deck not found: ${DECK}`);
  const slides = readDeck();
  console.log(`read ${slides.length} slides from ${path.basename(DECK)}`);

  if (wantPdf) buildPdf();
  if (wantScreens) await buildScreensPptx(slides);
}

main().catch((e) => {
  console.error(e.stack || e.message);
  process.exit(1);
});
