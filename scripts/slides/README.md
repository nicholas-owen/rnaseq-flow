# Deck builder

Builds the PDF and the PowerPoint rendering of the introduction deck from
[rnaseq-flow.html](../../rnaseq-flow.html), which is where the slides are
authored.

```bash
cd scripts/slides
npm install      # one-time: cheerio and pptxgenjs, 42 packages
npm run build    # both outputs
npm run pdf      # PDF only
npm run screens  # PowerPoint only, about a minute
```

Outputs, both written to the repository root and both gitignored:

| File | What it is |
|---|---|
| `rnaseq-flow_Slides.pdf` | Chrome printing the deck, one 1920x1080 landscape page per slide. The deck exactly as designed, SVG diagrams included |
| `rnaseq-flow_Slides_screens.pptx` | The same 29 slides as full-bleed 1920x1080 renders, pixel-identical to the deck, each carrying its speaker note. Nothing on it is editable |

## The editable deck is not built here

`rnaseq-flow-editable.pptx` at the repository root is hand-maintained, through
Claude Design, and is the copy to change when the deck needs design work.

A build that rebuilt each slide as PowerPoint text boxes, tables and shapes
existed briefly and was dropped on 2026-11-01. Reconstructing a designed slide
from its text gives you something that resembles the deck rather than the deck:
the workflow map, a few hundred positioned SVG elements, cannot survive the
conversion at all, and a bespoke composition flattens into a card grid. Slides
that need designing are better taken back to the tool that designed them.

What this builder is for is the other half of that problem. Three copies of the
deck were being handed out and only one was maintained: `rnaseq-flow-editable.pptx`
and `rnaseq-flow-screens.pptx` were built on 27 May 2026 and say **version
1.1.0** on their title slide, while the pipeline reached 1.8.0 and they had
missed `--counts`, `--validate_only`, `--reference_level`, the `reproduce/`
folders, the RO-Crate work, the GTF-derived transcriptome and the shared
significance thresholds. The PDF and the rendered PowerPoint now come from the
deck itself, so neither can go stale on its own.

## The PDF needs a browser

Both outputs shell out to Chrome or Edge, looking in the usual install
locations. Set `CHROME_PATH` if yours is somewhere else. For the PDF alone you
can skip the builder entirely: open `rnaseq-flow.html` and print it, since the
page carries the `@page` and `@media print` rules that put one slide on each
landscape page.

## Editing the deck

Edit `rnaseq-flow.html`, then rebuild. A new slide is a `<section>` inside
`<deck-stage>` with a `.frame` holding `.eyebrow`, `.slide-title`,
`.slide-subtitle` and its content.

Speaker notes live in the `#speaker-notes` JSON array at the top of the file and
are matched to slides **by position**. Add a slide, add its note at the same
index. The builder warns when the two counts disagree, which is how the gap from
1.6.0 was found: that slide never got a note, so every later note had been
showing against the wrong slide since September.
