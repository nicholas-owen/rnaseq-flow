# Deck builder

Builds the introductory talk from [deck/](deck/), where the slides are authored.

```bash
cd scripts/slides
npm install      # one-time: pptxgenjs, jszip and playwright-core
npm run deck
```

Writes three files to the repository root:

| File | What it is |
|---|---|
| `rnaseq-flow_slide_deck.html` | The talk in a browser, with speaker notes. Tracked |
| `rnaseq-flow_slide_deck.pptx` | Editable PowerPoint, every element placed where the browser draws it. Gitignored |
| `rnaseq-flow_slide_deck_images.pptx` | One picture per slide, identical to the HTML deck. Gitignored |

The build drives Microsoft Edge on Windows, Chrome elsewhere, or the browser named in `CHROMIUM_PATH`. How to change words, notes, slides and the look, what the build warns about, and how to keep the Slides artifact in step are in [deck/README.md](deck/README.md).

## The previous deck

`rnaseq-flow.html`, the deck this one replaces, and `build_slides.js`, which built its PDF and screens PowerPoint, were moved to `working_docs/archive-decks/` on 2026-11-16. Neither is built any more.
