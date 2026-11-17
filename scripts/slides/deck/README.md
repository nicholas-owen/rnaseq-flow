# The rnaseq-flow talk: how to update it

This folder is the source of the talk. The HTML deck and both PowerPoint decks are generated from it by `scripts/slides/build_deck.js`, so edit here, never the outputs.

## What is here

| File | What it holds |
|---|---|
| `deck.json` | The deck's title, the slide order (`order`), the outline (`sections`) and the typefaces (`faces`). |
| `deck.css` | The shared look, from the deck-design-flow design system: colours, type sizes, cards, code panels, route lines, footer. |
| `slides/<id>.html` | One slide each: exactly one `<section id="<id>">`, its id the same as the file name. The `<aside>` at the end of the section is the speaker notes. |

These are the same files the Slides artifact "rnaseq-flow: a FAIR approach" keeps, so the two can be copied either way.

`npm run deck` writes three files to the repository root:

| Output | What it is |
|---|---|
| `rnaseq-flow_slide_deck.html` | The talk in a browser. Arrow keys, Space or a click move on; N shows the speaker notes; F goes full screen; P prints one slide per page. `#5` at the end of the address opens slide 5. Tracked in git. |
| `rnaseq-flow_slide_deck.pptx` | Editable PowerPoint. Every text box, shape, route line and table is placed where the browser draws it, with the speaker notes. |
| `rnaseq-flow_slide_deck_images.pptx` | One picture per slide, with the speaker notes. It looks exactly like the HTML deck; nothing in it can be edited. |

Both `.pptx` files are gitignored by the existing `*.pptx` rule. Rebuild them from this folder.

## Setting up, once

1. Install [Node.js](https://nodejs.org/) 18 or later.
2. In `scripts/slides/`, run `npm install`. This brings `pptxgenjs`, `jszip` and `playwright-core`; it downloads no browser.
3. The build drives Microsoft Edge, which Windows already has. Elsewhere it uses Google Chrome, or the Chromium named in the `CHROMIUM_PATH` environment variable.
4. To present the editable deck, install [Roboto](https://fonts.google.com/specimen/Roboto) and [Roboto Mono](https://fonts.google.com/specimen/Roboto+Mono) on that computer. The decks do not carry the fonts, and without them PowerPoint substitutes its own. The images deck needs nothing installed; the HTML deck loads the fonts from Google Fonts, so it needs a connection.

## Changing the talk

### Words

Open `slides/<id>.html` and change the text between the tags. Leave the tags and their `class` attributes alone: the classes give each line its size and colour. `slides/aligners.html` is the slide whose id is `aligners`; the order is in `deck.json`.

No em dashes. The deck says nothing that ties rnaseq-flow to one institution's computers: it runs on laptops, desktops, clusters, the cloud and TREs.

### Speaker notes

Edit the plain text inside the slide's `<aside>...</aside>`. It must stay the last thing in the `<section>`.

### The look

A change to `deck.css` changes every slide that uses that class. The colour names at the top (`--brand-navy`, `--route-raspberry` and so on) are the design system's tokens; use them rather than new hex values. rnaseq-flow's own route colours are raspberry, violet, blue and amber. To change one slide only, put a `style="..."` on the element in that slide's file.

### Adding, removing or reordering slides

To add one, copy the slide most like it, give the file a new name (letters, digits, `-` and `_` only), set the `id` on its `<section>` to the same name, and add that name to `order` in `deck.json`.

To remove one, take its name out of `order`. A slide not in `order` is not built, so the file can stay.

To reorder, reorder `order`. Page numbers follow it. If a slide starts a new part of the talk, change the matching `start` in `sections` too.

The build stops if `order` names a file that does not exist, or a file's `<section>` id differs from its name.

## Building

```
cd scripts/slides
npm run deck
```

It prints the three files it wrote. Read anything else it prints:

- `WARNING: these faces did not load` means the browser could not reach Google Fonts and Roboto is not installed, so the slides were measured in another typeface and the editable deck will not match. Rebuild online, or install the fonts.
- `slide <id>: content runs past the bottom margin by up to N px` means the slide holds more than fits. Cut words or split the slide.
- `slide <id>: ... is not drawn in the editable deck` means the slide uses something the PowerPoint conversion does not handle (see Limits). The HTML and images decks still show it.
- `No browser found`: install Edge or Chrome, or set `CHROMIUM_PATH`.
- `EBUSY`, or a refusal to write a `.pptx`, means that deck is open in PowerPoint. Close it and build again.

Then open `rnaseq-flow_slide_deck.html` and page through it.

## Keeping the Slides artifact in step

The Slides artifact and this folder are copies of each other, not linked. After editing in one, bring the other up to date before you build or present.

If you edited in the Slides artifact, ask Claude to copy the deck's files from the artifact into `scripts/slides/deck/`, then run `npm run deck`. If you edited here, ask Claude to publish the changed files to the artifact.

Only `deck.json`, `deck.css` and `slides/*.html` move. The artifact's copy of the design system tokens (`ds/`) stays in the artifact.

## Limits

- PowerPoint wraps text itself. The boxes sit where the browser puts them, but now and then a line breaks one word earlier or later than in the HTML. The images deck has no such difference.
- PowerPoint offers regular and bold here, so text at weight 500 (leads, eyebrows) is regular in the editable deck, and weight 800 (the cover title) is bold.
- The cards that rise into place do so in the browser only.
- What converts: text, boxes with a fill, border or rounded corners, single-side rules, bulleted lists, tables, and inside an `<svg>` only `<line>` and `<circle>`, which is how the route lines and stations are drawn. Anything else (an image, a `<path>`, an `x-` element from the Slides editor) shows in the HTML and images decks, is left out of the editable deck, and is reported.

## The previous deck

The deck this one replaces, `rnaseq-flow.html`, and its builder `build_slides.js` are in `working_docs/archive-decks/` with the PDF and PowerPoint copies made from it. They are kept for reference and are not built.
