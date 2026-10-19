# Worked examples

Three examples: one in yeast and two in human. Each answers a different
question, so pick by what you want to find out.

| Example | Organism | Answers | Cost |
|---|---|---|---|
| [Four aligners on yeast](yeast/README.md) | *S. cerevisiae* | Does the pipeline run correctly on my machine, on every aligner? | ~1 GB subsampled, under an hour per route on a laptop |
| [Expression walkthrough](../rnaseq-flow_Example_Analysis.md) | Human | What does a realistic differential expression analysis look like end to end? | ~4 GB subsampled, a couple of hours |
| [Splicing and isoforms](human/README.md) | Human | Do `--dtu`, `--diffsplice` and `--isoform_switch` work, and do they find what the literature says? | ~10 GB at full depth, a few hours plus index builds |

They are listed in the order most people want them.

## Which to run first

**[Four aligners on yeast](yeast/README.md)** if you have just installed the
pipeline. It is the only example that exercises STAR, HISAT2, Salmon and
Kallisto on the same data and shows that all four agree, and it is sized for a
laptop. Nothing else proves an installation as quickly.

**[Expression walkthrough](../rnaseq-flow_Example_Analysis.md)** if you want to
see the shape of a real analysis. It takes GSE151251, six samples of human
hepatic stellate cells treated with TGF-beta, through HISAT2 to a fibrosis
expression signature you can check your output against. It lives at the
repository root rather than in this folder, because it is the document the
guide points newcomers at.

**[Splicing and isoforms](human/README.md)** if you care about
`--dtu`, `--diffsplice` or `--isoform_switch`. These cannot be exercised on
yeast at all: Ensembl R64-1-1 has 7127 genes and 7127 transcripts, so no gene
has the two or more isoforms those tests need, and they correctly produce empty
results there.

## Why the two human examples are separate

Both use GSE151251, and they are still two documents rather than one, because
they prepare the data differently. The expression example subsamples the reads;
splicing and transcript usage need full depth, because junction and isoform
evidence is a small fraction of gene-level counts and subsampling removes the
signal before the test sees it. One document cannot tell you to do both.
