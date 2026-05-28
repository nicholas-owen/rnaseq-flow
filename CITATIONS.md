# Citations

Every tool used by **rnaseq-flow**, with the publication to cite. Bibliographic
details were verified against PubMed. A [ready-to-paste methods
paragraph](#methods-paragraph) for manuscripts follows the reference list.

Cite only the tools that actually ran for your analysis — which tools run
depends on `--aligner` and the optional flags you set. The exact version of
every tool used in a run is recorded in `pipeline_info/software_versions.yml`.

---

## Workflow engine

- **Nextflow**
  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017;35(4):316-319. doi:[10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

## Read QC & trimming

- **FastQC** — raw-read quality control.
  > Andrews S. FastQC: A Quality Control Tool for High Throughput Sequence Data. Babraham Bioinformatics; 2010. <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>

- **fastp** — adapter and quality trimming.
  > Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34(17):i884-i890. doi:[10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)

## Alignment & quantification

- **STAR** — splice-aware genome alignment.
  > Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21. doi:[10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)

- **HISAT2** — graph-based genome alignment.
  > Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol. 2019;37(8):907-915. doi:[10.1038/s41587-019-0201-4](https://doi.org/10.1038/s41587-019-0201-4)

- **Salmon** — bias-aware transcript quantification.
  > Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods. 2017;14(4):417-419. doi:[10.1038/nmeth.4197](https://doi.org/10.1038/nmeth.4197)

- **kallisto** — pseudoalignment-based transcript quantification.
  > Bray NL, Pimentel H, Melsted P, Pachter L. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol. 2016;34(5):525-527. doi:[10.1038/nbt.3519](https://doi.org/10.1038/nbt.3519)

- **SAMtools** — BAM sorting and indexing.
  > Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008. doi:[10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)

## Post-alignment QC & coverage

- **RSeQC** — strandedness inference, read distribution, gene-body coverage.
  > Wang L, Wang S, Li W. RSeQC: quality control of RNA-seq experiments. Bioinformatics. 2012;28(16):2184-2185. doi:[10.1093/bioinformatics/bts356](https://doi.org/10.1093/bioinformatics/bts356)

- **deepTools** — CPM-normalised BigWig coverage tracks (`bamCoverage`).
  > Ramírez F, Ryan DP, Grüning B, Bhardwaj V, Kilpert F, Richter AS, Heyne S, Dündar F, Manke T. deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Res. 2016;44(W1):W160-W165. doi:[10.1093/nar/gkw257](https://doi.org/10.1093/nar/gkw257)

## Gene quantification & differential expression

- **featureCounts** (Subread) — gene-level and exon-level read counting.
  > Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014;30(7):923-930. doi:[10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656)

- **tximport** — summarises Salmon/kallisto transcript counts to gene level.
  > Soneson C, Love MI, Robinson MD. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research. 2015;4:1521. doi:[10.12688/f1000research.7563.2](https://doi.org/10.12688/f1000research.7563.2)

- **DESeq2** — differential gene expression.
  > Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550. doi:[10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

- **apeglm** — log2 fold-change shrinkage for the DESeq2 contrasts (`lfcShrink`).
  > Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. 2019;35(12):2084-2092. doi:[10.1093/bioinformatics/bty895](https://doi.org/10.1093/bioinformatics/bty895)

- **edgeR** — differential gene expression and, via `diffSpliceDGE`, differential exon/transcript splicing.
  > Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010;26(1):139-140. doi:[10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616)

## Alternative splicing, isoforms & fusions

- **rMATS** — event-based differential alternative splicing.
  > Shen S, Park JW, Lu ZX, Lin L, Henry MD, Wu YN, Zhou Q, Xing Y. rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. Proc Natl Acad Sci USA. 2014;111(51):E5593-E5601. doi:[10.1073/pnas.1419161111](https://doi.org/10.1073/pnas.1419161111)

- **STAR-Fusion** — gene-fusion detection from STAR chimeric junctions.
  > Haas BJ, Dobin A, Li B, Stransky N, Pochet N, Regev A. Accuracy assessment of fusion transcript detection via read-mapping and de novo fusion transcript assembly-based methods. Genome Biol. 2019;20(1):213. doi:[10.1186/s13059-019-1842-9](https://doi.org/10.1186/s13059-019-1842-9)

- **IsoformSwitchAnalyzeR** — transcript isoform-switch analysis.
  > Vitting-Seerup K, Sandelin A. IsoformSwitchAnalyzeR: analysis of changes in genome-wide patterns of alternative splicing and its functional consequences. Bioinformatics. 2019;35(21):4469-4471. doi:[10.1093/bioinformatics/btz247](https://doi.org/10.1093/bioinformatics/btz247)

- **DEXSeq** — differential transcript usage (the `--dtu` test).
  > Anders S, Reyes A, Huber W. Detecting differential usage of exons from RNA-seq data. Genome Res. 2012;22(10):2008-2017. doi:[10.1101/gr.133744.111](https://doi.org/10.1101/gr.133744.111)

## Functional enrichment

- **fgsea** — fast gene-set enrichment analysis (GSEA).
  > Korotkevich G, Sukhov V, Budin N, Shpak B, Artyomov MN, Sergushichev A. Fast gene set enrichment analysis. bioRxiv. 2021. doi:[10.1101/060012](https://doi.org/10.1101/060012)

  fgsea implements the GSEA method:
  > Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, Mesirov JP. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci USA. 2005;102(43):15545-15550. doi:[10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)

- **g:Profiler** (gprofiler2) — GO and pathway over-representation.
  > Raudvere U, Kolberg L, Kuzmin I, Arak T, Adler P, Peterson H, Vilo J. g:Profiler: a web server for functional enrichment analysis and conversions of gene lists (2019 update). Nucleic Acids Res. 2019;47(W1):W191-W198. doi:[10.1093/nar/gkz369](https://doi.org/10.1093/nar/gkz369)

## Reporting

- **MultiQC** — aggregated QC report.
  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-3048. doi:[10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)

- **Quarto** — rendered HTML analysis report.
  > Allaire JJ, Teague C, Scheidegger C, Xie Y, Dervieux C. Quarto. doi:[10.5281/zenodo.5960048](https://doi.org/10.5281/zenodo.5960048) — <https://quarto.org>

- **plotly** (R package) — interactive volcano plots in the analysis report.
  > Sievert C. Interactive Web-Based Data Visualization with R, plotly, and shiny. Chapman and Hall/CRC; 2020. ISBN 9781138331457 — <https://plotly-r.com>

- **DT** (R package) — searchable result tables in the analysis report.
  > Xie Y, Cheng J, Tan X. DT: A Wrapper of the JavaScript Library 'DataTables'. — <https://CRAN.R-project.org/package=DT>

## Reference data

- **Ensembl** — genome FASTA and GTF annotation (the `--download_refs` workflow). Cite the Ensembl release you used; the most recent reference is:
  > Dyer SC, Austine-Orimoloye O, Azov AG, et al. Ensembl 2025. Nucleic Acids Res. 2025;53(D1):D948-D957. doi:[10.1093/nar/gkae1071](https://doi.org/10.1093/nar/gkae1071)

- **MSigDB** — Hallmark / C2 / C5 gene sets for GSEA (the `--download_gmt` option).
  > Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015;1(6):417-425. doi:[10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)

## Software packaging & containers

- **Bioconda** — package distribution.
  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018;15(7):475-476. doi:[10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7)

- **BioContainers** — containerised tool images.
  > da Veiga Leprevost F, Grüning BA, Alves Aflitos S, et al. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017;33(16):2580-2582. doi:[10.1093/bioinformatics/btx192](https://doi.org/10.1093/bioinformatics/btx192)

- **Docker**
  > Merkel D. Docker: lightweight Linux containers for consistent development and deployment. Linux Journal. 2014;2014(239):2.

- **Singularity / Apptainer**
  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS ONE. 2017;12(5):e0177459. doi:[10.1371/journal.pone.0177459](https://doi.org/10.1371/journal.pone.0177459)

---

## Methods paragraph

A template for a manuscript Methods section. Edit it to match the run you
performed — delete the stages you did not use, and fill in the bracketed
tool **versions** from `pipeline_info/software_versions.yml`.

> RNA-seq data were processed with the rnaseq-flow pipeline (v1.1.0), implemented
> in Nextflow (Di Tommaso et al., 2017); each step ran in a Bioconda-based
> (Grüning et al., 2018) container (BioContainers; da Veiga Leprevost et al.,
> 2017). Raw reads were assessed with FastQC v[…] (Andrews, 2010) and
> adapter/quality trimmed with fastp v[…] (Chen et al., 2018). Trimmed reads
> were aligned to the [genome/transcriptome] with [STAR v[…] (Dobin et al.,
> 2013) | HISAT2 v[…] (Kim et al., 2019) | Salmon v[…] (Patro et al., 2017) |
> kallisto v[…] (Bray et al., 2016)]. [For the genome aligners, alignments were
> sorted and indexed with SAMtools v[…] (Danecek et al., 2021); library
> strandedness, read distribution and gene-body coverage were assessed with
> RSeQC v[…] (Wang et al., 2012), and CPM-normalised coverage tracks were
> generated with deepTools v[…] (Ramírez et al., 2016).] Gene-level counts were
> obtained with featureCounts v[…] (Liao et al., 2014) [or summarised from
> transcript-level estimates with tximport v[…] (Soneson et al., 2015)].
> Differential expression was tested with DESeq2 v[…] (Love et al., 2014),
> with log2 fold-change shrinkage via apeglm (Zhu et al., 2019), and edgeR v[…]
> (Robinson et al., 2010). [Differential alternative splicing was
> assessed with rMATS v[…] (Shen et al., 2014) and edgeR diffSpliceDGE (Robinson
> et al., 2010); differential transcript usage was tested with DEXSeq v[…]
> (Anders et al., 2012); isoform switches were identified with
> IsoformSwitchAnalyzeR v[…] (Vitting-Seerup & Sandelin, 2019); gene fusions
> were called with STAR-Fusion v[…] (Haas et al., 2019).] Functional enrichment
> used gene-set enrichment analysis via fgsea v[…] (Korotkevich et al., 2021;
> Subramanian et al., 2005) and over-representation analysis with g:Profiler
> (Raudvere et al., 2019). Quality-control metrics across all steps were
> aggregated with MultiQC v[…] (Ewels et al., 2016). Reference genome and
> annotation were obtained from Ensembl [release …] (Dyer et al., 2025) [and
> gene sets from MSigDB (Liberzon et al., 2015)].
