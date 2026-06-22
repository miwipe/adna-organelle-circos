#!/usr/bin/env Rscript

library(optparse)
library(circlize)
library(tidyverse)
library(ComplexHeatmap)

# ── Argument parsing ──────────────────────────────────────────────────────────

option_list <- list(
  make_option("--depth1",                   type = "character", help = "Depth file for alignment 1 (from samtools depth)"),
  make_option("--mafs1",                    type = "character", help = "Allele frequency file for alignment 1 (from ANGSD, .mafs)"),
  make_option("--depth2",                   type = "character", help = "Depth file for alignment 2 (from samtools depth)"),
  make_option("--mafs2",                    type = "character", help = "Allele frequency file for alignment 2 (from ANGSD, .mafs)"),
  make_option("--name1",                    type = "character", default = "Alignment 1", help = "Label for alignment 1 [default: 'Alignment 1']"),
  make_option("--name2",                    type = "character", default = "Alignment 2", help = "Label for alignment 2 [default: 'Alignment 2']"),
  make_option("--genes",                    type = "character", default = NULL,  help = "Optional NCBI reference gene table"),
  make_option("--no-genes",                 action = "store_true", default = FALSE, help = "Disable gene track even if --genes is supplied"),
  make_option("--plot-name",                type = "character", default = "plot", help = "Label shown in the centre of the plot [default: plot]"),
  make_option("--out",                      type = "character", help = "Output file path (.png or .pdf)"),
  make_option("--mutation-min-reads1",      type = "integer",   default = 10,   help = "Minimum reads supporting a mutation in alignment 1 [default: 10]"),
  make_option("--mutation-min-reads2",      type = "integer",   default = 10,   help = "Minimum reads supporting a mutation in alignment 2 [default: 10]"),
  make_option("--transition-min-frequency", type = "double",    default = 0.5,  help = "Minimum frequency for transitions (C>T, G>A, etc.) [default: 0.5]"),
  make_option("--segment-length",           type = "integer",   default = 5000, help = "Axis guide line spacing in bp [default: 5000]"),
  make_option("--point-size",               type = "double",    default = 0.5,  help = "Point size in mutation tracks [default: 0.5]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ── Input validation ──────────────────────────────────────────────────────────

required <- c("depth1", "mafs1", "depth2", "mafs2", "out")
missing  <- required[sapply(required, function(x) is.null(opt[[x]]))]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste0("--", missing, collapse = ", "))
}

for (f in c(opt$depth1, opt$mafs1, opt$depth2, opt$mafs2)) {
  if (!file.exists(f)) stop("File not found: ", f)
}

out_ext <- tolower(tools::file_ext(opt$out))
if (!out_ext %in% c("png", "pdf")) {
  stop("--out must end in .png or .pdf")
}

showgenes <- !isTRUE(opt[["no-genes"]]) && !is.null(opt$genes)
if (showgenes && !file.exists(opt$genes)) stop("Gene file not found: ", opt$genes)

plotname   <- opt[["plot-name"]]
alignment1 <- opt$name1
alignment2 <- opt$name2
mutatemin  <- opt[["mutation-min-reads1"]]
mutatemin2 <- opt[["mutation-min-reads2"]]
transifreq <- opt[["transition-min-frequency"]]
seglength  <- opt[["segment-length"]]
pointsize  <- opt[["point-size"]]

# ── Load data ─────────────────────────────────────────────────────────────────

depth   <- read.delim(opt$depth1, header = FALSE)
mutfreq <- read.delim(opt$mafs1,  header = TRUE)
depth2  <- read.delim(opt$depth2, header = FALSE)
mutfreq2 <- read.delim(opt$mafs2, header = TRUE)

if (showgenes) {
  genes <- read.delim(opt$genes, header = TRUE)
}

# ── Reference size and depth ──────────────────────────────────────────────────

refsize  <- max(depth$V2)
depthmax <- max(depth$V3)
depthmax2 <- max(depth2$V3)

# ── Filter mutations ──────────────────────────────────────────────────────────

transitions <- function(df) {
  (df$ref == "C" & df$minor == "T") |
  (df$ref == "T" & df$minor == "C") |
  (df$ref == "A" & df$minor == "G") |
  (df$ref == "G" & df$minor == "A")
}

filter_mutfreq <- function(mf, dep, minreads) {
  mf %>%
    filter(phat > 0 & (!transitions(.) | phat > transifreq)) %>%
    merge(dep, by.x = "position", by.y = "V2", all = FALSE) %>%
    mutate(Mutreads = round(phat * V3)) %>%
    filter(Mutreads >= minreads) %>%
    mutate(mutation = paste(ref, minor, sep = " to "))
}

mutfreq  <- filter_mutfreq(mutfreq,  depth,  mutatemin)
mutfreq2 <- filter_mutfreq(mutfreq2, depth2, mutatemin2)

# ── Compare mutations ─────────────────────────────────────────────────────────

uniqmut  <- anti_join(mutfreq,  mutfreq2, by = c("position", "mutation"))
uniqmut2 <- anti_join(mutfreq2, mutfreq,  by = c("position", "mutation"))
mutsame  <- inner_join(mutfreq, mutfreq2, by = c("position", "mutation", "minor", "major", "ref"))

# ── Gene annotation ───────────────────────────────────────────────────────────

if (showgenes) {
  genes <- genes %>%
    rename(
      start_position = start_position_on_the_genomic_accession,
      end_position   = end_position_on_the_genomic_accession
    ) %>%
    drop_na(start_position)

  genesmid <- (genes$start_position + genes$end_position) / 2

  filter_in_genes <- function(g, mf) {
    mf %>% filter(position >= g$start_position & position <= g$end_position)
  }
  mutgenes2 <- genes %>% split(1:nrow(.)) %>% map_df(~filter_in_genes(.x, mutsame))
}

# ── Color scale ───────────────────────────────────────────────────────────────

col_mut2 <- c(
  "A to C" = "chartreuse", "C to A" = "blue",
  "A to T" = "red3",       "T to A" = "purple",
  "G to T" = "gold",       "T to G" = "orangered",
  "G to C" = "cyan",       "C to G" = "magenta",
  "C to T" = "mediumaquamarine", "T to C" = "steelblue",
  "A to G" = "lightsalmon",     "G to A" = "plum"
)

lgnd_mut2 <- Legend(
  labels      = unique(mutfreq$mutation),
  legend_gp   = gpar(col = col_mut2[unique(mutfreq$mutation)]),
  type        = "points",
  title_position = "topleft",
  title       = "Mutation type"
)

# ── Open output device ────────────────────────────────────────────────────────

if (out_ext == "pdf") {
  pdf(file = opt$out, width = 10, height = 10, pointsize = 12)
} else {
  png(file = opt$out, width = 2800, height = 2800, res = 300)
}

# ── Plot ──────────────────────────────────────────────────────────────────────

axissegment <- seq(0, refsize, by = seglength)
df <- data.frame(name = " ", start = 0, end = refsize)

circos.par(
  "start.degree" = 90, "gap.degree" = 0,
  "track.margin" = c(0, 0),
  "cell.padding" = c(0.002, 0.002, 0.002, 0.002),
  "canvas.xlim" = c(-1.2, 1.2), "canvas.ylim" = c(-1.2, 1.2),
  "points.overflow.warning" = FALSE
)

circos.genomicInitialize(df, axis.labels.cex = 0.6 * par("cex"), major.by = seglength)

# Track 1: reference genes (optional)
if (showgenes) {
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.075,
    panel.fun = function(x, y) {
      circos.segments(axissegment, 0, axissegment, 1, col = "lightgrey")
      circos.rect(genes$start_position, 0.2, genes$end_position, 0.8)
      circos.segments(mutgenes2$position, 0.2, mutgenes2$position, 0.8,
                      col = col_mut2[mutgenes2$mutation])
      circos.text(x = genesmid, y = 1.75, labels = genes$description,
                  adj = c(0, 0), cex = 0.3, facing = "clockwise", niceFacing = TRUE)
    }, bg.border = NA)
}

# Track: shared mutations
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.15,
  panel.fun = function(x, y) {
    circos.segments(axissegment, 0, axissegment, 1, col = "lightgrey")
    circos.segments(0, 0, 0, 1)
    circos.segments(mutsame$position, 0, mutsame$position, 1,
                    col = add_transparency(col_mut2[mutsame$mutation], 0.3))
  })

# Track: mutations unique to alignment 1
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.10,
  panel.fun = function(x, y) {
    circos.segments(axissegment, 0, axissegment, 1, col = "lightgrey")
    circos.segments(0, 0, 0, 1)
    circos.segments(uniqmut$position, 0, uniqmut$position, 1,
                    col = add_transparency(col_mut2[uniqmut$mutation], 0.3))
  })

# Track: mutations unique to alignment 2
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.10,
  panel.fun = function(x, y) {
    circos.segments(axissegment, 0, axissegment, 1, col = "lightgrey")
    circos.segments(0, 0, 0, 1)
    circos.segments(uniqmut2$position, 0, uniqmut2$position, 1,
                    col = add_transparency(col_mut2[uniqmut2$mutation], 0.3))
  })

# Legend
draw(lgnd_mut2,
     x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))

# Centre label
text(0, 0, labels = plotname, cex = 1)

# Track legend (upper left)
if (showgenes) {
  text(-1.35, 1.1, adj = c(0, 0), cex = 0.8, labels = paste(
    "\n     Track 1: Reference genes",
    "\n     Track 2: Shared mutations",
    paste("\n     Track 3: Mutations unique to", alignment1),
    paste("\n     Track 4: Mutations unique to", alignment2)
  ))
} else {
  text(-1.35, 1.1, adj = c(0, 0), cex = 0.8, labels = paste(
    "\n     Track 1: Shared mutations",
    paste("\n     Track 2: Mutations unique to", alignment1),
    paste("\n     Track 3: Mutations unique to", alignment2)
  ))
}

# Stats (upper right)
if (showgenes) {
  text(1.30, 1.15, adj = 1, cex = 0.8, labels = paste(
    "Shared mutations:", nrow(mutsame),
    "\nShared mutations in genes:", nrow(mutgenes2),
    paste("\nUnique", alignment1, "mutations:"), nrow(uniqmut),
    paste("\nUnique", alignment2, "mutations:"), nrow(uniqmut2),
    paste("\n", alignment1, "minimum mutation depth:"), mutatemin,
    paste("\n", alignment2, "minimum mutation depth:"), mutatemin2,
    "\nMinimum transition frequency:", transifreq
  ))
} else {
  text(1.30, 1.15, adj = 1, cex = 0.8, labels = paste(
    "Shared mutations:", nrow(mutsame),
    paste("\nUnique", alignment1, "mutations:"), nrow(uniqmut),
    paste("\nUnique", alignment2, "mutations:"), nrow(uniqmut2),
    paste("\n", alignment1, "minimum mutation depth:"), mutatemin,
    paste("\n", alignment2, "minimum mutation depth:"), mutatemin2,
    "\nMinimum transition frequency:", transifreq
  ))
}

circos.clear()
dev.off()

message("Plot saved to: ", opt$out)
