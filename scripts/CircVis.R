#!/usr/bin/env Rscript

library(optparse)
library(circlize)
library(tidyverse)
library(ComplexHeatmap)

# ── Argument parsing ──────────────────────────────────────────────────────────

option_list <- list(
  make_option("--reads",                    type = "character", help = "BED file of aligned reads (from bedtools bamtobed)"),
  make_option("--pmd",                      type = "character", help = "PMD scores file (from metaDMG-cpp)"),
  make_option("--depth",                    type = "character", help = "Depth file (from samtools depth)"),
  make_option("--mafs",                     type = "character", help = "Allele frequency file (from ANGSD, .mafs)"),
  make_option("--genes",                    type = "character", default = NULL,  help = "Optional NCBI reference gene table"),
  make_option("--no-genes",                 action = "store_true", default = FALSE, help = "Disable gene track even if --genes is supplied"),
  make_option("--plot-name",                type = "character", default = "plot", help = "Label shown in the centre of the plot [default: plot]"),
  make_option("--out",                      type = "character", help = "Output file path (.png or .pdf)"),
  make_option("--mutation-min-reads",       type = "integer",   default = 10,    help = "Minimum reads supporting a mutation [default: 10]"),
  make_option("--transition-min-frequency", type = "double",    default = 0.5,   help = "Minimum frequency for transitions (C>T, G>A, etc.) [default: 0.5]"),
  make_option("--point-size",               type = "double",    default = 0.5,   help = "Point size in mutation tracks [default: 0.5]"),
  make_option("--segment-length",           type = "integer",   default = 5000,  help = "Axis guide line spacing in bp [default: 5000]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ── Input validation ──────────────────────────────────────────────────────────

required <- c("reads", "pmd", "depth", "mafs", "out")
missing  <- required[sapply(required, function(x) is.null(opt[[x]]))]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste0("--", missing, collapse = ", "))
}

for (f in c(opt$reads, opt$pmd, opt$depth, opt$mafs)) {
  if (!file.exists(f)) stop("File not found: ", f)
}

out_ext <- tolower(tools::file_ext(opt$out))
if (!out_ext %in% c("png", "pdf")) {
  stop("--out must end in .png or .pdf")
}

showgenes <- !isTRUE(opt[["no-genes"]]) && !is.null(opt$genes)
if (showgenes && !file.exists(opt$genes)) stop("Gene file not found: ", opt$genes)

plotname   <- opt[["plot-name"]]
mutatemin  <- opt[["mutation-min-reads"]]
transifreq <- opt[["transition-min-frequency"]]
pointsize  <- opt[["point-size"]]
seglength  <- opt[["segment-length"]]

# ── Load data ─────────────────────────────────────────────────────────────────

reads <- read.table(opt$reads, header = FALSE)
dmg   <- read.table(opt$pmd,   header = FALSE)
depth <- read.delim(opt$depth,  header = FALSE)
mutfreq <- read.delim(opt$mafs, header = TRUE)

if (showgenes) {
  genes <- read.delim(opt$genes, header = TRUE)
}

# ── Merge reads and PMD scores ────────────────────────────────────────────────

reads <- reads %>% rename(chrom = V1, readStart = V2, readEnd = V3, readID = V4)
dmg   <- dmg %>%
  subset(select = c(V1, V3)) %>%
  rename(readID = V1, PMD = V3)

reads <- merge(reads, dmg, by = "readID", all = TRUE)
reads <- reads[, c(2, 3, 4, 1, 7, 5, 6)]
reads <- reads %>% drop_na()

if (nrow(reads) == 0) {
  stop("No reads remain after merging BED and PMD files. Check that read IDs match between the two files.")
}

# ── Reference size and depth ──────────────────────────────────────────────────

refsize  <- max(depth$V2)
depthmax <- max(depth$V3)

# ── Filter mutations ──────────────────────────────────────────────────────────

transitions <- function(df) {
  (df$ref == "C" & df$minor == "T") |
  (df$ref == "T" & df$minor == "C") |
  (df$ref == "A" & df$minor == "G") |
  (df$ref == "G" & df$minor == "A")
}

mutfreq <- mutfreq %>%
  filter(phat > 0 & (!transitions(.) | phat > transifreq)) %>%
  merge(depth, by.x = "position", by.y = "V2", all = FALSE) %>%
  mutate(Mutreads = round(phat * V3)) %>%
  filter(Mutreads >= mutatemin) %>%
  mutate(mutation = paste(ref, minor, sep = " to "))

# ── Gene annotation ───────────────────────────────────────────────────────────

if (showgenes) {
  genes <- genes %>%
    rename(
      start_position = start_position_on_the_genomic_accession,
      end_position   = end_position_on_the_genomic_accession
    ) %>%
    drop_na(start_position)

  genesmid <- (genes$start_position + genes$end_position) / 2

  filter_mutations <- function(g, mf) {
    mf %>% filter(position >= g$start_position & position <= g$end_position)
  }
  mutgenes <- genes %>% split(1:nrow(.)) %>% map_df(~filter_mutations(.x, mutfreq))
}

# ── Color scales ──────────────────────────────────────────────────────────────

col_pmd  <- colorRamp2(
  c(min(reads$PMD), 1.999999, 2, max(reads$PMD)),
  c("darkgrey", "darkgrey", "yellow", "red")
)
col_mut  <- colorRamp2(c(0, 1), c(add_transparency("yellow", 0.3), "forestgreen"))
col_mut2 <- c(
  "A to C" = "chartreuse", "C to A" = "blue",
  "A to T" = "red3",       "T to A" = "purple",
  "G to T" = "gold",       "T to G" = "orangered",
  "G to C" = "cyan",       "C to G" = "magenta",
  "C to T" = "mediumaquamarine", "T to C" = "steelblue",
  "A to G" = "lightsalmon",     "G to A" = "plum"
)

# ── Legends ───────────────────────────────────────────────────────────────────

lgnd_pmd  <- Legend(at = c(min(reads$PMD), 2, max(reads$PMD)),
                    col_fun = col_pmd, title_position = "topleft", title = "PMD")
lgnd_mut  <- Legend(at = c(0, 0.25, 0.5, 0.75, 1),
                    col_fun = col_mut, title_position = "topleft", title = "Mutation frequency")
lgnd_mut2 <- Legend(labels = unique(mutfreq$mutation),
                    legend_gp = gpar(col = col_mut2[unique(mutfreq$mutation)]),
                    type = "points", title_position = "topleft", title = "Mutation type")

# ── Open output device ────────────────────────────────────────────────────────

if (out_ext == "pdf") {
  pdf(file = opt$out, width = 10, height = 10, pointsize = 12)
} else {
  png(file = opt$out, width = 2500, height = 2500, res = 300)
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
      circos.segments(mutgenes$position, 0.2, mutgenes$position, 0.8,
                      col = col_mut2[mutgenes$mutation])
      circos.text(x = genesmid, y = 1.75, labels = genes$description,
                  adj = c(0, 0), cex = 0.3, facing = "clockwise", niceFacing = TRUE)
    }, bg.border = NA)
}

# Track: aligned reads coloured by PMD
circos.trackPlotRegion(ylim = c(0, nrow(reads)), track.height = 0.2,
  panel.fun = function(x, y) {
    circos.segments(axissegment, 0, axissegment, nrow(reads), col = "lightgrey")
    circos.segments(0, 0, 0, nrow(reads))
    for (i in 1:nrow(reads)) {
      y_coord <- nrow(reads) - i
      circos.rect(
        xleft   = reads[i, "readStart"],
        xright  = reads[i, "readEnd"],
        ytop    = pmin(y_coord + nrow(reads) / 40, nrow(reads)),
        ybottom = pmax(y_coord - nrow(reads) / 40, 0),
        col     = col_pmd(reads[i, "PMD"]),
        border  = col_pmd(reads[i, "PMD"])
      )
    }
  })

# Track: depth of coverage
circos.trackPlotRegion(ylim = c(0, depthmax * 1.1), track.height = 0.2,
  panel.fun = function(x, y) {
    circos.segments(axissegment, 0, axissegment, depthmax * 1.1, col = "lightgrey")
    circos.lines(depth$V2, depth$V3, type = "l", straight = TRUE)
    circos.yaxis(at = depthmax, labels.cex = 0.5)
  })

# Track: mutation frequency
circos.trackPlotRegion(ylim = c(0, depthmax * 1.1), track.height = 0.15,
  panel.fun = function(x, y) {
    circos.segments(axissegment, 0, axissegment, depthmax * 1.1, col = "lightgrey")
    circos.points(mutfreq$position, mutfreq$V3, pch = 16, cex = pointsize,
                  col = add_transparency(col_mut(mutfreq$phat), 0.3))
    circos.yaxis(at = depthmax, labels.cex = 0.5)
  })

# Track: mutation type
circos.trackPlotRegion(ylim = c(0, depthmax * 1.1), track.height = 0.10,
  panel.fun = function(x, y) {
    circos.segments(axissegment, 0, axissegment, depthmax * 1.1, col = "lightgrey")
    circos.segments(0, 0, 0, depthmax * 1.1)
    circos.points(mutfreq$position, mutfreq$V3, pch = 16, cex = pointsize,
                  col = add_transparency(col_mut2[mutfreq$mutation], 0.3))
  })

# Legends
draw(packLegend(lgnd_pmd, lgnd_mut),
     x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
draw(lgnd_mut2,
     x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))

# Centre label
text(0, 0, label = plotname, cex = 0.8)

# Track legend (upper left)
if (showgenes) {
  text(-1.35, 1.1, "
     Track 1: Reference genes
     Track 2: Aligned reads with PMD score
     Track 3: Depth of coverage
     Track 4: Mutation frequency
     Track 5: Mutation type", adj = c(0, 0), cex = 0.8)
} else {
  text(-1.35, 1.1, "
     Track 1: Aligned reads with PMD score
     Track 2: Depth of coverage
     Track 3: Mutation frequency
     Track 4: Mutation type", adj = c(0, 0), cex = 0.8)
}

# Stats (upper right)
if (showgenes) {
  text(1.30, 1.2, adj = 1, cex = 0.8, labels = paste(
    "Read count:", nrow(reads),
    "\nMutation count:", nrow(mutfreq),
    "\nMutations in genes:", nrow(mutgenes),
    "\nMinimum mutation depth:", mutatemin,
    "\nMinimum transition frequency:", transifreq
  ))
} else {
  text(1.30, 1.2, adj = 1, cex = 0.8, labels = paste(
    "Read count:", nrow(reads),
    "\nMutation count:", nrow(mutfreq),
    "\nMinimum mutation depth:", mutatemin,
    "\nMinimum transition frequency:", transifreq
  ))
}

circos.clear()
dev.off()

message("Plot saved to: ", opt$out)
