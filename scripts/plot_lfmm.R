#!/usr/bin/env Rscript
library(argparser, quietly=TRUE)

ap = arg_parser("Plot a schnelLFMM result")

ap = add_argument(ap, "assoc_file", help=".tsv file from a schnelLFMM run")
ap = add_argument(ap, "--pval", help="Which phenotype p-value column to use", default="p_ath_sim_h50")
ap = add_argument(ap, "--output", help="Output file")
ap = add_argument(ap, "--highlight-snps", help="highlighted snp file")

args = parse_args(ap)

if (is.na(args$output)) {
    args$output = paste0(args$assoc_file, ".plot.png")
}

str(args)

library(readr)
library(dplyr)
library(ggplot2)

highlighted = c()
if (!is.na(args$highlight_snps)) {
    highlighted = read_tsv(args$highlight_snps) %>%
        pull(1)
    print(paste("Highlighting", length(highlighted), "SNPs"))
}

gwas = read_tsv(args$assoc_file) %>%
    rename(pval={!!args$pval}, rs=snp_id, ps=pos) %>%
    mutate(highlight=rs %in% highlighted) %>%
    arrange(highlight, chr, ps) %>%
    glimpse()
p = ggplot(gwas, aes(x=ps, y=-log10(pval))) +
    geom_point(aes(colour=highlight)) +
    facet_grid(~chr, scales="free_x", space="free_x") +
    theme_bw()
ggsave(args$output, plot=p, dpi=600, width=8, height=4)