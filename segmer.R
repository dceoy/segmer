#!/usr/bin/env Rscript

'Segment Detector for differentially methylated regions in methylome data

Usage:
  segmer bed [-v] [--platform=<str>] [--unfilter] [--out=<dir>]
  segmer dmr [-v] [--seed=<int>] [--sd-cutoff=<int>] [--ar=<ratio>]
             [--out=<dir>] <site_csv> <bv_csv>
  segmer cluster [-v] [--k=<int>] [--dist=<str>] [--hclust=<str>]
                 [--ar=<ratio>] [--out=<dir>] <dmrbv_csv>
  segmer plot [-v] [--ar=<ratio>] [--out=<dir>] <site_csv> <bv_csv> <seg_csv>
              <dmrbv_csv>
  segmer dmp [-v] [--seed=<int>] [--sd-cutoff=<dbl>] [--out=<dir>] <site_csv>
             <bv_csv>
  segmer --session
  segmer --version
  segmer -h|--help

Commands:
  bed               Download annotation data and write a target site CSV file
  segment           Segment target sites (Sites including NA are ignored.)
  cluster           Execute clustering and draw the heatmap
  plot              Visualize methylation data
  dmp               Determine differentially methylated positions

Options:
  -v                Run with debug logging
  --platform=<str>  Specify a methylation assay platform [default: EPIC]
                      choice: EPIC, hm450, hm27
  --unfilter        Skip recommended probe filtering
  --out=<dir>       Set an output directory [default: .]
  --seed=<int>      Set a random seed
  --k=<int>         Specify the number of clusters [default: 3]
  --dist=<str>      Specify the method of stats::dist [default: euclidean]
  --hclust=<str>    Specify the method of stats::hclust [default: ward.D2]
  --ar=<ratio>      Specify the aspect ratio of figures [default: 13:8]
  --sd-cutoff=<dbl> Specify the SD cutoff [default: 0.1]
  --session         Print session information and exit ({devtools} is required)
  --version         Print version and exit
  -h, --help        Print help and exit

Arguments:
  <site_csv>        Path to a target site CSV or BED file
                    (created by `segmer bed`)
  <bv_csv>          Path to a methylation CSV file
                      the 1st column: probe names
                      the other columns: beta-values
  <dmrbv_csv>       Path to a segmental beta-value CSV file of DMRs
                    (`*.seg.bv.mean.dmr.csv` created by `segmer segment`)
  <seg_csv>         Path to a segment CSV file
                    (`*.seg.csv` created by `segmer segment`)
' -> doc

command_version <- 'v0.0.9'


### controler

fetch_script_root <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  fa <- ca[grepl('^--file=', ca)]
  if (length(fa) == 1) {
    f <- sub('--file=', '', fa)
    l <- Sys.readlink(f)
    if (is.na(l)) {
      script <- normalizePath(f)
    } else if (startsWith(l, '/')) {
      script <- normalizePath(l)
    } else {
      script <- normalizePath(file.path(dirname(f), l))
    }
    return(dirname(script))
  } else {
    return(normalizePath(getwd()))
  }
}

main <- function(opts, root_dir = fetch_script_root()) {
  options(warn = 1)
  if (opts[['-v']]) {
    print(opts)
  }

  if (opts[['--session']]) {
    library('devtools', quietly = TRUE)
    print(devtools::session_info(pkgs = c('changepoint', 'GenomicRanges',
                                          'gplots', 'ggpubr', 'RColorBrewer',
                                          'tidyverse')))
  } else {
    if (opts[['bed']]) {
      add_pkgs <- 'GenomicRanges'
    } else if (opts[['dmr']]) {
      add_pkgs <- c('changepoint', 'ggpubr', 'parallel')
    } else if (opts[['cluster']]) {
      add_pkgs <- c('gplots', 'RColorBrewer')
    } else if (opts[['plot']]) {
      add_pkgs <- c('ggpubr', 'parallel')
    } else {
      add_pkgs <- NULL
    }
    load_packages(pkgs = c('tidyverse', add_pkgs))
    make_dir(path = opts[['--out']])
    dst_dir <- normalizePath(opts[['--out']])
    if (opts[['cluster']] | opts[['plot']]) {
      aspect_ratio <- as.integer(str_split(opts[['--ar']], pattern = ':',
                                           n = 2, simplify = TRUE))
    } else {
      aspect_ratio <- NULL
    }
    if (! is.null(opts[['--seed']])) {
      message('>>> Set a random seed')
      set.seed(opts[['--seed']])
      message(opts[['--seed']])
    }

    if (opts[['bed']]) {
      prepare_site_csv(dst_dir = dst_dir,
                       platform = opts[['--platform']],
                       unfilter = opts[['--unfilter']])
    } else if (opts[['dmr']]) {
      segment_sites(site_csv = normalizePath(opts[['<site_csv>']]),
                    bv_csv = normalizePath(opts[['<bv_csv>']]),
                    dst_dir = dst_dir,
                    sd_cutoff = as.numeric(opts[['--sd-cutoff']]))
    } else if (opts[['cluster']]) {
      cluster_segments(dmrbv_csv = normalizePath(opts[['<dmrbv_csv>']]),
                       dst_dir = dst_dir, k = opts[['--k']],
                       dist_method = opts[['--dist']],
                       hclust_method = opts[['--hclust']],
                       width = aspect_ratio[1], height = aspect_ratio[2])
    } else if (opts[['plot']]) {
      site_csv <- normalizePath(opts[['<site_csv>']])
      bv_csv <- normalizePath(opts[['<bv_csv>']])
      seg_csv <- normalizePath(opts[['<seg_csv>']])
      dmrbv_csv <- normalizePath(opts[['<dmrbv_csv>']])
      visualize_bv(site_csv = site_csv, bv_csv = bv_csv, dst_dir = dst_dir,
                   width = aspect_ratio[1], height = aspect_ratio[2])
      visualize_segments(seg_csv = seg_csv, site_csv = site_csv,
                         bv_csv = bv_csv, dmrbv_csv = dmrbv_csv,
                         dst_dir = dst_dir, width = aspect_ratio[1],
                         height = aspect_ratio[2])
    } else if (opts[['dmp']]) {
      write_dmp_bv(bv_csv = normalizePath(opts[['<bv_csv>']]),
                   site_csv = normalizePath(opts[['<site_csv>']]),
                   dst_dir = dst_dir,
                   sd_cutoff = as.numeric(opts[['--sd-cutoff']]))
    }
  }
}


### segmer bed

prepare_site_csv <- function(dst_dir, platform = 'EPIC', unfilter = FALSE,
                             hg_ver = 'hg38') {
  ann_df_list <- download_annotation_data(dst_dir = dst_dir)
  if (unfilter) {
    df_ann <- ann_df_list$gencode
    message('probes:\t', nrow(df_ann))
    bed_csv_name <- str_c(platform, hg_ver, 'bed.csv', sep = '.')
  } else {
    message('>>> Filter probes')
    walk(c('chrX', 'chrY', 'MASK_general', 'MASK_snp5_common'),
         function(i) message('- ', i))
    masked_ids <- filter(ann_df_list$mask,
                         seqnames %in% c('chrX', 'chrY')
                         | MASK_general | MASK_snp5_common)$name
    df_ann <- filter(ann_df_list$gencode,
                     ! seqnames %in% c('chrX', 'chrY'),
                     ! name %in% masked_ids)
    message('passing probes:\t', nrow(df_ann), ' / ',
            nrow(ann_df_list$gencode))
    bed_csv_name <- str_c(platform, hg_ver, 'f.bed.csv', sep = '.')
  }
  site_csv <- file.path(dst_dir, bed_csv_name)
  message('>>> Write a sorted target site BED:\t', site_csv)
  write_csv(select(mutate(dplyr::rename(df_ann,
                                        chrom = seqnames, chromEnd = end),
                          chromStart = start - 1),
                   chrom, chromStart, chromEnd, name, genesUniq, CGI,
                   CGIposition),
            path = site_csv)
}

download_annotation_data <- function(dst_dir = '.', platform = 'EPIC',
                                     hg_ver = 'hg38', gencode_ver = 'v22') {
  url_body <- file.path('http://zwdzwd.io/InfiniumAnnotation/current',
                        platform)
  return(lapply(c(mask = str_c(platform, hg_ver, 'manifest.rds', sep = '.'),
                  gencode = str_c(platform, hg_ver, 'manifest.gencode',
                                  gencode_ver, 'rds', sep = '.')),
                function(n) {
                  dst <- file.path(dst_dir, n)
                  download_file(src = file.path(url_body, n), dst = dst)
                  csv_gz <- sub('.rds$', '.csv.gz', dst)
                  message('>>> Convert RDS to CSV:\t', csv_gz)
                  d <- granges2tibble(readRDS(dst), var = 'name')
                  write_csv(d, path = csv_gz)
                  return(d)
                }))
}


### segmer segment

segment_sites <- function(site_csv, bv_csv, dst_dir, n_cpu = detectCores(),
                          cpt_method = 'PELT', mtc_method = 'bonferroni',
                          sd_cutoff = 0.1) {
  df_site <- read_csv_quietly(site_csv)
  df_bv <- filter(read_bv_csv(path = bv_csv), name %in% df_site$name)
  cl <- makeCluster(n_cpu)
  df_p <- mutate(arrange(inner_join(df_site,
                                    append_shapiro_wilk_p(cl, df_bv),
                                    by = 'name'),
                         chrom, chromStart),
                 mlogp = -log10(pval),
                 sid = row_number())
  stopCluster(cl)

  message('>>> Detect changepoints by ', cpt_method, ' method')
  cpt <- changepoint::cpt.mean(df_p$mlogp, method = cpt_method)
  summary(cpt)
  cp_sids <- c(1, changepoint::cpts(cpt))
  message('>>> Segment target sites')
  df_sid <- fill(select(left_join(dplyr::rename(df_p, cp_sid = sid),
                                  tibble(sid = cp_sids, cp_sid = cp_sids),
                                  by = 'cp_sid'),
                        -cp_sid),
                 sid)
  df_seg <- left_join(df_sid,
                      summarize(group_by(df_sid, chrom, sid),
                                segment = str_c(unique(chrom), ':',
                                                min(chromStart) + 1,
                                                '-', max(chromEnd)),
                                .groups = 'drop'),
                      by = c('chrom', 'sid'))
  print(summary(summarize(group_by(df_seg, segment),
                          sites_per_segment = n(),
                          .groups = 'drop')))
  message(nrow(df_seg), ' sites => ', length(unique(df_seg$segment)),
          ' segments')

  seg_csv <- file.path(dst_dir,
                       str_c(sub('.csv(|.gz)$', '', basename(bv_csv)),
                             sub('(|.bed).csv(|.gz)$', '', basename(site_csv)),
                             'seg.csv', sep = '.'))
  message('>>> Write a segment CSV:\t', seg_csv)
  write_csv(df_seg, path = seg_csv)

  message('>>> Calculate segmental median values')
  df_segbv <- summarize_all(group_by(select(left_join(select(df_seg,
                                                             name, segment),
                                                      df_bv, by = 'name'),
                                            -name),
                                     segment),
                            median,
                            .groups = 'drop')
  out_prefix <- file.path(dst_dir,
                          str_c(sub('.csv(|.gz)$', '', basename(seg_csv)),
                                'bv.median', sep = '.'))
  segbv_csv <- str_c(out_prefix, '.all.csv')
  message('>>> Write a segmental median CSV:\t', segbv_csv)
  write_csv(df_segbv, path = segbv_csv)

  df_dmrbv <- filter_bv_by_sd(df_segbv,
                              sd_cutoff = sd_cutoff)
  dmrbv_csv <- str_c(out_prefix, '.dmr.csv')
  message('>>> Write a DMR segmental median CSV:\t', dmrbv_csv)
  write_csv(df_dmrbv, path = dmrbv_csv)
  message('DMR segments:\t', nrow(df_dmrbv), ' / ', nrow(df_segbv))
}

append_shapiro_wilk_p <- function(cl, df_bv) {
   message('>>> Perform Shapiro-Wilk normality tests')
  return(mutate(select(df_bv, 1),
                pval = parApply(cl, select(df_bv, -1), 1,
                                function(v) return(shapiro.test(v)$p.value))))
}

read_site_csv <- function(path) {
  if (any(endsWith(path, c('.bed', '.bed.gz', '.bed.bz2')))) {
    return(read_bed(path))
  } else {
    return(read_csv_quietly(path))
  }
}

read_bv_csv <- function(path, drop_na = TRUE) {
  df_bv <- read_csv_quietly(path)
  if (drop_na) df_bv <- drop_na(df_bv)
  stopifnot(nrow(df_bv) > 1, ncol(df_bv) > 4)
  return(mutate(dplyr::rename(df_bv, name = 1), name = as.character(name)))
}

bv2var <- function(df_bv) {
  return(mutate(select(df_bv, 1), variance = apply(select(df_bv, -1), 1, var)))
}

filter_bv_by_sd <- function(df_bv, sd_cutoff = 0.1) {
  return(filter(df_bv,
                bv2var(df_bv)$variance > (sd_cutoff ^ 2)))
}

determine_bv_boundary <- function(df_bv, ...) {
  message('>>> Determine the boundary by K-means')
  bv <- df2num(df_bv)
  df_b <- summarize(group_by(tibble(bv = bv,
                                    label = kmeans(bv,
                                                   centers = 2,
                                                   ...)$cluster),
                             label),
                    min_bv = min(bv), max_bv = max(bv),
                    .groups = 'drop')
  return(mean(max(df_b$min_bv), min(df_b$max_bv)))
}

digitalize_bv <- function(df_bv, bv_boundary = NULL) {
  bvb <- ifelse(is.null(bv_boundary),
                determine_bv_boundary(df_bv = df_bv), bv_boundary)
  message('>>> Digitalize beta-values by the boundary:\t', bvb)
  return(cbind(select(df_bv, 1),
               mutate_all(as_tibble(select(df_bv, -1) > bvb),
                          as.integer)))
}

df2num <- function(df) {
  return(as.vector(as.matrix(select_if(df, is.numeric))))
}


### segmer cluster

cluster_segments <- function(dmrbv_csv, dst_dir, k = 3,
                             dist_method = 'euclidean',
                             hclust_method = 'ward.D2', width = 13,
                             height = 8) {
  stopifnot(k > 1)
  distfun <- function(...) return(dist(..., method = dist_method))
  hclustfun <- function(...) return(hclust(..., method = hclust_method))
  out_prefix <- str_c(sub('.csv(|.gz)$', '', dmrbv_csv),
                      dist_method, hclust_method,
                      str_c('k', k), sep = '.')
  df_dmrbv <- read_csv_quietly(dmrbv_csv)
  hclust_csv <- str_c(out_prefix, '.hclust.csv')
  message('>>> Write an observed cluster CSV:\t', hclust_csv)
  mt_top <- as.matrix(column_to_rownames(df_dmrbv,
                                         var = colnames(df_dmrbv)[1]))
  hclust_hclusts <- stats::cutree(hclustfun(distfun(t(mt_top))), k = k)
  write_csv(tibble(sample_id = names(hclust_hclusts),
                   observed_cluster = hclust_hclusts),
            path = hclust_csv)
  heatmap_pdf <- str_c(out_prefix, '.heatmap.pdf')
  message('>>> Plot a heatmap:\t', heatmap_pdf)
  to_pdf(heatmap_plot(mt = mt_top, col_labels = hclust_hclusts,
                      distfun = distfun, hclustfun = hclustfun,
                      margins = c(5, 10)),
         path = heatmap_pdf, w = width, h = height)
}

heatmap_plot <- function(mt, col_labels, col = rev(brewer.pal(9, 'RdBu')),
                         ...) {
  return(heatmap.2(mt,
                   ColSideColors = brewer.pal(length(unique(col_labels)),
                                              'Accent')[col_labels],
                   scale = 'none', trace = 'none', density.info = 'none',
                   col = col, ...))
}


### segmer plot

visualize_bv <- function(site_csv, bv_csv, dst_dir, n_cpu = detectCores(),
                         width = 13, height = 8) {
  df_site <- read_csv_quietly(site_csv)
  df_bv <- read_bv_csv(path = bv_csv)
  out_prefix <- file.path(dst_dir,
                          str_c(sub('.csv(|.gz)$', '', basename(bv_csv)),
                                sub('(|.bed).csv(|.gz)$', '',
                                    basename(site_csv)),
                                sep = '.'))

  hist_pdf <- str_c(out_prefix, '.hist.pdf')
  message('>>> Plot a beta-value histogram:\t', hist_pdf)
  to_pdf(plot(gghistogram(rbind(tibble(probe = 'unmasked',
                                       bv = df2num(filter(df_bv,
                                                          name %in%
                                                            df_site$name))),
                                tibble(probe = 'masked',
                                       bv = df2num(filter(df_bv,
                                                          ! name %in%
                                                            df_site$name)))),
                          x = 'bv', fill = 'probe', position = 'stack',
                          color = NA, alpha = 0.6,
                          palette = get_palette('nejm', 2),
                          xlab = 'beta-value', ylab = 'site count',
                          bins = 30)),
         path = hist_pdf, w = width, h = height)

  cl <- makeCluster(n_cpu)
  df_p <- inner_join(df_site,
                     mutate(append_shapiro_wilk_p(cl,
                                                  filter(df_bv,
                                                         name %in%
                                                           df_site$name)),
                            mlogp = -log10(pval)),
                     by = 'name')
  stopCluster(cl)
  ylim <- c(0, max(df_p$mlogp) * 1.1)
  for (n in unique(df_site$chrom)) {
    p <- str_c(out_prefix, 'shapiro_wilk', n, 'pdf', sep = '.')
    message('>>> Plot p-values from Shapiro-Wilk normality tests:\t', p)
    to_pdf(plot(ggscatter(filter(df_p, chrom == n),
                          x = 'chromStart', y = 'mlogp', color = 'navy',
                          alpha = 0.4, size = 0.1, ylim = ylim,
                          xlab = str_c('position on ', n),
                          ylab = '-log10(p-value)') +
                scale_x_continuous(labels = format_digit)),
           path = p, w = width, h = height)
  }
}

visualize_segments <- function(seg_csv, site_csv, bv_csv, dmrbv_csv, dst_dir,
                               width = 13, height = 8) {
  df_seg <- read_csv_quietly(path = seg_csv)
  hist_pdf <- file.path(dst_dir, sub('.csv(|.gz)$', '.hist.pdf',
                                     basename(seg_csv)))
  message('>>> Plot segmental histograms:\t', hist_pdf)
  segment2len <- function(s) return(1 - eval(parse(text = sub('^.*:', '', s))))
  gs <- list(gghistogram(summarize(group_by(df_seg, segment),
                                   n = n_distinct(name),
                                   .groups = 'drop'),
                         x = 'n', fill = 'navy', color = NA, alpha = 0.4,
                         xscale = 'log10', xlab = 'site count per segment',
                         ylab = 'segment count', bins = 30) +
             theme(plot.margin = margin(1, 2.5, 1, 1, 'lines')),
             gghistogram(mutate(distinct(df_seg, segment),
                                len = sapply(segment, segment2len)),
                         x = 'len', fill = 'navy', color = NA, alpha = 0.4,
                         xlab = 'segment length (bp)', ylab = 'segment count',
                         bins = 30) +
             scale_x_log10(labels = format_digit) +
             theme(plot.margin = margin(1, 2.5, 1, 1, 'lines')))
  to_pdf(plot(ggarrange(plotlist = gs, nrow = 2)),
         path = hist_pdf, w = width, h = height)

  feature_csv <- file.path(dst_dir,
                           sub('.csv(|.gz)$', '.feature.csv',
                               basename(dmrbv_csv)))
  message('>>> Write feature counts:\t', feature_csv)
  df_bv_raw <- read_bv_csv(bv_csv, drop_na = FALSE)
  df_feature <- tibble(feature = c('differentially methylated segments',
                                   'determined segments',
                                   'available filtered sites',
                                   'filtered sites', 'probe sites'),
                       count = c(nrow(read_csv_quietly(dmrbv_csv)),
                                 length(unique(df_seg$segment)), nrow(df_seg),
                                 sum(read_csv_quietly(site_csv)$name %in%
                                     df_bv_raw$name),
                                 nrow(df_bv_raw)))
  write_csv(df_feature, path = feature_csv)

  feature_pdf <- sub('.csv$', '.pdf', feature_csv)
  message('>>> Write feature counts:\t', feature_csv)
  to_pdf(plot(ggbarplot(df_feature, x = 'feature', y = 'count', fill = 'navy',
                        color = 'navy', alpha = 0.4, orientation = 'horiz',
                        label = TRUE, lab.hjust = -0.2) +
              scale_y_continuous(limits = c(0, max(df_feature$count) * 1.1)) +
              labs(x = element_blank(), y = 'feature count') +
              theme(plot.margin = margin(1, 2.5, 1, 1, 'lines'))),
         path = feature_pdf, w = width, h = height)
}

format_digit <- function(i) {
  return(format(i, scientific = FALSE))
}


### dmp

write_dmp_bv <- function(bv_csv, site_csv, dst_dir, sd_cutoff = 0.1) {
  dmp_csv <- file.path(dst_dir,
                       str_c(sub('.csv(|.gz)$', '', basename(bv_csv)),
                             sub('(|.bed).csv(|.gz)$', '',
                                 basename(site_csv)),
                             'bv.dmp.csv', sep = '.'))
  df_bv <- filter(read_bv_csv(path = bv_csv),
                  name %in% read_csv_quietly(path = site_csv)$name)
  message('>>> Write DMP beta-values:\t', dmp_csv)
  write_csv(filter_bv_by_sd(df_bv, sd_cutoff = sd_cutoff), path = dmp_csv)
}


### utility

download_file <- function(src, dst) {
  if (! file.exists(dst)) {
    message('>>> Download a file:\t', src)
    download.file(src, dst)
    message('saved:\t', dst)
  }
}

granges2tibble <- function(gr, var = 'name') {
  return(as_tibble(rownames_to_column(as(gr, 'data.frame'), var = var)))
}

load_packages <- function(pkgs) {
  message('>>> Load packages')
  print(suppressMessages(sapply(pkgs, library, character.only = TRUE)))
}

make_dir <- function(path) {
  if (! (is.null(path) | dir.exists(path))) {
    message('>>> Make a directory:\t', path)
    dir.create(path, mode = '0755')
  }
}

read_csv_quietly <- function(path, ...) {
  message('>>> Read a CSV file:\t', path)
  stopifnot(file.exists(path))
  return(suppressMessages(read_csv(file = path, ...)))
}

read_bed <- function(path) {
  message('>>> Read a BED file:\t', path)
  stopifnot(file.exists(path))
  return(arrange(as_tibble(setNames(read.table(path, header = FALSE)[, 1:4],
                                    c('chrom', 'chromStart', 'chromEnd',
                                      'name'))),
                 chrom, chromStart, chromEnd))
}

to_pdf <- function(graph, path, w = 10, h = 10) {
  if (file.exists(path)) {
    file.remove(path)
  }
  pdf(path, width = w, height = h, onefile = FALSE)
  graph
  dev.off()
}


### CLI

if (! interactive()) {
  library('docopt', quietly = TRUE)
  main(opts = docopt::docopt(gsub('>\n +<', '> <',
                                  gsub('\\]\n +\\[', '\\] \\[', doc)),
                             version = command_version))
}
