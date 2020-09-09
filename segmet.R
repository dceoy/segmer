#!/usr/bin/env Rscript

'Segment Detector for differentially methylated regions in methylome data

Usage:
  segmet.R bed [-v] [--platform=<str>] [--unfilter] [--out=<dir>]
  segmet.R segment [-v] [--seed=<int>] [--avgfun=<str>] [--ar=<ratio>]
                   [--out=<dir>] <site_csv> <met_csv>
  segmet.R cluster [-v] [--k=<int>] [--drop=<dbl>] [--dist=<str>]
                   [--hclust=<str>] [--ar=<ratio>] [--out=<dir>]
                   <segbv_csv>
  segmet.R --session
  segmet.R --version
  segmet.R -h|--help

Commands:
  bed               Download annotation data and write a target site CSV file
  segment           Segment target sites by sample variance
                    (Sites including NA are ignored.)
  cluster           Execute clustering and draw the heatmap

Options:
  -v                Run with debug logging
  --platform=<str>  Specify a methylation assay platform [default: EPIC]
                      choice: EPIC, hm450, hm27
  --unfilter        Skip recommended probe filtering
  --out=<dir>       Set an output directory [default: .]
  --seed=<int>      Set a random seed
  --avgfun=<str>    Specify the function for segment average [default: median]
  --k=<int>         Specify the number of clusters [default: 3]
  --drop=<dbl>      Specify the proportion of the ignored segments [default: 0]
  --dist=<str>      Specify the method of stats::dist [default: euclidean]
  --hclust=<str>    Specify the method of stats::hclust [default: ward.D2]
  --ar=<ratio>      Specify the aspect ratio of figures [default: 13:8]
  --session         Print session information and exit ({devtools} is required)
  --version         Print version and exit
  -h, --help        Print help and exit

Arguments:
  <site_csv>        Path to a target site CSV file (created by `segmet bed`)
  <met_csv>         Path to a methylation CSV file
                      the 1st column: probe names
                      the other columns: beta-values
  <segbv_csv>       Path to a segmental beta-value CSV file
                    (`*.seg.bv.median.dmr.csv` created by `segmet segment`)
' -> doc


command_version <- 'v0.0.7'

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
    } else if (opts[['segment']]) {
      add_pkgs <- c('changepoint', 'ggpubr')
    } else if (opts[['cluster']]) {
      add_pkgs <- c('gplots', 'RColorBrewer')
    } else {
      add_pkgs <- NULL
    }
    load_packages(pkgs = c('tidyverse', add_pkgs))
    make_dir(path = opts[['--out']])
    dst_dir <- normalizePath(opts[['--out']])
    if (opts[['segment']] | opts[['cluster']]) {
      aspect_ratio <- as.integer(str_split(opts[['--ar']], pattern = ':',
                                           n = 2, simplify = TRUE))
    } else {
      aspect_ratio <- NULL
    }

    if (opts[['bed']]) {
      prepare_site_csv(dst_dir = dst_dir,
                       platform = opts[['--platform']],
                       unfilter = opts[['--unfilter']])
    } else if (opts[['segment']]) {
      if (! is.null(opts[['--seed']])) {
        message('>>> Set a random seed')
        set.seed(opts[['--seed']])
        message(opts[['--seed']])
      }
      site_csv <- normalizePath(opts[['<site_csv>']])
      met_csv <- normalizePath(opts[['<met_csv>']])
      visualize_bv_var(site_csv = site_csv, met_csv = met_csv,
                       dst_dir = dst_dir, width = aspect_ratio[1],
                       height = aspect_ratio[2])
      csv_list <- segment_sites(site_csv = site_csv, met_csv = met_csv,
                                dst_dir = dst_dir)
      visualize_segments(seg_csv = csv_list$seg_csv, dst_dir = dst_dir,
                         width = aspect_ratio[1], height = aspect_ratio[2])
      calculate_segment_segbv(seg_csv = csv_list$seg_csv, met_csv = met_csv,
                              dst_dir = dst_dir,
                              bv_boundary = scan(csv_list$bv_boundary_txt),
                              avg = opts[['--avgfun']])
    } else if (opts[['cluster']]) {
      cluster_segments(segbv_csv = normalizePath(opts[['<segbv_csv>']]),
                       dst_dir = dst_dir, k = opts[['--k']],
                       drop_prop = as.numeric(opts[['--drop']]),
                       dist_method = opts[['--dist']],
                       hclust_method = opts[['--hclust']],
                       width = aspect_ratio[1], height = aspect_ratio[2])
    }
  }
}

cluster_segments <- function(segbv_csv, dst_dir, k = 3, drop_prop = 0,
                             dist_method = 'euclidean',
                             hclust_method = 'ward.D2', width = 13,
                             height = 8) {
  stopifnot((k > 1) & (drop_prop >= 0) & (drop_prop < 1))
  distfun <- function(...) return(dist(..., method = dist_method))
  hclustfun <- function(...) return(hclust(..., method = hclust_method))
  drop_percent <- as.integer(drop_prop * 100)
  out_prefix <- str_c(sub('.csv(|.gz)$', '', segbv_csv),
                      str_c('drop', drop_percent),
                      dist_method, hclust_method,
                      str_c('k', k), sep = '.')
  df_segbv <- read_csv_quietly(segbv_csv)
  n_seg <- nrow(df_segbv)
  if (drop_prop > 0) {
    message('>>> Drop low-variance segments (', drop_percent, '%)')
    df_segbv <- select(slice_max(mutate(df_segbv,
                                        s2 = apply(select(df_segbv, -segment),
                                                   1, var)),
                                 s2,
                                 prop = 1 - drop_prop),
                       -s2)
  }
  message('used segments:\t', nrow(df_segbv), ' / ', n_seg)
  used_csv <- str_c(out_prefix, '.used.csv')
  message('>>> Write used segments CSV:\t', used_csv)
  write_csv(df_segbv, path = used_csv)
  hclust_csv <- str_c(out_prefix, '.hclust.csv')
  message('>>> Write an observed cluster CSV:\t', hclust_csv)
  mt_top <- as.matrix(column_to_rownames(df_segbv, var = 'segment'))
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
  return(list(used_csv = used_csv, hclust_csv = hclust_csv))
}

heatmap_plot <- function(mt, col_labels, col = rev(brewer.pal(9, 'RdBu')),
                         ...) {
  return(heatmap.2(mt,
                   ColSideColors = brewer.pal(length(unique(col_labels)),
                                              'Accent')[col_labels],
                   scale = 'none', trace = 'none', density.info = 'none',
                   col = col, ...))
}

to_pdf <- function(graph, path, w = 10, h = 10) {
  if (file.exists(path)) {
    file.remove(path)
  }
  pdf(path, width = w, height = h, onefile = FALSE)
  graph
  dev.off()
}

calculate_segment_segbv <- function(seg_csv, met_csv, dst_dir,
                                    bv_boundary = NULL, avg = 'median') {
  df_seg <- read_csv_quietly(path = seg_csv)
  df_bv <- read_bv_csv(path = met_csv)
  message('>>> Calculate segmental ', avg, ' values')
  df_segbv <- summarize_all(group_by(select(left_join(select(df_seg,
                                                             name, segment),
                                                      df_bv, by = 'name'),
                                            -name),
                                     segment),
                            eval(parse(text = avg)),
                            .groups = 'drop')
  out_prefix <- file.path(dst_dir,
                          str_c(sub('.csv(|.gz)$', '', basename(seg_csv)),
                                'bv', avg, sep = '.'))
  segbv_csv <- str_c(out_prefix, '.all.csv')
  message('>>> Write a segmental ', avg, ' CSV:\t', segbv_csv)
  write_csv(df_segbv, path = segbv_csv)
  segd_csv <- str_c(out_prefix, '.digital.csv')
  message('>>> Write a digitalized segmental ', avg, ' CSV:\t', segd_csv)
  df_segd <- digitalize_bv(df_segbv, bv_boundary)
  write_csv(df_segd, path = segd_csv)
  dmr_csv <- str_c(out_prefix, '.dmr.csv')
  message('>>> Write a DMR segmental ', avg, ' CSV:\t', dmr_csv)
  df_dmr <- inner_join(select(filter(bv2var(df_segd), variance > 0),
                              -variance),
                       df_segbv, by = 'segment')
  message('DMR segments:\t', nrow(df_dmr), ' / ', nrow(df_segbv))
  write_csv(df_dmr, path = dmr_csv)
  return(dmr_csv)
}

determine_bv_boundary <- function(df_bv, ...) {
  message('>>> Determine the boundary by K-means')
  bv <- as.vector(as.matrix(select(df_bv, -1)))
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

bv2var <- function(df_bv) {
  return(mutate(select(df_bv, 1), variance = apply(select(df_bv, -1), 1, var)))
}

segment_sites <- function(site_csv, met_csv, dst_dir, method = 'PELT',
                          ...) {
  df_site <- read_csv_quietly(site_csv)
  df_bv <- read_bv_csv(path = met_csv)
  bv_boundary <- determine_bv_boundary(df_bv)
  bv_boundary_txt <- file.path(dst_dir,
                               sub('.csv(|.gz)$', '.boundary.txt',
                                   basename(met_csv)))
  message('>>> Write the beta-value boundary:\t', bv_boundary_txt)
  write(bv_boundary, file = bv_boundary_txt)
  df_var <- mutate(inner_join(df_site,
                              bv2var(digitalize_bv(df_bv, bv_boundary)),
                              by = 'name'),
                   sid = row_number())
  message(nrow(df_var), ' sites are used.')

  message('>>> Detect changepoints by ', method, ' method')
  cpt <- changepoint::cpt.meanvar(df_var$variance, method = method, ...)
  summary(cpt)
  cp_sids <- c(1, changepoint::cpts(cpt))
  message('>>> Segment target sites')
  df_reg <- summarize(group_by(fill(left_join(dplyr::rename(df_var,
                                                            cp_sid = sid),
                                              tibble(sid = cp_sids,
                                                     cp_sid = cp_sids),
                                              by = 'cp_sid'),
                                    sid),
                               CGI, CGIposition, sid),
                      segment = str_c(unique(chrom), ':',
                                      min(chromStart) + 1,
                                      '-', max(chromEnd)),
                      .groups = 'drop')
  df_seg <- fill(left_join(select(df_var, -chrom, -chromStart, -chromEnd),
                           df_reg,
                           by = c('CGI', 'CGIposition', 'sid')),
                 segment)
  print(summary(summarize(group_by(df_seg, segment),
                          sites_per_segment = n(),
                          .groups = 'drop')))

  seg_csv <- file.path(dst_dir,
                       str_c(sub('.csv(|.gz)$', '', basename(met_csv)),
                             sub('.csv(|.gz)$', '', basename(site_csv)),
                             'seg.csv', sep = '.'))
  message('>>> Write a segment CSV file:\t', seg_csv)
  write_csv(df_seg, path = seg_csv)
  return(list(seg_csv = seg_csv, bv_boundary_txt = bv_boundary_txt))
}

visualize_segments <- function(seg_csv, dst_dir, width = 13, height = 8) {
  df_seg <- read_csv_quietly(path = seg_csv)
  out_prefix <- file.path(dst_dir, sub('.csv(|.gz)$', '', basename(seg_csv)))
  pdf_list <- as.list(setNames(str_c(out_prefix, c('n_site', 'len'), 'pdf',
                                     sep = '.'),
                               c('n_site', 'len')))
  message('>>> Plot segmental site count:\t', pdf_list$n_site)
  to_pdf(plot(gghistogram(summarize(group_by(df_seg, segment),
                                    n = n_distinct(name),
                                    .groups = 'drop'),
                          x = 'n', fill = 'navy', color = NA, alpha = 0.4,
                          xlab = 'site count per segment',
                          ylab = 'segment count', bins = 30)),
         path = pdf_list$n_site, w = width, h = height)
  message('>>> Plot segmental lengths:\t', pdf_list$len)
  segment2len <- function(s) return(1 - eval(parse(text = sub('^.*:', '', s))))
  to_pdf(plot(gghistogram(mutate(distinct(df_seg, segment),
                                 len = sapply(segment, segment2len)),
                          x = 'len', fill = 'navy', color = NA, alpha = 0.4,
                          xlab = 'segment length (bp)', ylab = 'segment count',
                          bins = 30)),
         path = pdf_list$len, w = width, h = height)
}

read_bv_csv <- function(path) {
  df_bv <- mutate(drop_na(dplyr::rename(read_csv_quietly(path), name = 1)),
                  name = as.character(name))
  stopifnot(nrow(df_bv) > 1, ncol(df_bv) > 4)
  return(df_bv)
}

read_bed <- function(path) {
  message('>>> Read a BED file: ', path)
  stopifnot(file.exists(path))
  return(arrange(as_tibble(setNames(read.table(path, header = FALSE)[, 1:4],
                                    c('chrom', 'chromStart', 'chromEnd',
                                      'name'))),
                 chrom, chromStart, chromEnd))
}

visualize_bv_var <- function(site_csv, met_csv, dst_dir, width = 13,
                             height = 8) {
  df_site <- read_csv_quietly(site_csv)
  df_bv <- filter(read_bv_csv(path = met_csv), name %in% df_site$name)
  df_var <- inner_join(select(dplyr::rename(df_site, position = chromStart),
                              chrom, position, name, genesUniq, CGI,
                              CGIposition),
                       bv2var(df_bv),
                       by = 'name')
  ylim <- c(0, max(df_var$variance) * 1.2)
  out_prefix <- file.path(dst_dir, sub('.csv(|.gz)$', '', basename(met_csv)))

  hist_pdf <- str_c(out_prefix, '.hist.pdf')
  message('>>> Plot a beta-value histogram:\t', hist_pdf)
  to_pdf(plot(gghistogram(tibble(v = as.vector(as.matrix(select(df_bv,
                                                                -name)))),
                          x = 'v', fill = 'navy', color = NA, alpha = 0.4,
                          xlab = 'beta-value', ylab = 'site count',
                          bins = 30)),
         path = hist_pdf, w = width, h = height)

  digit2str <- function(i) return(format(i, scientific = FALSE))
  for (n in unique(df_site$chrom)) {
    p <- str_c(out_prefix, 'var', n, 'pdf', sep = '.')
    message('>>> Plot beta-value variances:\t', p)
    to_pdf(plot(ggscatter(filter(df_var, chrom == n),
                          x = 'position', y = 'variance', color = 'navy',
                          alpha = 0.1, size = 0.1, ylim = ylim,
                          xlab = str_c('position on ', n),
                          ylab = 'beta-value variance') +
                scale_x_continuous(labels = digit2str)),
           path = p, w = width, h = height)
  }
}

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
    bed_csv_name <- str_c(platform, hg_ver, 'filtered.bed.csv.gz', sep = '.')
  }
  site_csv <- file.path(dst_dir, bed_csv_name)
  message('>>> Write a sorted target site BED file:\t', site_csv)
  write_csv(select(mutate(dplyr::rename(df_ann,
                                        chrom = seqnames, chromEnd = end),
                          chromStart = start - 1),
                   chrom, chromStart, chromEnd, name, genesUniq, CGI,
                   CGIposition),
            path = site_csv)
  return(site_csv)
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
  message('>>> Read a CSV file: ', path)
  stopifnot(file.exists(path))
  return(suppressMessages(read_csv(file = path, ...)))
}

if (! interactive()) {
  library('docopt', quietly = TRUE)
  main(opts = docopt::docopt(gsub('\\]\n +([\\[<])', '] \\1', doc),
                             version = command_version))
}
