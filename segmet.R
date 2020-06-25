#!/usr/bin/env Rscript

'Segment Detector for differentially methylated regions in methylome data

Usage:
  segmet.R bed [-v] [--platform=<str>] [--unfilter] [--out=<dir>]
  segmet.R segment [-v] [--seed=<int>] [--out=<dir>] <probe_bed> <met_csv>
  segmet.R heatmap [-v] [--out=<dir>] [--k=<int>] [--cutoff=<dbl>] <segmet_csv>
  segmet.R --version
  segmet.R -h|--help

Commands:
  bed               Download annotation data and write a sorted probe BED file
  segment           Segment probes by sample variance
                    (Probes including NA are ignored.)
  heatmap           Draw a heatmap using detected segments

Options:
  -v                Run with debug logging
  --platform=<str>  Specify a methylation assay platform [default: EPIC]
                      choice: EPIC, hm450, hm27
  --unfilter        Skip recommended probe filtering
  --seed=<int>      Set a random seed
  --out=<dir>       Set an output directory [default: .]
  --k=<int>         Specify the number of clusters [default: 3]
  --cutoff=<dbl>    Specify the minimum cutoff for segmental values [default: 0]
  --version         Print version and exit
  -h, --help        Print help and exit

Arguments:
  <probe_bed>       Path to a probe BED file with IDs in the name column
  <met_csv>         Path to a methylation CSV file
                      the 1st column: probe names
                      the other columns: sample values (e.g., beta values)
  <segmet_csv>      Path to the CSV file created by `segmet segment`
' -> doc


script_version <- 'v0.0.1'

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
  if (! is.null(opts[['--seed']])) {
    message('>>> Set a random seed')
    set.seed(opts[['--seed']])
    message(opts[['--seed']])
  }

  if (opts[['bed']]) {
    add_pkgs <- 'GenomicRanges'
  } else if (opts[['segment']]) {
    add_pkgs <- 'changepoint'
  } else if (opts[['heatmap']]) {
    add_pkgs <- c('gplots', 'RColorBrewer')
  } else {
    add_pkgs <- NULL
  }
  load_packages(pkgs = c('tidyverse', add_pkgs))
  make_dir(path = opts[['--out']])

  if (opts[['bed']]) {
    write_probe_bed(dst_dir = normalizePath(opts[['--out']]),
                    platform = opts[['--platform']],
                    unfilter = opts[['--unfilter']])
  } else if (opts[['segment']]) {
    write_segment_files(probe_bed = normalizePath(opts[['<probe_bed>']]),
                        met_csv = normalizePath(opts[['<met_csv>']]),
                        dst_dir = normalizePath(opts[['--out']]))
  } else if (opts[['heatmap']]) {
    draw_segment_heatmap(segmet_csv =  normalizePath(opts[['<segmet_csv>']]),
                         dst_dir = normalizePath(opts[['--out']]),
                         k = opts[['--k']],
                         min_cutoff = as.numeric(opts[['--cutoff']]),
                         distfun = dist, hclustfun = ward_hclust)
  }
}

draw_segment_heatmap <- function(segmet_csv, dst_dir, k = 3, min_cutoff = 0,
                                 distfun = dist, hclustfun = hclust) {
  out_prefix <- sub('.csv(|.gz)$', '', segmet_csv)
  message('>>> Load a data frame')
  validate_input_files(paths = c('segmental methylation CSV' = segmet_csv))
  cluster_csv <- str_c(out_prefix, '.', k, 'clusters.csv')
  df_segmet <- column_to_rownames(read_csv_quietly(segmet_csv),
                                  var = 'segment')
  if (min_cutoff > 0) {
    mt_segmet <- as.matrix(filter(df_segmet,
                                  apply(df_segmet, 1, max) >=  min_cutoff))
  } else {
    mt_segmet <- as.matrix(df_segmet)
  }
  hclust_labels <- stats::cutree(hclustfun(distfun(t(mt_segmet))), k = k)
  message('>>> Write an observed cluster CSV:\t', cluster_csv)
  write_csv(tibble(sample_id = names(hclust_labels),
                   observed_cluster = hclust_labels),
            path = cluster_csv)
  heatmap_pdf <- str_c(out_prefix, '.', k, 'clusters.heatmap.pdf')
  message('>>> Draw a heatmap:\t', heatmap_pdf)
  to_pdf(heatmap_plot(mt = mt_segmet, col_labels = hclust_labels,
                      distfun = distfun, hclustfun = hclustfun,
                      margins = c(5, 10)),
         path = heatmap_pdf)
}

heatmap_plot <- function(mt, col_labels, col = rev(brewer.pal(9, 'RdBu')),
                         ...) {
  return(heatmap.2(mt,
                   ColSideColors = brewer.pal(length(unique(col_labels)),
                                              'Accent')[col_labels],
                   scale = 'none', trace = 'none', density.info = 'none',
                   col = col, ...))
}

ward_hclust <- function(...) {
  return(hclust(..., method = 'ward.D2'))
}

to_pdf <- function(graph, path, w = 10, h = 10) {
  if (file.exists(path)) {
    file.remove(path)
  }
  pdf(path, width = w, height = h, onefile = FALSE)
  graph
  dev.off()
}

write_segment_files <- function(probe_bed, met_csv, dst_dir, method = 'PELT',
                                avg = 'median', ...) {
  message('>>> Load data frames and calculate variance')
  validate_input_files(paths = c('probe BED' = probe_bed,
                                 'methylation CSV' = met_csv))
  df_var <- create_df_var(df_met = read_met_csv(path = met_csv),
                          df_probe = read_bed(path = probe_bed))
  message('>>> Detect changepoints by ', method, ' method')
  meanvar_pelt <- changepoint::cpt.meanvar(df_var$variance, method = 'PELT',
                                           ...)
  summary(meanvar_pelt)
  cp_sids <- c(1, changepoint::cpts(meanvar_pelt))
  message('>>> Segment probes')
  df_seg <- summarize(group_by(fill(left_join(rename(df_var,
                                                     cp_sid = sid),
                                              tibble(sid = cp_sids,
                                                     cp_sid = cp_sids),
                                              by = 'cp_sid'),
                                    sid),
                               chrom, sid),
                      segment = str_c(unique(chrom), ':',
                                      min(chromStart) + 1,
                                      '-', max(chromEnd)),
                      .groups = 'drop')
  df_ann <- fill(left_join(select(df_var, -chromStart, -chromEnd),
                           df_seg,
                           by = c('chrom', 'sid')),
                 segment)
  print(summary(summarize(group_by(df_ann, segment),
                          probes_per_segment = n(),
                          .groups = 'drop')))
  ann_csv <- file.path(dst_dir,
                       str_c(sub('.csv(|.gz)$', '', basename(met_csv)),
                             sub('.bed(|.gz)$', '', basename(probe_bed)),
                             'seg.csv', sep = '.'))
  message('>>> Write a segment CSV file:\t', ann_csv)
  write_csv(df_ann, path = ann_csv)
  message('>>> Calculate segmental ', avg, ' values by segment')
  df_segmet <- create_df_segmet(df_ann = df_ann,
                                df_met = read_met_csv(met_csv),
                                avgfun = eval(parse(text = avg)))
  segmet_csv <- sub('.csv$', str_c('.met.', avg, '.csv'), ann_csv)
  message('>>> Write an segmental methylation CSV:\t', segmet_csv)
  write_csv(df_segmet, path = segmet_csv)
}

create_df_segmet <- function(df_ann, df_met, avgfun = median) {
  return(summarize_all(group_by(select(left_join(df_ann, df_met, by = 'name'),
                                       -name, -chrom, -variance, -sid),
                                segment),
                       avgfun, .groups = 'drop'))
}

read_met_csv <- function(path) {
  df_met <- mutate(drop_na(rename(read_csv_quietly(path), name = 1)),
                   name = as.character(name))
  stopifnot(nrow(df_met) > 1, ncol(df_met) > 4)
  return(df_met)
}

read_bed <- function(path) {
  d <- read.table(path, header = FALSE)[, 1:4]
  return(arrange(rename(as_tibble(d),
                        set_names(names(d),
                                  c('chrom', 'chromStart', 'chromEnd',
                                    'name'))),
                 chrom, chromStart, chromEnd))
}

create_df_var <- function(df_met, df_probe) {
  df_var <- inner_join(df_probe,
                       mutate(select(df_met, name),
                              variance = apply(select(df_met, -name),
                                               1,
                                               function(v) {
                                                 return(var(v))
                                               })),
                       by = 'name')
  message(nrow(df_var), ' sites are used.')
  return(mutate(df_var, sid = row_number()))
}

write_probe_bed <- function(dst_dir, platform = 'EPIC', unfilter = FALSE,
                            hg_ver = 'hg38') {
  rds_files <- download_annotation_data(dst_dir = dst_dir)
  if (unfilter) {
    df_ann <- granges2tibble(readRDS(rds_files[2]), var = 'cg_id')
    bed_name <- str_c(platform, hg_ver, 'bed', sep = '.')
  } else {
    message('>>> Filter probes')
    walk(c('chrX', 'chrY', 'MASK_general', 'MASK_snp5_common'),
         function(i) message('- ', i))
    masked_ids <- filter(granges2tibble(readRDS(rds_files[1]), var = 'cg_id'),
                         seqnames %in% c('chrX', 'chrY')
                         | MASK_general | MASK_snp5_common)$cg_id
    df_ann <- filter(granges2tibble(readRDS(rds_files[2]), var = 'cg_id'),
                     ! seqnames %in% c('chrX', 'chrY'),
                     ! cg_id %in% masked_ids)
    bed_name <- str_c(platform, hg_ver, 'filtered.bed', sep = '.')
  }
  out_csv <- file.path(dst_dir, bed_name)
  message('>>> Write a sorted probe BED file:\t', out_csv)
  write_tsv(mutate(select(df_ann, seqnames, start, end, cg_id),
                   start = start - 1,
                   end = end),
            path = out_csv, col_names = FALSE)
}

download_annotation_data <- function(dst_dir = '.', platform = 'EPIC',
                                     hg_ver = 'hg38', gencode_ver = 'v22') {
  url_body <- file.path('http://zwdzwd.io/InfiniumAnnotation/current',
                        platform)
  return(map_chr(c(str_c(platform, hg_ver, 'manifest.rds', sep = '.'),
                   str_c(platform, hg_ver, 'manifest.gencode', gencode_ver,
                         'rds', sep = '.')),
                 function(n) {
                   src <- file.path(url_body, n)
                   dst <- file.path(dst_dir, n)
                   if (! file.exists(dst)) {
                     message('>>> Download:\t', src, ' => ', dst)
                     download.file(src, dst)
                   }
                   return(dst)
                 }))
}

granges2tibble <- function(gr, var) {
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

validate_input_files <- function(paths) {
  walk2(names(paths), paths,
        function(n, p) {
          message(n, ':\t', p)
          stopifnot(file.exists(p))
        })
}

read_csv_quietly <- function(...) {
  return(suppressMessages(read_csv(...)))
}

if (! interactive()) {
  library('docopt', quietly = TRUE)
  main(opts = docopt::docopt(doc, version = script_version))
}
