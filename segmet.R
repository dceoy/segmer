#!/usr/bin/env Rscript

'Segment Detector for differentially methylated regions in methylome data

Usage:
  segmet.R bed [--debug] [--platform=<str>] [--unfilter] [--out=<dir>]
  segmet.R segment [--debug] [--seed=<int>] [--out=<dir>] <probe_bed> <met_csv>
  segmet.R [--debug] [--seed=<int>] [--cpus=<int>] [--out=<dir>]
  segmet.R --version
  segmet.R -h|--help

Commands:
  bed               Download annotation data and write a sorted probe BED file
  segment           Segment probes by sample variance
                    (Probes including NA are ignored.)

Options:
  --debug           Run with debug logging
  --platform=<str>  Specify a methylation assay platform [default: EPIC]
                      choice: EPIC, hm450, hm27
  --unfilter        Skip recommended probe filtering
  --seed=<int>      Set a random seed
  --cpus=<int>      Limit CPU cores for multithreading
  --out=<dir>       Set an output directory [default: .]
  --version         Print version and exit
  -h, --help        Print help and exit

Arguments:
  <probe_bed>       Path to a probe BED file with IDs in the name column
  <met_csv>         Path to a methylation CSV file
                      the 1st column: probe names
                      the other columns: sample values (e.g., beta values)
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
  if (opts[['--debug']]) {
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
  } else {
    n_cpu <- ifelse(is.null(opts[['--cpus']]),
                    parallel::detectCores(), as.integer(opts[['--cpus']]))
    message('>>> Show configurations')
    message('Script directory path:       ', root_dir)
    message('CPU cores:                   ', n_cpu)
  }
}

write_segment_files <- function(probe_bed, met_csv, dst_dir, method = 'PELT') {
  message('>>> Load data frames and calculate variance')
  df_var <- mutate(create_df_var(met_csv = met_csv, probe_bed = probe_bed),
                   sid = row_number())
  message('>>> Detect changepoints by ', method, ' method')
  meanvar_pelt <- changepoint::cpt.meanvar(df_var$variance, method = 'PELT')
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
  out_csv <- file.path(dst_dir,
                       str_c(sub('.csv(|.gz)$', '', basename(met_csv)),
                             sub('.bed(|.gz)$', '', basename(probe_bed)),
                             'csv', sep = '.'))
  message('>>> Write a segment CSV file:\t', out_csv)
  write_csv(df_ann, path = out_csv)
}

create_df_var <- function(met_csv, probe_bed) {
  lapply(c(probe_bed, met_csv),
         function(p) stopifnot(file.exists(p)))
  df_met <- mutate(drop_na(rename(read_csv(met_csv), name = 1)),
                   name = as.character(name))
  df_var <- inner_join(read_bed(bed = probe_bed),
                       mutate(select(df_met, name),
                              variance = apply(select(df_met, -name),
                                               1, var)),
                       by = 'name')
  stopifnot(nrow(df_met) > 1, ncol(df_met) > 4)
  message(ncol(df_var), ' sites are used.')
  return(df_var)
}

read_bed <- function(bed) {
  d <- read.table(bed, header = FALSE)[, 1:4]
  return(arrange(rename(as_tibble(d),
                        setNames(names(d),
                                 c('chrom', 'chromStart', 'chromEnd',
                                   'name'))),
                 chrom, chromStart, chromEnd))
}

write_probe_bed <- function(dst_dir, platform = 'EPIC', unfilter = FALSE,
                            hg_ver = 'hg38') {
  rds_files <- download_annotation_data(dst_dir = dst_dir)
  if (unfilter) {
    df_ann <- granges2tibble(readRDS(rds_files[2]), var = 'cg_id')
    bed_name <- str_c(platform, hg_ver, 'bed', sep = '.')
  } else {
    message('>>> Filter probes')
    lapply(c('chrX', 'chrY', 'MASK_general', 'MASK_snp5_common'),
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
  urls <- sapply(c(str_c(platform, hg_ver, 'manifest.rds', sep = '.'),
                   str_c(platform, hg_ver, 'manifest.gencode',
                         gencode_ver,
                         'rds', sep = '.')),
                 function(f) {
                   return(file.path('http://zwdzwd.io/InfiniumAnnotation',
                                    'current', platform, f))
                 })
  dsts <- sapply(names(urls),
                 function(f) return(file.path(dst_dir, f)))
  lapply(names(urls),
         function(f) {
           src <- urls[f]
           dst <- dsts[f]
           if (! file.exists(dst)) {
             message('>>> Download:\t', src, ' => ', dst)
             download.file(src, dst)
           }
         })
  return(dsts)
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

if (! interactive()) {
  library('docopt', quietly = TRUE)
  main(opts = docopt::docopt(doc, version = script_version))
}
