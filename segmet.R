#!/usr/bin/env Rscript

'Segment Detector for differentially methylated regions in methylome data

Usage:
  segmet.R preproc [--debug] [--platform=<std>] [--unfilter] [--out=<dir>]
  segmet.R [--debug] [--seed=<int>] [--cpus=<int>] [--in=<dir>] [--out=<dir>]
  segmet.R --version
  segmet.R -h|--help

Commands:
  preproc           Download annotation data and write a sorted probe BED

Options:
  --debug           Run with debug logging
  --platform=<std>  Speciify a methylation assay platform [default: EPIC]
                    { EPIC, hm450, hm27 }
  --unfilter        Skip recommended probe filtering
  --seed=<int>      Set a random seed
  --cpus=<int>      Limit CPU cores for multithreading
  --in=<dir>        Set an input directory [default: .]
  --out=<dir>       Set an output directory [default: .]
  --version         Print version and exit
  -h, --help        Print help and exit' -> doc


script_version <- 'v0.0.1'

fetch_script_root <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  fa <- ca[grepl('^--file=', ca)]
  if (length(fa) == 1) {
    f <- sub('--file=', '', fa)
    l <- Sys.readlink(f)
    if (is.na(l)) {
      script_path <- normalizePath(f)
    } else if (startsWith(l, '/')) {
      script_path <- normalizePath(l)
    } else {
      script_path <- normalizePath(file.path(dirname(f), l))
    }
    return(dirname(script_path))
  } else {
    return(normalizePath(getwd()))
  }
}

main <- function(opts, root_dir = fetch_script_root()) {
  options(warn = 1)
  if (opts[['--debug']]) {
    print(opts)
  }

  if (opts[['preproc']]) {
    write_probe_bed(dir = normalizePath(opts[['--out']]),
                    platform = opts[['--platform']],
                    unfilter = opts[['--unfilter']])
  } else {
    if (! is.null(opts[['--seed']])) {
      message('>>> Set a random seed')
      set.seed(opts[['--seed']])
      message(opts[['--seed']])
    }

    load_packages(pkgs = c('stringr', 'tidyverse'))

    n_cpu <- ifelse(is.null(opts[['--cpus']]),
                    parallel::detectCores(), as.integer(opts[['--cpus']]))
    in_dir <- normalizePath(opts[['--in']])
    out_dir <- normalizePath(opts[['--out']])

    message('>>> Show configurations')
    message('Script directory path:       ', root_dir)
    message('Input directory path:        ', in_dir)
    message('Output directory path:       ', out_dir)
    message('CPU cores:                   ', n_cpu)
  }
}

write_probe_bed <- function(dir, platform = 'EPIC', unfilter = FALSE,
                            hg_ver = 'hg38') {
  message(str_c('>>> Download annotation data: ', platform))
  load_packages(pkgs = c('GenomicRanges', 'stringr', 'tidyverse'))
  rds_files <- fetch_epic_annotation(dir = dir)
  if (unfilter) {
    df_ann <- granges2tibble(readRDS(rds_files[2]), var = 'cg_id')
    bed_name <- str_c(platform, hg_ver, 'bed', sep = '.')
  } else {
    message('>>> Filter probes')
    masked_ids <- filter(granges2tibble(readRDS(rds_files[1]), var = 'cg_id'),
                         seqnames %in% c('chrX', 'chrY')
                         | MASK_general | MASK_snp5_common)$cg_id
    df_ann <- filter(granges2tibble(readRDS(rds_files[2]), var = 'cg_id'),
                     ! seqnames %in% c('chrX', 'chrY'),
                     ! cg_id %in% masked_ids)
    bed_name <- str_c(platform, hg_ver, 'filtered.bed', sep = '.')
  }
  message('>>> Write a sorted probe BED')
  write_csv(mutate(select(df_ann, seqnames, start, end, cg_id),
                   start = start - 1,
                   end = end - 1),
            path = file.path(dir, bed_name),
            col_names = FALSE)
}

load_packages <- function(pkgs) {
  message('>>> Load packages')
  print(suppressMessages(sapply(pkgs, library, character.only = TRUE)))
}

fetch_epic_annotation <- function(dir = '.', platform = 'EPIC',
                                  hg_ver = 'hg38', gencode_ver = 'v22',
                                  release_date = '20180909') {
  urls <- sapply(c(str_c(platform, hg_ver, 'manifest.rds', sep = '.'),
                   str_c(platform, hg_ver, 'manifest.gencode',
                         gencode_ver,
                         'rds', sep = '.')),
                 function(f) {
                   return(file.path('https://zwdzwd.io/InfiniumAnnotation',
                                    release_date, platform, f))
                 })
  dsts <- sapply(names(urls),
                 function(f) return(file.path(dir, f)))
  lapply(names(urls),
         function(f) {
           src <- urls[f]
           dst <- dsts[f]
           if (! file.exists(dst)) {
             message(str_c('Download:\t', src, ' => ', dst))
             download.file(src, dst)
           }
         })
  return(dsts)
}

granges2tibble <- function(gr, var) {
  return(as_tibble(rownames_to_column(as(gr, 'data.frame'), var = var)))
}

if (! interactive()) {
  library('docopt', quietly = TRUE)
  main(opts = docopt::docopt(doc, version = script_version))
}
