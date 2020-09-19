#!/usr/bin/env Rscript

'Segment Detector for differentially methylated regions in methylome data

Usage:
  segmer bed [-v] [--platform=<str>] [--unfilter] [--out=<dir>]
  segmer dmr [-v] [--seed=<int>] [--sd-cutoff=<dbl>] [--ar=<ratio>] [--out=<dir>] <site_csv>
             <mv_csv>
  segmer cluster [-v] [--k=<int>] [--dist=<str>] [--hclust=<str>] [--ar=<ratio>] [--out=<dir>]
                 <dmrmv_csv>
  segmer plot [-v] [--ar=<ratio>] [--out=<dir>] <site_csv> <mv_csv> <seg_csv> <dmrmv_csv>
  segmer dmp [-v] [--sd-cutoff=<dbl>] [--out=<dir>] <site_csv> <mv_csv>
  segmer idat2m [-v] [--offset=<dbl>] <idat_dir> <out_csv>
  segmer idat2beta [-v] [--offset=<dbl>] <idat_dir> <out_csv>
  segmer --session
  segmer --version
  segmer -h|--help

Commands:
  bed                 Download annotation data and write a target site CSV file
  dmr                 Segment target sites (Sites including NA are ignored.)
  cluster             Execute clustering and draw the heatmap
  plot                Visualize methylation data
  dmp                 Determine differentially methylated positions
  idat2m              Calculate M-values using IDAT files from Illumina methylation arrays
  idat2beta           Calculate beta-values using IDAT files from Illumina methylation arrays

Options:
  -v                  Run with debug logging
  --platform=<str>    Specify a methylation assay platform [default: EPIC]
                        choice: EPIC, hm450, hm27
  --unfilter          Skip recommended probe filtering
  --out=<dir>         Specify an output directory [default: .]
  --seed=<int>        Specify a random seed
  --sd-cutoff=<dbl>   Specify the SD cutoff
  --k=<int>           Specify the number of clusters [default: 3]
  --dist=<str>        Specify the method of stats::dist [default: euclidean]
  --hclust=<str>      Specify the method of stats::hclust [default: ward.D2]
  --ar=<ratio>        Specify the aspect ratio of figures [default: 13:8]
  --offset=<dbl>      Specify the offset for M-values [default: 1]
  --session           Print session information and exit (using {devtools})
  --version           Print version and exit
  -h, --help          Print help and exit

Arguments:
  <site_csv>          Path to a target site CSV or BED file (created by `segmer bed`)
  <mv_csv>            Path to a methylation CSV file
                        the 1st column: probe names
                        the other columns: M-values
  <dmrmv_csv>         Path to a M-value CSV file of DMRs or DMPs
                      (`*.seg.mv.dmr.csv` created by `segmer dmr`)
  <seg_csv>           Path to a segment CSV file (`*.seg.csv` created by `segmer dmr`)
  <idat_dir>          Path to an directory including IDAT files
  <out_csv>           Path to an output CSV file
' -> doc

command_version <- 'v0.1.0'


### controler

load_packages <- function(opts) {
  message('>>> Load packages')
  if (opts[['bed']]) {
    add_pkgs <- 'GenomicRanges'
  } else if (opts[['dmr']]) {
    add_pkgs <- c('changepoint', 'parallel')
  } else if (opts[['cluster']]) {
    add_pkgs <- c('gplots', 'RColorBrewer')
  } else if (opts[['plot']]) {
    add_pkgs <- c('ggpubr', 'parallel')
  } else if (opts[['dmp']]) {
    add_pkgs <- 'parallel'
  } else if (opts[['idat2m']] | opts[['idat2beta']]) {
    add_pkgs <- c('IlluminaHumanMethylationEPICmanifest', 'minfi')
  } else {
    add_pkgs <- NULL
  }
  suppressMessages(lapply(c('tidyverse', add_pkgs),
                          library, character.only = TRUE))
}

main <- function(opts) {
  options(warn = 1)
  if (opts[['-v']]) print(opts)

  if (opts[['--session']]) {
    library('devtools', quietly = TRUE)
    print(devtools::session_info(pkgs = c('changepoint', 'GenomicRanges',
                                          'gplots', 'ggpubr', 'minfi',
                                          'RColorBrewer', 'tidyverse')))
  } else {
    load_packages(opts = opts)
    make_dir(path = opts[['--out']])
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
      prepare_site_csv(dst_dir = opts[['--out']],
                       platform = opts[['--platform']],
                       unfilter = opts[['--unfilter']])
    } else if (opts[['dmr']]) {
      segment_sites(site_csv = opts[['<site_csv>']],
                    mv_csv = opts[['<mv_csv>']], dst_dir = opts[['--out']],
                    sd_cutoff = opts[['--sd-cutoff']])
    } else if (opts[['cluster']]) {
      cluster_segments(dmrmv_csv = opts[['<dmrmv_csv>']],
                       dst_dir = opts[['--out']], k = opts[['--k']],
                       dist_method = opts[['--dist']],
                       hclust_method = opts[['--hclust']],
                       width = aspect_ratio[1], height = aspect_ratio[2])
    } else if (opts[['plot']]) {
      visualize_segments(site_csv = opts[['<site_csv>']],
                         mv_csv = opts[['<mv_csv>']],
                         seg_csv = opts[['<seg_csv>']],
                         dmrmv_csv = opts[['<dmrmv_csv>']],
                         dst_dir = opts[['--out']], width = aspect_ratio[1],
                         height = aspect_ratio[2])
    } else if (opts[['dmp']]) {
      write_dmp_mv(mv_csv = opts[['<mv_csv>']],
                   site_csv = opts[['<site_csv>']], dst_dir = opts[['--out']],
                   sd_cutoff = opts[['--sd-cutoff']])
    } else if (opts[['idat2m']]) {
      convert_idat_to_mv(idat_dir = opts[['<idat_dir>']],
                         out_csv = opts[['<out_csv>']],
                         offset = opts[['--offset']])
    } else if (opts[['idat2beta']]) {
      convert_idat_to_bv(idat_dir = opts[['<idat_dir>']],
                         out_csv = opts[['<out_csv>']],
                         offset = opts[['--offset']])
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


### segmer dmr

segment_sites <- function(site_csv, mv_csv, dst_dir, sd_cutoff = NULL,
                          n_cpu = detectCores(), cpt_method = 'PELT') {
  df_site <- read_site_csv(site_csv)
  df_mv <- filter(read_met_csv(path = mv_csv), name %in% df_site$name)
  cl <- makeCluster(n_cpu)
  df_p <- mutate(arrange(inner_join(df_site,
                                    shapiro_wilk_test(cl, df_mv),
                                    by = 'name'),
                         chrom, chromStart, chromEnd),
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
                       str_c(sub('.csv(|.gz|.bz2)$', '', basename(mv_csv)),
                             sub('(|.bed).csv(|.gz|.bz2)$', '',
                                 basename(site_csv)),
                             'seg.csv', sep = '.'))
  message('>>> Write a segment CSV:\t', seg_csv)
  write_csv(df_seg, path = seg_csv)

  message('>>> Calculate segmental median values')
  df_segmv <- summarize_all(group_by(select(left_join(select(df_seg,
                                                             name, segment),
                                                      df_mv, by = 'name'),
                                            -name),
                                     segment),
                            median,
                            .groups = 'drop')
  out_prefix <- file.path(dst_dir,
                          str_c(sub('.csv(|.gz|.bz2)$', '', basename(seg_csv)),
                                'mv', sep = '.'))
  segmv_csv <- str_c(out_prefix, '.all.csv')
  message('>>> Write a segmental value CSV:\t', segmv_csv)
  write_csv(df_segmv, path = segmv_csv)

  df_dmrmv <- filter(df_segmv,
                     segment %in% filter_by_sd(df_segmv,
                                               sd_cutoff = sd_cutoff)$segment)
  dmrmv_csv <- str_c(out_prefix, '.dmr.csv')
  message('>>> Write a DMR segmental value CSV:\t', dmrmv_csv)
  write_csv(df_dmrmv, path = dmrmv_csv)
  message('DMR segments:\t', nrow(df_dmrmv), ' / ', nrow(df_segmv))
}

read_site_csv <- function(path) {
  if (any(endsWith(path, c('.bed', '.bed.gz', '.bed.bz2')))) {
    return(read_bed(path))
  } else {
    return(read_csv_quietly(path))
  }
}

shapiro_wilk_test <- function(cl, df_mv) {
  message('>>> Perform Shapiro-Wilk normality tests')
  return(mutate(select(df_mv, 1),
                pval = parApply(cl, select(df_mv, -1), 1,
                                function(v) return(shapiro.test(v)$p.value))))
}

read_met_csv <- function(path, dropna = TRUE) {
  df_v <- read_csv_quietly(path)
  if (dropna) df_v <- drop_na(df_v)
  stopifnot(nrow(df_v) > 1, ncol(df_v) > 4)
  return(mutate(dplyr::rename(df_v, name = 1), name = as.character(name)))
}

df2num <- function(df) {
  return(as.vector(as.matrix(select_if(df, is.numeric))))
}

determine_boundary <- function(values, ...) {
  message('>>> Determine the boundary using K-means')
  df_b <- summarize(group_by(tibble(v = values,
                                    label = kmeans(values,
                                                   centers = 2,
                                                   ...)$cluster),
                             label),
                    min_v = min(v), max_v = max(v), .groups = 'drop')
  return(mean(max(df_b$min_v), min(df_b$max_v)))
}

mv2var <- function(df_mv) {
  return(mutate(select(df_mv, 1), variance = apply(select(df_mv, -1), 1, var)))
}

filter_by_sd <- function(df_v, sd_cutoff = NULL, use_kmeans = TRUE) {
  df_var <- mv2var(df_v)
  if (! is.null(sd_cutoff)) {
    message('>>> Filter sites by the specified SD threshold')
    sd_co <- as.numeric(sd_cutoff)
    var_co <- sd_co ^ 2
  } else if (use_kmeans) {
    message('>>> Filter sites by the variance threshold from K-means')
    var_co <- determine_boundary(df_var$variance)
    sd_co <- sqrt(var_co)
  } else {
    message('>>> Filter sites by the variance of the whole data')
    var_co <- var(df2num(df_v))
    sd_co <- sqrt(var_co)
  }
  message('variance (SD) cutoff:\t', var_co, '\t(', sd_co, ')')
  return(filter(df_v, df_var$variance > var_co))
}


### segmer cluster

cluster_segments <- function(dmrmv_csv, dst_dir, k = 3,
                             dist_method = 'euclidean',
                             hclust_method = 'ward.D2', width = 13,
                             height = 8) {
  stopifnot(k > 1)
  distfun <- function(...) return(dist(..., method = dist_method))
  hclustfun <- function(...) return(hclust(..., method = hclust_method))
  out_prefix <- file.path(dst_dir,
                          str_c(sub('.csv(|.gz|.bz2)$', '',
                                    basename(dmrmv_csv)),
                                dist_method, hclust_method,
                                str_c('k', k), sep = '.'))
  df_dmrmv <- read_csv_quietly(dmrmv_csv)
  hclust_csv <- str_c(out_prefix, '.hclust.csv')
  message('>>> Write an observed cluster CSV:\t', hclust_csv)
  mt_dmr <- as.matrix(column_to_rownames(df_dmrmv,
                                         var = colnames(df_dmrmv)[1]))
  hclust_hclusts <- stats::cutree(hclustfun(distfun(t(mt_dmr))), k = k)
  write_csv(tibble(sample_id = names(hclust_hclusts),
                   observed_cluster = hclust_hclusts),
            path = hclust_csv)
  heatmap_pdf <- str_c(out_prefix, '.heatmap.pdf')
  message('>>> Plot a heatmap:\t', heatmap_pdf)
  to_pdf(heatmap_plot(mt = mt_dmr, col_labels = hclust_hclusts,
                      distfun = distfun, hclustfun = hclustfun,
                      key.title = NA,
                      key.xlab = ifelse(str_detect(dmrmv_csv, '.mv.'),
                                        'M-value', 'value'),
                      xlab = str_c('sample  (total: ', ncol(mt_dmr), ')'),
                      ylab = str_c(ifelse(str_detect(dmrmv_csv, '.seg.mv.'),
                                          ' segment', ' site'),
                                   '  (total: ', nrow(mt_dmr), ')'),
                      margins = c(max(nchar(colnames(mt_dmr))) * 0.3 + 4,
                                  max(nchar(rownames(mt_dmr))) * 0.2 + 4)),
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

visualize_segments <- function(site_csv, mv_csv, seg_csv, dmrmv_csv, dst_dir,
                               n_cpu = detectCores(), width = 13, height = 8) {
  df_site <- read_site_csv(site_csv)
  df_mv_raw <- read_met_csv(mv_csv, dropna = FALSE)
  df_mv <- drop_na(df_mv_raw)
  out_prefix <- file.path(dst_dir,
                          str_c(sub('.csv(|.gz|.bz2)$', '', basename(mv_csv)),
                                sub('(|.bed).csv(|.gz|.bz2)$', '',
                                    basename(site_csv)),
                                sep = '.'))

  hist_pdf <- str_c(out_prefix, '.hist.pdf')
  message('>>> Plot an M-value histogram:\t', hist_pdf)
  to_pdf(plot(gghistogram(rbind(tibble(probe = 'unmasked',
                                       mv = df2num(filter(df_mv,
                                                          name %in%
                                                            df_site$name))),
                                tibble(probe = 'masked',
                                       mv = df2num(filter(df_mv,
                                                          ! name %in%
                                                            df_site$name)))),
                          x = 'mv', fill = 'probe', position = 'stack',
                          color = NA, alpha = 0.6,
                          palette = get_palette('nejm', 2), xlab = 'M-value',
                          ylab = str_c('site count  (total: ', nrow(df_mv),
                                       ')'),
                          bins = 30, legend = 'right')),
         path = hist_pdf, w = width, h = height)

  cl <- makeCluster(n_cpu)
  df_p <- inner_join(df_site,
                     mutate(shapiro_wilk_test(cl,
                                              filter(df_mv,
                                                     name %in% df_site$name)),
                            mlogp = -log10(pval)),
                     by = 'name')
  stopCluster(cl)
  lapply(unique(df_site$chrom),
         function(n, ylim) {
           p <- str_c(out_prefix, 'shapiro_wilk', n, 'pdf', sep = '.')
           message('>>> Plot p-values from Shapiro-Wilk normality tests:\t', p)
           d <- filter(df_p, chrom == n)
           to_pdf(plot(ggscatter(d, x = 'chromStart', y = 'mlogp',
                                 color = 'navy', alpha = 0.4, size = 0.1,
                                 ylim = ylim,
                                 xlab = str_c('position on ', n, '  (',
                                              nrow(d), ' sites)'),
                                 ylab = '-log10(p-value)') +
                       scale_x_continuous(labels = format_digit)),
                  path = p, w = width, h = height)
         },
         ylim = c(0, max(df_p$mlogp) * 1.1))

  df_seg <- read_csv_quietly(path = seg_csv)
  hist_pdf <- file.path(dst_dir, sub('.csv(|.gz|.bz2)$', '.hist.pdf',
                                     basename(seg_csv)))
  message('>>> Plot segmental histograms:\t', hist_pdf)
  segment2len <- function(s) return(1 - eval(parse(text = sub('^.*:', '', s))))
  ylab <-  str_c('segment count  (total: ', n_distinct(df_seg$segment), ')')
  g1 <- list(gghistogram(summarize(group_by(df_seg, segment),
                                   n = n_distinct(name),
                                   .groups = 'drop'),
                         x = 'n', fill = 'navy', color = NA, alpha = 0.4,
                         xscale = 'log10',
                         xlab = str_c('site count per segment  (total: ',
                                      n_distinct(df_seg$name), ')'),
                         ylab = ylab, bins = 30) +
             theme(plot.margin = margin(1, 2.5, 1, 1, 'lines')),
             gghistogram(mutate(distinct(df_seg, segment),
                                len = sapply(segment, segment2len)),
                         x = 'len', fill = 'navy', color = NA, alpha = 0.4,
                         xlab = 'segment bp length', ylab = ylab, bins = 30) +
             scale_x_log10(labels = format_digit) +
             theme(plot.margin = margin(1, 2.5, 1, 1, 'lines')))
  to_pdf(plot(ggarrange(plotlist = g1, nrow = 2)),
         path = hist_pdf, w = width, h = height)

  feature_csv <- file.path(dst_dir,
                           sub('.csv(|.gz|.bz2)$', '.feature.csv',
                               basename(dmrmv_csv)))
  message('>>> Write feature counts:\t', feature_csv)
  df_feature <- tibble(feature = c('differentially methylated segments',
                                   'determined segments',
                                   'available filtered sites',
                                   'filtered sites', 'probe sites'),
                       count = c(nrow(read_csv_quietly(dmrmv_csv)),
                                 length(unique(df_seg$segment)), nrow(df_seg),
                                 sum(df_site$name %in% df_mv_raw$name),
                                 nrow(df_mv_raw)))
  write_csv(df_feature, path = feature_csv)

  feature_pdf <- sub('.csv$', '.pdf', feature_csv)
  message('>>> Plot feature counts:\t', feature_csv)
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

write_dmp_mv <- function(mv_csv, site_csv, dst_dir, sd_cutoff = NULL) {
  df_site <- read_site_csv(site_csv)
  df_mv <- filter(read_met_csv(path = mv_csv), name %in% df_site$name)
  df_dmpmv <- filter(df_mv,
                     name %in% filter_by_sd(df_mv,
                                            sd_cutoff = sd_cutoff,
                                            use_kmeans = FALSE)$name)
  message('DMP sites:\t', nrow(df_dmpmv), ' / ', nrow(df_mv))
  dmp_csv <- file.path(dst_dir,
                       str_c(sub('.csv(|.gz|.bz2)$', '', basename(mv_csv)),
                             sub('(|.bed).csv(|.gz|.bz2)$', '',
                                 basename(site_csv)),
                             'mv.dmp.csv', sep = '.'))
  message('>>> Write DMP value CSV:\t', dmp_csv)
  write_csv(df_dmpmv, path = dmp_csv)
}


### idat2m

convert_idat_to_mv <- function(idat_dir, out_csv, offset = 1) {
  message('>>> Read IDAT files and calculate M-values:\t', idat_dir)
  df_mv <- NULL
  for (i in sort(fetch_idat_names(idat_dir = idat_dir))) {
    message(i)
    d <- calculate_mv_from_idat(i, offset = as.numeric(offset))
    if (is.null(df_mv)) {
      df_mv <- d
    } else {
      df_mv <- left_join(df_mv, d, by = 'name')
    }
  }
  message('>>> Write M-values:\t', out_csv)
  write_csv(arrange(df_mv, name), path = out_csv)
}

fetch_idat_names <- function(idat_dir) {
  return(do.call(intersect,
                 lapply(c('_Grn.idat$', '_Red.idat$'),
                        function(s) {
                          return(str_replace(list.files(idat_dir,
                                                        pattern = s,
                                                        full.names = TRUE),
                                             s, ''))
                        })))
}

calculate_mv_from_idat <- function(idat_name, offset = 1) {
  rg <- preprocessIllumina(read.metharray(idat_name))
  return(rownames_to_column(as.data.frame(log2((getMeth(rg) + offset) /
                                               (getUnmeth(rg) + offset))),
                            var = 'name'))
}

### idat2beta

convert_idat_to_bv <- function(idat_dir, out_csv, offset = 100) {
  message('>>> Read IDAT files and calculate beta-values:\t', idat_dir)
  df_bv <- NULL
  for (i in sort(fetch_idat_names(idat_dir = idat_dir))) {
    message(i)
    d <- calculate_bv_from_idat(i, offset = as.numeric(offset))
    if (is.null(df_bv)) {
      df_bv <- d
    } else {
      df_bv <- left_join(df_bv, d, by = 'name')
    }
  }
  message('>>> Write beta-values:\t', out_csv)
  write_csv(arrange(df_bv, name), path = out_csv)
}

calculate_bv_from_idat <- function(idat_name, offset = 1) {
  return(rownames_to_column(as.data.frame(getBeta(read.metharray(idat_name),
                                                  offset = offset)),
                            var = 'name'))
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
                                  gsub('\\]\n +<', '\\] <', doc)),
                             version = command_version))
}
