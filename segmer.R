#!/usr/bin/env Rscript

'Segment Detector for differentially methylated regions in methylome data

Usage:
  segmer bed [-v] [--platform=<str>] [--unfilter] [--out=<dir>]
  segmer dmr [-v] [--kmeans-k=<int>] [--sd-cutoff=<dbl>] [--seed=<int>]
             [--cpus=<int>] [--ar=<ratio>] [--out=<dir>] <site_csv> <mv_csv>
  segmer cluster [-v] [--cutree-k=<int>] [--hclust=<str>] [--dist=<str>]
                 [--ar=<ratio>] [--out=<dir>] <dmrmv_csv>
  segmer plot [-v] [--cpus=<int>] [--ar=<ratio>] [--out=<dir>] <site_csv>
              <mv_csv> <seg_csv> <segmv_csv> <dmrmv_csv>
  segmer dmp [-v] [--kmeans-k=<int>] [--sd-cutoff=<dbl>] [--seed=<int>]
             [--out=<dir>] <site_csv> <mv_csv>
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
  --kmeans-k=<int>    Specify the number of clusters for K-means [default: 4]
  --sd-cutoff=<dbl>   Use a fixed SD cutoff
  --seed=<int>        Specify a random seed
  --cpus=<int>        Limit the CPU cores to be used
  --cutree-k=<int>    Specify the number of sample clusters [default: 3]
  --hclust=<str>      Specify the method of stats::hclust [default: ward.D2]
  --dist=<str>        Specify the method of stats::dist [default: euclidean]
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
  <segmv_csv>         Path to a segmental M-value CSV file
                      (`*.seg.mv.all.csv` created by `segmer dmr`)
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
    make_dir(opts[['--out']])
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
    if (opts[['dmr']] | opts[['plot']]) {
      message('>>> Set the number of CPU cores')
      n_cpu <- ifelse(is.null(opts[['--cpus']]), parallel::detectCores(),
                      as.integer(opts[['--cpus']]))
      message(n_cpu)
    } else {
      n_cpu <- 1
    }

    if (opts[['bed']]) {
      prepare_site_csv(dst_dir = opts[['--out']],
                       platform = opts[['--platform']],
                       unfilter = opts[['--unfilter']])
    } else if (opts[['dmr']]) {
      segment_sites(site_csv = opts[['<site_csv>']],
                    mv_csv = opts[['<mv_csv>']], dst_dir = opts[['--out']],
                    kmeans_k = opts[['--kmeans-k']],
                    sd_cutoff = opts[['--sd-cutoff']], n_cpu = n_cpu)
    } else if (opts[['cluster']]) {
      cluster_segments(dmrmv_csv = opts[['<dmrmv_csv>']],
                       dst_dir = opts[['--out']],
                       cutree_k = opts[['--cutree-k']],
                       dist_method = opts[['--dist']],
                       hclust_method = opts[['--hclust']],
                       width = aspect_ratio[1], height = aspect_ratio[2])
    } else if (opts[['plot']]) {
      visualize_segments(site_csv = opts[['<site_csv>']],
                         mv_csv = opts[['<mv_csv>']],
                         seg_csv = opts[['<seg_csv>']],
                         segmv_csv = opts[['<segmv_csv>']],
                         dmrmv_csv = opts[['<dmrmv_csv>']],
                         dst_dir = opts[['--out']], width = aspect_ratio[1],
                         height = aspect_ratio[2], n_cpu = n_cpu)
    } else if (opts[['dmp']]) {
      write_dmp_mv(mv_csv = opts[['<mv_csv>']],
                   site_csv = opts[['<site_csv>']], dst_dir = opts[['--out']],
                   kmeans_k = opts[['--kmeans-k']],
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
                  csv_gz <- str_replace(dst, '.rds$', '.csv.gz')
                  message('>>> Convert RDS to CSV:\t', csv_gz)
                  d <- granges2tibble(readRDS(dst), var = 'name')
                  write_csv(d, path = csv_gz)
                  return(d)
                }))
}


### segmer dmr

segment_sites <- function(site_csv, mv_csv, dst_dir, kmeans_k = 2,
                          sd_cutoff = NULL, n_cpu = detectCores(),
                          cpt_method = 'PELT') {
  df_site <- read_site_csv(site_csv)
  df_mv <- filter(read_met_csv(mv_csv), name %in% df_site$name)
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
  message(nrow(df_seg), ' sites => ', n_distinct(df_seg$segment),
          ' segments')

  seg_csv <- file.path(dst_dir,
                       str_c(fetch_csv_prefix(mv_csv),
                             fetch_csv_prefix(site_csv, bed = TRUE),
                             'seg.csv', sep = '.'))
  message('>>> Write a segment CSV:\t', seg_csv)
  write_csv(df_seg, path = seg_csv)

  message('>>> Calculate segmental mean values')
  df_segmv <- summarize_all(group_by(select(left_join(select(df_seg,
                                                             name, segment),
                                                      df_mv, by = 'name'),
                                            -name),
                                     segment),
                            mean,
                            .groups = 'drop')
  out_prefix <- file.path(dst_dir,
                          fetch_csv_prefix(seg_csv, replacement = '.mv'))
  segmv_csv <- str_c(out_prefix, '.all.csv')
  message('>>> Write a segmental value CSV:\t', segmv_csv)
  write_csv(df_segmv, path = segmv_csv)

  df_dmrmv <- filter(df_segmv,
                     segment %in% filter_by_sd(df_segmv,
                                               kmeans_k = kmeans_k,
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
  d <- read_csv_quietly(path)
  if (dropna) d <- drop_na(d)
  stopifnot(nrow(d) > 1, ncol(d) > 4)
  return(mutate(dplyr::rename(d, name = 1), name = as.character(name)))
}

filter_by_sd <- function(df_v, kmeans_k = 2, sd_cutoff = NULL,
                         kmeans_iter_max = 100) {
  df_var <- row2var(df_v)
  if (is.null(sd_cutoff)) {
    message('>>> Determine the variance threshold using K-means (k = ',
            kmeans_k, ')')
    km <- kmeans(df_var$variance, centers = kmeans_k,
                 iter.max = kmeans_iter_max)
    df_k <- tibble(v = df_var$variance,
                   label = km$cluster,
                   is_highest = (km$cluster ==
                                 which(km$center == max(km$center))))
    var_co <- mean(min(filter(df_k, is_highest)$v),
                   max(filter(df_k, ! is_highest)$v))
    sd_co <- sqrt(var_co)
    message('>>> Filter sites by the variance threshold')
  } else {
    sd_co <- as.numeric(sd_cutoff)
    var_co <- sd_co ^ 2
    message('>>> Filter sites by the specified SD threshold')
  }
  message('threshold variance (SD):\t', var_co, '  (', sd_co, ')')
  return(filter(df_v, df_var$variance > var_co))
}

row2var <- function(df_v) {
  return(mutate(select(df_v, 1), variance = apply(select(df_v, -1), 1, var)))
}

### segmer cluster

cluster_segments <- function(dmrmv_csv, dst_dir, cutree_k = 3,
                             dist_method = 'euclidean',
                             hclust_method = 'ward.D2', width = 13,
                             height = 8) {
  stopifnot(cutree_k > 1)
  distfun <- function(...) return(dist(..., method = dist_method))
  hclustfun <- function(...) return(hclust(..., method = hclust_method))
  out_prefix <- file.path(dst_dir,
                          str_c(fetch_csv_prefix(dmrmv_csv), dist_method,
                                hclust_method, str_c('k', cutree_k),
                                sep = '.'))
  df_dmrmv <- read_csv_quietly(dmrmv_csv)
  hclust_csv <- str_c(out_prefix, '.hclust.csv')
  message('>>> Write an observed cluster CSV:\t', hclust_csv)
  mt_dmr <- as.matrix(column_to_rownames(df_dmrmv,
                                         var = colnames(df_dmrmv)[1]))
  hclust_hclusts <- stats::cutree(hclustfun(distfun(t(mt_dmr))), k = cutree_k)
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
                      margins = c(6, 10)),
         path = heatmap_pdf, w = width, h = height)
}

heatmap_plot <- function(mt, col_labels, col = rev(brewer.pal(9, 'RdBu')),
                         ...) {
  return(heatmap.2(mt,
                   ColSideColors = brewer.pal(n_distinct(col_labels),
                                              'Accent')[col_labels],
                   scale = 'none', trace = 'none', density.info = 'none',
                   col = col, ...))
}


### segmer plot

fetch_csv_prefix <- function(path, replacement = '', bed = FALSE) {
  return(str_replace(basename(path),
                     str_c(ifelse(bed, '(|.bed)', ''), '.csv(|.gz|.bz2)$'),
                     replacement))
}

visualize_segments <- function(site_csv, mv_csv, seg_csv, segmv_csv, dmrmv_csv,
                               dst_dir, n_cpu = detectCores(), width = 13,
                               height = 8, palette = 'nejm') {
  df_site <- read_site_csv(site_csv)
  df_mv_raw <- read_met_csv(mv_csv, dropna = FALSE)
  df_mv <- drop_na(df_mv_raw)
  out_prefix <- file.path(dst_dir,
                          str_c(fetch_csv_prefix(mv_csv),
                                fetch_csv_prefix(site_csv, bed = TRUE),
                                sep = '.'))

  mv_hist_pdf <- str_c(out_prefix, '.hist.pdf')
  message('>>> Plot an M-value histogram:\t', mv_hist_pdf)
  to_pdf(plot(gghistogram(rbind(tibble(probe = 'used',
                                       mv = df2num(filter(df_mv,
                                                          name %in%
                                                            df_site$name))),
                                tibble(probe = 'masked',
                                       mv = df2num(filter(df_mv,
                                                          ! name %in%
                                                            df_site$name)))),
                          x = 'mv', fill = 'probe', position = 'stack',
                          color = NA, alpha = 0.6,
                          palette = get_palette(palette, 2), xlab = 'M-value',
                          ylab = str_c('site count  (total: ', nrow(df_mv),
                                       ')'),
                          bins = 30, legend = 'right')),
         path = mv_hist_pdf, w = width, h = height)

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

  df_dmrmv <- read_csv_quietly(dmrmv_csv)
  df_ndmrmv <- filter(read_csv_quietly(segmv_csv),
                      ! segment %in% df_dmrmv$segment)
  df_seg <- mutate(read_csv_quietly(seg_csv),
                   DMR = (segment %in% df_dmrmv$segment))
  seg_hist_pdf <- file.path(dst_dir,
                            fetch_csv_prefix(dmrmv_csv,
                                             replacement = '.hist.pdf'))
  message('>>> Plot segmental histograms:\t', seg_hist_pdf)
  n_site_used <- n_distinct(df_seg$name)
  n_seg <- n_distinct(df_seg$segment)
  seg_ylab <-  str_c('segment count  (total: ', n_seg, ')')
  gs <- lapply(list(gghistogram(summarize(group_by(df_seg, segment, DMR),
                                          n = n_distinct(name),
                                          .groups = 'drop'),
                                x = 'n', fill = 'DMR', color = NA, alpha = 0.6,
                                palette = get_palette(palette, 2),
                                position = 'stack', xscale = 'log10',
                                title = 'Probe site',
                                xlab = str_c('site count per segment',
                                             '  (total: ', n_site_used, ')'),
                                ylab = seg_ylab, bins = 30, legend = 'right'),
                    gghistogram(rbind(tibble(v = df2num(df_dmrmv), DMR = TRUE),
                                      tibble(v = df2num(df_ndmrmv),
                                             DMR = FALSE)),
                                x = 'v', fill = 'DMR', color = NA, alpha = 0.6,
                                palette = get_palette(palette, 2),
                                position = 'stack', title = 'Segment M-value',
                                xlab = 'mean M-value per segment',
                                ylab = str_c('count  (', n_site_used,
                                             ' sites x ', ncol(df_dmrmv) - 1,
                                             ' samples)'),
                                bins = 30, legend = 'right'),
                    gghistogram(mutate(distinct(df_seg, segment, DMR),
                                       len = sapply(segment, segment2len)),
                                x = 'len', fill = 'DMR', color = NA,
                                alpha = 0.6, palette = get_palette(palette, 2),
                                position = 'stack', xscale = 'log10',
                                title = 'Segment length',
                                xlab = 'segment bp length', ylab = seg_ylab,
                                bins = 30, legend = 'right'),
                    gghistogram(rbind(mutate(row2var(df_dmrmv), DMR = TRUE),
                                      mutate(row2var(df_ndmrmv), DMR = FALSE)),
                                x = 'variance', fill = 'DMR', color = NA,
                                alpha = 0.6, palette = get_palette(palette, 2),
                                position = 'stack',
                                title = 'Variance of segment M-value',
                                xlab = 'M-value variance per segment',
                                ylab = seg_ylab, bins = 30, legend = 'right')),
               function(g) {
                 return(g + theme(plot.margin = margin(1, 1, 1, 1, 'lines')))
               })
  to_pdf(plot(ggarrange(plotlist = gs, common.legend = TRUE,
                        legend = 'right')),
         path = seg_hist_pdf, w = width, h = height)

  feature_csv <- file.path(dst_dir,
                           fetch_csv_prefix(dmrmv_csv,
                                            replacement = '.feature.csv'))
  message('>>> Write feature counts:\t', feature_csv)
  df_feature <- tibble(feature = c('differentially methylated segments',
                                   'determined segments',
                                   'available filtered sites',
                                   'filtered sites', 'probe sites'),
                       count = c(nrow(df_dmrmv), n_seg, n_site_used,
                                 sum(df_mv_raw$name %in% df_site$name),
                                 nrow(df_mv_raw)))
  write_csv(df_feature, path = feature_csv)

  feature_pdf <- str_replace(feature_csv, '.csv$', '.pdf')
  message('>>> Plot feature counts:\t', feature_csv)
  to_pdf(plot(ggbarplot(df_feature, x = 'feature', y = 'count', fill = 'navy',
                        color = 'navy', alpha = 0.4, orientation = 'horiz',
                        label = TRUE, lab.hjust = -0.2) +
              scale_y_continuous(limits = c(0, max(df_feature$count) * 1.1)) +
              labs(x = element_blank(), y = 'feature count') +
              theme(plot.margin = margin(1, 2.5, 1, 1, 'lines'))),
         path = feature_pdf, w = width, h = height)
}

segment2len <- function(s) {
  return(1 - eval(parse(text = str_replace(s, '^.*:', ''))))
}

df2num <- function(df) {
  return(as.vector(as.matrix(select_if(df, is.numeric))))
}

format_digit <- function(i) {
  return(format(i, scientific = FALSE))
}


### dmp

write_dmp_mv <- function(mv_csv, site_csv, dst_dir, kmeans_k = 4,
                         sd_cutoff = NULL) {
  df_site <- read_site_csv(site_csv)
  df_mv <- filter(read_met_csv(mv_csv), name %in% df_site$name)
  df_dmpmv <- filter(df_mv,
                     name %in% filter_by_sd(df_mv,
                                            kmeans_k = kmeans_k,
                                            sd_cutoff = sd_cutoff)$name)
  message('DMP sites:\t', nrow(df_dmpmv), ' / ', nrow(df_mv))
  dmp_csv <- file.path(dst_dir,
                       str_c(fetch_csv_prefix(mv_csv),
                             fetch_csv_prefix(site_csv, bed = TRUE),
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
  main(opts = docopt::docopt(gsub('\n {8,}([\\[<])', ' \\1', doc),
                             version = command_version))
}
