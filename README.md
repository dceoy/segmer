segmer
======

Differentially methylated segment detector for methylation array data

Supporetd platforms:

- Illumina
  - EPIC (Infinium MethylationEPIC BeadChip)
  - hm450 (Infinium Human Methylation 450K BeadChip)
  - hm27 (Infinium HumanMethylation 27K BeadChip)

Docker image
------------

Pull the image from [Docker Hub](https://hub.docker.com/r/dceoy/segmer/).

```sh
$ docker image pull dceoy/segmer
```

Installation
------------

1.  Install [clir](https://github.com/dceoy/clir).

2.  Install dependencies.

    ```sh
    $ clir install --devt=cran changepoint gplots ggpubr RColorBrewer tidyverse
    $ clir install --devt=bioc GenomicRanges IlluminaHumanMethylationEPICmanifest minfi
    $ clir validate \
        changepoint GenomicRanges gplots ggpubr \
        IlluminaHumanMethylationEPICmanifest minfi RColorBrewer tidyverse
    ```

3.  Check out segmer.

    ```sh
    $ git clone https://github.com/dceoy/segmer.git
    $ cp -a segmer/segmer.R /path/to/bin/
    ```

Usage
-----

1.  Calculate M-values from IDAT files.

    ```sh
    $ segmer idat2m \
        /path/to/epic/idat/ \
        epic_mv.csv.gz
    ```

2.  Download annotation data and write a target site CSV file.

    ```sh
    $ segmer bed \
        --platform=EPIC \
        --out=.
    ```

3.  Segment target sites.

    ```sh
    $ segmer dmr \
        --out=. \
        EPIC.hg38.f.bed.csv \
        epic_mv.csv.gz
    ```

4.  Execute clustering and draw the heatmap.

    ```sh
    $ segmer clust \
        --out=. \
        epic_mv.EPIC.hg38.f.seg.mv.dmr.csv
    ```

5.  Visualize the result of segmentation

    ```sh
    $ segmer stats \
        --out=. \
        EPIC.hg38.f.bed.csv \
        epic_mv.csv.gz \
        epic_mv.EPIC.hg38.f.seg.csv \
        epic_mv.EPIC.hg38.f.seg.mv.all.csv
    ```

Annotation Data
---------------

This tool uses the following annotation data:

- [InfiniumMethylation BeadChips Annotation](https://zwdzwd.github.io/InfiniumAnnotation)
