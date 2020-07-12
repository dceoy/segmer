FROM dceoy/r-jupyter:latest

ADD https://raw.githubusercontent.com/dceoy/clir/master/install_clir.sh /tmp/install_clir.sh
ADD segmet.R /usr/local/bin/segmet.R

RUN set -e \
      && ln -s segmet.R /usr/local/bin/segmet

RUN set -e \
      && apt-get -y update \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        gnupg \
      && echo 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
        > /etc/apt/sources.list.d/r.list \
      && apt-key adv --keyserver keyserver.ubuntu.com \
        --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        libbz2-dev libjpeg-turbo8-dev libpng-dev \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -e \
      && bash /tmp/install_clir.sh --root --force --delete-r-lib \
      && rm -f /tmp/install_clir.sh

RUN set -e \
      && clir update \
      && clir install --devt=cran changepoint gplots RColorBrewer tidyverse \
      && clir install --devt=bioc GenomicRanges \
      && clir validate changepoint GenomicRanges gplots RColorBrewer tidyverse

RUN set -e \
      && clir install --devt=github IRkernel/IRkernel \
      && R -q -e 'IRkernel::installspec();'

ENTRYPOINT ["/usr/local/bin/segmet"]
