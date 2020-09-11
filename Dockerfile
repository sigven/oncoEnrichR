FROM ubuntu:18.04

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

RUN apt-get update && apt-get -y install apache2 apt-utils build-essential curl git unzip vim wget sudo libudunits2-dev libgeos-dev libgdal-dev
RUN apt-get update && apt-get install apt-transport-https

RUN sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN sudo apt-get update \
   && sudo apt-get -y install software-properties-common

ENV TZ=Europe/Minsk
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN sudo add-apt-repository 'deb [arch=amd64,i386] https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
RUN apt-get update && apt-get -y install --no-install-recommends r-base-core r-recommended r-base

USER root
WORKDIR /

ENV PACKAGE_BIO="libhts2"
ENV PACKAGE_DEV="gfortran gcc-multilib autoconf libmariadb-client-lgpl-dev liblzma-dev libncurses5-dev libblas-dev liblapack-dev libssh2-1-dev libxml2-dev vim libssl-dev libcairo2-dev libbz2-dev libcurl4-openssl-dev"
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		nano ed locales vim-tiny fonts-texgyre \
    $PACKAGE_DEV $PACKAGE_BIO \
	&& rm -rf /var/lib/apt/lists/*

RUN apt-get autoremove

RUN apt-get update && apt-get install -y --no-install-recommends gcc
# RUN R -e "install.packages(c('clipr', 'cli', 'crayon', 'curl', 'desc', 'fs', 'gh', 'git2r', 'glue', 'purrr', 'rematch2', 'rlang', 'rprojroot', 'rstudioapi', 'whisker', 'withr', 'yaml'))"
# RUN R -e "install.packages(c('callr', 'covr', 'DT', 'memoise', 'pkgbuild', 'pkgload', 'rcmdcheck', 'remotes', 'roxygen2', 'rversions', 'sessioninfo', 'testthat'))"
RUN R -e "install.packages('devtools')"
WORKDIR /

# RUN R -e "install.packages(c('reshape2','plyr'), dependencies = 'Depends')"
# RUN R -e "install.packages(c('zip','openxlsx'), dependencies = 'Depends')"
# RUN R -e "install.packages(c('generics'), dependencies = 'Depends')"
# RUN R -e "install.packages(c('rematch','hms'), dependencies = 'Depends')"
# RUN R -e "install.packages(c('cellranger','progress'), dependencies = 'Depends')"
# RUN R -e "install.packages(c('tinytex'), dependencies = 'Depends')"
# RUN R -e "install.packages(c('rmarkdown'), dependencies = 'Depends')"
# RUN R -e "install.packages(c('selectr'), dependencies = 'Depends')"
# RUN R -e "install.packages(c('DBI','tidyselect','blob'), dependencies = 'Depends')"
# RUN R -e "install.packages(c('dbplyr', 'dplyr', 'forcats', 'haven', 'hms', 'lubridate', 'readr', 'readxl', 'reprex', 'rvest'), dependencies = 'Depends')"
# RUN R -e "devtools::install_version('cpp11', version = '0.1.0', repos = 'http://cran.us.r-project.org')"
# RUN R -e "install.packages('tidyr', dependencies = F)"
# RUN R -e "install.packages('broom', dependencies = 'Depends')"
# RUN R -e "install.packages('modelr', dependencies = 'Depends')"



# RUN R -e "install.packages('colorspace', dependencies = 'Depends')"
# RUN R -e "install.packages(c('farver', 'labeling', 'munsell', 'RColorBrewer', 'viridisLite'), dependencies = 'Depends')"
# RUN R -e "install.packages(c('gtable', 'isoband', 'scales'), dependencies = 'Depends')"
# RUN R -e "install.packages(c('dichromat', 'colorspace', 'mapproj', 'maps'), dependencies = 'Depends')"
# RUN R -e "install.packages(c('hexbin','data.table'), dependencies = 'Depends')"
# RUN R -e "BioCManager::"
# RUN R -e "install.packages(c('tidyverse', 'visNetwork', 'plotly', 'pals', 'org.Hs.eg.db', 'clusterProfiler', 'igraph'), dependencies = 'Depends')"

RUN R -e "devtools::install_github('sigven/oncoEnrichR', force = T)"

## Install pandoc (for HTML report generation)
RUN wget https://github.com/jgm/pandoc/releases/download/2.9.2.1/pandoc-2.9.2.1-1-amd64.deb && \
  dpkg -i pandoc* && \
  rm pandoc* && \
  apt-get clean

#WORKDIR /workdir/output
