FROM ubuntu:18.04

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
ENV PACKAGE_DEV="gfortran gcc-multilib autoconf libcurl4-openssl-dev libmariadb-client-lgpl-dev liblzma-dev libncurses5-dev libblas-dev liblapack-dev libssh2-1-dev libxml2-dev vim libssl-dev libcairo2-dev libbz2-dev"
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		nano ed locales vim-tiny fonts-texgyre \
    $PACKAGE_DEV $PACKAGE_BIO \
	&& rm -rf /var/lib/apt/lists/*

RUN apt-get autoremove
RUN apt-get update && apt-get install -y --no-install-recommends iputils-ping gcc libgit2-dev

## oncoEnrichR dependencies: restore all R package versions and dependencies with renv package
## PINNED versions of all packages can be found in /renv.lock (JSON format)
RUN R -e "install.packages(repos = NULL, 'https://cloud.r-project.org/src/contrib/renv_0.13.2.tar.gz', lib = '/usr/lib/R/library/')"
COPY renv.lock /
COPY oncoEnrichR_0.9.6.tar.gz /

RUN R -e "library(renv); renv::restore(prompt = F, library = '/usr/lib/R/library/')"

WORKDIR /

## PINNED version - oncoEnrichR
RUN R -e "install.packages(repos = NULL, 'oncoEnrichR_0.9.6.tar.gz', lib = '/usr/lib/R/library/')"

## PINNED VERSION: pandoc (for HTML report generation)
RUN wget https://github.com/jgm/pandoc/releases/download/2.9.2/pandoc-2.9.2-1-amd64.deb && \
  dpkg -i pandoc* && \
  rm pandoc* && \
  apt-get clean
