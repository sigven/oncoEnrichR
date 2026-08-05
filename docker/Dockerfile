FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Minsk

RUN apt-get update && \
    apt-get install -y tzdata && \
    ln -fs /usr/share/zoneinfo/$TZ /etc/localtime && \
    dpkg-reconfigure --frontend noninteractive tzdata

RUN apt-get update && apt-get -y install apache2 apt-utils build-essential curl git unzip vim wget sudo libudunits2-dev libgeos-dev libgdal-dev
RUN apt-get update && apt-get install apt-transport-https

RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common dirmngr

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
    | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

RUN add-apt-repository \
   "deb [arch=amd64] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/"

RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base r-base-core r-recommended

USER root
WORKDIR /

ENV PACKAGE_BIO="libhts-dev"
ENV PACKAGE_DEV="gfortran gcc-multilib autoconf libmariadb-dev liblzma-dev libncurses5-dev libblas-dev liblapack-dev libssh2-1-dev libxml2-dev vim libssl-dev libcairo2-dev libbz2-dev"
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		nano ed locales vim-tiny fonts-texgyre \
    $PACKAGE_DEV $PACKAGE_BIO \
	&& rm -rf /var/lib/apt/lists/*

#RUN apt-get autoremove
RUN apt-get update && apt-get install -y --no-install-recommends iputils-ping gcc libgit2-dev libharfbuzz-dev librdf0-dev libsodium-dev libfribidi-dev libxt-dev

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        cmake \
        libcurl4-openssl-dev \
        libglpk-dev \
        pandoc

## oncoEnrichR dependencies: restore all R package versions and dependencies with renv package
## PINNED versions of all packages can be found in /renv.lock (JSON format)
RUN R -e "install.packages(repos = NULL, 'https://cloud.r-project.org/src/contrib/renv_1.1.5.tar.gz', lib = '/usr/lib/R/library/')"
COPY renv.lock /
COPY oncoEnrichR_1.6.2.tar.gz /

#RUN R -e "library(renv)"

RUN R -e "library(renv); renv::restore(prompt = F, library = '/usr/lib/R/library/')"
WORKDIR /

## PINNED version - oncoEnrichR
RUN R -e "install.packages(repos = NULL, 'oncoEnrichR_1.6.2.tar.gz', lib = '/usr/lib/R/library/')"

RUN set -e -x && \
    apt-get update && apt-get install -y --no-install-recommends \
    gdebi-core \
    librsvg2-bin && apt-get clean

RUN curl -LO https://quarto.org/download/latest/quarto-linux-amd64.deb
RUN gdebi --non-interactive quarto-linux-amd64.deb
RUN rm -f quarto-linux-amd64.deb

RUN rm -rf /root/.cache
RUN rm -rf /oncoEnrichR_1.6.2.tar.gz
