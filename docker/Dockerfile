FROM ubuntu:18.04

RUN apt-get update && apt-get -y install apache2 apt-utils build-essential curl git unzip vim wget sudo libudunits2-dev libgeos-dev libgdal-dev
RUN apt-get update && apt-get install apt-transport-https

#RUN sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
#RUN sudo apt-get update \
#   && sudo apt-get -y install software-properties-common

ENV TZ=Europe/Minsk
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN sudo apt update -qq
# install two helper packages we need
RUN sudo apt-get update && sudo apt-get -y install --no-install-recommends software-properties-common dirmngr
# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed

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
RUN apt-get update && apt-get install -y --no-install-recommends iputils-ping gcc libgit2-dev libharfbuzz-dev librdf0-dev libsodium-dev libfribidi-dev libxt-dev


RUN wget https://github.com/Kitware/CMake/releases/download/v3.20.0/cmake-3.20.0.tar.gz
RUN tar -zxvf cmake-3.20.0.tar.gz
RUN cd cmake-3.20.0 && ./bootstrap && make && make install

## oncoEnrichR dependencies: restore all R package versions and dependencies with renv package
## PINNED versions of all packages can be found in /renv.lock (JSON format)
RUN R -e "install.packages(repos = NULL, 'https://cloud.r-project.org/src/contrib/renv_1.0.7.tar.gz', lib = '/usr/lib/R/library/')"
COPY renv.lock /
COPY oncoEnrichR_1.5.1.tar.gz /

#RUN R -e "library(renv)"

RUN R -e "library(renv); renv::restore(prompt = F, library = '/usr/lib/R/library/')"
WORKDIR /

## PINNED version - oncoEnrichR
RUN R -e "install.packages(repos = NULL, 'oncoEnrichR_1.5.1.tar.gz', lib = '/usr/lib/R/library/')"

## PINNED VERSION: pandoc (for HTML report generation)
RUN wget https://github.com/jgm/pandoc/releases/download/2.13/pandoc-2.13-1-amd64.deb && \
  dpkg -i pandoc* && \
  rm pandoc* && \
  apt-get clean


RUN set -e -x && \
    apt-get update && apt-get install -y --no-install-recommends \
    gdebi-core \
    librsvg2-bin && apt-get clean

RUN curl -LO https://quarto.org/download/latest/quarto-linux-amd64.deb
RUN gdebi --non-interactive quarto-linux-amd64.deb
RUN rm -f quarto-linux-amd64.deb

RUN rm -rf /root/.cache
RUN rm -rf /oncoEnrichR_1.5.1.tar.gz
RUN rm -rf /cmake-3.20.0*
