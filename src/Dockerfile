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

RUN sudo add-apt-repository 'deb [arch=amd64,i386] https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
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

RUN R -e "install.packages('devtools', dependencies = T)"
WORKDIR /
RUN R -e "devtools::install_github('sigven/oncoEnrichR')"
RUN R -e "install.packages('openxlsx')"

## Install pandoc (for HTML report generation)
RUN wget https://github.com/jgm/pandoc/releases/download/2.9.2.1/pandoc-2.9.2.1-1-amd64.deb && \
  dpkg -i pandoc* && \
  rm pandoc* && \
  apt-get clean

#WORKDIR /workdir/output
