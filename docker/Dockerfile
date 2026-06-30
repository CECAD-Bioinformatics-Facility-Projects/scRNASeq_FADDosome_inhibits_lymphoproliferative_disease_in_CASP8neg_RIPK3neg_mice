# syntax=docker/dockerfile:1
###########################################################################################################################################################
# Extending the following:
# - Rocker RStudio container using new versioned paradigm: https://github.com/rocker-org/rocker-versioned2/blob/master/dockerfiles/rstudio_4.4.1.Dockerfile
# - license: GPLV2
###########################################################################################################################################################
FROM rocker/rstudio:4.4.1

### install key dependencies of certain packages that could be installed later
RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
    libcurl4-openssl-dev libfontconfig1-dev \
    libfreetype6-dev libicu-dev libx11-dev libxml2-dev

# install cmake and NLopt dev library so nloptr can build
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      cmake \
      libnlopt-dev && \
    rm -rf /var/lib/apt/lists/*
# install harfbuzz, fribidi, freetype2 & pkg-config so textshaping can build
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      pkg-config \
      libharfbuzz-dev \
      libfribidi-dev \
      libfreetype6-dev \
      libtiff5-dev \
      libjpeg-dev \
      libudunits2-dev \
      libpng-dev \
      libgdal-dev \
      libgeos-dev \
      libproj-dev \
      libmagick++-dev \
      libsqlite3-dev && \
    rm -rf /var/lib/apt/lists/*

# install R packages renv
RUN /rocker_scripts/bin/install2.r remotes renv 

RUN mkdir -p /home/rstudio/.cache/R
RUN mkdir -p /home/rstudio/project

RUN chown -R rstudio:rstudio /home/rstudio/.cache
RUN chown -R rstudio:rstudio /home/rstudio/project

USER rstudio

WORKDIR /home/rstudio/project

#RUN R -e "renv::install('tidyverse',dependencies=TRUE)"
COPY renv.lock /home/rstudio/project/
COPY renv/ /home/rstudio/project/renv/

# Make the image fully self-contained and safe to distribute, even when the whole
# project directory is bind-mounted at runtime:
#   * RENV_PATHS_LIBRARY puts the library ROOT OUTSIDE /home/rstudio/project so a
#     runtime bind mount of "." cannot shadow the installed packages. renv appends
#     a platform prefix (linux-ubuntu-jammy/R-4.4/x86_64-pc-linux-gnu) under it.
#   * RENV_CONFIG_CACHE_ENABLED=FALSE installs packages as REAL copies into that
#     library instead of symlinking them into ~/.cache/R/renv.
#   * NO BuildKit `type=cache` mount is used here: a cache mount is discarded from
#     the final image, so packages installed "through" it would not actually ship.
#     Installing into a normal layer bakes the real files into the image.
# The restore installs flat into RENV_PATHS_LIBRARY; a follow-up `mv` relocates
# the packages into the platform-prefixed subdirectory renv expects at runtime.
ENV RENV_PATHS_LIBRARY=/home/rstudio/renv-library
ENV RENV_CONFIG_CACHE_ENABLED=FALSE
RUN mkdir -p /home/rstudio/renv-library
RUN R -e "options(timeout = 1200); renv::restore(library = '/home/rstudio/renv-library', prompt = FALSE)"

# Relocate the flat library into the platform-prefixed path renv computes at runtime
# (RENV_PATHS_LIBRARY/<platform>/R-<ver>/<arch>). Without this, renv's autoloader
# sees an empty platform dir and packages fail to load.
RUN PLATFORM_DIR="/home/rstudio/renv-library/linux-ubuntu-jammy/R-4.4/x86_64-pc-linux-gnu" \
    && mkdir -p "$PLATFORM_DIR" \
    && find /home/rstudio/renv-library -maxdepth 1 -mindepth 1 -type d ! -name "linux-*" -exec mv -t "$PLATFORM_DIR" {} +

RUN echo 'setwd("/home/rstudio/project")\ngetwd()\n' > /home/rstudio/.Rprofile
RUN echo 'setHook("rstudio.sessionInit", function(newSession) { if (newSession && is.null(rstudioapi::getActiveProject())) rstudioapi::filesPaneNavigate("/home/rstudio/project") }, action = "append")' >> /home/rstudio/.Rprofile

RUN chown -R rstudio:rstudio /home/rstudio/.Rprofile

USER root

# Persist the external renv library location for every R/RStudio session so that
# renv's autoloader sets .libPaths() to the baked-in library (which lives outside
# the bind-mounted project directory) regardless of what the host's renv/ contains.
RUN echo 'RENV_PATHS_LIBRARY=/home/rstudio/renv-library' >> "$(R RHOME)/etc/Renviron.site"
RUN chown -R rstudio:rstudio /home/rstudio/renv-library

#RUN apt-get update && \
#    apt-get install -y --no-install-recommends \
#      cmake \
#      libnlopt-dev && \
#    rm -rf /var/lib/apt/lists/*
# install harfbuzz, fribidi, freetype2 & pkg-config so textshaping can build
#RUN apt-get update && \
#    apt-get install -y --no-install-recommends \
#      pkg-config \
#      libharfbuzz-dev \
#      libfribidi-dev \
#      libfreetype6-dev \
#      libtiff5-dev \
#      libjpeg-dev \
#      libudunits2-dev \
#      libpng-dev \
#      libgdal-dev \
#      libgeos-dev \
#      libproj-dev \
#      libmagick++-dev \
#      libsqlite3-dev && \
#    rm -rf /var/lib/apt/lists/*
