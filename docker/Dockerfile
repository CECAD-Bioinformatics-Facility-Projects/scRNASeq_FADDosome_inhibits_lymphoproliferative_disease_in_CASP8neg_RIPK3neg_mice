###########################################################################################################################################################
# Extending the following:
# - Rocker RStudio container using new versioned paradigm: https://github.com/rocker-org/rocker-versioned2/blob/master/dockerfiles/rstudio_4.4.3.Dockerfile
# - license: GPLV2
###########################################################################################################################################################
FROM rocker/rstudio:4.4.3

### install key dependencies of certain packages that could be installed later
RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
    libcurl4-openssl-dev libfontconfig1-dev \
    libfreetype6-dev libicu-dev libx11-dev libxml2-dev


# install R packages renv
RUN /rocker_scripts/bin/install2.r remotes renv 

RUN mkdir -p /home/rstudio/.cache/R
RUN mkdir -p /home/rstudio/project

RUN chown -R rstudio:rstudio /home/rstudio/.cache
RUN chown -R rstudio:rstudio /home/rstudio/project

USER rstudio

RUN R -e "renv::install('tidyverse',dependencies=TRUE)"

RUN echo 'setwd("/home/rstudio/project")\ngetwd()\n' > /home/rstudio/.Rprofile
RUN echo 'setHook("rstudio.sessionInit", function(newSession) { if (newSession && is.null(rstudioapi::getActiveProject())) rstudioapi::filesPaneNavigate("/home/rstudio/project") }, action = "append")' >> /home/rstudio/.Rprofile

RUN chown -R rstudio:rstudio /home/rstudio/.Rprofile

USER root
