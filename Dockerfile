FROM rocker/rstudio

ENV CRAN_MIRROR http://cran.rstudio.com

RUN apt-get update --fix-missing \
	&& apt-get install -y \
		ca-certificates \
    	        libglib2.0-0 \
	 	libxext6 \
	   	libsm6  \
	   	libxrender1 \
		libxml2-dev

# python3

RUN apt-get install -y \
		python3-pip \
		python3-dev \
	        && pip3 install virtualenv

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.4.10-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH /opt/conda/bin:$PATH


# Python packages from conda
RUN conda install -y numpy
RUN pip install --pre vaex

RUN apt-get install -y cmake

RUN apt-get install -y zlib1g-dev

RUN git clone https://github.com/mattgodbolt/zindex.git \
      && cd zindex \
      && make \
      && ln -s /zindex/build/Release/zq /usr/bin/zq \
      && ln -s /zindex/build/Release/zindex /usr/bin/zindex 

#Install gsutil
RUN apt-get install -y apt-transport-https \
                       curl 
ENV CLOUD_SDK_VERSION 193.0.0
RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk=${CLOUD_SDK_VERSION}-0 $INSTALL_COMPONENTS && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    gcloud --version
VOLUME ["/root/.config"]

#Install fuse
RUN export GCSFUSE_REPO="gcsfuse-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" > /etc/apt/sources.list.d/gcsfuse.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y gcsfuse

#Cairo
RUN apt-get install -y libgtk2.0-dev \
                       libcairo2-dev \
                       xvfb \
                       xauth \
                       xfonts-base \
                       libxt-dev \
#                       xdg-utils \
                       vim

#Install fonts
#RUN  mkdir ~/.fonts
COPY FONTS /usr/share/fonts/
RUN fc-cache -f -v /usr/share/fonts/
#Install R and R packages
RUN install2.r --repos ${CRAN_MIRROR}\
		Rcpp \
		devtools \
		roxygen2 \
		extrafont\
                Cairo \
                svglite \
                knitr \
		rmarkdown \
		yaml \
		reticulate 

RUN Rscript -e "library(extrafont);font_import('/usr/share/fonts',prompt = FALSE)"

#Install R dependencies
RUN install2.r --repos ${CRAN_MIRROR}\
                ggplot2 \
		dplyr \
                Matrix \
                colorRamps\
                RColorBrewer\
                dichromat\
                viridis\
                colourpicker

RUN install2.r --repos ${CRAN_MIRROR}\
                shiny\
                shinydashboard\
                shinythemes\
                scales\
                sp\
                plyr\
                ggrepel

COPY scSVA_0.1.0.tar.gz /home/

RUN install2.r --repos ${CRAN_MIRROR}\
                plotly\
                DT\
                data.table\
                tableHTML\
                && Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('rhdf5');  biocLite('BiocParallel')" \
                && Rscript -e "devtools::install_github(c('thomasp85/shinyFiles'))" \
                #&& Rscript -e "install.packages('/home/scSVA_0.1.0.tar.gz', repos = NULL, type='source')"\
                ## clean up
                && rm -rf /tmp/downloaded_packages/ /tmp/*.rds
#Install scSVA                 
RUN Rscript -e "install.packages('/home/scSVA_0.1.0.tar.gz', repos = NULL, type='source')" \
               && rm /home/scSVA_0.1.0.tar.gz

# Run rocker/rstudio init
CMD ["/init"]
