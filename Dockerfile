FROM rocker/rstudio

ENV CRAN_MIRROR http://cran.rstudio.com

RUN apt-get update --fix-missing \
	 && apt-get install -y \
	    ca-certificates \
    	    libglib2.0-0 \
	    libxext6 \
	    libsm6  \
	    libxrender1 \
	    libxml2-dev \
            fontocnfig \
            unzip

# Install Python
RUN apt-get install -y \
       	python3-pip \
	python3-dev \
        && pip3 install virtualenv

#Install numpy + vaex
RUN pip3 install --pre numpy
RUN pip3 install --pre vaex

#Install zindex
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

#Install Cairo
RUN apt-get install -y libgtk2.0-dev \
                       libcairo2-dev \
                       xvfb \
                       xauth \
                       xfonts-base \
                       libxt-dev \
                       vim

#Install google fonts
RUN wget https://github.com/google/fonts/archive/master.zip -O /usr/share/fonts/google_fonts.zip
RUN cd /usr/share/fonts; unzip google_fonts.zip
RUN fc-cache -f -v /usr/share/fonts/

#Install R and dependencies
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

#Import fonts                
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


RUN install2.r --repos ${CRAN_MIRROR}\
                plotly\
                DT\
                data.table\
                tableHTML\
                && Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('rhdf5');  biocLite('BiocParallel')" \
                && Rscript -e "devtools::install_github(c('thomasp85/shinyFiles'))" \
                && rm -rf /tmp/downloaded_packages/ /tmp/*.rds
  
#Install scsva
COPY scSVA_0.1.0.tar.gz /home/
RUN Rscript -e "install.packages('/home/scSVA_0.1.0.tar.gz', repos = NULL, type='source')" \
               && rm /home/scSVA_0.1.0.tar.gz

# Run rocker/rstudio init
CMD ["/init"]
