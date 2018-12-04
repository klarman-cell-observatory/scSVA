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
      unzip \
      cmake \
      zlib1g-dev \
      htop \
      libglu1-mesa-dev \
      freeglut3-dev \
      mesa-common-dev 
      
            
# Install Python
RUN apt-get -qq update && apt-get -qq -y install curl bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
#    && apt-get -qq -y remove curl bzip2 \
#    && apt-get -qq -y autoremove \
#    && apt-get autoclean \
#    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && conda clean --all --yes

ENV PATH /opt/conda/bin:$PATH
#Install numpy + vaex
RUN conda install -v -c anaconda pyopengl
RUN conda install -c anaconda numpy 
RUN conda install -c conda-forge vaex
RUN conda install -c plotly plotly-orca

#Install zindex
RUN git clone https://github.com/mattgodbolt/zindex.git \
      && cd zindex \
      && make \
      && ln -s /zindex/build/Release/zq /usr/bin/zq \
      && ln -s /zindex/build/Release/zindex /usr/bin/zindex 

#Install gsutil
RUN apt-get install -y apt-transport-https \
                       curl \
                       gnupg
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

#Install Image Magick
RUN apt-get update --fix-missing \
	          && apt-get install -y \ 
	          libmagick++-dev

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
                shinyBS\
                shinyjqui\
                shinyAce\
                scales\
                sp\
                plyr\
                ggrepel\
                googleComputeEngineR\
                magick


RUN install2.r --repos ${CRAN_MIRROR}\
                plotly\
                DT\
                data.table\
                tableHTML\
                future\
                future.apply\
                processx\
                plot3D\
                && Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('rhdf5')" \
                && Rscript -e "devtools::install_github(c('thomasp85/shinyFiles'))" \
                && rm -rf /tmp/downloaded_packages/ /tmp/*.rds
  
#Install scsva
COPY scSVA_0.2.0.tar.gz /home/
RUN Rscript -e "install.packages('/home/scSVA_0.2.0.tar.gz', repos = NULL, type='source')" \
               && rm -rf /home/scSVA_0.2.0.tar.gz
               
RUN apt-get install -y libxext-dev \
                       libxrender-dev \
                       libxtst-dev \
                       libxss1 \
                       libgconf-2-4 \
                       libnss3 \
                       libasound2 \
                       x11vnc \
                       xvfb \
                       desktop-file-utils

#RUN wget -q -O - https://dl-ssl.google.com/linux/linux_signing_key.pub | apt-key add - && \
#     sh -c 'echo "deb [arch=amd64] http://dl.google.com/linux/chrome/deb/ stable main" >> /etc/apt/sources.list.d/google-chrome.list' && \
#     apt-get update -y && \
#     apt-get install -y google-chrome-stable

RUN cp /usr/local/bin/orca /usr/bin/orca
RUN echo '#!/bin/bash\nxvfb-run --auto-servernum --server-args "-screen 0 1024x1024x24" /usr/bin/orca "$@" --enable-webgl ' > /usr/local/bin/orca && \
    chmod +x /usr/local/bin/orca

RUN sudo apt-get install ssh -y
RUN mkdir --mode=700 /root/.ssh

# Run rocker/rstudio init
CMD ["/init"]
#CMD ["R", "-e", "shiny::runApp('/usr/local/lib/R/site-library/scSVA/scSVA/',port=8787"]

