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
            libglu1-mesa-dev \
            freeglut3-dev \
            mesa-common-dev \ 
            unzip \
            cmake \
            zlib1g-dev \
            htop 
      
            
# Install Python
RUN apt-get -qq update && apt-get -qq -y install curl bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
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
ENV CLOUD_SDK_VERSION=226.0.0
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

#Remove local/openblas* 
#https://github.com/JuliaLang/julia/issues/11913
RUN rm /usr/local/lib/libopenblas*

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
                colourpicker\
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
                gridExtra\
                magick\
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
                
RUN apt-get update --fix-missing \
                       && apt-get install -y \
                       libxext-dev \
                       libxrender-dev \
                       libxtst-dev \
                       libxss1 \
                       libgconf-2-4 \
                       libnss3 \
                       libasound2 \
                       x11vnc \
                       xvfb \
                       desktop-file-utils 
                       
RUN cp /usr/local/bin/orca /usr/bin/orca
RUN echo '#!/bin/bash\nxvfb-run --auto-servernum --server-args "-screen 0 1024x1024x24" /usr/bin/orca "$@" --enable-webgl ' > /usr/local/bin/orca && \
    chmod +x /usr/local/bin/orca

#Install scsva
RUN echo '20181225845' >/dev/null && Rscript -e "devtools::install_github('klarman-cell-observatory/scSVA',dependencies=TRUE,repos=BiocInstaller::biocinstallRepos())"
 
# Run rocker/rstudio init
CMD ["/init"]

