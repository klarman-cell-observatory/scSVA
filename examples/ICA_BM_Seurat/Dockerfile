FROM mtabaka/scsvatools
ENV CRAN_MIRROR http://cran.rstudio.com

RUN apt-get update --fix-missing \
	 && apt-get install -y \
            libssl-dev \
            libcurl4-openssl-dev \
            libhdf5-dev \
            libpng-dev
          
RUN install2.r --repos ${CRAN_MIRROR}\
    Seurat
RUN Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('rhdf5')"

RUN pip3 install --upgrade umap-learn    


ADD http://www.fftw.org/fftw-3.3.8.tar.gz /tmp/
RUN tar -C /tmp -xzf /tmp/fftw-*.tar.gz
RUN cd /tmp/fftw*/ && ./configure && make && make install

#RUN git clone https://github.com/KlugerLab/FIt-SNE.git
#RUN cd /FIt-SNE && g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm  

RUN wget https://github.com/KlugerLab/FIt-SNE/archive/ec25f1b36598a2d21869d10a258ac366a12f0b05.zip \
    && unzip ec25f1b36598a2d21869d10a258ac366a12f0b05.zip \
    && mv FIt-SNE-ec25f1b36598a2d21869d10a258ac366a12f0b05 FIt-SNE \
    && rm  ec25f1b36598a2d21869d10a258ac366a12f0b05.zip \
    && cd /FIt-SNE && g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm
