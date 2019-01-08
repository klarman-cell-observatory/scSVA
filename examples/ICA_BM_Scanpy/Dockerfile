FROM mtabaka/scsvatools
ENV CRAN_MIRROR http://cran.rstudio.com

RUN apt-get update && apt-get install -y  libxml2-dev \
            lib32z1-dev 

RUN pip3 install scanpy python-igraph louvain 

