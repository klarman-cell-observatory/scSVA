### scSVA v0.2

---

The **scSVA** (single-cell Scalable Visualization and Analytics) is an R package for 
interactive visualization and exploratory analysis of massive 
single-cell omics datasets. The **scSVA** is built with [Shiny](https://shiny.rstudio.com) and is optimized for efficient 
visualization of cells on 2D/3D embedding and extracting cell features to visualize them 
from compressed big expression matrices stored on disk in [HDF5](https://www.hdfgroup.org/solutions/hdf5/), [Loom](http://loompy.org) and text file formats.
This reduces the memory resources needed to explore scRNA-Seq datasets by a factor of ~20,000.
The **scSVA** is able to visualize and explore interactively hundreds of millions of cells on a laptop 
or a billion cells on a moderate desktop computer. 
As a back-end it uses [VaeX](https://github.com/maartenbreddels/vaex), 
a fast python library for vector data processing on a grid. 

The **scSVA** simplifies a production of high-quality figures for scientific publications using [ggplot2](https://ggplot2.tidyverse.org) package and provides a comprehensive set of interactive tools for 2D/3D figure customization and annotation.
The **scSVA** allows for basic statistical analysis like computing cell counts and distributions of gene expression 
values across selected or provided groups of cells. 
In addition, users can run fast methods for diffusion maps and 3D force-directed layout embedding (FLE) interactively using 
[scSVAtools](https://github.com/broadinstitute/scSVAtools).
The full documentation is provided with the **scSVA** package in the "Help" tab.

A visualization of 100 Million cells (upsampled FLE created from 
the Human Bone Marrow single-cell [dataset](https://preview.data.humancellatlas.org)) 
on MacBook Pro (3.1 GHz i7, 16 GB) [![100 Million cells](inst/scSVA/www/movie.png?raw=t)](https://youtu.be/HoeIh7JnMoE):

---

### Installation

The **scSVA** package can be installed from GitHub as follows:

```
install.packages(devtools)
devtools::install_github("broadinstitute/scSVA",dependencies=TRUE,repos=BiocInstaller::biocinstallRepos())
```
You need R>=3.4.0 and Rstudio to be installed on your system to install and run **scSVA** package. 

### Prerequisites

The **scSVA** requires Python 3.6 to be installed on your machine. 
The most convenient way is to install Python through **Miniconda3** or **Anaconda3**.
Follow the instructions on https://conda.io/docs/user-guide/install/index.html. 

Then, install **numpy**, **pyopengl**, **vaex**, and **orca** packages:
```
conda install -c anaconda numpy
conda install -c anaconda pyopengl
conda install -c conda-forge vaex
conda install -c plotly plotly-orca
```

The **scSVA** uses [zindex](https://github.com/mattgodbolt/zindex) to create an index on gene names
from compressed expression matrices in a text format. 
If you plan using compressed text files to store the expression matrix, install **zindex** 
following the instruction on https://github.com/mattgodbolt/zindex.
Make sure that **zindex** is properly installed on your computer and
is visible as an executable file by your operating system. 
[ImageMagick](https://www.imagemagick.org/script/index.php) 
is required to combine multiple graphs, see 
R package [magick](https://github.com/ropensci/magick) for installation instructions.
To run [scSVAtools](https://github.com/broadinstitute/scSVAtools), install  [Google Cloud SDK](https://cloud.google.com/sdk/downloads) 
and  [googleComputeEngineR](https://cloudyr.github.io/googleComputeEngineR/index.html), see "Run scSVA on a cloud" section.


### Start quide

**scSVA** uses [reticulate](https://github.com/rstudio/reticulate) to run Python packages in R.
The first step is to configure the path to Python 3 installation with **numpy** and **vaex** libraries.
To get the default python version used by reticulate, use the following command
```
Sys.which("python")
```
To explore which Python versions are installed, use:
```
reticulate::py_discover_config()
```
To modify the default path e.g. if a Python executable file is in `/opt/anaconda3/bin/python`, use
```
Sys.setenv(PATH = paste("/opt/anaconda3/bin/", Sys.getenv("PATH"),sep=":"))
```
or 
```
reticulate::use_python("/opt/anaconda3/bin/")
```
To run **scSVA**, use the following command in R console
```
scSVA::scSVA()
```
This opens a **scSVA** Shiny app in a default browser. 
If the popup blocker is active in your browser click on "Try Again" 
to open **scSVA** Shiny App in a new tab.
The full documentation on how to use **scSVA** is provided with the package in the "Help" 
tab. 

### Install scSVA as a docker container

scSVA docker image has full installation of Rstudio server with R version 3.5.1 with openblass libraries 
and R dependencies preinstalled for **scSVA**.  
Python version 3.6.7 is preinstalled with Anaconda with **numpy** (version 1.15.4) and **vaex**
(version 1.0.0b7). It contains also **zindex** for creating an index on compressed files and 
gsutil tools for working with [Google Cloud](https://cloud.google.com).  

The first step is to download and install **docker** from https://docs.docker.com/install/ .
Then, pull the docker image using the 
following command in the terminal:

```
docker pull mtabaka/scsva
```

To pull [scSVAtools](https://github.com/broadinstitute/scSVAtools) use

```
docker pull mtabaka/scsvatools
```

Verify the docker image is available:

```
docker images
```

To run the **scSVA** docker container, use the command:

```
docker run -d -p 8787:8787 --rm -m 8g --cpus 4 -v PATH:/home/  -e USER=user_login -e PASSWORD=user_password mtabaka/scsva
```
The --cpus and -m flags specify the number of cores and memory available in the container, respectively.
PATH is a path to the directory with an expression matrix and metadata. 
The files will be available in the /home/ directory in the container.
To run R server, open web browser and visit http://localhost:8787 .
To run **scSVA**, use the command in R console
```
scSVA::scSVA()
```
If the popup blocker is active in your browser, click on "Try Again" to open **scSVA** Shiny App in a new tab.

Alternatively, scsva docker image can be installed using a Dockerfile provided with the **scSVA** package:
```
docker build -t scsva path_to_directory_with_Dockerfile
```

### Run scSVA on a cloud

We focus here on running scSVA using the Google Cloud Compute Engine (GCE):

1. Create and configure a Google Cloud Project, see [instructions](https://cloud.google.com/resource-manager/docs/creating-managing-projects)
2. Download GCE private key in JSON format, see [instructions](https://cloud.google.com/storage/docs/authentication#service_accounts})
3. Download and install [Google Cloud SDK](https://cloud.google.com/sdk/downloads)
4. Initialize and authorize Google Cloud SDK, see [instructions](https://cloud.google.com/sdk/docs/quickstarts)
5. Build a scSVA image from a directory containing Dockerfile provided with the scSVA package
    ```
    gcloud container builds submit --timeout=2h --tag gcr.io/project_id/scsva-image-name .
    ```
6. Create a storage and transfer files with an expression matrix and metadata, see [instructions](https://cloud.google.com/storage/docs/gsutil/commands/mb), e.g.
    ```
    gsutil mb gs://bucket-name/
    ```
    ```
    gsutil cp expression_matrix.h5 gs://bucket-name
    ```
7. Install [googleComputeEngineR](https://cloudyr.github.io/googleComputeEngineR/index.html) which provides a convenient R interface to the GCE
8. Initialize and run ***googleComputeEngineR*** library, for troubleshooting see this [link](https://cloudyr.github.io/googleComputeEngineR/)
    ```
    Sys.setenv("GCE_AUTH_FILE"          = "path_to_gce_privte_key_json_file",
               "GCE_DEFAULT_PROJECT_ID" = "project_id",
               "GCE_DEFAULT_ZONE"       = "us-east1-b")
    ```
    
    ```
    library(googleComputeEngineR)
    ```
  
    This will use GCE zone "us-east1-b", [see](https://cloud.google.com/compute/docs/regions-zones/) 
    for details or use `gce_list_zones("project_id")` for listing available zones. 
    To get information about your project use `gce_get_project()` 
9. Get the container tag
    ```
    tag <- gce_tag_container("scsva-image-name")
    ```
10. Run a virtual machine instance on GCE
    ```
    vm <- gce_vm(template        = "rstudio",
                 name            = "scsva",
                 username        = "user_login",
                 password        = "user_password",
                 dynamic_image   = tag,
                 predefined_type = "machine_type",
                 disk_size_gb    = "disk_size_in_gb")
    ```
    To check for available machine types use `gce_list_machinetype()`.
11. Get `externalIP` from R console returned after step 10. or use
    ```
    gce_list_instances()
    ```
    to list vm instances with `externalIP` address. 
    Alternatively, run `gce_get_external_ip("scsva")`.
12. Wait a few minutes for the `scsva-image-name` docker container to download and launch.
13. Open web browser and visit `http://localhost:externalIP`. 
    Provide `user_login` and `user_password` from the step 10. This will open Rstudio server.
14. To copy expression matrix with metadata from a bucket to your vm instance, go to Rstudio terminal and execute the command
    ```
    gsutil cp gs://bucket-name/* ./  
    ```
   This will copy the content of the bucket to the current local directory on vm instance. 
   
15. To start **scSVA** run 
    ```
    scSVA::scSVA()
    ```
  in the Rstudio console. If the popup blocker is active in your browser click on "Try Again" 
  to open **scSVA** Shiny App in a new tab.

16. Stop the running vm instance using
    ```
    gce_vm_stop("scsva")  
    ```
    User will not be charged for stopped instances (except for storage). 
    You can restart the vm instance by rerunning step 10.  
    User can also delete vm instance by running 
    ```
    gce_vm_delete("scsva")
    ```
  
### Authors

**scSVA** was developed by Marcin Tabaka and  **scSVAtools** by Marcin Tabaka and Joshua Gould from the [Regev Lab](https://www.broadinstitute.org/regev-lab) at the [Broad Institute of MIT and Harvard](https://www.broadinstitute.org).
Use this [link](mailto:mtabaka@broadinstitute.org) to contact me. 

### Bugs reports and feature requests
Use the issue tracker for bugs reporting and new features requesting.
 











