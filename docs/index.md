# **scSVA**

### **File Upload**

scSVA accepts tabular data stored in Comma-Separated Value (CSV) or in Hierarchical Data Format ( [HDF5](https://www.hdfgroup.org) ) files. The input data consists of (i) X- and Y-coordinates of your favorite embedding (e.g. tSNE or FLE), (ii) gene expression matrix or vectors that contain a gene quantification for each cell, and (iii) metadata information on each cell's group assignment The file loading depends on a file type and is described separately for each file type. The first thing to do is to open a File tab panel. This panel allows for loading all necessary files to scSVA.

#### **Selection of XY Coordinates**

First, select the "Input File Format" and then click on "File Select". This will open a new dialog box with the access to the local file system. Select the file that contains XY coordinates for your plot and click Select button in the dialog box. The output depends on the selected file type:

*   HDF5 - A new table opens below the panel. An HDF5 file includes two types of objects: datasets and groups that store datasets and other groups. The table lists both groups and datasets with detailed information like dataset/group name, storage class, dimension and size of datasets. Select the datasets with X and Y coordinates. scSVA supports HDF5-based [Loom](https://github.com/linnarsson-lab/loompy) files. The filename extension must be from the following set \[.h5,.hdf5,h5ad,.loom\] or from the same set but with all letters capitalized.
*   CSV - When you load a CSV file you make sure that the file has a header with column names. The column names (feature names) will appear in the new table below the panel. scSVA uses "fread" from [data.table](https://cran.r-project.org/web/packages/data.table/) package which reads "regular delimited files" i.e. each row in a file has the same number of columns. The field delimiter is recognized automatically and users should use one from the following set \[\\t ,;:|\]. The above field delimiters should not be used in column and row names. scSVA accepts both raw and compressed files (compressed with gzip) and the file type is recognized automatically. The file name extension must be from the following list \['.csv','.tsv','.txt'\], or with all letters capitalized or compressed equivalents (.gz).

#### **Selection of Expression Matrix**

*   HDF5 - scSVA accepts two dimensional datasets when a "Dataset" is selected or vectors of gene expression values when a "Group" is selected. If the "Dataset" is selected, specify if the rows or columns of the dataset are genes. If the "Group" is selected, scSVA will allow for loading expression vectors (one-dimensional datasets) within a group, with gene symbols as dataset's names. This option is appropriate for very big datasets.
*   CSV - We support only gene expression matrices, where the rows are genes and the columns are cells, and compressed with [Gzip](https://www.gnu.org/software/gzip/) . This layout allows for indexing and very fast retrieval of expression vectors. The first column must contain gene names with "Genes" as the column name. The first row contains Cell IDs. Below is an example of an expression matrix:
    
    | Genes | Cell.1 | Cell.2 | Cell.3 | ... |
    | --- | --- | --- | --- | --- |
    | Gene.1 | 0 | 1 | 0 | ... |
    | Gene.2 | 2 | 1 | 3 | ... |
    | Gene.3 | 1 | 0 | 0 | ... |
    | ... | ... | ... | ... | ... |
  
	To compress the expression matrix file use the following command in the terminal: 
	
	```gzip expression_matrix.csv```
	
	The next step is to upload a list with gene names that matches the gene names in the expression matrix. scSVA makes queries 
	on the compressed expression matrix based on that list. The list can be saved at the preprocessing step or extracted from the 
	uncompressed expression matrix file by the command:  
	```cat expression_matrix.csv |cut -d',' -f1 | awk '{if(NR>1) print} >Gene_Names.txt'```
    
    or from the compressed file:  
    ```gunzip -cd expression_matrix.csv.gz |cut -d',' -f1 | awk '{if(NR>1) print} >Gene_Names.txt'```  
        
    The index can be created on the compressed file by users using [Zindex](https://github.com/mattgodbolt/zindex) in the terminal :  
    ```zindex expression_matrix.csv.gz --field 1 --tab-delimiter --skip-first 1```  
    for tab-delimited files or  
    ```zindex expression_matrix.csv.gz --field 1 --delimiter , --skip-first 1```  
    for e.g. comma-delimited files. scSVA can create the index for you if a tick box "Create index on Compressed File" is checked. Make sure that [Zindex](https://github.com/mattgodbolt/zindex) is properly installed on your computer and is visible as an executable file by your operating system. Indexing allows for very efficient retrieval of the rows from compressed files without need for loading of the entire matrix into the memory. Thus, the memory usage is reduced by a factor of 10,000-40,000 for a typical expression matrix which makes possible a fast visualization of gene expression values on a 2D embedding or computing and plotting scores of a custom defined gene signature (set of genes).

  
#### **Selection of Groups**

Selection of a vector of group IDs is identical to the selection of XY coordinates. Group IDs must be a vector of integer (categorical data).

### **Plot Appearance**

The tab "Visualize" helps customizing the plot appearance. scSVA support visualization of 2D and 3D cell embeddings.

#### **2D**

The default plot is a density plot, where for a given grid resolution, it shows the number of cells within a square. The grid resolution can be modified by a scroll bar below the plot. The default resolution is 512 which means that a plot consists of 512x512=262144 square cells. Users can increase the resolution of the grid to produce higher quality figures but this results in increasing the plot refreshing time as more points needs to be plotted. The tradeoff between the figure quality and the interactive plot responsiveness needs to be exploited by users as it depends on the number of cells and machine-specific factors. We found that the default value of 512 is sufficient to produce high quality figures and the response is fast enough to interactively plot hundreds of millions of cells.  
The two sidebars in the "2D" tab simplify customization of the plot. Users can modify 1) size of the plotted points, 2) size and type of the font (x,y-axis labels, title, color bar), 3) transparency of the points by specifying the transparency value, 4) color palette of the colorbar, 5) colors of the panel and plot background. We included also a possibility of checking the colorschemes in terms of suitability for color blind readers ( [dichromat package](https://cran.r-project.org/web/packages/dichromat/index.html) ). Clicking on "Dichromatism" box will open options for various types of colorblindness (green, red, blue) and users can simulate the effects of the type of dichromatism which helps selecting a color palette suitable for people with different types of color perception deficiency. We provide also color palettes suitable for dichromats from [dichromat package](https://cran.r-project.org/web/packages/dichromat/index.html) . These color palette names begin with "dichromat::". The second sidebar modifies 1) the main title of the plot, 2) x- and y-axis label, 3) color bar title. Finally, the main plot can be saved in many file formats (vector or raster) with a custom resolution and size.  

##### **Plot Navigation**

Mouse left-click starts a selection of a rectangular region. Zooming in: double-click inside the selected area will zoom in the selected region with the specified resolution ("Grid Resolution" scrollbar). "Reset View" button below the main plot when clicked will restore the original plot xy range with the specified resolution.

#### **3D**

To initialize 3D plot, select X, Y, and Z coordinates and load datasets (see File section). Functions for 3D plot customization are similar to those present in 2D plots except the Camera Viewpoint which determines the camera view point about origin ( [see manual](https://plot.ly/python/3d-camera-controls/) ). A panel "Display Feature" changes the feature to be plotted. Users can visualize a density of cells or gene/signature expression values. 2D projection button opens a new window with 3D data projection onto 2D.

### **Analysis**

Analysis Menu consists of tools for exploratory analysis of uploaded datasets.

#### **Cell Abundance & Gene Expression**

Here users can explore number of cells in selected areas, plot gene expression profiles on 2D embedding, and compute gene expression distributions and statistics like mean, median, standard deviation in the selected areas.  
Users can use two selection tools: (i) rectangular selection (default), and (ii) polygonal selection. Mouse left-click on the plot and holding down the left button lets you create a rectangular selection. To start polygonal selection check "Polygonal Selection" box. From now, each mouse left-click on the plot will create a new polygon's vertex and connect it by edges with the first and the penultimate vertex.  
To compute number of cells, select a region of interest and click "Group A" button. This will mark cells in a selected area and return number of cells in that area. Users can select two groups of cells. "Reset Group" button will unselect cells. To plot gene expression values on 2D or 3D embedding, select a gene in right-hand panel. By default, the mean expression values are plotted on the grid but users can change the statistic by selecting "Select Statistic". Three options are supported: mean, median (the slowest), and maximal expression value. The mean is suitable for highly expressed genes and the maximum value option works well for lowly expressed genes (e.g. transcription factors or surface receptors in scRNA-Seq). For a chosen gene, users can select two groups of cells as described above and compute the gene expression statistics in these two groups by clicking on "Get Group Statistic" button.

#### **Gene Expression by Group**

Checking box "Compute Gene Expression by Group" starts computation of statistics in groups like number of cells, gene expression mean, median and standard deviation. Once completed, the statistics are listed in the table. The table is fully interactive: 1) groups can be sorted by a mean or number of cells for example, 2) searched by statistic values or group IDs, 3) filtered by statistic values by clicking on boxes below the table. Single-click on a row of the table (group) will show a gene expression distribution for the selected group. Users can mark many groups to compare the expression distributions among them. Each selection opens a new distribution plot on right-hand side. The plots are interactive and colored by the same colors as groups on the embedding. The scroll bar, on the bottom of the plots, changes the y-axis range for all plots. Moving a cursor to the plot bar will show fraction of cells falling in a given expression range. To compare all the distributions on one plot go to "ECDF" tab to open a plot with empirical cumulative distribution function (users can easily read dropout levels from the ecdf plot or check for bi- or multimodalities that may suggest underclustering) or "Violin" to visualize group expression values by violin plot. "Census" tab contains the fractions of cells in each group.

#### **Gene Signatures**

In order to compute gene signature scores (e.g. from Molecular Signature Database [MSigDB](http://software.broadinstitute.org/gsea/msigdb) ) for all cells and visualize them on 2D or 3D embedding, users need to provide a set of gene names to "Enter Gene List" field and specify "Gene Signature Name". The gene names in the list must be in the same format as in gene expression matrix (e.g. gene symbols). The gene expression values across all cells can be retrieved in parallel from an expression matrix by specifying the number of threads on local machine. Click on "Compute Gene Signature Score" box to start retrieving and computing the scores. scSVA computes gene expression z-score across all cells and average over the z-scores across the entire gene set.

### **Cloud**

scSVA supports interactive analytics and data storage mamagement on the cloud (Google Cloud Platform). Google Cloud Platform charges you for storage and running VMs, [see pricing](https://cloud.google.com/pricing/) . Stop or delete any GCP services once you are done with your analysis. Good practice is to often check the billing of your project. scSVA is an R package distributed as-is, without warranties of any kind. Make sure that you use the GCP services reasonably and follow best practices. To get started make sure that your system configuration is ready to use GCP and in particular to launch virtual machine instances in Google Cloud Compute Engine (GCE). Follow the instructions below to install and configure all required tools:  
(i) Create and configure a Google Cloud Project, see [instructions](https://cloud.google.com/resource-manager/docs/creating-managing-projects) .  
(ii) Download GCE private key in JSON format, see [instructions](https://cloud.google.com/storage/docs/authentication#service_accounts) .  
(iii) Download and install [Google Cloud SDK](https://cloud.google.com/sdk/downloads) .  
(iv) Initialize and authorize Google Cloud SDK, see [instructions](https://cloud.google.com/sdk/docs/quickstarts) .  
See [tutorial](https://cloud.google.com/getting-started/#quick-starts) for more information. scSVA uses [googleComputeEngineR](https://cloudyr.github.io/googleComputeEngineR/) , which provides R interface to the GCE API.

#### **Storage**

To list google buckets from the default project, choose "gs" and click "Explore FS". Select a bucket from the list and click right arrow below the list to show its content. The upper arrow below the list navigates up one directory. To remove the object from the bucket, select it and click "Remove Selected Object". Users can also create new buckets by specyfying the name of the bucket and clicking "Create New Bucket" button. Selection of "local" filesystem will list directories and files on your local computer. Users can easily copy objects between the local computer and the bucket by clicking left/right arrow between filesystem panels.

#### **Compute Engine**

scSVA uses [googleComputeEngineR](https://cloudyr.github.io/googleComputeEngineR/) to lunch and operate VMs. In the left corner, there is a link to Google Cloud Console. Clicking it opens the Console in new browser tab. Users can monitor there all services provided by GCP. To get started, go to "GCE Setup" and provide the path to GCE private key in JSON format (see step (ii)), your project name (see step (i)), and the zone you want to launch VM (see [instructions](https://cloud.google.com/compute/docs/regions-zones/) ) , After successful configuration, "GCE List Disks" in "GCE Info" allows to display all disks in the project, and "GCE List Instances" list all VMs (stopped and running) in the project. Right console displays all the information from GCE To launch a new VM, or restart a stopped VM, go to "RUN VM" window, and click "GCE List Machine Types". This will list all available Machines in your selected zone. Select container with computational tools you want to run on VM, and provide name for your virtual machine. To build a new image on Google Cloud, go to terminal, then to directory with a Dockerfile and run the following command:  
`gcloud container builds submit --timeout=2h --tag gcr.io/Project_Name/ContainerName`  
For further information, see instructions to [googleComputeEngineR](https://cloudyr.github.io/googleComputeEngineR/) package. Usually it takes less than one minute to lunch VM. You can check the progress in Rstudio console. If the VM instance is successfully launched you should see its name in the list of VM instances with status "RUNNING" after clicking the "GCE List Instances" button. "Get Running VMs" button displays all VM instances with status "RUNNING". Specify the VM on which you want to run jobs and click "Assign VM". Notice the change in the name of "Selected VM" above the right console. From now, you can run commands and containerized tools and copy files between local machine and VM instance. You can also stop or delete the VM instance and connect to VM instance through ssh. If you stop the VM instance you are charged only for storage, and the restarted VM instance will have access to all files from previous session. If you delete VM instance, it will permanently remove a disk and VM instance (you will not be able to restart it in the future). In order to list running processes on the VM instance click "Run Top". Users can easily transfer files from/to the VM instance using the panel below. Selected files from the list can be modified after clicking "Edit" button. Users can also create new files on the VM instances. If your container has tools installed that can run scripts (eg. bash, R, or python) selecting a script on VM and clicking "Run Script" will start running the script on the VM instances. If you have many running VM instances you can switch between them by selecting a running instance and assigning it. We added also the option of running commands on the VM instance. Write the command in the test field below the console (as you would do it in the bash terminal) and click "Run SSH Command". The output of the command appears in the console. To remove all text from the console, click "Clear" button.

### **Metadata**

#### **Groups**

Here, users can visualize and annotate cell groups on 2D embedding by selecting color palettes, modify colors of individual groups, show and modify cluster labels.  
The middle panel contains a table with coordinates of group labels and right-hand panel customizes the plot. First check "Show Groups". This will compute group positions (cell groups for example) on a grid. "Get Labels" check box, when selected, computes the position of labels by taking average x- and y-coordinates for each group. The table in the middle panel lists all groups (group IDs, x-, and y-coordinates, and group labels). Users can modify the coordinates of labels plotted on 2D and group names by double-click on the value or label in the table. The plot is updated after any modification. "Get Labels" check box opens a panel in which users can upload group names from a file. The size (number of rows in case of csv file or number of elements in hdf5 dataset) must be equal to the total number of groups. "Repel Labels" check box distribute labels in such a way that the labels don't overlap ( [ggrepel package](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html) ), which helps distributing the labels in such a way that they don't overlap. The labels can be further customized by removing the bounding box or changing font type, size, color, etc. in the "Annotation" tab.  
To show a legend with group annotations check "Show Legend". Users can modify the legend title, position, number of columns, order of groups, and displayed point size in the panel that opens after checking the box "Show Legend" Colors of cell groups on 2D embedding can be modified in "Color Palette" panel. We included support for most of [Crayola crayon colors](https://en.wikipedia.org/wiki/List_of_Crayola_crayon_colors) which generate distinctive and aesthetically-pleasing color sets to mark groups of cells. Clicking on Shuffle Colors will permute colors from a given set. Make sure that none of the groups has a color which is close to the background color. "Shuffle Cells" will permute the cells on the grid (important if some groups overlap for example). The individual group colors can be changed in two ways: (i) checking "Show color List" will open a table with group colors. Double-click on the column element will allow for color modification for the selected group. Users can provide a color in hex format or by giving a color name e.g. "navyblue"; (ii) Specify the "Group ID" whose color you wish to modify (the group can be highlighted by checking the "Highlight Group" box). Clicking "Change Group Color" field will open a color picker. Click on a color and the color of the cell group will change immediately.

#### **Multiplots**

Here, users can arrange multiple graphs on a grid. To start click on "Add Plot" button. This will add the main plot to a list of plots. Check the box "Show Multiplots" to make it visible. To add more graphs just click "Add Plot". To change the size of added plot or its resolution go to "Visualaize", "2D" tab and modify the plot properties on the right sidebar. scSVA uses [R bindings](https://cran.r-project.org/web/packages/magick/index.html) to [ImageMagick](https://www.imagemagick.org/script/index.php) to arrange graphs on a grid. Multiplots can be zoomed in/out by changing the "Scale" scroll bar, "Number of Columns" specify number of graphs per row.

### **Annotation**

To start annotating the main plot check the "Annotate" box. From now, any single mouse left-click on the main plot will create a new label. The coordinates and text labels appear in the table in the middle panel. The table is interactive and double click on its element allow for modification of this element. To remove an annotation label, specify the "Annotation ID" you wish to remove and click on "Remove Label". The labels can be further customized by specifying the background color of the label, its transparency, size, and font type of the text. The bounding box can be modified by changing the line width of its border or the radius of corners.

### **Extensions**

#### **Color Palette**

Users can customize their plots by creating their own color palettes. Click on "Select Colors For a New Palette". This will open a color picker widget that can be used to select colors for a custom palette (the colors can be added to a palette also by specifying the hex color code). Specifying the name of the color palette and clicking on "Add New Palette" will make it available in both "Color Palette" selection fields in "Visualize" and "Metadata" -> "Groups" tabs. The new custom color palette will appear at the end of the list of color palettes.

#### **Fonts**

To further customize plots, users can upload their own or open source fonts. scSVA accepts TrueType fonts (filename extension .ttf). This functionality is available on Mac or Linux operating systems. Windows users may consider installing scSVA as a docker container. The default installation path of fonts on linux is /usr/share/fonts, users will be asked to provide a password in the terminal to write the new fonts to this directory. In some cases, to see new fonts on the interactive or saved figure, it may be necessary to restart the current R session. For example, write in RStudio console: .rs.restartR() or use Command/Ctrl + Shift + F10 to restart the R session. We encourage users to import new fonts at the first step and restart R session before loading files with x-, and y- coordinates.  

#### Guided Video Tutorials:

Created by Marcin Tabaka  
Contact me at [mtabaka@broadinstitute.org](https://mail.google.com/mail/?view=cm&fs=1&tf=1&to=mtabaka@broadinstitute.org)
