ui <- fluidPage(theme = shinytheme("simplex"),
  navbarPage("scSVA: single cell Scalable Visualization and Analytics",
  tabPanel("File",
                 fluidRow(
                   column(width=4,
                     wellPanel(
                       h4("Select XY coordinates"),
                       selectInput(
                         inputId  =  "inputFileType_XY",
                         label    =  "Input File Format",
                         choices  =  c("HDF5", "CSV"),
                         selected =  "HDF5"),
                       shinyFilesButton(
                         id       = 'inputFile_XY',
                         label    =  'File Select',
                         title    =  'Select a file',
                         multiple =  FALSE),
                       style='border:32px;',
                       selectInput(
                         inputId  = "X_coord",
                         label    = "X-coordinate",
                         choices  = c("NULL"),
                         selected = "NULL"),
                       selectInput(
                         inputId  = "Y_coord",
                         label    = "Y-coordinate",
                         choices  = c("NULL"),
                         selected = "NULL"),
                       hr(),
                       h4("Select Expression Matrix"),
                       selectInput(
                         inputId  = "inputFileType_Genes",
                         label    = "Input File Format",
                         choices  = c("HDF5", "CSV"),
                         selected = "HDF5"),
                       shinyFilesButton(
                         id       = 'inputFile_Genes',
                         label    = 'File Select',
                         title    = 'Select a File',
                         multiple = FALSE),
                       conditionalPanel("input.inputFileType_Genes=='HDF5'",
                         radioButtons(
                           "h5_type",
                           "Select HDF5 type",
                           choices = c(Dataset   = "D",
                                       Group     = "G"),
                                        selected  = "D"),
                       selectInput(
                         inputId  = "inputGenes",
                         label    = "Expression Matrix",
                         choices  = c("NULL"),
                         selected = "NULL")),
                       conditionalPanel("input.inputFileType_Genes=='HDF5' & input.h5_type=='D'",
                         selectInput(
                           inputId  = "inputGeneVector",
                           label    = "Gene Names",
                           choices  = c("NULL"),
                           selected = "NULL"),
                         radioButtons("em_or", "Expression Matrix:",
                           choices = c(`Rows Are Genes`    = "R",
                                       `Columns Are Genes` = "C"),
                                       selected            = "R")),
                       conditionalPanel("input.inputFileType_Genes=='CSV'",
                         shinyFilesButton(
                         id       = 'inputFile_GenesVector',
                         label    = 'File select (Gene Names)',
                         title    = 'Select a file with Gene Names',
                         multiple = FALSE),
                       checkboxInput(
                         'computeIndex',
                         'Create Index on Compressed File',
                         FALSE)),
                       hr(),
                       h4("Select Clusters"),
                       selectInput(
                         inputId  = "inputFileType_Clusters",
                         label    = "Input File Format",
                         choices  = c("HDF5", "CSV"),
                         selected = "HDF5"),
                       shinyFilesButton(
                         id        =  'inputFile_Clusters',
                         label     =  'File Select',
                         title     =  'Select a File',
                         multiple  =  FALSE),
                       selectInput(
                         inputId  = "inputClusters",
                         label    = "Clusters",
                         choices  = c("NULL"),
                         selected = "NULL")),
                       actionButton("load_datasets", "Load Datasets")
                   ),
                   column(width=4,
                          DT::dataTableOutput("file_cont_XY",       height = 275),
                          DT::dataTableOutput("file_cont_Genes",    height = 275),
                          DT::dataTableOutput("file_cont_Clusters", height = 275)
                   )
                 )),
#View
  tabPanel("Visualize",
                   fluidRow(
                     column(width = 6,
                       plotOutput("plot1",
                         height   = 700,
                         click    = "plot1_click",
                         dblclick = "plot1_dblclick",
                         brush    = brushOpts(id = "plot1_brush")
                       ),
                       HTML("<br>"),
                       splitLayout(
                         sliderInput("grid",
                         "Grid Resolution:",
                         min   = 128,
                         max   = 5120,
                         value = 512,
                         step  = 1
                         ),
                         sliderInput("ColorbarRange",
                           "Colorbar Range:",
                           min   = 0,
                           max   = 1,
                           value = c(0,1),
                           step  = 0.1
                         )
                       ),
                       splitLayout(
                         actionButton("reset_view", "Reset View"),
                         actionButton("change_colorbar", "Change Colorbar")
                       )
                     ),
                     column(width=3,
                       offset=1,
                       wellPanel(style="max-height: 800px",
                         sliderInput("ps",
                           "Pointsize",
                           min   = 0.01,
                           max   = 3,
                           value = 0.1,
                           step  = 0.01),
                         sliderInput("fs",
                           "Fontsize",
                           min   = 6,
                           max   = 24,
                           value = 12,
                           step  = 1),
                         sliderInput("transparency",
                           "Transparency",
                           min   = 0.0,
                           max   = 1,
                           value = 1.0,
                           step  = 0.02),
                         selectInput(
                           inputId  = "pal",
                           label    = "Color Palette",
                           choices  = c("inferno", "viridis", "magma", "rainbow","matlab.like","matlab.like2","matlab.like3","gray.red"),
                           selected = "inferno"),
                         selectInput(
                           inputId  = "fonts",
                           label    = "Fonts",
                           choices  = fonts(),
                           selected = "Arial"),
                         colourInput("panel_background",
                           label="Panel Background",
                           "white"),
                         colourInput("plot_background",
                           label="Plot Background",
                           "white"),
                         selectInput(
                           inputId  = "dichromat",
                           label    = "Dichromatism",
                           choices  = c("No", "Deuteranopia (Green-Blind)"="deutan",
                                        "Protanopia (Red-Blind)"="protan", "Tritanopia (Blue-Blind)"="tritan"),
                           selected = "No"),
                         splitLayout(
                           checkboxInput("showlegend",
                             "Show Legend",
                             TRUE),
                           checkboxInput("showgrid",
                             "Show grid",
                             FALSE)),
                         splitLayout(
                           checkboxInput("showaxes",
                             "Show Axes",
                             TRUE),
                           checkboxInput("dens_scale",
                             "Scale Colorbar, log10",
                             TRUE)),
                         splitLayout(
                           checkboxInput("point_shape",
                             "Square Points",
                             FALSE))
                       )
                     ),
                     column(width=2,
                       wellPanel(style="min-height: 800px",
                         textInput("title",
                           "Title",
                           ""),
                         textInput("x_axis",
                           "X-axis Label",
                           "Dim 1"),
                         textInput("y_axis",
                           "Y-axis Label",
                           "Dim 2"),
                         textInput("clegend",
                           "Colorbar Title",
                           "Density"),
                         hr(),
                         selectInput("filetype",
                           "Save As:",
                           choices = c("pdf", "jpeg", "tiff", "png", "bmp", "svg","eps", "ps", "tex")),
                         numericInput("res",
                           "Resolution:", 300,
                           min = 1,
                           max = Inf),
                         selectInput("un",
                           "Unit:",
                           choices = c("in", "cm","mm")),
                         numericInput("width",
                           "Width:", 8.50,
                           min  = 1,
                           max  = Inf,
                           step = 0.05),
                         numericInput("height",
                           "Height:", 8.50,
                           min  = 1,
                           max  = Inf,
                           step = 0.05),
                         downloadButton("downloadPlot", "Download")
                       )
                     )
                   )
                 ),
#Analysis:Abundance & Expression
  navbarMenu("Analysis",
                   tabPanel("Cell Abundance & Gene Expression",
                     fluidRow(
                       column(width = 6,
                         plotOutput("plot3",
                           height   = 700,
                           click    = "plot1_click",
                           dblclick = "plot1_dblclick",
                           brush    = brushOpts(id = "plot1_brush")
                         ),
                         checkboxInput("select_polygon",
                           "Polygonal Selection",
                           FALSE),
                         actionButton("group_1",
                           "Group A",
                           style="color: #fff; background-color:  #337ab7; border-color: #2e6da4"),
                         actionButton("group_2",
                           "Group B",
                           style="color: #fff; background-color:  #177F75; border-color: #2e6da4"),
                         actionButton("reset_polygon",
                           "Reset Polygon"),
                         actionButton("reset_groups",
                           "Reset Groups"),
                         actionButton("reset_view2",
                           "Reset View")),
                       column(width = 3,
                         h6("Number of Cells in Group A:"),
                         verbatimTextOutput("count_cells_Group_A"),
                         h6("Number of Cells in Group B:"),
                         verbatimTextOutput("count_cells_Group_B"),
                         HTML("<br>"),
                         plotlyOutput("plot_GEDistGroupA",
                           height = 260,width = 380),
                         HTML("<br>"),
                         plotlyOutput("plot_GEDistGroupB",
                           height = 260,width=380),
                         tableOutput("table2")
                       ),
                       column(width = 3,
                         wellPanel(
                           selectInput(inputId  = "choose_gene",
                             label    = "Select a Gene:",
                             choices  = "",
                             selected = "",
                             multiple = F),
                           selectInput(inputId  = "select_method",
                             label    = "Select a Statistic:",
                             choices  = c("Mean","Median","Max"),
                             selected = "Mean",
                             multiple = F),
                           checkboxInput("plot_expression",
                             "Show Gene Expression",
                             FALSE),
                           hr(),
                           actionButton("get_statistic",
                             "Get Group Statistic")
                         )
                       )
                     )
                   ),
#Analysis:Gene Expression in Clusters
                   tabPanel("Gene Expression by Clusters",
                     fluidRow(
                       column(width = 5,
                         checkboxInput("get_exp_clusters",
                                      "Compute Gene Expression by Cluster",FALSE),
                         hr(),
                         DT::dataTableOutput('table_exp.in.clusters')
                         ),
                       column(width = 7,
                         checkboxInput("plot_ecdf",
                           "Compute ECDF",
                           value=FALSE),
                         plotlyOutput("plot_DistGEinCluster",
                           height = 600),
                         conditionalPanel(condition='plot_ecdf',
                           sliderInput("y_axisR",
                             "Y-Axis Range:",
                             min   = 0,
                             max   = 1,
                             value = 1,
                             step  = 0.001
                           )
                         )
                       )
                     )
                   ),
#Analysis:Gene Signatures
                   tabPanel("Gene Signatures",
                     fluidRow(
                       column(width=3,
                         wellPanel(style="max-height: 800px",
                           textInput("GeneSetName",
                             "Gene Signature Name",
                             ""),
                           textAreaInput("gene_sign",
                             "Enter Gene List",
                             "",
                             width = "345px",
                             height= "500px",
                             resize= "vertical"),
                           numericInput("n_workers",
                             "Number of Workers", 1,
                              min  = 1,
                              max  = Inf,
                              step = 1),
                           actionButton("get_gene_sign_score",
                             "Compute Gene Signature Score")
                           )
                         )
                       )
                     )
                   ),
#Metadata:Clusters
  navbarMenu("Metadata",
                   tabPanel("Clusters",
                     fluidRow(
                       column(width = 6,
                         plotOutput("plot4",
                           height   = 700,
                           click    = "plot1_click",
                           dblclick = "plot1_dblclick",
                           brush    = brushOpts(id = "plot1_brush")
                         ),
                         actionButton("reset_view3", "Reset View")
                       ),
                       column(width = 3,
                         DT::dataTableOutput('table_clusters')),
                       column(width = 3, wellPanel(
                         checkboxInput("show_clusters",
                           "Show Clusters",
                           FALSE),
                         hr(),
                         checkboxInput("plot_clusters",
                           "Get Cluster Label Positions",
                           FALSE),
                         checkboxInput("repel_labels",
                           "Repel Labels",
                           FALSE),
                         checkboxInput("show_text",
                           "Remove Bounding Box",
                           FALSE),
                         hr(),
                         selectInput(
                           inputId  = "pal.clusters",
                           label    = "Color Palette",
                           choices  = c("rainbow",
                                        "matlab",
                                        "matlab2",
                                        names(list.colors)
                             ),
                         selected = "rainbow"),
                         actionButton("shuffle_colors",
                           "Shuffle Colors"
                         ),
                         actionButton("shuffle_cells",
                           "Shuffle Cells"
                         ),
                         hr(),
                         numericInput("cluster_id",
                           "Cluster ID:", 1,
                           min  = 1,
                           max  = Inf,
                           step = 1),
                         checkboxInput("highlight_cluster",
                           "Highlight Cluster",
                           FALSE),
                         colourInput("colCluster",
                           label="Select Color for Cluster",
                           "white"),
                         checkboxInput("show_colors",
                           "Show Color List",
                           FALSE)
                         ),
                         conditionalPanel('show_colors == TRUE',
                         DT::dataTableOutput('table_colors'))
                       )
                     )
                   )
  ),
#Analysis:Annotation
  tabPanel("Annotation",
                   fluidRow(
                     column(width = 6,
                       plotOutput("plot2",
                         height   = 700,
                         click    = "plot1_click",
                         dblclick = "plot1_dblclick",
                         brush    = brushOpts(id = "plot1_brush")
                       )
                     ),
                     column(width = 3,
                       DT::dataTableOutput('table1')
                     ),
                     column(width =3,
                       wellPanel(
                         checkboxInput("annotate",
                           "Annotate",
                           FALSE),
                         numericInput("row_id",
                           "Annotation ID:", 1,
                           min  = 1,
                           max  = Inf,
                           step = 1),
                         actionButton("remove_label",
                           "Remove Label"),
                         hr(),
                         colourInput("annotation_background",
                           label  =  "Background Color",
                           "white"),
                         sliderInput("annotation_transparency",
                           "Transparency",
                           min   = 0.0,
                           max   = 1,
                           value = 1.0,
                           step  = 0.02),
                         sliderInput("annotation_size",
                           "Size",
                           min   = 1,
                           max   = 16,
                           value = 5,
                           step  = 1),
                         selectInput(inputId = "annotation_fonts",
                           label    = "Fonts",
                           choices  = fonts(),
                           selected = "Arial"),
                         selectInput(inputId = "annotation_fontface",
                           label    = "Fontface",
                           choices  = c("plain","italic","bold"),
                           selected = "plain"),
                         sliderInput("annotation_labelSize",
                           "Label Border Size",
                           min   = 0,
                           max   = 2,
                           value = 0.5,
                           step  = 0.05),
                         sliderInput("annotation_radiusCorners",
                           "Radius of Rounded Corners",
                           min   = 0,
                           max   = 1,
                           value = 0.15,
                           step  = 0.05)
                       )
                     )
                   )
  ),
#Extensions:Colors
  navbarMenu("Extensions",
                   tabPanel("Color Palette",
                     colourInput("ColPallete",
                     label="Select Colors For a New Pallete",
                     NULL),
                   textInput("NewPaletteName",
                     "New Palette Name",
                     ""),
                   actionButton("RemoveLastColor",
                     "Remove Last Color"),
                   actionButton("RemoveAllColors",
                     "Remove All Colors"),
                   plotOutput("PlotPalette",
                     width = "300px",
                     height = "250px"),
                   actionButton("AddNewPalette",
                     "Add New Palette ")
                   ),
#Extensions:Fonts
                   tabPanel("Fonts",
                     sidebarPanel(
                       fileInput('file_fonts',
                                 'Choose file to upload',
                       accept = c("font/ttf")
                       ),
                       h5("Accepting TrueType fonts only (.ttf)")
                       )
                     )
                   ),
  tabPanel("Documentation",
    style = "font-size: 12pt; line-height: 20pt; color: black;",
    column(width=2),
      column(width=8,

tags$br(),
tags$br(),
tags$br(),
tags$br(),
tags$br(),

h3(tags$b("Motivation")),

"The main aim of the scSVA (link to github) is to simplify exploratory data analysis in the field of single-cell transcriptomics.
The progress in the development of high-througput single-cell methods has allowed for study of cell types and states of
hundreds of thousands of cells. The interest of scientist in understanding biology at single cell level is not a fad but rather an ultimate goal -
e.g. the  ,",a("Human Cell Atlas Consortium", href="https://www.humancellatlas.org", target="_blank"), " is planning to profile at least 10 billion human cells
across all tissues. The massive growth of data volumes makes not only the data analysis harder but also put constraints on visualization of the data
and presentation of results to the scientific community. The motivation behind development of scSVA is to make visualization and exploration of data
quickly and easily. scSVA doesn't require coding or artistic skills to create impressive figures showing and explaining your data.
As a back-end it uses ", a("VaeX", href="http://vaex.astro.rug.nl", target="_blank") ,", an optimized python library for data binning.
scSVA provides a support for fast retrieval of one-dimensional arrays from compressed big martrices reducing the memory resources needed to visualize
and explore the data. scSVA is an R package distributed as-is, without warranties of any kind.
In case of any issues contact the author by email or use the issue tracker. Learn more about scSVA in sections below.",

h3(tags$b("File Upload Guide")),

"scSVA accepts tabular data stored in Comma-Separated Value (CSV) or in Hierarchical Data Format
(",a("HDF5", href="https://www.hdfgroup.org", target="_blank"),") files.",
"The input data consists of X- and Y-coordinates of your favorite embedding (e.g. tSNE or FLE),
a metadata information on each cell's cluster assignment and gene expression matrix or vectors that contain a gene quantification for each cell.",
"The file loading depends on a file type and is described separately for each file type.",

"The first thing to do is to open a File tab panel. The left sidebar allows for loading all necessary files for scSVA.",

h4(tags$b("Selection of XY coordinates")),

"First, select the \"Input File Format\" and then click on \"File Select\". This will open a new dialog box with the access to the file system.",
"Select the file that contains XY coordinates for your plot and click Select button in the dialog box.",
"The output depends on the selected file type:",

tags$ul(
  tags$li("HDF5 - A new table appears on the right-hand side of the side bar. An HDF5 file includes two types of objects:
          datasets and groups that store datasets and other groups. The table lists both groups and datasets with detailed information like
          dataset/group name, storage class, dimension and size of datasets. Select the datasets with X and Y cordinates in the left sidebar.
          scSVA supports HDF5-based ", a("Loom", href="https://github.com/linnarsson-lab/loompy", target="_blank"), "files.
          The filename extension must be from the following set [.h5,.hdf5,.loom] or from the same set but with all letters capitalized.
          "),
  tags$li("CSV - When you load a CSV file you make sure that the file has a header with column names.
          The column names (feature names) will appear in the new table on the right hand-side of the sidebar.
          scSVA uses \"fread\" from ",a("data.table", href="https://cran.r-project.org/web/packages/data.table/", target="_blank")," package
          which reads \"regular delimited files\" i.e. each row in a file has the same number of columns.
          The field delimiter is recognized automatically and users should use one from the folowing set [\\t ,;:|].
          The above field delimiters should not be used in column and row names.
          scSVA accepts both raw and compressed files (compressed with gzip) and the file type is recognized automatically.
          The file name extension must be from the following list ['.csv','.tsv','.txt'], or with all letters capitalized or compressed equivalents (.gz).
          ")
),

h4(tags$b("Selection of Expression Matrix")),

tags$ul(
  tags$li("HDF5 - scSVA accepts two dimensional datasets when a \"Dataset\" is selected or vectors of gene expression values when a \"Group\" is selected.
                If the \"Dataset\" is selected, specify if the rows or columns of the dataset are genes.
                If the \"Group\" is selected, scSVA will allow for loading expression vectors (one-dimensional datasets) within a group,
                with gene symbols as dataset's names.
          "),
  tags$li("CSV - We support only gene expression matrices, where the rows are genes and the columns are cells, and compressed with ",
          a("Gzip", href="https://www.gnu.org/software/gzip/", target="_blank"),
          ". This layout allows for indexing
          and very fast retrieval of expression vectors.
          The first column must contain gene names with \"Genes\" as the column name.
          The first row contains Cell IDs. Below is an example of an expression matrix:",
          tableHTML(data.frame(Genes=c("Gene.1","Gene.2", "Gene.3","..."),
                               "Cell.1"=c("0","2", "0","..."),
                               "Cell.2"=c("1","0", "1","..."),
                               'Cell.3'=c("0","3", "0","..."),
                               '...'=c("...","...", "...","...")),
                    rownames = F,widths = c(rep(60,5)),
                    theme= "default"),
          "The file saved as a CSV file with a comma as a delimiter:",
          tags$br(),
          "Genes,Cell.1,Cell.2,Cell.3,...", tags$br(),
          "Gene.1,0,1,0,...",tags$br(),
          "Gene.2,2,0,3,...",tags$br(),
          "Gene.3,0,1,0,...",tags$br(),
          "...",tags$br(),
          "To compress the expression matrix file use the following command in the terminal:",tags$br(),
          tags$code("gzip expression_matrix.csv"),tags$br(),
          "The next step is to upload a list with gene names that matches the gene names in the expression matrix.
          scSVA makes queries on the compressed expression matrix based on that list.
          The list can be saved at the preprocessing step or extracted from the uncompressed expression matrix file by the command:",tags$br(),
          tags$code("cat expression_matrix.csv |cut -d',' -f1 | awk '{if(NR>1) print} >Gene_Names.txt'"),tags$br(),
          "or from the compressed file:" ,tags$br(),
          tags$code("gunzip -cd expression_matrix.csv.gz |cut -d',' -f1 | awk '{if(NR>1) print} >Gene_Names.txt'"),tags$br(),
          "The index can be created on the compressed file by users using", a("Zindex", href="https://github.com/mattgodbolt/zindex", target="_blank"),
          " in the terminal",":",tags$br(),
          tags$code("zindex expression_matrix.csv.gz  --field 1  --tab-delimiter --skip-first 1 "),tags$br(),
          "for tab-delimited files or",tags$br(),
          tags$code("zindex expression_matrix.csv.gz  --field 1  --delimiter , --skip-first 1 "),tags$br(),
          "for e.g. comma-delimited files.",
          "scSVA can create the index for you if a tick box \"Create index on Compressed File\" is checked.",
          "Make sure that ",a("Zindex", href="https://github.com/mattgodbolt/zindex", target="_blank")," is properly installed on your computer and
          is visible as an executable file by your operating system.
          Indexing allows for very efficient retrieval of the rows from compressed files without need for loading of the entire matrix into the memory.
          Thus the memory usage is reduced by a factor of 10,000-40,000 for a typical expression matrix which makes possible a fast visualization
          of gene expression values on a 2D embedding or computing and plotting scores of a custom defined gene signature (set of genes)."
          )

),

h4(tags$b("Selection of Clusters")),

"Selection of a vector with cluster IDs is identical to the selection of XY coordinates. Cluster IDs must be an integer vector.",

h3(tags$b("Main Plot Apperance")),

"The tab \"View\" helps customizing the plot appearance. The defaul plot is a density plot where for a given grid resolution
it shows the number of cells within a square cell. The grid resolution can be modified by a
scroll bar below the plot. The default resolution is 512 which means that a plot consists of 512x512=262144 square cells.
Users can increase the resolution of the grid to produce higher quality figures but this results in increasing the plot refreshing time as more points needs to be plotted.
The tradeof between the figure quality and the interactive plot responsivness needs to be exploited by users as it depends on the number of cells and machine-specific factors.
We found that the default value of 512 is sufficient to produce high quality figures and the response is fast enough to interactively plot hundreds of millions of cells.",tags$br(),

"The two sidebars in the \"View\" tab simplify customization of the plot. Users can modify 1) size of the plotted points, 2) size and type of the font (x,y-axis labels, title, colorbar),
3) transparency of the points by specifying the alpha value, 4) color palette of the colorbar, 5) colors of the panel and plot background.
We included also a possibilty of checking the colorschemes in terms of suitability for color blind readers (",
a("dichromat package", href="https://cran.r-project.org/web/packages/dichromat/index.html", target="_blank"),"). Clicking on \"Dichromatism\" box
will open options for various types of colorblindness (green, red, blue) and users can simulte the effects of the type of dichromatism which helps selecting a color pallete suitable for
people with different types of color perception deficiency. We provide also color palletes suitable for dichromats from",
a("dichromat package", href="https://cran.r-project.org/web/packages/dichromat/index.html", target="_blank"),". These color palette names begin with \"dichromat::\".",

"The second sidbar modifies 1) the main title of the plot, 2) x- and y-axis label, 3) colorbar title.
Finaly, the main plot can be saved in many file formats with custom resolution and size.",tags$br(),

h4(tags$b("Plot Navigation")),

"Mouse left-click starts a selection of a rectangular region.",
"Zooming in: double-click inside the selected area will zoom in the selected region with a specified resolution (\"Grid Resolution\" scrollbar in the \"View\" tab).",
"\"Reset View\" button below the main plot when clicked will restore the original plot xy range.",

h3(tags$b("Analysis")),

"Analysis Menu consists of tools for exploratory analysis of uploaded datasets.",

h4(tags$b("Cell Abundance & Gene Expression")),

"Here users can explore number of cells in selected areas, plot gene expression profiles on 2d embedding,
and compute gene expression distributions and statistics like mean, median, standard deviation
in the selected areas.",tags$br(),

"Users can use two selection tools: 1) rectangular selection (default), and 2) polygonal selection.
Mouse left-click on the plot and holding down the left button lets you create a rectangular selection.
To start polygonal selection check \"Polygonal Selection\" box. From now, each mouse left-click on the plot
will create a new polygon's vertex and connect it by edges with the first and the penultimate vertex.
",tags$br(),

"To compute number of cells, select a region of interest and click \"Group A\" button.
This will mark cells in a selected area and return number of cells in that area.
Users can select two groups of cells. \"Reset Group\" button will unselect cells.
To plot gene expression values on 2D embedding, select a gene in right-hand panel.
By default the mean expression values are plotted on the grid but users can change the statistic by selecting \"Select a Statistic\".
Three options are supported: mean, median (the slowest), and maximal expression value.
The mean is suitable for highly expressed genes and the maximum value option works well for lowly expressed genes.
For a choosen gene, users can select two groups of cells as described above and compute the gene expression statistics in these two groups by clicking on
\"Get Group Statistic\" button.",

h4(tags$b("Gene Expression by Clusters")),

"Checking box \"Compute Gene Expression by Cluster\" starts computation of various statiscs in clusters like number of cells, gene expression mean, median and sd.
Once completed, the statistics are listed in the table. The table is fully interactive: 1) clusters can be sorted by a statistic,
2) searched by statistic values or cluster IDs, 3) filtered by statistic values by clicking on boxes below the table.
Single-click on a row of the table (cluster) will show a gene expression distribution for the selected cluster.
Users can mark many clusters to compare the distributions. Each selection opens a new distribution plot on right-hand side.
The plots are interactive and colored by the same colors as clusters on the 2D embedding. The scroll bar, on the top of the plots, changes the y-axis range for all plots.
Moving a coursor to the plot bar will show fraction of cells falling in a given expression range.
To compare all the distributions on one plot check \"Compute ECDF\" box.
This opens a new plot with empirical cumulative distribution functions for each selected cluster. Users can easily read dropout levels from the ecdf plot or
check for bi- or multimodalities that may suggest underclustering.",

h4(tags$b("Gene Signatures")),

"In order to compute gene signature scores (e.g. from Molecular Signature Database",a("MSigDB", href="http://software.broadinstitute.org/gsea/msigdb", target="_blank"),
") for all cells and visualize them on 2D embedding, users need to provide a set of gene names to \"Enter Gene List\" field and
specify \"Gene Signature Name\". The gene expression values across all cells can be retriewed in parallel from an expression matrix by specifying the number of workers.
Click on \"Compute Gene Signature Score\" box to start retrieving and computing the scores. scSVA computes gene expression z-score across all cells and average over
the z-scores across the entire gene set.",

h3(tags$b("Metadata")),

h4(tags$b("Clusters")),

"Here, users can plot clusters on 2D embedding using various color palletes, modify colors of individual clusters, show and modify cluster labels.",tags$br(),

"The middle panel contains a table with coordinates of cluster labels and right-hand panel customizes the plot.
First check \"Show Clusters\". This will compute cluster positions on a grid.
\"Get Cluster Label Positions\" check box, when selected, computes the position of labels by taking average x- and y-coordinates for each cluster.
The table in the middle panel lists all clusters (x-, and y-coordinates, and cluster IDs). Users can modify the coordinates of labels or cluster IDs by double-click
on the value or label in the table. The plot is updated after any modification.
\"Repel Labels\" check box distribute labels in such a way that the labels don't overlap
(",a("ggrepel package", href="https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html", target="_blank"),"),
which helps in case of many cluster labels. The labels can be furter customized by removing the bounding box or changing font type, size, color, etc. in the \"Annotation\" tab.",tags$br(),

"Colors of clusters on 2D embedding can be modified using various color palletes. We included support for many ",
a("Crayola crayon colors", href="https://en.wikipedia.org/wiki/List_of_Crayola_crayon_colors", target="_blank")," which generate distinctive and aesthetically-pleasing color sets to mark clusters.
Clicking on Shuffle Colors will permute colors from a given set. Make sure that none of the clusters has a color which is similar to the background color.
\"Shuffle Cells\" will permute the cells on the grid (important if some clusters overlap). The individual cluster colors can be changed in two ways:
1) checking \"Show color List\" will open a table will cluster colors. Double-click on the column element will allow for color modification for the selected cluster.
You can provide a color in hex format or by giving a color name e.g. \"navyblue\". 2) Specify the \"Cluster ID\" of the cluster color you want to modify
(the cluster can be highlighted by checking the \"Highlight Cluster\" box). Clicking \"Select Color for Cluster\" field will open a color picker.
Click on a color and the cluster color will change immediately.",

h3(tags$b("Annotation")),

"To start annotating the main figure check the \"Annotate\" box first. From now, any single mouse left-click on the main figure will create a new label.
The coordinates and text labels appear in the table in the middle panel. The table is interactive and double click on its element allow
for modification of this element. To remove an annotation label specify the \"Annotation ID\" you wish to remove and click on \"Remove Label\".
The labels can be further customized by specyfing the background color of the label, its transparency, size, and font type of the text.
The bounding box can be modified by changing the line width of its border or the radius of rounded corners.",

h3(tags$b("Extensions")),

h4(tags$b("Color Palette")),

"Users can customize their plots by creating their own color palletes. Click on \"Select Colors For a New Palette\".
This will open a color picker widget that can be used to select colors for a custom palette (the colors can be added to a palette also by specifying the hex color code).
Specifying the name of the color palette and clicking on \"Add New Palette\" will make it available in both \"Color Palette\" selection fields in \"Visualize\" and \"Metadata\" -> \"Clusters\" tabs.
The new custom color palette will appear at the end of the list of color palettes.",

h4(tags$b("Fonts")),

"To further customize plots, users can upload their own or open source fonts. scSVA accepts
TrueType fonts (filename extension .ttf). This functionality is available on Mac or Linux operating systems.
Windows users may consider installing scSVA as a docker container.
The default installation path of fonts on linux is /usr/share/fonts,
users will be asked to provide a password in the terminal to write the new fonts to this directory.
In some cases, to see new fonts on the interactive or saved figure, it may be necessary to restart the current R session.
For example, write in RStudio console: .rs.restartR() or use Command/Ctrl + Shift + F10 to restart the R session.
We encourage users to import new fonts at the first step and restart R session before loading files with x-, and y- coordinates.",

br(),

h4("Guided Video Tutorials:"),
"Created by Marcin Tabaka" ,
br(),
"Contact me at",
a("mtabaka@broadinstitute.org", href="mailto:mtabaka@broadinstitute.org"),
br(),
br(),
br(),
br()),
column(width=2)),

                 tabPanel("About",
tags$div(img(src = "scsva.logo.png", width = "400px", height = "150px",align="left")),
br(),
br(),
br(),
br(),
br(),
br(),
br(),
br(),
br(),
tags$div(img(src = "word_cloud.png", width = "450px", height = "450px",align="right")),
br(),
h4(
tags$ul(
tags$li("A graphical interface for fast visualization and exploratory analysis of single-cell data"),
br(),
tags$li("An interactive R package built with Shiny and VaeX as the back-end that allows for visualizing of a billion cells"),
br(),
tags$li("Applicable to any single-cell high-throughput genomic studies"),
br(),
tags$li("scSVA allows basic statistical analysis and help formulating hypotheses"),
br(),
tags$li("scSVA simplifies a production of high-quality figures for scientific publications"),
br(),
tags$li("scSVA allows for comprehensive figure customization by letting users to create their own color palletes or upload users' favorite fonts"),
br(),
tags$li("scSVA is optimized for efficient extracting of cell features
        from compressed big datasets in hdf5 and text file formats"),
br(),
tags$li("scSVA package can be installed from source or as a cloud-compatible image/container."),
br(),
tags$li("scSVA is open source, released under an MIT license, and is freely available for both commercial and academic research"))
),
br(),
br(),
"Developed by Marcin Tabaka, member of the ", a("Regev Lab", href="https://www.broadinstitute.org/regev-lab", target="_blank"),
", ", a("Broad Institute of MIT and Harvard", href="https://www.broadinstitute.org", target="_blank"),
br(),
"Contact me at",
a("mtabaka@broadinstitute.org", href="mailto:mtabaka@broadinstitute.org"),
br(),
h4("Acknowledgments:"),
br(),
h4("Citation:")
)
)
)
