ui <- fluidPage(theme = shinytheme("simplex"),
  navbarPage("scSVA: single cell Scalable Visualization and Analytics",
  tabPanel("File",
                 fluidRow(
                   column(width=4,
                     wellPanel(
                       h4("Select XY Coordinates"),
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
                       actionButton("load_dataset_XY", "Load Datasets")
                       ),
                       DT::dataTableOutput("file_cont_XY",       height = 275)
                       ),
                   column(width=4,
                          wellPanel(
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
                          actionButton("load_dataset_GE", "Load Datasets")
                          ),
                          DT::dataTableOutput("file_cont_Genes",    height = 275)
                   ),
                   column(width=4,
                          wellPanel(
                          h4("Select Groups"),
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
                            label    = "Groups",
                            choices  = c("NULL"),
                            selected = "NULL"),
                            actionButton("load_dataset_CL", "Load Datasets")),
                          DT::dataTableOutput("file_cont_Clusters", height = 275))
                   )
                 ),
#View
navbarMenu("Visualize",
  tabPanel("2D",
                   fluidRow(
                      column(width = 6,
                         plotOutput("plot1",
                         height   = 700,
                         click    = "plot1_click",
                         dblclick = "plot1_dblclick",
                         brush    = brushOpts(id = "plot1_brush",resetOnNew = TRUE)
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
                           "Color Bar Range:",
                           min   = 0,
                           max   = 1,
                           value = c(0,1),
                           step  = 0.1
                         )
                       ),
                       splitLayout(
                         actionButton("reset_view", "Refresh"),
                         actionButton("change_colorbar", "Apply")
                       )
                     ),
                       column(width=3,
                       offset=1,
                       wellPanel(style="max-height: 800px",
                         sliderInput("ps",
                           "Point Size",
                           min   = 0.01,
                           max   = 3,
                           value = 0.1,
                           step  = 0.01),
                         sliderInput("fs",
                           "Font Size",
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
                         selectizeInput(
                           inputId  = "fonts",
                           label    = "Font",
                           choices  = fonts(),
                           selected = NULL,
                           options=list(maxOptions = 3000)),
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
                             "Show Grid",
                             FALSE)),
                         splitLayout(
                           checkboxInput("showaxes",
                             "Show Axes",
                             TRUE),
                           checkboxInput("dens_scale",
                             "Scale Color Bar, log10",
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
                           "Color Bar Title",
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
                           "Units:",
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
  tabPanel("3D",
         fluidRow(
           column(width=6,
                  plotlyOutput("Plot3D",
                               height = 700,width = 700),
                  HTML("<br>"),
                  splitLayout(
                    sliderInput("grid_plot3D",
                                "Grid Resolution:",
                                min   = 100,
                                max   = 512,
                                value = 200,
                                step  = 1
                    ),
                    sliderInput("ColorbarRange3D",
                                "Color Bar Range:",
                                min   = 0,
                                max   = 1,
                                value = c(0,1),
                                step  = 0.1
                    )
                  ),
                  splitLayout(
                    actionButton("reset_view_3d", "Apply"),
                    actionButton("change_colorbar_3d", "Apply")
                  ),
                  DT::dataTableOutput("file_cont_XYZ",height = 400)
           ),
           column(width=2,
                  wellPanel(style="min-height: 800px",
                            h4("Select XYZ Coordinates"),
                            selectInput(
                              inputId  =  "inputFileType_XYZ",
                              label    =  "Input File Format",
                              choices  =  c("HDF5", "CSV"),
                              selected =  "HDF5"),
                            shinyFilesButton(
                              id       = 'inputFile_XYZ',
                              label    =  'File Select',
                              title    =  'Select a file',
                              multiple =  FALSE),
                            hr(),
                            style='border:32px;',
                            selectInput(
                              inputId  = "X_coord_3d",
                              label    = "X-coordinate",
                              choices  = c("NULL"),
                              selected = "NULL"),
                            selectInput(
                              inputId  = "Y_coord_3d",
                              label    = "Y-coordinate",
                              choices  = c("NULL"),
                              selected = "NULL"),
                            selectInput(
                              inputId  = "Z_coord_3d",
                              label    = "Z-coordinate",
                              choices  = c("NULL"),
                              selected = "NULL"),
                            actionButton("load_datasets_3d", "Load Datasets"),
                            hr(),
                            colourInput("panel_background_3d",
                                        label="Panel Background",
                                        "white"),
                            selectInput(
                              inputId  = "pal_3d",
                              label    = "Color Palette",
                              choices  = c("inferno", "viridis", "magma", "rainbow","matlab.like","matlab.like2","matlab.like3","gray.red"),
                              selected = "inferno"),
                            colourInput("axis_color_3d",
                                        label="Axis Color",
                                        "black"),
                            colourInput("font_color_3d",
                                        label="Axis Font Color",
                                        "black"),
                            colourInput("color_font_legend_3d",
                                        label="Legend Font Color",
                                        "black")
                  )),
           column(width=2,
                  wellPanel(style="min-height: 800px",
                            splitLayout(
                              checkboxInput("legend_3d", "Show Legend",value = TRUE),
                              checkboxInput("grid_3d", "Show Grid",value = TRUE)),
                            splitLayout(
                              checkboxInput("square_points_3d", "Square Points",value = FALSE),
                              checkboxInput("scale_colorbar_3d", "Scale Color Bar, log10",value = TRUE)
                            ),
                            hr(),
                            sliderInput("ps_3d",
                                        "Point Size",
                                        min   = 0.1,
                                        max   = 10,
                                        value = 2,
                                        step  = 0.1),
                            sliderInput("fs_3d",
                                        "Font Size",
                                        min   = 6,
                                        max   = 36,
                                        value = 14,
                                        step  = 1),
                            sliderInput("tick_3d",
                                        "Axis Font Size",
                                        min   = 6,
                                        max   = 36,
                                        value = 12,
                                        step  = 1),
                            sliderInput("points_transp_3d",
                                        "Transparency",
                                        min   = 0.0,
                                        max   = 1,
                                        value = 1.0,
                                        step  = 0.02),
                            hr(),
                            h5("Camera Viewpoint"),
                            sliderInput("cv_x",
                                        "X",
                                        min   = -5,
                                        max   = 5,
                                        value = 1.5,
                                        step  = 0.02),
                            sliderInput("cv_y",
                                        "Y",
                                        min   = -5,
                                        max   = 5,
                                        value = 1.5,
                                        step  = 0.02),
                            sliderInput("cv_z",
                                        "Z",
                                        min   = -5,
                                        max   = 5,
                                        value = 1.25,
                                        step  = 0.02)
                  )
           ),
           column(width=2,
                  wellPanel(
                  selectInput("display_feature",
                              "Display Feature",
                              choices = c("Density","Gene Expression/Signature")),
                  actionButton("project2d", "2D projection")),
                  bsModal("modal2Dprojection", "3D projection onto 2D", "project2d", size = "small",
                          tags$head(tags$style(HTML('
                                              .modal-sm {
                                                    width: 50%;
                                                    height: 1200px;
                                                    overflow-y: auto;
                                                    overflow-x: auto;
                                                    }
                                                    '))),
                          jqui_resizable(plotOutput("plot3dto2d",
                                                    height   = 600)),
                          splitLayout(
                          sliderInput("theta",
                                      "Theta:",
                                      value = 0,
                                      min = 1,
                                      max = 180),
                          sliderInput("phi",
                                      "Phi:",
                                      value = 0,
                                      min = 1,
                                      max = 360)
                          ),
                          actionButton("execute_project2d", "2D projection")
                          ),
                  wellPanel(
             textInput("title_3d",
                       "Title",
                       ""),
             textInput("x_axis_label_3d",
                       "X-axis Label",
                       "Dim 1"),
             textInput("y_axis_label_3d",
                       "Y-axis Label",
                       "Dim 2"),
             textInput("z_axis_label_3d",
                       "Z-axis Label",
                       "Dim 3"),
             textInput("colorbar_3d",
                       "Color Bar Title",
                       "Density"),
             hr(),
             selectInput("filetype_3d",
                         "Save As:",
                         choices = c("png", "pdf", "svg")),
             numericInput("width_3d",
                          "Width:", 800,
                          min  = 800,
                          max  = 5000,
                          step = 10),
             numericInput("height_3d",
                          "Height:", 800,
                          min  = 800,
                          max  = 5000,
                          step = 10),
             downloadButton("downloadPlot3D", "Download")
           )
           )
         ))
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
                           brush    = brushOpts(id = "plot1_brush",resetOnNew = TRUE)
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
                             label    = "Gene:",
                             choices  = "",
                             selected = "",
                             multiple = F),
                           checkboxInput("plot_expression",
                                         "Show Gene Expression",
                                         FALSE),
                           hr(),
                           selectInput(inputId  = "select_method",
                             label    = "Select Statistic:",
                             choices  = c("Mean","Median","Max"),
                             selected = "Mean",
                             multiple = F),
                           actionButton("get_statistic",
                             "Get Group Statistic")
                         )
                       )
                     )
                   ),
#Analysis:Gene Expression by Clusters
                   tabPanel("Gene Expression by Group",
                     fluidRow(
                       column(width = 4,
                         checkboxInput("get_exp_clusters",
                                      "Compute Gene Expression by Group",FALSE),
                         hr(),
                         DT::dataTableOutput('table_exp.in.clusters')
                         ),
                       column(width = 8,
                              tabsetPanel(type = "tabs",
                                          tabPanel("Density", value=1,plotlyOutput("plot_DistGEinCluster_Density",height = 600)),
                                          tabPanel("ECDF", plotlyOutput("plot_DistGEinCluster_ECDF",height = 600)),
                                          tabPanel("Violin", plotOutput("violin_plot",
                                                                        height   = 600)),
                                          tabPanel("Census", plotlyOutput("plot_DistGEinCluster_Census",height = 600)),
                                          id = "tabselected"
                              ),
                        splitLayout(
                      #   checkboxInput("rescale_y_axis",
                       #    "Rescale Y Axis",
                        #   value=FALSE),
                         conditionalPanel('input.tabselected==1',
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
                             "Number of Threads", 1,
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
navbarMenu("Cloud",
                      tabPanel("Storage",
                               fluidRow(  
                                 column(width=4,
                                        wellPanel(#style="min-height: 800px",
                                          selectInput("FS_type_1", "Filesystem: ", choices = c("local","gs"), selected = "local"),
                                          actionButton("FS_ls_1", "Explore FS",style='width:150px;'),
                                          DT::dataTableOutput("FS_1",       height = 500),
                                          HTML("<br>"),
                                          actionButton("FS_1_dirUP", label = icon("arrow-up")),
                                          actionButton("FS_1_dirRIGHT", label = icon("arrow-right")),
                                          HTML("<br><br>"),
                                          conditionalPanel('input.FS_type_1 == "local"',
                                                           splitLayout(
                                                             actionButton("FS_mkdir_1", "Create Directory",style='width:150px;'),
                                                             textInput('FS_mkdir_dirname_1',label=NULL,placeholder = "Dir Name")
                                                           ),
                                                           actionButton("FS_remove_local_1", "Remove Selected File",style='width:150px;')
                                          ),
                                          conditionalPanel('input.FS_type_1 == "gs"',
                                                           splitLayout(
                                                             actionButton("FS_mkdir_bucket_1", "Create New Bucket",style='width:180px;'),
                                                             textInput('FS_mkdir_bucket_name_1',label=NULL,value = "gs://",placeholder = "Dir Name")
                                                           ),
                                                           actionButton("FS_remove_object_1", "Remove Selected Object",style='width:180px;')
                                          )
                                        )
                                        
                                 ),
                                 column(width=1,
                                        HTML("<br><br><br><br><br><br><br><br><br><br><br>"),    
                                        HTML("Copy:<br>"),
                                        actionButton("FS_Right", label = icon("arrow-right")),
                                        actionButton("FS_Left", label = icon("arrow-left"))
                                 ),
                                 column(width=4,
                                        wellPanel(#style="min-height: 800px",
                                          selectInput("FS_type_2", "Filesystem: ", choices = c("local","gs"), selected = "local"),
                                          actionButton("FS_ls_2", "Explore FS",style='width:150px;'),
                                          DT::dataTableOutput("FS_2",       height = 500),
                                          HTML("<br>"),
                                          actionButton("FS_2_dirUP", label = icon("arrow-up")),
                                          actionButton("FS_2_dirRIGHT", label = icon("arrow-right")),
                                          HTML("<br><br>"),
                                          conditionalPanel('input.FS_type_2 == "local"',
                                                           splitLayout(
                                                             actionButton("FS_mkdir_2", "Create Directory",style='width:150px;'),
                                                             textInput('FS_mkdir_dirname_2',label=NULL,placeholder = "Dir Name")
                                                           ),
                                                           actionButton("FS_remove_local_2", "Remove Selected File",style='width:150px;')
                                          ),
                                          conditionalPanel('input.FS_type_2 == "gs"',
                                                           splitLayout(
                                                             actionButton("FS_mkdir_bucket_2", "Create New Bucket",style='width:180px;'),
                                                             textInput('FS_mkdir_bucket_name_2',label=NULL,placeholder = "Dir Name")
                                                           ),
                                                           actionButton("FS_remove_object_2", "Remove Selected Object",style='width:180px;')
                                          )                 
                                        )
                                 )
                               )
                      ),
                      
                      tabPanel("Compute Engine",
                               fluidRow(
                                 column(width=2,href="https://console.cloud.google.com",
                                        jqui_draggable(wellPanel(a("Open Google Cloud Console", href="https://console.cloud.google.com", target="_blank"))),
                                        jqui_draggable(wellPanel(#style="min-height: 800px",
                                          h5("GCE Setup"),
                                          textInput('GCE_AUTH_FILE', label='GCE Authorization File',placeholder = "PATH to json file"),
                                          textInput('GCE_DEFAULT_PROJECT_ID', label='GCE Project ID'),
                                          textInput('GCE_DEFAULT_ZONE', label='GCE Default Zone',value = "us-east1-b"),
                                          actionButton("gce_initialize", "GCE Initialize",style='width:150px;')
                                        )),
                                        jqui_draggable(wellPanel(#style="min-height: 800px",
                                          h5("GCE Info"),
                                          div(style="display:inline-block;width:10%;text-align: center;",actionButton("gce_list_disks", "GCE List Disks",style='width:150px;')),
                                          HTML("<br><br>"),
                                          div(style="display:inline-block;width:10%;text-align: center;",actionButton("gce_list_instances", "GCE List Instances",style='width:150px;'))
                                        ))
                                 ),
                                 column(width=2,
                                        jqui_draggable(wellPanel(#style="min-height: 800px",
                                          h5("RUN VM"),
                                          actionButton("gce_list_machinetypes", "GCE List Machine Types"),
                                          selectizeInput(inputId = "select_machine",
                                                         label    = "MachineTypes",
                                                         choices  = c("NA"),
                                                         selected = NULL),
                                          numericInput('disk_size', label='Disk Size (GB)', 10, min = 1, max = 1000),
                                          textInput('GCE_username', label='User Name'),
                                          textInput('GCE_password', label='Password Name'),
                                          textInput('container_name', label='Container Name'),
                                          selectizeInput(inputId = "select_template",
                                                         label    = "Template",
                                                         choices  = c("r-base"),
                                                         selected = "r-base"),
                                          textInput('vm_name', label='VM Name',placeholder = "(?:[a-z](?:[-a-z0-9]{0,61}[a-z0-9])?"),
                                          splitLayout(
                                            actionButton("run_vm", "Run VM")
                                          )
                                        ))
                                 ),
                                 column(width=2,
                                        jqui_draggable(wellPanel(
                                          actionButton("get_running_vms", "Get Running VMs"),
                                          HTML("<br><br>"),
                                          selectizeInput(inputId = "vm_name_sd",
                                                         label    = "Select VM",
                                                         choices  = NULL,
                                                         selected = NULL),
                                          actionButton("assign_vm", "Assign VM",style='width:100px;'),
                                          hr(),
                                          splitLayout(
                                            actionButton("stop_vm", "Stop VM  ",style='width:100px;'),
                                            actionButton("delete_vm", "Delete VM",style='width:100px;')),
                                          HTML("<br>"),
                                          splitLayout(
                                            actionButton("open_ssh", "Open SSH",style='width:100px;'),
                                            actionButton("show_processes", "Run Top  ",style='width:100px;')
                                          )
                                        )),
                                        jqui_draggable(wellPanel(
                                          h5("Local File System"),
                                          shinyFilesButton(
                                            id       =  'selectFiles_toVM',
                                            label    =  ' File Select ',
                                            title    =  ' File Select ',
                                            multiple =  TRUE),
                                          HTML("<br><br>"),
                                          splitLayout(                   
                                            actionButton("copyToVM", "Copy To VM",style='width:100px;'),
                                            textInput('dest_gcp',label=NULL,width = '100px',placeholder = "Insert PATH")
                                          ),
                                          hr(),
                                          h5("VM File System"),
                                          actionButton("selectFiles_fromVM", "File Select",style='width:100px;'),
                                          HTML("<br><br>"),
                                          splitLayout(
                                            actionButton("copyFromVM", "Copy To Local",style='width:100px;'),
                                            textInput('dest_local', label=NULL,width = '100px',placeholder = "Insert PATH")
                                          ),
                                          bsModal("modalVMFS", "VM File System:", "selectFiles_fromVM", size = "small",
                                                  dataTableOutput("dirTable"),
                                                  actionButton("dirUP", label = icon("arrow-up")),
                                                  actionButton("dirRIGHT", label = icon("arrow-right")),
                                                  dataTableOutput("fileTable"),
                                                  actionButton("selectFilesFromVM", "File Select")),
                                          actionButton("editselectefileFromVM", "Edit",style='width:100px;'),
                                          tags$head(tags$style(HTML('
                                                                    .modal-lg {
                                                                    width: 90%;
                                                                    max-height: 200vh;
                                                                    overflow-y: auto;
                                                                    overflow-x: auto;
                                                                    }
                                                                    '))),
                                          bsModal("modalEditFile", NULL, "editselectefileFromVM", size = "large",
                                                  jqui_resizable(sidebarPanel(width = 3,
                                                                              splitLayout(                          
                                                                                actionButton("EditSelectedFile", "Edit Selected File",style='width:140px;'),
                                                                                actionButton("SaveFile", "Save",style='width:140px;')
                                                                              ),
                                                                              HTML("<br><br>"),
                                                                              splitLayout(                   
                                                                                actionButton("SaveAs", "Save As",style='width:140px;'),
                                                                                textInput('SaveAsFile',label=NULL,width = '140px')
                                                                              ),
                                                                              HTML("<br>"),
                                                                              actionButton("CreateNew", "Create New",style='width:140px;'),
                                                                              hr(),
                                                                              selectInput("theme", "Theme: ", choices = themes, selected = "textmate"),
                                                                              selectInput("mode", "Mode: ", choices = modes, selected = "plain_text"),
                                                                              numericInput("size", "Tab Size:", 4),
                                                                              numericInput("FontSizeEditor", "Font Size:", 12),
                                                                              radioButtons("invisible", NULL, c("Hide Invisibles" = FALSE, "Show Invisibles" = TRUE), inline = TRUE)
                                                  )),
                                                  htmlOutput("EditedFile"),
                                                  #tags$style(HTML("
                                                  #ace {
                                                  #setShowPrintMargin:false;
                                                  #}")),
                                                  jqui_resizable(aceEditor(highlightActiveLine = T,
                                                                           height = '700px',
                                                                           vimKeyBinding = F,
                                                                           outputId = "ace", 
                                                                           value = ""
                                                  ))
                                                  
                                          )
                                          )),
                                        jqui_draggable(wellPanel(
                                          h5("RUN Script on VM"),
                                          h5("Selected script:"),
                                          verbatimTextOutput("selectedScript"),
                                          hr(),
                                          selectizeInput(inputId = "ContainerToRun",
                                                         label    = "Container to Run",
                                                         choices  = NULL,
                                                         selected = NULL),
                                          actionButton("runScriptOnVM", "Run Script",style='width:100px;')
                                        ))#,
                                        #    jqui_draggable(wellPanel(
                                        #      h5("Mount Bucket in VM"),
                                        #      actionButton("selectBuckettoVM", "Select Bucket",style='width:100px;'),
                                        #      actionButton("attachBuckettoVM", "Mount Bucket",style='width:100px;')
                                        #    ))
                                        ),
                                 column(width=6,
                                        wellPanel(
                                          htmlOutput("VM_assigned"),
                                          tags$style(HTML("
                                                          #terminalOutput {
                                                          height:600px;
                                                          overflow-y:scroll;
                                                          overflow-x:scroll;
                                                          background-color:#000000;
                                                          border-style:solid;
                                                          border-width: 1px;
                                                          border-color:#000000;
                                                          white-space:pre;
                                                          }
                                                          ")),
                                          tags$script(
                                            '
                                            Shiny.addCustomMessageHandler("scrollCallback",
                                            function(color) {
                                            var objDiv = document.getElementById("terminalOutput");
                                            objDiv.scrollTop = objDiv.scrollHeight;
                                            }
                                            );'
                    ),
                    
                    htmlOutput("terminalOutput")),
                    # DT::dataTableOutput("mytable_2",height = 700),
                    #verbatimTextOutput("summary", placeholder = TRUE),
                    textInput('command', label='Enter Command',width = '800px'),
                    splitLayout(
                      #actionButton("run_docker_command", "Run Docker Command"),
                      actionButton("run_ssh_command", "Run SSH Command"),
                      actionButton("clear_output", "Clear")
                    )
                                          )
                                 )),
           tabPanel("Diffusion Maps",
                    fluidRow(  
                      column(width=2,
                             wellPanel(
                               textInput('DMap_Path', 
                                         label='Path to Matrix',
                                         width = '800px'),
                               numericInput('DMap_nNN', 
                                            label='Number of NNs', 
                                            10, 
                                            min = 1, 
                                            max = Inf
                               ),
                               numericInput('DMap_k', 
                                            label='Number of Diffusion Components', 
                                            10, 
                                            min = 1, 
                                            max = Inf
                               ),
                               selectInput("DMap_EigDecompMethod",
                                           "Eigendecomposition Package",
                                           choices = c("Irlba","ARPACK")
                               ),
                               numericInput('DMap_ls', 
                                            label='Number of NNs to compute local sigma', 
                                            5, 
                                            min = 1, 
                                            max = Inf
                               ),
                               numericInput('DMap_nThreads', 
                                            label='Number of Threads', 
                                            1, 
                                            min = 1, 
                                            max = Inf
                               ),
                               selectInput("DMap_NNMethod",
                                           "ANNs method",
                                           choices = c("Annoy","Nmslib")
                                           ),
                               conditionalPanel(condition = "input.DMap_NNMethod == 'Annoy'",
                                                numericInput('DMap_Annoy_nTrees', 
                                                             label='Number of Trees', 
                                                             20, 
                                                             min = 1, 
                                                             max = Inf
                                                )),
                               conditionalPanel(condition = "input.DMap_NNMethod == 'Nmslib'",
                                                numericInput('DMap_Nmslib_M', 
                                                             label='M', 
                                                             10, 
                                                             min = 1, 
                                                             max = Inf
                                                ),
                                                numericInput('DMap_Nmslib_efC', 
                                                             label='efC', 
                                                             100, 
                                                             min = 1, 
                                                             max = Inf
                                                ),
                                                numericInput('DMap_Nmslib_efS', 
                                                             label='M', 
                                                             100, 
                                                             min = 1, 
                                                             max = Inf
                                                )
                                                ),
                                                actionButton("DMap_run", 
                                                             "Run DMaps",
                                                             style='width:100px;')
                               )
                             )
                      )
           ),
           tabPanel("Nearest Neighbor Graph",
                    fluidRow(  
                      column(width=2,
                             wellPanel(
                               textInput('NNG_Path', 
                                         label='Path to Matrix',
                                         width = '800px'),
                               numericInput('NNG_nNN', 
                                            label='Number of NNs', 
                                            10, 
                                            min = 1, 
                                            max = Inf
                               ),
                               numericInput('NNG_nThreads', 
                                            label='Number of Threads', 
                                            1, 
                                            min = 1, 
                                            max = Inf
                               ),
                               selectInput("NNG_NNMethod",
                                           "ANNs method",
                                           choices = c("Annoy","Nmslib")
                               ),
                               conditionalPanel(condition = "input.NNG_NNMethod == 'Annoy'",
                                                numericInput('NNG_Annoy_nTrees', 
                                                             label='Number of Trees', 
                                                             20, 
                                                             min = 1, 
                                                             max = Inf
                                                )),
                               conditionalPanel(condition = "input.NNG_NNMethod == 'Nmslib'",
                                                numericInput('NNG_Nmslib_M', 
                                                             label='M', 
                                                             10, 
                                                             min = 1, 
                                                             max = Inf
                                                ),
                                                numericInput('NNG_Nmslib_efC', 
                                                             label='efC', 
                                                             100, 
                                                             min = 1, 
                                                             max = Inf
                                                ),
                                                numericInput('NNG_Nmslib_efS', 
                                                             label='M', 
                                                             100, 
                                                             min = 1, 
                                                             max = Inf
                                                )
                               ),
                               actionButton("NNG_run", 
                                            "RUN NN Graph",
                                            style='width:100px;')
                             )
                      )
                    )
           ),
           tabPanel("3D FLE",
                    fluidRow(  
                      column(width=2,
                             wellPanel(
                               textInput('FLE_Path', 
                                         label='Path to Graph',
                                         width = '800px'),
                               selectInput("FLE_type",
                                           "FLE type",
                                           choices = c("fa_3d")),
                               numericInput('FLE_nsteps', 
                                            label='Number of Iterations', 
                                            1000, 
                                            min = 1, 
                                            max = Inf
                               ),
                               numericInput('FLE_nthreads', 
                                            label='Number of Threads', 
                                            1, 
                                            min = 1, 
                                            max = Inf
                               ),
                               numericInput('FLE_memmory', 
                                            label='Memmory (Gb)', 
                                            4, 
                                            min = 1, 
                                            max = Inf
                               ),
                               numericInput('FLE_scalingRatio', 
                                            label='Scaling Ratio', 
                                            1, 
                                            min = -Inf, 
                                            max = Inf
                               ),
                               numericInput('FLE_barnesHutTheta', 
                                            label='BarnesHut Theta', 
                                            1, 
                                            min = 1, 
                                            max = 10
                               ),
                               numericInput('FLE_barnesHutUpdateIter', 
                                            label='BarnesHut Update Every nth Iteration', 
                                            10, 
                                            min = 1, 
                                            max = Inf
                               ),
                               numericInput('FLE_barnesHutSplits', 
                                            label='BarnesHut Splits', 
                                            1, 
                                            min = 1, 
                                            max = 3
                               ),
                               checkboxInput("FLE_updateCenter",
                                             "Update Centers",
                                             FALSE),
                               checkboxInput("FLE_restart",
                                             "Restart",
                                             FALSE),
                               actionButton("FLE_run", 
                                            "Run",
                                            style='width:100px;'),
                               hr(),
                               actionButton("FLE_loadCoordinates", 
                                            "Load FLE Coordinates",
                                            style='width:150px;')
                             )),
                      column(width=5,
                             actionButton("FLE_getDistances", 
                                          "Get Distances",
                                          style='width:120px;'),     
                             plotlyOutput("plot_FLE_dist",height = 600),
                             actionButton("FLE_checkProgress", 
                                          "Check Progress",
                                          style='width:150px;'),
                             plotlyOutput("plot_FLE_checkProgress",height = 85,width = 500)
                      )
                    )         
           )
           ),
#Metadata:Groups
  navbarMenu("Metadata",
                   tabPanel("Groups",
                     fluidRow(
                       column(width = 6,
                         plotOutput("plot4",
                           height   = 700,
                           click    = "plot1_click",
                           dblclick = "plot1_dblclick",
                           brush    = brushOpts(id = "plot1_brush",resetOnNew = TRUE)
                         ),
                         actionButton("reset_view3", "Reset View")
                       ),
                       column(width = 3,
                         DT::dataTableOutput('table_clusters'),
                         HTML("<br><br>"),
                         conditionalPanel("input.plot_clusters == true",
                         wellPanel(
                         h5("Load Group Labels"),
                         selectInput(
                           inputId  =  "inputFileType_ClusterNames",
                           label    =  "Input File Format",
                           choices  =  c("HDF5", "CSV"),
                           selected =  "HDF5"),
                         shinyFilesButton(
                           id       = 'inputFile_ClusterNames',
                           label    =  'File Select',
                           title    =  'Load Group Labels',
                           multiple =  FALSE),
                         selectInput(
                           inputId  = "inputClusterNames",
                           label    = "Select Group Labels",
                           choices  = c("NULL"),
                           selected = "NULL"),
                       actionButton("load_dataset_ClusterNames", "Load")),
                       DT::dataTableOutput("file_cont_ClusterNames",       height = 275)
                       )),
                       column(width = 3, 
                        wellPanel(
                         checkboxInput("show_clusters",
                           "Show Groups",
                           FALSE),
                         conditionalPanel("input.show_clusters == true",
                         splitLayout(
                         checkboxInput("plot_clusters",
                           "Get Labels",
                           FALSE),
                         checkboxInput("show_label_clusters",
                                       "Show Labels",
                                       FALSE)),
                         splitLayout(
                         checkboxInput("repel_labels",
                           "Repel Labels",
                           FALSE),
                         checkboxInput("show_text",
                           "Remove Bounding Box",
                           FALSE)),
                         hr(),
                         checkboxInput("show_legend_clusters",
                                       "Show Legend",
                                       FALSE))),
                         conditionalPanel("input.show_legend_clusters == true",
                          wellPanel(                
                         textInput("legend_title_clusters",
                                   "Legend Title",
                                   ""),
                         selectInput(
                           inputId  = "legend_position",
                           label    = "Legend Position",
                           choices  = c("top","bottom","left","right"),
                           selected = "right"),
                         numericInput("legend_ncol_clusters",
                                      "Number of Columns", 1,
                                      min  = 1,
                                      max  = 10,
                                      step = 1),
                         checkboxInput("legend_byrow",
                                       "Group by Row",
                                       FALSE),
                         sliderInput("legend_point_size",
                                     "Legend Point Size",
                                     min   = 0.1,
                                     max   = 5,
                                     value = 1,
                                     step  = 0.1))),
                         wellPanel(
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
                           "Group ID:", 1,
                           min  = 1,
                           max  = Inf,
                           step = 1),
                         checkboxInput("highlight_cluster",
                           "Highlight Group",
                           FALSE),
                         colourInput("colCluster",
                           label="Change Group Color",
                           "white"),
                         checkboxInput("show_colors",
                           "Show Color List",
                           FALSE)
                         ),
                         conditionalPanel('input.show_colors == true',
                         DT::dataTableOutput('table_colors'))
                       )
                     )
                   ),
             tabPanel("Multiplots",
                      fluidRow(
                        column(width = 2,
                               wellPanel(
                               actionButton("collect_plots", "Add Plot"),
                               checkboxInput("show_multiplot",
                                             "Show Multiplots",
                                             FALSE)
                               ),
                               wellPanel(
                               sliderInput("multiplot_scale",
                                           "Scale:",
                                           min   = 0.1,
                                           max   = 1,
                                           value = 0.2,
                                           step  = 0.01
                               ),
                               numericInput("multiplot_ncols",
                                            "Number of Columns", 1,
                                            min  = 1,
                                            max  = Inf,
                                            step = 1),
                               hr(),
                               numericInput("multiplot_image_to_modify",
                                            "Image Number", 1,
                                            min  = 1,
                                            max  = Inf,
                                            step = 1),
                               actionButton("remove_image",
                                            "Remove Image"
                               ),
                               HTML("<br><br>"),
                               actionButton("remove_all_images",
                                            "Remove All Images"
                               )
                               )),
                        column(width = 10,
                               conditionalPanel('input.show_multiplot == true',
                               plotOutput("multiplot")
                               )
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
                         brush    = brushOpts(id = "plot1_brush",resetOnNew = TRUE)
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
                         selectizeInput(inputId = "annotation_fonts",
                           label    = "Font",
                           choices  = fonts(),
                           selected = NULL,
                           options=list(maxOptions = 3000)),
                         selectInput(inputId = "annotation_fontface",
                           label    = "Style",
                           choices  = c("plain","italic","bold"),
                           selected = "plain"),
                         sliderInput("annotation_labelSize",
                           "Stroke Weight",
                           min   = 0,
                           max   = 2,
                           value = 0.5,
                           step  = 0.05),
                         sliderInput("annotation_radiusCorners",
                           "Corner Radius",
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
                   tabPanel("Font",
                     sidebarPanel(
                       fileInput('file_fonts',
                                 'File',
                       accept = c("font/ttf")
                       ),
                       h5("Accepting TrueType fonts only (.ttf)")
                       )
                     )
                   ),
  tabPanel("Help",
    style = "font-size: 12pt; line-height: 20pt; color: black;",
    column(width=2),
      column(width=8,

tags$br(),
tags$br(),
tags$br(),
tags$br(),
tags$br(),


h3(tags$b("File Upload Guide")),

"scSVA accepts tabular data stored in Comma-Separated Value (CSV) or in Hierarchical Data Format
(",a("HDF5", href="https://www.hdfgroup.org", target="_blank"),") files.",
"The input data consists of (i) X- and Y-coordinates of your favorite embedding (e.g. tSNE or FLE),
(ii) gene expression matrix or vectors that contain a gene quantification for each cell, and (iii) metadata information on each cell's group assignment",
"The file loading depends on a file type and is described separately for each file type.",

"The first thing to do is to open a File tab panel. This panel allows for loading all necessary files to scSVA.",

h4(tags$b("Selection of XY Coordinates")),

"First, select the \"Input File Format\" and then click on \"File Select\". This will open a new dialog box with the access to the local file system.",
"Select the file that contains XY coordinates for your plot and click Select button in the dialog box.",
"The output depends on the selected file type:",

tags$ul(
  tags$li("HDF5 - A new table opens below the panel. An HDF5 file includes two types of objects:
          datasets and groups that store datasets and other groups. The table lists both groups and datasets with detailed information like
          dataset/group name, storage class, dimension and size of datasets. Select the datasets with X and Y coordinates.
          scSVA supports HDF5-based ", a("Loom", href="https://github.com/linnarsson-lab/loompy", target="_blank"), "files.
          The filename extension must be from the following set [.h5,.hdf5,h5ad,.loom] or from the same set but with all letters capitalized.
          "),
  tags$li("CSV - When you load a CSV file you make sure that the file has a header with column names.
          The column names (feature names) will appear in the new table below the panel.
          scSVA uses \"fread\" from ",a("data.table", href="https://cran.r-project.org/web/packages/data.table/", target="_blank")," package
          which reads \"regular delimited files\" i.e. each row in a file has the same number of columns.
          The field delimiter is recognized automatically and users should use one from the following set [\\t ,;:|].
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
                with gene symbols as dataset's names. This option is appropriate for very big datasets.
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
          Thus, the memory usage is reduced by a factor of 10,000-40,000 for a typical expression matrix which makes possible a fast visualization
          of gene expression values on a 2D embedding or computing and plotting scores of a custom defined gene signature (set of genes)."
          )

),

h4(tags$b("Selection of Groups")),

"Selection of a vector of group IDs is identical to the selection of XY coordinates. Group IDs must be a vector of integer (categorical data). ",

h3(tags$b("Plot Appearance")),

"The tab \"Visualize\" helps customizing the plot appearance. scSVA support visualization of 2D and 3D cell embeddings.",

h4(tags$b("2D")),

"The default plot is a density plot where for a given grid resolution
it shows the number of cells within a square. The grid resolution can be modified by a
scroll bar below the plot. The default resolution is 512 which means that a plot consists of 512x512=262144 square cells.
Users can increase the resolution of the grid to produce higher quality figures but this results in increasing the plot refreshing time as more points needs to be plotted.
The tradeoff between the figure quality and the interactive plot responsiveness needs to be exploited by users as it depends on the number of cells and machine-specific factors.
We found that the default value of 512 is sufficient to produce high quality figures and the response is fast enough to interactively plot hundreds of millions of cells.",tags$br(),

"The two sidebars in the \"2D\" tab simplify customization of the plot. Users can modify 1) size of the plotted points, 2) size and type of the font (x,y-axis labels, title, color bar),
3) transparency of the points by specifying the transparency value, 4) color palette of the colorbar, 5) colors of the panel and plot background.
We included also a possibility of checking the colorschemes in terms of suitability for color blind readers (",
a("dichromat package", href="https://cran.r-project.org/web/packages/dichromat/index.html", target="_blank"),"). Clicking on \"Dichromatism\" box
will open options for various types of colorblindness (green, red, blue) and users can simulate the effects of the type of dichromatism which helps selecting a color palette suitable for
people with different types of color perception deficiency. We provide also color palettes suitable for dichromats from",
a("dichromat package", href="https://cran.r-project.org/web/packages/dichromat/index.html", target="_blank"),". These color palette names begin with \"dichromat::\".",

"The second sidebar modifies 1) the main title of the plot, 2) x- and y-axis label, 3) color bar title.
Finally, the main plot can be saved in many file formats (vector or raster) with a custom resolution and size.",tags$br(),


h5(tags$b("Plot Navigation")),

"Mouse left-click starts a selection of a rectangular region.",
"Zooming in: double-click inside the selected area will zoom in the selected region with a specified resolution (\"Grid Resolution\" scrollbar).",
"\"Refresh\" button below the main plot when clicked will restore the original plot xy range.",

h4(tags$b("3D")),

"To initialize 3D plot, select X, Y, and Z coordinates and load datasets (see File section). 
Functions for 3D plot customization are similar to those present in 2D plots except the Camera Viewpoint 
which determines the camera view point about origin (", a("see manual", href="https://plot.ly/python/3d-camera-controls/", target="_blank"),
"). A panel \"Display Feature\" changes the feature to be plotted. Users can visualize a density of cells or gene/signature expression values.  
2D projection button opens a new window with 3D data projection onto 2D. 
",

h3(tags$b("Analysis")),

"Analysis Menu consists of tools for exploratory analysis of uploaded datasets.",

h4(tags$b("Cell Abundance & Gene Expression")),

"Here users can explore number of cells in selected areas, plot gene expression profiles on 2D embedding,
and compute gene expression distributions and statistics like mean, median, standard deviation
in the selected areas.",tags$br(),

"Users can use two selection tools: (i) rectangular selection (default), and (ii) polygonal selection.
Mouse left-click on the plot and holding down the left button lets you create a rectangular selection.
To start polygonal selection check \"Polygonal Selection\" box. From now, each mouse left-click on the plot
will create a new polygon's vertex and connect it by edges with the first and the penultimate vertex.
",tags$br(),

"To compute number of cells, select a region of interest and click \"Group A\" button.
This will mark cells in a selected area and return number of cells in that area.
Users can select two groups of cells. \"Reset Group\" button will unselect cells.
To plot gene expression values on 2D or 3D embedding, select a gene in right-hand panel.
By default, the mean expression values are plotted on the grid but users can change the statistic by selecting \"Select Statistic\".
Three options are supported: mean, median (the slowest), and maximal expression value.
The mean is suitable for highly expressed genes and the maximum value option works well for lowly expressed genes (e.g. transcription factors or surface receptors in scRNA-Seq).
For a chosen gene, users can select two groups of cells as described above and compute the gene expression statistics in these two groups by clicking on
\"Get Group Statistic\" button.",

h4(tags$b("Gene Expression by Group")),

"Checking box \"Compute Gene Expression by Group\" starts computation of statistics in groups like number of cells, gene expression mean, median and standard deviation.
Once completed, the statistics are listed in the table. The table is fully interactive: 1) groups can be sorted by a mean or number of cells for example,
2) searched by statistic values or group IDs, 3) filtered by statistic values by clicking on boxes below the table.
Single-click on a row of the table (group) will show a gene expression distribution for the selected group.
Users can mark many groups to compare the expression distributions among them. Each selection opens a new distribution plot on right-hand side.
The plots are interactive and colored by the same colors as groups on the embedding. The scroll bar, on the bottom of the plots, changes the y-axis range for all plots.
Moving a cursor to the plot bar will show fraction of cells falling in a given expression range.
To compare all the distributions on one plot go to \"ECDF\" tab to open a plot with empirical cumulative distribution function (users can easily read dropout levels from the ecdf plot or
check for bi- or multimodalities that may suggest underclustering) or \"Violin\" to visualize 
group expression values by violin plot. \"Census\" tab contains the fractions of cells in each group.",

h4(tags$b("Gene Signatures")),

"In order to compute gene signature scores (e.g. from Molecular Signature Database",a("MSigDB", href="http://software.broadinstitute.org/gsea/msigdb", target="_blank"),
") for all cells and visualize them on 2D or 3D embedding, users need to provide a set of gene names to \"Enter Gene List\" field and
specify \"Gene Signature Name\". The gene names in the list must be in the same format as in gene expression matrix (e.g. gene symbols).
The gene expression values across all cells can be retrieved in parallel from an expression matrix by specifying the number of threads on local machine.
Click on \"Compute Gene Signature Score\" box to start retrieving and computing the scores. scSVA computes gene expression z-score across all cells and average over
the z-scores across the entire gene set.",

h3(tags$b("Cloud")),

"scSVA supports interactive analytics and data storage mamagement on the cloud (Google Cloud Platform). 
Google Cloud Platform charges you for storage and running VMs, ", a("see pricing", href="https://cloud.google.com/pricing/", target="_blank"),".",
"Stop or delete any GCP services once you are done with your analysis. Good practice is to often check the billing of your project.
scSVA is an R package distributed as-is, without warranties of any kind. Make sure that you use the GCP services reasonably and follow best practices.",
"To get started make sure that your system configuration is ready to use GCP and in particular to launch virtual machine instances in Google Cloud Compute Engine (GCE).
Follow the instructions below to install and configure all required tools:",tags$br(),
"(i) Create and configure a Google Cloud Project, see ", a("instructions", href="https://cloud.google.com/resource-manager/docs/creating-managing-projects", target="_blank"),".",tags$br(),
"(ii) Download GCE private key in JSON format, see ", a("instructions", href="https://cloud.google.com/storage/docs/authentication#service_accounts", target="_blank"),".",tags$br(),
"(iii) Download and install ",a("Google Cloud SDK", href="https://cloud.google.com/sdk/downloads", target="_blank"),".",tags$br(),
"(iv) Initialize and authorize Google Cloud SDK, see ",a("instructions", href="https://cloud.google.com/sdk/docs/quickstarts", target="_blank"),".",tags$br(),
"See ",a("tutorial", href="https://cloud.google.com/getting-started/#quick-starts", target="_blank")," for more information.",
"scSVA uses ",a("googleComputeEngineR", href="https://cloudyr.github.io/googleComputeEngineR/", target="_blank"), ", which provides R interface to the GCE API.",



h4(tags$b("Storage")),

"To list google buckets from the default project, choose \"gs\" and click \"Explore FS\". Select a bucket from the list and click right arrow 
below the list to show its content. 
The upper arrow below the list navigates up one directory. To remove the object from the bucket, select it and click \"Remove Selected Object\". 
Users can also create new buckets by specyfying the name of the bucket and clicking \"Create New Bucket\" button. 
Selection of \"local\" filesystem will list directories and files on your local computer. Users can easily copy objects 
between the local computer and the bucket by clicking left/right arrow between filesystem panels.",

h4(tags$b("Compute Engine")),
"scSVA uses ", a("googleComputeEngineR", href="https://cloudyr.github.io/googleComputeEngineR/", target="_blank"), 
"to lunch and operate VMs. In the left corner, there is a link to Google Cloud Console. 
Clicking it opens the Console in new browser tab. Users can monitor there all services provided by GCP.

To get started, go to \"GCE Setup\" and provide the path to GCE private key in JSON format (see step (ii)), 
your project name (see step (i)), and the zone you want to launch VM (see ",
a("instructions", href="https://cloud.google.com/compute/docs/regions-zones/", target="_blank"),")",
", After successful configuration, \"GCE List Disks\" in \"GCE Info\" allows to display 
all disks in the project, and \"GCE List Instances\" list all VMs (stopped and running) 
in the project. Right console displays all the information from GCE",

"To launch a new VM, or restart a stopped VM, go to \"RUN VM\" window, and click
\"GCE List Machine Types\". This will list all available Machines in your selected zone.
Select container with computational tools you want to run on VM, and provide name for your virtual machine.
To build a new image on Google Cloud, go to terminal, then to directory with a Dockerfile and run the following command:", tags$br(),
tags$code("gcloud container builds submit --timeout=2h --tag gcr.io/Project_Name/ContainerName"),tags$br(),
"For further information, see instructions to ",a("googleComputeEngineR", href="https://cloudyr.github.io/googleComputeEngineR/", target="_blank"), "package.",
"Usually it takes less than one minute to lunch VM. You can check the progress in Rstudio console. 
If the VM instance is successfully launched you should see its name in the list of VM instances with status 
\"RUNNING\" after clicking the \"GCE List Instances\" button.",

"\"Get Running VMs\" button displays all VM instances with status \"RUNNING\". 
Specify the VM on which you want to run jobs and click \"Assign VM\". 
Notice the change in the name of \"Selected VM\" above the right console. 
From now, you can run commands and containerized tools and copy files between local machine and VM instance.
You can also stop or delete the VM instance and connect to VM instance through ssh. 
If you stop the VM instance you are charged only for storage, and the restarted VM instance will have access to all files from previous session.
If you delete VM instance, it will permanently remove a disk and VM instance (you will not be able to restart it in the future).  
In order to list running processes on the VM instance click \"Run Top\".
Users can easily transfer files from/to the VM instance using the panel below.
Selected files from the list can be modified after clicking \"Edit\" button. 
Users can also create new files on the VM instances. 
If your container has tools installed that can run scripts (eg. bash, R, or python) 
selecting a script on VM and clicking \"Run Script\" will start running the script on the VM instances.
If you have many running VM instances you can switch between them by selecting a running instance and assigning it.
We added also the option of running commands on the VM instance. Write the command in the test field below the console
(as you would do it in the bash terminal) and click \"Run SSH Command\".
The output of the command appears in the console. To remove all text from the console, click \"Clear\" button.
",

#h4(tags$b("Compute 3D FLE")),
#a("ForceAtlas2", href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0098679", target="_blank"),
#"number of iterations",
#"number of threads",
#"BarnesHut Theta",
#"BarnesHut Update Every nth Iteration",
#"BarnesHut Splits",
#"Update Centers",
#"Restart",
#"Get Distances",
#"Load FLE Coordinates",
#"Check Progress",

h3(tags$b("Metadata")),

h4(tags$b("Groups")),

"Here, users can visualize and annotate cell groups on 2D embedding by selecting color palettes, 
modify colors of individual groups, show and modify cluster labels.",tags$br(),

"The middle panel contains a table with coordinates of group labels and right-hand panel customizes the plot.
First check \"Show Groups\". This will compute group positions (cell groups for example) on a grid.
\"Get Labels\" check box, when selected, computes the position of labels by taking average x- and y-coordinates for each group.
The table in the middle panel lists all groups (group IDs, x-, and y-coordinates, and group labels). 
Users can modify the coordinates of labels plotted on 2D and group names by double-click
on the value or label in the table. The plot is updated after any modification.
\"Get Labels\" check box opens a panel in which users can upload group names from a file. 
The size (number of rows in case of csv file or number of elements in hdf5 dataset) must be equal to the total number of groups.  
\"Repel Labels\" check box distribute labels in such a way that the labels don't overlap
(",a("ggrepel package", href="https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html", target="_blank"),"),
which helps distributing the labels in such a way that they don't overlap. 
The labels can be further customized by removing the bounding box or changing font type, size, color, etc. in the \"Annotation\" tab.",tags$br(),
"To show a legend with group annotations check \"Show Legend\". Users can modify the legend title, position, number of columns, order of groups, 
and displayed point size in the panel that opens after checking the box \"Show Legend\"",

"Colors of cell groups on 2D embedding can be modified in \"Color Palette\" panel. We included support for most of ",
a("Crayola crayon colors", href="https://en.wikipedia.org/wiki/List_of_Crayola_crayon_colors", target="_blank")," 
which generate distinctive and aesthetically-pleasing color sets to mark groups of cells.
Clicking on Shuffle Colors will permute colors from a given set. Make sure that none of the groups has a color which is close to the background color.
\"Shuffle Cells\" will permute the cells on the grid (important if some groups overlap for example). The individual group colors can be changed in two ways:
(i) checking \"Show color List\" will open a table with group colors. Double-click on the column element will allow for color modification for the selected group.
Users can provide a color in hex format or by giving a color name e.g. \"navyblue\"; (ii) Specify the \"Group ID\" whose color you wish to modify
(the group can be highlighted by checking the \"Highlight Group\" box). Clicking \"Change Group Color\" field will open a color picker.
Click on a color and the color of the cell group will change immediately.",

h4(tags$b("Multiplots")),

"Here, users can arrange multiple graphs on a grid. To start click on \"Add Plot\" button. This will add the main plot to a list of plots.
Check the box \"Show Multiplots\" to make it visible. To add more graphs just click \"Add Plot\". To change the size of added plot or its resolution 
go to \"Visualaize\", \"2D\" tab and modify the plot properties on the right sidebar. scSVA uses ", 
a("R bindings", href="https://cran.r-project.org/web/packages/magick/index.html", target="_blank")," to",
a("ImageMagick", href="https://www.imagemagick.org/script/index.php", target="_blank"), "to arrange graphs on a grid.
Multiplots can be zoomed in/out by changing the \"Scale\" scroll bar, \"Number of Columns\" specify number of graphs per row.
",

h3(tags$b("Annotation")),

"To start annotating the main plot check the \"Annotate\" box. From now, any single mouse left-click on the main plot will create a new label.
The coordinates and text labels appear in the table in the middle panel. The table is interactive and double click on its element allow
for modification of this element. To remove an annotation label, specify the \"Annotation ID\" you wish to remove and click on \"Remove Label\".
The labels can be further customized by specifying the background color of the label, its transparency, size, and font type of the text.
The bounding box can be modified by changing the line width of its border or the radius of corners.",

h3(tags$b("Extensions")),

h4(tags$b("Color Palette")),

"Users can customize their plots by creating their own color palettes. Click on \"Select Colors For a New Palette\".
This will open a color picker widget that can be used to select colors for a custom palette (the colors can be added to a palette also by specifying the hex color code).
Specifying the name of the color palette and clicking on \"Add New Palette\" will make it available in both \"Color Palette\" selection fields in \"Visualize\" and \"Metadata\" -> \"Groups\" tabs.
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
"Created by Marcin Tabaka",
br(),
"Contact me at",
a("mtabaka@broadinstitute.org", href="mailto:mtabaka@broadinstitute.org"),
br(),
br(),
br(),
br()),
column(width=2)),

                 
tabPanel("About",

h3(tags$b("Motivation")),

h4("The main aim of the scSVA (link to github) is to simplify exploratory data analysis in the field of single-cell transcriptomics.
The progress in the development of high-throughput single-cell methods has allowed for study of cell types and states of
hundreds of thousands of cells. The interest of scientist in understanding biology at single cell level is not a fad but rather an ultimate goal -
e.g. the  ,",a("Human Cell Atlas Consortium", href="https://www.humancellatlas.org", target="_blank"), " is planning to profile at least 10 billion human cells
across all tissues. The massive growth of data volumes makes not only the data analysis harder but also put constraints on visualization of the data
and presentation of results to the scientific community. The motivation behind development of scSVA is to make visualization and exploration of data
quickly and easily. scSVA doesn't require coding or artistic skills to create impressive figures showing and explaining your data.
As a back-end it uses ", a("VaeX", href="http://vaex.astro.rug.nl", target="_blank") ,", an optimized python library for data binning.
scSVA provides a support for fast retrieval of one-dimensional arrays from compressed big matrices reducing the memory resources needed to visualize
and explore the data. scSVA is an R package distributed as-is, without warranties of any kind.
In case of any issues contact the author by email or use the issue tracker. Learn more about scSVA in sections below."),


br(),
br(),
br(),

h4(
tags$ul(
tags$li("A graphical interface for fast 2D & 3D visualization and exploratory analysis of single-cell data"),
br(),
tags$li("An interactive R package built with Shiny and VaeX as the back-end that allows for visualizing of a billion cells"),
br(),
tags$li("Applicable to any single-cell high-throughput genomic studies"),
br(),tags$div(img(src = "scsva.logo.png", width = "500px", height = "200px",align="right")),
tags$li("scSVA allows basic statistical analysis and help formulating hypotheses"),
br(),
tags$li("It simplifies a production of high-quality figures for scientific publications"),
br(),
tags$li("It allows for comprehensive figure customization by letting users to create their own color palettes or upload users' favorite fonts"),
br(),
tags$li("scSVA is optimized for efficient extracting of cell features
        from compressed big datasets in hdf5 and text file formats"),
br(),
tags$li("scSVA include fast algorithms for diffusion maps and 3D force-directed layout embedding (FLE)"),
br(),
tags$li("scSVA simplifies management of buckets and running scripts on the Google Cloud Platform"),
br(),
tags$li("scSVA package can be installed from source or as a cloud-compatible image/container."),
br(),
tags$li("scSVA is an open source, released under an MIT license, and is freely available for both commercial and academic research"))
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
