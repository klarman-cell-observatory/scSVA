server <- function(input, output, session) {
  vals <- reactiveValues(
    n.cols          =  20,
    X               =  c(py_to_r(ds$min('x')),py_to_r(ds$max('x'))),
    Y               =  c(py_to_r(ds$min('y')),py_to_r(ds$max('y'))),
    plot.ID         =  1,
    data            =  NULL,
    keeprows        =  NULL,
    GC_1            =  NULL,
    GC_2            =  NULL,
    GC_1.X          =  NULL,
    GC_1.Y          =  NULL,
    GC_2.X          =  NULL,
    GC_2.Y          =  NULL,
    trans           =  TRUE,
    p               =  NULL,
    Annotation.data =  data.frame(i = 0.0, j = 0.0,Annotation = NA,stringsAsFactors = F),
    f.annot         =  TRUE,
    group_Adist     =  NULL,
    group_Bdist     =  NULL,
    Ncells_Group_A  =  0,
    Ncells_Group_B  =  0,
    Summary_Stat    =  data.frame(Mean = 0.0, Median = 0.0,SD = 0.0,stringsAsFactors = F),
    gAa             =  FALSE,
    gBa             =  FALSE,
    rangeE          =  NULL,
    polygon         =  NULL,
    polygon_init    =  TRUE,
    polygon_groupA  =  NULL,
    polygon_groupB  =  NULL,
    polygon_gAa     =  FALSE,
    polygon_gBa     =  FALSE,
    cluster_labels  =  data.frame(i = 0.0, j = 0.0,ID = NA),
    cluster_df      =  data.frame(i = 0.0, j = 0.0,x = 0.0, Cluster_ID = 0 ),
    cluster_stat    =  data.frame(ID=0,Ncells=0,Mean=0.0,Median=0.0,SD=0.0, Q1=0.0, Q3=0.0,stringsAsFactors = F),
    loaded_feature  =  "",
    initialize      =  TRUE,
    cluster_colors  =  NULL,
    setup_colors    =  TRUE,
    N.clusters      =  NULL,
    cluster_h       =  NULL,
    CompClusFin     =  FALSE,
    color.df        =  NULL,
    df              =  NULL,
    GeneNames       =  NULL,
    modify_plot     =  NULL,
    inputFiles      =  list(XY=NULL,Genes=NULL,Clusters=NULL,Gene.Names=NULL,XYZ=NULL,Cluster.Names=NULL),
    EM              =  NULL,
    featureList     =  list(XY=NULL,Genes=NULL,Clusters=NULL,Gene.Names=NULL,XYZ=NULL,Cluster.Names=NULL),
    inputFiles_GCP  =  list(FILES=NULL),
    ClusterNames    =  NULL,
    ColPalette      =  NULL,
    InitColPalette  =  FALSE,
    Cols            =  FALSE,
    X_3d            =  NULL,
    Y_3d            =  NULL,
    Z_3d            =  NULL,
    trans_3d        = FALSE,
    data_3d         = NULL,
    p_3d            = NULL,
    load_3d         = FALSE,
    Cols3D          = FALSE,
    load_xy         = FALSE,  
    load_clust      = FALSE,
    FLE_dist        = FALSE,
    FLE_progress    = FALSE,
    gcp_dir                = "/",
    gcp_selected_dir       = "",
    gcp_files_indir        = "",
    gcp_selected_files     = "",
    local_path_edited_file = "",
    FS_ls_1         = data.frame(Content = NULL,stringsAsFactors = F),
    FS_ls_2         = data.frame(Content = NULL,stringsAsFactors = F),
    path_FS_1       = "/",
    path_FS_2       = "/",
    list_plots      = list()
  )
  
  #volumes = c(Computer = getwd())
  volumes = getVolumes()
  shinyFileChoose(input, 'inputFile_XY',               roots = volumes, filetypes=c('', 'h5','H5','hdf5','HDF5','h5ad','H5AD','loom','LOOM','csv','CSV','tsv','TSV','txt','TXT','gz'))
  shinyFileChoose(input, 'inputFile_Genes',            roots = volumes, filetypes=c('', 'h5','H5','hdf5','HDF5','h5ad','H5AD','loom','LOOM','gz'))
  shinyFileChoose(input, 'inputFile_Clusters',         roots = volumes, filetypes=c('', 'h5','H5','hdf5','HDF5','h5ad','H5AD','loom','LOOM','csv','CSV','tsv','TSV','txt','gz'))
  shinyFileChoose(input, 'inputFile_GenesVector',      roots = volumes, filetypes=c('', 'csv','CSV','tsv','TSV','txt','gz'))
  shinyFileChoose(input, 'inputFile_XYZ',              roots = volumes, filetypes=c('', 'h5','H5','hdf5','HDF5','h5ad','H5AD','loom','LOOM','csv','CSV','tsv','TSV','txt','TXT','gz'))
  shinyFileChoose(input, 'selectFiles_toVM',           roots = volumes)#, defaultPath = getwd())
  shinyFileChoose(input, 'inputFile_ClusterNames', roots = volumes, filetypes=c('', 'h5','H5','hdf5','HDF5','h5ad','H5AD','loom','LOOM','csv','CSV','tsv','TSV','txt','TXT','gz'))
  
  updateNumericInput(session = session,inputId = "n_workers",value = as.integer(availableCores()))
  updateSelectizeInput(session,"pal",choices=c(dens_cols))
  updateSelectizeInput(session,"pal_3d",choices=c(dens_cols))

  observeEvent(input$selectFiles_toVM, {
    vals$inputFiles_GCP$FILES  <- as.character(parseFilePaths(volumes, input$selectFiles_toVM)$datapath)
  })
  
  observeEvent(input$inputFile_ClusterNames, {
    vals$inputFiles$Cluster.Names  <- as.character(parseFilePaths(volumes, input$inputFile_ClusterNames)$datapath)
  })
  
  observeEvent(input$inputFile_XY, {
    vals$inputFiles$XY  <- as.character(parseFilePaths(volumes, input$inputFile_XY)$datapath)
  })
  
  observeEvent(input$inputFile_Genes, {
    vals$EM <- as.character(parseFilePaths(volumes, input$inputFile_Genes)$datapath)
  })

  observeEvent(input$inputFile_GenesVector, {
    vals$inputFiles$Gene.Names <- as.character(parseFilePaths(volumes, input$inputFile_GenesVector)$datapath)
  })

  observeEvent(input$inputFile_Clusters, {
    vals$inputFiles$Clusters <- as.character(parseFilePaths(volumes, input$inputFile_Clusters)$datapath)
  })
  
  observeEvent(input$inputFile_XYZ, {
    vals$inputFiles$XYZ  <- as.character(parseFilePaths(volumes, input$inputFile_XYZ)$datapath)
    print(vals$inputFiles$XYZ)
  })
  
  
   output$file_cont_XYZ <- renderDataTable({
    if(!is.null(vals$inputFiles$XYZ)){
      if(input$inputFileType_XYZ  ==  "HDF5"){
        choices  =  h5ls(vals$inputFiles$XYZ)
        ind      =  grep("H5I_DATASET",choices$otype)
        if(unique(choices$group[ind])=="/"){
          updateSelectizeInput(session,"X_coord_3d",choices=paste0(choices$group[ind],choices$name[ind]))
          updateSelectizeInput(session,"Y_coord_3d",choices=paste0(choices$group[ind],choices$name[ind]))
          updateSelectizeInput(session,"Z_coord_3d",choices=paste0(choices$group[ind],choices$name[ind]))
          } else {
          updateSelectizeInput(session,"X_coord_3d",choices=paste0(choices$group[ind],"/",choices$name[ind]))
          updateSelectizeInput(session,"Y_coord_3d",choices=paste0(choices$group[ind],"/",choices$name[ind]))
          updateSelectizeInput(session,"Z_coord_3d",choices=paste0(choices$group[ind],"/",choices$name[ind]))
        }
      DT::datatable(choices,options = list(scrollY = "150px",dom = 't'),rownames = FALSE,caption =paste0("File: ",vals$inputFiles$XYZ))
      }
      else if(input$inputFileType_XYZ  ==  "CSV"){
        if(summary(file(vals$inputFiles$XYZ))$class  ==  "gzfile"){
          vals$featureList$XYZ  <-  names(fread(paste0('gunzip -cq  ',gsub(" ","\\ ",vals$inputFiles$XYZ,fixed=T),' | head -n1'),nrows = 0))
          } else {
          vals$featureList$XYZ  <-  names(fread(input =vals$inputFiles$XYZ,nrows = 0))
          }
        updateSelectizeInput(session,"X_coord_3d",choices=vals$featureList$XYZ)
        updateSelectizeInput(session,"Y_coord_3d",choices=vals$featureList$XYZ)
        updateSelectizeInput(session,"Z_coord_3d",choices=vals$featureList$XYZ)
        
        DT::datatable(data.frame(Features=vals$featureList$XYZ),
                      options = list(scrollY = "150px",dom = 't'),
                      rownames = FALSE,caption =paste0("File: ",vals$inputFiles$XYZ))
      }
    }
  }) 
  output$file_cont_XY <- renderDataTable({
    if(!is.null(vals$inputFiles$XY)){
      if(input$inputFileType_XY  ==  "HDF5"){
        choices  =  h5ls(vals$inputFiles$XY)
        ind      =  grep("H5I_DATASET",choices$otype)
        if(unique(choices$group[ind])=="/"){
          updateSelectizeInput(session,"X_coord",choices=paste0(choices$group[ind],choices$name[ind]))
          updateSelectizeInput(session,"Y_coord",choices=paste0(choices$group[ind],choices$name[ind]))
          } else {
          updateSelectizeInput(session,"X_coord",choices=paste0(choices$group[ind],"/",choices$name[ind]))
          updateSelectizeInput(session,"Y_coord",choices=paste0(choices$group[ind],"/",choices$name[ind]))
        }
      DT::datatable(choices,options = list(scrollY = "150px",dom = 't'),rownames = FALSE,caption =paste0("File: ",vals$inputFiles$XY))
      }
      else if(input$inputFileType_XY  ==  "CSV"){
        if(summary(file(vals$inputFiles$XY))$class  ==  "gzfile"){
          vals$featureList$XY  <-  names(fread(paste0('gunzip -cq  ',gsub(" ","\\ ",vals$inputFiles$XY,fixed=T),' | head -n1'),nrows = 0))
          } else {
          vals$featureList$XY  <-  names(fread(input =vals$inputFiles$XY,nrows = 0))
          }
        updateSelectizeInput(session,"X_coord",choices=vals$featureList$XY)
        updateSelectizeInput(session,"Y_coord",choices=vals$featureList$XY)
        DT::datatable(data.frame(Features=vals$featureList$XY),
                      options = list(scrollY = "150px",dom = 't'),
                      rownames = FALSE,caption =paste0("File: ",vals$inputFiles$XY))
      }
    }
  })

  output$file_cont_Genes <- renderDataTable({
    if(!is.null(vals$EM)){
      if(input$inputFileType_Genes  ==  "HDF5"){
        choices  =  h5ls(vals$EM)

        if(input$h5_type  ==  "D"){
          ind  =  grep("H5I_DATASET",choices$otype)
          if(unique(choices$group[ind])  ==  "/"){
            updateSelectizeInput(session,"inputGenes",choices=paste0(choices$group[ind],choices$name[ind]))
            updateSelectizeInput(session,"inputGeneVector",choices=paste0(choices$group[ind],choices$name[ind]))
          } else {
            updateSelectizeInput(session,"inputGenes",choices=paste0(choices$group[ind],"/",choices$name[ind]))
            updateSelectizeInput(session,"inputGeneVector",choices=paste0(choices$group[ind],"/",choices$name[ind]))
          }
        }
        else if(input$h5_type  ==  "G"){
          un_choices=unique(choices$group)
          if(unique(choices$group)  ==  "/"){
            updateSelectizeInput(session,"inputGenes",choices=un_choices)
          } else {
            updateSelectizeInput(session,"inputGenes",choices=un_choices)
          }
        }

        DT::datatable(choices,
                      options  = list(scrollY = "150px",dom = 't'),
                      rownames = FALSE,
                      caption  = paste0("File: ",vals$EM))
      }
      else if(input$inputFileType_Genes  ==  "CSV"){
        DT::datatable(data.frame(data.frame(Files=c(paste0("Expression matrix: ",vals$EM),
                                                    paste0("Gene Names: ",vals$inputFiles$Gene.Names)))),
                      options  = list(scrollY = "150px",dom = 't'),
                      rownames = FALSE)
      }
    }
  })

  output$file_cont_Clusters <- renderDataTable({
    if(!is.null(vals$inputFiles$Clusters)){
      if(input$inputFileType_Clusters  ==  "HDF5"){
        choices  =  h5ls(vals$inputFiles$Clusters)
        ind      =  grep("H5I_DATASET",choices$otype)
        if(unique(choices$group[ind])  ==  "/"){
          updateSelectizeInput(session,"inputClusters",choices=paste0(choices$group[ind],choices$name[ind]))
        } else {
          updateSelectizeInput(session,"inputClusters",choices=paste0(choices$group[ind],"/",choices$name[ind]))
        }
        DT::datatable(choices,options = list(scrollY = "150px",dom = 't'),rownames = FALSE,caption =paste0("File: ",vals$inputFiles$Clusters))
      }
      else if(input$inputFileType_Clusters  ==  "CSV"){
        if(summary(file(vals$inputFiles$Clusters))$class  ==  "gzfile"){
          vals$featureList$Clusters<-names(fread(paste0('gunzip -cq  ',gsub(" ","\\ ",vals$inputFiles$Clusters,fixed=T),' | head -n1'),nrows = 0))
        } else {
          vals$featureList$Clusters  <-  names(fread(input =vals$inputFiles$Clusters,nrows = 0))
        }
        updateSelectizeInput(session,"inputClusters",choices=vals$featureList$Clusters)
        DT::datatable(data.frame(Features=vals$featureList$Clusters),
                      options   = list(scrollY = "150px",dom = 't'),
                      rownames  = FALSE,
                      caption   = paste0("File: ",vals$inputFiles$Clusters))
      }
    }
  })

output$file_cont_ClusterNames <- renderDataTable({
    if(!is.null(vals$inputFiles$Cluster.Names)){
      print(vals$inputFiles$Cluster.Names)
      if(input$inputFileType_ClusterNames  ==  "HDF5"){
        choices  =  h5ls(vals$inputFiles$Cluster.Names)
        ind      =  grep("H5I_DATASET",choices$otype)
        if(unique(choices$group[ind])  ==  "/"){
          updateSelectizeInput(session,"inputClusterNames",choices=paste0(choices$group[ind],choices$name[ind]))
        } else {
          updateSelectizeInput(session,"inputClusterNames",choices=paste0(choices$group[ind],"/",choices$name[ind]))
        }
        DT::datatable(choices,options = list(scrollY = "150px",dom = 't'),rownames = FALSE,caption =paste0("File: ",vals$inputFiles$Cluster.Names))
      }
      else if(input$inputFileType_ClusterNames  ==  "CSV"){
        if(summary(file(vals$inputFiles$Cluster.Names))$class  ==  "gzfile"){
          vals$featureList$Cluster.Names<-names(fread(paste0('gunzip -cq  ',gsub(" ","\\ ",vals$inputFiles$Cluster.Names,fixed=T),' | head -n1'),nrows = 0))
        } else {
          vals$featureList$Cluster.Names  <-  names(fread(input =vals$inputFiles$Cluster.Names,nrows = 0))
        }
        updateSelectizeInput(session,"inputClusterNames",choices=vals$featureList$Cluster.Names)
        DT::datatable(data.frame(Features=vals$featureList$Cluster.Names),
                      options   = list(scrollY = "150px",dom = 't'),
                      rownames  = FALSE,
                      caption   = paste0("File: ",vals$inputFiles$Cluster.Names))
      }
    }
  })
  
  
  observeEvent(input$colCluster, {
    vals$cluster_colors[input$cluster_id]  <-  input$colCluster
  })

  observeEvent(input$load_dataset_XY, {
#XY coordinates
    if(input$inputFileType_XY  ==  "HDF5"){
      withProgress(message = 'Loading X-coordinate',value = 0.25, {
      ds <<- vaex$from_arrays(x=numpy$array(h5read(vals$inputFiles$XY, input$X_coord)))})
      withProgress(message = 'Loading Y-coordinate',value = 0.5, {
      ds$add_column('y',data = numpy$array(h5read(vals$inputFiles$XY, input$Y_coord)))})
    }
    else if(input$inputFileType_XY  ==  "CSV"){
      withProgress(message = 'Loading X-coordinate',value = 0.25, {
      if(summary(file(vals$inputFiles$XY))$class  ==  "gzfile"){
        ds <<- vaex$from_arrays(x=numpy$array(fread(paste0('gunzip -cq  ',gsub(" ","\\ ",vals$inputFiles$XY,fixed=T)),
                                                    header = F,
                                                    select = c(which(vals$featureList$XY==input$X_coord)),skip=1)[[1]]))
      } else {
        ds <<- vaex$from_arrays(x=numpy$array(fread(vals$inputFiles$XY,
                                                    header = F,
                                                    select = c(which(vals$featureList$XY==input$X_coord)),skip=1)[[1]]))
      }
      })
      withProgress(message = 'Loading Y-coordinate',value = 0.25, {
      if(summary(file(vals$inputFiles$XY))$class  ==  "gzfile"){
        ds$add_column('y',data = numpy$array(fread(paste0('gunzip -cq  ',gsub(" ","\\ ",vals$inputFiles$XY,fixed=T)),
                                                    header = F,
                                                    select = c(which(vals$featureList$XY==input$Y_coord)),skip=1)[[1]]))
      } else {
        ds$add_column('y',data = numpy$array(fread(vals$inputFiles$XY,
                                                    header = F,
                                                    select = c(which(vals$featureList$XY==input$Y_coord)),skip=1)[[1]]))
      }})}
        gc()
        gc_python$collect()
        vals$X    <- c(py_to_r(ds$min('x')),py_to_r(ds$max('x')))
        vals$Y    <- c(py_to_r(ds$min('y')),py_to_r(ds$max('y')))
        vals$trans       =  TRUE
        vals$initialize  =  TRUE  
        vals$load_xy<-TRUE
        if(input$plot_expression){
        updateTextInput(session     = session,inputId = "choose_gene",    value ="")
        updateCheckboxInput(session = session,inputId = "plot_expression",value = FALSE)
        }
  })
  
observeEvent(input$load_dataset_GE, {  
#Cell Features
    vals$inputFiles$Genes<-vals$EM
    if(input$inputFileType_Genes  ==  "HDF5"){

      if(input$h5_type  ==  "D"){
        withProgress(message = 'Loading Gene List',value = 0.75, {
        vals$GeneNames<-h5read(vals$inputFiles$Genes, input$inputGeneVector)
        updateSelectizeInput(session,"choose_gene",choices=c("",vals$GeneNames),selected = "")})
      }
      else if(input$h5_type  ==  "G"){
        withProgress(message = 'Loading Gene List',value = 0.75, {
        choices=h5ls(vals$inputFiles$Genes)
        ind<-which(choices$group  ==  input$inputGenes & choices$otype  ==  "H5I_DATASET")
        updateSelectizeInput(session,"choose_gene",choices=c("",choices$name[ind]),selected = "")})
      }
    }
    else if(input$inputFileType_Genes  ==  "CSV"){
      withProgress(message = 'Loading Gene List',value = 0.75, {
        if(summary(file(vals$inputFiles$Gene.Names))$class  ==  "gzfile"){
          vals$GeneNames<-fread(paste0('gunzip -cq  ',gsub(" ","\\ ",vals$inputFiles$Gene.Names,fixed=T)),header=F)[[1]]
        } else {
          vals$GeneNames<-fread(vals$inputFiles$Gene.Names,header=F)[[1]]
        }
        updateSelectizeInput(session,"choose_gene",choices=c("",vals$GeneNames),selected = "")
      })
      if(input$computeIndex){
       withProgress(message = 'Creating an Index of the Selected Compressed File',value = 0.8, {
         delim=c(",","\t"," ","\\|",";",":")
         delim_in_matrix<-delim[which.max(sapply(1:length(delim),
                                                 function(x)  length(unlist(strsplit(
                                                                                     system(paste0('gunzip -cq ',
                                                                                                   gsub(" ",
                                                                                                        "\\ ",
                                                                                                        vals$inputFiles$Genes,
                                                                                                        fixed=T),
                                                                                                   ' | head -c 10000'),
                                                                                            intern = T),
                                                                                     delim[x])))))]
         if(delim_in_matrix  ==  "\t"){
          system(paste0("zindex ",gsub(" ","\\ ",vals$inputFiles$Genes,fixed=T)," --field 1  --tab-delimiter --skip-first 1"))
         }
         else{
          system(paste0("zindex ",gsub(" ","\\ ",vals$inputFiles$Genes,fixed=T)," --field 1 --delimiter ",delim_in_matrix," --skip-first 1"))
         }
      })
      }
    }
}) 
  
observeEvent(input$load_dataset_CL, {  
    if(input$inputFileType_Clusters  ==  "HDF5"){
      withProgress(message = 'Loading Clusters',value = 1, {
      ds$add_column('clusters',data = numpy$array(h5read(vals$inputFiles$Clusters, input$inputClusters)))})
    }
    else if(input$inputFileType_Clusters  ==  "CSV"){
      withProgress(message = 'Loading Clusters',value = 1, {
      if(summary(file(vals$inputFiles$Clusters))$class  ==  "gzfile"){
        ds$add_column('clusters',data = numpy$array(fread(paste0('gunzip -cq  ',gsub(" ","\\ ",vals$inputFiles$Clusters,fixed=T)),
                                                    header = F,
                                                    select = c(which(vals$featureList$Clusters==input$inputClusters)),skip=1)[[1]]))
      } else {
        ds$add_column('clusters',data = numpy$array(fread(vals$inputFiles$Clusters,
                                                    header = F,
                                                    select = c(which(vals$featureList$Clusters==input$inputClusters)),skip=1)[[1]]))
      }
      })
    }
   vals$load_clust<-TRUE
})
 
observeEvent(input$load_dataset_ClusterNames, {  
  if(input$inputFileType_ClusterNames  ==  "HDF5"){
    withProgress(message = 'Loading Cluster Names',value = 1, {
     vals$ClusterNames<- h5read(vals$inputFiles$Cluster.Names, input$inputClusterNames)})
  }
  else if(input$inputFileType_ClusterNames  ==  "CSV"){
    withProgress(message = 'Loading Cluster Names',value = 1, {
      if(summary(file(vals$inputFiles$Cluster.Names))$class  ==  "gzfile"){
        vals$ClusterNames<-fread(paste0('gunzip -cq  ',gsub(" ","\\ ",vals$inputFiles$Cluster.Names,fixed=T)),
                                                          header = F,
                                                          select = c(which(vals$featureList$Cluster.Names==input$inputClusterNames)),skip=1)[[1]]
      } else {
        vals$ClusterNames<-fread(vals$inputFiles$Clusters,
                                                          header = F,
                                                          select = c(which(vals$featureList$Cluster.Names==input$inputClusterNames)),skip=1)[[1]]
      }
    })
  }
  if(input$plot_clusters & input$show_clusters) {vals$cluster_labels$ID<-vals$ClusterNames}
})


   observeEvent(input$load_datasets_3d, {
#XYZ coordinates
    if(input$inputFileType_XYZ  ==  "HDF5"){
      tmp=h5read(vals$inputFiles$XYZ, input$X_coord_3d)
      if(ds$length()!=length(tmp)){
        if(input$show_clusters | input$plot_clusters | input$plot_expression) {
          updateCheckboxInput(session = session,inputId = "show_clusters",value = FALSE)
          updateCheckboxInput(session = session,inputId = "plot_clusters",value = FALSE)
          updateCheckboxInput(session = session,inputId = "plot_expression",value = FALSE)
        }
        ds <<- vaex$from_arrays(x=numpy$array(tmp))
        ds$add_virtual_column("x3d", "x")
        ds$add_column('y',data = numpy$array(h5read(vals$inputFiles$XYZ, input$Y_coord_3d)))
        ds$add_virtual_column("y3d", "y")
        ds$add_column('z3d',data = numpy$array(h5read(vals$inputFiles$XYZ, input$Z_coord_3d)))
      } else {
      withProgress(message = 'Loading X-coordinate',value = 0.25, {
      ds$add_column('x3d',data = numpy$array(tmp))})
      withProgress(message = 'Loading Y-coordinate',value = 0.5, {
      ds$add_column('y3d',data = numpy$array(h5read(vals$inputFiles$XYZ, input$Y_coord_3d)))})
      withProgress(message = 'Loading Z-coordinate',value = 0.75, {
      ds$add_column('z3d',data = numpy$array(h5read(vals$inputFiles$XYZ, input$Z_coord_3d)))})
      }
    }
    else if(input$inputFileType_XYZ  ==  "CSV"){
      withProgress(message = 'Loading X-coordinate',value = 0.25, {
      if(summary(file(vals$inputFiles$XYZ))$class  ==  "gzfile"){
        tmp=fread(paste0('gunzip -cq  ',gsub(" ","\\ ",vals$inputFiles$XYZ,fixed=T)),
                                                      header = F,
                                                      select = c(which(vals$featureList$XYZ==input$X_coord_3d)),skip=1)[[1]] 
        if(ds$length()!=length(tmp)){
          if(input$show_clusters | input$plot_clusters | input$plot_expression) {
            updateCheckboxInput(session = session,inputId = "show_clusters",value = FALSE)
            updateCheckboxInput(session = session,inputId = "plot_clusters",value = FALSE)
            updateCheckboxInput(session = session,inputId = "plot_expression",value = FALSE)
          }
          ds <<- vaex$from_arrays(x=numpy$array(tmp))
          ds$add_virtual_column("x3d", "x")
        } else {
          ds$add_column('x3d',data = numpy$array(tmp)) 
        }
      } else {
        tmp=fread(vals$inputFiles$XYZ,
                  header = F,
                  select = c(which(vals$featureList$XYZ==input$X_coord_3d)),skip=1)[[1]]
        if(ds$length()!=length(tmp)){
          if(input$show_clusters | input$plot_clusters | input$plot_expression) {
            updateCheckboxInput(session = session,inputId = "show_clusters",value = FALSE)
            updateCheckboxInput(session = session,inputId = "plot_clusters",value = FALSE)
            updateCheckboxInput(session = session,inputId = "plot_expression",value = FALSE)
          }
          ds <<- vaex$from_arrays(x=numpy$array(tmp))
          ds$add_virtual_column("x3d", "x")
        } else {
          ds$add_column('x3d',data = numpy$array(tmp)) 
        }
      }
      })
      withProgress(message = 'Loading Y-coordinate',value = 0.25, {
      if(summary(file(vals$inputFiles$XYZ))$class  ==  "gzfile"){
        ds$add_column('y',data = numpy$array(fread(paste0('gunzip -cq  ',gsub(" ","\\ ",vals$inputFiles$XYZ,fixed=T)),
                                                    header = F,
                                                    select = c(which(vals$featureList$XYZ==input$Y_coord_3d)),skip=1)[[1]]))
        ds$add_virtual_column("y3d", "y")
      } else {
        ds$add_column('y',data = numpy$array(fread(vals$inputFiles$XYZ,
                                                    header = F,
                                                    select = c(which(vals$featureList$XYZ==input$Y_coord_3d)),skip=1)[[1]]))
        ds$add_virtual_column("y3d", "y")
      }
      })
      withProgress(message = 'Loading Z-coordinate',value = 0.25, {
        if(summary(file(vals$inputFiles$XYZ))$class  ==  "gzfile"){
          ds$add_column('z3d',data = numpy$array(fread(paste0('gunzip -cq  ',gsub(" ","\\ ",vals$inputFiles$XYZ,fixed=T)),
                                                     header = F,
                                                     select = c(which(vals$featureList$XYZ==input$Z_coord_3d)),skip=1)[[1]]))
        } else {
          ds$add_column('z3d',data = numpy$array(fread(vals$inputFiles$XYZ,
                                                     header = F,
                                                     select = c(which(vals$featureList$XYZ==input$Z_coord_3d)),skip=1)[[1]]))
        }
      })
    }
    gc()
    gc_python$collect()
    vals$X_3d    <- c(py_to_r(ds$min('x3d')),py_to_r(ds$max('x3d')))
    vals$Y_3d    <- c(py_to_r(ds$min('y3d')),py_to_r(ds$max('y3d')))
    vals$Z_3d    <- c(py_to_r(ds$min('z3d')),py_to_r(ds$max('z3d')))
    
    vals$load_3d        =  TRUE
    vals$trans_3d       =  TRUE
    vals$plot.ID        =  1
    vals$trans           <- TRUE
    vals$initialize      <- TRUE
    vals$X               <- c(py_to_r(ds$min('x')),py_to_r(ds$max('x')))
    vals$Y               <- c(py_to_r(ds$min('y')),py_to_r(ds$max('y')))
    vals$group_Adist     =  NULL
    vals$group_Bdist     =  NULL
    vals$Ncells_Group_A  =  0
    vals$Ncells_Group_B  =  0
    vals$gAa             =  FALSE
    vals$gBa             =  FALSE
    vals$polygon         <- NULL
    Summary_Stat         =  data.frame(Mean = 0.0, Median = 0.0,SD = 0.0,stringsAsFactors = F)
    #if(input$plot_expression){
    #  updateTextInput(session     = session,inputId = "choose_gene",    value ="")
    #  updateCheckboxInput(session = session,inputId = "plot_expression",value = FALSE)
    #}

  })
  
  observeEvent(input$choose_gene, {
    if(input$choose_gene != ""){
      if(input$inputFileType_Genes  ==  "HDF5"){
        if(input$h5_type  ==  "D"){
          if(input$em_or  ==  "R"){
            ds$add_column('Gene',
                          data = numpy$array(drop(rhdf5::h5read(vals$inputFiles$Genes,
                                                                input$inputGenes,
                                                                index=list(which(vals$GeneNames==input$choose_gene),NULL)))))
          }
          else if(input$em_or  ==  "C"){
            ds$add_column('Gene',
                          data = numpy$array(drop(rhdf5::h5read(vals$inputFiles$Genes,
                                                                input$inputGenes,
                                                                index=list(NULL,which(vals$GeneNames==input$choose_gene))))))
          }
        }
        else if(input$h5_type=="G"){
          withProgress(message = paste0('Loading Expression Vector - ',input$choose_gene),value = 0.5, {
          ds$add_column('Gene',
                        data = numpy$array(h5read(vals$inputFiles$Genes,
                                                  input$choose_gene)))})
          }
      }
      else if(input$inputFileType_Genes=="CSV"){
        withProgress(message = paste0('Loading Expression Vector - ',input$choose_gene),value = 0.5, {
        ds$add_column('Gene',
                      data = numpy$array(unlist(fread(
                                                      paste0("zq ",
                                                             gsub(" ",
                                                                  "\\ ",
                                                                  vals$inputFiles$Genes,
                                                                  fixed=T),
                                                             " ",
                                                              input$choose_gene," |cut -f2-"),
                                                      stringsAsFactors = F))))
        })
      }
    updateCheckboxInput(session, "get_exp_clusters", value = FALSE)
    updateTextInput(session,     "title",            value = input$choose_gene)
    if(!is.null(vals$group_Adist) | !is.null(vals$group_Bdist)){vals$group_Adist    <- NULL
                                                                vals$group_Bdist    <- NULL
                                                                vals$gAa            =  FALSE
                                                                vals$gBa            =  FALSE
                                                                vals$Ncells_Group_A =  0
                                                                vals$Ncells_Group_B =  0
                                                                vals$Summary_Stat   <- data.frame(Mean   = 0.0,
                                                                                                  Median = 0.0,
                                                                                                  SD     = 0.0,
                                                                                                  stringsAsFactors = F)
                                                                }
    vals$rangeE          <-  c(py_to_r(ds$min('Gene')),py_to_r(ds$max('Gene')))
    vals$trans           =   TRUE
    vals$loaded_feature  <-  input$choose_gene
    if(input$plot_expression == FALSE) {
      updateTextInput(session,"plot_expression",value=TRUE)
    } else {
      updateTextInput(session,"clegend",value="Expression")}
    }
  })

  observeEvent(input$select_method, {
    vals$trans        =  TRUE
    vals$modify_plot  =  paste0("Selected method ",input$select_method)
  })

  update_slider<-function(){
    if(input$dens_scale){
      min_slider  =  log10(min(vals$data$x[vals$data$x>0]))
      max_slider  =  log10(max(vals$data$x[vals$data$x>0]))
    } else {
      min_slider  =  min(vals$data$x)
      max_slider  =  max(vals$data$x)
    }
    updateSliderInput(session,
                      "ColorbarRange",
                      value = round(c(min_slider,max_slider),3),
                      min   = round(min_slider,3),
                      max   = round(max_slider,3),
                      step  = round((max_slider-min_slider)/100,3))
    return(list(min=min_slider,max=max_slider))
  }

  update_slider_3d<-function(){
    if(input$scale_colorbar_3d){
      min_slider  =  log10(min(vals$data_3d$x[vals$data_3d$x>0]))
      max_slider  =  log10(max(vals$data_3d$x[vals$data_3d$x>0]))
    } else {
      min_slider  =  min(vals$data_3d$x)
      max_slider  =  max(vals$data_3d$x)
    }
    updateSliderInput(session,
                      "ColorbarRange3D",
                      value = round(c(min_slider,max_slider),3),
                      min   = round(min_slider,3),
                      max   = round(max_slider,3),
                      step  = round((max_slider-min_slider)/100,3))
    return(list(min=min_slider,max=max_slider))
  }
  
  observeEvent(input$dens_scale,{
    tmp       = update_slider()
    vals$Cols <- c(tmp$min,tmp$max)
  })

  observeEvent(input$change_colorbar,{
    vals$Cols  <-  input$ColorbarRange
  })

  plotInput <- reactive({
    isolate({if(vals$initialize){
        vals$X[1]       <-  vals$X[1]-(vals$X[2]-vals$X[1])*0.005
        vals$X[2]       <-  vals$X[2]+(vals$X[2]-vals$X[1])*0.005
        vals$Y[1]       <-  vals$Y[1]-(vals$Y[2]-vals$Y[1])*0.005
        vals$Y[2]       <-  vals$Y[2]+(vals$Y[2]-vals$Y[1])*0.005
        vals$initialize = FALSE
    }
    if(vals$trans){
      if(vals$plot.ID  ==  1){
        xycounts  =  ds$count(binby  =  tuple(ds$x,ds$y),
                              limits =  tuple(tuple(vals$X[1],vals$X[2]),tuple(vals$Y[1],vals$Y[2])),
                              shape  =  tuple(input$grid,input$grid),selection=FALSE)
        vals$data <- summary(as(py_to_r(xycounts),'dgTMatrix'))
      }
      else if(vals$plot.ID  ==  2){
        xycounts   =  ds$count(binby=tuple(ds$x,ds$y),
                            limits=tuple(tuple(vals$X[1],vals$X[2]),tuple(vals$Y[1],vals$Y[2])),
                            shape=tuple(input$grid,input$grid),selection=FALSE)
        allCells   <- summary(as(py_to_r(xycounts),'dgTMatrix'))
        if(input$select_method  ==  "Mean"){
        xycounts = ds$mean('Gene',
                           binby  =  tuple(ds$x,ds$y),
                           limits =  tuple(tuple(vals$X[1],vals$X[2]),tuple(vals$Y[1],vals$Y[2])),
                           shape  =  tuple(input$grid,input$grid),selection=FALSE)
        }
        else if(input$select_method=="Median"){
          xycounts = ds$median_approx('Gene',
                                      binby=tuple(ds$x,ds$y),
                                      limits=tuple(tuple(vals$X[1],vals$X[2]),tuple(vals$Y[1],vals$Y[2])),
                                      shape=tuple(input$grid,input$grid),selection=FALSE)
        }
        else if(input$select_method=="Max"){
          xycounts = ds$max('Gene',
                            binby  =  tuple('x','y'),
                            limits =  tuple(tuple(vals$X[1],vals$X[2]),tuple(vals$Y[1],vals$Y[2])),
                            shape  =  tuple(input$grid,input$grid),selection=FALSE)
        }
        vals$data<-summary(as(py_to_r(xycounts),'dgTMatrix'))
        if(input$select_method  ==  "Max") {
          vals$data$x[!is.finite(vals$data$x)]  <-  NA
        }
        vals$data  <-  vals$data[complete.cases(vals$data),]
        zeros      <-  suppressMessages(dplyr::anti_join(allCells[,c(1,2)],vals$data[,c(1,2)]))
        zeros$x    <-  rep(0,nrow(zeros))
        vals$data  <-  rbind(zeros,vals$data)
      }
      else if(vals$plot.ID  ==  3){
        vals$N.clusters  =  py_to_r(ds$max('clusters'))
        if(vals$setup_colors){
          vals$cluster_colors  <-  sample(rainbow(vals$N.clusters))
          vals$setup_colors    <-  FALSE
        }
        vals$cluster_df  <-  data.frame(i = 0.0, j = 0.0,x = 0, Cluster_ID = 0 )
        withProgress(message = 'Getting groups', value = 0, {
        for(l in 1:vals$N.clusters){
          ds$select(paste0("clusters==",l,""))
          cluster_grids<-py_to_r(ds$count('clusters',
                                          selection=T,
                                          binby   =  tuple(ds$x,ds$y),
                                          limits  =  tuple(tuple(vals$X[1],vals$X[2]),tuple(vals$Y[1],vals$Y[2])),
                                          shape   =  tuple(input$grid,input$grid)))
          cluster_grids   <-  summary(as(cluster_grids,'dgTMatrix'))
          cluster_grids   <-  cluster_grids[complete.cases(cluster_grids),]
          cluster_grids   <-  cbind(cluster_grids,Cluster_ID=rep(l,nrow(cluster_grids)))
          vals$cluster_df <- rbind(vals$cluster_df,cluster_grids)
          incProgress(1/vals$N.clusters)
        }
        })
        vals$CompClusFin   =  TRUE
        range_X            =  (vals$X[2]-vals$X[1])/input$grid
        range_Y            =  (vals$Y[2]-vals$Y[1])/input$grid
        vals$cluster_df$i  <- vals$X[1]+range_X*(vals$cluster_df$i-.5)
        vals$cluster_df$j  <- vals$Y[1]+range_Y*(vals$cluster_df$j-.5)
        vals$cluster_df    <- vals$cluster_df[-1,]
        rownames(vals$cluster_df) <- NULL
      }
      range_X       =  (vals$X[2]-vals$X[1])/input$grid
      range_Y       =  (vals$Y[2]-vals$Y[1])/input$grid
      vals$data$i   <- vals$X[1]+range_X*(vals$data$i-.5)
      vals$data$j   <- vals$Y[1]+range_Y*(vals$data$j-.5)
      vals$keeprows =  rep(TRUE, nrow(vals$data))
      vals$GC_1     =  rep(TRUE, nrow(vals$data))
      vals$GC_2     =  rep(TRUE, nrow(vals$data))
      vals$trans    =  FALSE
      tmp           =  update_slider()
      vals$Cols     <- c(tmp$min,tmp$max)
    }})

    if(input$plot_expression) {print(vals$modify_plot)}
      keep  <-  vals$data[ !vals$keeprows, , drop = FALSE]
      GC_1  <-  vals$data[ !vals$GC_1,     , drop = FALSE]
      GC_2  <-  vals$data[ !vals$GC_2,     , drop = FALSE]
    if(input$highlight_cluster){
      vals$cluster_h  <-  subset(vals$cluster_df,
                                 Cluster_ID  ==  input$cluster_id)
    }

    ggplot(vals$data, aes(x = i, y = j,colour = if(input$show_legend_clusters){as.factor(vals$cluster_df$Cluster_ID)}else{if(input$dens_scale){log10(x)} else{x}}))+
    {if(input$show_clusters & vals$CompClusFin) {
      if(input$show_legend_clusters){
      geom_point(data    =  vals$cluster_df,
                 size    =  input$ps,
                 alpha   =  input$transparency,
                 shape   =  ifelse(input$point_shape,15,19)#,
      )
      } else {
        geom_point(data    =  vals$cluster_df,
                   size    =  input$ps,
                   alpha   =  input$transparency,
                   shape   =  ifelse(input$point_shape,15,19),
                   colour  =  if(input$dichromat  !=  "No"){
                                  dichromat(vals$cluster_colors[vals$cluster_df$Cluster_ID],
                                            input$dichromat)
                              }
                              else{vals$cluster_colors[vals$cluster_df$Cluster_ID]})
        }
     }else{
      geom_point(size   =  input$ps,
                 alpha  =  input$transparency,
                 shape  =  ifelse(input$point_shape,15,19))
    }}+
    {if(input$highlight_cluster) geom_point(data    =  vals$cluster_h,
                                            size    =  input$ps*5,
                                            alpha   =  input$transparency,
                                            shape   =  ifelse(input$point_shape,15,19),
                                            colour  =  vals$cluster_colors[input$cluster_id])
    }+
    geom_point(data   =  GC_1,
               shape  =  ifelse(input$point_shape,15,19),
               fill   =  "#FF7A00",
               color  =  "#FF7A00",
               size   =  0.75)+
    geom_point(data   =  GC_2,
               shape  =  ifelse(input$point_shape,15,19),
               fill   =  "#353839",
               color  =  "#353839",
               size   =  0.75)+
    geom_point(data   = keep,
               shape  = ifelse(input$point_shape,15,19),
               fill   = "#FF7A00",
               color  = "#FF7A00",
               size   = 0.75)+
    {if(input$select_polygon & !is.null(vals$polygon)){
      geom_point(data  =  vals$polygon,
                 color =  "grey",
                 size  =  0.4)
    }}+
    {if(input$select_polygon & !is.null(vals$polygon)){
      geom_polygon(data  =  vals$polygon,
                   alpha =  0.1,
                   color =  "grey",
                   fill  =  "darkcyan",
                   size  =  0.4)
    }}+
    xlab(input$x_axis) +
    ylab(input$y_axis) +
    labs(title = input$title) +
    theme(plot.background    =  element_rect(fill  =  input$plot_background),
          legend.background  =  element_rect(fill  =  input$plot_background),
          panel.background   =  element_rect(fill  =  input$panel_background, colour = "black"),
          panel.grid.major   =  if(input$showgrid) {element_line(colour = "grey")} else {element_blank()},
          panel.grid.minor   =  element_blank(),
          axis.line          =  element_line(colour = "black",size = 0.5),
          text               =  element_text(family = input$fonts, face = "bold",size = input$fs),
          axis.ticks         =  element_line(size   = 0.5,colour = 'black'),
          axis.text          =  element_text(size   = input$fs,colour = 'black'))+
    {if(input$showlegend  ==  FALSE) theme(plot.background    =  element_rect(fill = input$plot_background),
                                           panel.background   =  element_rect(fill = input$panel_background),
                                           legend.background  =  element_rect(fill = input$plot_background),
                                           legend.position    =  "none",
                                           text               =  element_text(family  =  input$fonts,
                                                                              face    =  "bold",
                                                                              size    =  input$fs),
                                           legend.box         =  NULL,
                                           panel.grid.major   =  if(input$showgrid) {
                                                                   element_line(colour  =  "grey")
                                                                 } else {
                                                                   element_blank()
                                                                 })
    }+
    {if(!input$showaxes) theme(axis.line  = element_blank(),
                               axis.ticks = element_blank(),
                               axis.text  = element_blank(),
                               axis.title = element_blank())+
                         theme(plot.background  =  element_rect(fill    =  input$plot_background,
                                                                colour  =  input$plot_background),
                               panel.background =  element_rect(fill    =  input$panel_background,
                                                                colour  =  input$panel_background),
                               line             =  element_blank(),
                               rect             =  element_blank(),
                               legend.box       =  NULL,
                               panel.grid.major =  if(input$showgrid) {
                                                     element_line(colour = "grey")
                                                   } else {
                                                     element_blank()
                                                   })
    }+
   {if(input$annotate) geom_label(data             =  vals$Annotation.data, aes(x = i, y=j, label=Annotation),
                                  fill             =  input$annotation_background,
                                  inherit.aes      =  FALSE,
                                  size             =  input$annotation_size,
                                  family           =  input$annotation_fonts,
                                  fontface         =  input$annotation_fontface,
                                  show.legend      =  FALSE,
                                  alpha            =  input$annotation_transparency,
                                  label.size       =  input$annotation_labelSize,
                                  label.r          =  unit(input$annotation_radiusCorners, "lines"))
   }+
   {if(input$show_label_clusters){
      if(input$repel_labels){
        subset_labels  =  subset(vals$cluster_labels, i>=vals$X[1] & i<=vals$X[2] & j>=vals$Y[1] & j<=vals$Y[2])
          if(input$show_text) {
            geom_text_repel( data         =  subset_labels, aes(x = i, y=j, label=ID),
                             fill         =  input$annotation_background,
                             inherit.aes  =  FALSE,
                             size         =  input$annotation_size,
                             family       =  input$annotation_fonts,
                             fontface     =  input$annotation_fontface,
                             show.legend  =  FALSE,
                             alpha        =  input$annotation_transparency,
                             label.size   =  input$annotation_labelSize,
                             label.r      =  unit(input$annotation_radiusCorners, "lines"))
          } else {
            geom_label_repel(data         =  subset_labels, aes(x = i, y=j, label=ID),
                             fill         =  input$annotation_background,
                             inherit.aes  =  FALSE,
                             size         =  input$annotation_size,
                             family       =  input$annotation_fonts,
                             fontface     =  input$annotation_fontface,
                             show.legend  =  FALSE,
                             alpha        =  input$annotation_transparency,
                             label.size   =  input$annotation_labelSize,
                             label.r      =  unit(input$annotation_radiusCorners, "lines"))
          }
      } else {
      if(input$show_text) {
        geom_text(data         = vals$cluster_labels, aes(x = i, y=j, label=ID),
                  fill         = input$annotation_background,
                  inherit.aes  = FALSE,
                  size         = input$annotation_size,
                  family       = input$annotation_fonts,
                  fontface     = input$annotation_fontface,
                  show.legend  = FALSE,
                  alpha        = input$annotation_transparency,
                  label.size   = input$annotation_labelSize,
                  label.r      = unit(input$annotation_radiusCorners, "lines"))
          } else {
          geom_label(data         = vals$cluster_labels, aes(x = i, y=j, label=ID),
                     fill         = input$annotation_background,
                     inherit.aes  = FALSE,
                     size         = input$annotation_size,
                     family       = input$annotation_fonts,
                     fontface     = input$annotation_fontface,
                     show.legend  = FALSE,
                     alpha        = input$annotation_transparency,
                     label.size   = input$annotation_labelSize,
                     label.r      = unit(input$annotation_radiusCorners, "lines"))
          }
      }
    }}+
    {if(input$show_legend_clusters) {
      ind<-unique(vals$cluster_df$Cluster_ID)
      #subset_labels  = subset(vals$cluster_labels, i>=vals$X[1] & i<=vals$X[2] & j>=vals$Y[1] & j<=vals$Y[2])
      scale_colour_manual(name=input$legend_title_clusters,
                          values=if(input$dichromat  !=  "No"){
                                     dichromat(vals$cluster_colors,
                                              input$dichromat)
                                   }
                                 else{vals$cluster_colors[ind]},
                          labels=vals$cluster_labels$ID[ind])
    } else {scale_colour_gradientn(name     =  input$clegend,
                            limits   =  vals$Cols,
                            oob      =  squish,
                            colours  =  get_colors(input$pal,
                                                   input$dichromat  !=  "No",
                                                   input$dichromat,
                                                   vals$n.cols))}}+
    {if(input$show_legend_clusters) {guides(color = guide_legend(ncol  = input$legend_ncol_clusters,
                                                                 byrow = input$legend_byrow,
                                                                 override.aes = list(size=input$legend_point_size)))}}+
    {if(input$show_legend_clusters) {theme(legend.position = input$legend_position, legend.key = element_blank())}}+
    coord_cartesian(xlim = vals$X, ylim = vals$Y)
  })

  output$plot1 <- output$plot2 <- output$plot3 <- output$plot4 <- renderPlot(res = 72,{
  plotInput()
  })

  output$count_cells_Group_A <- renderPrint({
    cat(comma(vals$Ncells_Group_A, sep = "\n"))
  })

  output$count_cells_Group_B <- renderPrint({
    cat(comma(vals$Ncells_Group_B, sep = "\n"))
  })

  # Reset all points
  observeEvent(input$reset_view, {
    vals$trans           <- TRUE
    vals$initialize      <- TRUE
    vals$X               <- c(py_to_r(ds$min('x')),py_to_r(ds$max('x')))
    vals$Y               <- c(py_to_r(ds$min('y')),py_to_r(ds$max('y')))
    vals$group_Adist     =  NULL
    vals$group_Bdist     =  NULL
    vals$Ncells_Group_A  =  0
    vals$Ncells_Group_B  =  0
    vals$gAa             =  FALSE
    vals$gBa             =  FALSE
    vals$polygon         <- NULL
    gc_python$collect()
    Summary_Stat         =  data.frame(Mean = 0.0, Median = 0.0,SD = 0.0,stringsAsFactors = F)
  })

  observeEvent(input$reset_view2, {
    vals$trans           <- TRUE
    vals$initialize      <- TRUE
    vals$X               <- c(py_to_r(ds$min('x')),py_to_r(ds$max('x')))
    vals$Y               <- c(py_to_r(ds$min('y')),py_to_r(ds$max('y')))
    vals$group_Adist     =  NULL
    vals$group_Bdist     =  NULL
    vals$Ncells_Group_A  =  0
    vals$Ncells_Group_B  =  0
    vals$gAa             =  FALSE
    vals$gBa             =  FALSE
    vals$polygon_gAa     =  FALSE
    vals$polygon_gBa     =  FALSE
    polygon_groupA       =  NULL
    polygon_groupB       =  NULL
    vals$polygon         <- NULL
    gc_python$collect()
    vals$Summary_Stat    =  data.frame(Mean = 0.0, Median = 0.0,SD = 0.0,stringsAsFactors = F)
  })

   observeEvent(input$reset_view3, {
    vals$trans           <- TRUE
    vals$initialize      <- TRUE
    vals$X               <- c(py_to_r(ds$min('x')),py_to_r(ds$max('x')))
    vals$Y               <- c(py_to_r(ds$min('y')),py_to_r(ds$max('y')))
    vals$group_Adist     =  NULL
    vals$group_Bdist     =  NULL
    vals$Ncells_Group_A  =  0
    vals$Ncells_Group_B  =  0
    vals$gAa             =  FALSE
    vals$gBa             =  FALSE
    vals$polygon         <- NULL
    gc_python$collect()
    vals$Summary_Stat    =  data.frame(Mean = 0.0, Median = 0.0,SD = 0.0,stringsAsFactors = F)
  })

   observeEvent(input$reset_polygon,{
     vals$polygon        <-  NULL
   })

  observeEvent(input$reset_groups, {
    vals$GC_1            =  rep(TRUE, nrow(vals$data))
    vals$GC_2            =  rep(TRUE, nrow(vals$data))
    vals$GC_1.X          =  NULL
    vals$GC_1.Y          =  NULL
    vals$GC_2.X          =  NULL
    vals$GC_2.Y          =  NULL
    vals$group_Adist     =  NULL
    vals$group_Bdist     =  NULL
    vals$Ncells_Group_A  =  0
    vals$Ncells_Group_B  =  0
    vals$gAa             =  FALSE
    vals$gBa             =  FALSE
    vals$polygon         <- NULL
    vals$Summary_Stat    =  data.frame(Mean = 0.0, Median = 0.0,SD = 0.0,stringsAsFactors = F)
  })

  observeEvent(input$select, {
    res           <-  brushedPoints(vals$data, input$plot1_brush, allRows = TRUE)
    vals$keeprows <-  xor(vals$keeprows, res$selected_)
  })

  observeEvent(input$plot1_dblclick, {
    brush  <-  input$plot1_brush
    if (!is.null(brush)) {
      vals$X      <-  c(brush$xmin, brush$xmax)
      vals$Y      <-  c(brush$ymin, brush$ymax)
      vals$trans  <-  TRUE
    } else {
      vals$X      <-  NULL
      vals$Y      <-  NULL
    }
  })

  observeEvent(input$plot1_click, {
    if(input$annotate){
      vals$Annotation.data <- rbind(vals$Annotation.data,
                                    data.frame(i                = input$plot1_click$x,
                                               j                = input$plot1_click$y,
                                               Annotation       = "Add text",
                                               stringsAsFactors = F))
      if(vals$f.annot){
        vals$Annotation.data  <-  vals$Annotation.data[-1,]
        vals$f.annot          =   FALSE
      }
      n.rows                         =  nrow(vals$Annotation.data)
      rownames(vals$Annotation.data) <- 1:n.rows
    }
    else if(input$select_polygon){
      if(vals$polygon_init){
        vals$polygon      <- data.frame(i = input$plot1_click$x,
                                        j = input$plot1_click$y)
        vals$polygon_init = FALSE
      } else {
        vals$polygon  <-  rbind(vals$polygon,
                                data.frame(i  =  input$plot1_click$x,
                                           j  =  input$plot1_click$y))
      }
    }
  })

  observeEvent(input$select_polygon,{
    if(!input$select_polygon) {  vals$polygon  <-  NULL  }
  })

  observeEvent(input$remove_label, {
    if(input$annotate){
      if(nrow(vals$Annotation.data)  ==  1){
        vals$Annotation.data  =  data.frame(i                 =  0.0,
                                            j                 =  0.0,
                                            Annotation        =  NA,
                                            stringsAsFactors  =  F)
        n.rows                =  nrow(vals$Annotation.data)
        rownames(vals$Annotation.data)  <-  1:n.rows
      } else {
        vals$Annotation.data  <-  vals$Annotation.data[-input$row_id,]
        n.rows                =   nrow(vals$Annotation.data)
        rownames(vals$Annotation.data)  <-  1:n.rows
      }
    }
  })

  output$table_clusters <- DT::renderDataTable(datatable({plyr::rename(vals$cluster_labels,
                                                                       c("i"  =  input$x_axis,"j"  =  input$y_axis,"ID" = "Label"))},
                                               editable  =  TRUE,
                                               options   =  list(paging = FALSE,scrollY = "250px")) %>% formatRound(c(input$x_axis,input$y_axis), 3))

  proxy1 = DT::dataTableProxy('table_clusters')
  observeEvent(input$table_clusters_cell_edit, {
    info  =  input$table_clusters_cell_edit
    i     =  info$row
    j     =  info$col
    v     =  info$value
    vals$cluster_labels[i,j]  <<-  DT::coerceValue(v, vals$cluster_labels[i, j])
    DT::replaceData(proxy1, vals$cluster_labels, resetPaging = FALSE)
  })

  output$table1 <- DT::renderDataTable(datatable({plyr::rename(vals$Annotation.data,
                                                               c("i"  =  input$x_axis,"j"  =  input$y_axis, "Annotation" = "Label"))},
                                       editable  =  TRUE,
                                       options   =  list(pageLength = 100,scrollY = "250px")) %>% formatRound(c(input$x_axis,input$y_axis), 3))

  proxy2 = DT::dataTableProxy('table1')
  observeEvent(input$table1_cell_edit, {
    info = input$table1_cell_edit
    i    = info$row
    j    = info$col
    v    = info$value
    vals$Annotation.data[i,j] <<- DT::coerceValue(v, vals$Annotation.data[i, j])
    DT::replaceData(proxy2, vals$Annotation.data, resetPaging = FALSE)
  })

  output$table_colors <- DT::renderDataTable({
    if(input$show_colors & !input$plot_clusters){
      vals$color.df<-data.frame(ID    =  as.character(1:vals$N.clusters),
                                Color =  vals$cluster_colors[1:vals$N.clusters])
      datatable(vals$color.df,
                options  = list(paging   = FALSE,
                                scrollY  = "150px",
                                dom      = 't'),
                editable = TRUE)  %>% DT::formatStyle('ID', backgroundColor = styleEqual(vals$color.df$ID, vals$cluster_colors[1:vals$N.clusters]))
      } else if(input$show_colors & input$plot_clusters){
        vals$color.df<- data.frame(ID     =  as.character(vals$cluster_labels$ID),
                                   Color  =  vals$cluster_colors[1:vals$N.clusters],
                                   stringsAsFactors = F)
     datatable(vals$color.df,
               options   = list(paging = FALSE,
                                scrollY = "150px",
                                dom = 't'),
               editable  = TRUE)  %>%  DT::formatStyle('ID',backgroundColor = styleEqual(vals$color.df$ID, vals$cluster_colors[1:vals$N.clusters]))
      } else {}
  })

  proxy3 = DT::dataTableProxy('table_colors')
  observeEvent(input$table_colors_cell_edit, {
    info   =   input$table_colors_cell_edit
    i      =   info$row
    j      =   info$col
    v      =   info$value
    color.conv=function(x){
      tryCatch(rgb(t(col2rgb(x)),maxColorValue = 255),
               error = function(e) {rgb(t(col2rgb("black")),maxColorValue = 255)})
    }
    vals$cluster_colors[i]  <-  color.conv(v)
  })

  output$downloadPlot <- downloadHandler(
    filename = function() { paste('plot.',input$filetype, sep='') },
    content  = function(file) {
    if(input$filetype %in% c("png","bmp","jpeg","tiff")){
      ggsave(file,
             plotInput(),
             type    = "cairo",
             width   = input$width,
             height  = input$height,
             units   = input$un,
             dpi     = input$res)
    }
    else if(input$filetype %in% c("pdf")){
      ggsave(file,
             plotInput(),
             device  = cairo_pdf,
             width   = input$width,
             height  = input$height,
             units   = input$un,
             dpi     = input$res)
    }
    else if(input$filetype %in% c("eps","ps")){
      ggsave(file,
             plotInput(),
             device  = cairo_ps,
             width   = input$width,
             height  = input$height,
             units   = input$un,
             dpi     = input$res)

    }
    else if(input$filetype %in% c("svg","tex")){
      ggsave(file,
             plotInput(),
             device  = input$filetype,
             width   = input$width,
             height  = input$height,
             units   = input$un,
             dpi     = input$res)
    }
    })

  observeEvent(input$group_1, {
    if(input$select_polygon & is.null(input$plot1_brush)){
      vals$polygon_groupA  <-  vals$polygon
      vals$polygon_gAa     <-  TRUE
      ds$select_lasso('x', 'y', numpy$array(vals$polygon$i), numpy$array(vals$polygon$j))
      vals$GC_1            <-  (point.in.polygon(vals$data$i, vals$data$j, vals$polygon$i, vals$polygon$j)  !=  1)
      vals$Ncells_Group_A  <-  py_to_r(ds$count(selection=T))
      vals$polygon         <-  NULL
    } else {
      vals$polygon_gAa    <- FALSE
      vals$GC_1           <- rep(TRUE, nrow(vals$data))
      vals$GC_1.X         <- c(input$plot1_brush$xmin,input$plot1_brush$xmax)
      vals$GC_1.Y         <- c(input$plot1_brush$ymin,input$plot1_brush$ymax)
      res                 <- brushedPoints(vals$data, input$plot1_brush, allRows = TRUE)
      vals$GC_1           <- xor(vals$GC_1, res$selected_)
      vals$Ncells_Group_A <- sum(py_to_r(ds$count(binby=tuple(ds$x,ds$y),
                                                  shape=1L,
                                                  limits=tuple(tuple(vals$GC_1.X[1],vals$GC_1.X[2]),
                                                               tuple(vals$GC_1.Y[1],vals$GC_1.Y[2]))
                                )))
    }
    if(vals$gAa    == TRUE) {vals$group_Adist  <-  NULL}
    else{vals$gAa  <- TRUE}
  })
  observeEvent(input$group_2, {
    if(input$select_polygon & is.null(input$plot1_brush)){
      vals$polygon_groupB  <-  vals$polygon
      vals$polygon_gBa     <- TRUE
      ds$select_lasso('x', 'y', numpy$array(vals$polygon$i), numpy$array(vals$polygon$j))
      vals$GC_2            <- (point.in.polygon(vals$data$i, vals$data$j, vals$polygon$i, vals$polygon$j)!=1)
      vals$Ncells_Group_B  <- py_to_r(ds$count(selection=T))
      vals$polygon         <- NULL
    } else {
    vals$polygon_gBa    <- FALSE
    vals$GC_2           <- rep(TRUE, nrow(vals$data))
    vals$GC_2.X         <- c(input$plot1_brush$xmin,input$plot1_brush$xmax)
    vals$GC_2.Y         <- c(input$plot1_brush$ymin,input$plot1_brush$ymax)
    res                 <- brushedPoints(vals$data, input$plot1_brush, allRows = TRUE)
    vals$GC_2           <- xor(vals$GC_2, res$selected_)
    vals$Ncells_Group_B <- sum(py_to_r(ds$count(binby=tuple(ds$x,ds$y),
                                                shape=1L,
                                                limits=tuple(tuple(vals$GC_2.X[1],vals$GC_2.X[2]),tuple(vals$GC_2.Y[1],vals$GC_2.Y[2]))
                            )))
    }
    if(vals$gBa    == TRUE) {vals$group_Bdist<-NULL}
    else{vals$gBa  <- TRUE}
  })

  observeEvent(input$shuffle_colors, {
    vals$cluster_colors<-sample(vals$cluster_colors)
  })

  observeEvent(input$shuffle_cells, {
    vals$cluster_df<-vals$cluster_df[sample(1:nrow(vals$cluster_df)),]
  })

  observeEvent(input$pal.clusters, {
    if(input$show_clusters){
    if(input$pal.clusters       ==  "rainbow") {vals$cluster_colors   <-  sample(rainbow(vals$N.clusters))}
    else if(input$pal.clusters  ==  "matlab")  {vals$cluster_colors   <-  sample(matlab.like(vals$N.clusters))}
    else if(input$pal.clusters  ==  "matlab2") {vals$cluster_colors   <-  sample(matlab.like2(vals$N.clusters))}
      else{
        if(vals$N.clusters  <=  length(list.colors[[input$pal.clusters]])){vals$cluster_colors  <-  sample(list.colors[[input$pal.clusters]])}
        else{vals$cluster_colors  <-  sample(colorRampPalette(list.colors[[input$pal.clusters]])(vals$N.clusters))}
      }
    }
  })

  observeEvent(input$plot_expression, {
    if(input$plot_expression & input$choose_gene  ==  "") {
      showNotification("Select a Gene",type = "warning",duration = 10)
      updateCheckboxInput(session = session,inputId = "plot_expression",value = FALSE)
      return()
    }
    if(input$plot_expression){
    vals$plot.ID  =  2
    vals$rangeE   <- c(py_to_r(ds$min('Gene')),py_to_r(ds$max('Gene')))
    vals$trans    <- TRUE
    updateTextInput(session,     "clegend",          value = "Expression")
    updateCheckboxInput(session, "get_exp_clusters", value = FALSE)
    }
    else
    {
      vals$plot.ID  =  1
      updateTextInput(session, "clegend", value = "Density")
      updateTextInput(session, "title",   value="")
      vals$trans    <- TRUE
    }
  })

  observeEvent(input$plot_clusters, {
    if(input$plot_clusters){
      vals$N.clusters  =  py_to_r(ds$max('clusters'))
      withProgress(message = 'Getting group label positions', value = 0, {
      for(l in 1:vals$N.clusters){
        ds$select(paste0("clusters==",l,""))
        vals$cluster_labels <- rbind(vals$cluster_labels,
                                     data.frame(i  =  py_to_r(ds$mean("x",selection=TRUE)),
                                                j  =  py_to_r(ds$mean("y",selection=TRUE)),
                                                ID =  as.character(l),
                                                stringsAsFactors = F))
        incProgress(1/vals$N.clusters)
      }})
      vals$cluster_labels            <-  vals$cluster_labels[-1,]
      rownames(vals$cluster_labels)  <-  NULL
    } else {
    vals$cluster_labels=data.frame(i = 0.0, j = 0.0,ID = NA)
    }
  })
  
  observeEvent(input$show_legend_clusters, {
      updateTextInput(session, "showlegend", value = T)  
    })
  
  observeEvent(input$show_clusters, {
    if(input$show_clusters){
      if(!vals$load_clust & vals$load_xy){
        showNotification("Load Group Data",type = "warning",duration = 10)
        updateCheckboxInput(session,inputId = "show_clusters",value = FALSE) 
        return()
        }
      updateTextInput(session, "showlegend", value = FALSE)
      updateTextInput(session, "title",      value="")
      vals$plot.ID  =  3
      vals$trans    <- TRUE
    } else {
    updateTextInput(session, "showlegend", value = TRUE)
    updateCheckboxInput(session,inputId = "show_legend_clusters",value = FALSE)  
    updateCheckboxInput(session,inputId = "show_label_clusters",value = FALSE)  
    if(input$plot_expression){
       vals$plot.ID  =  2
       updateTextInput(session,"title",value=vals$loaded_feature)
     } else {
       vals$plot.ID  =  1
     }
     vals$trans      <- TRUE
    }
  })

  observeEvent(input$get_statistic, {
    if(vals$loaded_feature  ==  ""){showNotification("Select a feature first. E.g. Select a Gene or Compute a Gene Signature",type = "warning",duration = 10);return()}
    if(vals$gAa){
      if(vals$polygon_gAa){
        ds$select_lasso('x', 'y', numpy$array(vals$polygon_groupA$i), numpy$array(vals$polygon_groupA$j))
        vals$group_Adist  <-  py_to_r(ds$count(binby      =  tuple(ds$Gene),
                                               limits     =  tuple(vals$rangeE[1],vals$rangeE[2]),
                                               shape      =  tuple(40L),
                                               selection  =  TRUE))

        } else {
        ds$select(paste0("(x>",vals$GC_1.X[1],") & (x<",vals$GC_1.X[2],") & (y>",vals$GC_1.Y[1],") & (y<",vals$GC_1.Y[2],")"))
        vals$group_Adist  <-  py_to_r(ds$count(binby      =  tuple(ds$Gene),
                                               limits     =  tuple(vals$rangeE[1],vals$rangeE[2]),
                                               shape      =  tuple(40L),
                                               selection  =  TRUE))
        }
        Summary_Stat.A<-data.frame(Mean   = py_to_r(ds$mean('Gene',selection=TRUE)),
                                   Median = py_to_r(ds$median_approx('Gene',selection=TRUE)),
                                   SD     = py_to_r(ds$std('Gene',selection=TRUE)))
    }
    if(vals$gBa){
      if(vals$polygon_gBa){
        ds$select_lasso('x', 'y', numpy$array(vals$polygon_groupB$i), numpy$array(vals$polygon_groupB$j))
        vals$group_Bdist  <-  py_to_r(ds$count(binby      =  tuple(ds$Gene),
                                               limits     =  tuple(vals$rangeE[1],vals$rangeE[2]),
                                               shape      =  tuple(40L),
                                               selection  =  TRUE))
      } else {
        ds$select(paste0("(x>",vals$GC_2.X[1],") & (x<",vals$GC_2.X[2],") & (y>",vals$GC_2.Y[1],") & (y<",vals$GC_2.Y[2],")"))
        vals$group_Bdist  <-  py_to_r(ds$count(binby      =  tuple(ds$Gene),
                                               limits     =  tuple(vals$rangeE[1],vals$rangeE[2]),
                                               shape      =  tuple(40L),
                                               selection  =  TRUE))
      }
        Summary_Stat.B<-data.frame(Mean   = py_to_r(ds$mean('Gene',selection=TRUE)),
                                   Median = py_to_r(ds$median_approx('Gene',selection=TRUE)),
                                   SD     = py_to_r(ds$std('Gene',selection=TRUE)))
    }
    if(vals$gAa & vals$gBa){
      vals$Summary_Stat  <-  rbind(Summary_Stat.A,
                                   Summary_Stat.B)
      rownames(vals$Summary_Stat)  <-  c("Group A","Group B")
    }
    else if(vals$gAa & !vals$gBa){
      vals$Summary_Stat               <-  Summary_Stat.A
      rownames(vals$Summary_Stat)[1]  <-  "Group A"
    }
    else if(vals$gBa & !vals$gAa){
      vals$Summary_Stat               <-  Summary_Stat.B
      rownames(vals$Summary_Stat)[1]  <-  "Group B"
    }
  })

  observeEvent(input$get_exp_clusters, {
    if(input$get_exp_clusters  & !input$plot_clusters) {
      showNotification("Load group labels:  Metadata -> Groups -> Get Group Labels",type = "warning",duration = 10)
      updateCheckboxInput(session, "get_exp_clusters", value = FALSE)
      return()
    }
    if(input$get_exp_clusters  & !input$plot_expression) {
      showNotification("Select a Gene:  Analysis -> Cell Abundance & Gene Expression -> Select a Gene",type = "warning",duration = 10)
      updateCheckboxInput(session, "get_exp_clusters", value = FALSE)
      return()
    }
    if(input$get_exp_clusters){
    withProgress(message = 'Getting Gene Expression by Group', value = 0, {
      for(l in 1:vals$N.clusters){
        ds$select(paste0("clusters==",l,""))
        vals$cluster_stat  <-  rbind(vals$cluster_stat,
                                     data.frame(ID      =  vals$cluster_labels$ID[l],
                                                Ncells  =  py_to_r(ds$count(selection=TRUE)),
                                                Mean    =  py_to_r(ds$mean("Gene",selection=TRUE)),
                                                Median  =  py_to_r(ds$median_approx('Gene',selection=TRUE)),
                                                SD      =  py_to_r(ds$std('Gene',selection=TRUE)),
                                                Q1      =  py_to_r(ds$percentile_approx("Gene", 25,selection=TRUE)),
                                                Q3      =  py_to_r(ds$percentile_approx("Gene", 75,selection=TRUE)),
                                                stringsAsFactors = F))
          incProgress(1/vals$N.clusters)
      }
      vals$cluster_stat        <-  vals$cluster_stat[-1,]
      vals$cluster_stat[,3:7]  <-  round(vals$cluster_stat[,3:7],digits=5)
    })
    } else {
      vals$cluster_stat  =  data.frame(ID     =  0,
                                       Ncells =  0,
                                       Mean   =  0.0,
                                       Median =  0.0,
                                       SD     =  0.0,
                                       Q1     =  0.0,
                                       Q3     =  0.0,
                                       stringsAsFactors = F)
    }
  })

  output$table_exp.in.clusters <- DT::renderDataTable({DT::datatable(vals$cluster_stat,
                                                         options  =  list(paging = FALSE,scrollY = "500px",columnDefs =list(list(visible=FALSE, targets=c(5,6)))),
                                                         caption  =  vals$loaded_feature,
                                                         filter   =  'bottom',
                                                         rownames =  FALSE) %>%
                                                         formatCurrency('Ncells',
                                                                        currency  = "",
                                                                        digits    = 0,
                                                                        interval  = 3,
                                                                        mark      = ","
                                                                        ) %>%
                                                         formatRound(c("Mean","Median","SD","Q1","Q3"), 3)})

  output$plot_GEDistGroupA <- renderPlotly({
    if(input$plot_expression & !is.null(vals$group_Adist)){
      df  =  data.frame(Expression  =  seq(vals$rangeE[1], vals$rangeE[2], length.out  =  40),
                        Cell_Number =  vals$group_Adist)
      plot_ly(df,
              x      = ~Expression,
              y      = ~Cell_Number,
              type   = "bar",
              marker = list(color = "#FF7A00")) %>%
      layout(title = paste0(vals$loaded_feature,": Group A"),
             yaxis = list(title = "Cell Number"),
             xaxis = list(title = input$clegend))
    } else {
    plotly_empty()
    }
  })

  output$plot_GEDistGroupB <- renderPlotly({
    if(input$plot_expression & !is.null(vals$group_Bdist)){
      df  =  data.frame(Expression  =  seq(vals$rangeE[1], vals$rangeE[2], length.out=40),
                        Cell_Number =  vals$group_Bdist)
      plot_ly(df,
              x      = ~Expression,
              y      = ~Cell_Number,
              type   = "bar",
              marker = list(color = "#353839")) %>%
      layout(title = paste0(vals$loaded_feature,": Group B"),
             yaxis = list(title = "Cell Number"),
             xaxis = list(title = input$clegend))
    } else {
      plotly_empty()
    }
  })

  output$plot_DistGEinCluster_Census = renderPlotly({
      if(input$get_exp_clusters){
        cols=vals$cluster_colors
        colList=list()
        for(i in 1:length(cols)){colList[[i]]<-c((i-1)/(length(cols)-1),cols[i])}  
        plot_ly(vals$cluster_stat, 
                labels = ~ID, 
                values = ~Ncells,
                marker = list(colors = vals$cluster_colors,
                              line = list(color = '#FFFFFF', width = 1)),
                type = 'pie') %>%
          layout(title = 'Population by Group',
                 xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                 yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
      } else {
               plotly_empty()
      }
  })
  
  
  output$plot_DistGEinCluster_Density = renderPlotly({
    if(input$get_exp_clusters){
      s  =  input$table_exp.in.clusters_rows_selected
      if(length(s)){
        plotList  <-  list()
        distCl    <-  list()
        for(i in 1:length(s)){
          ds$select(paste0("clusters==",s[i],""))
          distCl[[i]]  =  py_to_r(ds$count(binby      =  tuple(ds$Gene),
                                           limits     =  tuple(vals$rangeE[1],vals$rangeE[2]),
                                           shape      =  tuple(40L),
                                           selection  =  TRUE))
        }
          for(i in 1:length(s)){
            df=data.frame(Expression   =  seq(vals$rangeE[1],vals$rangeE[2],length.out=40),
                          Cell_Number  =  distCl[[i]]/sum(distCl[[i]]))
            plotList[i]  <-  plot_ly(df,
                                     x      = ~Expression,
                                     y      = ~Cell_Number,
                                     type   = "bar",
                                     name   = paste0("Group: ",vals$cluster_labels$ID[s[i]]),
                                     marker = list(color = vals$cluster_colors[s[i]])) %>%
              layout(xaxis   = list(title = input$clegend),
                     yaxis   = list(title = "Fraction of Cells",
                                    range = c(0,input$y_axisR)),
                     title   = input$clegend)
          }
          subplot(plotList,
                  shareX  =  T,
                  shareY  =  T,
                  nrows   =  ceiling(length(s)/4)) %>%
            layout(title    =  vals$loaded_feature)
      } else {
        plotly_empty()
      }
    } else {
      plotly_empty()
    }
  })
  
  
  
  output$plot_DistGEinCluster_ECDF = renderPlotly({
    if(input$get_exp_clusters){
      s  =  input$table_exp.in.clusters_rows_selected
      if(length(s)){
        plotList  <-  list()
        distCl    <-  list()
        for(i in 1:length(s)){
          ds$select(paste0("clusters==",s[i],""))
          distCl[[i]]  =  py_to_r(ds$count(binby      =  tuple(ds$Gene),
                                           limits     =  tuple(vals$rangeE[1],vals$rangeE[2]),
                                           shape      =  tuple(40L),
                                           selection  =  TRUE))
        }
            p  <-  plot_ly()
          for(i in 1:length(s)){
            p  <-  add_lines(p,
                             x      =  seq(vals$rangeE[1],vals$rangeE[2],length.out=40),
                             y      =  cumsum(distCl[[i]])/sum(distCl[[i]]),
                             name   =  paste0("Group: ",vals$cluster_labels$ID[s[i]]),
                             marker =  list(color = vals$cluster_colors[s[i]]),
                             line   =  list(color = vals$cluster_colors[s[i]]),
                             mode   = 'lines')
          }
          p <- p %>% layout(yaxis  =  list(title = "Cumulative Probability"),
                            title  =  vals$loaded_feature,
                            xaxis  =  list(title = input$clegend))
          p
      } else {
        plotly_empty()
      }
    } else {
      plotly_empty()
    }
  })
 
  output$violin_plot = renderPlot(res = 72,{
    if(input$get_exp_clusters){
      s  =  input$table_exp.in.clusters_rows_selected
      if(length(s)){
        plotList  <-  list()
        distCl    <-  list()
        for(i in 1:length(s)){
          ds$select(paste0("clusters==",s[i],""))
          distCl[[i]]  =  py_to_r(ds$count(binby      =  tuple(ds$Gene),
                                           limits     =  tuple(vals$rangeE[1],vals$rangeE[2]),
                                           shape      =  tuple(40L),
                                           selection  =  TRUE))
        }
        ListPlots<-list()
        ListDFs<-list()
        maxX<-0
        maxY<-0
        rotate=FALSE
        for(i in 1:length(s)){
          ListDFs[[i]]<-data.frame(X=seq(vals$rangeE[1],vals$rangeE[2],length.out=40),Y=smooth.spline(distCl[[i]]/sum(distCl[[i]]))$y)
          for(j in 40:2){if(distCl[[i]][j]==0 & distCl[[i]][j-1]==0) {ListDFs[[i]]<-ListDFs[[i]][-j,]} else{break}}  
          if(maxX<max(ListDFs[[i]]$X)) maxX=max(ListDFs[[i]]$X)
          if(maxY<max(ListDFs[[i]]$Y)) maxY=max(ListDFs[[i]]$Y)
        }
        maxX<-maxX+maxX*0.25
        for(i in 1:length(s)){  
        if(i==1){
          ListPlots[[i]]<-ggplot(ListDFs[[i]],aes(x=X,y=Y))+
            stat_chull(fill=vals$cluster_colors[s[i]],colour = "black")+
            annotate("rect", xmin=-max(ListDFs[[i]]$Y)*0.05, xmax=max(ListDFs[[i]]$Y)*0.05, 
                     ymin=vals$cluster_stat$Q1[s[i]], ymax=vals$cluster_stat$Q3[s[i]], alpha=1, fill="grey10",color='black') +
            annotate("point", x = 0, y = vals$cluster_stat$Median[s[i]], colour = "white") +
            xlab(vals$cluster_labels$ID[s[i]])+
            ylab("Expression")+ 
            geom_hline(aes(yintercept=0))+
            theme(panel.background = element_rect(fill=NA),
                  #panel.grid.major.x = element_blank(),
                  #panel.grid.minor.x = element_blank(),
                  #axis.line.x = element_blank(),
                  axis.title.x=if(rotate){element_text(size=10,vjust = 0.5, hjust=1,angle=90)}else{element_text(size=14)},
                  axis.ticks.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title.y = element_text(size=14),
                  axis.line.y = element_line(),
                  axis.text.y = element_text(size=14))+
            ylim(0, maxX)
        } else{
          ListPlots[[i]]<-ggplot(ListDFs[[i]],aes(x=X,y=Y))+
            stat_chull(fill=vals$cluster_colors[s[i]],colour = "black")+
            annotate("rect", xmin=-max(ListDFs[[i]]$Y)*0.05, xmax=max(ListDFs[[i]]$Y)*0.05, 
                     ymin=vals$cluster_stat$Q1[s[i]], ymax=vals$cluster_stat$Q3[s[i]], alpha=1, fill="grey10",color='black') +
            annotate("point", x = 0, y = vals$cluster_stat$Median[s[i]], colour = "white") +
            xlab(vals$cluster_labels$ID[s[i]])+
            ylab("Expression")+  
            geom_hline(aes(yintercept=0))+ 
             theme(panel.background = element_rect(fill=NA),
                 #panel.grid.major.x = element_blank(),
                 #panel.grid.minor.x = element_blank(),
                 #axis.line.x = element_blank(),
                 text=element_text(size=14),
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank(),
                axis.title.y = element_blank(),
             axis.line.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank())+
            ylim(0, maxX) 
        }}
        p<-do.call(grid.arrange, c(ListPlots, list(nrow=1)))
        print(p)
      }
      }
  })
  
  
  output$table2 <- renderTable({
    if(input$plot_expression){
      vals$Summary_Stat
    }},
    include.rownames=TRUE)

  observeEvent(input$get_gene_sign_score, {
    GeneSignature  <-  unlist(strsplit(input$gene_sign, "[,;\n ]+"))
    GeneSignature  <-  vals$GeneNames[vals$GeneNames %in% GeneSignature]
    updateTextInput(session, "clegend", value = "Score")
    if(length(GeneSignature)  ==  0) {
      return()
    }
###CSV
    if(input$inputFileType_Genes  ==  "CSV"){
      if(input$n_workers  >  1){
        withProgress(message = 'Computing Gene Signature Scores',value = 0.5, {
        fun<-function(x,inputFilesGenes){
          scale(unlist(fread(paste0("zq ",
                                    gsub(" ","\\ ",inputFilesGenes,fixed=T),
                                    " ",
                                    x,
                                    " |cut -f2-"),
                             stringsAsFactors = F)))
        }
        inputFilesGenes   = vals$inputFiles$Genes
        VectSignScores  <-  rowMeans(sapply(future_lapply(1:length(GeneSignature), function(x) fun(GeneSignature[x],inputFilesGenes)),unlist))
        })
      } else {
        withProgress(message = 'Computing Gene Signature',value = 0.5, {
        VectSignScores  <-  scale(unlist(fread(paste0("zq ",
                                                      gsub(" ","\\ ",vals$inputFiles$Genes,fixed=T),
                                                      " ",
                                                      GeneSignature[1],
                                                      " |cut -f2-"),
                                               stringsAsFactors = F)))
        for(i in 2:length(GeneSignature)){
          VectSignScores  <-  VectSignScores  +  scale(unlist(fread(paste0("zq ",
                                                                           gsub(" ","\\ ",vals$inputFiles$Genes,fixed=T),
                                                                           " ",
                                                                           GeneSignature[i],
                                                                           " |cut -f2-"),
                                                                           stringsAsFactors = F)))
          incProgress(1/length(GeneSignature))
        }
        VectSignScores<-VectSignScores/length(GeneSignature)
        })
      }
      ds$add_column('Gene',data = numpy$array(drop(VectSignScores)))
      rm(VectSignScores)
    }
#HDF5
    if(input$inputFileType_Genes  ==  "HDF5"){
      if(input$n_workers  >  1){
        withProgress(message = 'Computing Gene Signature',value = 0.5, {
          fun<-function(inputFilesGenes,inputGenes,Gene_ind,Or){
            if(Or  ==  "R"){
            return(scale(drop(rhdf5::h5read(inputFilesGenes,
                                            inputGenes,
                                            index  =  list(Gene_ind,NULL)))))
            }
            else if(Or  ==  "C"){
            return(scale(drop(rhdf5::h5read(inputFilesGenes,
                                            inputGenes,
                                            index  =  list(NULL,Gene_ind)))))
            }
          }
          ind               = match(GeneSignature,vals$GeneNames)
          inputFilesGenes   = vals$inputFiles$Genes
          inputGenes        = input$inputGenes
          Or                = input$em_or
          VectSignScores  <-  rowMeans(sapply(future_lapply(1:length(GeneSignature), function(x) fun(inputFilesGenes,inputGenes,ind[x],Or)),unlist))
        })
      } else {
        withProgress(message = 'Computing Gene Signature',value = 0.5, {
          if(input$em_or  ==  "R"){
            VectSignScores  <-  scale(drop(rhdf5::h5read(vals$inputFiles$Genes,
                                                         input$inputGenes,
                                                         index  =  list(which(vals$GeneNames  ==  GeneSignature[1]),NULL))))
          }
          else if(input$em_or  ==  "C"){
            VectSignScores  <-  scale(drop(rhdf5::h5read(vals$inputFiles$Genes,
                                                         input$inputGenes,
                                                         index  =  list(NULL,which(vals$GeneNames  ==  GeneSignature[1])))))
          }
          for(i in 2:length(GeneSignature)){
            if(input$em_or  ==  "R"){
              VectSignScores  <-  VectSignScores  +  scale(drop(rhdf5::h5read(vals$inputFiles$Genes,
                                                                              input$inputGenes,
                                                                              index  =  list(which(vals$GeneNames  ==  GeneSignature[i]),NULL))))
            }
            else if(input$em_or  ==  "C"){
              VectSignScores  <-  VectSignScores  +  scale(drop(rhdf5::h5read(vals$inputFiles$Genes,
                                                                              input$inputGenes,
                                                                              index  =  list(NULL,which(vals$GeneNames  ==  GeneSignature[i])))))
            }
            incProgress(1/length(GeneSignature))
          }
          VectSignScores  <-  VectSignScores/length(GeneSignature)
        })
      }
      ds$add_column('Gene',data = numpy$array(drop(VectSignScores)))
      rm(VectSignScores)
    }
############End H5
    vals$plot.ID  =  2
    vals$trans    =  TRUE
    vals$rangeE   <- c(py_to_r(ds$min('Gene')),py_to_r(ds$max('Gene')))
    updateTextInput(session,"title",value=input$GeneSetName)
    vals$loaded_feature  <-  input$GeneSetName
    updateTextInput(session,     "choose_gene",      value  =  "")
    updateCheckboxInput(session, "get_exp_clusters", value  =  FALSE)
#    if(input$plot_expression  ==  FALSE) {
      updateTextInput(session,"plot_expression",value=TRUE)
 #   }
  })

  output$PlotPalette <- renderPlot({
    plot_palette(vals$ColPalette)
  })

  observeEvent(input$ColPallete, {
    if(!vals$InitColPalette) {  vals$InitColPalette  <-  TRUE  }
    else {  vals$ColPalette  <-  c(vals$ColPalette,input$ColPallete)  }
  })

  observeEvent(input$RemoveLastColor, {
    if(length(vals$ColPalette)  >  0) {
      vals$ColPalette  <-  vals$ColPalette[-length(vals$ColPalette)]
    }
  })

  observeEvent(input$AddNewPalette, {
    if(input$NewPaletteName  !=  "") {
      list.colors[[length(list.colors)+1]]     <<-  vals$ColPalette
      names(list.colors)[length(list.colors)]  <<-  input$NewPaletteName
    }
    updateSelectizeInput(session,
                         "pal.clusters",
                         choices=c("rainbow",
                                   "matlab",
                                   "matlab2",
                                   names(list.colors)))
    updateSelectizeInput(session,
                         "pal",
                         choices=c(dens_cols,input$NewPaletteName))
    updateSelectizeInput(session,
                         "pal_3d",
                         choices=c(dens_cols,input$NewPaletteName))
    dens_cols  <<-  c(dens_cols,input$NewPaletteName)
  })

  observeEvent(input$RemoveAllColors, {
    vals$InitColPalette  <-  FALSE
    vals$ColPalette      <-  NULL
  })

  observeEvent(input$file_fonts,{
    inFile <- input$file_fonts
    if (is.null(inFile)) return(NULL)
    if(Sys.info()[['sysname']]  ==  "Darwin"){
    file.copy(inFile$datapath,
              paste0("~/Library/Fonts/",inFile$name))
    extrafont::font_import("~/Library/Fonts/",prompt = FALSE)
    }
    else if(Sys.info()[['sysname']]  ==  "Linux"){
    system(paste0("sudo -S sh -c 'cp ",
                  inFile$datapath,
                  " /usr/share/fonts/",
                  inFile$name,
                  "; Rscript -e \"extrafont::font_import(prompt = F)\"'"),
                  input  =  readline("Enter your password: "))
    }
    extrafont::loadfonts()
    updateSelectInput(session,inputId = "fonts",           choices = fonts(),selected = "Arial")
    updateSelectInput(session,inputId = "annotation_fonts",choices = fonts(),selected = "Arial")
  })
  
  output$Plot3D <- renderPlotly({
    if(vals$trans_3d){
      vals$X_3d<-c(py_to_r(ds$min('x3d')),py_to_r(ds$max('x3d')));
      vals$Y_3d<-c(py_to_r(ds$min('y3d')),py_to_r(ds$max('y3d')));
      vals$Z_3d<-c(py_to_r(ds$min('z3d')),py_to_r(ds$max('z3d')))
      N=input$grid_plot3D
      withProgress(message = 'Processing on the grid',value = 0.5, {
      if(input$display_feature=="Density"){
      xycounts = ds$count(binby=tuple(ds$x3d,ds$y3d,ds$z3d), 
                          limits=tuple(tuple(vals$X_3d[1],vals$X_3d[2]),
                                       tuple(vals$Y_3d[1],vals$Y_3d[2]),
                                       tuple(vals$Z_3d[1],vals$Z_3d[2])),
                          shape=tuple(N,N,N),selection = FALSE)
      tmp<-py_to_r(xycounts)
      pl<-future_lapply(1:N, function(x) {sp<-summary(as(tmp[x,,],'dgTMatrix'));
      sp<-sp[complete.cases(sp),];
      sp<-cbind(h=rep(x,nrow(sp)),sp)})
      pl<-do.call("rbind",pl)
      vals$data_3d<-pl
      } else if (input$display_feature=="Gene Expression/Signature" & input$choose_gene != ""){
        xycounts = ds$count(binby=tuple(ds$x3d,ds$y3d,ds$z3d), 
                            limits=tuple(tuple(vals$X_3d[1],vals$X_3d[2]),
                                         tuple(vals$Y_3d[1],vals$Y_3d[2]),
                                         tuple(vals$Z_3d[1],vals$Z_3d[2])),
                            shape=tuple(N,N,N),selection = FALSE)
        tmp<-py_to_r(xycounts)
        pl<-future_lapply(1:N, function(x) {sp<-summary(as(tmp[x,,],'dgTMatrix'));
        sp<-sp[complete.cases(sp),];
        sp<-cbind(h=rep(x,nrow(sp)),sp)})
        allCells<-do.call("rbind",pl)
        if(input$select_method  ==  "Mean"){
          xycounts = ds$mean('Gene',
                             binby  =  tuple(ds$x3d,ds$y3d,ds$z3d),
                             limits =  tuple(tuple(vals$X_3d[1],vals$X_3d[2]),
                                             tuple(vals$Y_3d[1],vals$Y_3d[2]),
                                             tuple(vals$Z_3d[1],vals$Z_3d[2])),
                             shape  =  tuple(N,N,N),selection=FALSE)
        }
        else if(input$select_method=="Median"){
          xycounts = ds$median_approx('Gene',
                                      binby=tuple(ds$x3d,ds$y3d,ds$z3d),
                                      limits =  tuple(tuple(vals$X_3d[1],vals$X_3d[2]),
                                                      tuple(vals$Y_3d[1],vals$Y_3d[2]),
                                                      tuple(vals$Z_3d[1],vals$Z_3d[2])),
                                      shape  =  tuple(N,N,N),selection=FALSE)
        }
        else if(input$select_method=="Max"){
          xycounts = ds$max('Gene',
                            binby=tuple(ds$x3d,ds$y3d,ds$z3d),
                            limits =  tuple(tuple(vals$X_3d[1],vals$X_3d[2]),
                                            tuple(vals$Y_3d[1],vals$Y_3d[2]),
                                            tuple(vals$Z_3d[1],vals$Z_3d[2])),
                            shape  =  tuple(N,N,N),selection=FALSE)
        }
        tmp<-py_to_r(xycounts)
        pl<-future_lapply(1:N, function(x) {sp<-summary(as(tmp[x,,],'dgTMatrix'));
                                                     sp<-sp[complete.cases(sp),];
                                                     sp<-cbind(h=rep(x,nrow(sp)),sp)})
        pl<-do.call("rbind",pl)
        vals$data_3d<-pl
        
        if(input$select_method  ==  "Max") {
          vals$data_3d$x[!is.finite(vals$data_3d$x)]  <-  NA
        }
        
        vals$data_3d  <-  vals$data_3d[complete.cases(vals$data_3d),]
        zeros         <-  suppressMessages(dplyr::anti_join(allCells[,c(1,2,3)],vals$data_3d[,c(1,2,3)]))
        zeros$x       <-  rep(0,nrow(zeros))
        vals$data_3d  <-  rbind(zeros,vals$data_3d)
      }
      
      
      range_X       =  (vals$X_3d[2]-vals$X_3d[1])/input$grid_plot3D
      range_Y       =  (vals$Y_3d[2]-vals$Y_3d[1])/input$grid_plot3D
      range_Z       =  (vals$Z_3d[2]-vals$Z_3d[1])/input$grid_plot3D
      
      vals$data_3d$h   <- vals$X_3d[1]+range_X*(vals$data_3d$h-.5)
      vals$data_3d$i   <- vals$Y_3d[1]+range_Y*(vals$data_3d$i-.5)
      vals$data_3d$j   <- vals$Z_3d[1]+range_Z*(vals$data_3d$j-.5)
      
      rm(pl)
      tmp           =  update_slider_3d()
      vals$Cols3D   =  c(tmp$min,tmp$max)
      vals$trans_3d=FALSE
      
    })
    } 
    f1 <- list(
      family = "Arial, sans-serif",
      size = input$fs_3d,
      color = input$font_color_3d
    )
    f1t <- list(
      family = "Arial, sans-serif",
      size = input$tick_3d,
      color = input$font_color_3d
    )
    f2 <- list(
      family = "Arial, sans-serif",
      size = input$fs_3d,
      color = input$color_font_legend_3d
    )
    axx <- list(
      title = input$x_axis_label_3d,
      color = input$axis_color_3d,
      titlefont = f1,
      tickfont = f1t,
      showgrid = input$grid_3d
    )
    
    axy <- list(
      title = input$y_axis_label_3d,
      color = input$axis_color_3d,
      titlefont = f1,
      tickfont = f1t,
      showgrid = input$grid_3d
    )
    
    axz <- list(
      title = input$z_axis_label_3d,
      color = input$axis_color_3d,
      titlefont = f1,
      tickfont = f1t,
      showgrid = input$grid_3d
    )
    m <- list(
      l = 0,
      r = 0,
      b = 0,
      t = 30,
      pad = 0
    )
    if(vals$load_3d){
    cols=get_colors(input$pal_3d,
                    input$dichromat  !=  "No",
                    input$dichromat,
                    vals$n.cols)
    colList=list()
    for(i in 1:length(cols)){colList[[i]]<-c((i-1)/(length(cols)-1),cols[i])}  

    vals$p_3d<-plot_ly(vals$data_3d,
                       x      = ~h,
                       y      = ~i,
                       z      = ~j, 
                       marker = list(size = input$ps_3d, 
                                     color = if(input$scale_colorbar_3d){~log10(x)} else {~x}, 
                                     cmin=vals$Cols3D[1],
                                     cmax=vals$Cols3D[2],
                                     colorscale = colList,
                                     showscale = input$legend_3d,
                                     colorbar=list(title=input$colorbar_3d,
                                                   len=0.2,
                                                   tickcolor = input$color_font_legend_3d,
                                                   outlinecolor = input$color_font_legend_3d,
                                                   titlefont = f2, 
                                                   tickfont = f2,
                                                   ticks = 'outside'),
                                     symbol=if(input$square_points_3d) {'square'} else{'circle'},
                                     opacity=input$points_transp_3d), 
                       type = "scatter3d") %>%
      layout(title=input$title_3d,
             margin = m,
             paper_bgcolor =input$panel_background_3d,
             plot_bgcolor  =input$panel_background_3d,
             scene = list(xaxis=axx,
                          yaxis=axy,
                          zaxis=axz,
                          camera = list(eye = list(x = input$cv_x, 
                                                   y = input$cv_y, 
                                                   z = input$cv_z)
                          )
             ))
    }
    else{
      vals$p_3d<-plotly_empty()
    }
    vals$p_3d
    
  })
  observeEvent(input$change_colorbar_3d,{
    vals$Cols3D  <-  input$ColorbarRange3D
  })
  observeEvent(input$scale_colorbar_3d,{
    tmp       = update_slider_3d()
    vals$Cols3D <- c(tmp$min,tmp$max)
  })
  observeEvent(input$reset_view_3d, {
    vals$trans_3d<-TRUE
  })
  
  observeEvent(input$display_feature, {
    if(vals$load_3d){
      if(input$display_feature == "Density"){
        updateTextInput(session     = session,inputId = "colorbar_3d",    value ="Density")
        updateTextInput(session     = session,inputId = "title_3d",       value ="")
        } else if (input$display_feature=="Gene Expression/Signature" & input$choose_gene != ""){
        updateTextInput(session     = session,inputId = "colorbar_3d",    value ="Expression")
        updateTextInput(session     = session,inputId = "title_3d",       value =input$choose_gene)
        }
    vals$trans_3d<-TRUE  
    }
  })
  
  output$downloadPlot3D <- downloadHandler(
    filename <- function() { paste('plot.',input$filetype_3d, sep='') },
    content <- function(file) {
      b        <- vals$p_3d
      b$x$layout$width  <- input$width_3d 
      b$x$layout$height <- input$height_3d
      b        <- plotly_build(b)$x[c("data", "layout")]
      json     <- plotly:::to_JSON(b)
      filename_tmp<-tempfile(fileext = ".json")
      write(json, file=filename_tmp)
      OrcaPath<-system('which orca',intern = T)
      if(input$filetype_3d %in% c("png")){
        args <- c(
          "graph", filename_tmp,
          "-o", filename())
        invisible(processx::run(OrcaPath, args, echo = TRUE, spinner = TRUE,echo_cmd = T))
        file.copy(filename(), file)
      }
      if(input$filetype_3d %in% c("pdf")){
        args <- c(
          "graph", filename_tmp,
          "-o", filename(),
          "--format", "pdf")
        invisible(processx::run(OrcaPath, args, echo = TRUE, spinner = TRUE,echo_cmd = T))
        file.copy(filename(), file)
      }
      if(input$filetype_3d %in% c("svg")){
        args <- c(
          "graph", filename_tmp,
          "-o", filename(),
          "--format", "svg")
        invisible(processx::run(OrcaPath, args, echo = TRUE, spinner = TRUE,echo_cmd = T))
        file.copy(filename(), file)
      }
    }
  )
  getTerminalOutputCommand<-function(command){
    output$terminalOutput <- renderText({
      termOut<<-paste0(termOut,"<font size = '3' face = 'arial' color='white'><b>",":>",command,"</b></font>",sep = '<br/>')
      # Call custom javascript to scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      return(termOut)
    })
  }
  getTerminalOutputResult<-function(){
    output$terminalOutput <- renderText({
      session$sendCustomMessage(type = "scrollCallback", 1)
      return(termOut)
    })
  }
 # observeEvent(input$run_docker_command,{
#    isolate(command<-input$command)
#    getTerminalOutputCommand(command)
#  })
  observeEvent(input$clear_output,{
    updateTextInput(session = session,inputId = "command",value = "")
    command<-input$command
    termOut<<-NULL
    getTerminalOutputCommand(command)
  })
  observeEvent(input$gce_list_disks,{
    out<-gce_list_disks()$items[,1:6]
    tmp<-capture.output(out)
    termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>","List Disks:","</b></font>",sep = '<br/>')
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",tmp,"</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
  })
  observeEvent(input$gce_list_instances,{
    out<-gce_list_instances()
    tmp<-capture.output(out)
    termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>","List Instances:","</b></font>",sep = '<br/>')
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",tmp,"</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
  })
  observeEvent(input$gce_list_machinetypes,{
    out<-gce_list_machinetype()
    tmp<-capture.output(out)
    termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>","List Machine Types:","</b></font>",sep = '<br/>')
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",tmp,"</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
    updateSelectizeInput(session = session,inputId = "select_machine",choices = out$items$name)
  })
 
  observeEvent(input$run_vm,{
    vm  <<- gce_vm(template        = input$select_template,
                   name            = input$vm_name,
                   username        = input$GCE_username, 
                   password        = input$GCE_password,
                   dynamic_image   = input$container_name,
                   predefined_type = input$select_machine,
                   disk_size_gb    = input$disk_size)
    updateTextInput(session  = session,
                    inputId  = "ContainerToRun",
                    value    = input$container_name)
    termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>","Running VM: ",input$vm_name,"</b></font>",sep = '<br/>')
    #termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",tmp,"</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
  })
  
  observeEvent(input$get_running_vms,{
    list_all_instances<-gce_list_instances()
    running_machines<-which(list_all_instances$items$status=="RUNNING")
    if(length(running_machines)  ==  0) {
      showNotification("No running instances",type = "warning",duration = 10)
      updateSelectizeInput(session = session,inputId = "vm_name_sd",choices = NULL,selected = NULL,server = T)
      return()
    } else {
      updateSelectizeInput(session = session,inputId = "vm_name_sd",choices = list_all_instances$items$name[running_machines])
    }
  })
  observeEvent(input$ContainerToRun,{
    tag<<-input$ContainerToRun 
  })
  observeEvent(input$assign_vm,{
    if(is.null(input$vm_name_sd) | input$vm_name_sd=="") {
      showNotification("Select VM Name",type = "warning",duration = 10)
      return()
    } 
    vm   <<- gce_get_instance(input$vm_name_sd)
    list_containers <- gce_ssh(vm, "docker images|cut -f1 -d' '",capture_text = T)[-1] 
    updateSelectizeInput(session = session,
                         inputId = "ContainerToRun",
                         choices = list_containers,
                         selected = list_containers[1])
    tag <<- list_containers[1]
      
    termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>","Selected VM: ",input$vm_name,"</b></font>",sep = '<br/>')
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",capture.output(vm),"</font>",collapse = '<br/>'),sep = '<br/>')
    
    output$VM_assigned <- renderText(
      paste0("<font size = '2'>","Selected VM: ","<font size = '3' color=\"#FD0E35\"><b>", vm$name , "</b></font>")
      ) 
    vals$gcp_selected_files==""
    suppressWarnings(file.remove(paste0(tempdir(),"/hosts")));
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>", gce_ssh(vm, "uname -a",capture_text = T),"</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
  })
  
  observeEvent(input$stop_vm,{
    if(is.null(input$vm_name_sd ) | input$vm_name_sd=="") {
      showNotification("Select VM Name",type = "warning",duration = 10)
      return()
    }
    tmp<-capture.output(gce_vm_stop(input$vm_name_sd))
    termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>","Stopping VM: ",input$vm_name_sd,"</b></font>",sep = '<br/>')
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",tmp,"</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
    vm<<-NULL
    updateSelectizeInput(session = session,inputId = "vm_name_sd",choices = NULL,selected = NULL,server = T)
    output$VM_assigned <- renderText(
      paste0("<font size = '2'>","Selected VM: ","None", "</b></font>")
    ) 
    vals$gcp_selected_files==""
  })
  
  observeEvent(input$delete_vm,{
    if(is.null(input$vm_name_sd) | input$vm_name_sd=="") {
      showNotification("Select VM Name",type = "warning",duration = 10)
      return()
    }
    tmp<-capture.output(gce_vm_delete(input$vm_name_sd))
    termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>","Deleting VM: ",input$vm_name_sd,"</b></font>",sep = '<br/>')
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",tmp,"</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
    vm<<-NULL
    updateSelectizeInput(session = session,inputId = "vm_name_sd",choices = NULL,selected = NULL,server = T)
    output$VM_assigned <- renderText(
      paste0("<font size = '2'>","Selected VM: ","None", "</b></font>")
    ) 
    vals$gcp_selected_files==""
  })
   
  observeEvent(input$open_ssh,{
    if(is.null(vm)) {
      showNotification("Select and assign VM",type = "warning",duration = 10)
      return()
    }
    gce_ssh_browser(input$vm_name_sd)
  })
  
  observeEvent(input$show_processes,{
    if(is.null(vm)) {
      showNotification("Select and Assign VM",type = "warning",duration = 10)
      return()
    }
    tmp<-gce_ssh(vm, "top -n1 -b|head -n 20",capture_text = T)
    termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>","List processes: ",input$vm_name_sd,"</b></font>",sep = '<br/>')
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",tmp,"</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
  })
  observeEvent(input$run_ssh_command,{
    if(input$command=="" | is.null(vm)) {
      showNotification("Write a command",type = "warning",duration = 10)
      return()
    }
    isolate(command<-input$command)
    tmp<-gce_ssh(vm, input$command,capture_text = T)
    termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>",input$command,"</b></font>",sep = '<br/>')
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",tmp,"</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
    updateTextInput(session = session,inputId = "command",value = "")
  })
  output$VM_assigned <- renderText({ 
    if(is.null(vm)){
      paste0("<font size = '2'>","Selected VM: ",  "None","</b></font>")  
    } else {
      paste0("<font size = '2'>","Selected VM: ", "<font size = '3' color=\"#FD0E35\"><b>",vm$name,"</b></font>") 
    }
  })

  
  output$dirTable <- renderDataTable({
    vals$gcp_selected_dir<-gce_ssh(vm, paste0("ls -d -1 ",vals$gcp_dir,"*/"),capture_text = T)
    if(length(grep("cannot access",vals$gcp_selected_dir))!=0) {
      return(data.frame(".."))
    } else {
    return(data.frame(vals$gcp_selected_dir))
    }
  }, options  = list(paging   = F,
                     scrollY  = "200px"),rownames = NULL,colnames=NULL,caption =paste0("Folders: "),selection=list(mode="single"))
  
  output$fileTable <- renderDataTable({
    tmp<-gce_ssh(vm, paste0("ls -d -1 ",vals$gcp_dir,"*.*"),capture_text = T)
    if(length(grep("cannot access",tmp))!=0) tmp<-".."
    vals$gcp_files_indir<-tmp
    return(data.frame(tmp))
  }, options  = list(paging   = FALSE,
                     scrollY  = "200px"),rownames = NULL,colnames=NULL,caption =paste0("Files: "),selection=list(mode="multiple"))
  
  observeEvent(input$dirRIGHT, {
    s = input$dirTable_rows_selected
    if(length(s)>0){
      vals$gcp_dir<-vals$gcp_selected_dir[s]
    } else {
      return()
    }
  })
  
  observeEvent(input$dirUP, {
    dirs=unlist(strsplit(vals$gcp_dir,"/"))
    dirs=dirs[-length(dirs)]
    vals$gcp_dir<-paste0(paste0(dirs,collapse = "/"),"/")
  })
  
  
  observeEvent(input$copyToVM, {
    if(is.null(vm)) {
      showNotification("Select and Assign VM",type = "warning",duration = 10)
      return()
    }
    termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>","Copying Files to VM:","</b></font>",sep = '<br/>')
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",vals$inputFiles_GCP$FILES,"</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
    dest<-if(input$dest_gcp=="") {"~/"} else {input$dest_gcp} 
    gce_ssh_upload(vm,gsub(" ","\\ ",vals$inputFiles_GCP$FILES,fixed=T),dest,verbose = TRUE)
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>","Done","</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
  }
  )
  
  observeEvent(input$selectFilesFromVM, {
    s = input$fileTable_rows_selected
    if(length(s)>0){
     vals$gcp_selected_files<-vals$gcp_files_indir[s]
     termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>","Selectd Files on VM:","</b></font>",sep = '<br/>')
     termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",vals$gcp_selected_files,"</font>",collapse = '<br/>'),sep = '<br/>')
     getTerminalOutputResult()
    output$selectedScript <- renderText({ vals$gcp_selected_files[1] }) 
    } else {
      showNotification("Select Files",type = "warning",duration = 10)
      return()
    }
  })
  observeEvent(input$copyFromVM, {
    if(is.null(vm)) {
      showNotification("Select and Assign VM",type = "warning",duration = 10)
      return()
    }
    termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>","Copying Files from VM:","</b></font>",sep = '<br/>')
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",vals$gcp_selected_files,"</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
    dest<-if(input$dest_local=="") {"~/"} else {input$dest_local} 
    gce_ssh_download(vm,vals$gcp_selected_files,dest,verbose = TRUE)
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>","Done","</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
  }
  )
  observe({
    updateAceEditor(
      session, "ace", 
      theme = input$theme, 
      mode = input$mode,
      tabSize = input$size, 
      fontSize = input$FontSizeEditor,
      showInvisibles = as.logical(input$invisible)
    )
  })
  observeEvent(input$EditSelectedFile, {
    if(vals$gcp_selected_files=="") {
      showNotification("Select Files",type = "warning",duration = 10)
      return()
    }
    tmpdir<-tempdir()
    if(length(vals$gcp_selected_files)>0) vals$gcp_selected_files<-vals$gcp_selected_files[1] 
    gce_ssh_download(vm,vals$gcp_selected_files,tmpdir,verbose = TRUE,overwrite = T)
    vals$local_path_edited_file<-paste0(tmpdir,"/",basename(vals$gcp_selected_files))
    output$EditedFile <- renderText(
      paste0("<font size = '2' color=\"#FD0E35\"><b>",paste0(vals$gcp_selected_files," (VM: ",vm$name,")"), "</b></font>")
    ) 
    #print(vals$local_path_edited_file)
    updateAceEditor(session, "ace", value = paste0(readLines(vals$local_path_edited_file),collapse = "\r"))
  })
  
  observeEvent(input$CreateNew, {
    output$EditedFile <- renderText(
      paste0("")
    ) 
    vals$local_path_edited_file<-tempfile()
    updateAceEditor(session, "ace", value = "")
  })
  
  observeEvent(input$SaveFile, {
    file_name<-vals$local_path_edited_file
    writeLines(text = input$ace,file_name)
    system(paste0("awk -v RS='\r' '{print $0}' ",file_name," > ",file_name,"tmp ; mv ",file_name,"tmp ",file_name))
    gce_ssh_upload(vm,file_name,dirname(vals$gcp_selected_files),verbose = TRUE)
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>","File Saved:","</font>",collapse = '<br/>'),sep = '<br/>')
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",paste0(vals$gcp_selected_files," (VM: ",vm$name,")"),"</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
  })
  
  observeEvent(input$SaveAs, {
    if(vals$local_path_edited_file==""){
    file_name<-paste0(tempdir(),"/",basename(input$SaveAsFile))
    } else {
    file_name<-paste0(dirname(vals$local_path_edited_file),"/",basename(input$SaveAsFile))
    }
    writeLines(text = input$ace,file_name)
    system(paste0("awk -v RS='\r' '{print $0}' ",file_name," > ",file_name,"tmp ; mv ",file_name,"tmp ",file_name))
    dest<-dirname(input$SaveAsFile)
    gce_ssh_upload(vm,file_name,dest,verbose = TRUE)
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>","File Saved:","</font>",collapse = '<br/>'),sep = '<br/>')
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",paste0(input$SaveAsFile," (VM: ",vm$name,")"),"</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
  })

  observeEvent(input$runScriptOnVM, {
    script<-basename(vals$gcp_selected_files[1])
    dir<-dirname(vals$gcp_selected_files[1])
    cmd=paste0("run -v ",dir,":/home --rm ",tag," /bin/bash /home/",script,"  &>logs.txt &")
    print(cmd)
    docker_cmd(vm,cmd = cmd, capture_text = TRUE)
  })
  
  observeEvent(input$gce_initialize, {
    Sys.setenv("GCE_AUTH_FILE" = input$GCE_AUTH_FILE,
               "GCE_DEFAULT_PROJECT_ID" = input$GCE_DEFAULT_PROJECT_ID,
               "GCE_DEFAULT_ZONE" = input$GCE_DEFAULT_ZONE)
    print(input$GCE_AUTH_FILE)
    print(input$GCE_DEFAULT_PROJECT_ID)
    print(input$GCE_DEFAULT_ZONE)
    require(googleComputeEngineR)
  })
  output$plot3dto2d <- renderPlot({
    tmp<-scatter3D( vals$data_3d$h, 
                     vals$data_3d$i, 
                     vals$data_3d$j,
                     phi = input$phi,
                     theta = input$theta)
    tmp<-tmp %*% rbind(vals$data_3d$h,
                       vals$data_3d$i,
                       vals$data_3d$j,
                      rep(1,length(vals$data_3d$h)))
    tmp<-t(tmp)
    colnames(tmp)<-c("x","y","z","w")
    ggplot(as.data.frame(tmp), aes(x=x, y=y, colour=if(input$scale_colorbar_3d){log10(vals$data_3d$x)} else {vals$data_3d$x})) +
      geom_point(size=0.1) +
      theme(plot.background    =  element_rect(fill  =  input$plot_background),
            legend.background  =  element_rect(fill  =  input$plot_background),
            panel.background   =  element_rect(fill  =  input$panel_background, colour = "black"),
            panel.grid.major   =  if(input$showgrid) {element_line(colour = "grey")} else {element_blank()},
            panel.grid.minor   =   element_blank(),
            legend.position    =  "none",
            axis.line          =  element_line(colour = "black",size = 0.5),
            text               =  element_text(family = input$fonts, face = "bold",size = input$fs),
            axis.ticks         =  element_line(size   = 0.5,colour = 'black'),
            axis.text          =  element_text(size   = input$fs,colour = 'black'))+
      scale_colour_gradientn(colours  =  get_colors(input$pal_3d,
                                                    input$dichromat  !=  "No",
                                                    input$dichromat,
                                                    vals$n.cols))
    
  })
  
  observeEvent(input$DMap_run, {
    Parameters.DMap<<-data.frame( path=paste0("/home/",basename(input$DMap_Path)),
                                  nNN=input$DMap_nNN,
                                  k=input$DMap_k,
                                  EigDecompMethod=input$DMap_EigDecompMethod,
                                  nLocalsigma=input$DMap_ls,
                                  nThreads=input$DMap_nThreads,
                                  AnnMethod=input$DMap_NNMethod,
                                  nTrees=input$DMap_Annoy_nTrees,
                                  M=input$DMap_Nmslib_M,
                                  efC=input$DMap_Nmslib_efC,
                                  efS =input$DMap_Nmslib_efS
                                  )
    DMap_filename=paste0(tempdir(),"/DiffusionMap_InputParameters.json")
    jsonlite::write_json(Parameters.DMap,DMap_filename,pretty=T)
    if(is.null(vm)) {
      showNotification("Select and Assign VM",type = "warning",duration = 10)
      return()
    }
    termOut<<-paste0(termOut,
                     "<font size = '2' face = 'arial' color='white'><b>",":>",
                     "Copying Files to VM:",
                     "</b></font>",
                     sep = '<br/>')
    termOut<<-paste0(termOut,
                     paste0("<font size = '2' face = 'arial' color='white'>",
                            "DiffusionMap_InputParameters.json",
                            "</font>",
                            collapse = '<br/>'),
                     sep = '<br/>')
    getTerminalOutputResult()
    dest="~/"
    gce_ssh_upload(vm,DMap_filename,dest,verbose = TRUE)
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",
                                    "Done",
                                    "</font>",
                                    collapse = '<br/>'),
                     sep = '<br/>')
    getTerminalOutputResult()
    cmd=paste0("run -v ",dirname(input$DMap_Path),":/home --rm ",tag," R -e 'source(\"/scsvatools/DMap.R\")' &>DMap.logs &")
    print(cmd)
    docker_cmd(vm,cmd = cmd, capture_text = TRUE)
    
  })
  
  observeEvent(input$NNG_run, {
    Parameters.NNG<<-data.frame( path= if(input$NNG_Path==""){input$NNG_Path}else{paste0("/home/",basename(input$NNG_Path))},
                                 nNN=input$NNG_nNN,
                                 nThreads=input$NNG_nThreads,
                                 AnnMethod=input$NNG_NNMethod,
                                 nTrees=input$NNG_Annoy_nTrees,
                                 M=input$NNG_Nmslib_M,
                                 efC=input$NNG_Nmslib_efC,
                                 efS =input$NNG_Nmslib_efS 
    )
    NNG_filename=paste0(tempdir(),"/NNG_InputParameters.json")
    jsonlite::write_json(Parameters.NNG,NNG_filename,pretty=T)
    if(is.null(vm)) {
      showNotification("Select and Assign VM",type = "warning",duration = 10)
      return()
    }
    termOut<<-paste0(termOut,
                     "<font size = '2' face = 'arial' color='white'><b>",":>",
                     "Copying Files to VM:",
                     "</b></font>",
                     sep = '<br/>')
    termOut<<-paste0(termOut,
                     paste0("<font size = '2' face = 'arial' color='white'>",
                            "DiffusionMap_InputParameters.json",
                            "</font>",
                            collapse = '<br/>'),
                     sep = '<br/>')
    getTerminalOutputResult()
    dest="~/"
    gce_ssh_upload(vm,NNG_filename,dest,verbose = TRUE)
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",
                                    "Done",
                                    "</font>",
                                    collapse = '<br/>'),
                     sep = '<br/>')
    getTerminalOutputResult()
    cmd=paste0("run -v ",if(input$NNG_Path==""){"~/"}else{dirname(input$NNG_Path)},":/home --rm ",tag," R -e 'source(\"/scsvatools/NNG.R\")' &>NNG.logs &")
    print(cmd)
    docker_cmd(vm,cmd = cmd, capture_text = TRUE)
    
  })
  
  observeEvent(input$FLE_run, {
  Parameters.FA_3D<<-data.frame( path= if(input$FLE_Path==""){input$FLE_Path}else{paste0("/home/",basename(input$FLE_Path))},
                                 memmory = paste0(input$FLE_memmory,"g"),
                                 nsteps=input$FLE_nsteps,
                                 nThreads=input$FLE_nthreads,
                                 scalingRatio=input$FLE_scalingRatio,
                                 seed=1,
                                 barnesHutTheta=input$FLE_barnesHutTheta,
                                 barnesHutUpdateIter=input$FLE_barnesHutUpdateIter,
                                 updateCenter=input$FLE_updateCenter,
                                 barnesHutSplits=input$FLE_barnesHutSplits,
                                 restart=input$FLE_restart)
    FLE_filename=paste0(tempdir(),"/ForceDirectedLayout_InputParameters.json")
    jsonlite::write_json(Parameters.FA_3D,FLE_filename,pretty=T)
    if(is.null(vm)) {
      showNotification("Select and Assign VM",type = "warning",duration = 10)
      return()
    }
    termOut<<-paste0(termOut,
                     "<font size = '2' face = 'arial' color='white'><b>",":>",
                     "Copying Files to VM:",
                     "</b></font>",
                     sep = '<br/>')
    termOut<<-paste0(termOut,
                     paste0("<font size = '2' face = 'arial' color='white'>",
                            "ForceDirectedLayout_InputParameters.json",
                            "</font>",
                            collapse = '<br/>'),
                     sep = '<br/>')
    getTerminalOutputResult()
    dest="~/"
    gce_ssh_upload(vm,FLE_filename,dest,verbose = TRUE)
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>",
                                    "Done",
                                    "</font>",
                                    collapse = '<br/>'),
                     sep = '<br/>')
    getTerminalOutputResult()
    cmd=paste0("run -v ",if(input$FLE_Path==""){"~/"}else{dirname(input$FLE_Path)},":/home --rm ",tag," R -e 'source(\"/scsvatools/FLE.R\")' &> FLE.logs &")
    print(cmd)
    docker_cmd(vm,cmd = cmd, capture_text = TRUE)
    
  })
  
  observeEvent(input$FLE_getDistances, {
    if(is.null(vm)) {
      showNotification("Select and Assign VM",type = "warning",duration = 10)
      return()
    }
    vals$FLE_dist<-FALSE
    termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>","Copying Files from VM:","</b></font>",sep = '<br/>')
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'> FLE_InputData_FLE.distances.txt</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
    dest=tempdir()
    gce_ssh_download(vm,"~/NNG_output_FLE.distances.txt",dest,verbose = TRUE)
    termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'>","Done","</font>",collapse = '<br/>'),sep = '<br/>')
    getTerminalOutputResult()
    FLE_dist<<-as.data.frame(fread(paste0(dest,"/NNG_output_FLE.distances.txt")))
    vals$FLE_dist<-TRUE
  })
  
  output$plot_FLE_dist = renderPlotly({
    if(vals$FLE_dist){
    p  =  plot_ly()
    p  =  plot_ly(FLE_dist,
                   x      =  ~step,
                   y      =  ~log10(distance),
                   type = 'scatter', 
                   mode = 'lines',
                   line = list(color = '#004C71', width = 0.5)) %>% 
      add_lines(y = ~fitted(loess(log10(distance) ~ step, data = FLE_dist, span=0.2)),
                line = list(color = 'black',width=1.5)) %>%
      layout(title = "",
             yaxis = list(title = "Total distance, log10"),
             xaxis = list (title = "Iteration"),
             showlegend = FALSE)
    p
    }
    else {
      plotly_empty()
    }
  })
 
observeEvent(input$FLE_checkProgress, {
    
    val= as.integer(gce_ssh(vm, "cat ~/FLE.logs |grep -A4 -B3 'Edges loaded' |fgrep -o '*'  | wc -l",capture_text = T))
    output$plot_FLE_checkProgress<- renderPlotly({
      plot_ly(x = val, type = 'bar', orientation = 'h',marker=list(color='#004C71')) %>% 
        layout(xaxis = list(range = c(0, 100),showline = TRUE,mirror = "ticks",
                            ticksuffix = "%"),
               yaxis=list(
                 title = "",
                 zeroline = FALSE,
                 showline = TRUE,
                 mirror = "ticks",
                 showticklabels = FALSE,
                 showgrid = FALSE
               ),
               margin = list(r = 50)) %>% config(displayModeBar = F)
    })
  })

observeEvent(input$FLE_loadCoordinates, {
  termOut<<-paste0(termOut,"<font size = '2' face = 'arial' color='white'><b>",":>","Copying Files from VM:","</b></font>",sep = '<br/>')
  termOut<<-paste0(termOut,paste0("<font size = '2' face = 'arial' color='white'> NNG_output_FLE.txt</font>",collapse = '<br/>'),sep = '<br/>')
  getTerminalOutputResult()
  gce_ssh_download(vm,"~/NNG_output_FLE.txt","~/",verbose = TRUE)
  file_name="~/NNG_output_FLE.txt"
  tmp=fread(file_name)
  tmp=tmp[order(tmp$id),]
  if(ds$length()!=nrow(tmp)){
          if(input$show_clusters | input$plot_clusters | input$plot_expression) {
            updateCheckboxInput(session = session,inputId = "show_clusters",value = FALSE)
            updateCheckboxInput(session = session,inputId = "plot_clusters",value = FALSE)
            updateCheckboxInput(session = session,inputId = "plot_expression",value = FALSE)
          }
          ds <<- vaex$from_arrays(x=numpy$array(tmp$x))
          ds$add_virtual_column("x3d", "x")
        } else {
          ds$add_column('x3d',data = numpy$array(tmp$x)) 
        }
  ds$add_column('y',data = numpy$array(tmp$y))
  ds$add_virtual_column("y3d", "y")
  ds$add_column('z3d',data = numpy$array(tmp$z))
  gc()
  gc_python$collect()
  vals$X_3d    <- c(py_to_r(ds$min('x3d')),py_to_r(ds$max('x3d')))
  vals$Y_3d    <- c(py_to_r(ds$min('y3d')),py_to_r(ds$max('y3d')))
  vals$Z_3d    <- c(py_to_r(ds$min('z3d')),py_to_r(ds$max('z3d')))
  
  vals$load_3d        =  TRUE
  vals$trans_3d       =  TRUE
  vals$plot.ID        =  1
  vals$trans           <- TRUE
  vals$initialize      <- TRUE
  vals$X               <- c(py_to_r(ds$min('x')),py_to_r(ds$max('x')))
  vals$Y               <- c(py_to_r(ds$min('y')),py_to_r(ds$max('y')))
  vals$group_Adist     =  NULL
  vals$group_Bdist     =  NULL
  vals$Ncells_Group_A  =  0
  vals$Ncells_Group_B  =  0
  vals$gAa             =  FALSE
  vals$gBa             =  FALSE
  vals$polygon         <- NULL
  Summary_Stat         =  data.frame(Mean = 0.0, Median = 0.0,SD = 0.0,stringsAsFactors = F)
    
})


observeEvent(input$FS_ls_1, {
  if(input$FS_type_1=="local"){
    vals$path_FS_1<-"/"
    all.files=list.files(vals$path_FS_1,rec=F,full.names = T)
    vals$FS_ls_1<-data.frame(Content=substring(c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),2),stringsAsFactors = F)
  } else if(input$FS_type_1=="gs"){
    vals$path_FS_1<-"gs://"
    vals$FS_ls_1<-data.frame(Content=system(paste0(GCU,' ls ',vals$path_FS_1),intern = T),stringsAsFactors = F)
  }
})

output$FS_1 <- renderDataTable({
      DT::datatable(vals$FS_ls_1,
                    options  = list(paging   = FALSE,
                                    scrollY  = "400px"),rownames = NULL,colnames=NULL,caption =vals$path_FS_1,selection=list(mode="multiple"))
    })



observeEvent(input$FS_1_dirRIGHT, {
  s = input$FS_1_rows_selected
  if(length(s)>0){
    if(input$FS_type_1=="local"){
    print(vals$FS_ls_1$Content[s[1]])  
    vals$path_FS_1<-vals$FS_ls_1$Content[s[1]]
    all.files=list.files(vals$path_FS_1,rec=F,full.names = T)
    vals$FS_ls_1<-data.frame(Content=c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),stringsAsFactors = F)
    } else if(input$FS_type_1=="gs"){
      vals$path_FS_1<-vals$FS_ls_1$Content[s[1]]
      vals$FS_ls_1<-data.frame(Content=setdiff(system(paste0(GCU,' ls ',vals$path_FS_1),intern = T),vals$path_FS_1),stringsAsFactors = F)
    }
  } else {
    return()
  }
})

observeEvent(input$FS_1_dirUP, {
  if(input$FS_type_1=="local"){
  vals$path_FS_1<-dirname(vals$path_FS_1)
  all.files=list.files(vals$path_FS_1,rec=F,full.names = T)
  if(vals$path_FS_1=="/"){
    vals$FS_ls_1<-data.frame(Content=substring(c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),2),stringsAsFactors = F)
  } else {
  vals$FS_ls_1<-data.frame(Content=c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),stringsAsFactors = F)}
  }
  else if(input$FS_type_1=="gs"){
    vals$path_FS_1<-dirname(vals$path_FS_1)
    print(vals$path_FS_1)
    vals$FS_ls_1<-data.frame(Content=setdiff(system(paste0(GCU,' ls ',if(vals$path_FS_1=="gs:"){"gs://"}else{vals$path_FS_1}),intern = T),vals$path_FS_1),stringsAsFactors = F)
  }
  
})

observeEvent(input$FS_mkdir_1, {
  dir.create(path = paste0(vals$path_FS_1,"/",input$FS_mkdir_dirname_1),showWarnings = T)
  all.files=list.files(vals$path_FS_1,rec=F,full.names = T)
  vals$FS_ls_1<-data.frame(Content=c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),stringsAsFactors = F)
}
)

observeEvent(input$FS_ls_2, {
  if(input$FS_type_2=="local"){
    vals$path_FS_2<-"/"
    all.files=list.files(vals$path_FS_2,rec=F,full.names = T)
    vals$FS_ls_2<-data.frame(Content=substring(c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),2),stringsAsFactors = F)
  } else if(input$FS_type_2=="gs"){
    vals$path_FS_2<-"gs://"
    vals$FS_ls_2<-data.frame(Content=system(paste0(GCU,' ls ',vals$path_FS_2),intern = T),stringsAsFactors = F)
  }
})

output$FS_2 <- renderDataTable({
  DT::datatable(vals$FS_ls_2,
                options  = list(paging   = FALSE,
                                scrollY  = "400px"),
                rownames = NULL,
                colnames=NULL,
                caption =vals$path_FS_2,
                selection=list(mode="multiple"))
})
#Copy
observeEvent(input$FS_Right, {
  s = input$FS_1_rows_selected
  if(length(s)>0){
   if(input$FS_type_1=="local" & input$FS_type_2=="local"){
    for(i in 1:length(s)){
        file.copy(from = as.character(vals$FS_ls_1$Content[s[i]]),
                to = vals$path_FS_2,
                recursive = T,
                overwrite = T
                 )
    }
     all.files=list.files(vals$path_FS_2,rec=F,full.names = T)
     vals$FS_ls_2<-data.frame(Content=c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),stringsAsFactors = F)
   } else {
    system(paste0(GCU,' -m cp -r ',paste(vals$FS_ls_1$Content[s],collapse = " "),' ',vals$path_FS_2),intern = T) 
    vals$FS_ls_2<-data.frame(Content=system(paste0(GCU,' ls ',vals$path_FS_2),intern = T),stringsAsFactors = F) 
   }
  } else {
    return()
  }
})
#Copy
observeEvent(input$FS_Left, {
  s = input$FS_2_rows_selected
  if(length(s)>0){
    if(input$FS_type_1=="local" & input$FS_type_2=="local"){
      for(i in 1:length(s)){
        file.copy(from = as.character(vals$FS_ls_2$Content[s[i]]),
                  to = vals$path_FS_1,
                  recursive = T,
                  overwrite = T
        )
      }
      all.files=list.files(vals$path_FS_1,rec=F,full.names = T)
      vals$FS_ls_1<-data.frame(Content=c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),stringsAsFactors = F)
    } else {
      system(paste0(GCU,' -m cp -r ',paste(vals$FS_ls_2$Content[s],collapse = " "),' ',vals$path_FS_1),intern = T) 
      vals$FS_ls_1<-data.frame(Content=system(paste0(GCU,' ls ',vals$path_FS_1),intern = T),stringsAsFactors = F)
    }
  } else {
    return()
  }
})

observeEvent(input$FS_2_dirRIGHT, {
  s = input$FS_2_rows_selected
  if(length(s)>0){
    if(input$FS_type_2=="local"){
      print(vals$FS_ls_2$Content[s[1]])  
      vals$path_FS_2<-vals$FS_ls_2$Content[s[1]]
      all.files=list.files(vals$path_FS_2,rec=F,full.names = T)
      vals$FS_ls_2<-data.frame(Content=c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),stringsAsFactors = F)
    } else if(input$FS_type_2=="gs"){
      vals$path_FS_2<-vals$FS_ls_2$Content[s[1]]
      vals$FS_ls_2<-data.frame(Content=setdiff(system(paste0(GCU,' ls ',vals$path_FS_2),intern = T),vals$path_FS_2),stringsAsFactors = F)
    }
  } else {
    return()
  }
})

observeEvent(input$FS_2_dirUP, {
  if(input$FS_type_2=="local"){
    vals$path_FS_2<-dirname(vals$path_FS_2)
    all.files=list.files(vals$path_FS_2,rec=F,full.names = T)
    if(vals$path_FS_2=="/"){
      vals$FS_ls_2<-data.frame(Content=substring(c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),2),stringsAsFactors = F)
    } else {
      vals$FS_ls_2<-data.frame(Content=c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),stringsAsFactors = F)}
  }
  else if(input$FS_type_2=="gs"){
    vals$path_FS_2<-dirname(vals$path_FS_2)
    vals$FS_ls_2<-data.frame(Content=setdiff(system(paste0(GCU,' ls ',if(vals$path_FS_2=="gs:"){"gs://"}else{vals$path_FS_2}),intern = T),vals$path_FS_1),stringsAsFactors = F)
  }
})

observeEvent(input$FS_mkdir_2, {
  dir.create(path = paste0(vals$path_FS_2,"/",input$FS_mkdir_dirname_2),showWarnings = T)
  all.files=list.files(vals$path_FS_2,rec=F,full.names = T)
  vals$FS_ls_2<-data.frame(Content=c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),stringsAsFactors = F)
}
)

observeEvent(input$FS_mkdir_bucket_1, {
  system(paste0(GCU,' mb ',input$FS_mkdir_bucket_name_1),intern = T)
  vals$path_FS_1<-"gs://"
  vals$FS_ls_1<-data.frame(Content=system(paste0(GCU,' ls ',vals$path_FS_1),intern = T),stringsAsFactors = F)
  }
)

observeEvent(input$FS_mkdir_bucket_2, {
  system(paste0(GCU,' mb ',input$FS_mkdir_bucket_name_2),intern = T)
  vals$path_FS_2<-"gs://"
  vals$FS_ls_2<-data.frame(Content=system(paste0(GCU,' ls ',vals$path_FS_2),intern = T),stringsAsFactors = F)
}
)

observeEvent(input$FS_remove_local_1, {
    s = input$FS_1_rows_selected
    if(length(s)>0){
        for(i in 1:length(s)){
          unlink(vals$FS_ls_1$Content[s[i]],recursive = T)
        }
        all.files=list.files(vals$path_FS_1,rec=F,full.names = T)
        vals$FS_ls_1<-data.frame(Content=c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),stringsAsFactors = F)
    } else {
      return()
    }
  })

observeEvent(input$FS_remove_local_2, {
  s = input$FS_2_rows_selected
  if(length(s)>0){
    for(i in 1:length(s)){
      unlink(vals$FS_ls_2$Content[s[i]],recursive = T)
    }
    all.files=list.files(vals$path_FS_2,rec=F,full.names = T)
    vals$FS_ls_2<-data.frame(Content=c(all.files[file.info(all.files)$isdir],all.files[!file.info(all.files)$isdir]),stringsAsFactors = F)
  } else {
    return()
  }
})

observeEvent(input$FS_remove_object_1, {
  s = input$FS_1_rows_selected
  if(length(s)>0){
    system(paste0(GCU,' -m rm -r ',paste(vals$FS_ls_1$Content[s],collapse = " ")),intern = T) 
    vals$FS_ls_1<-data.frame(Content=system(paste0(GCU,' ls ',vals$path_FS_1),intern = T),stringsAsFactors = F)
  } else {
    return()
  }
})

observeEvent(input$FS_remove_object_2, {
  s = input$FS_2_rows_selected
  if(length(s)>0){
    system(paste0(GCU,' -m rm -r ',paste(vals$FS_ls_2$Content[s],collapse = " ")),intern = T) 
    vals$FS_ls_2<-data.frame(Content=system(paste0(GCU,' ls ',vals$path_FS_2),intern = T),stringsAsFactors = F)
  } else {
    return()
  }
})

observeEvent(input$show_multiplot,{
  if(length(vals$list_plots)==0 & input$show_multiplot==TRUE){ 
    updateCheckboxInput(session = session,inputId = "show_multiplot",value = FALSE)
    isolate(showNotification("Nothing to plot",type = "warning",duration = 10))
  }
})

observeEvent(input$collect_plots,{
  file=paste0(tempfile(tmpdir = tmpdir),".png")
  ggsave(file,
         plotInput(),
         type    = "cairo",
         width   = input$width,
         height  = input$height,
         units   = input$un,
         dpi     = input$res)
  vals$list_plots[[length(vals$list_plots)+1]]<<-file
})


output$multiplot <- renderImage({
  input$collect_plots
  ImageVec<-magick::image_read(vals$list_plots[[1]])
  if(length(vals$list_plots)>1){
    for(i in 2:length(vals$list_plots)){
      ImageVec<-c(ImageVec,magick::image_read(vals$list_plots[[i]]))
    }
  }

  d=1:length(ImageVec);n=input$multiplot_ncols;chunks<-split(d, ceiling(seq_along(d)/n))
  images_byRow<-magick::image_append(ImageVec[chunks[[1]]],stack = F)
  if(length(chunks)>1){
    for(i in 2:length(chunks)){
      images_byRow<-c(images_byRow,magick::image_append(ImageVec[chunks[[i]]],stack = F))
    }
  }
  images<-magick::image_append(images_byRow,stack = T)
  tmpfile <- images  %>% 
    magick::image_write(tempfile(fileext='png'), format = 'png')
  tmp=magick::image_info(images)
  list(src = tmpfile,
       contentType = 'image/png',
       width = tmp$width*input$multiplot_scale,
       height = tmp$height*input$multiplot_scale)
}, deleteFile = F)

observeEvent(input$remove_image,{
  vals$list_plots[[input$multiplot_image_to_modify]]<-NULL
  if(length(vals$list_plots)==0) {updateCheckboxInput(session,inputId = "show_multiplot",value = FALSE)}
  })

observeEvent(input$remove_all_images,{
  vals$list_plots<-list()
  updateCheckboxInput(session,inputId = "show_multiplot",value = FALSE)
})

}

