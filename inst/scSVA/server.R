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
    cluster_stat    =  data.frame(ID=0,Ncells=0,Mean=0.0,Median=0.0,SD=0.0, stringsAsFactors = F),
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
    inputFiles      =  list(XY=NULL,Genes=NULL,Clusters=NULL,Gene.Names=NULL),
    featureList     =  list(XY=NULL,Genes=NULL,Clusters=NULL,Gene.Names=NULL),
    ColPalette      =  NULL,
    InitColPalette  =  FALSE,
    Cols            =  FALSE
  )

  volumes = getVolumes()
  shinyFileChoose(input, 'inputFile_XY',          roots = volumes, filetypes=c('', 'h5','H5','hdf5','HDF5','loom','LOOM','csv','CSV','tsv','TSV','txt','TXT','gz'))
  shinyFileChoose(input, 'inputFile_Genes',       roots = volumes, filetypes=c('', 'h5','H5','hdf5','HDF5','loom','LOOM','gz'))
  shinyFileChoose(input, 'inputFile_Clusters',    roots = volumes, filetypes=c('', 'h5','H5','hdf5','HDF5','loom','LOOM','csv','CSV','tsv','TSV','txt','gz'))
  shinyFileChoose(input, 'inputFile_GenesVector', roots = volumes, filetypes=c('', 'csv','CSV','tsv','TSV','txt','gz'))

  updateNumericInput(session = session,inputId = "n_workers",value = parallel::detectCores())
  updateSelectizeInput(session,"pal",choices=c(dens_cols))

  observeEvent(input$inputFile_XY, {
    vals$inputFiles$XY  <- as.character(parseFilePaths(volumes, input$inputFile_XY)$datapath)
  })

  observeEvent(input$inputFile_Genes, {
    vals$inputFiles$Genes <- as.character(parseFilePaths(volumes, input$inputFile_Genes)$datapath)
  })

  observeEvent(input$inputFile_GenesVector, {
    vals$inputFiles$Gene.Names <- as.character(parseFilePaths(volumes, input$inputFile_GenesVector)$datapath)
  })

  observeEvent(input$inputFile_Clusters, {
    vals$inputFiles$Clusters <- as.character(parseFilePaths(volumes, input$inputFile_Clusters)$datapath)
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
    if(!is.null(vals$inputFiles$Genes)){
      if(input$inputFileType_Genes  ==  "HDF5"){
        choices  =  h5ls(vals$inputFiles$Genes)

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
                      caption  = paste0("File: ",vals$inputFiles$Genes))
      }
      else if(input$inputFileType_Genes  ==  "CSV"){
        DT::datatable(data.frame(data.frame(Files=c(paste0("Expression matrix: ",vals$inputFiles$Genes),
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

  observeEvent(input$colCluster, {
    vals$cluster_colors[input$cluster_id]  <-  input$colCluster
  })

  observeEvent(input$load_datasets, {
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
      }
      })
    }
#Cell Features
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
    gc()
    gc_python$collect()
    vals$X    <- c(py_to_r(ds$min('x')),py_to_r(ds$max('x')))
    vals$Y    <- c(py_to_r(ds$min('y')),py_to_r(ds$max('y')))
    vals$trans       =  TRUE
    vals$initialize  =  TRUE
    if(input$plot_expression){
      updateTextInput(session     = session,inputId = "choose_gene",    value ="")
      updateCheckboxInput(session = session,inputId = "plot_expression",value = FALSE)
    }

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
        withProgress(message = 'Getting clusters', value = 0, {
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

    ggplot(vals$data, aes(x = i, y = j,colour = if(input$dens_scale){log10(x)} else{x}))+
    {if(input$show_clusters & vals$CompClusFin) {
      geom_point(data    =  vals$cluster_df,
                 size    =  input$ps,
                 alpha   =  input$transparency,
                 shape   =  ifelse(input$point_shape,15,19),
                 colour  =  if(input$dichromat  !=  "No"){
                              dichromat(vals$cluster_colors[vals$cluster_df$Cluster_ID],
                                        input$dichromat)
                            }
                 else{vals$cluster_colors[vals$cluster_df$Cluster_ID]})
     }else{
      geom_point(size   =  input$ps,
                 alpha  =  input$transparency,
                 shape  =  ifelse(input$point_shape,15,19))
    }}+
    {if(input$highlight_cluster) geom_point(data    =  vals$cluster_h,
                                            size    =  input$ps*15,
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
   {if(input$plot_clusters){
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
    {scale_colour_gradientn(name     =  input$clegend,
                            limits   =  vals$Cols,
                            oob      =  squish,
                            colours  =  get_colors(input$pal,
                                                   input$dichromat  !=  "No",
                                                   input$dichromat,
                                                   vals$n.cols))}+
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
                                                                       c("i"  =  input$x_axis,"j"  =  input$y_axis))},
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
                                                               c("i"  =  input$x_axis,"j"  =  input$y_axis))},
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
      withProgress(message = 'Getting cluster label positions', value = 0, {
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

  observeEvent(input$show_clusters, {
    if(input$show_clusters){
      updateTextInput(session, "showlegend", value = FALSE)
      updateTextInput(session, "title",      value="")
      vals$plot.ID  =  3
      vals$trans    <- TRUE
    } else {
    updateTextInput(session, "showlegend", value = TRUE)
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
      showNotification("Load cluster labels first:  Metadata -> Clusters -> Get Cluster Label Positions",type = "warning",duration = 10)
      updateCheckboxInput(session, "get_exp_clusters", value = FALSE)
      return()
    }
    if(input$get_exp_clusters){
    withProgress(message = 'Getting Gene Expression by Cluster', value = 0, {
      for(l in 1:vals$N.clusters){
        ds$select(paste0("clusters==",l,""))
        vals$cluster_stat  <-  rbind(vals$cluster_stat,
                                     data.frame(ID      =  vals$cluster_labels$ID[l],
                                                Ncells  =  py_to_r(ds$count(selection=TRUE)),
                                                Mean    =  py_to_r(ds$mean("Gene",selection=TRUE)),
                                                Median  =  py_to_r(ds$median_approx('Gene',selection=TRUE)),
                                                SD      =  py_to_r(ds$std('Gene',selection=TRUE)),
                                                stringsAsFactors = F))
          incProgress(1/vals$N.clusters)
      }
      vals$cluster_stat        <-  vals$cluster_stat[-1,]
      vals$cluster_stat[,3:5]  <-  round(vals$cluster_stat[,3:5],digits=5)
    })
    } else {
      vals$cluster_stat  =  data.frame(ID     =  0,
                                       Ncells =  0,
                                       Mean   =  0.0,
                                       Median =  0.0,
                                       SD     =  0.0,
                                       stringsAsFactors = F)
    }
  })

  output$table_exp.in.clusters <- DT::renderDataTable({DT::datatable(vals$cluster_stat,
                                                         options  =  list(paging = FALSE,scrollY = "500px"),
                                                         caption  =  vals$loaded_feature,
                                                         filter   =  'bottom',
                                                         rownames =  FALSE) %>%
                                                         formatCurrency('Ncells',
                                                                        currency  = "",
                                                                        digits    = 0,
                                                                        interval  = 3,
                                                                        mark      = ","
                                                                        ) %>%
                                                         formatRound(c("Mean","Median","SD"), 3)})

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

  output$plot_DistGEinCluster = renderPlotly({
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
        if(input$plot_ecdf){
          p  <-  plot_ly()
          for(i in 1:length(s)){
            p  <-  add_lines(p,
                             x      =  seq(vals$rangeE[1],vals$rangeE[2],length.out=40),
                             y      =  cumsum(distCl[[i]])/sum(distCl[[i]]),
                             name   =  paste0("Cluster: ",vals$cluster_labels$ID[s[i]]),
                             marker =  list(color = vals$cluster_colors[s[i]]),
                             line   =  list(color = vals$cluster_colors[s[i]]),
                             mode   = 'lines')
          }
          p <- p %>% layout(yaxis  =  list(title = "Cumulative Probability"),
                            title  =  vals$loaded_feature,
                            xaxis  =  list(title = input$clegend))
          p
        } else {
          for(i in 1:length(s)){
            df=data.frame(Expression   =  seq(vals$rangeE[1],vals$rangeE[2],length.out=40),
                          Cell_Number  =  distCl[[i]]/sum(distCl[[i]]))
            plotList[i]  <-  plot_ly(df,
                                     x      = ~Expression,
                                     y      = ~Cell_Number,
                                     type   = "bar",
                                     name   = paste0("Cluster: ",vals$cluster_labels$ID[s[i]]),
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
        }
      } else {
        plotly_empty()
      }
    } else {
      plotly_empty()
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
        BiocParallel::register(BiocParallel::MulticoreParam(progressbar  =  T,
                                                            workers      =  input$n_workers))
        fun<-function(x){
          scale(unlist(fread(paste0("zq ",
                                    gsub(" ","\\ ",vals$inputFiles$Genes,fixed=T),
                                    " ",
                                    x,
                                    " |cut -f2-"),
                             stringsAsFactors = F)))
        }
        VectSignScores  <-  rowMeans(sapply(bplapply(GeneSignature,fun),unlist))
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
          BiocParallel::register(BiocParallel::MulticoreParam(progressbar = T,
                                                              workers     = input$n_workers))
          fun<-function(x){
            if(input$em_or  ==  "R"){
            return(scale(drop(rhdf5::h5read(vals$inputFiles$Genes,
                                            input$inputGenes,
                                            index  =  list(which(vals$GeneNames  ==  x),NULL)))))
            }
            else if(input$em_or  ==  "C"){
            return(scale(drop(rhdf5::h5read(vals$inputFiles$Genes,
                                            input$inputGenes,
                                            index  =  list(NULL,which(vals$GeneNames  ==  x))))))
            }
          }

          VectSignScores  <-  rowMeans(sapply(bplapply(GeneSignature,fun),unlist))
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
    if(input$plot_expression  ==  FALSE) {
      updateTextInput(session,"plot_expression",value=TRUE)
    }
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
}


