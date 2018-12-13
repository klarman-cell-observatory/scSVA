#shiny::runApp("inst/scSVA",launch.browser =TRUE)
require(reticulate)
require(rhdf5)
require(Matrix)
require(ggplot2)
require(Cairo)
require(MASS)
library(RColorBrewer)
require(colorRamps)
require(dichromat)
require(viridis)
require(colourpicker)
require(shiny)
require(shinydashboard)
require(shinythemes)
require(shinyFiles)
require(scales)
require(sp)
require(extrafont)
require(plyr)
require(dplyr)
require(ggrepel)
require(plotly)
require(plot3D)
require(DT)
require(data.table)
require(tableHTML)
require(future)
require(future.apply)
require(shinyBS)
require(shinyjqui)
require(shinyAce)
require(gridExtra)
#require(magick)
#require(googleComputeEngineR)

suppressMessages(loadfonts())

if(Sys.info()["sysname"]=="Darwin" & file.exists("~/.bash_profile")){
  reticulate::use_python(dirname(system("source ~/.bash_profile; which python3",intern = T)))
} else if(Sys.info()["sysname"]=="Linux"){
  reticulate::use_python(dirname(system("which python3",intern = T)))
}

if(file.exists("/.dockerenv")){Sys.setenv(PATH = paste("/opt/conda/bin/", Sys.getenv("PATH"),sep=":"))}

if(Sys.getenv("GOOGLE_CLOUD_SDK") == ""){
GCU<-if(Sys.info()["sysname"]=="Darwin" & file.exists("~/.bash_profile")){
  system("source ~/.bash_profile; which gsutil",intern = T)
} else if(Sys.info()["sysname"]=="Linux"){
    system("which gsutil",intern = T)
}
} else {
  GCU<-file.path(Sys.getenv("GOOGLE_CLOUD_SDK"),"gsutil")
}

options(future.globals.maxSize=Inf)
plan(multiprocess)

options(warn=-1)
options(shiny.usecairo=FALSE)

numpy     <-  NULL
vaex      <-  NULL
gc_python <-  NULL


vaex       <-reticulate::import("vaex", convert  =  FALSE)
numpy      <-reticulate::import("numpy", convert =  FALSE)
gc_python  <-reticulate::import("gc", convert  =  FALSE)

     
InitPlot<-read.table(gzfile(paste0(system.file("scSVA",package="scSVA"),"/scSVA.txt.gz")),
                     header = T,
                     sep    = ",")

ds <-  vaex$from_arrays(x        = numpy$array(InitPlot$X),
                        y        = numpy$array(InitPlot$Y),
                        Gene     = numpy$array(InitPlot$Gene),
                        clusters = numpy$array(InitPlot$Cluster))

termOut <- NULL
vm      <- NULL
tag     <- NULL

Parameters.FA_3D <- NULL
FLE_dist         <- NULL
Parameters.DMap  <- NULL
Parameters.NNG   <- NULL

modes <- getAceModes()
themes <- getAceThemes()
tmpdir<-tempdir()
#https://en.wikipedia.org/wiki/List_of_Crayola_crayon_colors
list.colors<-list(
  Crayola.12  =  c("#000000","#0066FF","#AF593E","#01A368","#FF861F","#ED0A3F","#FF3F34","#76D7EA","#8359A3","#FFFFFF","#FBE870","#C5E17A"),
  Crayola.24  =  c("#000000","#0066FF","#AF593E","#01A368","#FF861F","#ED0A3F","#FF3F34","#76D7EA","#8359A3","#FFFFFF","#FBE870","#C5E17A",
                   "#03BB85","#FFDF00","#8B8680","#0A6B0D","#8FD8D8","#A36F40","#F653A6","#CA3435","#FFCBA4","#FF99CC","#FA9D5A","#FFAE42"),
  Crayola.36  =  c("#000000","#0066FF","#AF593E","#01A368","#FF861F","#ED0A3F","#FF3F34","#76D7EA","#8359A3","#FFFFFF","#FBE870","#C5E17A",
                   "#03BB85","#FFDF00","#8B8680","#0A6B0D","#8FD8D8","#A36F40","#F653A6","#CA3435","#FFCBA4","#FF99CC","#FA9D5A","#FFAE42",
                   "#A78B00","#788193","#514E49","#1164B4","#F4FA9F","#FED8B1","#C32148","#01796F","#E90067","#FF91A4","#404E5A","#6CDAE7"),
  Crayola.50  =  c("#000000","#0066FF","#AF593E","#01A368","#FF861F","#ED0A3F","#FF3F34","#76D7EA","#8359A3","#FFFFFF","#FBE870","#C5E17A",
                   "#03BB85","#FFDF00","#8B8680","#0A6B0D","#8FD8D8","#A36F40","#F653A6","#CA3435","#FFCBA4","#FF99CC","#FA9D5A","#FFAE42",
                   "#A78B00","#788193","#514E49","#1164B4","#F4FA9F","#FED8B1","#C32148","#01796F","#E90067","#FF91A4","#404E5A","#6CDAE7",
                   "#FFC1CC","#006A93","#867200","#E2B631","#6EEB6E","#FFC800","#CC99BA","#00003B","#BC6CAC","#DCCCD7","#EBE1C2","#A6AAAE",
                   "#B99685","#0086A7"),
  Crayola.60  =  c("#000000","#0066FF","#AF593E","#01A368","#FF861F","#ED0A3F","#FF3F34","#76D7EA","#8359A3","#FFFFFF","#FBE870","#C5E17A",
                   "#03BB85","#FFDF00","#8B8680","#0A6B0D","#8FD8D8","#A36F40","#F653A6","#CA3435","#FFCBA4","#FF99CC","#FA9D5A","#FFAE42",
                   "#A78B00","#788193","#514E49","#1164B4","#F4FA9F","#FED8B1","#C32148","#01796F","#E90067","#FF91A4","#404E5A","#6CDAE7",
                   "#FFC1CC","#006A93","#867200","#E2B631","#6EEB6E","#FFC800","#CC99BA","#00003B","#BC6CAC","#DCCCD7","#EBE1C2","#A6AAAE",
                   "#B99685","#0086A7","#5E4330","#C8A2C8","#708EB3","#BC8777","#B2592D","#497E48","#6A2963","#E6335F","#03228E","#B5A895"),
  Crayola.105 =  c("#000000","#0066FF","#AF593E","#01A368","#FF861F","#ED0A3F","#FF3F34","#76D7EA","#8359A3","#FFFFFF","#FBE870","#C5E17A",
                   "#03BB85","#FFDF00","#8B8680","#0A6B0D","#8FD8D8","#A36F40","#F653A6","#CA3435","#FFCBA4","#FF99CC","#FA9D5A","#FFAE42",
                   "#A78B00","#788193","#514E49","#1164B4","#F4FA9F","#FED8B1","#C32148","#01796F","#E90067","#FF91A4","#404E5A","#6CDAE7",
                   "#FFC1CC","#006A93","#867200","#E2B631","#6EEB6E","#FFC800","#CC99BA","#00003B","#BC6CAC","#DCCCD7","#EBE1C2","#A6AAAE",
                   "#B99685","#0086A7","#5E4330","#C8A2C8","#708EB3","#BC8777","#B2592D","#497E48","#6A2963","#E6335F","#03228E","#B5A895",
                   "#0048ba","#EED9C4","#C88A65","#FDD5B1","#D0FF14","#B2BEB5","#926F5B","#00B9FB","#6456B7","#DB5079",
                   "#C62D42","#FA9C44","#DA8A67","#FD7C6E","#93CCEA","#FCF686","#503E32","#63B76C","#FF5470","#87421F",
                   "#9DE093","#FF7A00","#4F69C6","#A50B5E","#F0E68C","#FDFF00","#F091A9","#FF6E4A","#2D383A","#6F9940",
                   "#FC74FD","#652DC1","#D6AEDD","#EE34D2","#BB3385","#6B3FA0","#33CC99","#FFDB00","#87FF2A","#FC80A5",
                   "#D9D6CF","#00755E","#FFFF66","#7A89B8","#FF007C"),
  Crayola_Mix =  c("#000000","#95918C","#FDD9B5","#1F75FE","#FFAE42","#FFAACC","#1CAC78","#1DACD6","#7366BD","#FC2847",
                   "#50BFE6","#00468C","#C3CDE6","#0066FF","#91E351","#C0E7F1","#6ADA8E","#FFF7CC","#FFB97B","#FDFF00",
                   "#0048BA","#D0C6C6","#6E7FE7","#FF7A00","#6D3834","#0D98BA","#E7C697","#FF355E","#54626F","#66FF66",
                   "#AAF0D1","#FF00CC","#FFCC33"),
Crayola_Fluorescent    =  c("#FF355E","#FD5B78","#FF6037","#FF9966","#FF9933","#FFCC33","#FFFF66","#CCFF00","#66FF66","#AAF0D1","#50BFE6","#FF6EFF","#EE34D2","#FF00CC"),
Crayola_Silver_Swirls  =  c("#C39953","#AD6F69","#6EAEA1","#B768A2","#9E5E6F","#5F8A8B","#A17A74","#2E2D88","#AE98AA","#8BA8B7","#DA2C43","#914E75",
                            "#6D9BC3","#AB92B3","#BBB477","#5DA493","#778BA5","#8A496B","#CD607E","#676767","#AD4379","#A6A6A6","#5FA778","#56887D"),
Crayola_Magic_Scent    =  c("#FFFFFF","#CA3435","#FBE870","#8359A3","#000000","#0066FF","#ED0A3F","#FED85D","#C32148","#9E5B40","#FF8833","#C9A0DC",
                            "#FFA6C9","#FF3399","#4570E6","#AF593E","#29AB87","#C5E17A","#FFCBA4","#8B8680","#FC80A5","#76D7EA","#FDD5B1","#01786F","#C3CDE6"),
Crayola_Changeables    =  c("#C0E7F1","#91E351","#FF8071","#FF8ABA","#F4405D","#FDFD07","#EB58DD","#963D7F","#000000","#FFF7CC","#131391","#4F7948","#FFE9D1"),
Crayola_Pearl_Brite    =  c("#5FBED7","#E8F48C","#4F42B5","#F1444A","#54626F","#F37A48","#48BF91","#F2F27A","#6ADA8E","#702670","#7B4259","#F1CC79","#F5F5F5",
                            "#D65282","#F03865","#3BBCD0"),
Crayola_Glitter        =  c("#000000","#0D98BA","#1CAC78","#FF7538","#EE204D","#7851A9","#FFFFFF","#FCE883","#1F75FE","#FFAACC","#C8385A","#E6A8D7","#C0448F",
                            "#80DAEB","#C5E384"),
Crayola_Metallic_FX    =  c("#C46210","#2E5894","#9C2542","#BF4F51","#A57164","#58427C","#4A646C","#85754E","#319177","#0A7E8C","#9C7C38","#8D4E85","#8FD400",
                            "#D98695","#757575","#0081AB"),
Crayola_Gel_FX         =  c("#FFBF7F","#00D0B5","#7853A8","#63A5C3","#0081FF","#FF3399","#CF0060","#FFFFFF","#FF6699","#8F5873","#F26D7D","#6666CC","#F58345","#FFFF66",
                            "#99FF99","#12E3DB"),
Crayola_Silly_Scents   =  c("#C5E17A","#FED85D","#FF681F","#00CCCC","#D99A6C","#B94E48","#6456B7","#FF8833","#ECEBBD","#ED0A3F","#FC80A5","#E77200",
                            "#C32148","#8B8680","#F7468A","#76D7EA"),
Crayola_Heads_n_Tails  =  c("#FF3855","#FFAA1D","#2243B6","#A83731","#FD3A4A","#FFF700","#5DADEC","#AF6E4D","#FB4D46","#299617","#5946B2","#1B1B1B","#FA5B3D",
                            "#A7F432","#9C51B6","#BFAFB2"),
Crayola_Mini_Twistables=  c("#FDD9B5","#000000","#1F75FE","#0D98BA","#7366BD","#B4674D","#FFAACC","#1DACD6","#FDDB6D","#95918C","#1CAC78","#F0E891","#5D76CB",
                            "#FF7538","#EE204D","#FF5349","#C0448F","#FC2847","#926EAE","#F75394","#FFFFFF","#FCE883","#C5E384","#FFAE42")
)
tmp           <-  dichromat::colorschemes
names(tmp)    <-  paste0("dichromat::",names(dichromat::colorschemes))
tmp.2         <-  lapply(1:nrow(brewer.pal.info), function(x) brewer.pal(name=rownames(brewer.pal.info)[x],n=brewer.pal.info$maxcolors[x]))
names(tmp.2)  <-  rownames(brewer.pal.info)
list.colors   <-  c(list.colors,tmp,tmp.2)
dens_cols     <-  c("inferno", "viridis", "magma", "rainbow","matlab.like","matlab.like2","matlab.like3","gray.red",names(tmp),names(tmp.2))
rm(tmp,tmp.2)
gray.red      <-  function(n.colors){colorRampPalette(c("grey80","gray70","gray60","red","red2","red4","black"))(n.colors)}
matlab.like3  <-  function(n.colors){c("black","gray20",matlab.like2(n.colors-2))}
get_colors    <-  function(palette,if.dichromat,d.type,n.colors){
  if(!is.null(tryCatch(get(palette), error  =  function(e) {}))){
    if(if.dichromat) {dichromat(do.call(palette,list(n.colors)),d.type)} else
    {do.call(palette,list(n.colors))}
  } else {
    if(if.dichromat) {dichromat(colorRampPalette(list.colors[[palette]])(n.colors),d.type)} else
      {colorRampPalette(list.colors[[palette]])(n.colors)}
  }
}

plot_palette  <-  function(list.cols,
                           border = "gray80"){
    plot(0, 0,
         axes =  FALSE,
         xlab =  "",
         ylab =  "",
         xlim =  c(0L, 1L),
         ylim =  c(0L, 1L),
         type =  "n")
  N  =  length(list.cols)
  rect(0:(N-1)/N, 0, 1:N/N, 1.0,
       col     =  list.cols,
       border  =  border)
}

#Violin plots
StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                       tmp<-rbind(data.frame(x=data$y,y=data$x),
                                  data.frame(x=-rev(data$y),y=rev(data$x)))
                       #   tmp<-rbind(data.frame(x=tmp$x,y=tmp$y,col=rep(0,length(tmp$x))),
                       #           data.frame(x=tmp$x,y=tmp$y,col=rep(1,length(tmp$x))),
                       #      )
                       #   tmp$col=factor(tmp$col)
                       tmp
                     },
                     default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                                       alpha = NA, linetype = "solid"),
                     
                     required_aes = c("x", "y")
)

stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}










