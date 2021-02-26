## Functions for Pocillopora microfragment physiology data analysis 
library("plyr")
library("tidyverse")

###

### Ana's import functions
Import_YII <- function (dir)  {
# 1. Read all the .csv files (; delimited) and consolidate them in a single data frame
  # dir: folder containing the .csv files. The data in the dir can be separated in subfolders, recursive = TRUE will find the .csv files inside the subfolders
  # You can use the subfolders names to later extract information with the option full.names=TRUE, in this 
  # example the subfolders are the differents sampling dates, but it could be treatments, locations...

   IPAM.data <- plyr::ldply (list.files(path = dir, pattern = "*.csv", 
                                full.names = TRUE, recursive = TRUE),
                     function(filename) {
                       dum = (read.csv(filename, sep = ";", blank.lines.skip = T))
                       dum$filename = filename
                       return(dum)
                     })  

# 2. Select the columns that contain the file path and the YII values (remove F0 and Fm values)
  IPAM.data.Y.II <- IPAM.data[,grep("filename|Y.II.", names(IPAM.data))]
  IPAM.data.F0 <- IPAM.data[,grep("filename|FO", names(IPAM.data))]
  IPAM.data.Fm <- IPAM.data[,grep("filename|Fm", names(IPAM.data))]
  
  # Change to long format
  IPAM.data.Y.II <- na.omit(reshape2::melt(IPAM.data.Y.II, id=c("filename")))
  IPAM.data.F0 <- na.omit(reshape2::melt(IPAM.data.F0, id=c("filename")))
  IPAM.data.Fm <- na.omit(reshape2::melt(IPAM.data.Fm, id=c("filename")))
  
# 3.If you have subfolders separate folder and file information contained in "filename"
  # In this example each subfolder correspond to different dates, but it could be treatments, locations...
  # Followed for the specific file name
  file.picture <- plyr::rbind.fill(lapply(strsplit(as.character(IPAM.data$filename), split="/"), 
                                    function(X) data.frame(t(X))))
  colnames(file.picture) <- c("Folder", "Date", "File")
  
  YII <- cbind(file.picture[,-1], IPAM.data.Y.II[,-1])
  F0 <- cbind(file.picture[,-1], IPAM.data.F0[,-1])
  Fm <- cbind(file.picture[,-1], IPAM.data.Fm[,-1])
  
  return(YII)
}

Import_F0 <- function (dir)  {
  # 1. Read all the .csv files (; delimited) and consolidate them in a single data frame
  # dir: folder containing the .csv files. The data in the dir can be separated in subfolders, recursive = TRUE will find the .csv files inside the subfolders
  # You can use the subfolders names to later extract information with the option full.names=TRUE, in this 
  # example the subfolders are the differents sampling dates, but it could be treatments, locations...
  
  IPAM.data <- plyr::ldply (list.files(path = dir, pattern = "*.csv", 
                                       full.names = TRUE, recursive = TRUE),
                            function(filename) {
                              dum=(read.csv(filename, sep = ";", blank.lines.skip=T))
                              dum$filename=filename
                              return(dum)
                            })  
  
  # 2. Select the columns that contain the file path and the YII values (remove F0 and Fm values)
  IPAM.data.F0 <- IPAM.data[,grep("filename|FO", names(IPAM.data))]
  
  # Change to long format
  IPAM.data.F0 <- na.omit(reshape2::melt(IPAM.data.F0, id=c("filename")))
  
  # 3.If you have subfolders separate folder and file information contained in "filename"
  # In this example each subfolder correspond to different dates, but it could be treatments, locations...
  # Followed for the specific file name
  file.picture <- plyr::rbind.fill(lapply(strsplit(as.character(IPAM.data$filename), split="/"), 
                                          function(X) data.frame(t(X))))
  colnames(file.picture) <- c("Folder", "Date", "File")
  
  F0 <- cbind(file.picture[,-1], IPAM.data.F0[,-1])
  
  return(F0)
}
  