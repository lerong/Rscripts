wdpath = function(){
  HOME_drive = ifelse(.Platform$OS.type=="windows","//mymdafiles/usersdqs3/lli11","/home/lli11") 
  HN_drive = ifelse(.Platform$OS.type=="windows","//q1prphn/projects/lli11","/projects/lli11")
  Bioinfo2 = ifelse(.Platform$OS.type=="windows", "//d1prpccifs/bcb/bioinfo/bioinfo2", "/data/bioinfo2")
  assign("HOME_drive",HOME_drive,envir = .GlobalEnv)
  assign("HN_drive",HN_drive,envir = .GlobalEnv)
  assign("Bioinfo2",Bioinfo2,envir = .GlobalEnv)
}

createProj = function(proj_pi,proj_name){
  HN_proj_path = file.path(HN_drive,"Projects")
  create_date = gsub("-", "",Sys.Date(), fixed = TRUE)
  proj_name = paste(create_date,proj_name,sep="_")
  proj_path = file.path(HN_proj_path,proj_pi,proj_name)
  if(!file.exists(proj_path)){ 
    dir.create(proj_path,recursive = T) 
    folders = c("dataInput","dataOutput","Rcodes","refs","reports","figures","outputs","caches")
    folder_list = paste(proj_path,folders,sep=.Platform$file.sep)
    lapply(folder_list,dir.create)
    # copy templates files
    template_folder = file.path(HN_drive,"backups/templates/")
    files_list = paste(template_folder,dir(template_folder),sep=.Platform$file.sep)
    reports_folder = file.path(proj_path,"reports")
    file.copy(files_list,reports_folder)
    cat('New Project folder created for: \n',
        "-------------------------------------\n",
        proj_pi,'\n',proj_name,"\n",
        "@",proj_path ,"\n",basename(dir(proj_path)),"\n",
        "-------------------------------------\n")
  } else{ 
    cat('Folders already existed \n')  
  }
  list(pi = proj_pi,project = proj_name )
  
  assign("proj_pi",proj_pi,envir = .GlobalEnv)
  assign("proj_name",proj_name,envir = .GlobalEnv)
  assign("proj_path",proj_pi,envir = .GlobalEnv)
}

makeSubAnalysisFolders = function(projHome,subAnalysis){
  wdpath()
  dataInput_dir = file.path(projHome,"dataInput")
  dataOutput_dir = file.path(projHome,"dataOutput",subAnalysis)
  figOut_dir = file.path(projHome,"figures",subAnalysis)
  tabOut_dir = file.path(projHome,"outputs",subAnalysis)
  caches_dir = file.path(projHome,"caches",subAnalysis)
  createFolder(dataOutput_dir)
  createFolder(figOut_dir,echo = TRUE)
  createFolder(tabOut_dir,echo = TRUE)
  createFolder(caches_dir,echo = TRUE)
}






