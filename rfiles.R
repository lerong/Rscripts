#-------------------------------------------------------------------------------------------------------------------------
# rfunctions for project organizations
# ---projCreate.R
#--------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------
# create the project directories
#--------------------------------------------------------------------------------------------------------------------------
projCreate = function(proj_pi,proj_name){
  # create project folder at current folder
  HOME_drive = ifelse(.Platform$OS.type=="windows","//mymdafiles/usersdqs3/lli11","/home/lli11") 
  HN_drive = ifelse(.Platform$OS.type=="windows","//q1prphn/projects/lli11","/projects/lli11")
  HN_proj_path = file.path(HN_drive,"Projects")
  # pwd = getwd()
  create_date = gsub("-", "",Sys.Date(), fixed = TRUE)
  proj_name = paste(create_date,proj_name,sep="_")
  proj_path = file.path(HN_proj_path,proj_pi,proj_name)
  if(!file.exists(proj_path)){ 
    dir.create(proj_path,recursive = T) 
    #setwd(proj_path)  
    folders = c("data","Rcodes","refs","reports","figures","outputs")
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
    warning('Folder already existed \n')  
  }
  
  tmp = list()
  tmp$PID = proj_pi
  tmp$PNAME = proj_name
  tmp
}

