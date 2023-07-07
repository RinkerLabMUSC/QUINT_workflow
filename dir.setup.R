#' @title QUINT Workflow Automated Directory Setup
#' @description  Creates a set of directories for the QUINT workflow based on a vector of unique IDs like animal numbers/names.
#' @param IDs (required) Vector, set of unique identifiers assigned to a given brain, usually animal IDs.
#' @param project.path (required) String, file directory where you would like the new directories to be created. Use forward slashes.
#' @param fluorophore (required) Vector, the name of the fluorophore(s) for your signal channel.
#'          If you have multiple signal channels, a vector of fluorphore names equal to the number of signal channels is required.
#'          If you are performing colocalization, each overlay or merge counts as a signal channel of its own and must be included.
#' @export

dir.setup <- function(IDs, project.path, fluorophore){

  #check if there is a '/' at the end of the project filepath and add one if there is not.
  if(substr(project.path, nchar(project.path), nchar(project.path)) != "/"){
    project.path <- paste0(project.path,"/")
  }

  if(dir.exists(project.path) == T &&
     dir.exists(paste0(project.path,"QUINT/")) != T){

    dir.create(paste0(project.path,"QUINT/"))
    dir.create(paste0(project.path,"QUINT/QUINTegrate Input/"))
    dir.create(paste0(project.path,"QUINT/QUINTegrate Output/"))
    dir.create(paste0(project.path,"QUINT/Subjects/"))
      for(i in IDs){
        ID.path <- paste0(project.path, "/QUINT/Subjects/", i, "/") #Create path for unique ID
        #Check to see if this directory already exists so you dont overwrite data
        if(dir.exists(ID.path) != T){
          #Create directories for this ID
          dir.create(paste0(ID.path))
          dir.create(paste0(ID.path,"AtlasMap/"))
          dir.create(paste0(ID.path,"Output/"))
          dir.create(paste0(ID.path,"Segmentations/"))
          dir.create(paste0(ID.path,"Segmentations/Simple/"))
          dir.create(paste0(ID.path,"Segmentations/Probabilities/"))
          dir.create(paste0(ID.path,"Images/"))
          dir.create(paste0(ID.path,"Images/JPGs/"))
          dir.create(paste0(ID.path,"Images/TIFFs/"))
          dir.create(paste0(ID.path,"Images/TIFFs/Registration Channel/"))
          for(j in fluorophore){
            dir.create(paste0(ID.path,"Images/TIFFs/", j, "/"))
          }
        } else {cat(paste0("For: ", ID.path, "\nDirectory already exists. There is nothing to create."))}
      }

  } else {stop(paste0("Error either: ", project.path, " does not exist. OR\n", project.path,"QUINT/ ", "already exists in ", project.path))}

}
