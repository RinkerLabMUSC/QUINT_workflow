#' @title Integrate 'QUINT' Workflow Datasets Respective to Brain Structure Hierarchy
#' @description  Integrates datasets from the 'QUINT' workflow with additional details from the 'Allen CCFv3 Ontogeny'.
#    The main function is to sum Object count, Object area, Object pixels, Region pixels and Region area according
#    to each region's position in the hierarchy. All sub-nuclei data is summed to produce a total for each parent
#    region. E.G. 10 cells in the nucleus accumbens come from summing the child structures 2 cells in the nucleus
#    accumbens core and 8 in the nucleus accumbens shell. This "summing up the hierarchy" is applied to all structures
#    and sub-structures and then output into a list object and saved as a .csv file.
#' @param project.name (required) String, Name or abbreviation for your project for use in naming outputs
#' @param input.path (required) String, file directory containing only unaltered Nutil output files in html format
#' @param save.path (required) String, file directory where outputs should be saved.
#' @return A list object including the fully summed and integrated data frame.
#' @export
QUINTegrate <- function(project.name,input.path,save.path){
  # Check to see if directories exist:
  if(dir.exists(input.path) == T &&
     dir.exists(save.path) == T &&
     endsWith(input.path, "/") == T &&
     endsWith(save.path, "/") == T){
    # Create a list object to store relevant data
    main <- list()

    # Store directories
    # Make a list of directories for the raw cell count data of only files with the .html file extension at the end in the specified directory
    main$dir$html <- dir(input.path, pattern = "\\.html$", full.names = T)
    # Write save path for workspace and data saving into list
    main$dir$save <- save.path
    # Write project name into list
    main$project <- project.name

    #loop to import and rename raw data
    for(i in 1:length(main$dir$html)){
      html<- rio::import(main$dir$html[[i]])[,1:10] #Keeps only the first 10 columns
      html[1,9] <- "Area.unit.1"
      colnames(html) <- gsub(" ", ".", html[1,]) #ads "." for every space and uses first row as column names
      html <- html[-1,]
      #change the data types of columns
      html[,1] <- as.integer(html[,1])
      html[,3] <- as.integer(html[,3])
      html[,4] <- as.integer(html[,4])
      html[,6] <- as.integer(html[,6])
      html[,7] <- as.integer(html[,7])
      html[,8] <- as.integer(html[,8])

      cat(paste0("Reading in ", basename(main$dir$html[[i]]), "\n"))
      # Check that you have the correct columns in number and in names. Remove extraneous columns and data registered outside the atlas.
      # Check that all the necessary columns are present.
      if("Region.ID" %in% colnames(html)
         && "Region.Name" %in% colnames(html)
         && "Region.pixels" %in% colnames(html)
         && "Region.area" %in% colnames(html)
         && "Area.unit" %in% colnames(html)
         && "Object.count" %in% colnames(html)
         && "Object.pixels" %in% colnames(html)
         && "Object.area" %in% colnames(html)
         && "Area.unit.1" %in% colnames(html)){
        cat("Column names checked. All necessary columns present.\n")
      } else {
        stop("Essential columns or column names are missing from the animal's data file.
       \n**These files should NOT be altered after being created by Nutil**
       \nEach animal's data file MUST include the following column names:
       \nRegion ID, Region Name, Region pixels, Region area, Area unit, Object count, Object pixels, Object area, Area unit.
       \nCheck the column names before trying again.\n")
      }
      # If the first row contains "Clear Label" which is all cells registered outside the atlas, then remove this column from the dataset.
      if(html$Region.Name[1] == "Clear Label"){
        html<-html[-1,]
        cat(paste0("Clear Label row removed from ", basename(main$dir$html[[i]]),"\n"))
      }
      # If there are any extra junk columns that dont have names and are thus named by r with X1, X2, etc... remove these columns from the data
      # Check the vector of column names to see if any contain "X" and if they do, return a TRUE at that position in a vector.
      if(all(stringr::str_detect(colnames(html),"X")!=F)){
        # Get a vector of numeric column positons of all those column names that contain X and subset them out to remove them from html
        html <- html[,-which(stringr::str_detect(colnames(html),"X") %in% T)]
        cat(paste0("Junk columns removed from ", basename(main$dir$html[[i]]),"\n"))
      }

      # Add animalID info column to df
      html[,length(colnames(html))+1] <- tools::file_path_sans_ext(basename(main$dir$html[[i]]))
      colnames(html)[[length(colnames(html))]] <- "animalID"

      # Store html data in list object
      main$html[[i]] <- html
      names(main$html)[[i]] <- gsub("-",".",tools::file_path_sans_ext(basename(main$dir$html[[i]]))) #substitues the - for . in the filename without the extension.
      remove(html)
    }

    # Merge Atlas Ontogeny with cell count data

    for(i in 1:length(main$html)){
      merged <- merge(x = main$html[[i]], y = Allen.CCFv3.Ontog,
                      by.x = c("Region.ID", "Region.Name"),
                      by.y = c("id", "name"),
                      all = T)

      main$merged[[i]] <- merged
      names(main$merged)[[i]] <- gsub("-",".",tools::file_path_sans_ext(basename(main$dir$html[[i]]))) #substitues the - for . in the filename without the extension.
      remove(merged)
    }

    cat("Allen CCFv3 Ontogeny merged with user data.\n")



    ##Sum cell counts or other variables for cases where the parent structure is the same, then input that value into the appropriate column in the parent structure's row.
    cat("Child region areas, pixels, and counts have been summed and input into respective parent region cells for:\n")
    for(i in seq_along(main$merged)){ #mouse level

      #rearrange the data to start with the lowest depth so no higher depths are summed before their lower child depths.
      main$merged[[i]] <- dplyr::arrange(main$merged[[i]], -depth) #-depth organizes from highest to lowest value here

      main$summed[[i]] <- main$merged[[i]] #write to new object where it will be manipulated.
      names(main$summed)[[i]] <- gsub("-",".",tools::file_path_sans_ext(basename(main$dir$html[[i]]))) #substitutes the - for . in the filename without the extension.

      for(j in main$summed[[i]]$Region.ID){

        # Object Count
        OCsum <- sum(main$summed[[i]][which(main$summed[[i]]$parent.structure.id == j), "Object.count"])
        if(main$summed[[i]]$Object.count[which(main$summed[[i]]$Region.ID == j)] == 0){
          main$summed[[i]]$Object.count[which(main$summed[[i]]$Region.ID == j)] <- OCsum
        }
        # Object Pixels
        OPsum <- sum(main$summed[[i]][which(main$summed[[i]]$parent.structure.id == j), "Object.pixels"])
        if(main$summed[[i]]$Object.pixels[which(main$summed[[i]]$Region.ID == j)] == 0){
          main$summed[[i]]$Object.pixels[which(main$summed[[i]]$Region.ID == j)] <- OPsum
        }
        # Object Area
        OAsum <- sum(main$summed[[i]][which(main$summed[[i]]$parent.structure.id == j), "Object.area"])
        if(main$summed[[i]]$Object.area[which(main$summed[[i]]$Region.ID == j)] == 0){
          main$summed[[i]]$Object.area[which(main$summed[[i]]$Region.ID == j)] <- OAsum
        }
        # Region Pixels
        RPsum <- sum(main$summed[[i]][which(main$summed[[i]]$parent.structure.id == j), "Region.pixels"])
        if(main$summed[[i]]$Region.pixels[which(main$summed[[i]]$Region.ID == j)] == 0){
          main$summed[[i]]$Region.pixels[which(main$summed[[i]]$Region.ID == j)] <- RPsum
        }
        # Region Area
        RAsum <- sum(main$summed[[i]][which(main$summed[[i]]$parent.structure.id == j), "Region.area"])
        if(main$summed[[i]]$Region.area[which(main$summed[[i]]$Region.ID == j)] == 0){
          main$summed[[i]]$Region.area[which(main$summed[[i]]$Region.ID == j)] <- RAsum
        }
      }

      main$summed[[i]] <- dplyr::arrange(main$summed[[i]], graph.order)
      #print progress
      cat(paste0(names(main$merged)[[i]]))

    }


    #use rbind to append each dataframe by rows
    #write the first main$summed dataframe to a temporary object.
    df <- main$summed[[1]]

    #loop along each dataframe in main$summed excluding the first one which is already in the df at the start of the loop.
    for(i in seq_along(main$summed)[-1]){

      df <- rbind(df, main$summed[[i]])

    }

    #write completed df to object in our list structure.
    main$data <- df
    #add column with row labels to arrange the data back into its original order later.
    main$data <- dplyr::mutate(main$data, order = row_number())

    cat("Individual data bound into a single main dataset. Call $data for the data frame.","\n")

    # Generate .csv file for main dataset
    write.table(main$data, file = paste0(main$dir$save, main$project, " Cell Count Data.csv"), sep = ",", row.names = F)

    cat(paste0("Dataset saved to: ", main$dir$save, main$project, " Cell Count Data.csv","\n"))

  } else {
    stop("File path does not exist.")
  }
  return(main)
}#End of QUINTegrate function


