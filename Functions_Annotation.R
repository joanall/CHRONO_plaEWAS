#------------------------------------------------------------------------------#
#                            Define functions                                  #----
#------------------------------------------------------------------------------#

get_annot_EPIraction <- function(gzip_file, out_folder) {
  
  # Go to EPIraction Tool and Get the form 
  url <- "https://epiraction.crg.eu/"
  epiraction_url <- read_html(url)
  
  forms <- html_form(epiraction_url)
  
  # We have 6 forms, if we go to the web browser we'll get to see the matching so 
  # 1 - interval --> genomic insights
  # 2 - rss_gene --> gene variants
  # 3 - overlap  --> gwas annotation (file) --> this is the one we want
  # 4 - enrichment  --> gwas annotation (gene name)
  # 5 - regions_list --> regions 
  # 6 - gene --> enhancer-gene pairs
  
  # Prepare form for file upload
  form_index <- 3  # The actual form index for "gwas annotation (file)"
  file_input_name <- "overlap"  # The name of the file input field
  file_path <- gzip_file
  form <- forms[[form_index]]
  
  # Check if the file exists
  if (!file.exists(file_path)) {
    stop("The file does not exist. Please check the file path.")
  } else { cat("File exists: ", file_path, "\n") }
  
  # Step 5: Create a session and set the form values using html_form_set
  session <- session(url)
  form <- html_form_set(form, !!file_input_name := file_path)
  print(form)
  
  # Prepare the form data for submission manually
  form_data <- list(overlap = upload_file(file_path))
  cat("Form data to be submitted:\n")
  print(form_data)
  
  # Step 6: Submit the form using POST request with multipart encoding
  response <- POST(url = form$action,
                   body = form_data,
                   encode = "multipart")
  
  # Check the response status
  if (response$status_code == 200) {
    cat("Form submission successful.\n")
    
    # Extract the HTML content
    content_text <- rawToChar(response$content)
    content_html <- read_html(content_text)
    
    # Extract tables
    tables <- html_table(content_html, fill = TRUE)
    
    if (length(tables) > 0) {
      print(tables[[1]]) # Parameters used (1st table it appears in the results web)
      
      # Extract the second table -> table with links to download 
      table2 <- tables[[2]]
      
      # Extract download links
      download_links <- content_html %>%
        html_nodes(xpath = "//table[2]//td/a[text()='Download']") %>%
        html_attr("href")
      
      # We want the first link -> all tissues as row with QTls 
      link <- download_links[1]
      download_url <- paste0("https://epiraction.crg.eu/", link)
      
      response <- GET(download_url)
      
      # Check if the request was successful
      if (response$status_code == 200) {
        destfile <- paste0(out_folder, "/", basename(link), ".txt")
        writeBin(response$content, destfile)
        cat("Downloaded and saved as: ", destfile, "\n")
      } else {
        cat("Failed to download: ", download_url, "\n")
      }
      
      return(destfile)
      
    } else {
      cat("No tables found in the response.\n")
    }
  } else {
    cat("Form submission failed with status code:", response$status_code, "\n")
    
    # Print the response content for debugging
    content_text <- rawToChar(response$content)
    cat("Response content:\n")
    cat(content_text)
  }
}

# Find Overlaps between two genomics ranges and returns the coincident values. 
# Search overlaps between CpGs (selected data) and Fetal Placenta (complete data)

findoverlapvalues <- function(selected, complete) {
  
  # Find overlaps
  ranges <- subsetByOverlaps(selected, complete)
  hits <- findOverlaps(selected, complete)
  
  # Find the names of states associated to each position
  idx <- subjectHits(hits)
  values <- DataFrame(States  = complete$name[idx])
  qh <- queryHits(hits)
  
  ans <- list("ranges" = ranges,
              "hits" = hits,
              "shits" = idx,
              "qhits" = qh,
              "values" = values)
  
  return(ans)}
