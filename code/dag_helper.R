
#!/usr/bin/env Rscript

# getting arguments from terminal command

args <- commandArgs(trailingOnly = TRUE)

cluster_id <- args[1]
num_jobs <- as.numeric(args[2])
email <- args[3]

# function to write dag file

write_dag <- function(cluster_id, num_jobs, email) {
  
  dag_txt <- character(0)
  
  for(i in 0:(num_jobs-1)) {dag_txt <- c(dag_txt, paste0("JOB job_", i, " ", cluster_id, ".", i))}
  
  parent_line <- paste("PARENT", paste(paste0("job_", 0:(num_jobs-1)), collapse = " "), "-> all_done")
  sub_line <- "JOB all_done notify.sub"
  notif_line <- "NOTIFICATION = Complete"
  email_line <- paste0("NOTIFY_USER = ", email)

  dag_txt <- c(dag_txt, "", parent_line, "", sub_line, notif_line, email_line)
  
  return(dag_txt)
}

# write text to dag file

output_file <- "notify_batch.dag"

writeLines(write_dag(cluster_id, num_jobs, email), output_file)

