# Imputation 

dataset_og <- readRDS('./raw_data_AFAR.rds')

# At least one observation of every variable
dataset <- dataset_og[which(rowSums(is.na(dataset_og)) < 373),]

print(nrow(dataset))

# Random forest imputation 
start.time <- Sys.time()
imp_rf <- mice(dataset, method = 'rf', m = 1, maxit = 40, ntree = 10) 
end.time <- Sys.time()

time.taken <- round(end.time - start.time,2)
time.taken

imp_rf$loggedEvents

saveRDS(imp_rf, file.path('./', paste0('impset_AFAR_',as.character(Sys.Date()),'.rds')))

pdf('QC-imputation_AFAR.pdf')
for (v in names(dataset)) { if (nrow(imp_rf$imp[[v]]) > 1) {
  message(v)
  
  nmiss <- sum(is.na(imp_rf$data[v]))
  nmiss <- paste0('\n n missing = ', nmiss, ' (',round(nmiss/nrow(imp_rf$data)*100,1),'%)')
  
  try(print(mice::densityplot(imp_rf, as.formula(paste('~',v)), main=paste(v, nmiss))))
}
}
dev.off()

# ------------------------------------------------------------------------------
# imp_rf <- readRDS(list.files(path = datapath, pattern = 'impset_AFAR_', full.names = TRUE))
