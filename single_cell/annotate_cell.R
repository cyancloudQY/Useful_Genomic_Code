##### SELF annotation ######

new.cluster.ids <- c(
  "0" = "Stromal Cells",
  "1" = "Supporting Cells",
  "2" = "Dark Cells",
  "3" = "NA",
  '4' = 'Macrophage', 
  '5' = 'Supporting Cells',
  '6' = 'Progenitor Cells',
  '7' = 'Hair Cells',
  '8' = 'Vascular Cells', 
  '9' = 'Roof Cells', 
  '10' = 'Schwann Cells',
  '11' = 'Supporting Cells',
  '12' = 'NA', 
  '13' = 'Pericytes',
  '14' = 'Melanocytes',
  '15' = 'NA',
  '16' = 'T cells'
)


Idents(dc_subject) <- recode(Idents(dc_subject), !!!new.cluster.ids)
dc_subject$cell_type <- Idents(dc_subject)


####
