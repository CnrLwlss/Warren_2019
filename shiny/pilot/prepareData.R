library(data.table)

pdat = fread("../pilot_data.csv",sep=",",stringsAsFactors=FALSE)
grp = ifelse(substring(pdat$Folder,1,1)=="C","Control","Patient")
draw = data.frame(value=pdat$Value,id=pdat$ID,channel=pdat$Channel,patient_type=grp,replicate=1,subject_group=grp,patrep_id = paste(pdat$Folder,sprintf("R%02d",1),sep=""),stringsAsFactors=FALSE)
draw$cell_id = paste(draw$patrep_id,sprintf("%04d",draw$id),sep="_")
draw$cluster = 1

chans = c("133Cs_133Cs", "134Xe_134X", "138Ba_138Ba", "153Eu_SDHA", "155Gd_TOMM20", 
"158Gd_158Gd", "160Gd_NDUFB8", "161Dy_OSCP", "164Dy_GRIM19", 
"166Er_VDAC1", "168Er_Cox4", "172Yb_MTCO1", "174Yb_UqCRC2", "176Yb_Dystrophin", 
"191Ir_191Ir-DNA1", "193Ir_193Ir-DNA2", "195Pt_195P", "208Pb_208Pb", 
"80ArAr_80ArAr", "89Y_89Y")

chans = c("Cs", "Xe", "Ba", "SDHA", "TOM20","Gd", "NDUFB8", "OSCP", "NDUFA13", 
"VDAC1", "COX4+4L2", "MTCO1", "UqCRC2", "Dystrophin","DNA1", "DNA2", "Pt", "Pb","Ar", "Y")

names(chans) = c("133Cs_133Cs", "134Xe_134X", "138Ba_138Ba", "153Eu_SDHA", "155Gd_TOMM20", 
"158Gd_158Gd", "160Gd_NDUFB8", "161Dy_OSCP", "164Dy_GRIM19", 
"166Er_VDAC1", "168Er_Cox4", "172Yb_MTCO1", "174Yb_UqCRC2", "176Yb_Dystrophin", 
"191Ir_191Ir-DNA1", "193Ir_193Ir-DNA2", "195Pt_195P", "208Pb_208Pb", 
"80ArAr_80ArAr", "89Y_89Y")

draw$channel=chans[draw$channel] 

write.table(draw,"../pilot_data.txt",sep="\t",quote=FALSE,row.names=FALSE)