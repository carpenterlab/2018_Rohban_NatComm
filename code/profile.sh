#
# proj.		name								batch				plate list
# TA ORF	2011_07_13_TargetAccelerator_CancerProgram_MPG 			SIGMA2_Pilot_2013_10_11		processed_plates_TA.txt
# CDRP		2015_Bray_GigaScience						CDRP				processed_plates_CDRP_bio.txt
# Repurposing	2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad	2016_04_01_a549_48hr_batch1	processed_plates_repurposing.txt

parallel -j 1 './profile.R --name=2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad --batch=2016_04_01_a549_48hr_batch1 --plate={1} --dim=3000 --rdensity=0.1 --core=2 --col=Metadata_broad_sample --value=DMSO --feats="../input/feature_list.txt"' :::: ../input/processed_plates_repurposing.txt
