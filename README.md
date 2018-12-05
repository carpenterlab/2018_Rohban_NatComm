# Capturing single-cell heterogeneity via data fusion improves morphological profiling #

## Abstract ##
Single-cell resolution technologies warrant computational methods that capture cell heterogeneity while allowing efficient comparisons of populations. Here, we summarize cell populations by adding features’ measures of dispersion and covariances to population averages, in the context of morphological profiling. We find that data fusion is critical for these metrics to improve results over the prior state-of-the-art, providing ~30% better performance in tasks including predicting a compound’s mechanism of action (MoA) and a gene’s pathway.

## Prerequisites ##
* Mac OS X
* R Ver. 3.3.3 
* Following R packages: dplyr 0.7.4
magrittr 1.5
foreach 1.4.4
stringr 1.2.0
readr 1.1.1
doParallel 1.0.11
SNFtool 2.2
ggplot2 2.2.1
Matrix 1.2-8
htmlTable 1.6
readbulk 1.1.0
cytominer 0.1.0.9000 (https://github.com/cytomining/cytominer)
* aws command line interface (https://docs.aws.amazon.com/cli/latest/userguide/cli-install-macos.html) configured to access `cellpainting-datasets` s3 bucket. 
* Package installation time is about an hour on a typical PC. 
* Note : For each dataset, create a separate clone of the repository. Then, `cd code`.

* Note : TA-ORF-BBBC037-Rohban is the smallest dataset consisting of only around 5 plates, so can also be used for the demo purposes. Each plate takes on average between 2 to 3 hours to get processed on a normal PC. Bioactives-BBBC022-Gustafsdottir and CDRPBIO-BBBC036-Bray consist of 20 and 55 plates, respectively. 

## Creating median+MAD profiles ##
* Bioactives-BBBC022-Gustafsdottir : `parallel -j 1 './profile_trad.R --name=Bioactives-BBBC022-Gustafsdottir --batch=BBBC022_2013 --plate={1} --operation="median+mad" --col="Metadata_broad_sample" --value="DMSO" --cores=2 --feats="../input/feature_list_BBBC022.txt"' :::: ../input/processed_plates_BBBC022.txt`
* TA-ORF-BBBC037-Rohban : `parallel -j 1 './profile_trad.R --name=TA-ORF-BBBC037-Rohban --batch=SIGMA2_Pilot_2013_10_11 --plate={1} --operation="median+mad" --col="Metadata_ASSAY_WELL_ROLE" --value="Untreated" --cores=2 --feats="../input/feature_list.txt"' :::: ../input/processed_plates_TA.txt`
* CDRPBIO-BBBC036-Bray : `parallel -j 1 './profile_trad.R --name=CDRPBIO-BBBC036-Bray --batch=CDRP --plate={1} --operation="median+mad" --col="Metadata_broad_sample" --value="DMSO" --cores=2 --feats="../input/feature_list.txt"' :::: ../input/processed_plates_CDRP_bio.txt`

## Creating cov. profiles ##
* Bioactives-BBBC022-Gustafsdottir : 
``` 
rm ../input/random_projection_unified.rds

mv ../input/random_projection_unified_BBBC022.rds ../input/random_projection_unified.rds

parallel -j 1 './profile.R --name=Bioactives-BBBC022-Gustafsdottir --batch=BBBC022_2013 --plate={1} --dim=3000 --rdensity=0.1 --core=2 --col=Metadata_broad_sample --value="DMSO" --feats="../input/feature_list_BBBC022.txt"' :::: ../input/processed_plates_BBBC022.txt 

```
* TA-ORF-BBBC037-Rohban : `parallel -j 1 './profile.R --name=TA-ORF-BBBC037-Rohban --batch=SIGMA2_Pilot_2013_10_11 --plate={1} --dim=3000 --rdensity=0.1 --core=2 --col=Metadata_ASSAY_WELL_ROLE --value="Untreated" --feats="../input/feature_list.txt"' :::: ../input/processed_plates_TA.txt`
* CDRPBIO-BBBC036-Bray : `parallel -j 1 './profile.R --name=CDRPBIO-BBBC036-Bray --batch=CDRP --plate={1} --dim=3000 --rdensity=0.1 --core=2 --col=Metadata_broad_sample --value="DMSO" --feats="../input/feature_list.txt"' :::: ../input/processed_plates_CDRP_bio.txt`

## Creating the treatment correlation matrices ##
* Bioactives-BBBC022-Gustafsdottir :
``` 
./evaluate.R -m "median" -p "../input/processed_plates_BBBC022.txt" -e ../input/metadata_BBBC022.csv -f "../input/feature_list_BBBC022.txt"

./evaluate.R -m "mad" -p "../input/processed_plates_BBBC022.txt" -e ../input/metadata_BBBC022.csv -f "../input/feature_list_BBBC022.txt"

./evaluate.R -m "cov" -p "../input/processed_plates_BBBC022.txt" -e ../input/metadata_BBBC022.csv -f "../input/feature_list_BBBC022.txt"

./evaluate.R -m "median+mad" -p "../input/processed_plates_BBBC022.txt" -e ../input/metadata_BBBC022.csv -f "../input/feature_list_BBBC022.txt"

```
* TA-ORF-BBBC037-Rohban :
``` 
./evaluate.R -m "median" -p "../input/processed_plates_TA.txt" -e ../input/metadata_TA.csv -f "../input/feature_list.txt"

./evaluate.R -m "mad" -p "../input/processed_plates_TA.txt" -e ../input/metadata_TA.csv -f "../input/feature_list.txt"

./evaluate.R -m "cov" -p "../input/processed_plates_TA.txt" -e ../input/metadata_TA.csv -f "../input/feature_list.txt"

./evaluate.R -m "median+mad" -p "../input/processed_plates_TA.txt" -e ../input/metadata_TA.csv -f "../input/feature_list.txt"

```
* CDRPBIO-BBBC036-Bray :
``` 
./evaluate.R -m "median" -p "../input/processed_plates_CDRP_bio.txt" -e ../input/metadata_CDRP.csv -f "../input/feature_list.txt"

./evaluate.R -m "mad" -p "../input/processed_plates_CDRP_bio.txt" -e ../input/metadata_CDRP.csv -f "../input/feature_list.txt"

./evaluate.R -m "cov" -p "../input/processed_plates_CDRP_bio.txt" -e ../input/metadata_CDRP.csv -f "../input/feature_list.txt"

./evaluate.R -m "median+mad" -p "../input/processed_plates_CDRP_bio.txt" -e ../input/metadata_CDRP.csv -f "../input/feature_list.txt"

```

## Generating Fig. 1A (enrichment comparison plot) ##
* Run `./compare_mean_cov.R -p chemical` for Bioactives-BBBC022-Gustafsdottir and CDRPBIO-BBBC036-Bray
* Run `./compare_mean_cov.R -p genetic` for TA-ORF-BBBC037-Rohban

## Generating Fig. 1B (similarity graphs for an MOA) ##
* Run `sub_corr_plot.R` for CDRPBIO-BBBC036-Bray
