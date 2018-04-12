# Capturing single-cell heterogeneity via data fusion improves morphological profiling #

## Abstract ##
Recent advances in single-cell resolution technologies warrants designing of computational methods that can capture cell heterogeneity beyond population averages. In this paper, we aim at adding measure of dispersion and feature covariances to the population average through data fusion techniques in the context of morphological profiling. In comparison to state-of-the-art, the proposed method provides substantial improvement (typically around 30%) in enrichment of treatment pairs with most similar profiles in having same Mechanisms of Action (MOA)/pathways.


## Prerequisites ##
* Mac OS X
* R Ver. 3.3.3 
* cytominer package (http://github.com/cytomining/cytominer)
* aws command line interface (https://docs.aws.amazon.com/cli/latest/userguide/cli-install-macos.html) 
* Note : For each dataset, create a separate clone of the repository. Then, `cd code`

## Creating median+MAD profiles ##
* BBBC022 : `parallel -j 1 './profile_trad.R --name=2016_12_13_Cytominer_Janssen --batch=BBBC022_2013 --plate={1} --operation="median+mad" --col="Metadata_broad_sample" --value="DMSO" --cores=2 --feats="../input/feature_list_BBBC022.txt"' :::: processed_plates_BBBC022.txt`
* TA ORF : `parallel -j 1 './profile_trad.R --name=2011_07_13_TargetAccelerator_CancerProgram_MPG --batch=SIGMA2_Pilot_2013_10_11 --plate={1} --operation="median+mad" --col="Metadata_ASSAY_WELL_ROLE" --value="Untreated" --cores=2 --feats="../input/feature_list.txt"' :::: processed_plates_TA.txt`
* CDRP : `parallel -j 1 './profile_trad.R --name=2015_Bray_GigaScience --batch=CDRP --plate={1} --operation="median+mad" --col="Metadata_broad_sample" --value="DMSO" --cores=2 --feats="../input/feature_list.txt"' :::: processed_plates_CDRP_bio.txt`

## Creating cov. profiles ##
* BBBC022 : 
``` 
rm ../input/random_projection_unified.rds

mv ../input/random_projection_unified_BBBC022.rds ../input/random_projection_unified.rds

parallel -j 1 './profile.R --name=2016_12_13_Cytominer_Janssen --batch=BBBC022_2013 --plate={1} --dim=3000 --rdensity=0.1 --core=2 --col=Metadata_broad_sample --value="DMSO" --feats="../input/feature_list_BBBC022.txt"' :::: processed_plates_BBBC022.txt 

```
* TA ORF : `parallel -j 1 './profile.R --name=2011_07_13_TargetAccelerator_CancerProgram_MPG --batch=SIGMA2_Pilot_2013_10_11 --plate={1} --dim=3000 --rdensity=0.1 --core=2 --col=Metadata_ASSAY_WELL_ROLE --value="Untreated" --feats="../input/feature_list.txt"' :::: processed_plates_TA.txt`
* CDRP : `parallel -j 1 './profile.R --name=2015_Bray_GigaScience --batch=CDRP --plate={1} --dim=3000 --rdensity=0.1 --core=2 --col=Metadata_broad_sample --value="DMSO" --feats="../input/feature_list.txt"' :::: processed_plates_CDRP_bio.txt`

## Creating the treatment correlation matrices ##
* BBBC022 :
``` 
./evaluate.R -m "median" -p "../input/processed_plates_BBBC022.txt" -e ../input/metadata_BBBC022.csv -f "../input/feature_list_BBBC022.txt"

./evaluate.R -m "mad" -p "../input/processed_plates_BBBC022.txt" -e ../input/metadata_BBBC022.csv -f "../input/feature_list_BBBC022.txt"

./evaluate.R -m "cov" -p "../input/processed_plates_BBBC022.txt" -e ../input/metadata_BBBC022.csv -f "../input/feature_list_BBBC022.txt"

./evaluate.R -m "median+mad" -p "../input/processed_plates_BBBC022.txt" -e ../input/metadata_BBBC022.csv -f "../input/feature_list_BBBC022.txt"
```
* TA ORF :
``` 
./evaluate.R -m "median" -p "../input/processed_plates_TA.txt" -e ../input/metadata_TA.csv -f "../input/feature_list.txt"

./evaluate.R -m "mad" -p "../input/processed_plates_TA.txt" -e ../input/metadata_TA.csv -f "../input/feature_list.txt"

./evaluate.R -m "cov" -p "../input/processed_plates_TA.txt" -e ../input/metadata_TA.csv -f "../input/feature_list.txt"

./evaluate.R -m "median+mad" -p "../input/processed_plates_TA.txt" -e ../input/metadata_TA.csv -f "../input/feature_list.txt"
```
* CDRP :
``` 
./evaluate.R -m "median" -p "../input/processed_plates_CDRP_bio.txt" -e ../input/metadata_CDRP.csv -f "../input/feature_list.txt"

./evaluate.R -m "mad" -p "../input/processed_plates_CDRP_bio.txt" -e ../input/metadata_CDRP.csv -f "../input/feature_list.txt"

./evaluate.R -m "cov" -p "../input/processed_plates_CDRP_bio.txt" -e ../input/metadata_CDRP.csv -f "../input/feature_list.txt"

./evaluate.R -m "median+mad" -p "../input/processed_plates_CDRP_bio.txt" -e ../input/metadata_CDRP.csv -f "../input/feature_list.txt"
```

## Generating Fig. 1 (enrichment comparison plot) ##
* Run `compare_mean_cov.R` 

