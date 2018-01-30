parallel -j 1 './profile.R --plate={1} --dim=3000 --rdensity=0.1 --core=2' :::: ../input/processed_plates.txt
