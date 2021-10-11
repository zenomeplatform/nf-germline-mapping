## Documentation for NF-Germline-SNV pipeline (In progress)



### Temporary

To run a test with differnt values of LB for same samples:
```
nextflow run . -profile test --multiqc_prealignment_by_sample false --multiqc_prealignment_all false --input testdata/TESTX/TESTX_test_data_paths_copied_for_LB.csv
```

To run a test using full hg38 reference data:
```
nextflow run . -profile test --genome hg38 --max_memory 6.GB --max_cpus 4
```
