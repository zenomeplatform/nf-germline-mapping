# nf-germline-snv

Nextflow pipeline for germline variant calling.

Test examples:

Using public test data:
```
nextflow run . -profile test
```

Using private test data:
```
nextflow run . -profile test2,yandex --accessKey '<accessKey>' --secretKey '<secretKey>'
```
where `<accessKey>` and `<secretKey>` are credentials to private s3 bucket s3://zenome-ngs-data .

