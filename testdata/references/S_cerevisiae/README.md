## File sources

### saccharomyces_cerevisiae.vcf.gz
File containing known S. cerevisiae variants
Initially obtained from http://ftp.ebi.ac.uk/ensemblgenomes/pub/current/fungi/variation/vcf/saccharomyces_cerevisiae/
Then bgzipped and tabix-indexed. Commands:
```
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/current/fungi/variation/vcf/saccharomyces_cerevisiae/saccharomyces_cerevisiae.vcf.gz
gunzip saccharomyces_cerevisiae.vcf.gz
bgzip saccharomyces_cerevisiae.vcf
tabix saccharomyces_cerevisiae.vcf.gz
```
