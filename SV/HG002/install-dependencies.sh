## Update R version
## Added "deb https://cloud.r-project.org/bin/linux/debian stretch-cran35/" to /etc/apt/sources.list
sudo apt-get update
sudo apt-get install r-base

## Typical dependencies used by useful R packages
sudo apt-get install libxml2-dev libssl-dev libmariadbclient-dev libcurl4-openssl-dev

## Install R packages. In R:
# install.packages(c('devtools', "BiocManager"))
# BiocManager::install("jmonlong/sveval")
