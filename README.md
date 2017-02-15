# AppScanGEO
Source code for ScanGEO Shiny App

Requirements:
15 GB disk space
8 GB RAM
R installed (https://www.r-project.org/)
Internet connection

ScanGEO installation instructions:
1) Download files from this repository into your local target directory.
2) In your local target directory, create two additional subdirectories named "softfiles" and "results".
3) unzip www.zip (this will create a www subdirectory with the files for the app documentation and logo).
4) Download the GEOmetadb data base into your local directory by running the following commands in R:
> library(GEOmetadb)
> getSQLiteFile(destdir = getwd(), destfile = "GEOmetadb.sqlite.gz")
5) When scanning a particular GDS for the first time, the scan will take longer, as the soft file needs to be downloaded once.
