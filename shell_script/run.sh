input=$1
Rscript -e "rmarkdown::render('$input')"
