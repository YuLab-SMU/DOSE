md5:
	Rscript -e 'source("generate-md5.R")'


hpo:
	Rscript -e 'source("create-HPO.R")'

html:
	Rscript -e 'rmarkdown::render("index.rmd")'
