build:
	(export R_LIBS=/usr/lib/R/library:/usr/local/lib/R/site-library; R CMD build --no-build-vignettes ~/github/trena)
	docker build -t pshannon/trenabase ./

bash:
	docker run -ti --rm pshannon/trenabase bash

