build:
	docker build -t trena ./

build.1:
	docker build -t trena -f ./Dockerfile.1 .

build.2:
	docker build -t trena -f ./Dockerfile.2 .

run:
	docker run -ti --rm trena

bash:
	docker run -ti --rm trena bash

buildTrena:
	docker build -t trenarunnable -f ./Dockerfile.trena .

runTrena:
	docker run -ti --rm trenarunnable bash


