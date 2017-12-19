build:
	docker build -t trena.mef2c .

bash:
	docker run -ti --rm trena.mef2c bash

shared:
	docker volume create myTest
	docker run -ti --rm -v myTest:/home/trena/sharedData testgelinastrena bash
