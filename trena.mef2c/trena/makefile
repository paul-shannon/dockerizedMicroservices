build:
	docker build -t testgelinastrena .

bash:
	docker run -ti --rm testgelinastrena bash

shared:
	docker volume create myTest
	docker run -ti --rm -v myTest:/home/trena/sharedData testgelinastrena bash
