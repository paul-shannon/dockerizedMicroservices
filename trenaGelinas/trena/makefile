build:
	docker build -t testgelinastrena .

bash:
	docker run -ti --rm testgelinastrena bash

bashWithVolume:
	docker volume create myTest
	docker run -ti --rm -v myTest:/home/trena/work testgelinastrena bash
