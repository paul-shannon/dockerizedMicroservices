build:
	docker build -t testgelinasnbserver .

bash:
	docker run -ti --rm testgelinasnbserver bash

shared:
	docker volume create myTest
	docker run -ti --rm -v myTest:/home/jovyan/sharedData testgelinasnbserver bash

