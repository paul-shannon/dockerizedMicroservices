build:
	docker-compose build
up:
	docker-compose --no-ansi up

clean:
	docker-compose rm --force --stop --no-ansi
	- docker rmi `docker images | grep pypy | awk '{print $1}'` --force
