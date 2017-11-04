build:
	docker-compose build

up:
	docker-compose --no-ansi up

clean:
	- docker-compose rm --force --stop --no-ansi
	- docker images | grep pyr | awk '{print $$3}' | xargs docker rmi --force

