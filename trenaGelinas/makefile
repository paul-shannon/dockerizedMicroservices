build:
	docker-compose build

up:
	docker-compose --no-ansi up

clean:
	- docker-compose rm --force --stop
	- docker images | grep trena.gelinas | awk '{print $$3}' | xargs docker rmi --force

