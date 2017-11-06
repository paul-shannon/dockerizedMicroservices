build:
	docker-compose build

up:
	docker-compose --no-ansi up

clean:
	- docker-compose rm --force --stop
	- docker images | grep jupyrwebsocketcustomwidgetsdemocompose | awk '{print $$3}' | xargs docker rmi --force

