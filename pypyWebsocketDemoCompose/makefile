buildServer:
	docker build -t zmqserver -f ./Dockerfile.server .

buildClient:
	docker build -t zmqclient -f ./Dockerfile.client .

runServer:
	docker run -ti --rm zmqserver

runClient:
	docker run -ti --rm zmqclient

bashServer:
	docker run -ti --rm zmqserver bash

bashClient:
	docker run -ti --rm zmqclient bash




