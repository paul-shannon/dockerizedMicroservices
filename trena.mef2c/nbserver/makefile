build:
	docker build -t trenaeqtlserver .

bash:
	docker run -ti --rm trenaeqtlserver bash

shared:
	docker volume create myTest
	docker run -ti --rm -v myTest:/home/jovyan/sharedData trenaeqtlserver bash

