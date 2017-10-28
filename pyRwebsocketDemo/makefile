build.py:
	docker build -t pywebsocketclient -f ./Dockerfile.pyClient .

run.py:
	docker run --add-host="dockerhost:192.168.65.1" -it -P pywebsocketclient bash

