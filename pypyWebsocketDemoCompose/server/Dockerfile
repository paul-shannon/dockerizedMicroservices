FROM python:3
RUN mkdir /code
WORKDIR /code
ADD requirements.txt /code/
ADD server.py /code/
RUN pip install -r requirements.txt
ADD . /code/
CMD ["python", "-i", "server.py"]
