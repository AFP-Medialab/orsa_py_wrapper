FROM gcc:latest as builder

LABEL Maintainer="Medialab AFP"

RUN apt-get update && apt-get install -y cmake
COPY ./src/ /app/src
WORKDIR /app/build

RUN cmake ../src && make


FROM python:3.7-slim as compile
RUN apt-get update && apt-get install -y  build-essential libjpeg-dev libpng-dev libtiff-dev 
COPY ./requirements.txt /requirements.txt
RUN pip install  -r requirements.txt

#FROM python:3.7-slim
RUN mkdir -p app/out
COPY ./*.py /app
COPY ./data/ /app/data
#COPY --from=compile /root/.cache /root/.cache
COPY --from=builder /app/build/demo/*.so /app/build/demo/
COPY --from=builder /app/build/demo/all.h /app/build/demo/
#COPY ./requirements.txt /requirements.txt
#RUN pip install  -r requirements.txt
ENTRYPOINT ["/bin/sh", "-c"]
CMD ["/bin/bash"]



####
#RUN apt-get update && apt-get install -y build-essential libjpeg-dev libpng-dev libtiff-dev 