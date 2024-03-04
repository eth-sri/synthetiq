FROM python:3.10

WORKDIR /usr/synthetiq

RUN apt-get update && apt-get install -y \
    gcc-11 \
    make

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 100

COPY . .
RUN rm -rf data/input data/gates

RUN make clean
RUN make all

RUN pip install -r requirements.txt

COPY data/input data/input
COPY data/gates data/gates