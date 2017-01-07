FROM ubuntu:latest
RUN apt-get update && apt-get install -y \
  cmake \
  make \
  g++ \
  git \
  openmpi-bin libopenmpi-dev libhdf5-openmpi-dev \
  gcc
COPY . /usr/src/myapp
WORKDIR /usr/src/myapp/build
RUN cmake ..
RUN cmake --build .
RUN ctest -VV
