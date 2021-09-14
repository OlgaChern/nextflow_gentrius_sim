# Set the base image to ubuntu
FROM ubuntu

# File Author / Maintainer
MAINTAINER Olga Chernomor <o.chernomor@gmail.com>

# Install basics
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential \
        zlib1g-dev \
	clang-tidy \
#       libhdf5-dev \
	llvm \
	gcc \
	g++ \
        cmake \
	make \
#       autoconf \
	git \
	wget \
	ca-certificates

# Install llvm packages
#RUN bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"

# Software folder
RUN mkdir /software

# Install eigen library
WORKDIR /usr/include/
RUN git clone https://gitlab.com/libeigen/eigen.git


WORKDIR /usr/include/
# Install boost
#RUN wget https://boostorg.jfrog.io/ui/native/main/release/1.77.0/source/boost_1_77_0.tar.gz/download
RUN wget https://sourceforge.net/projects/boost/files/boost/1.77.0/boost_1_77_0.tar.gz
RUN mkdir -p /usr/include/boost && tar xvf boost_1_77_0.tar.gz -C /usr/include/boost

# Install IQ-TREE, branch with Gentrius
WORKDIR /software/
RUN git clone https://github.com/iqtree/iqtree2.git
WORKDIR /software/iqtree2/
RUN git checkout terragen
RUN mkdir /software/iqtree2/build
WORKDIR /software/iqtree2/build
RUN cmake -DEIGEN3_INCLUDE_DIR=/usr/include/eigen -DBoost_INCLUDE_DIR=/usr/include/boost ..
#RUN cmake ..
RUN make 

WORKDIR /

