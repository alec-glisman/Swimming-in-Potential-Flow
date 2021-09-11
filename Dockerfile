# Base image
FROM ubuntu:20.04

# Label directives
LABEL maintainer=alec.glisman@gmail.com
LABEL project=bodies-in-potential-flow
LABEL version=0.1
LABEL environment=dev

# General environment directives
ENV DEBIAN_FRONTEND=noninteractive 


# Update packages on base image and install sudo
RUN apt-get update --fix-missing && \
    apt-get upgrade -y && \
    apt-get install -y sudo build-essential software-properties-common

# command line packages
RUN apt-get update && apt-get install -y \
    git wget curl file \
    zsh fonts-powerline


# Add newer gcc package repositories
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
# Add intel package repositories

# Add intel repositories
WORKDIR "/tmp"
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list


# coding languages
RUN apt-get update && apt-get install -y \
    cmake \
    gcc-11 g++-11 \
    perl cpanminus \
    python3 python3-dev python3-pip

# library depdencencies
RUN apt-get update && apt-get install -y \
    libboost-all-dev \
    locales tzdata


# Encoding
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Change terminal timezone to PDT 
ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN sudo dpkg-reconfigure -f noninteractive tzdata

# Install Intel MKL via Intel OneAPI
# RUN apt-get update && apt-get install -y intel-basekit
RUN pip3 install mkl

# Reduce image size
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get autoclean && apt-get autoremove


# Copy the folders that are not ignored into the Docker image
WORKDIR "/"
COPY . /bodies-in-potential-flow/

# Python requirements
WORKDIR "/bodies-in-potential-flow/requirements/Python"
RUN pip3 install --upgrade pip && \
    pip3 install --upgrade virtualenv
RUN pip3 install -r requirements.txt

# Perl requirements
WORKDIR "/bodies-in-potential-flow/requirements/Perl"
RUN cpanm --installdeps .


# Install eigen3 (commit 8ce341caf2947e4b5ac4580c20254ae7d828b009 from git repo)
WORKDIR "/bodies-in-potential-flow/include/Eigen"
RUN mkdir build
WORKDIR "/bodies-in-potential-flow/include/Eigen/build"
RUN cmake .. && make install

# Install spdlog (v1.9.1 from git repo)
WORKDIR "/bodies-in-potential-flow/include/spdlog"
RUN mkdir build
WORKDIR "/bodies-in-potential-flow/include/spdlog/build"
RUN cmake .. && make install

# Install catch2 (v2.13.7 from git repo)
WORKDIR "/bodies-in-potential-flow/include/Catch2"
RUN mkdir build
WORKDIR "/bodies-in-potential-flow/include/Catch2/build"
RUN cmake .. && make install


# REVIEW: Other options of commands to run in Docker container
WORKDIR "/bodies-in-potential-flow"
ENTRYPOINT [ "perl", "scripts/collinear-swimmer-wall.pl" ]   # starts simulations from perl script
# CMD ["zsh"] # launches zsh terminal to test run commands after build
