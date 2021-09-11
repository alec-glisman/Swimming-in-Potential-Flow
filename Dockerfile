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

# Add newer gcc package repositories
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
# Add intel package repositories

WORKDIR "/tmp"
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list


# command line packages
RUN apt-get update && apt-get install -y \
    git wget curl file \
    zsh fonts-powerline

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

# Python requirements
WORKDIR "/bodies-in-potential-flow/requirements/Python"
RUN pip3 install --upgrade pip && \
    pip3 install --upgrade virtualenv
RUN pip3 install -r requirements.txt

# Perl requirements
WORKDIR "/bodies-in-potential-flow/requirements/Perl"
RUN cpanm --installdeps .


# Install Intel MKL via Intel OneAPI
RUN apt-get update && apt-get install -y intel-basekit

# Reduce image size
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get autoclean && apt-get autoremove


# Copy the folders that are not ignored into the Docker image
WORKDIR "/"
COPY . /bodies-in-potential-flow/

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

# Intel OneAPI setvars.sh environment variables
ENV ONEAPI_ROOT='/opt/intel/oneapi'
ENV MKLROOT='/opt/intel/oneapi/mkl/latest'
ENV TBBROOT='/opt/intel/oneapi/tbb/latest/env/..'
ENV CMPLR_ROOT='/opt/intel/oneapi/compiler/latest'
ENV SETVARS_VARS_PATH='/opt/intel/oneapi/vtune/latest/env/vars.sh'
ENV CMAKE_PREFIX_PATH='/opt/intel/oneapi/vpl/latest:/opt/intel/oneapi/tbb/latest/env/..:/opt/intel/oneapi/dal/latest'
ENV CPATH='/opt/intel/oneapi/vpl/latest/include:/opt/intel/oneapi/tbb/latest/env/../include:/opt/intel/oneapi/mpi/latest//include:/opt/intel/oneapi/mkl/latest/include:/opt/intel/oneapi/ippcp/latest/include:/opt/intel/oneapi/ipp/latest/include:/opt/intel/oneapi/dpl/latest/linux/include:/opt/intel/oneapi/dnnl/latest/cpu_dpcpp_gpu_dpcpp/lib:/opt/intel/oneapi/dev-utilities/latest/include:/opt/intel/oneapi/dal/latest/include:/opt/intel/oneapi/compiler/latest/linux/include:/opt/intel/oneapi/ccl/latest/include/cpu_gpu_dpcpp'
ENV LD_LIBRARY_PATH='/opt/intel/oneapi/vpl/latest/lib:/opt/intel/oneapi/tbb/latest/env/../lib/intel64/gcc4.8:/opt/intel/oneapi/mpi/latest//libfabric/lib:/opt/intel/oneapi/mpi/latest//lib/release:/opt/intel/oneapi/mpi/latest//lib:/opt/intel/oneapi/mkl/latest/lib/intel64:/opt/intel/oneapi/ippcp/latest/lib/intel64:/opt/intel/oneapi/ipp/latest/lib/intel64:/opt/intel/oneapi/dnnl/latest/cpu_dpcpp_gpu_dpcpp/lib:/opt/intel/oneapi/debugger/10.1.2/gdb/intel64/lib:/opt/intel/oneapi/debugger/10.1.2/libipt/intel64/lib:/opt/intel/oneapi/debugger/10.1.2/dep/lib:/opt/intel/oneapi/dal/latest/lib/intel64:/opt/intel/oneapi/compiler/latest/linux/lib:/opt/intel/oneapi/compiler/latest/linux/lib/x64:/opt/intel/oneapi/compiler/latest/linux/lib/emu:/opt/intel/oneapi/compiler/latest/linux/lib/oclfpga/host/linux64/lib:/opt/intel/oneapi/compiler/latest/linux/lib/oclfpga/linux64/lib:/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/ccl/latest/lib/cpu_gpu_dpcpp'
ENV LIBRARY_PATH='/opt/intel/oneapi/vpl/latest/lib:/opt/intel/oneapi/tbb/latest/env/../lib/intel64/gcc4.8:/opt/intel/oneapi/mpi/latest//libfabric/lib:/opt/intel/oneapi/mpi/latest//lib/release:/opt/intel/oneapi/mpi/latest//lib:/opt/intel/oneapi/mkl/latest/lib/intel64:/opt/intel/oneapi/ippcp/latest/lib/intel64:/opt/intel/oneapi/ipp/latest/lib/intel64:/opt/intel/oneapi/dnnl/latest/cpu_dpcpp_gpu_dpcpp/lib:/opt/intel/oneapi/dal/latest/lib/intel64:/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/compiler/latest/linux/lib:/opt/intel/oneapi/ccl/latest/lib/cpu_gpu_dpcpp'
ENV PATH='/opt/intel/oneapi/vtune/2021.5.0/bin64:/opt/intel/oneapi/vpl/latest/bin:/opt/intel/oneapi/mpi/latest//libfabric/bin:/opt/intel/oneapi/mpi/latest//bin:/opt/intel/oneapi/mkl/latest/bin/intel64:/opt/intel/oneapi/intelpython/latest/bin:/opt/intel/oneapi/intelpython/latest/condabin:/opt/intel/oneapi/dev-utilities/latest/bin:/opt/intel/oneapi/debugger/10.1.2/gdb/intel64/bin:/opt/intel/oneapi/compiler/latest/linux/lib/oclfpga/llvm/aocl-bin:/opt/intel/oneapi/compiler/latest/linux/lib/oclfpga/bin:/opt/intel/oneapi/compiler/latest/linux/bin/intel64:/opt/intel/oneapi/compiler/latest/linux/bin:/opt/intel/oneapi/advisor/latest/bin64:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin'


# REVIEW: Other options of commands to run in Docker container
WORKDIR "/bodies-in-potential-flow"
ENTRYPOINT [ "perl", "scripts/collinear-swimmer-wall.pl" ]   # starts simulations from perl script
# CMD ["zsh"] # launches zsh terminal to test run commands after build