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

# Add additional package repositories
RUN add-apt-repository ppa:ubuntu-toolchain-r/test

# Install APT packages
RUN apt-get update && apt-get install -y \
    git wget curl file \
    zsh fonts-powerline
# command line packages
RUN apt-get update && apt-get install -y \
    cmake \
    gcc-11 g++-11 \
    perl cpanminus \
    python3 python3-dev python3-pip
# coding languages
RUN apt-get update && apt-get install -y \
    libboost-all-dev libspdlog-dev \
    locales tzdata
# library depdencencies

# Install Catch2 (v2.13.7-1)
WORKDIR "/tmp"
RUN wget -O catch2_2.13.7-1_amd64.deb http://ftp.us.debian.org/debian/pool/main/c/catch2/catch2_2.13.7-1_amd64.deb
RUN apt-get install -y ./catch2_2.13.7-1_amd64.deb
RUN rm catch2_2.13.7-1_amd64.deb

# Install Intel MKL via Intel OneAPI
WORKDIR "/tmp"
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
RUN apt-get update && apt-get install -y intel-basekit

# Download oh-my-zsh and make it the default terminal
RUN wget https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh -O - | zsh || true
RUN chsh -s $(which zsh)


# Copy the folders that are not ignored into the Docker image
WORKDIR "/"
COPY . /bodies-in-potential-flow/

# Install eigen3 (git-repo)
WORKDIR "/bodies-in-potential-flow/include/Eigen"
RUN mkdir build && cd build
RUN cmake .. && make install


# Encoding
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Change terminal timezone to PDT 
ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN sudo dpkg-reconfigure -f noninteractive tzdata

# Reduce image size
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get autoclean && apt-get autoremove


# Python requirements
WORKDIR "/bodies-in-potential-flow/requirements/Python"
RUN pip3 install --upgrade pip && \
    pip3 install --upgrade virtualenv
RUN pip3 install -r requirements.txt

# Perl requirements
WORKDIR "/bodies-in-potential-flow/requirements/Perl"
RUN cpanm --installdeps .


# Intel OneAPI setvars.sh environment variables
ENV ONEAPI_ROOT='/opt/intel/oneapi'
ENV MKLROOT='/opt/intel/oneapi/mkl/2021.3.0'
ENV TBBROOT='/opt/intel/oneapi/tbb/2021.3.0/env/..'
ENV CMPLR_ROOT='/opt/intel/oneapi/compiler/2021.3.0'
ENV SETVARS_VARS_PATH='/opt/intel/oneapi/vtune/latest/env/vars.sh'
ENV CMAKE_PREFIX_PATH='/opt/intel/oneapi/vpl/2021.4.0:/opt/intel/oneapi/tbb/2021.3.0/env/..:/opt/intel/oneapi/dal/2021.3.0'
ENV CPATH='/opt/intel/oneapi/vpl/2021.4.0/include:/opt/intel/oneapi/tbb/2021.3.0/env/../include:/opt/intel/oneapi/mpi/2021.3.0//include:/opt/intel/oneapi/mkl/2021.3.0/include:/opt/intel/oneapi/ippcp/2021.3.0/include:/opt/intel/oneapi/ipp/2021.3.0/include:/opt/intel/oneapi/dpl/2021.4.0/linux/include:/opt/intel/oneapi/dnnl/2021.3.0/cpu_dpcpp_gpu_dpcpp/lib:/opt/intel/oneapi/dev-utilities/2021.3.0/include:/opt/intel/oneapi/dal/2021.3.0/include:/opt/intel/oneapi/compiler/2021.3.0/linux/include:/opt/intel/oneapi/ccl/2021.3.0/include/cpu_gpu_dpcpp'
ENV LD_LIBRARY_PATH='/opt/intel/oneapi/vpl/2021.4.0/lib:/opt/intel/oneapi/tbb/2021.3.0/env/../lib/intel64/gcc4.8:/opt/intel/oneapi/mpi/2021.3.0//libfabric/lib:/opt/intel/oneapi/mpi/2021.3.0//lib/release:/opt/intel/oneapi/mpi/2021.3.0//lib:/opt/intel/oneapi/mkl/2021.3.0/lib/intel64:/opt/intel/oneapi/ippcp/2021.3.0/lib/intel64:/opt/intel/oneapi/ipp/2021.3.0/lib/intel64:/opt/intel/oneapi/dnnl/2021.3.0/cpu_dpcpp_gpu_dpcpp/lib:/opt/intel/oneapi/debugger/10.1.2/gdb/intel64/lib:/opt/intel/oneapi/debugger/10.1.2/libipt/intel64/lib:/opt/intel/oneapi/debugger/10.1.2/dep/lib:/opt/intel/oneapi/dal/2021.3.0/lib/intel64:/opt/intel/oneapi/compiler/2021.3.0/linux/lib:/opt/intel/oneapi/compiler/2021.3.0/linux/lib/x64:/opt/intel/oneapi/compiler/2021.3.0/linux/lib/emu:/opt/intel/oneapi/compiler/2021.3.0/linux/lib/oclfpga/host/linux64/lib:/opt/intel/oneapi/compiler/2021.3.0/linux/lib/oclfpga/linux64/lib:/opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/ccl/2021.3.0/lib/cpu_gpu_dpcpp'
ENV LIBRARY_PATH='/opt/intel/oneapi/vpl/2021.4.0/lib:/opt/intel/oneapi/tbb/2021.3.0/env/../lib/intel64/gcc4.8:/opt/intel/oneapi/mpi/2021.3.0//libfabric/lib:/opt/intel/oneapi/mpi/2021.3.0//lib/release:/opt/intel/oneapi/mpi/2021.3.0//lib:/opt/intel/oneapi/mkl/2021.3.0/lib/intel64:/opt/intel/oneapi/ippcp/2021.3.0/lib/intel64:/opt/intel/oneapi/ipp/2021.3.0/lib/intel64:/opt/intel/oneapi/dnnl/2021.3.0/cpu_dpcpp_gpu_dpcpp/lib:/opt/intel/oneapi/dal/2021.3.0/lib/intel64:/opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/compiler/2021.3.0/linux/lib:/opt/intel/oneapi/ccl/2021.3.0/lib/cpu_gpu_dpcpp'
ENV PATH='/opt/intel/oneapi/vtune/2021.5.0/bin64:/opt/intel/oneapi/vpl/2021.4.0/bin:/opt/intel/oneapi/mpi/2021.3.0//libfabric/bin:/opt/intel/oneapi/mpi/2021.3.0//bin:/opt/intel/oneapi/mkl/2021.3.0/bin/intel64:/opt/intel/oneapi/intelpython/latest/bin:/opt/intel/oneapi/intelpython/latest/condabin:/opt/intel/oneapi/dev-utilities/2021.3.0/bin:/opt/intel/oneapi/debugger/10.1.2/gdb/intel64/bin:/opt/intel/oneapi/compiler/2021.3.0/linux/lib/oclfpga/llvm/aocl-bin:/opt/intel/oneapi/compiler/2021.3.0/linux/lib/oclfpga/bin:/opt/intel/oneapi/compiler/2021.3.0/linux/bin/intel64:/opt/intel/oneapi/compiler/2021.3.0/linux/bin:/opt/intel/oneapi/advisor/2021.3.0/bin64:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin'


# REVIEW: Other options of commands to run in Docker container
WORKDIR "/bodies-in-potential-flow"
ENTRYPOINT [ "perl", "scripts/collinear-swimmer-wall.pl" ]
# CMD ["zsh"] # launches zsh terminal to test run commands after build