# Base image
FROM ubuntu:20.04

# Label directives
LABEL maintainer=alec.glisman@gmail.com
LABEL project=bodies-in-potential-flow
LABEL version=0.1
LABEL environment=dev

# Environment directives
ENV DEBIAN_FRONTEND=noninteractive 

# Update packages on base image
RUN apt install -y sudo
RUN sudo apt update
RUN sudo apt upgrade -y

# Install APT packages
RUN sudo apt install -y \
    build-essential software-properties-common \
    git wget \
    zsh cmake gcc cpanminus \
    python3-pip python-dev \
    libboost-all-dev \
    fonts-powerline \
    tzdata


# Install Python requirements
RUN pip3 install -r requirements/Python/requirements.txt

# Install Perl requirements
RUN cpanm File::Slurp
RUN cpanm File::Tee
RUN cpanm Path::Tiny
RUN cpanm Switch


# Install Intel MKL, Source: https://software.intel.com/content/www/us/en/develop/articles/installing-intel-free-libs-and-python-apt-repo.html
# now get the key:
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
# now install that key
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
# now remove the public key file exit the root shell
RUN rm GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
# Add all libraries
RUN wget https://apt.repos.intel.com/setup/intelproducts.list -O /etc/apt/sources.list.d/intelproducts.list
# Update packages
RUN sudo apt-get -y update
# Install MKL
RUN sudo apt-get install intel-mkl-2020.0-088 -y

# Install Boost, Source: https://stackoverflow.com/questions/12578499/how-to-install-boost-on-ubuntu
RUN wget -O boost_1_72_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.72.0/boost_1_72_0.tar.gz/download
RUN tar xzvf boost_1_72_0.tar.gz
WORKDIR "/boost_1_72_0"
RUN sudo ./bootstrap.sh
RUN sudo ./b2 toolset=gcc cxxflags="-std=c++11" linkflags="-std=c++11" -j 10 || true
RUN sudo ./b2 install toolset=gcc cxxflags="-std=c++11" linkflags="-std=c++11" -j 10 || true
WORKDIR "/"

# Download oh-my-zsh and make it the default terminal
RUN wget https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh -O - | zsh || true
RUN chsh -s $(which zsh)


# Change terminal timezone to PDT 
ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN sudo dpkg-reconfigure -f noninteractive tzdata


# Copy the folders that are not ignored into the Docker image
COPY . /potential-swimmer-dynamics/
WORKDIR "/potential-swimmer-dynamics"


# Run the program to build from CMAKE
ENTRYPOINT [ "perl", "scripts/collinearAsymmetric.pl" ]  

# REVIEW: Other options of commands to run in Docker container
# ENTRYPOINT [ "perl", "scripts/collinearConstrained.pl"]  
# CMD ["zsh"] # launches zsh terminal to test run commands after build
