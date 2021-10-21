#!/usr/bin/env bash

# Install dependnencies via vcpkg
vcpkg/./bootstrap-vcpkg.sh -disableMetrics
vcpkg/./vcpkg install

# NOTE: Using Assuming DEBIAN linux environment
# Install OneAPI
# @SOURCE: https://software.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-linux/top/installation/install-using-package-managers/apt.html

# use wget to fetch the Intel repository public key
cd /tmp
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
# add to your apt sources keyring so that archives signed with this key will be trusted.
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
# remove the public key
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB

# configure the APT client to use Intel's repository
echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list

# Install new packages
sudo apt update
sudo apt install -y intel-basekit
sudo apt install -y intel-hpckit
