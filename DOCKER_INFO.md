# Information on Docker Integration and Data Syncing

## Build and run this project

[Source](https://unix.stackexchange.com/a/531575)

Build the project:

```[bash]
docker build --tag aglisman/bodies-in-potential-flow:latest .
```

Run the container and have volume for data output shared between VM and host:

```[bash]
docker run -it --volume "$(pwd)/Docker_output/:/bodies-in-potential-flow/data" aglisman/bodies-in-potential-flow:latest
```

**All together:**

```[bash]
docker build --tag aglisman/bodies-in-potential-flow:latest . && docker run -it --volume "$(pwd)/Docker_output/:/bodies-in-potential-flow/data" aglisman/bodies-in-potential-flow:latest
```

### Save to version control and image

```[bash]
git add .
git commit -m "Message"
git push origin master

docker push aglisman/bodies-in-potential-flow:latest
```

### Free space on local machine

```[bash]
# Remove all resources not associated with a container
docker system prune

# Remove any stopped containers and all unused images
docker system prune -a
```

## `rsync` Use between server and local machine

- Final slash in folder path will sync ALL FILE CONTENTS, NOT FILE ITSELF
- Run `rsync` as a daemon with `--daemon`

```[bash]
# rsync options source destination
# -r: recursive, recurse into directories
# -a: archive mode, archive mode allows copying files recursively and it also preserves symbolic links, file permissions, user & group ownerships and timestamps
# -v: verbose
# -h, --human-readable, output numbers in a human-readable format
# -e: specify protocol over which to transfer data (use ssh)
# --progress: display file transfer progress
# --delete: delete extraneous files from dest dirs
# --exclude: do not sync certain directories,
#       notably any C++ CMake build directories with '*build*' pattern
#       also all hidden files, potentially: .git


# Sync data FROM local TO remote, excluding any build directories
rsync -ravh --progress --delete                    \
-e "ssh -i KEY.pem.txt"                            \
/home/localuser/testfolder                         \
remoteuser@X.X.X.X:/home/remoteuser/testfolder


# Sync data FROM remote TO local (DO NOT DELETE AT LOCAL)
rsync -ravh --progress                            \
-e "ssh -i KEY.pem.txt"                           \
remoteuser@X.X.X.X:/home/remoteuser/testfolder    \
/home/localuser/testfolder



# Current build commands

# FROM local TO remote
rsync -ravh --progress --exclude '*build*' --exclude '.*' --delete     \
-e "ssh -i ~/keys/cpp-explorations-keypair.pem.txt"     \
"CLionProjects/cpp-explorations/"                       \
"ubuntu@ec2-18-144-15-61.us-west-1.compute.amazonaws.com:project"

# FROM remote TO local
rsync -ravh --progress --exclude '*build*' --exclude '.*'    \
-e "ssh -i ~/keys/cpp-explorations-keypair.pem.txt"                  \
"ubuntu@ec2-18-144-15-61.us-west-1.compute.amazonaws.com:project/"   \
"CLionProjects/cpp-explorations/"
```

## Set-up Docker and Git in a new Repository

```[bash]
# Start repo
git init
touch README.md

#  Add all project files, .gitignore, .dockerignore, and Dockerfile

# log-in to docker
docker login --username=yourhubusername --password=yourpassword



# Add docker image
docker tag local-image:tagname new-repo:tag-name
docker push new-repo:tagname

# Run docker container (-v for linking volume for data I/O)
docker build --rm -f "Dockerfile" -t repo:tagname .
docker run -it --rm -v $(pwd)/file:/file repo:tagname

# Save to version control
git add .
git commit -m "Message"

# Add remote server for git
git remote add origin {remote repository URL}
git push -u origin master

# Synchronize Changes with Docker
docker push repo:tagname



# Pull changes down from git and docker
docker pull repo
git pull

# Locally save and load docker images
docker save repo > repo.tar
docker load --input repo.tar
```
