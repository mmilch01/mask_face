#!/bin/bash

echo "building Docker image for Facemasking"
echo docker build . -t registry.nrg.wustl.edu/docker/nrg-repo/facemasking 
docker build . -t registry.nrg.wustl.edu/docker/nrg-repo/facemasking 
