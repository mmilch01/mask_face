docker tag registry.nrg.wustl.edu/docker/nrg-repo/facemasking:1.21 \
    mmilch01/facemasking:1.21
docker login docker.io
docker push mmilch01/facemasking:1.21

