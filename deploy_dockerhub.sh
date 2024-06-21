docker tag registry.nrg.wustl.edu/docker/nrg-repo/facemasking:latest \
    mmilch01/facemasking:1.1
docker login docker.io
docker push mmilch01/facemasking:1.1

