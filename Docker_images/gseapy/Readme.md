## Build the docker image from within folder
```bash
docker build --rm -t gseapy:0.10.4 .
```

## This changes the name of the image
```bash
docker tag 218edd5fcf10 evafast1/gseapy:0.10.4
```

## this does the login
```bash
docker login -u evafast1 -p <password>
```

## this pushes to the repo
```bash
docker push evafast1/gseapy:0.10.4
```

## Launch the docker image
```bash
docker run \
--rm \
-d \
--name gseapy \
-p 8883:8888 \
-e JUPYTER_ENABLE_LAB=YES \
-v /Users/efast/:/home/jovyan/work \
gseapy:0.10.4
```