## Build the docker image from within folder
```bash
docker build --rm -t signac:0.2.5 .
```

## This changes the name of the image
```bash
docker tag aae86fcb0170 evafast1/signac:0.2.5
```

## this does the login
```bash
docker login -u evafast1 -p <password>
```

## this pushes to the repo
```bash
docker push evafast1/signac:0.2.5
```

## Launch the docker image
```bash
docker run \
--rm \
-d \
--name signac_old \
-p 8883:8888 \
-e JUPYTER_ENABLE_LAB=YES \
-v /Users/efast/Documents/:/home/jovyan/work \
signac:0.2.5
```