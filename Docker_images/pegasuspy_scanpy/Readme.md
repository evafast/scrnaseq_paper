## Build the docker image from within folder
```bash
docker build --rm -t pegasuspy_scanpy:vs1 .
```

## This changes the name of the image
```bash
docker tag 23ce8738897e evafast1/pegasuspy_scanpy:vs1
```

## this does the login
```bash
docker login -u evafast1 -p <password>
```

## this pushes to the repo
```bash
docker push evafast1/pegasuspy_scanpy:vs1
```

## Launch the docker image
```bash
docker run \
--rm \
-d \
--name demuxEM \
-p 8880:8888 \
-e JUPYTER_ENABLE_LAB=YES \
-v /Users/efast/Documents/:/home/jovyan/work \
pegasuspy_scanpy:vs1
```