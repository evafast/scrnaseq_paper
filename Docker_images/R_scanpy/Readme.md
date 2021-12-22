## Build the docker image from within folder
```bash
docker build --rm -t r_scanpy:vs5 .
```

## This changes the name of the image
```bash
docker tag adb321099d9f evafast1/r_scanpy:vs5
```

## this does the login
```bash
docker login -u evafast1 -p <password>
```

## this pushes to the repo
```bash
docker push evafast1/r_scanpy:vs5
```

## Launch the docker image
```bash
docker run \
--rm \
-d \
--name r_scanpy2 \
-p 8885:8888 \
-e JUPYTER_ENABLE_LAB=YES \
-v /Users/efast/Documents/:/home/jovyan/work \
r_scanpy:vs5
```