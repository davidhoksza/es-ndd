# es-ndd
Essential Sites in NDD proteins

## Deploying using Docker

Clone the ES-NDD repository
```
git clone https://github.com/davidhoksza/es-ndd
```

Enter the directory with the repo
```
cd es-ndd
```

Build the ES-NDD image 
```
docker build . -t es-ndd
```

Start a container with the ES-NDD image and make it run in the background exposing port 8080 in the container to the 8080 port on the localhost
```
docker run -dp 8080:8080 es-ndd
```

The ES-NDD app should now be available at http://localhost:8080/.
