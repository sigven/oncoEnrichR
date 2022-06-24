echo "Build the Docker Image"
TAG=`date "+%Y%m%d"`
docker build --no-cache -t sigven/oncoenrichr:$TAG --rm=true .

