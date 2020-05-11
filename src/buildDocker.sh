echo "Build the Docker Image"
TAG=`date "+%Y%m%d"`
docker build -t sigven/oncoenrichr:$TAG --rm=true .

