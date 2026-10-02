FROM gcr.io/oss-fuzz-base/base-builder
RUN apt-get update && apt-get install -y cmake ninja-build python3
RUN git clone --recurse-submodules https://github.com/loumalouomega/meshioplusplus.git $SRC/meshioplusplus
WORKDIR $SRC/meshioplusplus
COPY build.sh $SRC/build.sh
