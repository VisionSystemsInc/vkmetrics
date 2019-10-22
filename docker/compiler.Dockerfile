FROM vsiri/recipe:gosu as gosu
FROM vsiri/recipe:tini as tini
FROM vsiri/recipe:ninja as ninja
FROM vsiri/recipe:cmake as cmake
FROM vsiri/recipe:vsi as vsi

# This image is responsible for setting up the install folder. This is accomplished
# in two pieces:
# 1. Compile and install source code - This is achieved by installing the libraries
#    and compilers needed to compile all of the core3d source code. The entrypoint
#    will call the actual build script to use at run time
# 2. Building the virtualenv - This is achieved by building the virtualenv and
#    installing packages at docker build time. And the install directory is
#    populated at docker load time, as long as the install directory is empty
#    to start with. For this requirement to be achieved, the install directory
#    is wiped at sync time

FROM debian:stretch

SHELL ["/usr/bin/env", "bash", "-euxvc"]

RUN apt-get update; \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python python3-dev python-pip libgeotiff-dev \
        gcc g++ curl bzip2 rsync unzip ca-certificates; \
    rm -rf /var/lib/apt/lists/*

# copy from recipes
COPY --from=gosu /usr/local/bin/gosu /usr/local/bin/gosu
COPY --from=tini /usr/local/bin/tini /usr/local/bin/tini
COPY --from=ninja /usr/local/bin/ninja /usr/local/bin/ninja
COPY --from=cmake /cmake /usr/local/
COPY --from=vsi /vsi /vsi

# process recipes
RUN chmod u+s /usr/local/bin/gosu

# entrypoint setup
COPY docker/compiler.Justfile /
ENTRYPOINT ["/usr/local/bin/tini", "--", "/usr/bin/env", "bash", "/vsi/linux/just_entrypoint.sh"]

# default command
CMD ["bash"]