FROM vsiri/recipe:gosu as gosu
FROM vsiri/recipe:tini as tini
FROM vsiri/recipe:vsi as vsi
FROM vsiri/recipe:pipenv as pipenv

FROM debian:stretch

SHELL ["/usr/bin/env", "bash", "-euxvc"]

# install base packages
RUN build_deps="wget ca-certificates"; \
    apt-get update; \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        ${build_deps} python3; \
    DEBIAN_FRONTEND=noninteractive apt-get purge -y --autoremove ${build_deps}; \
    rm -rf /var/lib/apt/lists/*

# copy from recipes
COPY --from=gosu /usr/local/bin/gosu /usr/local/bin/gosu
COPY --from=tini /usr/local/bin/tini /usr/local/bin/tini
COPY --from=vsi /vsi /vsi
COPY --from=pipenv /tmp/pipenv /tmp/pipenv

# process recipes
RUN chmod u+s /usr/local/bin/gosu
RUN python3 /tmp/pipenv/get-pipenv; rm -r /tmp/pipenv

# PIPENV setup
ENV WORKON_HOME=/venv \
    PIPENV_PIPFILE=/src/Pipfile \
    # Needed for pipenv shell \
    PYENV_SHELL=/bin/bash \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8
COPY Pipfile Pipfile.lock /src/

# GDAL specific hacks
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal \
    C_INCLUDE_PATH=/usr/include/gdal

# install pipenv packages (numpy, GDAL with python bindings)
RUN build_deps="libgdal-dev python3-dev g++"; \
    apt-get update; \
    # Install build dependencies and runtime dependencies
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      gdal-bin ${build_deps}; \
    cd /src; \
    # Build up venv
    if [ ! -s Pipfile.lock ]; then \
      pipenv lock; \
    fi; \
    # install numpy first (using version specified in Pipfile.lock)
    pipenv install --ignore-pipfile numpy; \
    # install all remaining items from Pipfile.lock
    pipenv sync; \
    # Cleanup
    apt-get purge -y --auto-remove ${build_deps}; \
    rm -r /var/lib/apt/lists/* /src/*

ENV LD_LIBRARY_PATH=/install/lib \
    PATH="/install/bin:${PATH}"

# entrypoint setup
COPY docker/common_entrypoint.bsh docker/vkm_entrypoint.bsh /
RUN chmod 755 /common_entrypoint.bsh /vkm_entrypoint.bsh
ENTRYPOINT ["/usr/local/bin/tini", "/usr/bin/env", "bash", "/vkm_entrypoint.bsh"]

# default command
CMD ["enter"]