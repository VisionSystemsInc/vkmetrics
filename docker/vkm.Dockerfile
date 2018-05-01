FROM vsiri/recipe:gosu as gosu

FROM debian:stretch

SHELL ["bash", "-euxvc"]

RUN apt-get update; \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        gdal-bin libgdal-dev python3-dev g++ \
        curl ca-certificates; \
    curl -L https://bootstrap.pypa.io/get-pip.py | python3; \
    # Do no add more pip packages here, use the Pipfile instead
    pip install pipenv==11.9.0; \
    # Cleanup
    apt-get purge -y --auto-remove curl ca-certificates; \
    rm -r /var/lib/apt/lists/*

# PIPENV specific settings
ENV WORKON_HOME=/venv \
    PIPENV_PIPFILE=/src/Pipfile \
    # Needed for pipenv shell \
    PYENV_SHELL=/bin/bash \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    # Gdal specific hacks
    CPLUS_INCLUDE_PATH=/usr/include/gdal \
    C_INCLUDE_PATH=/usr/include/gdal

# Comment out this ADD the first time for a project, and run "pipenv install"
# in /install to generate the Pipfiles the first time, then uncomment these
# lines
ADD Pipfile Pipfile.lock /src/

# This bootstraps the pipenv environment to pre-include packages in the image.
# The project dir is the dirname of PIPENV_PIPFILE
RUN cd /src; \
    # Build up venv
    pipenv install --skip-lock "numpy$(python3 -c "import json; print(json.load(open('Pipfile.lock', 'r'))['default']['numpy']['version'])" 2>/dev/null)"; \
    pipenv sync; \
    # Cleanup
    rm -rv /src/*

# install gdb
RUN apt-get update; \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends gdb; \
    rm -r /var/lib/apt/lists/*

COPY --from=gosu /usr/local/bin/gosu /usr/local/bin/gosu
RUN chmod 4111 /usr/local/bin/gosu

ENV LD_LIBRARY_PATH=/install/lib \
    PATH="/install/bin:${PATH}"

ADD docker/pyrc /etc/pyrc
ENV PYTHONSTARTUP=/etc/pyrc

ADD docker/vkm_entrypoint.bsh /
RUN chmod 755 /vkm_entrypoint.bsh
ENTRYPOINT ["/vkm_entrypoint.bsh"]

# Snapshot virtualenv so a fresh install contains the virtual env
VOLUME /venv

CMD ["enter"]