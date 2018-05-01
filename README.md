This README is incomplete, but will have additions made to it as the project develops.

# Core3D Bottom up processing

This project aims to extract information from 3D models created from high resolution Digital Elevation Models. The goal is compact high level explanation of the underlying geometry. This is accomplished through a set of techniques to fit primitive geometric shapes to point clouds.

## Dependencies

This project is designed to be run within a docker container, meaning the dependencies should be mitigated. There are some requirements however.

[git](https://git-scm.com/downloads) version control (how you should be downloading this repo)

[Docker](https://www.docker.com/get-docker) which provides a platform independent container for everything we build

[Docker compose 1.19 or newer](https://docs.docker.com/compose/install/) which manages running multiple docker applications

### submodules

This project consists of a number of submodules that are themselves separate repositories. Make sure to get all the submodules by using

```
git clone --recursive git@bitbucket.org:visionsystemsinc/vkm.git
```

### Install (first time instructions):

```
source setup.env

just install

```

Install may take some time (possibly 10 minutes or more depending on hardware).


## Building

We recommend using bash. This is the default shell for macOS and many flavors of Linux. On windows installing git installs git bash by default.


## Elevated permissions

`gosu` is installed and accessible by the `user` account. A `sudo` function is
available to simulate basic sudo functionality.

## Updating

Updating the entire core3d system can involve many steps. These steps are all run
for you when you run

```
just sync
```

After a pull or merge, don't forget to run `just sync`

## Module Breakdown

### External modules

These modules come from external projects and may be updated periodically from those external sources.

#### Openmesh

Library for generic and efficient representation and manipulation of polygonal meshes. More info at https://www.openmesh.org/

#### VSI_common

A set of common tools for use across VSI applications. This includes the 'just' command and other scripts

#### VXL

The Vision something Library. contains a wide array of different computer vision, 3D processing and mathematics libraries for use in computer vision and related applications. See http://vxl.sourceforge.net/ for more information.

## Adding 3rd party software

### Python packages

Python packages used by the containers, are maintained by `Pipfile` and
`Pipfile.lock`, described [here](https://github.com/pypa/pipfile). This is the
new replacement for `requirements.txt` files allowing us to track our python
dependencies.

The python virtual environment is stored in `/venv` in a directory starting
with `src-`. By default, the virtual environment is stored in a docker volume
to remain persistent, allowing you to make local changes for development. If
for any reason, you need to reset this directory, removing the docker volume
will cause it to be reset to factory settings upon starting the core3d container.
This can be done by `just clean venv`

- **Adding a python package for everyone** - `pipenv install scipy`. This will
add scipy as a tracker package in the `Pipfile`, and all of its dependencies
(*numpy*, in this case) will be tracked in the `Pipfile.lock`. The exact
versions of all the packages are store in `Pipfile.lock`, which is why it is
important to `git add` and `commit` these back to the repository.
- **Manually adding a python package for everyone** - Packages can be added
to `Pipfile` by editing the `[packages]` section. Once you have made your changes
run `pipenv update` to apply those changes to your virtual env.
- **Locally adding a python package** - (not recommended) `pip install scipy`.
This will install the package **scipy** and all of its dependencies (*numpy*).
These packages are not tracked and will not be there if you clean your virtual
env, nor will anyone else see that you are using that packages.
- **Checking for the latest versions of tracked python packages** - Run
`pipenv update`, and `Pipfile.lock` is updated with all the newest versions.
- **Updating your virtualenv** - Use `pipenv sync` to updated your packages to
match what is in the core3d repository, `Pipfile`. This is done for you
every time you run `just sync`
- **Check to see what local changes are untracked** - Not a perfect solution,
but running `pipenv clean --dry-run` will list packages you have locally install
that are not tracked.
- **Remove local packages** - While `just clean venv` will do this, in order to
not completely restart the container, `pipenv clean` will remove all local
packages. NOTE: This does not run `pipenv sync` also, you will have to call that
manually.

Using `pipenv` and these tools, you should be able to both customize your local
virtual environment, and shared python packages with everyone on the team.

## TODOs

Integrate preprocessing and DEM generation

Full pipeline automation

More complete unit tests
