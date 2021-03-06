This README is incomplete, but will have additions made to it as the project develops.  

# TL;DR

Install dependencies:

- **git** [instructions](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- **docker** [instructions](https://docs.docker.com/compose/install/)
- **docker-compose** at least v1.22, [instructions](https://docs.docker.com/compose/install/)
- **JUST** at least v0.0.8, [instructions](https://github.com/VisionSystemsInc/just)

Clone repo & setup VK-metrics system
```
cd <PARENT_DIRECTORY> # where <PARENT_DIRECTORY> is a directory of your choice
git clone --recursive https://github.com/VisionSystemsInc/vkmetrics.git
cd vkmetrics
just setup
```

# VK-Metrics

This project aims to provide metrics for the IARPA CORE3D program.  Functionality includes:

* Constructing 3D ground truth models from gridded LiDAR and manual building annotation

* Measuring similiarty between 3D ground truth models and performer models

## Dependencies

This project is designed to be run within a docker container, meaning the dependencies should be mitigated. Minimal requirements:

- [git](https://git-scm.com/downloads) version control (how you should be downloading this repo)
- [Docker](https://www.docker.com/get-docker) platform independent containers
- [Docker compose](https://docs.docker.com/compose/install/) manages multiple docker applications (v1.20 or newer)
- [J.U.S.T.](https://github.com/VisionSystemsInc/just/releases) simplified tasking (v0.0.8 or newer)

### Submodules

This project consists of a number of submodules that are themselves separate repositories. Make sure to get all the submodules by using the `--recursive` flag.

```
git clone --recursive git@bitbucket.org:visionsystemsinc/vkm.git
```

Submodules include:

- **Openmesh** Library for generic and efficient representation and manipulation of polygonal meshes. More info at https://www.openmesh.org/
- **VXL** The Vision something Library. contains a wide array of different computer vision, 3D processing and mathematics libraries for use in computer vision and related applications. See http://vxl.sourceforge.net/ for more information.
- **pybind11** pybind11 is a lightweight header-only library that exposes C++ types in Python and vice versa, mainly to create Python bindings of existing C++ code. See https://pybind11.readthedocs.io/en/stable/ for more information.  VXL and VK-Metric functionality are accessible in python using pybind11 bindings.

### Install (first time instructions):

Install may take some time (possibly 10 minutes or more depending on hardware).

```
just setup
```

### Update

Updating the entire system can involve many steps. After a pull or merge, don't forget to update.

```
just sync
```

### Ground-truth model construction

Run the following command to construct a ground truth model.
```
just truth -g <TRUTH_INPUT_JSON> -o <TRUTH_MODEL_DIR>
```

**TRUTH_MODEL_DIR** is the directory for constructed ground truth.

**TRUTH_INPUT_JSON** is a JSON file describing the ground truth inputs.
All support files must reside in the same directory as the TRUTH_JSON file.
The JSON file is expected to contain the following fields:

- **name**: ground truth name
- **region**: "regions" file generated by Kitware annotation
- **type**: "types" file generated by Kitware annotation
- **dsm**: ground truth LiDAR geotiff
- **permim** optional perimeter files (names/files)

```
{
  "name": "site_ground_truth",
  "region": "CORE3D-Phase1a.kw18.regions",
  "type": "CORE3D-Phase1a.kw18.types",
  "dsm": "LiDAR.tif",
  "perim": [
    {
      "name": "building",
      "file": "building_perimeter.txt"
    }
  ]
}
```


### Metrics

Run the following command to compare performer model to ground truth model.
```
just metrics -i <INPUT_JSON> -g <TRUTH_MODEL_JSON> -o <METRIC_OUTPUT_DIR>
```

**TRUTH_MODEL_JSON** is a JSON file created by the `just truth ...` command in the TRUTH_MODEL_DIR directory.

**METRIC_OUTPUT_DIR** is the directory to store metric results.

**INPUT_JSON** is a JSON file describing the performer model.
All support files must reside in the same directory as the INPUT_JSON file.
The model is expected to consist only of top-level componenets, such as building
roofs or top-most bridge decks.
The JSON file is expected to contain the following fields:

- **name**: model name
- **coordinate_system**: model file coordinate system (currently only accepts LVCS)
    - **type**: coordinate system type (e.g., LVCS)
    - **parameters**: coordinate system description
- **models**: performer model names/files

```
{
  "name": "site_performer_model",
  "coordinate_system": {
    "type": "LVCS",
    "parameters": "wgs84 meters degrees 10.0 11.0 120.0 0.0 0.0 0 0 0"
  },
  "models": [
    {
      "name": "building",
      "file": "building.obj"
    }
  ]
}
```


## Misc.

### Building

We recommend using bash. This is the default shell for macOS and many flavors of Linux.
On windows installing git installs git bash by default.

### Elevated permissions

`gosu` is installed and accessible by the `user` account. A `sudo` function is
available to simulate basic sudo functionality.

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

