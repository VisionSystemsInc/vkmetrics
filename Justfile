#!/usr/bin/env bash
set -euE

export CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd)"
source_environment_files "${CWD}/vkm.env"

#Import things like Docker-compose helper function
source "${VSI_COMMON_DIR}/linux/just_docker_functions.bsh"

function caseify()
{
  local just_arg=$1
  shift 1

  case ${just_arg} in

    ### Docker build tasks
    build) #Build all necessary docker images
      (justify build recipes gosu ninja cmake)
      (justify build main)
      ;;
    build_main) # Build maintenance docker images
      Docker-compose -f "${VKM_CWD}/docker-compose.yml" build ${@+"${@}"}
      extra_args+=$#
      ;;

    ### Code compilation
    compile) #Build the source code
      Docker-compose -f "${VKM_CWD}/docker-compose.yml" run compiler compile
      ;;

    ### Run tasks
    run) # Start a container using the vkm image
      Just-docker-compose -f "${VKM_CWD}/docker-compose.yml" run vkm ${@+"${@}"}
      extra_args+=$#
      ;;
    run_compiler) # Start a container using the compiler image
      Just-docker-compose run compiler ${@+"${@}"}
      extra_args+=$#
      ;;

    truth) #Run ground truth routine
      (source "${VKM_CWD}/scripts/task.bsh" truth "${@}")
      extra_args+=$#
      ;;

    metrics) #Run metrics routine
      (source "${VKM_CWD}/scripts/task.bsh" metrics "${@}")
      extra_args+=$#
      ;;

    ### Environment synchronize
    sync) # Sync the environment so that everyone is on the same version. This
          # includes building the docker images, building the source code, etc.
      (justify git_submodule-update) # For those users who don't remember!
      (justify build)
      (justify clean install) || :
      (justify compile)
      (justify run pipenv sync)
      ;;

    sync_quiet) #Same as sync, but quieter
      echo "==Updating submodules=="
      (justify git_submodule-update) > /dev/null
      echo "==Rebuilding docker images=="
      (justify build) > /dev/null
      echo "==Cleaning up old install directory=="
      (justify clean install) &>/dev/null || :
      echo "==Incremental build and install=="
      (justify compile) > /dev/null
      echo "==Incremental build and install=="
      (justify run pipenv sync)
      echo "Sync done. Database still running."
      ;;

    ### Testing tasks
    test) #Run unit tests
      (justify run test)
      ;;

    ### Clean tasks (e.g., remove volumes)
    clean_all) # Delete all local volumes
      ask_question "Are you sure you want to remove all local volumes?" n
      (justify clean_compile clean_install clean_venv)
      ;;
    clean_compile) # Delete only build artifacts volume
      Docker volume rm "${COMPOSE_PROJECT_NAME}_vkm-build"
      ;;
    clean_install) # Delete only install volume
      Docker volume rm "${COMPOSE_PROJECT_NAME}_vkm-bin"
      ;;
    clean_venv) # Delete only virtual env volume. This provides a \
                # clean slate for the virtual env.
      Docker volume rm "${COMPOSE_PROJECT_NAME}_vkm-venv"
      ;;

    ### Debugging tasks
    d) # Execute arbitrary docker command in Just environment
      Docker ${@+"${@}"}
      extra_args+=$#
      ;;
    dc) # Execute arbitrary docker-compose command in Just environment
      Docker-compose ${@+"${@}"}
      extra_args+=$#
      ;;
    list_volumes) # List docker volumes relevant to your environment
      Docker volume ls --filter "label=com.docker.compose.project=${COMPOSE_PROJECT_NAME}"
      ;;
    printenv) #print local environment variables
      printenv | grep -E "VKM|COMPOSE|MSYS"
      ;;
    *)
      defaultify "${just_arg}" ${@+"${@}"}
      ;;
  esac
}
