#!/usr/bin/env false

function caseify()
{
  local cmd="${1}"
  shift 1
  case "${cmd}" in

    compile) # compile
      exec bash "${VKM_SOURCE_DIR}/scripts/compile.bsh"
      ;;

    *) # Default: run command
      exec ${cmd} "${@}"
      ;;
  esac
}
