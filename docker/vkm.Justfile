#!/usr/bin/env false

for s in /install/lib/python3.5/site-packages/*.so; do
  gosu user ln -sf "${s}" /venv/src-*/lib/python3.5/site-packages/
done

# check for pipenv lock file
if [ ! -s "${VKM_SOURCE_DIR}/Pipfile.lock" ]; then
  pipenv lock
fi

function caseify()
{
  local cmd="${1}"
  shift 1
  case "${cmd}" in

    test) # test
      pipenv run bash /src/scripts/test.bsh
      ;;

    enter) # enter
      pipenv shell
      ;;

    *) # Default: run command
      exec ${cmd} "${@}"
      ;;
  esac
}
