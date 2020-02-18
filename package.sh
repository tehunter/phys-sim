#!/usr/bin/env bash
# TODO: Ignore git files properly.
NAME=$1

git archive --format=zip -o "${NAME}" HEAD
git submodule --quiet foreach "cd \$toplevel; echo $path; zip -ru ${NAME} \$path -x \*.git\*"
