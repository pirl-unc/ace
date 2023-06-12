#!/bin/bash
set -o errexit

find src test -name '*.py' \
  | xargs pylint \
  --errors-only
  --disable=print-statement,unsubscriptable-object,not-an-iterable,no-member


echo 'Passes pylint check'