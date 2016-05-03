#!/bin/bash
# Remove any unneaded files from stacks_workflow
# Launch from the base `stacks_workflow` directory

rm -rf .git
rm .gitignore
rm -r 99-documentation
rm -r 00-scripts/deprecated
rm -r 00-scripts/test_stacks_workflow
rm CHANGELOG LICENSE MANUAL.html README.md TODO.md
find . | grep \.gitignore | xargs rm
