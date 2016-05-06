#!/bin/bash
# Remove any unneaded files from stacks_workflow
# Launch from the base `stacks_workflow` directory

rm -r 99-documentation
rm -r 00-scripts/utility_scripts/deprecated
rm -r 00-scripts/utility_scripts/test_stacks_workflow
rm CHANGELOG.txt LICENSE MANUAL.html README.md TODO.md
