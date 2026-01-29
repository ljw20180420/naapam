#!/bin/bash

uv version --bump patch
git tag -a $(uv version) -m $(uv version)
git push --tags
