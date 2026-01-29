#!/bin/bash

uv version --dump patch
git tag -a $(uv version) -m $(uv version)
git push --tags
