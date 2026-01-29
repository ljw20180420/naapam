#!/bin/bash

uv version --bump patch
git tag -a "v$(uv version | cut -d" " -f2)" -m "v$(uv version | cut -d" " -f2)"
git push --tags
