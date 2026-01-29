#!/bin/bash

uv version --frozen --bump patch
version="v$(uv version | cut -d" " -f2)"
gh release create "${version}" --notes "release ${version}"
git pull
