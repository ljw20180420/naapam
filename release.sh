#!/bin/bash

uv version --frozen --bump patch
version="v$(uv version | cut -d" " -f2)"
git commit -am "release ${version}"
git push
gh release create "${version}" --notes "release ${version}"
git pull
