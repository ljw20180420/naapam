#!/bin/bash

uv version --frozen --bump patch
version="v$(uv version | cut -d" " -f2)"
git add pyproject.toml
git commit -m "release ${version}"
gh release create "${version}" --notes "release ${version}"
git pull
