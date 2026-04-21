#!/bin/bash

uv sync
uv build
uv run cibuildwheel
