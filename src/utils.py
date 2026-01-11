import os
import pathlib
import re


def infer_cas(file: os.PathLike) -> str:
    name = pathlib.Path(os.fspath(file)).name.lower()
    if (
        name.startswith("a2")
        or name.startswith("a7")
        or name.startswith("d2")
        or name.startswith("h73")
        or name.startswith("h80")
        or name.startswith("p79")
        or name.startswith("p125")
        or name.startswith("q101")
        or name.startswith("q102")
        or name.startswith("q192")
    ):
        return "spycas9"
    if name.startswith("36t") or name.startswith("b2") or name.startswith("x"):
        return "spymac"
    if name.startswith("i10t") or name.startswith("i83"):
        return "ispymac"
    if name.startswith("wt"):
        return "control"
    raise Exception("unknown cas")


def infer_polq(file: os.PathLike) -> str:
    name = pathlib.Path(os.fspath(file)).name.lower()
    if (
        name.startswith("a2")
        or name.startswith("a7")
        or name.startswith("d2")
        or name.startswith("36t")
        or name.startswith("b2")
        or name.startswith("x")
        or name.startswith("i10t")
        or name.startswith("i83")
        or name.startswith("wt")
    ):
        return "wt"
    if name.startswith("h73") or name.startswith("h80"):
        return "hel"
    if name.startswith("p79") or name.startswith("p125"):
        return "pol"
    if name.startswith("q101") or name.startswith("q102") or name.startswith("q192"):
        return "ko"
    raise Exception("unknown polq")


def infer_chip(file: os.PathLike) -> str:
    name = pathlib.Path(os.fspath(file)).name.lower()
    return re.search(r"-(a1|a2|a3|g1n|g2n|g3n)-", name).group(1)


def infer_time(file: os.PathLike) -> str:
    name = pathlib.Path(os.fspath(file)).name.lower()
    return re.search(r"-(?:a1|a2|a3|g1n|g2n|g3n)-(\d)", name).group(1)


def infer_sample(file: os.PathLike) -> str:
    name = pathlib.Path(os.fspath(file)).name.lower()
    return re.search(r"^(.+?)-", name).group(1)
