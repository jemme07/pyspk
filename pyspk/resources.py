"""Loading of packaged resources (data files)."""

from __future__ import annotations

import importlib.resources as _resources


def read_text(package: str, name: str, *, encoding: str = "utf-8") -> str:
    """Read a text resource bundled in the package.

    Args:
        package: Package name that owns the resource (e.g. "pyspk").
        name: Resource file name.
        encoding: Text encoding.

    Returns:
        The decoded text contents.

    Raises:
        FileNotFoundError: If the resource does not exist.
    """

    return _resources.files(package).joinpath(name).read_text(encoding=encoding)
