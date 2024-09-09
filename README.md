# stjames

[![pypi](https://img.shields.io/pypi/v/stjames.svg)](https://pypi.python.org/pypi/stjames)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v1.json)](https://github.com/charliermarsh/ruff)

*STructured JSON Atom/Molecule Encoding Scheme*

<img src='img/james_icon.jpg' width=350>

This is the Rowan schema for passing molecule/calculation data back and forth between different parts of the software.

This is not intended to be run as a standalone library: it's basically just a big composite Pydantic model which does some validation and intelligent default selection.
(A benefit of doing validation on the client side is that it's transparent to the end user—you can see all of the settings that the calculation will use.)

## Installation

To install, ensure you have Python 3.11 or newer. Then run:

```
pip install stjames
```

For bug reports, please use the Issues tab (above).

## New Releases

To release a new version, update the version number in ``pyproject.toml`` (using semantic versioning) and then issue a release with a matching version number. Github Actions will take care of the rest!

*Corin Wagen, 2023*
