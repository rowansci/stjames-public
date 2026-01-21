# stjames

[![pypi](https://img.shields.io/pypi/v/stjames.svg)](https://pypi.python.org/pypi/stjames)
[![License](https://img.shields.io/github/license/rowansci/stjames-public)](https://github.com/rowansci/stjames-public/blob/master/LICENSE)
[![Powered by: uv](https://img.shields.io/badge/-uv-purple)](https://docs.astral.sh/uv)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

*STructured JSON Atom/Molecule Encoding Scheme*

<img src='img/james_icon.jpg' width=350>

This is the Rowan schema for passing molecule/calculation data back and forth between different parts of the software.

This is not intended to be run as a standalone library: it's basically just a big composite Pydantic model which does some validation and intelligent default selection.
(A benefit of doing validation on the client side is that it's transparent to the end user—you can see all of the settings that the calculation will use.)
