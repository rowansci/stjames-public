# stjames

*STructured JSON Atom/Molecule Encoding Scheme*

<img src='img/james_icon.jpg' width=350>

This is the Rowan schema for passing molecule/calculation data back and forth between different parts of the software.

This is not intended to be run as a standalone library: it's basically just a big composite Pydantic model which does some validation and intelligent default selection.
(A benefit of doing validation on the client side is that it's transparent to the end user—you can see all of the settings that the calculation will use.)

## Installation

To install, ensure you have Python 3.7 or newer. Then run:

```
pip install stjames
```

For bug reports, please use the Issues tab (above).

*Corin Wagen, 2023*
