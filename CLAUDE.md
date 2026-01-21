# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## AI Agent Development Guide

This document provides essential guidance for AI agents working on this repository. It covers tooling, conventions, and workflows needed to contribute effectively.

## How to use this document

### When to read this file
- First time working on this repository
- Before making any code changes or commits
- When unsure about code conventions or tooling

## Before every commit
- Ensure all code has type annotations
- Add sphinx-style docstrings (NO types, NO leading articles)
- Run checks: `prek -a`
- Pre-commit hooks will run automatically and must pass

## When in doubt
- Check code conventions section below
- Run individual tools to identify issues
- Ask user for clarification on ambiguous requirements

## Code conventions

### Docstrings

Required for: all public modules, classes, functions, and methods

#### Format: Sphinx-style
1. Do not place type information in docstrings - use type annotations only
2. Do not use leading articles in parameter, return, and error descriptions "a", "an", or "the"

##### Section order
1. Args
2. Returns
3. Raises
4. Yields
5. Examples
6. Note
7. Warning

#### Example
```python
def process_data(input_data: list[str], threshold: int = 10) -> dict[str, int]:
    """
    Process input data and return summary statistics.

    :param input_data: strings to process
    :param threshold: minimum count threshold for inclusion
    :return: Mapping of categories to counts

    >>> process_data(["a", "b"], 5)
    {"valid": 2}
    """
```

#### Incorrect example (do not do this!)
```python
def process_data(input_data: list[str], threshold: int = 10) -> dict[str, int]:
    """
    Process input data and return summary statistics.

    :param input_data: (list[str]): A list of strings to process.  # ❌ Has type and article
    :param threshold: (int): A minimum count threshold.            # ❌ Has type and article
    :return: dict[str, int]: A dictionary mapping categories.      # ❌ Has type and article
    """
```

### Type annotations

Requirements
- All functions must have complete type annotations
- Use modern syntax: `list[str]`, `dict[str, int]` (not `List[str]`, `Dict[str, int]`)
- Use `|` for union types (Python 3.10+): `str | None`
- Import types from `typing` only when necessary (prefer built-ins)

#### Verification
```bash
mypy .
```

### Code formatting

Via ruff
- Line length: 160
- Indentation: 4 spaces (no tabs except Makefiles)
- Encoding: UTF-8
- Line endings: LF (Unix-style)
- Trailing newline: required
- No trailing whitespace

### Naming conventions
- Functions/methods: snake_case
- Variables: snake_case
- Constants: UPPER_SNAKE_CASE
- Classes: PascalCase
- Modules: snake_case
- Private attributes/methods: _leading_underscore

### Imports
- Absolute imports preferred
- Group imports: standard library, third-party, local
- No wildcard imports (`from module import *`) except in `__init__.py`
- Import sorting handled by ruff (isort)

### Error handling
- Use specific exceptions (ValueError, OSError, etc.) rather than generic Exception
- Avoid bare `except:` clauses; catch specific exceptions
- Use context managers (`with` statements) for resource management
- Log errors appropriately using the `logging` module
- Raise custom exceptions for domain-specific errors

### General style
- Use f-strings for string formatting (Python 3.6+)
- Prefer list/dict comprehensions over loops when appropriate
- Use `pathlib.Path` for file operations instead of `os.path`
- Avoid global variables; use dependency injection
- Write readable code; prefer explicit over implicit

## Repository overview

Purpose: STructured JSON Atom/Molecule Encoding Scheme (stjames) - a Pydantic-based schema library for passing molecular and calculation data between components of the Rowan computational chemistry platform. Provides validation and intelligent default selection for computational chemistry settings.

Structure:
- `stjames/` - source code (flat layout)
- `stjames/workflows/` - high-level computational chemistry workflow definitions
- `stjames/data/` - element data and isotope information
- `stjames/atomium_stjames/` - PDB/mmCIF parsing utilities
- `tests/` - test suite
- `.github/workflows/` - CI/CD configuration

Python Version: >=3.11

Key configuration files:
- `pyproject.toml` - Project metadata, dependencies, all tool configuration
- `.pre-commit-config.yaml` - Prek hook configuration
- `.coveragerc` - Test coverage settings
- `.editorconfig` - Editor formatting settings

## Architecture

### Core Model Hierarchy

All models inherit from `Base` (in `base.py`), which extends `pydantic.BaseModel` with numpy array coercion.

**Key base classes:**
- `Base` - Pydantic BaseModel with numpy handling
- `LowercaseStrEnum` - Case/hyphen/underscore-insensitive string enum
- `UniqueList` - Annotated list type that validates uniqueness

### Main Data Models

- **`Molecule`** (`molecule.py`) - Molecular structure with atoms, charge, multiplicity, and computed properties (energy, gradient, vibrational modes, thermochemistry)
- **`Atom`** (`atom.py`) - Single atom with atomic number and position
- **`Calculation`** (`calculation.py`) - Wraps molecules, tasks, and settings for a computation
- **`Settings`** (`settings.py`) - Computational settings (method, basis set, engine, SCF/optimization parameters)

### Workflows (`stjames/workflows/`)

Workflows define high-level computational chemistry tasks. They inherit from workflow base classes in `workflow.py`:
- `Workflow` - Base for all workflows
- `MoleculeWorkflow` - Operations on a single molecule (takes `initial_molecule`, `mode`)
- `SMILESWorkflow` - Operations on a SMILES string (takes `initial_smiles`)
- `BatchSMILESWorkflow` - Operations on multiple SMILES
- `FASTAWorkflow` - Operations on biological sequences
- `ProteinStructureWorkflow` - Operations on protein PDB structures

`WORKFLOW_MAPPING` in `workflows/__init__.py` maps workflow names to classes.

### Enums and Settings

- **`Method`** (`method.py`) - Computational methods (HF, DFT functionals, NNPs, semiempirical, force fields)
- **`Engine`** (`engine.py`) - Computation backends (Psi4, PySCF, xTB, MACE, etc.)
- **`Task`** (`task.py`) - Calculation tasks (ENERGY, OPTIMIZE, GRADIENT, CHARGE, etc.)
- **`Mode`** (`mode.py`) - Accuracy levels (RECKLESS, RAPID, CAREFUL, METICULOUS, DEBUG)
- **`BasisSet`** (`basis_set.py`) - Atomic basis sets
- **`Correction`** (`correction.py`) - Dispersion corrections (D3BJ, D4, etc.)
- **`Solvent`/`SolventSettings`** (`solvent.py`) - Implicit solvent models

### Data Flow Pattern

```
Molecule/SMILES → Workflow → Settings → Calculation → [external engine] → Molecule (with results)
```

Workflows set up appropriate `Settings` with method, basis set, and parameters. `Calculation` wraps the input molecules, settings, and tasks. Results are stored back in `Molecule` properties (energy, gradient, charges, etc.).

## Essential commands

```bash
# Setup
uv sync                         # Install dependencies
prek install                    # Install git hooks

# Code quality
ruff format .                   # Format code
ruff check .                    # Lint code
mypy .                          # Type check (uses pydantic plugin)

# Testing
pytest                          # Run tests
pytest --cov                    # Run tests with coverage
pytest -k "pattern"             # Run tests matching pattern
pytest -v                       # Verbose output

# Prek (pre-commit replacement)
prek -a                         # Run all hooks
prek run <hook-id>              # Run specific hook

# Git workflow
git add .
git commit -m "feat: message"   # Hooks run automatically

# Package management
uv add <package>                # Add a package
uv add --dev <package>          # Add a package to dev
uv lock                         # Check the lockfile matches the pyproject.toml (and update if different)
uv update                       # Update all pacakges in the lockfile
uv tree                         # Print the dependencies tree
```

## Development environment

Install dependencies:
```bash
uv sync
```

Install prek hooks:
```bash
prek install
```

Verify environment:
```bash
python --version
prek -a
```

## Claude Code integration

This project includes Claude Code configuration for AI-assisted development.

### Configuration files

- `.claude/settings.local.json` - Claude permissions and hooks

### Auto-formatting hooks

Claude is configured with PostToolUse hooks that run automatically after Edit or Write operations:
1. `ruff format .` - Formats all code
2. `ruff check . --fix` - Applies auto-fixable linting corrections

Effect: When Claude modifies files, formatting and basic lints are applied immediately, reducing pre-commit friction.

To modify permissions, edit `.claude/settings.local.json`.

### Workflow impact

1. File edits trigger automatic formatting - no manual `ruff format` needed
2. Pre-commit checks still run on commit - hooks are complementary, not redundant
3. Permissions reduce interruptions for common development commands

## Direnv integration

This project uses direnv for automatic virtual environment activation.

### Behavior

When entering the project directory:
1. Checks if `uv.lock` has changed
2. Runs `uv sync --frozen` to sync dependencies
3. Activates the virtual environment automatically

### Setup

Direnv was configured during project initialization. If you encounter "direnv: error" messages:

```bash
direnv allow .
```

This grants permission for direnv to run the `.envrc` script in this directory.

## Code quality tools

### Ruff (formatting & linting)

Enabled rule categories:
- `B` - bugbear (common bugs and design problems)
- `E`/`W` - pycodestyle (PEP 8 style errors and warnings)
- `F` - pyflakes (logical errors)
- `I` - isort (import sorting)
- `N` - pep8-naming (naming conventions)
- `C4` - comprehensions (list/dict/set comprehension improvements)
- `PL` - pylint (code quality and error detection)
- `PT` - pytest-style (pytest best practices)
- `PIE` - misc lints (miscellaneous improvements)
- `PYI` - flake8-pyi (stub file best practices)
- `TID` - tidy imports (import hygiene)
- `TCH` - type-checking imports (TYPE_CHECKING block enforcement)
- `RUF` - Ruff-specific rules
- `RSE` - flake8-raise (exception raising improvements)
- `ICN001` - unconventional import aliases

Note: `D` (pydocstyle) is currently disabled

Ignored rules (globally):
- `N801` - CapWords convention
- `N802` - Arguments should be lowercase
- `N803` - Invalid argument name
- `N806` - Non-lowercase variable in function (allows PascalCase variables)
- `N815` - Mixed case attribute in class
- `E741` - Ambiguous variable name
- `PLR0911` - Too many return statements
- `PLR0912` - Too many branches
- `PLR0913` - Too many arguments to function call
- `PLR0914` - Too many local variables
- `PLR0915` - Too many statements
- `PLR1702` - Too many nested blocks
- `PLR2004` - Magic value used in comparison
- `PT006` - Pytest parameterize names wrong type
- `PT011` - Pytest raises too broad
- `PT013` - Incorrect import of pytest
- `TID252` - Prefer absolute imports from parent
- `RUF002` - Docstring contains ambiguous dash
- `RUF003` - Ambiguous dash
- `RUF012` - Mutable class defaults
- `RUF028` - Invalid suppression comment

Per-file ignores:
- `__init__.py`:
  - `F401` - Unused imports (common for `__all__` exports)
  - `F403` - `from module import *` (acceptable in `__init__.py`)

Configuration: `pyproject.toml` under `[tool.ruff]` and `[tool.ruff.lint]`

### Mypy (type checking)

Configuration:
- Uses pydantic plugin (`pydantic.mypy`)
- `strict = true` enabled
- `warn_unused_ignores = true`

Requirements:
- All functions must have complete type annotations
- Modern syntax required (e.g., `list[str]` not `List[str]`)

### Pytest (testing)

Configuration:
- Test paths:
  - `tests/` - Unit and integration tests
  - `stjames/` - Source code (for doctests)
- Doctests: Enabled automatically via `--doctest-modules` flag
- Doctest normalization: `NORMALIZE_WHITESPACE` applied to all doctests (allows flexible spacing in examples)
- Default addopts: `--doctest-modules --durations=5 -m 'not regression'`
- Markers:
  - `smoke` - sanity tests to reveal simple failures
  - `regression` - tests to make sure bugs stay closed (excluded by default)

Both test files in `tests/` and docstring examples in source code are automatically discovered and run

## Testing

### Test structure

- Location: `tests/` directory
- Naming: test files must start with `test_`
- Doctests: automatically discovered in source code

### Doctests

Include examples in docstrings:
```python
def add(a: int, b: int) -> int:
    """
    Add two integers.

    :param a: first integer
    :param b: second integer
    :return: sum of a and b

    >>> add(2, 3)
    5
    >>> add(-1, 1)
    0
    """
    return a + b
```

### Coverage configuration

- Source: `stjames`
- Omitted: `stjames/__main__.py`
- Excluded lines: `pragma: no cover`, `__repr__`, `if self.debug`, `raise AssertionError`, `raise NotImplementedError`, `if 0:`, `if __name__ == "__main__":`, ellipsis-only lines (`...`)

## Prek hooks

Configuration: `.pre-commit-config.yaml`

Hooks enabled (from pre-commit-hooks v6.0.0 and local):

1. check-yaml - Validate YAML syntax
2. check-toml - Validate TOML syntax
3. end-of-file-fixer - Ensure single newline at EOF
4. trailing-whitespace - Remove trailing whitespace
5. ruff-format - Format Python code (`uv run ruff format .`)
6. ruff-check - Lint Python code (`uv run ruff check .`)
7. mypy - Type check Python code (`uv run mypy .`)
8. pytest - Run test suite (`uv run pytest`)

Behavior:
- Stages: `pre-commit` and `pre-push`
- All hooks use `pass_filenames: false` (process entire project)
- Hooks must pass before commit succeeds
- Some hooks auto-fix issues (formatting, whitespace)

Manual execution:
```bash
prek -a                 # All hooks
prek run ruff-check     # Specific hook
```

## CI/CD

File: `.github/workflows/test.yml`

Triggers:
- All pull requests
- Pushes to `master` branch

Test matrix:

- Python versions: 3.11, 3.12, 3.13
- Operating systems: ubuntu-latest
- Total combinations: 3

Checks performed:
1. `uv run ruff format .` - format check
2. `uv run ruff check .` - lint check
3. `uv run mypy .` - type check
4. `uv run pytest --cov --cov-report=xml` - tests with coverage

Ensure local checks pass before pushing to avoid CI failures.

## Troubleshooting

### Pre-commit hook failures

Formatting issues:
- Usually auto-fixed by ruff
- Re-stage files: `git add .`
- Try committing again

Lint issues:
- Read error message for specific rule
- Fix; if `# noqa: <rule>` if absolutely necessary, ask before adding

Type issues:
- Add missing type annotations
- Fix type mismatches
- Use `mypy .` to verify locally

Test failures:
- Fix failing tests or code
- Run `pytest -v` for detailed output
- Run specific test: `pytest tests/test_file.py::test_name`

## Commit guidelines

Format: conventional commits recommended
- `feat:` - new features
- `fix:` - bug fixes
- `docs:` - documentation changes
- `test:` - test changes
- `refactor:` - code refactoring
- `chore:` - maintenance tasks

### Example
```bash
git commit -m "feat: add user authentication

- Implement JWT token generation
- Add login/logout endpoints
- Include comprehensive tests"
```

## Additional resources

- uv documentation: https://docs.astral.sh/uv
- ruff documentation: https://docs.astral.sh/ruff
- mypy documentation: https://mypy.readthedocs.io
- pydantic documentation: https://docs.pydantic.dev
- pytest documentation: https://docs.pytest.org
- prek documentation: https://prek.j178.dev
