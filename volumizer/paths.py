"""
Define project paths.
"""

from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]

C_CODE_DIR = ROOT_DIR / "src"
PACKAGE_DIR = ROOT_DIR / "volumizer"
TEST_DIR = ROOT_DIR / "tests"
