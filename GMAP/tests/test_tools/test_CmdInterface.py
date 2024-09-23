"""
Tests all the functions/classes/methods in the file:
src/tools/CLibLoader.py.

Missing tests:

(@ August 2nd '24):
  113  (1 missed statements)

- We don't actually invoke a program here (test_GEM does do that) (113)
"""

# Standard library imports
import pytest
import subprocess

# Local imports
import GMAP
from GMAP.src.tools.Exceptions import GmapAttributeError
import GMAP.src.tools.CmdInterface as GM_CI
import GMAP.src.tools.PrintTools as GM_PT


def test_cmd_interface(capsys):
    callcommand = ["GMAP", "nothing", "nothing"]
    with pytest.raises(GmapAttributeError, match="SU_GM_1$"):
        GM_CI.cmd_interface(callcommand)

    callcommand = ["GMAP"]  # == GMAP help help
    GM_CI.cmd_interface(callcommand)
    captured = capsys.readouterr()
    assert captured.out.endswith(GM_PT.prettifier(GMAP.__doc__) + "\n")

    callcommand = ["GMAP", "HeLp"]  # == GMAP help help
    GM_CI.cmd_interface(callcommand)
    captured = capsys.readouterr()
    assert captured.out.endswith(GM_PT.prettifier(GMAP.__doc__) + "\n")

    callcommand = ["GMAP", "GEM"]  # == GMAP GEM help
    GM_CI.cmd_interface(callcommand)
    captured = capsys.readouterr()
    assert captured.out.endswith(GM_PT.prettifier(GMAP.GEM.__doc__) + "\n")

    callcommand = ["GMAP", "help", "GEM"]
    GM_CI.cmd_interface(callcommand)
    captured = capsys.readouterr()
    assert captured.out.endswith(GM_PT.prettifier(GMAP.GEM.__doc__) + "\n")

    callcommand = ["GMAP", "help", "nothing"]
    with pytest.raises(GmapAttributeError, match="SU_GM_1$"):
        GM_CI.cmd_interface(callcommand)


def test_main():
    callcommand = ["GMAP", "nothing", "nothing"]
    captured = subprocess.run(callcommand, capture_output=True)
    stderr = captured.stderr
    # On windows the \r will be generated, not on linux/mac
    assert stderr.endswith(b"SU_GM_1\r\n") or stderr.endswith(b"SU_GM_1\n")

    callcommand = ["GMAP"]  # == GMAP help help
    captured = subprocess.run(callcommand, capture_output=True)
    # On windows
    outbytes = (GM_PT.prettifier(GMAP.__doc__) + "\n").replace(
        "\n", "\r\n").encode('utf-8')
    # On linux/mac
    outbyteslinux = (GM_PT.prettifier(GMAP.__doc__) + "\n").encode('utf-8')
    assert captured.stdout.endswith(outbytes) or captured.stdout.endswith(outbyteslinux)
