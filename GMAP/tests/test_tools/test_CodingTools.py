"""
Tests all the functions/classes/methods in the file:
src/tools/CmdInterface.py.

Missing tests:

(@ Sept 20th '24):
  (0 missed statements)

- None!

"""

# Standard library imports
import pytest

# Local imports
import GMAP.src.tools.CodingTools as GM_CT


class TestFrozenDict:
    def test_retrieve(self):
        mydict = {2: 3, 4: 5}
        frozendict = GM_CT.FrozenDict(mydict)
        assert mydict[2] == frozendict[2]

    def test_add(self):
        mydict = {2: 3, 4: 5}
        frozendict = GM_CT.FrozenDict(mydict)
        mydict = {1: 0, 0: 1}
        with pytest.raises(TypeError):
            _ = frozendict + mydict

    def test_radd(self):
        mydict = {2: 3, 4: 5}
        frozendict = GM_CT.FrozenDict(mydict)
        mydict = {1: 0, 0: 1}
        with pytest.raises(TypeError):
            _ = mydict + frozendict

    def test_setitem(self):
        mydict = {2: 3, 4: 5}
        frozendict = GM_CT.FrozenDict(mydict)
        with pytest.raises(TypeError):
            frozendict[1] = 1

    def test_update(self):
        mydict = {2: 3, 4: 5}
        frozendict = GM_CT.FrozenDict(mydict)
        mydict = {1: 0, 0: 1}
        with pytest.raises(TypeError):
            frozendict.update(mydict)
