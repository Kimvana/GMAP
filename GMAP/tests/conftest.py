"""
A basis for all tests. Things here are available to all tests!
"""

# 3rd party imports
import pytest

# local imports
import GMAP.src.tools.CodingTools as GM_CT
import GMAP.src.tools.FileHandler as GM_FH


# To report to the user that this file is present and active
@pytest.fixture(scope="session", autouse=True)
def my_fixture(request):
    capmanager = request.config.pluginmanager.getplugin("capturemanager")
    with capmanager.global_and_fixture_disabled():
        print(
            "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            "\n!!!!! Running tests using tests/conftest.py !!!!!"
            "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        )


# To prevent any created singletons leaking between tests.
# https://github.com/pytest-dev/pytest-mock/issues/100
@pytest.fixture()
def reset_singletons():
    GM_CT.Singleton._instances = {}


@pytest.fixture()
def init_files():
    GM_FH.FileLocations()


@pytest.fixture(autouse=True)
def pre_test(reset_singletons, init_files):
    pass
