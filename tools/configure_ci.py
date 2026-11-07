"""Configure native CI with the active Python interpreter and external GoogleTest."""
import pathlib,subprocess,sys
subprocess.run(['cmake','-S','.','-B','build','-DCMAKE_BUILD_TYPE=Release',
    '-DSCE_GTEST_SOURCE='+str(pathlib.Path('.deps/googletest').resolve()),
    '-DSCE_REFERENCE_PYTHON='+sys.executable],check=True)
