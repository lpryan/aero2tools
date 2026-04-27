import os
import subprocess
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext

class MakeBuild(build_ext):
    
    def run(self):
        # build C++ core
        subprocess.check_call(["make"])
        
        # tell setuptools where .so files are
        for filename in os.listdir("aero2tools/cpp"):
            if filename.endswith(".so"):
                self.copy_file(
                    os.path.join("aero2tools/cpp", filename),
                    os.path.join(self.build_lib, "cpp", filename),
                )

setup(
    name = "aero2tools",
    version = "4.0.1b",
    packages=find_packages(),
    cmdclass={"build_ext": MakeBuild}
)