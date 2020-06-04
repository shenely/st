from setuptools import setup, find_packages

setup(name="st",
      version="0.1",
      author="shenely",
      description="A star tracker algorithm",
      python_requires=">=3.6",
      install_requires=["numpy",
                        "scipy",
                        "matplotlib"],
      packages=find_packages())


